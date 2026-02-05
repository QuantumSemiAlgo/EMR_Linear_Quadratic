/**
 * @file petsc_diagonalizer.cpp
 * @brief Linear Solver and Matrix Analysis Interface.
 *
 * This module manages the solution of the linear system Ax = b using PETSc's
 * KSP (Krylov Subspace Methods) solvers. It also includes comprehensive tools
 * for analyzing the condition number and stability of the system matrix.
 *
 * Key Features:
 * - Configurable KSP solver setup.
 * - Condition number estimation using SLEPc (Eigenvalue analysis).
 * - Matrix norm and sparsity analysis.
 * - Diagonal scaling checks.
 * - Post-solve verification (residual check).
 * - Solution scattering and output.
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#include "global_params_1.h"
#include "prototypes_1.h"
#include <petscksp.h>
#include <slepceps.h>

// Forward declaration of the condition number analysis function
PetscErrorCode compute_condition_number(Mat A, data &dat);

/**
 * @brief Solves the linear system and performs post-processing.
 *
 * @param gmat Reference to the global matrices structure (A, rhs).
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode petsc_diagonalizer(global_matrices &gmat, data &dat) {
  PetscErrorCode ierr;
  KSP ksp;        // Linear solver context
  Vec x, tx;      // Solution vector (parallel) and gathered vector (seq)
  VecScatter ctx; // Scatter context

  // -------------------------------------------------------------------------
  // 1. Create Solution Vector
  // -------------------------------------------------------------------------
  ierr = MatCreateVecs(gmat.A, NULL, &x);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);
  CHKERRQ(ierr);

  // Check matrix and RHS norms for sanity
  PetscReal anorm, bnorm;
  ierr = MatNorm(gmat.A, NORM_INFINITY, &anorm);
  CHKERRQ(ierr);
  ierr = VecNorm(gmat.rhs, NORM_INFINITY, &bnorm);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Checking assembly: ||A||_∞ = %g, ||b||_∞ = %g\n", anorm,
                     bnorm);
  CHKERRQ(ierr);

  // -------------------------------------------------------------------------
  // 2. Condition Number Analysis
  // -------------------------------------------------------------------------
  // Perform a detailed analysis of the matrix properties before solving.
  // This helps diagnose convergence issues or accuracy loss.
  // ierr = compute_condition_number(gmat.A, dat);
  // CHKERRQ(ierr);

  // -------------------------------------------------------------------------
  // 3. Setup KSP Solver
  // -------------------------------------------------------------------------
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, gmat.A, gmat.A);
  CHKERRQ(ierr);

  // Set solver type (e.g., KSPCG, KSPGMRES) from input file
  ierr = KSPSetType(ksp, dat.kspType);
  CHKERRQ(ierr);

  // Set convergence tolerances: relative, absolute, divergence, max iterations
  ierr = KSPSetTolerances(ksp, dat.rtol, dat.atol, dat.divtol, dat.maxit);
  CHKERRQ(ierr);

  // Allow runtime customization via command line (e.g., -ksp_monitor)
  ierr = KSPSetFromOptions(ksp);
  CHKERRQ(ierr);

  // -------------------------------------------------------------------------
  // 4. Solve System
  // -------------------------------------------------------------------------
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Solving system...\n");
  CHKERRQ(ierr);
  ierr = KSPSolve(ksp, gmat.rhs, x);
  CHKERRQ(ierr);

  // -------------------------------------------------------------------------
  // 5. Check Convergence
  // -------------------------------------------------------------------------
  KSPConvergedReason reason;
  PetscInt its;
  PetscReal rnorm;

  ierr = KSPGetConvergedReason(ksp, &reason);
  CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp, &its);
  CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp, &rnorm);
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     ">>> KSP convergence:\n"
                     "    Reason:     %d (%s)\n"
                     "    Iterations: %d\n"
                     "    Residual:   %.6e\n",
                     reason, reason > 0 ? "CONVERGED" : "DIVERGED", its, rnorm);
  CHKERRQ(ierr);

  if (reason < 0) {
    ierr =
        PetscPrintf(PETSC_COMM_WORLD, "⚠ WARNING: Solver did not converge!\n");
    CHKERRQ(ierr);
  }

  // -------------------------------------------------------------------------
  // 6. Post-Solve: Enforce Periodicity in Solution
  // -------------------------------------------------------------------------
  // For annular domains, we must copy the solution from the "master" nodes
  // to the "slave" (overlaid) nodes to ensure the output vector is complete.
  if (dat.use_annular) {
    ierr =
        PetscPrintf(PETSC_COMM_WORLD, "Copying solution to overlaid DOFs...\n");
    CHKERRQ(ierr);

    int copy_count = 0;
    for (int i = 0; i < dat.nglobal; ++i) {
      if (dat.Convertnode[i] != i) {
        // This is an overlaid DOF - copy value from canonical DOF
        PetscInt overlaid_dof = i;
        PetscInt canonical_dof = dat.Convertnode[i];

        PetscScalar canonical_value;
        ierr = VecGetValues(x, 1, &canonical_dof, &canonical_value);
        CHKERRQ(ierr);
        ierr = VecSetValue(x, overlaid_dof, canonical_value, INSERT_VALUES);
        CHKERRQ(ierr);

        copy_count++;
      }
    }

    ierr = VecAssemblyBegin(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "Copied solution to %d overlaid DOFs\n", copy_count);
    CHKERRQ(ierr);
  }

  // -------------------------------------------------------------------------
  // 7. Verify Residual (Ax - b)
  // -------------------------------------------------------------------------
  Vec r;
  PetscReal res_norm;
  ierr = MatCreateVecs(gmat.A, NULL, &r);
  CHKERRQ(ierr);
  ierr = MatMult(gmat.A, x, r);
  CHKERRQ(ierr); // r = A * x
  ierr = VecAXPY(r, -1.0, gmat.rhs);
  CHKERRQ(ierr); // r = r - b
  ierr = VecNorm(r, NORM_INFINITY, &res_norm);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "||A x - b||_∞ = %g\n", res_norm);
  CHKERRQ(ierr);
  ierr = VecDestroy(&r);
  CHKERRQ(ierr);

  // -------------------------------------------------------------------------
  // 8. Output Solution
  // -------------------------------------------------------------------------
  // Gather the distributed solution vector to process 0 for writing to disk.
  ierr = VecScatterCreateToAll(x, &ctx, &tx);
  CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx, x, tx, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx, x, tx, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRQ(ierr);
  ierr = VecScatterDestroy(&ctx);
  CHKERRQ(ierr);

  // Call the solution output function (interpolates to grid)
  // ierr = solution(gmat, dat, tx);
  // CHKERRQ(ierr);

  // -------------------------------------------------------------------------
  // 9. EMR Post-Processing (Resistance Calculation)
  // -------------------------------------------------------------------------
  if (dat.Nsample_R2 > 0) { // EMR Mode Check
    PetscScalar *sol_array;
    ierr = VecGetArray(tx, &sol_array);
    CHKERRQ(ierr);

    int port4_nodes = 0;
    double v4_sum = 0.0;

    // Calculate Average Voltage on Port 4 (Probe)
    // Rank 0 has the full solution 'tx' due to scatter above
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (rank == 0) {
      for (int i = 0; i < dat.ngnodes; ++i) {
        if (dat.port_code[i] == 4) {
          v4_sum += PetscRealPart(sol_array[i]);
          port4_nodes++;
        }
      }

      double v4_avg = 0.0;
      if (port4_nodes > 0)
        v4_avg = v4_sum / port4_nodes;

      double Resistance = 0.0;
      if (dat.Io != 0.0) {
        // Resistance (Ohms) = (Volt / Amps)
        // Note: Stiffness matrix already includes dat.t, so V is already scaled
        // correctly.
        Resistance = (v4_avg / dat.Io);
      }
      dat.latest_resistance = Resistance;

      // ============================================================================
      // DIAGNOSTIC OUTPUT - Resistance Calculation
      // ============================================================================
      {
        double V4_volts = v4_avg;
        double V3_volts = 0.0; // Port 3 is reference (0V)
        double I_amps = dat.Io;
        double R_val = Resistance;

        ierr = PetscPrintf(PETSC_COMM_WORLD,
                           "\n"
                           "───────────────────────────────────────────────────"
                           "───────────────\n"
                           "RESISTANCE CALCULATION DETAILS:\n"
                           "───────────────────────────────────────────────────"
                           "───────────────\n"
                           "  Voltage at Port 4 (avg): V4 = %.12e V\n"
                           "  Voltage at Port 3 (ref): V3 = %.12e V\n"
                           "  Voltage difference:     ΔV = %.12e V\n"
                           "  Input current:          Io = %.12e A\n"
                           "  Resistance:         R(H=%.2f) = %.12e Ω\n",
                           V4_volts, V3_volts, V4_volts - V3_volts, I_amps,
                           dat.H_current, R_val);
        CHKERRQ(ierr);

        if (dat.Ro > 0) {
          double EMR_percent = 100.0 * (R_val - dat.Ro) / dat.Ro;
          ierr = PetscPrintf(PETSC_COMM_WORLD,
                             "  Reference resistance:   Ro = %.12e Ω\n"
                             "  EMR:                    = %.6f %%\n"
                             "  Resistance ratio:       R(H)/R(0) = %.6f\n"
                             "─────────────────────────────────────────────────"
                             "─────────────────\n"
                             "\n",
                             dat.Ro, EMR_percent, R_val / dat.Ro);
          CHKERRQ(ierr);
        }
      }
      // ============================================================================

      // Compute EMR Ratio (needs R0 from previous run or calc on fly)
      // Usually we just output R(H) and process later.
      // BUT simple EMR% = 100 * (R(H) - R(0)) / R(0).
      // We don't have R(0) stored here easily across calls unless we pass it.
      // Let's just write R(H) to file.

      /*
      // FILE WRITING MOVED TO main.cpp TO ENSURE R0 CONSISTENCY
      FILE *fp = fopen("EMRdata.out", "a");
      if (fp) {
        // Header if new file? (Ideally handle in main, but safe to append)
        // Format: H  R(H)  R2  ...
        fprintf(fp, "%.6e  %.8e  %.6e  %.6e  %.6e\n", dat.H_current, Resistance,
                dat.R2, dat.sigma1, dat.sigma2);
        fclose(fp);
      } else {
        PetscPrintf(PETSC_COMM_WORLD, "Error opening EMRdata.out for append\n");
      }
      */

      PetscPrintf(PETSC_COMM_WORLD,
                  "  >> EMR Calculation: H=%.4g T, V4_avg=%.6g V, R=%.8g Ohm\n",
                  dat.H_current, v4_avg, Resistance);
    }

    ierr = VecRestoreArray(tx, &sol_array);
    CHKERRQ(ierr);
  }

  // Clean up
  ierr = KSPDestroy(&ksp);
  CHKERRQ(ierr);
  ierr = VecDestroy(&x);
  CHKERRQ(ierr);
  ierr = VecDestroy(&tx);
  CHKERRQ(ierr);
  ierr = MatDestroy(&gmat.A);
  CHKERRQ(ierr);
  ierr = VecDestroy(&gmat.rhs);
  CHKERRQ(ierr);
  return 0;
}

// =============================================================================
// CONDITION NUMBER COMPUTATION
// =============================================================================

/**
 * @brief Computes or estimates the condition number of the matrix A.
 *
 * Uses SLEPc to find the largest and smallest eigenvalues to compute κ(A).
 * Also provides fallback estimates using matrix norms and diagonal analysis.
 *
 * @param A The system matrix.
 * @param dat Global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode compute_condition_number(Mat A, data &dat) {
  PetscErrorCode ierr;
  EPS eps;
  PetscReal lambda_max, lambda_min, cond_num;
  PetscInt nconv;
  PetscScalar eigval;

  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "\n"
      "==============================================================\n"
      "CONDITION NUMBER ANALYSIS\n"
      "==============================================================\n");
  CHKERRQ(ierr);

  PetscInt m, n;
  ierr = MatGetSize(A, &m, &n);
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Matrix size: %d × %d\n", m, n);
  CHKERRQ(ierr);

  // -------------------------------------------------------------------------
  // Method 1: Eigenvalue-based Condition Number (κ = λ_max / λ_min)
  // -------------------------------------------------------------------------
  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "\nMethod 1: Eigenvalue-based condition number\n"
      "--------------------------------------------------------------\n");
  CHKERRQ(ierr);

  // Create Eigenvalue Problem Solver (EPS) context
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
  CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, A, NULL);
  CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_HEP);
  CHKERRQ(ierr); // Hermitian Eigenproblem

  // --- Step 1: Find Largest Eigenvalue ---
  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Computing largest eigenvalue...\n");
  CHKERRQ(ierr);

  ierr = EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);
  CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps, 1, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps, 1e-6, 300);
  CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);
  CHKERRQ(ierr);
  ierr = EPSSolve(eps);
  CHKERRQ(ierr);

  ierr = EPSGetConverged(eps, &nconv);
  CHKERRQ(ierr);
  if (nconv > 0) {
    ierr = EPSGetEigenvalue(eps, 0, &eigval, NULL);
    CHKERRQ(ierr);
    lambda_max = PetscRealPart(eigval);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "    λ_max = %.6e\n", lambda_max);
    CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "     Warning: Largest eigenvalue did not converge\n");
    CHKERRQ(ierr);
    lambda_max = -1.0;
  }

  // --- Step 2: Find Smallest Eigenvalue ---
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "  Computing smallest eigenvalue (basic method)...\n");
  CHKERRQ(ierr);

  // Use SMALLEST_REAL. Note: Shift-and-invert is better but requires LU.
  ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
  CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps, 10, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps, 1e-4, 1000);
  CHKERRQ(ierr);

  ierr = EPSSolve(eps);
  CHKERRQ(ierr);
  ierr = EPSGetConverged(eps, &nconv);
  CHKERRQ(ierr);

  if (nconv > 0) {
    // Filter out zero eigenvalues (rigid body modes) to find smallest positive
    lambda_min = 1e100;
    int found_positive = 0;

    for (int i = 0; i < nconv; ++i) {
      ierr = EPSGetEigenvalue(eps, i, &eigval, NULL);
      CHKERRQ(ierr);
      PetscReal mag = PetscRealPart(eigval);

      if (mag > 1e-10 && mag < lambda_min) {
        lambda_min = mag;
        found_positive = 1;
      }
    }

    if (found_positive) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "    λ_min = %.6e (smallest positive)\n", lambda_min);
      CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "     No positive eigenvalues found (converged: %d)\n",
                         nconv);
      CHKERRQ(ierr);
      lambda_min = -1.0;
    }
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "    ⚠ Smallest eigenvalue did not converge\n");
    CHKERRQ(ierr);
    lambda_min = -1.0;
  }

  ierr = EPSDestroy(&eps);
  CHKERRQ(ierr);

  // --- Compute Condition Number ---
  if (lambda_max > 0 && lambda_min > 0) {
    cond_num = lambda_max / lambda_min;
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "\n  ★ Condition number (2-norm):\n"
                       "    κ(A) = λ_max / λ_min = %.6e\n\n",
                       cond_num);
    CHKERRQ(ierr);

    // Interpret result
    if (cond_num < 1e3) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "    ✓✓✓ EXCELLENT: Well-conditioned\n"
                         "    → Expected accuracy: ~14 digits\n");
      CHKERRQ(ierr);
    } else if (cond_num < 1e6) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "    ✓✓ GOOD: Moderately conditioned\n"
                         "    → Expected accuracy: ~10-13 digits\n");
      CHKERRQ(ierr);
    } else if (cond_num < 1e9) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "    ✓ FAIR: Some conditioning issues\n"
                         "    → Expected accuracy: ~7-9 digits\n");
      CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "    ⚠ POOR: Ill-conditioned\n"
                         "    → Expected accuracy: ~4-6 digits\n");
      CHKERRQ(ierr);
    }

    PetscReal precision_loss = PetscLog10Real(cond_num);
    ierr =
        PetscPrintf(PETSC_COMM_WORLD,
                    "    → Loss of decimal digits: ~%.1f\n"
                    "    → Target residual: < %.1e\n",
                    precision_loss, PetscPowReal(10.0, -16.0 + precision_loss));
    CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "\n  ⚠ Cannot compute exact condition number\n");
    CHKERRQ(ierr);
  }

  // -------------------------------------------------------------------------
  // Method 2: Matrix Norm Estimates
  // -------------------------------------------------------------------------
  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "\nMethod 2: Matrix norm-based estimates\n"
      "--------------------------------------------------------------\n");
  CHKERRQ(ierr);

  PetscReal norm_1, norm_inf, norm_fro;

  ierr = MatNorm(A, NORM_1, &norm_1);
  CHKERRQ(ierr);
  ierr = MatNorm(A, NORM_INFINITY, &norm_inf);
  CHKERRQ(ierr);
  ierr = MatNorm(A, NORM_FROBENIUS, &norm_fro);
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "  Matrix norms:\n"
                     "    ||A||_1   = %.6e (max column sum)\n"
                     "    ||A||_∞   = %.6e (max row sum)\n"
                     "    ||A||_F   = %.6e (Frobenius norm)\n",
                     norm_1, norm_inf, norm_fro);
  CHKERRQ(ierr);

  // Rough estimate: κ(A) ~ O(N²) for 2D Laplacian
  PetscInt N_estimate = (PetscInt)PetscSqrtReal((PetscReal)m);
  PetscReal kappa_estimate = PetscPowReal((PetscReal)N_estimate, 2.0);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n  Rough estimate (for FEM on ~%d×%d mesh):\n"
                     "    κ(A) ≈ O(%d²) ≈ %.1e\n",
                     N_estimate, N_estimate, N_estimate, kappa_estimate);
  CHKERRQ(ierr);

  // Matrix sparsity stats
  MatInfo info;
  ierr = MatGetInfo(A, MAT_GLOBAL_SUM, &info);
  CHKERRQ(ierr);

  double avg_nnz = info.nz_used / (double)m;
  double fill = info.nz_used / (double)(m * n);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n  Sparsity:\n"
                     "    Nonzeros:     %.0f\n"
                     "    Nonzeros/row: %.1f (avg)\n"
                     "    Fill:         %.2f%%\n",
                     info.nz_used, avg_nnz, fill * 100.0);
  CHKERRQ(ierr);

  // -------------------------------------------------------------------------
  // Method 3: Diagonal Scaling Analysis
  // -------------------------------------------------------------------------
  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "\nMethod 3: Diagonal analysis\n"
      "--------------------------------------------------------------\n");
  CHKERRQ(ierr);

  Vec diag;
  ierr = MatCreateVecs(A, NULL, &diag);
  CHKERRQ(ierr);
  ierr = MatGetDiagonal(A, diag);
  CHKERRQ(ierr);

  PetscReal diag_min, diag_max;
  ierr = VecMin(diag, NULL, &diag_min);
  CHKERRQ(ierr);
  ierr = VecMax(diag, NULL, &diag_max);
  CHKERRQ(ierr);

  PetscReal diag_ratio = diag_max / PetscAbsReal(diag_min);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "  Diagonal entries:\n"
                     "    min(diag) = %.6e\n"
                     "    max(diag) = %.6e\n"
                     "    ratio     = %.6e\n",
                     diag_min, diag_max, diag_ratio);
  CHKERRQ(ierr);

  if (diag_ratio < 1e6) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "    ✓ Diagonal is well-scaled\n");
    CHKERRQ(ierr);
  } else {
    ierr =
        PetscPrintf(PETSC_COMM_WORLD,
                    "    ⚠ Large diagonal ratio (consider diagonal scaling)\n");
    CHKERRQ(ierr);
  }

  // Check for zero or near-zero diagonals
  const PetscScalar *diag_array;
  PetscInt n_local;
  ierr = VecGetLocalSize(diag, &n_local);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(diag, &diag_array);
  CHKERRQ(ierr);

  int zero_diag = 0, small_diag = 0;
  for (int i = 0; i < n_local; ++i) {
    PetscReal abs_val = PetscAbsScalar(diag_array[i]);
    if (abs_val < 1e-14)
      zero_diag++;
    else if (abs_val < 1e-8)
      small_diag++;
  }

  ierr = VecRestoreArrayRead(diag, &diag_array);
  CHKERRQ(ierr);
  ierr = VecDestroy(&diag);
  CHKERRQ(ierr);

  int zero_global, small_global;
  ierr = MPI_Allreduce(&zero_diag, &zero_global, 1, MPI_INT, MPI_SUM,
                       PETSC_COMM_WORLD);
  CHKERRQ(ierr);
  ierr = MPI_Allreduce(&small_diag, &small_global, 1, MPI_INT, MPI_SUM,
                       PETSC_COMM_WORLD);
  CHKERRQ(ierr);

  if (zero_global > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "     WARNING: %d near-zero diagonal entries!\n",
                       zero_global);
    CHKERRQ(ierr);
  } else if (small_global > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "    Note: %d small diagonal entries (< 1e-8)\n",
                       small_global);
    CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "     All diagonals well-bounded\n");
    CHKERRQ(ierr);
  }

  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "==============================================================\n\n");
  CHKERRQ(ierr);

  return 0;
}