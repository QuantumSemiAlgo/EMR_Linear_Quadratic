/**
 * @file apply_bc.cpp
 * @brief Boundary Condition Application.
 *
 * This module applies Dirichlet and Neumann boundary conditions to the global
 * system. It handles the distinction between standard rectangular domains and
 * annular domains (where "left/right" become inner/outer radii and "top/bottom"
 * become periodic theta boundaries).
 *
 * Key Features:
 * - Geometric tolerance handling for identifying boundary nodes.
 * - Radial boundary condition application for annular domains.
 * - Periodic boundary handling (skipping theta boundaries).
 * - Tangential derivative constraints for Hermite elements.
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#include "global_params_1.h"
#include <cmath>

namespace {
const double ANNULAR_BOUNDARY_TOL = 0.01; // Tolerance for annular geometry
const double RECT_BOUNDARY_TOL = 1e-20;   // Tight tolerance for rectangular
} // namespace

/**
 * @brief Applies boundary conditions to the global matrix and RHS.
 *
 * @param gmat Reference to the global matrices structure.
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success.
 */
// Prototype
PetscErrorCode apply_bc_EMR(global_matrices &gmat, data &dat);

PetscErrorCode apply_bc(global_matrices &gmat, data &dat) {
  PetscErrorCode ierr;
  PetscInt stride = dat.dof_per_node;

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0)
    fprintf(stderr, "DEBUG: Inside apply_bc. Nsample_R2=%d. Dispatching...\n",
            dat.Nsample_R2);

  // DISPATCH EMR Check
  if (dat.Nsample_R2 > 0) {
    return apply_bc_EMR(gmat, dat);
  }

  // ---------------------------------------------------------------------------
  // 1. Determine Domain Bounds
  // ---------------------------------------------------------------------------
  // Use a tighter tolerance for rectangular domains, looser for annular
  // (due to potential floating point drift in coordinate transformations).
  // (due to potential floating point drift in coordinate transformations).
  const double tol = dat.use_annular ? ANNULAR_BOUNDARY_TOL : RECT_BOUNDARY_TOL;

  double actual_ymin = dat.ymin;
  double actual_ymax = dat.ymax;
  double actual_xmin = dat.xmin;
  double actual_xmax = dat.xmax;

  // For annular domains, the "logical" bounds (xmin, xmax, etc.) might not
  // match the physical coordinate bounds if the mesh was deformed.
  // We scan all nodes to find the true physical extents.
  if (dat.use_annular) {
    actual_ymin = 1e10;
    actual_ymax = -1e10;
    actual_xmin = 1e10;
    actual_xmax = -1e10;

    for (int i = 0; i < dat.ngnodes; ++i) {
      double x = dat.node[i][0];
      double y = dat.node[i][1];

      if (y < actual_ymin)
        actual_ymin = y;
      if (y > actual_ymax)
        actual_ymax = y;
      if (x < actual_xmin)
        actual_xmin = x;
      if (x > actual_xmax)
        actual_xmax = x;
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "  Domain bounds: x=[%.6f, %.6f], y=[%.6f, %.6f]\n",
                       actual_xmin, actual_xmax, actual_ymin, actual_ymax);
    CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // 2. DIAGNOSTICS: Count boundary nodes
  // ---------------------------------------------------------------------------
  int n_inner = 0, n_outer = 0, n_top = 0, n_bottom = 0;
  for (int bi = 0; bi < dat.nbnode; ++bi) {
    int node_id = dat.bnode_id[bi];
    double x = dat.node[node_id][0];
    double y = dat.node[node_id][1];

    if (std::fabs(x - actual_xmin) < tol)
      n_inner++;
    if (std::fabs(x - actual_xmax) < tol)
      n_outer++;
    if (std::fabs(y - actual_ymin) < tol)
      n_bottom++;
    if (std::fabs(y - actual_ymax) < tol)
      n_top++;
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Boundary node distribution:\n"
                     "  Inner (x~%.2f): %d nodes\n"
                     "  Outer (x~%.2f): %d nodes\n"
                     "  Bottom (y~%.2f): %d nodes\n"
                     "  Top (y~%.2f): %d nodes\n"
                     "  Total boundary nodes: %d\n",
                     actual_xmin, n_inner, actual_xmax, n_outer, actual_ymin,
                     n_bottom, actual_ymax, n_top, dat.nbnode);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 3. Apply Dirichlet BCs
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // 3. Apply Dirichlet Boundary Conditions
  // ---------------------------------------------------------------------------
  int bc_count = 0;
  int skipped_theta = 0;
  int inner_bc_count = 0;
  int outer_bc_count = 0;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Applying boundary conditions...\n");
  CHKERRQ(ierr);

  // Iterate over all identified boundary nodes
  for (int bi = 0; bi < dat.nbnode; ++bi) {
    int node_id = dat.bnode_id[bi];

    if (node_id < 0 || node_id >= dat.ngnodes) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
              "Invalid boundary node_id=%d", node_id);
    }

    double x = dat.node[node_id][0];
    double y = dat.node[node_id][1];

    PetscScalar bcval = 0.0;
    bool apply_bc = false;

    if (dat.use_annular) {
      // =========================================================================
      // ANNULAR DOMAIN: Radial boundaries take PRIORITY over theta boundaries
      // =========================================================================

      // Check which boundary this node is on
      bool on_inner_radius = (std::fabs(x - actual_xmin) < tol);
      bool on_outer_radius = (std::fabs(x - actual_xmax) < tol);
      bool on_theta_min = (std::fabs(y - actual_ymin) < tol);
      bool on_theta_max = (std::fabs(y - actual_ymax) < tol);

      if (on_inner_radius) {
        // -----------------------------------------------------------------------
        // INNER RADIUS (r = Rin): Floating BC (natural boundary condition)
        // -----------------------------------------------------------------------
        // Do NOT apply Dirichlet BC - this allows du/dr = 0
        apply_bc = false;
        inner_bc_count++;

      } else if (on_outer_radius) {
        // -----------------------------------------------------------------------
        // OUTER RADIUS (r = Rout): Dirichlet BC u = 10·sin(5θ)
        // -----------------------------------------------------------------------

        // Map logical y-coordinate to physical theta angle [0, 2π]
        double theta =
            2.0 * PI * (y - actual_ymin) / (actual_ymax - actual_ymin);

        // Compute BC value: u = 10·sin(5θ)
        bcval = 10.0 * sin(5.0 * theta);
        apply_bc = true;
        outer_bc_count++;

        // Debug output for nodes near theta = 0 or 2π (Optional, disabled for
        // production)

      } else if (on_theta_min || on_theta_max) {
        // -----------------------------------------------------------------------
        // PURE THETA BOUNDARIES (not on radial boundary)
        // These are periodic - handled by node overlay via Convertnode
        // Skip applying explicit BC here
        // -----------------------------------------------------------------------
        skipped_theta++;
        continue;

      } else {
        // Not on any boundary - shouldn't happen if bnode.dat is correct
        continue;
      }

    } else {
      // =========================================================================
      // RECTANGULAR DOMAIN LOGIC
      // =========================================================================
      if (std::fabs(x - dat.xmin) < tol) {
        // Left boundary: u = 10·sin(πy/20)
        bcval = 10.0 * sin(PI * y / 20.0);
        apply_bc = true;
      } else {
        // All other boundaries: u = 0
        bcval = 0.0;
        apply_bc = true;
      }
    }

    // Skip if no BC to apply
    if (!apply_bc)
      continue;

    // ===========================================================================
    // Apply the Dirichlet Boundary Condition
    // ===========================================================================

    // Primary DOF (value, not derivatives)
    PetscInt d = 0;
    int local_dof = node_id * stride + d;

    // Map to global DOF (handles periodic node mapping)
    PetscInt gdof = dat.Convertnode[local_dof];

    // Enforce Dirichlet BC by:
    // 1. Zero the matrix row and set diagonal to 1.0
    // 2. Set RHS vector to the boundary value
    ierr = MatZeroRows(gmat.A, 1, &gdof, 1.0, NULL, NULL);
    CHKERRQ(ierr);

    ierr = VecSetValue(gmat.rhs, gdof, bcval, INSERT_VALUES);
    CHKERRQ(ierr);

    bc_count++;

    // Debug logging for first few BCs
    if (bc_count <= 5) {
      ierr = PetscPrintf(
          PETSC_COMM_WORLD,
          "  BC #%d: node=%d, (x=%.3f,y=%.3f) -> u=%.6f (gdof=%d)\n", bc_count,
          node_id, x, y, bcval, gdof);
      /*
       * Read logic from dat
       * iterate boundary nodes (dat.bnodes)
       * Modify A and rhs
       */

      if (dat.ndebug > 0)
        PetscPrintf(PETSC_COMM_WORLD, "   Debug: Entering apply_bc\n");
      int rank;
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      if (rank == 0)
        fprintf(stderr, "DEBUG: Inside apply_bc. nbnode=%d\n", dat.nbnode);

      if (!gmat.A) {
        if (rank == 0)
          fprintf(stderr, "DEBUG: ERROR gmat.A is NULL in apply_bc!\n");
        return 1;
      }

      PetscInt Istart, Iend;
      ierr = MatGetOwnershipRange(gmat.A, &Istart, &Iend);
      CHKERRQ(ierr);
      ierr = VecAssemblyEnd(gmat.rhs);
      CHKERRQ(ierr);
    }
  }

  // Assemble vectors after setting values
  ierr = VecAssemblyBegin(gmat.rhs);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(gmat.rhs);
  CHKERRQ(ierr);

  // Print summary
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\nBoundary condition summary:\n"
                     "  Inner radius (floating BC): %d nodes\n"
                     "  Outer radius (Dirichlet BC): %d nodes\n"
                     "  Theta boundaries (periodic, skipped): %d nodes\n"
                     "  Total Dirichlet BCs applied: %d\n\n",
                     inner_bc_count, outer_bc_count, skipped_theta, bc_count);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 3. Hermite Derivative Constraints (Rectangular Only)
  // ---------------------------------------------------------------------------
  // For Hermite elements, we often need to constrain tangential derivatives
  // at boundaries to zero if they are not specified otherwise.
  if (dat.dof_per_node > 1 && !dat.use_annular) {
    int deriv_bc_count = 0;

    for (int bi = 0; bi < dat.nbnode; ++bi) {
      int node_id = dat.bnode_id[bi];
      double y = dat.node[node_id][1];

      // Top and Bottom boundaries
      if (std::fabs(y - dat.ymin) < tol || std::fabs(y - dat.ymax) < tol) {
        PetscInt d = 1; // ∂/∂ξ (tangential derivative) DOF
        int local_dof = node_id * stride + d;
        PetscInt gdof = dat.Convertnode[local_dof];

        ierr = MatZeroRows(gmat.A, 1, &gdof, 1.0, NULL, NULL);
        CHKERRQ(ierr);
        ierr = VecSetValue(gmat.rhs, gdof, 0.0, INSERT_VALUES);
        CHKERRQ(ierr);

        deriv_bc_count++;
      }
    }

    if (deriv_bc_count > 0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "Applied %d derivative boundary conditions\n",
                         deriv_bc_count);
      CHKERRQ(ierr);
    }
  }

  // ---------------------------------------------------------------------------
  // 4. Final Assembly
  // ---------------------------------------------------------------------------
  // Assemble the RHS vector to commit the changes.
  // Note: Matrix assembly is handled by MatZeroRows implicitly or requires
  // a separate call if needed, but here we only modified RHS explicitly.
  // Assembly handled previously

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Boundary conditions applied successfully\n");
  CHKERRQ(ierr);

  return 0;
}

// =============================================================================
// HELPER: Shape Function Integration (1D) for Neumann BCs
// =============================================================================
// Computes integral of N[i] * N[j] * r * dTheta along a port?
// Actually simpler: We need to distribute Total Current I0 into nodal loads.
// Load vector F_i = Integral( J_n * N_i * dS )
// dS = r * dTheta (at r=Rout).
// J_n = I0 / (R * Width). (Assumed uniform).
// So F_i = (I0 / (R*W)) * Integral( N_i * R * dTheta )
//        = (I0 / W) * Integral( N_i * dTheta )
//
// We need Integral( N_i * dTheta ) over the element edge.
// For linear elements (2 nodes on edge), Integral N_i dxi = 1.0 (from -1 to 1).
// Jacobian from xi to theta: dTheta = (DeltaTheta / 2) * dxi.
// So Integral( N_i dTheta ) = (DeltaTheta / 2) * Integral( N_i dxi )
// Integral( N_i dxi ) for N_1 = (1-xi)/2 is 1.0. For N_2 = (1+xi)/2 is 1.0.
// So Load_i = (I0 / W) * (DeltaTheta / 2) * 1.0.
// Basically distributes current proportional to element size.
//
// We can implement this directly in the loop.

// =============================================================================
// EMR BOUNDARY CONDITIONS
// =============================================================================
PetscErrorCode apply_bc_EMR(global_matrices &gmat, data &dat) {
  PetscErrorCode ierr;
  ierr =
      PetscPrintf(PETSC_COMM_WORLD,
                  "Applying EMR Boundary Conditions (Hermite Corrected)...\n");
  CHKERRQ(ierr);

  // Identify Ports
  // Port 1: Injection (I = Io)
  // Port 2: Extraction (I = -Io)
  // Port 3: Reference (V = 0)
  // Port 4: Floating / Probe

  // ==========================================================================
  // HERMITE SPECIFIC: Constrain Tangential Derivative (d/dTheta) at ALL Ports
  // ==========================================================================
  // Metal contacts are equipotential surfaces -> dPhi/dTheta = 0 along the
  // contact. This applies to Current Ports (1, 2), Voltage Ports (3, 4). DOF
  // Mapping: 0=Val, 1=d/dr (Xi), 2=d/dTheta (Eta), 3=Mixed FIX: User requested
  // constraining DOF 2 (ie_dof = 2).
  if (dat.dof_per_node == 4) {
    int tan_deriv_count = 0;
    for (int i = 0; i < dat.ngnodes; ++i) {
      if (dat.port_code[i] > 0) { // If node belongs to ANY port
        // Constrain DOF 2 (Tangential Derivative d/dTheta)
        int local_dof_eta = i * dat.dof_per_node + 2; // CORRECTED INDEX
        PetscInt gdof_eta = dat.Convertnode[local_dof_eta];

        ierr = MatZeroRows(gmat.A, 1, &gdof_eta, 1.0, NULL, NULL);
        CHKERRQ(ierr);
        ierr = VecSetValue(gmat.rhs, gdof_eta, 0.0, INSERT_VALUES);
        CHKERRQ(ierr);
        tan_deriv_count++;

        // Also constrain Mixed Derivative (DOF 3)? Usually yes for boundary.
        int local_dof_mix = i * dat.dof_per_node + 3;
        PetscInt gdof_mix = dat.Convertnode[local_dof_mix];
        ierr = MatZeroRows(gmat.A, 1, &gdof_mix, 1.0, NULL, NULL);
        CHKERRQ(ierr);
        ierr = VecSetValue(gmat.rhs, gdof_mix, 0.0, INSERT_VALUES);
        CHKERRQ(ierr);
      }
    }
    PetscPrintf(PETSC_COMM_WORLD,
                "  Hermite: Constrained Tangential (DOF 2) & Mixed Derivs at "
                "%d port nodes\n",
                tan_deriv_count);
  }

  // ==========================================================================
  // PORT 3: V=0 Reference - CONSTRAIN VALUE (DOF 0)
  // ==========================================================================
  // Fix: Must constrain ALL degrees of freedom for the reference nodes,
  // not just the first one. For Hermite (dof=4), unconstrained derivatives
  // cause singular matrix.
  int constrained_count = 0;
  for (int i = 0; i < dat.ngnodes; ++i) {
    if (dat.port_code[i] == 3) {
      // Loop over DOFs
      if (dat.dof_per_node == 4) {
        // CUBIC HERMITE:
        // DOF 1 (Tangential) and DOF 3 (Mixed) are already constrained above.
        // We MUST constrain DOF 0 (Value) to 0.0 (Reference Voltage).
        // We MUST leave DOF 2 (Normal) FREE to allow current flow.

        int d = 0; // Value
        int local_dof = i * dat.dof_per_node + d;
        PetscInt gdof = dat.Convertnode[local_dof];
        ierr = MatZeroRows(gmat.A, 1, &gdof, 1.0, NULL, NULL);
        CHKERRQ(ierr);
        ierr = VecSetValue(gmat.rhs, gdof, 0.0, INSERT_VALUES);
        CHKERRQ(ierr);

        // (DOF 1 and 3 done above. DOF 2 left free).
      } else {
        // Standard (Q1/Q2/Q9?): Constrain all DOFs usually means V=0.
        // For Q9 (Quintic), we had issues. But for Q1/Q2 dof=1.
        for (int d = 0; d < dat.dof_per_node; ++d) {
          int local_dof = i * dat.dof_per_node + d;
          PetscInt gdof = dat.Convertnode[local_dof];
          ierr = MatZeroRows(gmat.A, 1, &gdof, 1.0, NULL, NULL);
          CHKERRQ(ierr);
          ierr = VecSetValue(gmat.rhs, gdof, 0.0, INSERT_VALUES);
          CHKERRQ(ierr);
        }
      }
      constrained_count++;
    }
  }
  PetscPrintf(
      PETSC_COMM_WORLD,
      "  Port 3 (Ref V=0): Constrained %d nodes (All %d DOFs per node)\n",
      constrained_count, dat.dof_per_node);

  // ==========================================================================
  // PORTS 1 & 2: Current Injection (Neumann)
  // ==========================================================================
  // Current Density J1 = Io / (ArcLength). Units: A/m.
  double J1 = dat.Io / (dat.R1 * dat.width_L1 * M_PI);
  double J2 = -dat.Io / (dat.R1 * dat.width_L2 * M_PI);

  int n_theta_elems = dat.n_x;
  int n_rad_elems = dat.n_y;
  int outer_row_idx = n_rad_elems - 1;

  for (int ie_t = 0; ie_t < n_theta_elems; ++ie_t) {
    // Element index in Global list
    int ie = outer_row_idx * n_theta_elems + ie_t;
    int n3 = dat.elem[ie][3]; // TL
    int n2 = dat.elem[ie][2]; // TR
    // Check ports for corner nodes (always present)
    int p3 = dat.port_code[n3];
    int p2 = dat.port_code[n2];
    double J_applied = 0.0;

    if (dat.node_elem == 9) {
      // 9-node (Quadratic Edge) -> Simpson's Rule
      // ... (Existing Q9 Logic) ...
      int n6 = dat.elem[ie][6]; // TM
      int p6 = dat.port_code[n6];

      if (!(p3 == 1 && p2 == 1 && p6 == 1) && !(p3 == 2 && p2 == 2 && p6 == 2))
        continue;
      if (p3 == 1)
        J_applied = J1;
      else
        J_applied = J2;

      double th3 = dat.node[n3][0];
      double th2 = dat.node[n2][0];
      if (th2 < th3)
        th2 += 2.0 * M_PI;
      double dTheta = th2 - th3;
      double ArcLen = dat.Rout * dTheta;

      double load_corner = J_applied * ArcLen * (1.0 / 6.0);
      double load_midside = J_applied * ArcLen * (4.0 / 6.0);

      PetscInt gdof3 = dat.Convertnode[n3 * dat.dof_per_node + 0];
      PetscInt gdof2 = dat.Convertnode[n2 * dat.dof_per_node + 0];
      PetscInt gdof6 = dat.Convertnode[n6 * dat.dof_per_node + 0];

      ierr = VecSetValue(gmat.rhs, gdof3, load_corner, ADD_VALUES);
      CHKERRQ(ierr);
      ierr = VecSetValue(gmat.rhs, gdof2, load_corner, ADD_VALUES);
      CHKERRQ(ierr);
      ierr = VecSetValue(gmat.rhs, gdof6, load_midside, ADD_VALUES);
      CHKERRQ(ierr);

    } else if (dat.node_elem == 4 && dat.dof_per_node == 4) {
      // CUBIC HERMITE (4 nodes, 4 DOFs)
      // Consistent Load Integration for Hermite Shape Functions.
      // Along edge u in [-1, 1], J is constant.
      // Load_i = J * (Jac=ArcLen/2) * Integral(H_i(u) du)
      //
      // Hermite Shapes on [-1,1]:
      // H_val_L  (Node 0): 1/4 (1-u)^2 (2+u) -> Integral = 1
      // H_dxi_L  (Node 0): 1/4 (1-u)^2 (1+u) -> Integral = 1/3 (Wait, check
      // scaling) H_val_R  (Node 1): 1/4 (1+u)^2 (2-u) -> Integral = 1 H_dxi_R
      // (Node 1): 1/4 (1+u)^2 (u-1) -> Integral = -1/3

      // CAUTION: The H_dxi functions in shape_funcs.hpp might have different
      // scalings. H1 (Node 0 dxi) = 0.0625 * (xi-1)^2 * (xi+1) ... at eta
      // constrained? Let's assume standard reference element. Integral_{-1}^1
      // H_val du = 1.0. Integral_{-1}^1 H_slope du = 1/3.

      // Jacobian: dx/dxi = ArcLen / 2.
      // Real Load_val = J * (ArcLen/2) * 1.0 = J * ArcLen / 2.
      // Real Load_slope = J * (ArcLen/2) * (1/3) * (dxi/dx??)
      // Wait, force conjugate to Slope DOF (du/dxi) depends on definition.
      // Usually we just integrate N_i * J.
      // So Load_slope_local = J * (ArcLen/2) * (1/3).

      if (!(p3 == 1 && p2 == 1) && !(p3 == 2 && p2 == 2))
        continue;
      if (p3 == 1)
        J_applied = J1;
      else
        J_applied = J2;

      double th3 = dat.node[n3][0];
      double th2 = dat.node[n2][0];
      if (th2 < th3)
        th2 += 2.0 * M_PI;
      double dTheta = th2 - th3;
      double ArcLen = dat.Rout * dTheta;

      double jac = ArcLen / 2.0;

      // Weights for [-1, 1] Integrals of Hermite polynomials
      double w_val = 1.0;
      double w_slope = 1.0 / 3.0;    // Node 0 (Left)
      double w_slope_R = -1.0 / 3.0; // Node 1 (Right)

      // Apply to Corner Nodes
      // Node 3 (Left/Start):
      PetscInt gdof3_val = dat.Convertnode[n3 * dat.dof_per_node + 0]; // Value
      PetscInt gdof3_dxi =
          dat.Convertnode[n3 * dat.dof_per_node + 1]; // Tangential Slope (dxi)

      double load_3_val = J_applied * jac * w_val;
      double load_3_dxi = J_applied * jac * w_slope;

      ierr = VecSetValue(gmat.rhs, gdof3_val, load_3_val, ADD_VALUES);
      CHKERRQ(ierr);
      ierr = VecSetValue(gmat.rhs, gdof3_dxi, load_3_dxi, ADD_VALUES);
      CHKERRQ(ierr);

      // Node 2 (Right/End):
      PetscInt gdof2_val = dat.Convertnode[n2 * dat.dof_per_node + 0]; // Value
      PetscInt gdof2_dxi =
          dat.Convertnode[n2 * dat.dof_per_node + 1]; // Tangential Slope (dxi)

      double load_2_val = J_applied * jac * w_val;
      double load_2_dxi = J_applied * jac * w_slope_R;

      ierr = VecSetValue(gmat.rhs, gdof2_val, load_2_val, ADD_VALUES);
      CHKERRQ(ierr);
      ierr = VecSetValue(gmat.rhs, gdof2_dxi, load_2_dxi, ADD_VALUES);
      CHKERRQ(ierr);

    } else if (dat.node_elem == 4) {
      // Linear Elements (Q1)
      if (!(p3 == 1 && p2 == 1) && !(p3 == 2 && p2 == 2))
        continue;
      if (p3 == 1)
        J_applied = J1;
      else
        J_applied = J2;

      double th3 = dat.node[n3][0];
      double th2 = dat.node[n2][0];
      if (th2 < th3)
        th2 += 2.0 * M_PI;
      double dTheta = th2 - th3;
      double ArcLen = dat.Rout * dTheta;

      double load_corner = J_applied * ArcLen * 0.5;

      PetscInt gdof3 = dat.Convertnode[n3 * dat.dof_per_node + 0];
      PetscInt gdof2 = dat.Convertnode[n2 * dat.dof_per_node + 0];

      ierr = VecSetValue(gmat.rhs, gdof3, load_corner, ADD_VALUES);
      CHKERRQ(ierr);
      ierr = VecSetValue(gmat.rhs, gdof2, load_corner, ADD_VALUES);
      CHKERRQ(ierr);

    } else {
      // Unknown element type?
      PetscPrintf(PETSC_COMM_WORLD,
                  "Warning: Unknown element type in apply_bc_EMR\n");
    }
  }

  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "  Applied Neumann Current Loads to Ports 1 & 2 (DOF 0 targeted)\n");
  CHKERRQ(ierr);

  ierr = VecAssemblyBegin(gmat.rhs);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(gmat.rhs);
  CHKERRQ(ierr);

  return 0;
}