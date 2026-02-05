/**
 * @file make_global_EMR_FIXED.cpp
 * @brief CORRECTED VERSION - Global Stiffness Matrix Assembly for EMR
 *
 * CRITICAL FIX: Added proper Fermi function implementation
 *
 * BUG FOUND: The original code had only a declaration:
 *     double fermi(double r, double R2, double delta);
 *
 * This meant either:
 * 1. No material transition was happening (f = constant)
 * 2. Linker was using a stub that returned wrong values
 * 3. Material properties were uniform across the disk
 *
 * RESULT: No Hall parameter variation → No EMR effect!
 *
 * @author Arya Retheeshan (with critical bugfix)
 * @date 2025-02-02
 */

#include "global_params_1.h"
#include "prototypes_1.h"
#include "shape_funcs_1.hpp"
#include "utilrect_1.h"
#include <cmath>
#include <petsc.h>
#include <petscksp.h>
#include <vector>

// =============================================================================
// FERMI FUNCTION - PROPER IMPLEMENTATION
// =============================================================================
/**
 * @brief Smooth material transition function (Fermi-Dirac form).
 *
 * Creates a smooth "S-curve" transition between metal (inner) and semiconductor
 * (outer) regions. This avoids sharp discontinuities that cause numerical
 * issues.
 *
 * Mathematical Form:
 *   f(r) = 1 / (1 + exp((r - R2) / delta))
 *
 * Behavior:
 *   - r << R2 (inside metal):  f → 1.0
 *   - r = R2  (interface):     f = 0.5
 *   - r >> R2 (outside, semi): f → 0.0
 *
 * @param r Current radial position [meters]
 * @param R2 Metal-semiconductor interface radius [meters]
 * @param delta Transition width parameter [dimensionless, typically 0.1-10]
 * @return Fermi factor (0 = pure semiconductor, 1 = pure metal)
 */
double fermi(double r, double R2, double delta) {
  // Normalized distance from interface
  double arg = (r - R2) / delta;

  // Prevent overflow in exp() for large |arg|
  if (arg > 50.0) {
    return 0.0; // Far outside R2 → pure semiconductor
  }
  if (arg < -50.0) {
    return 1.0; // Far inside R2 → pure metal
  }

  // Standard Fermi function
  return 1.0 / (1.0 + std::exp(arg));
}

// =============================================================================
// EMR ELEMENT STIFFNESS MATRIX
// =============================================================================
/**
 * @brief Assembles local element stiffness matrix with EMR physics.
 *
 * Implements the weak form of the EMR equation in polar coordinates:
 *
 *   ∫∫ [∇φ]ᵀ · [σ_eff] · [∇φ] · r · dr · dθ
 *
 * Where the effective conductivity tensor includes Hall effect:
 *
 *   [σ_eff] = (σ / (1 + β²)) · [ 1    -β  ]
 *                                [ β     1  ]
 *
 * With spatially-varying material properties:
 *   σ(r) = fermi(r) · σ_metal + (1 - fermi(r)) · σ_semi
 *   β(r) = μ(r) · |H|
 *
 * @param dat Global data structure (contains mesh, physics params)
 * @param ie Element index
 * @param Ae Output local stiffness matrix [nloc × nloc], flattened row-major
 */
void element_stiffness_EMR(const data &dat, int ie, PetscScalar *Ae) {

  // -------------------------------------------------------------------------
  // 1. Initialize and Extract Element Data
  // -------------------------------------------------------------------------
  int nloc = dat.nele_mat;
  for (int i = 0; i < nloc * nloc; ++i) {
    Ae[i] = 0.0;
  }

  // Get element node connectivity
  std::vector<int> nodes(dat.node_elem);
  for (int k = 0; k < dat.node_elem; ++k) {
    nodes[k] = dat.elem[ie][k];
  }

  // Extract nodal coordinates (x = theta, y = r)
  std::vector<double> theta_n(dat.node_elem);
  std::vector<double> r_n(dat.node_elem);
  for (int k = 0; k < dat.node_elem; ++k) {
    int nid = nodes[k];
    theta_n[k] = dat.node[nid][0]; // x-coordinate is theta
    r_n[k] = dat.node[nid][1];     // y-coordinate is r
  }

  // -------------------------------------------------------------------------
  // 2. Setup Quadrature & Shape Functions
  // -------------------------------------------------------------------------
  // Use 3×3 Gauss for Quintic Hermite/Q9 integration
  int n_q_1d = 3;
  double xi_q[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
  double w_q[3] = {0.5555555555555556, 0.8888888888888889, 0.5555555555555556};

  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------
  // HERMITE SHAPE FUNCTION POINTERS
  // -------------------------------------------------------------------------
  typedef double (*ShapeFunc)(double, double);

  // Quintic Hermite (36 DOFs) - Q9
  static const ShapeFunc ptr_N_Q9[] = {
      N0,  N1,  N2,  N3,  N4,  N5,  N6,  N7,  N8,  N9,  N10, N11,
      N12, N13, N14, N15, N16, N17, N18, N19, N20, N21, N22, N23,
      N24, N25, N26, N27, N28, N29, N30, N31, N32, N33, N34, N35};

  static const ShapeFunc ptr_dNdxi_Q9[] = {
      dN0_dxi,  dN1_dxi,  dN2_dxi,  dN3_dxi,  dN4_dxi,  dN5_dxi,
      dN6_dxi,  dN7_dxi,  dN8_dxi,  dN9_dxi,  dN10_dxi, dN11_dxi,
      dN12_dxi, dN13_dxi, dN14_dxi, dN15_dxi, dN16_dxi, dN17_dxi,
      dN18_dxi, dN19_dxi, dN20_dxi, dN21_dxi, dN22_dxi, dN23_dxi,
      dN24_dxi, dN25_dxi, dN26_dxi, dN27_dxi, dN28_dxi, dN29_dxi,
      dN30_dxi, dN31_dxi, dN32_dxi, dN33_dxi, dN34_dxi, dN35_dxi};

  static const ShapeFunc ptr_dNdeta_Q9[] = {
      dN0_deta,  dN1_deta,  dN2_deta,  dN3_deta,  dN4_deta,  dN5_deta,
      dN6_deta,  dN7_deta,  dN8_deta,  dN9_deta,  dN10_deta, dN11_deta,
      dN12_deta, dN13_deta, dN14_deta, dN15_deta, dN16_deta, dN17_deta,
      dN18_deta, dN19_deta, dN20_deta, dN21_deta, dN22_deta, dN23_deta,
      dN24_deta, dN25_deta, dN26_deta, dN27_deta, dN28_deta, dN29_deta,
      dN30_deta, dN31_deta, dN32_deta, dN33_deta, dN34_deta, dN35_deta};

  // Bicubic Hermite (16 DOFs) - Q4
  static const ShapeFunc ptr_N_Q4[] = {H0, H1, H2,  H3,  H4,  H5,  H6,  H7,
                                       H8, H9, H10, H11, H12, H13, H14, H15};

  static const ShapeFunc ptr_dNdxi_Q4[] = {
      dH0_dxi,  dH1_dxi,  dH2_dxi,  dH3_dxi, dH4_dxi,  dH5_dxi,
      dH6_dxi,  dH7_dxi,  dH8_dxi,  dH9_dxi, dH10_dxi, dH11_dxi,
      dH12_dxi, dH13_dxi, dH14_dxi, dH15_dxi};

  static const ShapeFunc ptr_dNdeta_Q4[] = {
      dH0_deta,  dH1_deta,  dH2_deta,  dH3_deta, dH4_deta,  dH5_deta,
      dH6_deta,  dH7_deta,  dH8_deta,  dH9_deta, dH10_deta, dH11_deta,
      dH12_deta, dH13_deta, dH14_deta, dH15_deta};

  // -------------------------------------------------------------------------
  // 3. Quadrature Loop
  // -------------------------------------------------------------------------
  for (int iq = 0; iq < n_q_1d; ++iq) {
    for (int jq = 0; jq < n_q_1d; ++jq) {
      double xi = xi_q[iq];
      double eta = xi_q[jq];
      double w = w_q[iq] * w_q[jq];

      // -----------------------------------------------------------------
      // 3a. Geometry Shape Functions
      // -----------------------------------------------------------------
      // Used ONLY for Jacobian and coordinate transformation
      std::vector<double> N_geom(dat.node_elem);
      std::vector<double> dNdxi_geom(dat.node_elem);
      std::vector<double> dNdeta_geom(dat.node_elem);

      if (dat.node_elem == 9) {
        shapeQ2(N_geom.data(), xi, eta);
        deriv1Q2(dNdxi_geom.data(), dNdeta_geom.data(), xi, eta);
      } else if (dat.node_elem == 4) {
        // Linear Elements (Wait! If using CUBIC HERMITE, geometry can still be
        // Q1 Linear) If dat.node_elem == 4, we use Q1 shape for geometry
        // mapping.
        shapeQ1(N_geom.data(), xi, eta);
        deriv1Q1(dNdxi_geom.data(), dNdeta_geom.data(), xi, eta);
      } else {
        // Fallback default (Q2)
        shapeQ2(N_geom.data(), xi, eta);
        deriv1Q2(dNdxi_geom.data(), dNdeta_geom.data(), xi, eta);
      }

      // -----------------------------------------------------------------
      // 3b. Jacobian Calculation (Parameter → Physical)
      // -----------------------------------------------------------------
      double J11 = 0.0, J12 = 0.0;
      double J21 = 0.0, J22 = 0.0;
      double theta_phys = 0.0;
      double r_phys = 0.0;

      for (int k = 0; k < dat.node_elem; ++k) {
        J11 += dNdxi_geom[k] * theta_n[k];
        J12 += dNdxi_geom[k] * r_n[k];
        J21 += dNdeta_geom[k] * theta_n[k];
        J22 += dNdeta_geom[k] * r_n[k];

        theta_phys += N_geom[k] * theta_n[k];
        r_phys += N_geom[k] * r_n[k];
      }

      double detJ = J11 * J22 - J12 * J21;
      if (detJ <= 0.0)
        continue;
      double invDet = 1.0 / detJ;

      // -----------------------------------------------------------------
      // 3c. Field Shape Functions (HERMITE vs LAGRANGE)
      // -----------------------------------------------------------------
      int num_field_dofs = dat.node_elem * dat.dof_per_node; //
      std::vector<double> N_field(num_field_dofs);
      std::vector<double> dNdxi_field(num_field_dofs);
      std::vector<double> dNdeta_field(num_field_dofs);

      // Populate Hermite shapes using function pointers
      if (dat.dof_per_node == 4 && num_field_dofs == 36) {
        // QUINTIC HERMITE
        for (int k = 0; k < 36; ++k) {
          N_field[k] = ptr_N_Q9[k](xi, eta);
          dNdxi_field[k] = ptr_dNdxi_Q9[k](xi, eta);
          dNdeta_field[k] = ptr_dNdeta_Q9[k](xi, eta);
        }
      } else if (dat.dof_per_node == 4 && num_field_dofs == 16) {
        // BICUBIC HERMITE (NEW)
        for (int k = 0; k < 16; ++k) {
          N_field[k] = ptr_N_Q4[k](xi, eta);
          dNdxi_field[k] = ptr_dNdxi_Q4[k](xi, eta);
          dNdeta_field[k] = ptr_dNdeta_Q4[k](xi, eta);
        }
      } else {
        // Fallback to Geometry shapes if not Hermite (e.g. Lagrangian Q2 or Q1)
        // This allows code to technically run if user sets ndof=1
        for (int k = 0; k < dat.node_elem; ++k) {
          N_field[k] = N_geom[k];
          dNdxi_field[k] = dNdxi_geom[k];
          dNdeta_field[k] = dNdeta_geom[k];
        }
        num_field_dofs = dat.node_elem; // Reduces to 9 or 4
      }

      // -----------------------------------------------------------------
      // 3d. Transform Derivatives (Chain Rule)
      // -----------------------------------------------------------------
      std::vector<double> dNdtheta(num_field_dofs);
      std::vector<double> dNdr(num_field_dofs);

      for (int k = 0; k < num_field_dofs; ++k) {
        dNdtheta[k] = invDet * (J22 * dNdxi_field[k] - J12 * dNdeta_field[k]);
        dNdr[k] = invDet * (-J21 * dNdxi_field[k] + J11 * dNdeta_field[k]);
      }

      // -----------------------------------------------------------------
      // 3e. Material Properties via Fermi Function
      // -----------------------------------------------------------------
      double f = fermi(r_phys, dat.R2, dat.delta);
      double sigma = f * dat.sigma2 + (1.0 - f) * dat.sigma1;
      double mu = f * dat.mu2 + (1.0 - f) * dat.mu1;
      double beta = mu * std::fabs(dat.H_current);

      // -----------------------------------------------------------------
      // 3f. Effective Conductivity Tensor
      // -----------------------------------------------------------------
      double D = 1.0 + beta * beta;
      double sigma_rr = sigma / D;
      double sigma_rt = -sigma * beta / D;
      double sigma_tr = sigma * beta / D;
      double sigma_tt = sigma / D;

      // -----------------------------------------------------------------
      // 3g. Assemble Element Matrix
      // -----------------------------------------------------------------
      double r_inv = 1.0 / r_phys;
      double weight = r_phys * detJ * w * dat.t;

      for (int a = 0; a < num_field_dofs; ++a) {
        double grad_a_r = dNdr[a];
        double grad_a_t = r_inv * dNdtheta[a];

        for (int b = 0; b < num_field_dofs; ++b) {
          double grad_b_r = dNdr[b];
          double grad_b_t = r_inv * dNdtheta[b];

          // Weak Form:  (∇u) · σ · (∇v)
          // Includes Hall Terms:
          //   grad_a_r * sigma_rt * grad_b_t corresponds to:
          //   (∂N_A/∂r) * (-σ·β/D) * (1/r ∂N_B/∂θ)
          // This correctly captures the mixed derivative term essential for
          // EMR.

          double term = grad_a_r * (sigma_rr * grad_b_r + sigma_rt * grad_b_t) +
                        grad_a_t * (sigma_tr * grad_b_r + sigma_tt * grad_b_t);

          Ae[a * nloc + b] += term * weight;
        }
      }
    }
  }
}

// =============================================================================
// MAIN ASSEMBLY FUNCTION (Unchanged from your version)
// =============================================================================
PetscErrorCode make_global(global_matrices &gmat, data &dat) {
  PetscErrorCode ierr;

  // Setup node overlay for periodic BCs
  if (dat.use_annular) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "Initializing node overlay for annular domain...\n");
    CHKERRQ(ierr);
    Convertnode(dat);
  } else {
    dat.Convertnode = new int[dat.nglobal];
    for (int i = 0; i < dat.nglobal; ++i) {
      dat.Convertnode[i] = i;
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "Rectangular domain: using identity mapping\n");
    CHKERRQ(ierr);
  }

  // Create matrix and vectors
  // Use MATAIJ (auto-selects SEQ/MPI) and call MatSetUp before MatSetOption
  // to avoid segfault on Linux PETSc 3.12+ when setting options on
  // an uninitialized MPIAIJ internal structure.
  ierr = MatCreate(PETSC_COMM_WORLD, &gmat.A);
  CHKERRQ(ierr);
  ierr =
      MatSetSizes(gmat.A, PETSC_DECIDE, PETSC_DECIDE, dat.nglobal, dat.nglobal);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(gmat.A);
  CHKERRQ(ierr);
  ierr = MatSetType(gmat.A, MATAIJ);
  CHKERRQ(ierr);
  ierr = MatSetUp(gmat.A);
  CHKERRQ(ierr);
  ierr = MatSetOption(gmat.A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  CHKERRQ(ierr);
  ierr = MatSetOption(gmat.A, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = MatCreateVecs(gmat.A, NULL, &gmat.rhs);
  CHKERRQ(ierr);
  ierr = VecSet(gmat.rhs, 0.0);
  CHKERRQ(ierr);

  // Parallel element assembly
  int rank, size;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  int start = (dat.nelem * rank) / size;
  int end = (dat.nelem * (rank + 1)) / size;

  int nloc = dat.nele_mat;
  PetscScalar *Ke = new PetscScalar[nloc * nloc];

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Assembling elements %d to %d on rank %d...\n", start,
                     end - 1, rank);
  CHKERRQ(ierr);

  // Precompute constants or just enter loop

  for (int ie = start; ie < end; ++ie) {
    // Choose physics kernel
    if (dat.Nsample_R2 > 0) {
      // EMR MODE - use corrected function
      if (!dat.elem[ie]) {
        printf("FATAL: dat.elem[%d] is NULL\n", ie);
        fflush(stdout);
        exit(1);
      }
      element_stiffness_EMR(dat, ie, Ke);
    } else {
      // Standard Laplace
      int lay = dat.layer_id_for_elem[ie];
      int matID = dat.layer_material[lay];
      double kappa = (matID == 1) ? 2.0 : 1.0;
      element_stiffness(dat, ie, kappa, Ke);
    }

    // Scatter to global matrix
    for (int a = 0; a < dat.node_elem; ++a) {
      for (int c = 0; c < dat.dof_per_node; ++c) {
        int local_dof_row = dat.elem[ie][a] * dat.dof_per_node + c;
        PetscInt row = dat.Convertnode[local_dof_row];

        for (int b = 0; b < dat.node_elem; ++b) {
          for (int d = 0; d < dat.dof_per_node; ++d) {
            int local_dof_col = dat.elem[ie][b] * dat.dof_per_node + d;
            PetscInt col = dat.Convertnode[local_dof_col];

            PetscScalar val = Ke[(a * dat.dof_per_node + c) * nloc +
                                 (b * dat.dof_per_node + d)];

            ierr = MatSetValue(gmat.A, row, col, val, ADD_VALUES);
            CHKERRQ(ierr);
          }
        }

        // =============================================================
        // HERMITE FIX: Hall effect at PORT nodes
        // =============================================================
        // Apply inside scatter loop where 'row' and node info are available.
        if (dat.dof_per_node == 4 && c == 2) {
          int node_idx = dat.elem[ie][a];
          if (dat.port_code[node_idx] != 0 &&
              std::abs(dat.H_current) > 1.0e-6) {
            // Option A: Set dPhi/dTheta = 0 approx at all port nodes (simplest,
            // robust)
            ierr = MatSetValue(gmat.A, row, row, 1.0, ADD_VALUES);
            CHKERRQ(ierr);
          }
        }
      }
    }
  }

  delete[] Ke;

  // Insert placeholder diagonals for overlaid DOFs
  // Insert placeholder diagonals for ALL DOFs to prevent singular matrix
  // from orphan nodes or boundary issues.

  // Insert placeholder diagonals for ALL DOFs to prevent singular matrix
  int placeholder_count = 0;
  for (int i = 0; i < dat.nglobal; ++i) {
    PetscInt row = i;
    PetscScalar placeholder = 1.0e-9;
    ierr = MatSetValue(gmat.A, row, row, placeholder, ADD_VALUES);
    CHKERRQ(ierr);
    placeholder_count++;
  }

  // First assembly
  ierr = MatAssemblyBegin(gmat.A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(gmat.A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(gmat.rhs);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(gmat.rhs);
  CHKERRQ(ierr);

  // Apply periodic constraints
  std::vector<PetscInt> overlaid_dofs;
  std::vector<PetscInt> canonical_dofs;

  for (int i = 0; i < dat.nglobal; ++i) {
    if (dat.Convertnode[i] != i) {
      overlaid_dofs.push_back(i);
      canonical_dofs.push_back(dat.Convertnode[i]);
    }
  }

  int nconstraints = overlaid_dofs.size();
  if (nconstraints > 0) {
    ierr = MatZeroRows(gmat.A, nconstraints, overlaid_dofs.data(), 1.0, NULL,
                       NULL);
    CHKERRQ(ierr);

    for (int k = 0; k < nconstraints; ++k) {
      PetscInt top_dof = overlaid_dofs[k];
      PetscInt bottom_dof = canonical_dofs[k];

      ierr = MatSetValue(gmat.A, top_dof, bottom_dof, -1.0, INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = VecSetValue(gmat.rhs, top_dof, 0.0, INSERT_VALUES);
      CHKERRQ(ierr);
    }
  }

  // Final assembly
  ierr = MatAssemblyBegin(gmat.A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(gmat.A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(gmat.rhs);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(gmat.rhs);
  CHKERRQ(ierr);

  return 0;
}