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
PetscErrorCode apply_bc(global_matrices &gmat, data &dat) {
  PetscErrorCode ierr;
  PetscInt stride = dat.dof_per_node;

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