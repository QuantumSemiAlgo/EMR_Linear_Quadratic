/**
 * @file convertnode.cpp
 * @brief Node Overlay for Periodic Boundary Conditions.
 *
 * This module implements the "node overlay" technique used to enforce periodic
 * boundary conditions in annular domains. It identifies nodes on the "top" edge
 * (theta = 2pi) and maps them to the corresponding nodes on the "bottom" edge
 * (theta = 0).
 *
 * The `Convertnode` array stores this mapping:
 * - Convertnode[i] = i  (for standard nodes)
 * - Convertnode[top_node] = bottom_node (for periodic nodes)
 *
 * This allows the solver to treat the two boundaries as physically connected.
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#include "global_params_1.h"
#include <algorithm>
#include <cmath>
#include <vector>

/**
 * @brief Generates the node mapping array for periodic boundaries.
 *
 * Scans the mesh to identify nodes on the top and bottom boundaries (y-limits).
 * Pairs them based on their x-coordinates and populates the `Convertnode`
 * array.
 *
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Convertnode(data &dat) {
  PetscErrorCode ierr;

  // ---------------------------------------------------------------------------
  // 1. Initialization
  // ---------------------------------------------------------------------------
  // Allocate conversion array
  dat.Convertnode = new int[dat.nglobal];

  // Initialize: each DOF maps to itself by default (Identity mapping)
  for (int i = 0; i < dat.nglobal; ++i) {
    dat.Convertnode[i] = i;
  }

  // If not annular, we don't need overlay (rectangular domain uses different
  // BCs)
  if (!dat.use_annular) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "Rectangular domain: no node overlay needed\n");
    CHKERRQ(ierr);
    return 0;
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Setting up node overlay for annular domain...\n");
  CHKERRQ(ierr);

  int stride = dat.dof_per_node;

  // ---------------------------------------------------------------------------
  // 2. Determine Tolerances
  // ---------------------------------------------------------------------------
  // Use mesh spacing to set robust geometric tolerances.
  double tol_x = 1e-6;
  double tol_y = 1e-6;

  // Use mesh density if available to set reasonable tolerances
  if (dat.n_x > 0)
    tol_x = (dat.xmax - dat.xmin) / (double)dat.n_x * 0.1;
  if (dat.n_y > 0)
    tol_y = (dat.ymax - dat.ymin) / (double)dat.n_y * 0.1;

  // Fallback if spacing info is missing
  // Scale to 0.1% of domain size for robustness
  if (tol_x <= 0.0)
    tol_x = (dat.xmax - dat.xmin) * 0.001;
  if (tol_y <= 0.0)
    tol_y = (dat.ymax - dat.ymin) * 0.001;

  // ADD THIS DEBUG BLOCK
  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Tolerance: tol_x=%.6e, tol_y=%.6e\n",
                     tol_x, tol_y);
  CHKERRQ(ierr);
  ierr =
      PetscPrintf(PETSC_COMM_WORLD,
                  "  Domain bounds from dat: x=[%.6f, %.6f], y=[%.6f, %.6f]\n",
                  dat.xmin, dat.xmax, dat.ymin, dat.ymax);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 3. Find Domain Extents
  // ---------------------------------------------------------------------------
  // Scan all nodes to find the actual physical bounds of the mesh.
  double actual_ymin = 1e10;
  double actual_ymax = -1e10;
  double actual_xmin = 1e10;
  double actual_xmax = -1e10;

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
                     "  Actual domain: x=[%.6f, %.6f], y=[%.6f, %.6f]\n",
                     actual_xmin, actual_xmax, actual_ymin, actual_ymax);
  CHKERRQ(ierr);

  // ===========================================================================
  // 4. Identify Candidate Nodes
  // ===========================================================================
  // Build lists of nodes that lie on the top (y ~ ymax) and bottom (y ~ ymin)
  // edges.

  int overlay_count = 0;
  std::vector<int> bottom_nodes;
  std::vector<int> top_nodes;

  for (int n = 0; n < dat.ngnodes; ++n) {
    double y_n = dat.node[n][1];
    if (std::fabs(y_n - actual_ymin) < tol_y)
      bottom_nodes.push_back(n);
    if (std::fabs(y_n - actual_ymax) < tol_y)
      top_nodes.push_back(n);
  }

  // Copy to candidate lists (we might filter them later if needed, but for now
  // keep all)
  std::vector<int> top_candidates = top_nodes;
  std::vector<int> bottom_candidates = bottom_nodes;

  // ===========================================================================
  // 5. Sort Candidates
  // ===========================================================================
  // Sort both lists by x-coordinate. This ensures that the i-th node on the
  // bottom corresponds to the i-th node on the top, assuming a structured mesh.

  struct XCmp {
    double **nodeptr;
    XCmp(double **n) : nodeptr(n) {}
    bool operator()(int a, int b) const {
      return nodeptr[a][0] < nodeptr[b][0];
    }
  };

  std::sort(bottom_candidates.begin(), bottom_candidates.end(), XCmp(dat.node));
  std::sort(top_candidates.begin(), top_candidates.end(), XCmp(dat.node));

  // Check for mismatch in node counts
  int npair =
      std::min((int)bottom_candidates.size(), (int)top_candidates.size());
  if (npair == 0) {
    ierr = PetscPrintf(
        PETSC_COMM_WORLD,
        "  WARNING: no top/bottom edge node pairs found (bottom=%d, top=%d)\n",
        (int)bottom_candidates.size(), (int)top_candidates.size());
    CHKERRQ(ierr);
  }

  // ===========================================================================
  // 6. Pair Nodes and Update Mapping
  // ===========================================================================
  for (int i = 0; i < npair; ++i) {
    int top_n = top_candidates[i];
    int bot_n = bottom_candidates[i];

    double x_top = dat.node[top_n][0];
    double y_top = dat.node[top_n][1];
    double x_bot = dat.node[bot_n][0];
    double y_bot = dat.node[bot_n][1];

    // Verify x-alignment (should match for periodic boundaries)
    double diff = std::fabs(x_top - x_bot);
    bool within_tol = (diff <= tol_x);

    if (!within_tol) {
      ierr =
          PetscPrintf(PETSC_COMM_WORLD,
                      "  WARNING: Pair %d has large x-diff=%.6e (tol=%.6e)\n",
                      i, diff, tol_x);
      CHKERRQ(ierr);
    }

    // Map all DOFs of the top node to the bottom node
    for (int d = 0; d < stride; ++d) {
      int top_dof = top_n * stride + d;
      int bottom_dof = bot_n * stride + d;

      if (top_dof >= 0 && top_dof < dat.nglobal && bottom_dof >= 0 &&
          bottom_dof < dat.nglobal) {
        dat.Convertnode[top_dof] = bottom_dof;
      }
    }
    overlay_count++;

    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "  Overlay: top node %d (x=%.6f,y=%.6f) → bottom node "
                       "%d (x=%.6f,y=%.6f) diff=%.6e %s\n",
                       top_n, x_top, y_top, bot_n, x_bot, y_bot, diff,
                       within_tol ? "" : "(outside tol)");
    CHKERRQ(ierr);
  }

  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "Node overlay complete: %d nodes overlaid (of %d top-edge candidates)\n",
      overlay_count, (int)top_candidates.size());
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 7. Validation
  // ---------------------------------------------------------------------------
  // Check for any out-of-bounds mappings.
  int error_count = 0;
  for (int i = 0; i < dat.nglobal; ++i) {
    if (dat.Convertnode[i] < 0 || dat.Convertnode[i] >= dat.nglobal) {
      error_count++;
      if (error_count <= 3) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                           "  ERROR: Invalid mapping Convertnode[%d] = %d\n", i,
                           dat.Convertnode[i]);
        CHKERRQ(ierr);
      }
    }
  }

  if (error_count > 0) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
            "Found %d invalid Convertnode mappings", error_count);
  }

  return 0;
}