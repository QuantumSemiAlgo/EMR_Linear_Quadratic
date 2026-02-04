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

  // If not annular, we don't need overlay (unless specified otherwise)
  if (!dat.use_annular) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "Rectangular domain: no node overlay needed\n");
    CHKERRQ(ierr);
    return 0;
  }

  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "Setting up node overlay for annular domain (Theta Periodicity)...\n");
  CHKERRQ(ierr);

  int stride = dat.dof_per_node;

  // ---------------------------------------------------------------------------
  // 2. Determine Tolerances
  // ---------------------------------------------------------------------------
  double tol_x = 1e-6;
  double tol_y = 1e-6;

  // Set tolerances based on domain size if available
  if (dat.n_x > 0)
    tol_x = (dat.xmax - dat.xmin) / (double)dat.n_x * 0.1;
  if (dat.n_y > 0)
    tol_y = (dat.ymax - dat.ymin) / (double)dat.n_y * 0.1;

  if (tol_x <= 0.0)
    tol_x = 1e-5;
  if (tol_y <= 0.0)
    tol_y = 1e-5;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Tolerance: tol_x=%.6e, tol_y=%.6e\n",
                     tol_x, tol_y);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 3. Find Domain Extents
  // ---------------------------------------------------------------------------
  double actual_xmin = 1e10, actual_xmax = -1e10;
  double actual_ymin = 1e10, actual_ymax = -1e10;

  for (int i = 0; i < dat.ngnodes; ++i) {
    double x = dat.node[i][0]; // Theta
    double y = dat.node[i][1]; // Radius

    if (x < actual_xmin)
      actual_xmin = x;
    if (x > actual_xmax)
      actual_xmax = x;
    if (y < actual_ymin)
      actual_ymin = y;
    if (y > actual_ymax)
      actual_ymax = y;
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "  Actual domain: Theta=[%.6f, %.6f], R=[%.6f, %.6f]\n",
                     actual_xmin, actual_xmax, actual_ymin, actual_ymax);
  CHKERRQ(ierr);

  // ===========================================================================
  // 4. Identify Periodic Nodes (Theta Boundaries)
  // ===========================================================================
  // In Annular mode: x is Theta, y is Radius.
  // Periodicity maps x_max (2pi) -> x_min (0).
  // We look for nodes where x ~ x_max vs x ~ x_min.

  std::vector<int> left_nodes;  // Theta ~ Min
  std::vector<int> right_nodes; // Theta ~ Max

  // Use slightly larger tolerance for identifying boundary membership
  double boundary_tol = std::max(tol_x, 1e-4);

  for (int n = 0; n < dat.ngnodes; ++n) {
    double x_n = dat.node[n][0];
    if (std::fabs(x_n - actual_xmin) < boundary_tol)
      left_nodes.push_back(n);
    if (std::fabs(x_n - actual_xmax) < boundary_tol)
      right_nodes.push_back(n);
  }

  // ===========================================================================
  // 5. Match Nodes by Y-Coordinate (Radius)
  // ===========================================================================
  // Sort by Y (Radius) to facilitate matching
  struct YCmp {
    double **nodeptr;
    YCmp(double **n) : nodeptr(n) {}
    bool operator()(int a, int b) const {
      return nodeptr[a][1] < nodeptr[b][1];
    }
  };

  std::sort(left_nodes.begin(), left_nodes.end(), YCmp(dat.node));
  std::sort(right_nodes.begin(), right_nodes.end(), YCmp(dat.node));

  int match_count = 0;
  size_t i_left = 0;
  size_t i_right = 0;

  // Linear scan to match (since sorted by Y)
  while (i_left < left_nodes.size() && i_right < right_nodes.size()) {
    int n_l = left_nodes[i_left];
    int n_r = right_nodes[i_right];

    double y_l = dat.node[n_l][1];
    double y_r = dat.node[n_r][1];

    if (std::fabs(y_l - y_r) < tol_y) {
      // MATCH FOUND: Map RIGHT node (source) to LEFT node (canonical)
      // Convention: Map High Theta -> Low Theta

      // Update Convertnode for all DOFs
      for (int d = 0; d < stride; ++d) {
        int gdof_r = n_r * stride + d;
        int gdof_l = n_l * stride + d;
        dat.Convertnode[gdof_r] = gdof_l;
      }

      // Debug print for first few
      if (match_count < 5) {
        PetscPrintf(PETSC_COMM_WORLD,
                    "  Overlay: Node %d (Theta=2pi, R=%.5f) -> Node %d "
                    "(Theta=0, R=%.5f)\n",
                    n_r, y_r, n_l, y_l);
      }

      match_count++;
      i_left++;
      i_right++;
    } else if (y_l < y_r) {
      i_left++;
    } else {
      i_right++;
    }
  }

  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "Node overlay complete: %d periodic pairs mapped (Left=%lu, Right=%lu)\n",
      match_count, left_nodes.size(), right_nodes.size());
  CHKERRQ(ierr);

  return 0;
}