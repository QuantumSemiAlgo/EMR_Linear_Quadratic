/**
 * @file solution_clean.cpp
 * @brief Solution Post-Processing and Output.
 *
 * This module handles the reconstruction of the solution field from the
 * finite element coefficients. It interpolates the solution onto a regular
 * grid for visualization and handles the coordinate transformation for
 * annular domains (mapping logical (x,y) to physical (r,θ) to Cartesian (X,Y)).
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#include "global_params_1.h"
#include "prototypes_1.h"
#include "utilrect_1.h"
#include <fstream>
#include <iomanip>
#include <petsc.h>
#include <vector>

/**
 * @brief Interpolates the solution onto a regular grid and writes to file.
 *
 * This function takes the solution vector (psi), scatters it to the root
 * process (if not already done), and then iterates over a dense grid of points.
 * For each point, it:
 * 1. Locates the containing element.
 * 2. Computes the shape functions at the local coordinates.
 * 3. Interpolates the solution value.
 * 4. Transforms coordinates if the domain is annular.
 * 5. Writes the result to `solution_recon.dat`.
 *
 * @param gmat Global matrices structure.
 * @param dat Global data structure.
 * @param psi Solution vector (gathered on root).
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode solution(global_matrices &gmat, data &dat, Vec psi) {
  PetscErrorCode ierr;
  const PetscScalar *u;

  // Access the solution array (read-only)
  ierr = VecGetArrayRead(psi, &u);
  CHKERRQ(ierr);

  // Only rank 0 performs the file I/O
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank != 0) {
    ierr = VecRestoreArrayRead(psi, &u);
    CHKERRQ(ierr);
    return 0;
  }

  // Open output file
  std::string reconfile = std::string(dat.outpath) + "solution_recon.dat";
  std::ofstream fout(reconfile);

  // Write header
  if (dat.use_annular) {
    fout << "#   x_cart    y_cart    r    theta    value\n";
  } else {
    fout << "#   zx        zy        value\n";
  }
  fout << std::scientific << std::setprecision(6);

  // Pre-allocate working vectors
  int nele_mat = dat.node_elem * dat.dof_per_node;
  std::vector<double> phi(nele_mat);
  std::vector<PetscScalar> uel(nele_mat);

  // Store the first row of values to enforce periodicity in the output
  // (i.e., ensure the value at θ=2π exactly matches θ=0)
  std::vector<double> first_row_uh;
  if (dat.use_annular)
    first_row_uh.resize(dat.ndz_x + 1, 0.0);

  // Loop limits
  int j_start = 0;
  int j_end = dat.ndz_y;
  // For annular, we stop one step short and manually duplicate the first row
  // to the last row to guarantee perfect visual periodicity.
  // Note: Periodicity is enforced in-place at j==dat.ndz_y (see line 125)

  // ---------------------------------------------------------------------------
  // Grid Loop
  // ---------------------------------------------------------------------------
  for (int j = j_start; j <= j_end; ++j) {
    double zy = dat.ymin + j * dat.dzy;

    for (int i = 0; i <= dat.ndz_x; ++i) {
      double zx = dat.xmin + i * dat.dzx;

      // 1. Find the element containing (zx, zy)
      int ieFound;
      ierr = locelem(zx, zy, dat, ieFound);
      CHKERRQ(ierr);

      // 2. Compute local coordinates (xi, eta)
      // Note: This assumes a rectangular regular mesh structure.
      int nn1 = dat.elem[ieFound][0]; // Bottom-Left
      int nn2 = dat.elem[ieFound][1]; // Bottom-Right
      int nn4 = dat.elem[ieFound][3]; // Top-Left

      double x1 = dat.node[nn1][0], y1 = dat.node[nn1][1];
      double x2 = dat.node[nn2][0], y2 = dat.node[nn2][1];
      double x4 = dat.node[nn4][0], y4 = dat.node[nn4][1];

      double xi = 2.0 * (zx - x1) / (x2 - x1) - 1.0;
      double eta = 2.0 * (zy - y1) / (y4 - y1) - 1.0;

      // 3. Evaluate shape functions
      shapeRect(phi.data(), xi, eta, dat);

      // 4. Gather element nodal values
      for (int a = 0; a < dat.node_elem; ++a) {
        for (int c = 0; c < dat.dof_per_node; ++c) {
          int local_dof = dat.elem[ieFound][a] * dat.dof_per_node + c;
          // Map to global ID using Convertnode (handles periodicity)
          int gid = dat.Convertnode[local_dof];
          uel[a * dat.dof_per_node + c] = u[gid];
        }
      }

      // 5. Interpolate solution value u_h
      double uh = 0.0;
      for (int k = 0; k < nele_mat; ++k)
        uh += phi[k] * uel[k];

      // PERIODICITY FIX: Enforce exact periodicity at theta=2pi boundary
      if (dat.use_annular) {
        if (j == 0) {

          // First row: store for later
          first_row_uh[i] = uh;
        } else if (j == dat.ndz_y) {
          // Last row: force to match first row exactly
          uh = first_row_uh[i];
        }
      }

      // 6. Write Output
      if (dat.use_annular) {
        // Coordinate Transformation:
        // x (logical) -> r (physical)
        // y (logical) -> theta (physical)
        double r_param = zx;
        double theta_param = zy;

        double r_phys = dat.Rin + (r_param - dat.xmin) / (dat.xmax - dat.xmin) *
                                      (dat.Rout - dat.Rin);
        double theta_phys =
            2.0 * PI * (theta_param - dat.ymin) / (dat.ymax - dat.ymin);

        double x_cart = r_phys * std::cos(theta_phys);
        double y_cart = r_phys * std::sin(theta_phys);

        fout << x_cart << "  " << y_cart << "  " << r_phys << "  " << theta_phys
             << "  " << uh << "\n";

        if (j == 0) {
          // first_row_uh already stored in interpolation loop
        }
      } else {
        fout << zx << "  " << zy << "  " << uh << "\n";
      }
    }
    fout << "\n"; // Blank line for Gnuplot grid structure
  }

  // ---------------------------------------------------------------------------
  // Append Periodic Boundary (Annular Only)
  // ---------------------------------------------------------------------------
  // Output block removed (integrated into main loop)

  fout.close();
  ierr = VecRestoreArrayRead(psi, &u);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Wrote reconstructed solution to %s\n",
                     reconfile.c_str());
  CHKERRQ(ierr);
  return 0;
}
