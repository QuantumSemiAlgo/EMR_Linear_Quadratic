/**
 * @file make_global.cpp
 * @brief Global Stiffness Matrix Assembly.
 *
 * This module handles the assembly of the global stiffness matrix (A) and the
 * right-hand side vector (rhs) for the 2D Laplace equation. It supports both
 * standard rectangular domains and annular domains using a node overlay
 * technique for periodic boundary conditions.
 *
 * Key Features:
 * - Parallel assembly using MPI and PETSc.
 * - Support for periodic boundary conditions via `Convertnode` mapping.
 * - Stabilization of the matrix for periodic constraints using placeholder
 * diagonals.
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#include "global_params_1.h"
#include "prototypes_1.h"
#include "utilrect_1.h"
#include <petsc.h>
#include <petscksp.h>
#include <vector>

/**
 * @brief Assembles the global stiffness matrix and RHS vector.
 *
 * @param gmat Reference to the global matrices structure (A, rhs).
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode make_global(global_matrices &gmat, data &dat) {
  PetscErrorCode ierr;

  // ===========================================================
  // 1. Setup Node Overlay (Periodic BCs)
  // ===========================================================
  // The "Convertnode" array maps a local degree of freedom (DOF) to its
  // canonical global DOF. For standard nodes, Convertnode[i] = i.
  // For periodic nodes (e.g., θ=0 and θ=2π), they map to the same ID.
  if (dat.use_annular) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "Initializing node overlay for annular domain...\n");
    CHKERRQ(ierr);
    // Generates the mapping for periodic boundaries (θ direction).
    Convertnode(dat);
  } else {
    // Rectangular domain: Identity mapping (no overlay).
    dat.Convertnode = new int[dat.nglobal];
    for (int i = 0; i < dat.nglobal; ++i) {
      dat.Convertnode[i] = i;
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "Rectangular domain: using identity mapping\n");
    CHKERRQ(ierr);
  }

  // ===========================================================
  // 2. Create PETSc Matrix and Vectors
  // ===========================================================
  ierr = MatCreate(PETSC_COMM_WORLD, &gmat.A);
  CHKERRQ(ierr);
  ierr =
      MatSetSizes(gmat.A, PETSC_DECIDE, PETSC_DECIDE, dat.nglobal, dat.nglobal);
  CHKERRQ(ierr);
  ierr = MatSetFromOptions(gmat.A);
  CHKERRQ(ierr);
  ierr = MatSetType(gmat.A, MATMPIAIJ);
  CHKERRQ(ierr);

  // IMPORTANT: We need to modify the matrix structure later (e.g., for BCs).
  // These options prevent PETSc from throwing errors when we add new nonzeros.
  ierr = MatSetOption(gmat.A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  CHKERRQ(ierr);
  ierr = MatSetOption(gmat.A, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
  CHKERRQ(ierr);

  ierr = MatSetUp(gmat.A);
  CHKERRQ(ierr);

  ierr = MatCreateVecs(gmat.A, NULL, &gmat.rhs);
  CHKERRQ(ierr);
  ierr = VecSet(gmat.rhs, 0.0);
  CHKERRQ(ierr);

  // ===========================================================
  // 3. Parallel Element Assembly
  // ===========================================================
  int rank, size;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  CHKERRQ(ierr);

  // Determine the range of elements this process is responsible for.
  int start = (dat.nelem * rank) / size;
  int end = (dat.nelem * (rank + 1)) / size;

  int nloc = dat.nele_mat; // Local DOFs per element
  PetscScalar *Ke = new PetscScalar[nloc * nloc];

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Assembling elements %d to %d on rank %d...\n", start,
                     end - 1, rank);
  CHKERRQ(ierr);

  // Loop over local elements assigned to this process
  for (int ie = start; ie < end; ++ie) {

    // Retrieve material properties (currently uniform diffusion)
    int lay = dat.layer_id_for_elem[ie];
    int matID = dat.layer_material[lay];

    // Future-proofing: Lookup material property (C++98 safe)
    double kappa;
    switch (matID) {
    case 0:
      kappa = 1.0;
      break;
    case 1:
      kappa = 2.0;
      break;
    case 2:
      kappa = 0.5;
      break;
    default:
      kappa = 1.0;
      break;
    }

    // Compute the local element stiffness matrix (Ke)
    element_stiffness(dat, ie, kappa, Ke);

    // -------------------------------------------------------
    // Scatter to Global Matrix with Overlay Mapping
    // -------------------------------------------------------
    // We iterate over the local element matrix and add values to the
    // global matrix. The `Convertnode` array ensures that contributions
    // from periodic nodes are accumulated into the same global row/col.
    for (int a = 0; a < dat.node_elem; ++a) {
      for (int c = 0; c < dat.dof_per_node; ++c) {

        // Local DOF index (before overlay)
        int local_dof_row = dat.elem[ie][a] * dat.dof_per_node + c;

        // Global DOF index (after overlay mapping)
        PetscInt row = dat.Convertnode[local_dof_row];

        for (int b = 0; b < dat.node_elem; ++b) {
          for (int d = 0; d < dat.dof_per_node; ++d) {

            int local_dof_col = dat.elem[ie][b] * dat.dof_per_node + d;
            PetscInt col = dat.Convertnode[local_dof_col];

            PetscScalar val = Ke[(a * dat.dof_per_node + c) * nloc +
                                 (b * dat.dof_per_node + d)];

            // ADD_VALUES: accumulates contributions from overlaid nodes
            ierr = MatSetValue(gmat.A, row, col, val, ADD_VALUES);
            CHKERRQ(ierr);
          }
        }
      }
    }
  }

  delete[] Ke;

  // ===========================================================
  // 4. Insert Placeholder Diagonals
  // ===========================================================
  // CRITICAL: For overlaid DOFs (the "slave" nodes in a periodic pair),
  // no element contributions are assembled into their rows because they
  // are remapped to the "master" node. This leaves their rows empty.
  // PETSc's MatZeroRows (used later for BCs) requires diagonal entries to
  // exist. We insert a tiny non-zero value to prevent "missing diagonal"
  // errors.
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Inserting placeholder diagonals for overlaid DOFs...\n");
  CHKERRQ(ierr);

  int placeholder_count = 0;
  for (int i = 0; i < dat.nglobal; ++i) {
    if (dat.Convertnode[i] != i) {
      // This DOF is overlaid (slave) - insert placeholder
      PetscInt row = i;
      PetscScalar placeholder = 1.0e-30; // Arbitrary small value
      ierr = MatSetValue(gmat.A, row, row, placeholder, ADD_VALUES);
      CHKERRQ(ierr);
      placeholder_count++;
    }
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Inserted %d placeholder diagonals\n",
                     placeholder_count);
  CHKERRQ(ierr);

  // ===========================================================
  // 5. First Assembly
  // ===========================================================
  // Commit all element contributions and placeholders to the matrix.
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "First assembly (committing element contributions)...\n");
  CHKERRQ(ierr);

  ierr = MatAssemblyBegin(gmat.A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(gmat.A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  ierr = VecAssemblyBegin(gmat.rhs);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(gmat.rhs);
  CHKERRQ(ierr);

  // ===========================================================
  // 6. Apply Periodic Constraints (Explicitly)
  // ===========================================================
  // While the overlay accumulates stiffness into the master node, we still
  // need to constrain the slave node to equal the master node:
  // u_slave - u_master = 0

  std::vector<PetscInt> overlaid_dofs;
  std::vector<PetscInt> canonical_dofs;

  for (int i = 0; i < dat.nglobal; ++i) {
    if (dat.Convertnode[i] != i) {
      overlaid_dofs.push_back(i);
      canonical_dofs.push_back(dat.Convertnode[i]);
    }
  }

  int nconstraints = overlaid_dofs.size();
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Found %d overlaid DOFs to constrain\n",
                     nconstraints);
  CHKERRQ(ierr);

  if (nconstraints > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Applying periodic constraints...\n");
    CHKERRQ(ierr);

    // Step A: Zero out the rows corresponding to slave nodes.
    // This sets the diagonal to 1.0 and RHS to 0.0 by default.
    ierr = MatZeroRows(gmat.A, nconstraints, overlaid_dofs.data(), 1.0, NULL,
                       NULL);
    CHKERRQ(ierr);

    // Step B: Add the coupling term to enforce u_slave = u_master.
    // Equation: 1.0 * u_slave - 1.0 * u_master = 0
    // We use INSERT_VALUES because MatZeroRows cleared the row.
    for (int k = 0; k < nconstraints; ++k) {
      PetscInt top_dof = overlaid_dofs[k];     // Slave
      PetscInt bottom_dof = canonical_dofs[k]; // Master

      ierr = MatSetValue(gmat.A, top_dof, bottom_dof, -1.0, INSERT_VALUES);
      CHKERRQ(ierr);

      ierr = VecSetValue(gmat.rhs, top_dof, 0.0, INSERT_VALUES);
      CHKERRQ(ierr);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD, "Applied %d periodic constraints\n",
                       nconstraints);
    CHKERRQ(ierr);
  }

  // ===========================================================
  // 7. Final Assembly
  // ===========================================================
  // Commit the constraint changes.
  ierr = MatAssemblyBegin(gmat.A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(gmat.A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  ierr = VecAssemblyBegin(gmat.rhs);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(gmat.rhs);
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Assembly complete (final)\n");
  CHKERRQ(ierr);

  return 0;
}