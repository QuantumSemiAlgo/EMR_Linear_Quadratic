/**
 * @file global_params.h
 * @brief Global Data Structures.
 *
 * This header defines the core data structures used to store the simulation
 * state, mesh data, solver configuration, and global matrices. It serves as the
 * central repository for all runtime parameters.
 *
 * Key Structures:
 * - `data`: Holds all mesh topology, physical parameters, and solver settings.
 * - `global_matrices`: Holds the PETSc matrix (A) and RHS vector (b).
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#ifndef GLOBAL_PARAMS_H
#define GLOBAL_PARAMS_H

#include "constants_1.h"
#include "prototypes_1.h" // Forward declarations
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <petsc.h>
#include <slepc.h>
#include <vector>

/**
 * @struct data
 * @brief Main data structure for the FEM simulation.
 *
 * Stores everything from MPI rank info to mesh connectivity, boundary
 * conditions, and solver tolerances. This struct is passed by reference to
 * almost all functions.
 */
struct data {
  // =========================================================================
  // SYSTEM & MPI
  // =========================================================================
  int ndebug;                     ///< Debug level (0=none, 1=info, 2=verbose)
  char xdate[PETSC_MAX_PATH_LEN]; ///< Date string
  char xtime[PETSC_MAX_PATH_LEN]; ///< Time string
  PetscMPIInt rank;               ///< MPI rank of this process
  PetscMPIInt size;               ///< Total number of MPI processes

  // =========================================================================
  // PROBLEM SIZE
  // =========================================================================
  int nglobal;   ///< Total number of global DOFs
  int ngnodes;   ///< Total number of mesh nodes
  int nelem;     ///< Total number of elements
  int nbelem;    ///< Total number of boundary elements
  int *bnode_id; ///< Array of boundary-node indices (0-based)

  // =========================================================================
  // BOUNDARY CONDITIONS
  // =========================================================================
  int nbnode;      ///< Number of boundary nodes
  double **bnode;  ///< [nbnode][2] Physical coordinates of each boundary node
  double **bvalue; ///< [nbnode][dof_per_node] Prescribed Dirichlet values
  int *btob; ///< [nbnode] Boundary type flags (0 = Dirichlet, others custom)
  int nboundary; ///< Number of distinct boundaries (unused)

  // =========================================================================
  // ELEMENT TOPOLOGY & DOFS
  // =========================================================================
  int node_elem;    ///< Nodes per element (4 for Q1, 9 for Q2)
  int dof_per_node; ///< DOFs per node (1 for scalar, 4 for Hermite)
  int nband;        ///< Bandwidth parameter (unused for sparse solvers)
  int nele_mat;     ///< Local DOFs per element = node_elem * dof_per_node
  int ndeg;         ///< Polynomial degree index (nele_mat - 1)
  int ndof;         ///< Total DOFs (redundant with nglobal in some contexts)

  // =========================================================================
  // PHYSICAL DOMAIN EXTENTS
  // =========================================================================
  double xmin, xmax; ///< Domain limits in X (logical)
  double ymin, ymax; ///< Domain limits in Y (logical)

  // =========================================================================
  // ANNULAR DOMAIN SUPPORT
  // =========================================================================
  bool use_annular; ///< Flag: true for circular/annular domains
  double Rin;       ///< Inner radius (physical units)
  double Rout;      ///< Outer radius (physical units)
  int *Convertnode; ///< [nglobal] Mapping for periodic BCs (Overlay:
                    ///< slave->master)

  // =========================================================================
  // MESH GENERATION PARAMETERS
  // =========================================================================
  int n_x;                ///< Number of elements in x-direction
  int n_y;                ///< Number of elements in y-direction
  int n_layers;           ///< Number of material layers
  int auto_generate_mesh; ///< Flag: 1=generate mesh, 0=load from files
  double stretch_factor;  ///< Mesh stretching factor (1.0 = uniform)

  // =========================================================================
  // MESH DATA ARRAYS
  // =========================================================================
  double **node;   ///< [ngnodes][2] Nodal coordinates
  int **elem;      ///< [nelem][node_elem] Element connectivity table
  int *material;   ///< [nelem] Material ID for each element
  int **belem;     ///< [nbelem][2] Boundary element connectivity
  int **bmaterial; ///< [nbelem][2] Boundary material/tag

  // =========================================================================
  // LAYER MAPPING (Optional)
  // =========================================================================
  int *layer_id_for_elem; ///< [nelem] Layer index for each element
  int nlayer;             ///< Number of defined layers
  int *layer_id;          ///< [nlayer] Layer IDs
  double *layer_ymin;     ///< [nlayer] Minimum Y for layer
  double *layer_ymax;     ///< [nlayer] Maximum Y for layer
  int *layer_material;    ///< [nlayer] Material ID for layer

  // =========================================================================
  // QUADRATURE
  // =========================================================================
  int ngaus;       ///< Total number of Gauss points per element
  double *xigaus;  ///< [ngaus] Gaussian integration points (xi)
  double *etagaus; ///< [ngaus] Gaussian integration points (eta)
  double *wgaus;   ///< [ngaus] Gaussian weights

  // =========================================================================
  // OUTPUT GRID
  // =========================================================================
  int ndz_x;  ///< Number of interpolation points along X
  int ndz_y;  ///< Number of interpolation points along Y
  double dzx; ///< Grid spacing in X
  double dzy; ///< Grid spacing in Y

  // =========================================================================
  // SOLVER SETTINGS
  // =========================================================================
  char kspType[PETSC_MAX_PATH_LEN]; ///< KSP solver type (e.g., "gmres", "cg")
  double rtol;                      ///< Relative convergence tolerance
  double atol;                      ///< Absolute convergence tolerance
  double divtol;                    ///< Divergence tolerance
  int maxit;                        ///< Maximum number of iterations

  // =========================================================================
  // I/O PATHS
  // =========================================================================
  char outpath[PETSC_MAX_PATH_LEN]; ///< Directory for output files
};

/**
 * @struct global_matrices
 * @brief Container for PETSc global linear algebra objects.
 */
struct global_matrices {
  Mat A;   ///< Global stiffness matrix (distributed sparse matrix)
  Vec rhs; ///< Global right-hand side vector (distributed)
};

#endif // GLOBAL_PARAMS_H