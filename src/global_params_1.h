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

  // =========================================================================
  // EMR-SPECIFIC GEOMETRY
  // =========================================================================
  double R1;           ///< Outer semiconductor radius (fixed)
  double R2;           ///< Inner metal radius (variable sweep parameter)
  double R2min, R2max; ///< Range for R2 sweep
  double dR2;          ///< Step size for R2 sweep
  int Nsample_R2;      ///< Number of R2 values to simulate

  // =========================================================================
  // MATERIAL PROPERTIES
  // =========================================================================
  double sigma1;     ///< Conductivity of semiconductor (outer)
  double sigma2;     ///< Conductivity of metal (inner)
  double mu1;        ///< Mobility of semiconductor
  double mu2;        ///< Mobility of metal
  double n1, n2;     ///< Carrier density (semiconductor, metal)
  double tau1, tau2; ///< Collision time
  double m1, m2;     ///< Effective mass
  double delta;      ///< Width parameter for Fermi function transition

  // =========================================================================
  // MAGNETIC FIELD SWEEP
  // =========================================================================
  double Hmin;      ///< Minimum magnetic field (Tesla)
  double Hmax;      ///< Maximum magnetic field (Tesla)
  double dH;        ///< Step size for H sweep
  int Nsample_H;    ///< Number of H field values
  double H_current; ///< Current H value in simulation

  // =========================================================================
  // ELECTRICAL PARAMETERS
  // =========================================================================
  double Io;                ///< Input current (Amperes)
  double t;                 ///< Disk thickness (meters)
  double Ro = 0.0;                ///< Reference resistance (at H=0)
  double latest_resistance = 0.0; ///< Most recently calculated resistance

  // =========================================================================
  // PORT GEOMETRY
  // =========================================================================
  double theta1, theta2, theta3, theta4; ///< Port angular positions (x PI)
  double width_L1, width_L2, width_L3,
      width_L4; ///< Port angular widths (x PI)
  int iela1, ielb1, iela2, ielb2, iela3, ielb3, iela4,
      ielb4;                                          ///< Port element indices
  int ina1, inb1, ina2, inb2, ina3, inb3, ina4, inb4; ///< Port node indices
  int *port_code; ///< Array[ngnodes_t]: identifies which port each theta node
                  ///< belongs to
                  ///< 0=none, 1=port1, 2=port2, 3=port3, 4=port4

  // =========================================================================
  // MESH REFINEMENT AROUND R2
  // =========================================================================
  double width_R2;  ///< Radial width of refined zone around R2
  double width_Rs2; ///< Width of side transition zones
  int nelem_R2;     ///< Elements in main R2 zone
  int nelem_Rs2;    ///< Elements in side zones
  int nelem_L;      ///< Elements per port (angular)
  int nelem_otherR; ///< Elements in other radial zones
  int nelem_otherT; ///< Elements in other angular zones
};

/**
 * @struct global_matrices
 * @brief Container for PETSc global linear algebra objects.
 */
struct global_matrices {
  Mat A;   ///< Global stiffness matrix (distributed sparse matrix)
  Vec rhs; ///< Global right-hand side vector (distributed)
};

/**
 * @brief Get element type string based on node_elem and dof_per_node
 * @param node_elem Number of nodes per element
 * @param dof_per_node Degrees of freedom per node
 * @return String describing the element type
 */
inline const char* get_element_type_string(int node_elem, int dof_per_node) {
  if (dof_per_node == 1) {
    // Lagrangian elements
    if (node_elem == 4) return "lagrangian_linear";   // Q1
    if (node_elem == 9) return "lagrangian_quadratic"; // Q2
  } else if (dof_per_node == 4) {
    // Hermite elements
    if (node_elem == 4) return "hermite_cubic";   // Bicubic Hermite
    if (node_elem == 9) return "hermite_quintic"; // Biquintic Hermite
  }
  return "unknown"; // Fallback for unsupported combinations
}

#endif // GLOBAL_PARAMS_H