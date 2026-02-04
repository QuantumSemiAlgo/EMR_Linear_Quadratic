/**
 * @file main.cpp
 * @brief Parallel Finite Element Solver for the 2D Laplace Equation.
 *
 * This is the main entry point for the application. It orchestrates the entire
 * Finite Element Analysis (FEA) workflow, including:
 * 1.  Initialization of PETSc and MPI environments.
 * 2.  Reading runtime parameters from an input file.
 * 3.  Reading mesh data (nodes, elements, boundary conditions).
 * 4.  Assembling the global stiffness matrix and load vector.
 * 5.  Applying boundary conditions (Dirichlet/Neumann).
 * 6.  Solving the linear system using PETSc/SLEPc solvers.
 * 7.  Post-processing and cleanup.
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#include "global_params_1.h"
#include "mesh_generator_1.h"
#include "prototypes_1.h"
#include <ctime>
#include <slepceps.h>

// Add these for directory creation
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

// -----------------------------------------------------------------------------
// Global Help Message
// -----------------------------------------------------------------------------
// This string is displayed when the user runs the program with the -help
// option. It provides a brief description of the program's purpose.
// -----------------------------------------------------------------------------
// Global Help Message
// -----------------------------------------------------------------------------
static char help[] = "Parallel FEM for 2D Laplace / EMR Simulation\n\n";

// Helper to free dynamic mesh memory
void FreeMesh(data &dat) {
  if (dat.Convertnode) {
    delete[] dat.Convertnode;
    dat.Convertnode = nullptr;
  }
  if (dat.bnode_id) {
    delete[] dat.bnode_id;
    dat.bnode_id = nullptr;
  }
  if (dat.node) {
    for (int i = 0; i < dat.ngnodes; ++i)
      delete[] dat.node[i];
    delete[] dat.node;
    dat.node = nullptr;
  }
  if (dat.elem) {
    for (int i = 0; i < dat.nelem; ++i)
      delete[] dat.elem[i];
    delete[] dat.elem;
    dat.elem = nullptr;
  }
  if (dat.bnode) {
    for (int i = 0; i < dat.nbnode; ++i)
      delete[] dat.bnode[i];
    delete[] dat.bnode;
    dat.bnode = nullptr;
  }
  if (dat.bvalue) {
    for (int i = 0; i < dat.nbnode; ++i)
      delete[] dat.bvalue[i];
    delete[] dat.bvalue;
    dat.bvalue = nullptr;
  }
  if (dat.port_code) {
    delete[] dat.port_code;
    dat.port_code = nullptr;
  }
  // Note: material arrays, layer arrays etc might be mesh dependent too.
  // For safety, clear layer mappings.
  if (dat.layer_id_for_elem) {
    delete[] dat.layer_id_for_elem;
    dat.layer_id_for_elem = nullptr;
  }

  // Reset counts
  dat.ngnodes = 0;
  dat.nelem = 0;
  dat.nbnode = 0;
}

// =============================================================================
// PORT IDENTIFICATION
// =============================================================================
PetscErrorCode identify_ports(data &dat) {
  PetscErrorCode ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Identifying Ports for EMR BCs...\n");
  CHKERRQ(ierr);

  // Allocate port_code
  if (dat.port_code)
    delete[] dat.port_code;
  dat.port_code = new int[dat.ngnodes];
  for (int i = 0; i < dat.ngnodes; ++i)
    dat.port_code[i] = 0;

  double w_ang = (10.0 * M_PI / 180.0); // 10 degrees half-width
  if (dat.width_L1 > 1e-9)
    w_ang = dat.width_L1 * M_PI; // If read

  // Angles
  double ang1 = 0.0;
  double ang2 = M_PI;
  double ang3 = M_PI / 2.0;
  double ang4 = 3.0 * M_PI / 2.0;

  if (dat.theta1 > 1e-9 || std::abs(dat.theta1) > 1e-9)
    ang1 = dat.theta1 * M_PI;
  if (dat.theta2 > 1e-9)
    ang2 = dat.theta2 * M_PI;
  if (dat.theta3 > 1e-9)
    ang3 = dat.theta3 * M_PI;

  int c1 = 0, c2 = 0, c3 = 0, c4 = 0;

  // Loop nodes
  // Only check OUTER boundary nodes for ports?
  // Normally contacts are on periphery.
  double tol = 1e-5;

  for (int i = 0; i < dat.ngnodes; ++i) {
    double t = dat.node[i][0]; // theta
    double r = dat.node[i][1]; // r

    // Normalize t to 0..2pi
    while (t < 0)
      t += 2 * M_PI;
    while (t >= 2 * M_PI)
      t -= 2 * M_PI;

    // Check if on Outer Boundary
    if (std::abs(r - dat.Rout) < tol) {
      // Check angles
      // Distance in angle
      // Lambda requires C++11. makefile standard?
      // Inline logic.
      auto check_port = [&](double center, int code) {
        double d = std::abs(t - center);
        if (d > M_PI)
          d = 2 * M_PI - d; // Wrap around
        if (d < w_ang) {
          dat.port_code[i] = code;
          return true;
        }
        return false;
      };

      if (check_port(ang1, 1))
        c1++;
      else if (check_port(ang2, 2))
        c2++;
      else if (check_port(ang3, 3))
        c3++;
      else if (check_port(ang4, 4))
        c4++;
    }
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "  Port Nodes identified: P1=%d, P2=%d, P3=%d, P4=%d\n",
                     c1, c2, c3, c4);
  CHKERRQ(ierr);

  return 0;
}

int main(int argc, char **argv) {
  PetscErrorCode ierr;
  int rank, size;
  std::time_t t0, t1;

  ierr = SlepcInitialize(&argc, &argv, nullptr, help);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "=== PETSc initialized ===\n");
  CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  CHKERRQ(ierr);

  if (argc < 2) {
    PetscPrintf(PETSC_COMM_WORLD,
                "Error: please provide FEMstruct input file\n");
    SlepcFinalize();
    return EXIT_FAILURE;
  }

  std::time(&t0);
  data dat;
  dat.rank = rank;
  dat.size = size;

  // Initialize pointers to null
  dat.node = nullptr;
  dat.elem = nullptr;
  dat.Convertnode = nullptr;
  dat.bnode_id = nullptr;
  dat.bnode = nullptr;
  dat.bvalue = nullptr;
  dat.port_code = nullptr;

  ierr = input_reader(argv[1], dat);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // EMR SWEEP LOGIC
  // ---------------------------------------------------------------------------
  if (dat.Nsample_R2 > 0) {
    if (rank == 0) {
      PetscPrintf(PETSC_COMM_WORLD,
                  "\n=== STARTING EMR SIMULATION (R2 & H Sweep) ===\n");
      // Initialize Output File
      FILE *fp = fopen("../output/EMRdata.out", "w");
      if (fp) {
        fprintf(fp, "# H(T)  Ro(Ohm)  R(H)(Ohm)  EMR(%%)\n");
        fclose(fp);
      }
    }

    // R2 Sweep Loop
    for (int i_R2 = 0; i_R2 < dat.Nsample_R2; ++i_R2) {
      double R2_curr = dat.R2min;
      if (dat.Nsample_R2 > 1) {
        R2_curr = dat.R2min + i_R2 * dat.dR2;
      }
      dat.R2 = R2_curr;

      PetscPrintf(PETSC_COMM_WORLD, "\n>> Step %d/%d: R2 = %.6e m\n", i_R2 + 1,
                  dat.Nsample_R2, dat.R2);

      // 1. Generate Mesh for current R2
      // Rank 0 generates, then all reload.
      if (dat.auto_generate_mesh) {
        if (rank == 0) {
          fprintf(stderr, "DEBUG: Calling generate_mesh...\n");
          // generate_mesh calls generate_emr_mesh which uses dat.R2
          ierr = generate_mesh(dat, dat.node_elem);
          CHKERRQ(ierr);
          fprintf(stderr, "DEBUG: generate_mesh done.\n");
        }
        MPI_Barrier(PETSC_COMM_WORLD);
      }

      // Load Mesh Data
      if (rank == 0)
        fprintf(stderr, "DEBUG: Calling mesh_input...\n");
      ierr = mesh_input(dat);
      CHKERRQ(ierr);
      if (rank == 0)
        fprintf(stderr, "DEBUG: mesh_input done.\n");

      // Initialize Ports
      ierr = identify_ports(dat);
      CHKERRQ(ierr);

      dat.nglobal = dat.ngnodes * dat.dof_per_node;

      // =========================================================================
      // PRE-CALCULATE REFERENCE RESISTANCE (R0) at H=0
      // =========================================================================
      // We calculate R0 once per geometry (R2) to ensure consistent EMR
      // baseline
      {
        if (rank == 0)
          PetscPrintf(PETSC_COMM_WORLD,
                      "   Computing Reference Resistance R0 at H=0.0...\n");

        double H_save = dat.H_current; // Just in case
        dat.H_current = 0.0;

        global_matrices gmat_ref = {NULL, NULL};

        if (rank == 0)
          fprintf(stderr, "DEBUG: Calling make_global...\n");
        ierr = make_global(gmat_ref, dat);
        CHKERRQ(ierr);

        if (rank == 0)
          fprintf(stderr, "DEBUG: Calling apply_bc...\n");
        ierr = apply_bc(gmat_ref, dat);
        CHKERRQ(ierr);

        if (rank == 0)
          fprintf(stderr, "DEBUG: Calling petsc_diagonalizer...\n");
        ierr = petsc_diagonalizer(gmat_ref, dat);
        CHKERRQ(ierr);

        if (rank == 0)
          fprintf(stderr, "DEBUG: R0 calculation done.\n");

        dat.Ro = dat.latest_resistance;

        if (rank == 0)
          PetscPrintf(PETSC_COMM_WORLD, "   >> R0 fixed at %.6e Ohm\n", dat.Ro);

        // Restore or Reset H seems unnecessary as loop overrides it, but good
        // practice
      }

      // H Sweep Loop
      for (int i_H = 0; i_H < dat.Nsample_H; ++i_H) {
        double H_curr = dat.Hmin;
        if (dat.Nsample_H > 1) {
          H_curr = dat.Hmin + i_H * dat.dH;
        }
        dat.H_current = H_curr;

        if (rank == 0 && i_H % 5 == 0) {
          PetscPrintf(PETSC_COMM_WORLD, "   .. Solving for H = %.4f T\n",
                      dat.H_current);
        }

        // Create fresh Matrix structures
        global_matrices gmat = {NULL, NULL};

        ierr = make_global(gmat, dat);
        CHKERRQ(ierr);
        ierr = apply_bc(gmat, dat);
        CHKERRQ(ierr);
        ierr = petsc_diagonalizer(gmat, dat);
        CHKERRQ(ierr);

        // gmat cleaned up by petsc_diagonalizer

        // =======================================================================
        // DATA LOGGING (Moved from petsc_diagonalizer to ensure consistency)
        // =======================================================================
        double R_H = dat.latest_resistance;
        double EMR = 0.0;
        if (dat.Ro != 0.0) {
          EMR = 100.0 * (R_H - dat.Ro) / dat.Ro;
        }

        if (rank == 0) {
          FILE *fp = fopen("../output/EMRdata.out", "a");
          if (fp) {
            // Format: H  R0  R(H)  EMR
            fprintf(fp, "%.4f  %.6e  %.6e  %.4f\n", dat.H_current, dat.Ro, R_H,
                    EMR);
            fclose(fp);
          }
          PetscPrintf(PETSC_COMM_WORLD,
                      "   >> H=%+.4f T: R=%.6e Ohm, EMR=%+.4f%%\n",
                      dat.H_current, R_H, EMR);
        }
      }

      // Cleanup Mesh before regenerating
      FreeMesh(dat);
    }

  } else {
    // -----------------------------------------------------------------------
    // STANDARD SIMULATION (Single Run)
    // -----------------------------------------------------------------------
    if (dat.auto_generate_mesh) {
      if (rank == 0) {
        ierr = generate_mesh(dat, dat.node_elem);
        CHKERRQ(ierr);
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }
    ierr = mesh_input(dat);
    CHKERRQ(ierr);
    dat.nglobal = dat.ngnodes * dat.dof_per_node;

    global_matrices gmat;
    ierr = make_global(gmat, dat);
    CHKERRQ(ierr);
    ierr = apply_bc(gmat, dat);
    CHKERRQ(ierr);
    ierr = petsc_diagonalizer(gmat, dat);
    CHKERRQ(ierr);

    FreeMesh(dat);
  }

  // Final Cleanup
  // Free other static arrays if any

  std::time(&t1);
  double minutes = std::difftime(t1, t0) / 60.0;
  PetscPrintf(PETSC_COMM_WORLD, "\nTotal execution time: %.2f mins\n", minutes);

  SlepcFinalize();
  return EXIT_SUCCESS;
}