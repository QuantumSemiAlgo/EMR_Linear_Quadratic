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
static char help[] = "Parallel FEM for 2D Laplace on rectangular mesh\n\n";

/**
 * @brief Main function of the FEM solver.
 *
 * @param argc Argument count.
 * @param argv Argument vector. argv[1] is expected to be the input file path.
 * @return int Exit status (EXIT_SUCCESS or EXIT_FAILURE).
 */
int main(int argc, char **argv) {
  // ---------------------------------------------------------------------------
  // Variable Declarations
  // ---------------------------------------------------------------------------
  // ierr: Used to capture and check the return codes of PETSc functions.
  PetscErrorCode ierr;

  // rank: The ID of the current MPI process (0 to size-1).
  // size: The total number of MPI processes in the communicator.
  int rank, size;

  // Timers for measuring total execution time.
  std::time_t t0, t1;

  // ---------------------------------------------------------------------------
  // 1. Initialize PETSc/SLEPc Environment
  // ---------------------------------------------------------------------------
  // SlepcInitialize handles the initialization of both SLEPc and PETSc.
  // It sets up memory management, parses command-line arguments, and
  // initializes MPI.
  ierr = SlepcInitialize(&argc, &argv, nullptr, help);
  CHKERRQ(ierr);

  // Print a welcome message to the standard output (only from the root
  // process).
  ierr = PetscPrintf(PETSC_COMM_WORLD, "=== PETSc initialized ===\n");
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 2. MPI Setup
  // ---------------------------------------------------------------------------
  // Determine the rank of the current process and the total number of
  // processes.
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 3. Input Validation
  // ---------------------------------------------------------------------------
  // Ensure that the user has provided an input file argument.
  if (argc < 2) {
    ierr = PetscPrintf(
        PETSC_COMM_WORLD,
        "Error: please provide FEMstruct input file on command line\n");
    CHKERRQ(ierr);
    // Finalize PETSc before exiting to ensure proper cleanup.
    PetscFinalize();
    return EXIT_FAILURE;
  }

  // ---------------------------------------------------------------------------
  // 4. Print Execution Banner
  // ---------------------------------------------------------------------------
  // Display current time, date, and user information for logging purposes.
  {
    std::time_t now = std::time(nullptr);
    std::tm *ltm = std::localtime(&now);

    char timebuf[32], datebuf[32];
    std::strftime(timebuf, sizeof(timebuf), "%H:%M:%S", ltm);
    std::strftime(datebuf, sizeof(datebuf), "%Y-%m-%d", ltm);

    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "\n"
                       "===============================\n"
                       ":      FEM calculation         :\n"
                       ":     Laplace Equation         :\n"
                       ":     Arya Retheeshan          :\n"
                       ":                              :\n"
                       ":   Time: %-10s           :\n"
                       ":   Date: %-10s           :\n"
                       "===============================\n",
                       timebuf, datebuf);
    CHKERRQ(ierr);
  }

  // Store the input filename.
  const char *fname = argv[1];
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Input file: %s\n", fname);
  CHKERRQ(ierr);

  // Start the execution timer.
  std::time(&t0);

  // ---------------------------------------------------------------------------
  // 5. Read Configuration Input
  // ---------------------------------------------------------------------------
  // Initialize the main data structure 'dat' which holds all simulation
  // parameters.
  data dat;
  dat.rank = rank;
  dat.size = size;

  // Parse the input file to populate 'dat'.
  ierr = input_reader(fname, dat);
  CHKERRQ(ierr);

// Create output directory if it doesn't exist
#ifdef _WIN32
  _mkdir(dat.outpath);
#else
  mkdir(dat.outpath, 0755);
#endif

  // LOG CONFIGURATION SUMMARY
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n=== Configuration Summary ===\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Mesh generation:\n"
                     "  n_x = %d elements\n"
                     "  n_y = %d elements\n"
                     "  auto_generate = %d\n"
                     "  use_annular = %d\n",
                     dat.n_x, dat.n_y, dat.auto_generate_mesh, dat.use_annular);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Output grid:\n"
                     "  ndz_x = %d points\n"
                     "  ndz_y = %d points\n",
                     dat.ndz_x, dat.ndz_y);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Element type:\n"
                     "  node_elem = %d\n"
                     "  dof_per_node = %d\n",
                     dat.node_elem, dat.dof_per_node);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "=============================\n\n");
  CHKERRQ(ierr);

  // Log key configuration parameters.
  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "FEMstruct read: ndz_x=%d, ndz_y=%d, node_elem=%d, dof_per_node=%d\n",
      dat.ndz_x, dat.ndz_y, dat.node_elem, dat.dof_per_node);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 6. Compute Grid Spacings
  // ---------------------------------------------------------------------------
  // Calculate the spatial step sizes (dx, dy) based on domain size and grid
  // resolution.
  dat.dzx = (dat.xmax - dat.xmin) / double(dat.ndz_x);
  dat.dzy = (dat.ymax - dat.ymin) / double(dat.ndz_y);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Grid spacings: dzx=%.6f, dzy=%.6f\n",
                     dat.dzx, dat.dzy);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 7. Generate Mesh Files (if requested)
  // ---------------------------------------------------------------------------
  if (dat.auto_generate_mesh) {
    // CRITICAL STABILITY FIX:
    // Only Rank 0 generates the mesh files to prevent race conditions and file
    // corruption. All other ranks wait at the barrier.
    if (dat.rank == 0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "Auto-generating mesh files (Rank 0)...\n");
      CHKERRQ(ierr);

      // Generate for 4-node elements
      ierr = generate_mesh(dat, 4);
      CHKERRQ(ierr);

      // Generate for 9-node elements
      ierr = generate_mesh(dat, 9);
      CHKERRQ(ierr);

      ierr = PetscPrintf(PETSC_COMM_WORLD, "Mesh generation complete.\n\n");
      CHKERRQ(ierr);
    }

    // Ensure all ranks wait until mesh files are fully written and closed
    ierr = MPI_Barrier(PETSC_COMM_WORLD);
    CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // 7. Read Mesh Data
  // ---------------------------------------------------------------------------
  // Load nodes, elements, and boundary definitions from external mesh files.
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Reading mesh files...\n");
  CHKERRQ(ierr);
  ierr = mesh_input(dat);
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Mesh input complete: ngnodes=%d, nelem=%d\n", dat.ngnodes,
                     dat.nelem);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 8. Compute Global Degrees of Freedom (DOFs)
  // ---------------------------------------------------------------------------
  // Total DOFs = Number of Nodes * DOFs per Node.
  dat.nglobal = dat.ngnodes * dat.dof_per_node;

  ierr =
      PetscPrintf(PETSC_COMM_WORLD, "Global DOFs: nglobal=%d\n", dat.nglobal);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Optional: Hermite Shape Function Testing
  // ---------------------------------------------------------------------------
  // This block runs only for specific element types (Hermite Cubic) to verify
  // shape function derivatives. It serves as a sanity check for the basis
  // functions.
  if (dat.node_elem == 4 && dat.dof_per_node == 4 && dat.rank == 0) {
    PetscPrintf(
        PETSC_COMM_WORLD,
        "── Testing Hermite nodal-value partition & derivative sums ──\n");

    double N[16], dNdxi[16], dNdeta[16];
    double xi_vals[5] = {-1.0, -0.5, 0.0, 0.5, 1.0};
    double eta_vals[5] = {-1.0, -0.5, 0.0, 0.5, 1.0};

    // Iterate over test points in the reference element (xi, eta).
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
        // Evaluate shape functions and derivatives.
        shapeHermiteCubic(N, xi_vals[i], eta_vals[j]);
        deriv1HermiteCubic(dNdxi, dNdeta, xi_vals[i], eta_vals[j]);

        // Sum the shape functions associated with nodal values (u).
        // For a partition of unity, these should sum to 1 (if only u-DOFs are
        // considered).
        double sumN = N[0] + N[4] + N[8] + N[12];
        double sumXi = dNdxi[0] + dNdxi[4] + dNdxi[8] + dNdxi[12];
        double sumEta = dNdeta[0] + dNdeta[4] + dNdeta[8] + dNdeta[12];

        PetscPrintf(PETSC_COMM_WORLD,
                    " at (ξ,η)=(% .1f,% .1f): ΣN^u=% .6f, Σ∂N^u/∂ξ=% .6f, "
                    "Σ∂N^u/∂η=% .6f\n",
                    xi_vals[i], eta_vals[j], sumN, sumXi, sumEta);
      }
    }
    PetscPrintf(PETSC_COMM_WORLD,
                "──────────────────────────────────────────────────\n");
  }

  // ---------------------------------------------------------------------------
  // 9. Assemble Global System
  // ---------------------------------------------------------------------------
  // Construct the global stiffness matrix and load vector.
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Assembling global system...\n");
  CHKERRQ(ierr);

  global_matrices gmat;
  ierr = make_global(gmat, dat);
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Assembly complete\n");
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 10. Apply Boundary Conditions
  // ---------------------------------------------------------------------------
  // Apply Dirichlet and Neumann boundary conditions to the system.
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Applying boundary conditions...\n");
  CHKERRQ(ierr);
  ierr = apply_bc(gmat, dat);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "BCs applied\n");
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 11. Solve System
  // ---------------------------------------------------------------------------
  // Solve the linear system Ax = b using the configured PETSc solver.
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Solving system...\n");
  CHKERRQ(ierr);
  ierr = petsc_diagonalizer(gmat, dat);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Solve + postprocess done\n");
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 12. Finalize and Report
  // ---------------------------------------------------------------------------
  // Stop the timer and report total execution time.
  std::time(&t1);
  double minutes = std::difftime(t1, t0) / 60.0;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "===========================================\n"
                     "Total time: %.2f mins\n"
                     "===========================================\n",
                     minutes);
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // 13. Cleanup Dynamic Memory
  // ---------------------------------------------------------------------------
  // Free dynamically allocated arrays to prevent memory leaks
  if (dat.Convertnode != nullptr) {
    delete[] dat.Convertnode;
    dat.Convertnode = nullptr;
  }

  if (dat.bnode_id != nullptr) {
    delete[] dat.bnode_id;
    dat.bnode_id = nullptr;
  }

  if (dat.node != nullptr) {
    for (int i = 0; i < dat.ngnodes; ++i) {
      delete[] dat.node[i];
    }
    delete[] dat.node;
    dat.node = nullptr;
  }

  if (dat.elem != nullptr) {
    for (int i = 0; i < dat.nelem; ++i) {
      delete[] dat.elem[i];
    }
    delete[] dat.elem;
    dat.elem = nullptr;
  }

  if (dat.bnode != nullptr) {
    for (int i = 0; i < dat.nbnode; ++i) {
      delete[] dat.bnode[i];
    }
    delete[] dat.bnode;
    dat.bnode = nullptr;
  }

  if (dat.bvalue != nullptr) {
    for (int i = 0; i < dat.nbnode; ++i) {
      delete[] dat.bvalue[i];
    }
    delete[] dat.bvalue;
    dat.bvalue = nullptr;
  }

  if (dat.btob != nullptr) {
    delete[] dat.btob;
    dat.btob = nullptr;
  }

  if (dat.material != nullptr) {
    delete[] dat.material;
    dat.material = nullptr;
  }

  if (dat.layer_id_for_elem != nullptr) {
    delete[] dat.layer_id_for_elem;
    dat.layer_id_for_elem = nullptr;
  }

  if (dat.layer_id != nullptr) {
    delete[] dat.layer_id;
    dat.layer_id = nullptr;
  }

  if (dat.layer_ymin != nullptr) {
    delete[] dat.layer_ymin;
    dat.layer_ymin = nullptr;
  }

  if (dat.layer_ymax != nullptr) {
    delete[] dat.layer_ymax;
    dat.layer_ymax = nullptr;
  }

  if (dat.layer_material != nullptr) {
    delete[] dat.layer_material;
    dat.layer_material = nullptr;
  }

  if (dat.xigaus != nullptr) {
    delete[] dat.xigaus;
    dat.xigaus = nullptr;
  }

  if (dat.etagaus != nullptr) {
    delete[] dat.etagaus;
    dat.etagaus = nullptr;
  }

  if (dat.wgaus != nullptr) {
    delete[] dat.wgaus;
    dat.wgaus = nullptr;
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Memory cleanup complete\n");
  CHKERRQ(ierr);

  // Clean up PETSc/SLEPc internal structures.
  ierr = SlepcFinalize();
  CHKERRQ(ierr);

  return EXIT_SUCCESS;
}