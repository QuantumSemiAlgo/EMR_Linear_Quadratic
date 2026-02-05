# EMR Simulation: Linear & Quadratic Elements

[![Docker](https://img.shields.io/badge/Docker-ready-blue.svg)](Dockerfile)
[![C++14](https://img.shields.io/badge/C%2B%2B-14-blue.svg)](src/)
[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](requirements.txt)

Finite Element Method (FEM) simulation of **Extraordinary Magnetoresistance (EMR)** in semiconductor-metal hybrid structures. Solves a modified Laplace equation with Hall effect terms on a 2D annular domain.

## Features

- **Physics**: Solves ∇·(σ̂∇φ) = 0 where σ̂ is the conductivity tensor with magnetic field dependence
- **Geometry**: Annular semiconductor disk with concentric metallic inclusion
- **Elements**: Supports Linear (Q1) and Quadratic (Q2) Lagrange elements
- **Solver**: PETSc/SLEPc sparse matrix solver with KSP linear solver
- **Visualization**: Python scripts for plotting resistance, EMR vs magnetic field

## Quick Start

### Option 1: Docker (Recommended)

No dependencies needed! Just Docker.

```bash
# Build the container
docker build -t emr-sim .

# Run the simulation
docker run --rm -v $(pwd)/output:/app/output emr-sim

# Results will appear in ./output/
```

### Option 2: Local Build

#### Prerequisites

1. **C++ Dependencies:**
   - PETSc (3.x+) - Portable Extensible Toolkit for Scientific Computation
   - SLEPc (3.x+) - Scalable Library for Eigenvalue Problem Computations
   - MPI implementation (OpenMPI or MPICH)
   - C++14 compatible compiler

2. **Python Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

#### Installation

**On Ubuntu/Debian:**
```bash
# Install PETSc and SLEPc
sudo apt-get update
sudo apt-get install petsc-dev slepc-dev libopenmpi-dev

# Set environment variables (add to ~/.bashrc)
export PETSC_DIR=/usr/lib/petsc
export SLEPC_DIR=/usr/lib/slepc
```

**On macOS:**
```bash
# Using Homebrew
brew install petsc slepc open-mpi

# Set environment variables
export PETSC_DIR=/opt/homebrew/opt/petsc
export SLEPC_DIR=/opt/homebrew/opt/slepc
```

**Manual Installation:**
If not using package managers, set `PETSC_DIR` and `SLEPC_DIR` to your installation paths.

#### Build & Run

```bash
# Navigate to source directory
cd src

# Build the solver
make

# Run a single simulation
./Laplace ../config/FEMstruct2d.inp

# Or run full magnetic field sweep
./run_emr_sweep.sh
```

## Project Structure

```
EMR_Linear_Quadratic/
├── config/                    # Configuration files
│   ├── FEMstruct2d.inp       # Main input file (annular geometry)
│   ├── FEMstruct_hermite.inp # Hermite element config
│   └── FEMstruct_test.inp    # Test configuration
├── src/                       # C++ source code
│   ├── main_1.cpp            # Main entry point
│   ├── mesh_generator_1.cpp  # Mesh generation
│   ├── petsc_diagonalizer_1.cpp # PETSc solver interface
│   ├── apply_bc_1.cpp        # Boundary conditions
│   ├── Makefile              # Build configuration
│   ├── run_emr_sweep.sh      # Automated sweep script
│   └── *.py                  # Python analysis scripts
├── mesh/                      # Mesh data files
├── output/                    # Simulation outputs
├── requirements.txt           # Python dependencies
├── Dockerfile                 # Docker container definition
└── README.md                  # This file
```

## Usage

### Configuration

Edit `config/FEMstruct2d.inp` to modify:
- Geometry (inner/outer radii R1, R2)
- Material properties (conductivity, mobility)
- Magnetic field sweep range (H_min, H_max)
- Mesh density (number of elements)
- Element type (Linear or Quadratic)

### Running Simulations

**Single run:**
```bash
cd src
./Laplace ../config/FEMstruct2d.inp
```

**Magnetic field sweep:**
```bash
cd src
./run_emr_sweep.sh
# Outputs: ../output/EMRdata.out and plots
```

**Diagnostic mode:**
```bash
cd src
./run_emr_diagnostics.sh
```

### Visualization

```bash
# Plot EMR vs H
python3 src/plot_emr_vs_h.py output/EMRdata.out

# Visualize mesh
python3 src/visualize_mesh.py

# Verify solution
python3 src/verify_solution.py
```

## Output Files

After running, check `output/` directory:
- `EMRdata.out` - Resistance and EMR data vs magnetic field
- `EMR_Plot_EMR_vs_H.png` - EMR vs H plot
- `EMR_Plot_R_vs_H.png` - Resistance vs H plot
- `emr_sweep.log` - Detailed simulation log

## Troubleshooting

**"PETSC_NOT_FOUND" error:**
```bash
# Set environment variables
export PETSC_DIR=/path/to/petsc
export SLEPC_DIR=/path/to/slepc
```

**MPI errors:**
```bash
# Use fewer processes or run serially
mpiexec -n 1 ./Laplace ../config/FEMstruct2d.inp
```

**Python import errors:**
```bash
pip install -r requirements.txt
```

## Technical Details

### Physics Model

The simulation solves:
```
∇·(σ̂∇φ) = 0
```

where the conductivity tensor with Hall effect is:
```
σ̂ = σ / (1 + (μB)²) * [[1, -μB], [μB, 1]]
```

- σ: scalar conductivity
- μ: carrier mobility
- B: magnetic field (perpendicular to plane)

### Boundary Conditions

Four contact ports on outer boundary:
1. Current input (I_in)
2. Current output (I_out)
3. Voltage reference (V_ref)
4. Voltage measurement (V_meas)

Periodic BCs on angular boundaries.

### EMR Calculation

```
EMR(H) = [R(H) - R(0)] / R(0) × 100%
```

where R = V_meas / I_in

## References

For more detailed technical documentation, see:
- `CODE_EXPLANATION.md` - Architecture and code walkthrough
- `README_DOCKER.md` - Docker-specific instructions

## License

[Add your license here]

## Contact

[Add contact information]
