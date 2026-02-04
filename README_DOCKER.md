# EMR Simulation Docker Container

## Overview
This container runs the standard Extraordinary Magnetoresistance (EMR) simulation pipeline, including:
1.  **FEM Solver**: C++ Finite Element Method solver using PETSc/SLEPc.
2.  **Simulation Sweep**: Shell script (`run_emr_sweep.sh`) orchestrating magnetic field (H) sweeps.
3.  **Visualization**: Python script (`plot_emr_vs_h.py`) generating R vs H and EMR vs H plots.

## Prerequisites
-   Docker installed on your system.

## Usage

### 1. Build the Image
Navigate to the project root (where `Dockerfile` is located) and run:
```bash
docker build -t emr-sim .
```

### 2. Run the Simulation
To run the simulation and retrieve the output plots, mount the `output` directory:
```bash
docker run --rm -v $(pwd)/output:/app/output emr-sim
```
*   The simulation will run inside the container.
*   Results (logs, data files, and PNG plots) will be saved to your local `output/` directory.

### 3. Interactive Mode
To run manually or debug inside the container:
```bash
docker run --rm -it -v $(pwd)/output:/app/output emr-sim /bin/bash
```

## Structure
-   `/app/src`: Source code (compiled inside container).
-   `/app/output`: Output directory (should be mounted).
