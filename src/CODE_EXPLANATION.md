# EMR Simulation Codebase Explanation

## Overview
This project simulates **Extraordinary Magnetoresistance (EMR)** in a semiconductor-metal hybrid structure using the **Finite Element Method (FEM)**. It solves a modified Laplace equation with Hall effect terms on a 2D annulus domain.

**Key Features:**
*   **Physics**: Solves $\nabla \cdot (\hat{\sigma} \nabla \phi) = 0$ where $\hat{\sigma}$ is the conductivity tensor including magnetic field ($B$) dependence.
*   **Geometry**: Annular semiconductor disk with a concentric metallic inclusion (shunts current).
*   **Numerics**: Uses PETSc/SLEPc for sparse matrix assembly and linear solving (KSP).
*   **Elements**: Supports Linear (Q1) and Quadratic (Q2) Lagrange elements.

---

## 1. Core Architecture

### `global_params_1.h`
*   **Purpose**: Defines the global `data` structure passed to all functions.
*   **Key Fields**:
    *   `R1`, `R2`: Geometry radii.
    *   `sigma1`, `sigma2`, `mu1`, `mu2`: Material properties (Conductivity, Mobility).
    *   `H_current`: Current magnetic field value.
    *   `Convertnode`: Mapping array for handling periodic boundary conditions.
    *   `elem`, `coor`: Mesh connectivity and node coordinates.

### `main_1.cpp`
*   **Purpose**: The entry point and orchestration loop.
*   **Workflow**:
    1.  Initializes MPI and PETSc.
    2.  Reads configuration via `read_input()`.
    3.  **Outer Loop**: Sweeps Inner Radius `R2` (if configured).
    4.  **Inner Loop**: Sweeps Magnetic Field `H` (from `-Hmax` to `+Hmax`).
    5.  For each step:
        *   Regenerates mesh (if geometry changed).
        *   Calls `Solve_system()` (assembled in `petsc_diagonalizer`).
        *   Calculates Resistance and EMR.
        *   Writes output to `EMRdata.out`.

---

## 2. Mesh Generation

### `mesh_generator_1.cpp`
*   **Function**: `generate_emr_mesh()`
*   **Logic**: Creates a structured grid mapped to an annulus.
    *   **Radial Zoning**: Divides radius into 5 zones to concentrate element density near the Semiconductor-Metal interface ($R_2$) where physics is most active.
    *   **Angular Zoning**: Ensures nodes align exactly with the 4 contact ports (Current I/O, Voltage Sense).
    *   **Port Codes**: Tags nodes on the boundary with `port_code` (1=I_in, 2=I_out, 3=V_ref, 4=V_meas) for BC application.

### `convertnode_1.cpp`
*   **Purpose**: Handles **Periodic Boundary Conditions**.
*   **Logic**: Identifies nodes at $\theta=0$ and $\theta=2\pi$. Maps the "overlap" nodes (at $2\pi$) back to their canonical counterparts (at $0$), ensuring continuity of the solution around the ring.

---

## 3. Physics Implementation

### `make_global_1.cpp` (The Engine)
*   **Purpose**: Assembles the Global Stiffness Matrix ($A$) and Load Vector ($b$).
*   **Key Logic**:
    *   Loops over all elements.
    *   Computes **Material Properties** at each quadrature point using `fermiF.cpp` (smooth transition between semiconductor and metal).
    *   Calculates the **Conductivity Tensor** $\sigma$ with Hall terms:
        $$ \sigma_{rr} = \sigma / (1+\beta^2) $$
        $$ \sigma_{r\theta} = -\sigma \beta / (1+\beta^2) $$
    *   **Weak Form**: Assembles $\int (\nabla v)^T \cdot \sigma \cdot (\nabla u) d\Omega$.
    *   **Hall Effect**: The off-diagonal term $\sigma_{r\theta}$ couples the radial and tangential fields, creating the Lorentz force effect.

### `fermiF.cpp`
*   **Purpose**: Provides a smooth step function (Fermi-Dirac like).
*   **Usage**: Transitions conductivity $\sigma(r)$ and mobility $\mu(r)$ sharpness from Metal to Semiconductor to avoid numerical instability at the sharp interface.

---

## 4. Boundary Conditions

### `apply_bc_1.cpp`
*   **Purpose**: Enforces physical constraints on the system.
*   **Applied BCs**:
    1.  **Periodic BCs**: Enforced via matrix row elimination (linking $\theta=0$ to $\theta=2\pi$).
    2.  **Reference Voltage**: Sets Potential $V=0$ at Port 3 (`MatZeroRows`).
    3.  **Current Injection**: Applies Neumann BC (current flux) at Ports 1 (Source) and 2 (Sink) by integrating current density $J$ into the RHS vector $b$.
    4.  **Port Constraints**: Constrains tangential derivatives at metal ports (making them equipotential surfaces).

---

## 5. Solver & Post-Processing

### `petsc_diagonalizer_1.cpp`
*   **Purpose**: Solves the linear system $Ax=b$.
*   **Logic**:
    *   Uses **KSP** (Krylov Subspace Methods) typically GMRES or MUMPS (direct solver).
    *   Extracts solution potential $\phi$ at Port 4.
    *   Computes **Resistance**: $R = (V_4 - V_3) / I_{source}$.
    *   Computes **EMR Ratio**: $EMR = (R(H) - R(0)) / R(0) \times 100\%$.

### `analyze_emr_results.py` & `plot_emr_vs_h.py`
*   **Purpose**: Visualization.
*   **Output**: Reads `EMRdata.out` and generates plots:
    *   `EMR_Plot_R_vs_H.png`: Resistance vs Magnetic Field (Parabolic).
    *   `EMR_Plot_EMR_vs_H.png`: EMR % vs Magnetic Field.

---

## 6. How to Run

### Manual Run
```bash
cd src
make clean && make
./run_emr_sweep.sh
python3 plot_emr_vs_h.py
```

### Docker Run
```bash
docker build -t emr-sim .
docker run --rm -v $(pwd)/output:/app/output emr-sim
```

---

## 7. Configuration (`FEMstruct2d.inp`)
*   **Element Type**:
    *   `4 1 0` = Linear Lagrange (Fast, Robust).
    *   `9 1 0` = Quadratic Lagrange (Higher Precision).
*   **Geometry**: $R_1$ (Outer), $R_2$ (Inner Metal).
*   **Physics**: $\mu$ (Mobility) - Key driver of EMR magnitude.
*   **Sweep**: $H_{min}, H_{max}$ - Range of magnetic field.
