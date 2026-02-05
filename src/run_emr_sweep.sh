#!/bin/bash

# Clean previous output
rm -f ../output/EMRdata.out
rm -f EMRdata.out

# Run FEM Solver
echo "Starting EMR Mesh and Simulation Sweep..."
# Ensure we use the correct input file
# Binary is named Laplace
./Laplace ../config/FEMstruct2d.inp -ksp_rtol 1e-4 2>&1 | tee emr_sweep.log

# Check if simulation finished correctly
if [ $? -ne 0 ]; then
    echo "Error: Simulation failed!"
    exit 1
fi

echo "Simulation complete. Generating plot..."

# Run Plotter
python3 plot_emr_vs_h.py ../output/EMRdata.out

if [ -f "../output/EMR_Plot_EMR_vs_H.png" ]; then
    echo "Success! Plots generated in ../output/"
else
    echo "Error: Plots not generated."
fi
