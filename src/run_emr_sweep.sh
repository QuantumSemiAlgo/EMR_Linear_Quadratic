#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

OUTDIR="../output"
mkdir -p "$OUTDIR"

LOG="$OUTDIR/emr_sweep.log"

cleanup() {
  echo "Cancelled."
  # kill whole process group
  kill -- -$$ 2>/dev/null || true
  exit 130
}
trap cleanup INT TERM
set -m

# Clean previous output
rm -f ../output/EMRdata.out
rm -f EMRdata.out

# Run FEM Solver
echo "Starting EMR Mesh and Simulation Sweep..."
# Ensure we use the correct input file
# Binary is named Laplace
./Laplace ../config/FEMstruct2d.inp -ksp_rtol 1e-4 2>&1 | tee -a "$LOG"

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
