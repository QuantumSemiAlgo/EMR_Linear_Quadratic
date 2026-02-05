#!/bin/bash
# Fast smoke test for EMR simulation
# Uses reduced mesh and minimal H sweep

set -e  # Exit on any error

# Create smoke output directory
mkdir -p ../output/smoke

# Clean previous smoke test output
rm -f ../output/smoke/EMRdata.out
rm -f EMRdata_smoke.out

echo "=========================================="
echo "  EMR SMOKE TEST (Fast Validation)"
echo "=========================================="
echo "Using config: ../config/FEMstruct2d_smoke.inp"
echo "Output dir: ../output/smoke/"
echo ""

# Run FEM Solver with smoke config
./Laplace ../config/FEMstruct2d_smoke.inp -ksp_rtol 1e-4 > smoke.log 2>&1

# Move output to smoke directory
if [ -f "../output/EMRdata.out" ]; then
    mv ../output/EMRdata.out ../output/smoke/EMRdata.out
    echo "✓ Simulation completed"
else
    echo "✗ ERROR: No output file generated"
    cat smoke.log
    exit 1
fi

echo ""
echo "Smoke test simulation complete."
echo ""

# Validate output
if [ -f "./validate_output.sh" ]; then
    echo "Running output validation..."
    ./validate_output.sh "../output/smoke"
else
    echo "Note: validate_output.sh not found, skipping validation"
fi
