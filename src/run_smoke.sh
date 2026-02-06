#!/bin/bash
# Fast smoke test for EMR simulation
# Uses reduced mesh and minimal H sweep
# Usage: ./run_smoke.sh [q1|q2|hermite_cubic|hermite_quintic]

set -e  # Exit on any error

# Determine which config to use
ELEMENT_TYPE="${1:-q2}"
case "$ELEMENT_TYPE" in
    q1)
        CONFIG="../config/FEMstruct2d_smoke.inp"
        LABEL="Q1 Lagrangian"
        OUTPUT_SUFFIX="lagrangian_linear"
        ;;
    q2)
        CONFIG="../config/FEMstruct2d_smoke.inp"
        LABEL="Q2 Lagrangian"
        OUTPUT_SUFFIX="lagrangian_quadratic"
        ;;
    hermite_cubic)
        CONFIG="../config/FEMstruct2d_hermite_cubic_smoke.inp"
        LABEL="Bicubic Hermite (4 nodes, 4 DOF/node)"
        OUTPUT_SUFFIX="hermite_cubic"
        ;;
    hermite_quintic)
        CONFIG="../config/FEMstruct2d_hermite_quintic_smoke.inp"
        LABEL="Biquintic Hermite (9 nodes, 4 DOF/node)"
        OUTPUT_SUFFIX="hermite_quintic"
        ;;
    *)
        echo "Unknown element type: $ELEMENT_TYPE"
        echo "Usage: $0 [q1|q2|hermite_cubic|hermite_quintic]"
        exit 1
        ;;
esac

# Create smoke output directory
mkdir -p ../output/smoke

# Clean previous smoke test output
rm -f ../output/smoke/EMRdata*.out
rm -f ../output/EMRdata*.out
rm -f EMRdata_smoke.out

echo "=========================================="
echo "  EMR SMOKE TEST (Fast Validation)"
echo "=========================================="
echo "Element Type: $LABEL"
echo "Using config: $CONFIG"
echo "Output dir: ../output/smoke/"
echo ""

# Run FEM Solver with smoke config
./Laplace "$CONFIG" -ksp_rtol 1e-4 > ../output/smoke/smoke.log 2>&1

# Move output to smoke directory
TYPE_SPECIFIC_FILE="../output/EMRdata_${OUTPUT_SUFFIX}.out"
if [ -f "$TYPE_SPECIFIC_FILE" ]; then
    # Move type-specific file
    mv "$TYPE_SPECIFIC_FILE" "../output/smoke/EMRdata_${OUTPUT_SUFFIX}.out"
    # Also move generic file for backwards compatibility
    if [ -f "../output/EMRdata.out" ]; then
        mv ../output/EMRdata.out ../output/smoke/EMRdata.out
    fi
    echo "✓ Simulation completed"
    echo "  Output: ../output/smoke/EMRdata_${OUTPUT_SUFFIX}.out"
else
    echo "✗ ERROR: No output file generated (expected: $TYPE_SPECIFIC_FILE)"
    cat ../output/smoke/smoke.log
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
