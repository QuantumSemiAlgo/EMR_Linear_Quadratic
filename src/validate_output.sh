#!/bin/bash
# Validates EMR simulation output
# Usage: ./validate_output.sh [output_dir]
#   Default output_dir: ../output

OUTPUT_DIR="${1:-../output}"
ERRORS=0
WARNINGS=0

echo "=========================================="
echo "  OUTPUT VALIDATION"
echo "=========================================="
echo "Checking: $OUTPUT_DIR"
echo ""

# Check 1: EMRdata.out exists
if [ ! -f "$OUTPUT_DIR/EMRdata.out" ]; then
    echo "✗ FAIL: EMRdata.out not found in $OUTPUT_DIR"
    ERRORS=$((ERRORS + 1))
else
    echo "✓ PASS: EMRdata.out exists"

    # Check 2: File is not empty
    if [ ! -s "$OUTPUT_DIR/EMRdata.out" ]; then
        echo "✗ FAIL: EMRdata.out is empty"
        ERRORS=$((ERRORS + 1))
    else
        echo "✓ PASS: EMRdata.out is not empty"

        # Check 3: Has at least 2 lines (header + 1 data point)
        LINE_COUNT=$(wc -l < "$OUTPUT_DIR/EMRdata.out")
        if [ "$LINE_COUNT" -lt 2 ]; then
            echo "✗ FAIL: EMRdata.out has only $LINE_COUNT lines (need at least 2)"
            ERRORS=$((ERRORS + 1))
        else
            echo "✓ PASS: EMRdata.out has $LINE_COUNT lines"
        fi

        # Check 4: Contains expected header
        if ! head -1 "$OUTPUT_DIR/EMRdata.out" | grep -q "H(T)"; then
            echo "⚠ WARN: EMRdata.out header doesn't match expected format"
            WARNINGS=$((WARNINGS + 1))
        else
            echo "✓ PASS: EMRdata.out header format correct"
        fi
    fi
fi

# Check 5: Plot files (warnings only, not failures)
PLOT_FILES=(
    "EMR_Plot_EMR_vs_H.png"
    "EMR_Plot_R_vs_H.png"
)

for PLOT in "${PLOT_FILES[@]}"; do
    if [ -f "$OUTPUT_DIR/$PLOT" ]; then
        echo "✓ PASS: $PLOT exists"
    else
        echo "⚠ WARN: $PLOT not found (may not have been generated)"
        WARNINGS=$((WARNINGS + 1))
    fi
done

echo ""
echo "=========================================="
if [ $ERRORS -eq 0 ]; then
    echo "✓ OK: Output validated in $OUTPUT_DIR"
    echo "   Errors: $ERRORS, Warnings: $WARNINGS"
    exit 0
else
    echo "✗ FAILED: Output validation failed"
    echo "   Errors: $ERRORS, Warnings: $WARNINGS"
    exit 1
fi
