#!/bin/bash
# EMR Diagnostic Runner
# Compiles code with diagnostics and runs

echo "════════════════════════════════════════════════════════════════"
echo "           EMR Diagnostic Mode                                  "
echo "════════════════════════════════════════════════════════════════"

# Compile
echo "Compiling with diagnostic output..."
make clean
make

if [ $? -ne 0 ]; then
    echo "❌ Compilation failed!"
    exit 1
fi

# Run simulation
echo ""
echo "Running simulation..."
echo ""

# Run and capture output
mpiexec -n 1 ./Laplace ../config/FEMstruct2d.inp 2>&1 | tee emr_diagnostic_output.txt

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "Diagnostic output saved to: emr_diagnostic_output.txt"
echo "EMR data saved to: ../output/EMRdata.out"
echo ""
echo "QUICK CHECKS:"
echo "════════════════════════════════════════════════════════════════"

# Extract key values
if [ -f "EMRdata.out" ]; then
    echo "Latest EMR data:"
    tail -5 EMRdata.out
    echo ""
    
    # Calculate EMR from last line
    last_line=$(tail -1 EMRdata.out)
    H=$(echo $last_line | awk '{print $1}')
    R0=$(echo $last_line | awk '{print $2}') # Wait, column 2 is R(H), where is R0? 
    # Usually EMRdata.out format is H, R(H), R2, Sig1, Sig2.
    # We don't have R0 in the file.
    RH=$(echo $last_line | awk '{print $2}')
    
    echo "At H = $H T:"
    echo "  R(H) = $RH Ω"
    echo ""
    
    # Note: R0 is not column 2. R(H) is column 2.
    # The shell script logic in prompt assumed format H R0 RH EMR ?
    # But petsc_diagonalizer outputs: fprintf(fp, "%.6e  %.8e  %.6e  %.6e  %.6e\n", dat.H_current, Resistance, dat.R2, dat.sigma1, dat.sigma2);
    # So Columns: 1:H, 2:R, 3:R2, 4:Sig1, 5:Sig2.
    
    echo "Use analyze_emr_results.py for full analysis."
fi

echo "════════════════════════════════════════════════════════════════"
