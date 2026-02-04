#!/usr/bin/env python3
"""
Analyze EMR Results and Compare to Expected Values
"""

import numpy as np
import sys
import os

def analyze_emr_data(filename='EMRdata.out'):
    """Parse EMRdata.out and check for issues"""
    
    print("=" * 70)
    print("                    EMR ANALYSIS REPORT")
    print("=" * 70)
    
    if not os.path.exists(filename):
        # Try checking ../output/EMRdata.out
        if os.path.exists('../output/EMRdata.out'):
            filename = '../output/EMRdata.out'
        elif os.path.exists('src/EMRdata.out'):
             filename = 'src/EMRdata.out'

    print(f"Reading file: {filename}")

    try:
        # Check if file is empty
        if os.path.getsize(filename) == 0:
            print("❌ File is empty")
            return

        data = np.loadtxt(filename, comments='#', skiprows=1)
    except Exception as e:
        print(f"❌ Could not read {filename}: {e}")
        return
    
    if data.ndim == 1:
        data = data.reshape(1, -1)
        
    # Columns: H, R(H), R2, Sig1, Sig2
    H = data[:, 0]
    R_vals = data[:, 1]
    
    # Find zero-field point (closest to 0)
    idx_zero = np.argmin(np.abs(H))
    R_zero = R_vals[idx_zero]
    
    # Calculate EMR % = 100 * (R - R0) / R0
    EMR = 100.0 * (R_vals - R_zero) / R_zero
    
    # Find maximum field points
    idx_max_pos = np.argmax(H)
    idx_max_neg = np.argmin(H)
    
    print(f"\nDATA SUMMARY:")
    print(f"  Number of H values: {len(H)}")
    print(f"  H range: [{H[0]:.2f}, {H[-1]:.2f}] Tesla")
    print(f"  R(H={H[idx_zero]:.4f}) = {R_zero:.6e} Ω")
    
    print(f"\nKEY RESISTANCE VALUES:")
    print(f"  H = {H[idx_max_neg]:.2f} T: R = {R_vals[idx_max_neg]:.6e} Ω, EMR = {EMR[idx_max_neg]:.2f}%")
    print(f"  H = {H[idx_zero]:.2f} T: R = {R_zero:.6e} Ω, EMR = {EMR[idx_zero]:.2f}%")
    print(f"  H = {H[idx_max_pos]:.2f} T: R = {R_vals[idx_max_pos]:.6e} Ω, EMR = {EMR[idx_max_pos]:.2f}%")
    
    print(f"\nRESISTANCE RATIOS:")
    ratio_pos = R_vals[idx_max_pos] / R_zero
    ratio_neg = R_vals[idx_max_neg] / R_zero
    print(f"  R(+H_max) / R(0) = {ratio_pos:.2f}×")
    print(f"  R(-H_max) / R(0) = {ratio_neg:.2f}×")
    
    if abs(ratio_pos - ratio_neg) / ratio_pos > 0.1:
        print(f"  ⚠️  WARNING: Asymmetry {abs(ratio_pos-ratio_neg)/ratio_pos*100:.1f}% (should be symmetric)")
    
    print(f"\nEMR MAGNITUDE:")
    EMR_max = max(abs(EMR[idx_max_pos]), abs(EMR[idx_max_neg]))
    print(f"  Peak EMR = {EMR_max:.2f}%")
    
    # Diagnosis
    print(f"\n" + "=" * 70)
    print("DIAGNOSIS:")
    print("=" * 70)
    
    if EMR_max < 1.0:
        print("❌ CRITICAL: EMR < 1%")
        print("   → Material properties likely uniform (no Hall effect)")
        print("   → Check: Fermi function returning constant?")
        print("   → Check: mu1 value in input file")
        
    elif EMR_max < 10.0:
        print("⚠️  EMR < 10% - ISSUE DETECTED")
        print("   Possible causes:")
        print("   1. Units error: Check R1, R2, t are in METERS")
        print("   2. Conductivity too low: sigma1 should be ~1e4 S/m")
        print("   3. Mobility too low: mu1 should be ~0.5 m²/(V·s)")
        print("   4. Current too small: Check Io value")
        print("   5. Hall cross-terms missing in stiffness matrix")
        
        # Additional checks
        if R_zero < 0.001:
            print(f"   → R(0) = {R_zero:.2e} < 1 mΩ suggests geometry issue (disk too large or thickness issue?)")
        if R_zero > 1000:
            print(f"   → R(0) = {R_zero:.2e} > 1 kΩ suggests geometry issue (disk too small?)")
        if ratio_pos < 1.1:
            print("   → R barely changes with H (Hall parameter too small)")
            
    elif EMR_max < 100.0:
        print(f"⚠️  EMR = {EMR_max:.1f}% - Low but possibly correct for poor geometry")
        print("   Check: Is R2/R1 ratio optimal? (Should be 0.2-0.4)")
        
    elif EMR_max < 1000.0:
        print(f"✓ EMR = {EMR_max:.1f}% - Reasonable for some geometries")
        print("  Consider: May be correct if R2/R1 is not optimal")
        
    else:
        print(f"✓ EMR = {EMR_max:.0f}% - EXCELLENT!")
        print("  Physics appears correct.")
    
    # Parabolic fit check
    print(f"\n" + "=" * 70)
    print("SHAPE ANALYSIS:")
    print("=" * 70)
    
    # Fit parabola: EMR = a·H²
    try:
        coeffs = np.polyfit(H, EMR, 2)
        EMR_fit = np.polyval(coeffs, H)
        residual = np.std(EMR - EMR_fit)
        
        print(f"  Parabolic fit: EMR ≈ {coeffs[0]:.2f}·H² + {coeffs[1]:.2e}·H + {coeffs[2]:.2e}")
        print(f"  Fit quality (RMSE): {residual:.4f}%")
        
        if residual < 0.1:
            print("  ✓ Excellent parabolic shape (Hall physics working correctly)")
        elif abs(coeffs[1]) > 0.1:
            print(f"  ⚠️  Linear term present ({coeffs[1]:.2e}) - Check for asymmetry")
    except:
        print("Could not fit parabola (too few points?)")
    
    print("=" * 70)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        analyze_emr_data(sys.argv[1])
    else:
        analyze_emr_data()
