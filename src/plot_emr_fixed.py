#!/usr/bin/env python3
"""
Plot EMR vs H with proper error checking
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

def load_and_verify_data(filename='../output/EMRdata.out'):
    """Load data and verify R0 consistency"""
    
    print("=" * 70)
    print("LOADING AND VERIFYING EMR DATA")
    print("=" * 70)
    
    try:
        data = np.loadtxt(filename, comments='#')
    except Exception as e:
        # Try local path if ../output fails
        try:
             data = np.loadtxt('EMRdata.out', comments='#')
        except:
            print(f"ERROR: Could not load {filename}")
            print(f"Error: {e}")
            sys.exit(1)
    
    # Handle single row
    if data.ndim == 1:
        data = data.reshape(1, -1)
    
    # Extract columns
    if data.shape[1] < 4:
        print(f"ERROR: Expected 4 columns, got {data.shape[1]}")
        print("File format should be: H  R0  R(H)  EMR")
        sys.exit(1)
    
    H = data[:, 0]
    R0 = data[:, 1]
    RH = data[:, 2]
    EMR = data[:, 3]
    
    print(f"Loaded {len(H)} data points")
    print(f"H range: [{H.min():.4f}, {H.max():.4f}] T")
    
    # CRITICAL CHECK 1: Is R0 constant?
    R0_mean = np.mean(R0)
    R0_std = np.std(R0)
    R0_variation_pct = 100.0 * R0_std / R0_mean if R0_mean != 0 else 0
    
    print(f"\nR0 Statistics:")
    print(f"  Mean:      {R0_mean:.6e} Ω")
    print(f"  Std Dev:   {R0_std:.6e} Ω")
    print(f"  Variation: {R0_variation_pct:.6f}%")
    
    if R0_variation_pct > 0.001:
        print("  ⚠️  WARNING: R0 is not constant!")
        print("  R0 values:")
        for i in range(min(10, len(H))):
            print(f"    H = {H[i]:+.4f} T: R0 = {R0[i]:.6e} Ω")
        if len(H) > 10:
            print(f"    ... ({len(H)-10} more rows)")
    else:
        print("  ✓ R0 is constant (good!)")
    
    # CRITICAL CHECK 2: Is EMR = 0 at H = 0?
    idx_zero = np.argmin(np.abs(H))
    H_zero = H[idx_zero]
    EMR_zero = EMR[idx_zero]
    RH_zero = RH[idx_zero]
    
    print(f"\nAt H ≈ {H_zero:.6f} T:")
    print(f"  R(H=0) = {RH_zero:.6e} Ω")
    print(f"  R0     = {R0[idx_zero]:.6e} Ω")
    print(f"  EMR    = {EMR_zero:.6f}%")
    
    if abs(EMR_zero) < 1e-6:
        print("  ✓ EMR = 0 at H = 0 (good!)")
    else:
        print(f"  ⚠️  EMR ≠ 0 at H = 0 (should be 0, got {EMR_zero:.6f}%)")
    
    # CRITICAL CHECK 3: Symmetry
    # Check max value magnitude in pos vs neg H
    EMR_pos_vals = EMR[H > 1e-6]
    EMR_neg_vals = EMR[H < -1e-6]
    
    EMR_pos_max = np.max(EMR_pos_vals) if len(EMR_pos_vals) > 0 else 0
    EMR_neg_max = np.max(EMR_neg_vals) if len(EMR_neg_vals) > 0 else 0
    
    asymmetry_pct = 100 * abs(EMR_pos_max - EMR_neg_max) / max(EMR_pos_max, 1e-10)
    
    print(f"\nSymmetry Check:")
    print(f"  EMR(+H_max) = {EMR_pos_max:.4f}%")
    print(f"  EMR(-H_max) = {EMR_neg_max:.4f}%")
    print(f"  Asymmetry:  {asymmetry_pct:.2f}%")
    
    if asymmetry_pct < 5.0:
        print("  ✓ Symmetric (good!)")
    else:
        print("  ⚠️  Asymmetric!")
    
    print("=" * 70)
    
    return H, R0, RH, EMR, R0_mean

def plot_emr(H, R0, RH, EMR, R0_mean):
    """Create EMR plot with diagnostics"""
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot data
    ax.plot(H, EMR, 'o', color='red', markersize=8, 
            label='Simulation Data', zorder=3)
    
    H_fit = np.linspace(H.min(), H.max(), 200)

    # Model 1: Parabolic (Classic Low Field)
    # EMR = a * H^2
    H_nonzero = H[np.abs(H) > 1e-6]
    EMR_nonzero = EMR[np.abs(H) > 1e-6]
    if len(H_nonzero) > 0:
        a_quad = np.mean(EMR_nonzero / (H_nonzero**2))
        EMR_quad = a_quad * H_fit**2
        ax.plot(H_fit, EMR_quad, '--', color='blue', linewidth=2,
                label=f'Parabolic Fit (aH²): a={a_quad:.2f}')
                
        # Model 2: Linear (High Field / Saturation)
        # EMR = b * |H|
        b_lin = np.mean(EMR_nonzero / np.abs(H_nonzero))
        EMR_lin = b_lin * np.abs(H_fit)
        ax.plot(H_fit, EMR_lin, ':', color='green', linewidth=2,
                label=f'Linear Fit (b|H|): b={b_lin:.2f}')

    # Model 3: User's Polynomial (to show offset)
    coeffs = np.polyfit(H, EMR, 2)
    EMR_poly = np.polyval(coeffs, H_fit)
    # ax.plot(H_fit, EMR_poly, '-', color='gray', alpha=0.5, label='PolyFit (Offset)')
    
    ax.set_xlabel('Magnetic Field H (T)', fontsize=14, fontweight='bold')
    ax.set_ylabel('EMR (%)', fontsize=14, fontweight='bold')
    ax.set_title(f'Extraordinary Magnetoresistance\nR0 = {R0_mean*1e3:.4f} mΩ', 
                 fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11, loc='upper center')
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    
    # Save
    output_file = '../output/EMR_vs_H_DIAGNOSTIC.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Plot saved to: {output_file}\n")
    
    # Diagnostics
    print("=" * 70)
    print("PHYSICS MODEL CHECK")
    print("=" * 70)
    
    # Check R^2 for both
    # Parabolic
    EMR_pred_quad = a_quad * H**2
    r2_quad = 1 - np.sum((EMR - EMR_pred_quad)**2) / np.sum((EMR - np.mean(EMR))**2)
    
    # Linear
    EMR_pred_lin = b_lin * np.abs(H)
    r2_lin = 1 - np.sum((EMR - EMR_pred_lin)**2) / np.sum((EMR - np.mean(EMR))**2)
    
    print(f"Parabolic Fit R²: {r2_quad:.4f}")
    print(f"Linear Fit R²:    {r2_lin:.4f}")
    
    if r2_lin > r2_quad:
         print("\nCONCLUSION: Data is LINEAR (|H|), not Parabolic (H²).")
         print("This explains why a parabolic fit forces an offset!")
    else:
         print("\nCONCLUSION: Data is Parabolic.")

    print("=" * 70)

if __name__ == '__main__':
    H, R0, RH, EMR, R0_mean = load_and_verify_data()
    plot_emr(H, R0, RH, EMR, R0_mean)
