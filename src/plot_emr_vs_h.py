
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def parse_emr_data(filepath):
    """
    Parses EMRdata.out.
    Format: H(T)  Resistance(Ohm)  R2(m)  Sigma1(S/m)  Sigma2(S/m)
    """
    try:
        data = np.loadtxt(filepath, skiprows=1, comments='#')
        if data.size == 0:
            return None, None, None
        if data.ndim == 1:
            data = data.reshape(1, -1)
        
        # Check dimensions to determine format
        cols = data.shape[1]
        
        if cols >= 4:
            # New format: H(T)  Ro(Ohm)  R(H)(Ohm)  EMR(%)
            print("Detected 4-column format (H, Ro, R, EMR)")
            H = data[:, 0]
            R = data[:, 2] # Column 2 is R(H)
            EMR = data[:, 3]
            return H, R, EMR
        else:
            # Old format: H(T)  Resistance(Ohm) ...
            print("Detected old format (H, R, ...)")
            H = data[:, 0]
            R = data[:, 1]
            return H, R, None
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None, None, None

def calculate_emr(H, R):
    # Find R0 (R at minimal |H|, ideally H=0)
    # Assuming the sweep includes H=0 or close to it.
    idx_min_H = np.argmin(np.abs(H))
    R0 = R[idx_min_H]
    
    if R0 == 0:
        print("Warning: R0 is zero. EMR calculation might vary.")
        R0 = 1e-9 # Avoid div by zero
        
    emr_percent = ((R - R0) / R0) * 100.0
    return R0, emr_percent

def plot_emr(H, R, EMR_file=None, output_prefix="EMR_Plot"):
    if EMR_file is not None:
        emr = EMR_file
        # Use first Ro from file if possible, else calc
        R0 = R[np.argmin(np.abs(H))]
        print("Using pre-calculated EMR from file.")
    else:
        R0, emr = calculate_emr(H, R)
        print("Calculated EMR from R(H).")
    
    # 1. Resistance vs H
    plt.figure(figsize=(10, 6))
    plt.plot(H, R * 1000.0, 'bo-', linewidth=2, markersize=6) # R in mOhm
    plt.xlabel('Magnetic Field B (T)', fontsize=14)
    plt.ylabel('Resistance (mΩ)', fontsize=14)
    plt.title('Resistance vs Magnetic Field', fontsize=16)
    plt.grid(True, which='both', ls='-', alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"../output/{output_prefix}_R_vs_H.png")
    print(f"Saved ../output/{output_prefix}_R_vs_H.png")
    
    # 2. EMR vs H
    plt.figure(figsize=(10, 6))
    plt.plot(H, emr, 'rd-', linewidth=2, markersize=6)
    plt.xlabel('Magnetic Field B (T)', fontsize=14)
    plt.ylabel('EMR (%)', fontsize=14)
    if EMR_file is not None:
         plt.title(f'Extraordinary Magnetoresistance (EMR)', fontsize=16)
    else:
         plt.title(f'Extraordinary Magnetoresistance (EMR)\nR0 = {R0*1000:.4f} mΩ', fontsize=16)
    plt.grid(True, which='both', ls='-', alpha=0.5)
    
    # Fit parabola (Removed by request)
    # try:
    #     p = np.polyfit(H, emr, 2)
    #     H_fit = np.linspace(min(H), max(H), 100)
    #     emr_fit = np.polyval(p, H_fit)
    #     plt.plot(H_fit, emr_fit, 'k--', label=f'Fit: {p[0]:.2f}H² + {p[1]:.2f}H + {p[2]:.2f}')
    #     plt.legend()
    # except:
    #     pass

    plt.tight_layout()
    plt.savefig(f"../output/{output_prefix}_EMR_vs_H.png")
    print(f"Saved ../output/{output_prefix}_EMR_vs_H.png")

if __name__ == "__main__":
    filepath = "../output/EMRdata.out"
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        sys.exit(1)
        
    print(f"Analyzing {filepath}...")
    H, R, EMR_val = parse_emr_data(filepath)
    
    if H is not None:
        plot_emr(H, R, EMR_file=EMR_val)
        
        # Print summary
        peak_emr = 0.0
        if EMR_val is not None:
            peak_emr = np.max(EMR_val)
        else:
             _, emr_calc = calculate_emr(H, R)
             peak_emr = np.max(emr_calc)

        idx_peak = np.argmax(R)
        print(f"\nSummary:")
        print(f"  Field Range: {np.min(H):.2f} T to {np.max(H):.2f} T")
        print(f"  Peak Resistance: {np.max(R):.6e} Ohm at {H[idx_peak]:.2f} T")
        print(f"  Peak EMR: {peak_emr:.2f} %")
