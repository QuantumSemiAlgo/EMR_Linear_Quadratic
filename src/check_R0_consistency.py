#!/usr/bin/env python3
import numpy as np
import sys

try:
    data = np.loadtxt('../output/EMRdata.out', comments='#')
except OSError:
    try:
        data = np.loadtxt('EMRdata.out', comments='#')
    except OSError:
        print("Error: Could not find EMRdata.out")
        sys.exit(1)

if data.ndim == 1:
    data = data.reshape(1, -1)

# Assuming columns: H, R0, R(H), EMR
H = data[:, 0]
R0 = data[:, 1]
RH = data[:, 2]
EMR = data[:, 3]

print("=" * 60)
print("R0 CONSISTENCY CHECK")
print("=" * 60)

# Check if R0 is constant
R0_mean = np.mean(R0)
R0_std = np.std(R0)
R0_variation = 100 * R0_std / R0_mean if R0_mean != 0 else 0

print(f"R0 mean:  {R0_mean:.6e} Ohm")
print(f"R0 std:   {R0_std:.6e} Ohm")
print(f"R0 variation: {R0_variation:.4f}%")

if R0_variation < 0.001:
    print("✓ R0 is CONSTANT (good!)")
else:
    print("✗ R0 is VARYING (BUG!)")
    print(f"\nR0 values (first 5 and last 5):")
    for i in range(min(5, len(H))):
        print(f"  H = {H[i]:+.4f} T: R0 = {R0[i]:.6e} Ohm")
    if len(H) > 5:
        print("  ...")
        for i in range(max(5, len(H)-5), len(H)):
            print(f"  H = {H[i]:+.4f} T: R0 = {R0[i]:.6e} Ohm")

# Check EMR at H=0
idx_zero = np.argmin(np.abs(H))
H_at_zero = H[idx_zero]
EMR_at_zero = EMR[idx_zero]
R_at_zero = RH[idx_zero]

print(f"\nAt H ~= {H_at_zero:.4f} T:")
print(f"  R(H={H_at_zero:.4f}) = {R_at_zero:.6e} Ohm")
print(f"  R0          = {R0[idx_zero]:.6e} Ohm")
print(f"  EMR         = {EMR_at_zero:.4f}%")

if abs(EMR_at_zero) < 0.01:
    print("✓ EMR = 0 at H = 0 (good!)")
else:
    print(f"✗ EMR != 0 at H = 0 (BUG! Should be 0, got {EMR_at_zero:.4f}%)")

# Fit parabola
try:
    coeffs = np.polyfit(H, EMR, 2)
    print(f"\nParabolic fit: EMR = {coeffs[0]:.2f}*H^2 + {coeffs[1]:.2e}*H + {coeffs[2]:.2e}")

    if abs(coeffs[2]) < 0.01:
        print("✓ Zero offset (good!)")
    else:
        print(f"✗ Non-zero offset {coeffs[2]:.2f}% (BUG!)")
except:
    print("\nCould not fit parabola (not enough points?)")

print("=" * 60)
