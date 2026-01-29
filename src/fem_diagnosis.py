#!/usr/bin/env python3
"""
Quick diagnostic to check FEM solution data
"""
import numpy as np
import matplotlib.pyplot as plt
import os

# Read FEM solution
print("Reading FEM solution...")
try:
    data = np.loadtxt('../output/solution_recon.dat', comments='#')
except Exception as e:
    print(f"Error reading file: {e}")
    exit(1)

print(f"Data shape: {data.shape}")
try:
    print(f"Columns: {data.shape[1]}")
except IndexError:
    print("Error: Data has no columns (empty or scalar?)")
    exit(1)

if data.shape[1] == 5:
    x_cart, y_cart, r, theta, u_fem = data.T
    print("\nFormat: x_cart, y_cart, r, theta, u_fem")
elif data.shape[1] == 3:
    x, y, u_fem = data.T
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    theta[theta < 0] += 2*np.pi
    print("\nFormat: x, y, u_fem (converted to polar)")
else:
    print(f"ERROR: Unexpected data format with {data.shape[1]} columns")
    exit(1)

# Check radial range
print(f"\nRadial coordinate statistics:")
print(f"  Min r: {r.min():.6f}")
print(f"  Max r: {r.max():.6f}")
print(f"  Mean r: {r.mean():.6f}")
print(f"  Median r: {np.median(r):.6f}")

# Count points near r=0
n_near_zero = np.sum(r < 0.2)
n_at_zero = np.sum(r < 0.01)
print(f"\nPoints near origin:")
print(f"  r < 0.2: {n_near_zero} points")
print(f"  r < 0.01: {n_at_zero} points")

# Check solution values
print(f"\nSolution statistics:")
print(f"  Min u: {u_fem.min():.6f}")
print(f"  Max u: {u_fem.max():.6f}")
print(f"  Mean u: {u_fem.mean():.6f}")
print(f"  u at min(r): {u_fem[np.argmin(r)]:.6f}")
print(f"  u at max(r): {u_fem[np.argmax(r)]:.6f}")

# Check for NaN or Inf
n_nan = np.sum(np.isnan(u_fem))
n_inf = np.sum(np.isinf(u_fem))
print(f"\nData quality:")
print(f"  NaN values: {n_nan}")
print(f"  Inf values: {n_inf}")

# Plot radial distribution
plt.figure(figsize=(12, 4))

plt.subplot(131)
plt.hist(r, bins=50)
plt.xlabel('r')
plt.ylabel('Count')
plt.title('Distribution of radial points')
plt.grid(True)

plt.subplot(132)
plt.scatter(r, u_fem, alpha=0.3, s=1)
plt.xlabel('r')
plt.ylabel('u_fem')
plt.title('FEM solution vs radius')
plt.grid(True)

plt.subplot(133)
# Get points at theta ≈ 90° (where sin(5θ) has interesting behavior)
theta_target = np.pi/2
mask = np.abs(theta - theta_target) < 0.1
if np.sum(mask) > 10:
    r_slice = r[mask]
    u_slice = u_fem[mask]
    idx = np.argsort(r_slice)
    plt.plot(r_slice[idx], u_slice[idx], 'b.-', markersize=2)
    plt.xlabel('r')
    plt.ylabel('u_fem')
    plt.title(f'Radial profile at θ≈{90}°')
    plt.grid(True)
    
    # Add expected analytical for comparison
    r_analytical = r_slice[idx]
    u_analytical = 10.0 * (r_analytical / 10.0)**5
    plt.plot(r_analytical, u_analytical, 'r--', label='Analytical')
    plt.legend()
else:
    plt.text(0.5, 0.5, 'Not enough points\nat θ≈90°', 
             ha='center', va='center', transform=plt.gca().transAxes)

plt.tight_layout()
plt.tight_layout()

# Ensure figures directory exists
figures_dir = '../output/figures'
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

save_path = os.path.join(figures_dir, 'diagnostic_plot.png')
plt.savefig(save_path, dpi=150)
print(f"\nSaved diagnostic plot to {save_path}")

# Print first few and last few lines
print("\n" + "="*70)
print("FIRST 10 DATA POINTS:")
print("="*70)
if data.shape[1] == 5:
    print(f"{'x_cart':>10} {'y_cart':>10} {'r':>10} {'theta':>10} {'u_fem':>10}")
    for i in range(min(10, len(data))):
        print(f"{data[i,0]:10.4f} {data[i,1]:10.4f} {data[i,2]:10.4f} "
              f"{data[i,3]:10.4f} {data[i,4]:10.4f}")
else:
    print(f"{'x':>10} {'y':>10} {'u_fem':>10}")
    for i in range(min(10, len(data))):
        print(f"{data[i,0]:10.4f} {data[i,1]:10.4f} {data[i,2]:10.4f}")

print("\n" + "="*70)
print("LAST 10 DATA POINTS:")
print("="*70)
if data.shape[1] == 5:
    print(f"{'x_cart':>10} {'y_cart':>10} {'r':>10} {'theta':>10} {'u_fem':>10}")
    for i in range(max(0, len(data)-10), len(data)):
        print(f"{data[i,0]:10.4f} {data[i,1]:10.4f} {data[i,2]:10.4f} "
              f"{data[i,3]:10.4f} {data[i,4]:10.4f}")
else:
    print(f"{'x':>10} {'y':>10} {'u_fem':>10}")
    for i in range(max(0, len(data)-10), len(data)):
        print(f"{data[i,0]:10.4f} {data[i,1]:10.4f} {data[i,2]:10.4f}")

print("\n" + "="*70)
print("POINTS WITH SMALLEST r:")
print("="*70)
idx_smallest = np.argsort(r)[:10]
if data.shape[1] == 5:
    print(f"{'x_cart':>10} {'y_cart':>10} {'r':>10} {'theta':>10} {'u_fem':>10}")
    for i in idx_smallest:
        print(f"{data[i,0]:10.4f} {data[i,1]:10.4f} {data[i,2]:10.4f} "
              f"{data[i,3]:10.4f} {data[i,4]:10.4f}")
else:
    print(f"{'r':>10} {'theta':>10} {'u_fem':>10}")
    for i in idx_smallest:
        print(f"{r[i]:10.4f} {theta[i]:10.4f} {u_fem[i]:10.4f}")

print("\n" + "="*70)
print("DIAGNOSIS COMPLETE")
print("="*70)
