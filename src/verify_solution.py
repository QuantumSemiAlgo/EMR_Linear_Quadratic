import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.special import jv, jn_zeros  # Bessel functions

# ============================================================================
# ANALYTICAL SOLUTION FOR LAPLACE EQUATION IN ANNULUS
# ============================================================================
# Problem: ∇²u = 0 in annular domain (Rin < r < Rout, 0 ≤ θ < 2π)
# BC: u(Rin, θ) = 0 (or floating/natural)
#     u(Rout, θ) = 10·sin(5θ)
#     Periodic in θ
#
# For BC u = A·sin(mθ) on outer radius, the solution is:
# u(r,θ) = C·r^m·sin(mθ) + D·r^(-m)·sin(mθ)
#
# For m=5, with floating inner BC (∂u/∂r = 0 at r=Rin):
# u(r,θ) = A·[f(r)]·sin(5θ)
# where f(r) satisfies the radial equation and boundary conditions
# ============================================================================

def analytical_solution(r, theta, Rin, Rout, m=5, A=10.0):
    """
    Analytical solution for Laplace equation in polar coordinates.
    """
    # Normalize theta
    theta = np.mod(theta, 2*np.pi)
    
    # Snap small angles to exactly zero
    theta[np.abs(theta) < 1e-8] = 0.0
    theta[np.abs(theta - 2*np.pi) < 1e-8] = 0.0
    
    # Treat Rin < 1e-10 as disk
    if Rin < 1e-10:
        # Disk case: C2 must be 0 for regularity at origin
        C1 = A / (Rout ** m)
        u = C1 * (r ** m) * np.sin(m * theta)
        
        # Force u=0 at zeros of sin(mθ)
        u[np.abs(np.sin(m * theta)) < 1e-10] = 0.0
        
    else:
        # Annulus case (original logic)
        denominator = Rout**m + (Rin**(2*m)) * (Rout**(-m))
        C1 = A / denominator
        C2 = C1 * (Rin**(2*m))
        
        # R(r) = C₁·r^m + C₂·r^(-m)
        R_r = C1 * (r**m) + C2 * (r**(-m))
        
        # u(r,θ) = R(r)·sin(mθ)
        u = R_r * np.sin(m * theta)
        
    return u


def read_solution_file(filename):
    """
    Read the FEM solution from solution_recon.dat
    
    Expected format for annular domain:
    # x_cart  y_cart  r  theta  value
    """
    data = np.loadtxt(filename, comments='#')
    
    if data.shape[1] == 5:
        # Annular format: x_cart, y_cart, r, theta, value
        x_cart = data[:, 0]
        y_cart = data[:, 1]
        r = data[:, 2]
        theta = data[:, 3]
        u_fem = data[:, 4]
        # Normalize theta
        theta = np.mod(theta, 2*np.pi)
        
        # Snap to exact values for small angles
        theta_threshold = 1e-8
        mask_zero = np.abs(theta) < theta_threshold
        mask_2pi = np.abs(theta - 2*np.pi) < theta_threshold
        theta[mask_zero] = 0.0
        theta[mask_2pi] = 0.0
        
        # Recompute Cartesian with corrected theta
        x_cart = r * np.cos(theta)
        y_cart = r * np.sin(theta)

        return r, theta, u_fem, x_cart, y_cart
    elif data.shape[1] == 3:
        # Rectangular format: x, y, value
        x = data[:, 0]
        y = data[:, 1]
        u_fem = data[:, 2]
        # Convert to polar
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        theta[theta < 0] += 2*np.pi
        return r, theta, u_fem, x, y
    else:
        raise ValueError(f"Unexpected data format with {data.shape[1]} columns")


def compute_errors(u_fem, u_analytical):
    """
    Compute various error metrics
    Returns errors and a valid mask
    """
    # Remove any NaN or Inf values
    mask = np.isfinite(u_fem) & np.isfinite(u_analytical)
    u_fem_clean = u_fem[mask]
    u_analytical_clean = u_analytical[mask]
    
    # Absolute error
    abs_error_full = np.zeros_like(u_fem)
    abs_error_full[:] = np.nan
    abs_error_clean = np.abs(u_fem_clean - u_analytical_clean)
    abs_error_full[mask] = abs_error_clean
    
    # Relative error (avoid division by zero)
    rel_error_full = np.zeros_like(u_fem)
    rel_error_full[:] = np.nan
    denom = np.abs(u_analytical_clean)
    denom[denom < 1e-10] = 1.0  # Avoid division by very small numbers
    rel_error_clean = abs_error_clean / denom
    rel_error_full[mask] = rel_error_clean
    
    # L2 norm error
    l2_error = np.sqrt(np.mean(abs_error_clean**2))
    l2_norm_analytical = np.sqrt(np.mean(u_analytical_clean**2))
    relative_l2_error = l2_error / l2_norm_analytical if l2_norm_analytical > 0 else 0
    
    # L_infinity (max) error
    linf_error = np.max(abs_error_clean)
    
    return {
        'abs_error': abs_error_full,
        'rel_error': rel_error_full,
        'l2_error': l2_error,
        'relative_l2_error': relative_l2_error,
        'linf_error': linf_error,
        'mean_abs_error': np.mean(abs_error_clean),
        'max_rel_error': np.max(rel_error_clean),
        'mask': mask
    }


def plot_comparison(r, theta, u_fem, u_analytical, x_cart, y_cart, errors, 
                   Rin, Rout, output_prefix='verification'):
    """
    Create comprehensive comparison plots
    """
    # Use only valid points for plotting
    mask = errors['mask']
    
    fig = plt.figure(figsize=(18, 12))
    
    # 1. FEM Solution (Cartesian)
    ax1 = fig.add_subplot(2, 3, 1)
    scatter1 = ax1.scatter(x_cart[mask], y_cart[mask], c=u_fem[mask], 
                          cmap='RdBu_r', s=2)
    ax1.set_aspect('equal')
    ax1.set_title('FEM Solution', fontsize=14, fontweight='bold')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.colorbar(scatter1, ax=ax1, label='u')
    
    # 2. Analytical Solution (Cartesian)
    ax2 = fig.add_subplot(2, 3, 2)
    scatter2 = ax2.scatter(x_cart[mask], y_cart[mask], c=u_analytical[mask], 
                          cmap='RdBu_r', s=2)
    ax2.set_aspect('equal')
    ax2.set_title('Analytical Solution', fontsize=14, fontweight='bold')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.colorbar(scatter2, ax=ax2, label='u')
    
    # 3. Absolute Error (Cartesian)
    ax3 = fig.add_subplot(2, 3, 3)
    abs_error_valid = errors['abs_error'][mask]
    scatter3 = ax3.scatter(x_cart[mask], y_cart[mask], c=abs_error_valid, 
                          cmap='hot', s=2, vmin=0)
    ax3.set_aspect('equal')
    ax3.set_title(f'Absolute Error\nMax: {errors["linf_error"]:.2e}', 
                  fontsize=14, fontweight='bold')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    plt.colorbar(scatter3, ax=ax3, label='|u_FEM - u_exact|')
    
    # 4. Solution comparison along outer radius
    ax4 = fig.add_subplot(2, 3, 4)
    # Find points near outer radius
    outer_mask = (np.abs(r - Rout) < 0.01 * (Rout - Rin)) & mask
    if np.sum(outer_mask) > 0:
        theta_outer = theta[outer_mask]
        u_fem_outer = u_fem[outer_mask]
        u_analytical_outer = u_analytical[outer_mask]
        
        # Sort by theta for plotting
        sort_idx = np.argsort(theta_outer)
        ax4.plot(theta_outer[sort_idx], u_fem_outer[sort_idx], 'b.-', 
                 label='FEM', linewidth=2, markersize=4)
        ax4.plot(theta_outer[sort_idx], u_analytical_outer[sort_idx], 'r--', 
                 label='Analytical', linewidth=2)
    ax4.set_xlabel('θ (radians)', fontsize=12)
    ax4.set_ylabel('u', fontsize=12)
    ax4.set_title(f'Solution at Outer Radius (r={Rout:.3f})', 
                  fontsize=14, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. Solution comparison along inner radius
    ax5 = fig.add_subplot(2, 3, 5)
    inner_mask = (np.abs(r - Rin) < 0.01 * (Rout - Rin)) & mask
    if np.sum(inner_mask) > 0:
        theta_inner = theta[inner_mask]
        u_fem_inner = u_fem[inner_mask]
        u_analytical_inner = u_analytical[inner_mask]
        
        sort_idx = np.argsort(theta_inner)
        ax5.plot(theta_inner[sort_idx], u_fem_inner[sort_idx], 'b.-', 
                 label='FEM', linewidth=2, markersize=4)
        ax5.plot(theta_inner[sort_idx], u_analytical_inner[sort_idx], 'r--', 
                 label='Analytical', linewidth=2)
        ax5.set_xlabel('θ (radians)', fontsize=12)
        ax5.set_ylabel('u', fontsize=12)
        ax5.set_title(f'Solution at Inner Radius (r={Rin:.3f})', 
                      fontsize=14, fontweight='bold')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
    
    # 6. Error statistics
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.axis('off')
    
    stats_text = f"""
    ERROR STATISTICS
    ════════════════════════════════════
    
    L² Error:          {errors['l2_error']:.4e}
    Relative L² Error: {errors['relative_l2_error']:.4e}
    
    L∞ Error (max):    {errors['linf_error']:.4e}
    Mean Abs Error:    {errors['mean_abs_error']:.4e}
    Max Rel Error:     {errors['max_rel_error']:.4e}
    
    ────────────────────────────────────
    Domain Parameters:
    Inner Radius:  {Rin:.4f}
    Outer Radius:  {Rout:.4f}
    
    Number of Points: {np.sum(mask)}
    """
    
    ax6.text(0.1, 0.5, stats_text, fontsize=11, family='monospace',
             verticalalignment='center', transform=ax6.transAxes)
    
    plt.tight_layout()
    
    # Ensure figures directory exists
    figures_dir = '../output/figures'
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)
        
    save_path = os.path.join(figures_dir, f'{output_prefix}_comparison.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved comparison plot to {save_path}")
    
    return fig


def plot_radial_profile(r, theta, u_fem, u_analytical, mask, Rin, Rout, 
                       output_prefix='verification'):
    """
    Plot solution along radial lines at different theta values
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    theta_values = [0, np.pi/6, np.pi/3, np.pi/2, 2*np.pi/3, 5*np.pi/6]
    
    for idx, theta_target in enumerate(theta_values):
        ax = axes[idx]
        
        # Find points near this theta value
        theta_tol = 0.1
        line_mask = (np.abs(theta - theta_target) < theta_tol) & mask
        
        if np.sum(line_mask) > 10:
            r_line = r[line_mask]
            u_fem_line = u_fem[line_mask]
            u_analytical_line = u_analytical[line_mask]
            
            # Sort by r
            sort_idx = np.argsort(r_line)
            
            ax.plot(r_line[sort_idx], u_fem_line[sort_idx], 'b.-', 
                   label='FEM', linewidth=2, markersize=4)
            ax.plot(r_line[sort_idx], u_analytical_line[sort_idx], 'r--', 
                   label='Analytical', linewidth=2)
            ax.set_xlabel('r', fontsize=11)
            ax.set_ylabel('u', fontsize=11)
            ax.set_title(f'θ = {theta_target:.3f} rad ({np.degrees(theta_target):.1f}°)', 
                        fontsize=12, fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_xlim([Rin*0.95, Rout*1.05])
    
    plt.tight_layout()
    
    # Ensure figures directory exists
    figures_dir = '../output/figures'
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    save_path = os.path.join(figures_dir, f'{output_prefix}_radial_profiles.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved radial profiles to {save_path}")
    
    return fig


def main():
    """
    Main verification routine
    """
    print("=" * 70)
    print("FEM SOLUTION VERIFICATION AGAINST ANALYTICAL SOLUTION")
    print("=" * 70)
    
    # ========================================================================
    # CONFIGURATION - ADJUST THESE VALUES TO MATCH YOUR SIMULATION
    # ========================================================================
    solution_file = '../output/solution_recon.dat'  # Path to FEM solution
    Rin = 0.0    # Inner radius - MUST MATCH YOUR INPUT FILE
    Rout = 10.0   # Outer radius - MUST MATCH YOUR INPUT FILE
    m = 5        # Mode number (from BC: 10*sin(5*theta))
    A = 10.0     # Amplitude (from BC: 10*sin(5*theta))
    
    # ========================================================================
    # READ FEM SOLUTION
    # ========================================================================
    print(f"\nReading FEM solution from: {solution_file}")
    try:
        r, theta, u_fem, x_cart, y_cart = read_solution_file(solution_file)
        print(f"  ✓ Successfully loaded {len(u_fem)} data points")
        print(f"  ✓ r range: [{np.min(r):.4f}, {np.max(r):.4f}]")
        print(f"  ✓ θ range: [{np.min(theta):.4f}, {np.max(theta):.4f}]")
    except FileNotFoundError:
        print(f"  ✗ ERROR: File '{solution_file}' not found!")
        print(f"    Please ensure the file exists in the current directory.")
        return
    except Exception as e:
        print(f"  ✗ ERROR: Failed to read file: {e}")
        return
    
    # -------------------------------------------------------------------------
    # ========================================================================
    # COMPUTE ANALYTICAL SOLUTION
    # ========================================================================
    print(f"\nComputing analytical solution...")
    print(f"  Domain: Annulus with Rin={Rin:.4f}, Rout={Rout:.4f}")
    print(f"  BC: u(Rout, θ) = {A}·sin({m}θ), floating inner BC")
    
    u_analytical = analytical_solution(r, theta, Rin, Rout, m=m, A=A)
    print(f"  ✓ Analytical solution computed")
    
    # ========================================================================
    # COMPUTE ERRORS
    # ========================================================================
    print(f"\nComputing error metrics...")
    errors = compute_errors(u_fem, u_analytical)
    
    print(f"\n" + "=" * 70)
    print("ERROR ANALYSIS RESULTS")
    print("=" * 70)
    print(f"L² Error:               {errors['l2_error']:.6e}")
    print(f"Relative L² Error:      {errors['relative_l2_error']:.6e}")
    print(f"L∞ Error (maximum):     {errors['linf_error']:.6e}")
    print(f"Mean Absolute Error:    {errors['mean_abs_error']:.6e}")
    print(f"Maximum Relative Error: {errors['max_rel_error']:.6e}")
    
    # Find location of max error
    diff = np.abs(u_fem - u_analytical)
    max_err_idx = np.argmax(diff)
    r_flat = r.flatten()
    theta_flat = theta.flatten()
    print(f"Location of Max Error:  r={r_flat[max_err_idx]:.4f}, theta={theta_flat[max_err_idx]:.4f}")
    print(f"                        u_fem={u_fem.flatten()[max_err_idx]:.4f}, u_ana={u_analytical.flatten()[max_err_idx]:.4f}")
    
    # ========================================================================
    # CONVERGENCE ASSESSMENT
    # ========================================================================
    print(f"\n" + "=" * 70)
    print("CONVERGENCE ASSESSMENT")
    print("=" * 70)
    
    if errors['relative_l2_error'] < 1e-2:
        print("✓✓✓ EXCELLENT: Relative L² error < 1%")
    elif errors['relative_l2_error'] < 5e-2:
        print("✓✓  GOOD: Relative L² error < 5%")
    elif errors['relative_l2_error'] < 1e-1:
        print("✓   ACCEPTABLE: Relative L² error < 10%")
    else:
        print("✗   WARNING: Relative L² error > 10%")
        print("    Consider refining the mesh or checking boundary conditions")
    
    # ========================================================================
    # GENERATE PLOTS
    # ========================================================================
    print(f"\nGenerating comparison plots...")
    
    fig1 = plot_comparison(r, theta, u_fem, u_analytical, x_cart, y_cart, 
                          errors, Rin, Rout)
    
    fig2 = plot_radial_profile(r, theta, u_fem, u_analytical, errors['mask'],
                               Rin, Rout)
    
    print(f"\n" + "=" * 70)
    print("VERIFICATION COMPLETE")
    print("=" * 70)
    print(f"Plots saved:")
    print(f"  - verification_comparison.png")
    print(f"  - verification_radial_profiles.png")
    
    # Show plots
    plt.show()


if __name__ == "__main__":
    main()