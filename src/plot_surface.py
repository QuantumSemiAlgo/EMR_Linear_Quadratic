#!/usr/bin/env python3
"""
Surface Plotter for FEM Solution
Reads solution_recon.dat and generates 2D and 3D surface visualizations.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import sys

def read_solution_file(filename):
    """
    Read the FEM solution from solution_recon.dat
    Handles both Annular (5 cols) and Rectangular (3 cols) formats.
    """
    print(f"Reading solution from: {filename}")
    try:
        data = np.loadtxt(filename, comments='#')
    except OSError:
        print(f"Error: Could not find file {filename}")
        sys.exit(1)

    if data.shape[1] == 5:
        # Annular format: x_cart, y_cart, r, theta, value
        x = data[:, 0]
        y = data[:, 1]
        u = data[:, 4]
        return x, y, u, "Annular"
    elif data.shape[1] == 3:
        # Rectangular format: x, y, value
        x = data[:, 0]
        y = data[:, 1]
        u = data[:, 2]
        return x, y, u, "Rectangular"
    else:
        print(f"Error: Unexpected data shape {data.shape}")
        sys.exit(1)

def plot_3d_surface(x, y, u, output_dir):
    """
    Creates a 3D surface plot of the solution.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create the surface plot
    # utilizing triangulation since data might not be on a perfect grid for plotting directly
    surf = ax.plot_trisurf(x, y, u, cmap='viridis', edgecolor='none', alpha=0.9, linewidth=0, antialiased=True)
    
    ax.set_title("3D Surface Plot of Solution u(x,y)", fontsize=14, fontweight='bold')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("u")
    
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label='Solution u')
    
    # Save
    save_path = os.path.join(output_dir, "surface_plot_3d.png")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Saved 3D plot to: {save_path}")
    plt.close(fig)

def plot_2d_contour(x, y, u, output_dir):
    """
    Creates a 2D filled contour plot (heatmap) of the solution.
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create filled contour plot using triangulation
    # standard scatter is essentially points, tricontourf interpolates to fill
    contour = ax.tricontourf(x, y, u, levels=50, cmap='viridis')
    
    ax.set_title("2D Contour Plot of Solution u(x,y)", fontsize=14, fontweight='bold')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect('equal')
    
    # Add colorbar
    cbar = fig.colorbar(contour, ax=ax, label='Solution u')
    
    # Add points (optional, maybe too dense)
    # ax.scatter(x, y, c='k', s=0.1, alpha=0.1) 
    
    # Save
    save_path = os.path.join(output_dir, "surface_plot_2d.png")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Saved 2D plot to: {save_path}")
    plt.close(fig)

def main():
    # Configuration
    # Assumes script is in src/ and output is in ../output/
    solution_path = "../output/solution_recon.dat"
    figures_dir = "../output/figures"

    # Ensure output directory exists
    if not os.path.exists(figures_dir):
        try:
            os.makedirs(figures_dir)
            print(f"Created directory: {figures_dir}")
        except OSError as e:
            print(f"Error creating directory {figures_dir}: {e}")
            sys.exit(1)

    # 1. Read Data
    x, y, u, geom_type = read_solution_file(solution_path)
    print(f"Loaded {len(u)} points (Geometry: {geom_type})")
    print(f"Range: u in [{u.min():.4f}, {u.max():.4f}]")

    # 2. Generate 3D Surface Plot
    print("Generating 3D surface plot...")
    plot_3d_surface(x, y, u, figures_dir)

    # 3. Generate 2D Contour Plot
    print("Generating 2D contour plot...")
    plot_2d_contour(x, y, u, figures_dir)
    
    print("\nVisualization complete.")

if __name__ == "__main__":
    main()
