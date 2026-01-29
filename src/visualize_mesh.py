#!/usr/bin/env python3
"""
Mesh Visualization Tool for Annular FEM Domains

This script reads mesh files (node.dat, elem.dat, bnode.dat) and creates
comprehensive visualization plots of the mesh geometry, including:
- Full mesh with element edges
- Boundary nodes highlighted
- Element quality analysis
- Mesh statistics

Author: Auto-generated
Date: 2026-01-29
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import LineCollection
import os
import sys
from pathlib import Path
from io import StringIO

# Set matplotlib style for publication-quality plots
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 9


class MeshVisualizer:
    """Class to handle mesh visualization for annular FEM domains"""
    
    def __init__(self, mesh_dir='mesh', node_elem=9):
        """
        Initialize mesh visualizer
        
        Parameters:
        -----------
        mesh_dir : str
            Directory containing mesh files
        node_elem : int
            Number of nodes per element (4 or 9)
        """
        self.mesh_dir = Path(mesh_dir)
        self.node_elem = node_elem
        self.nodes = None
        self.elements = None
        self.boundary_nodes = None
        self.Rin = 0.0
        self.Rout = 1.0
        self.use_annular = False
        self.domain_params = {}
        
    def read_input_parameters(self):
        """Read parameters from FEMstruct2d.inp if it exists"""
        inp_file = Path('FEMstruct2d.inp')
        if not inp_file.exists():
            # Try parent directory
            inp_file = Path('../FEMstruct2d.inp')
            
        if not inp_file.exists():
            print("  ⚠ Input file FEMstruct2d.inp not found. Assuming generic mesh.")
            return

        print(f"Reading parameters from {inp_file}...")
        try:
            with open(inp_file, 'r') as f:
                lines = [l.strip() for l in f if l.strip() and not l.strip().startswith('#')]
            
            # This is a simple parser based on the known file structure
            # It relies on the order of parameters matching the file exactly
            # 1. Debug level
            # 2. xmin, xmax
            vals = lines[1].split()
            self.domain_params['xmin'] = float(vals[0])
            self.domain_params['xmax'] = float(vals[1])
            
            # 3. ymin, ymax
            vals = lines[2].split()
            self.domain_params['ymin'] = float(vals[0])
            self.domain_params['ymax'] = float(vals[1])
            
            # 4. use_annular
            self.use_annular = int(lines[3]) == 1
            
            # 5. Rin
            self.Rin = float(lines[4])
            
            # 6. Rout
            self.Rout = float(lines[5])
            
            print(f"  ✓ Domain: Annular={self.use_annular}, R=[{self.Rin}, {self.Rout}]")
            
        except Exception as e:
            print(f"  ⚠ Error parsing input file: {e}")
            print("  Using default/detected parameters.")

    def read_mesh_files(self):
        """Read node.dat, elem.dat, and bnode.dat"""
        
        # First try to read input parameters to know about geometry
        self.read_input_parameters()
        
        print(f"Reading mesh files from {self.mesh_dir}...")
        
        # Read node.dat
        node_file = self.mesh_dir / 'node.dat'
        if not node_file.exists():
            raise FileNotFoundError(f"Node file not found: {node_file}")
            
        self.nodes = self._read_data_file(node_file)
        print(f"  ✓ Loaded {len(self.nodes)} nodes")
        
        # Transform coordinates if annular
        if self.use_annular and 'xmax' in self.domain_params:
            print("  ⟳ Transforming logical coordinates to physical annular mesh...")
            self._transform_coordinates()
        
        # Read elem.dat
        elem_file = self.mesh_dir / 'elem.dat'
        if not elem_file.exists():
            raise FileNotFoundError(f"Element file not found: {elem_file}")
            
        self.elements = self._read_data_file(elem_file, dtype=int)
        print(f"  ✓ Loaded {len(self.elements)} elements")
        
        # Read bnode.dat
        bnode_file = self.mesh_dir / 'bnode.dat'
        if not bnode_file.exists():
            raise FileNotFoundError(f"Boundary node file not found: {bnode_file}")
            
        self.boundary_nodes = self._read_data_file(bnode_file, dtype=float)
        print(f"  ✓ Loaded {len(self.boundary_nodes)} boundary nodes")
        
        # Detect geometry type and compute radii (verification)
        self._detect_geometry()

    def _transform_coordinates(self):
        """Transform logical (x,y) to physical (X,Y) for annular domains"""
        u = self.nodes[:, 1].copy()
        v = self.nodes[:, 2].copy()
        
        xmin, xmax = self.domain_params['xmin'], self.domain_params['xmax']
        ymin, ymax = self.domain_params['ymin'], self.domain_params['ymax']
        
        # Map u -> r
        # r = Rin + (u - xmin) / (xmax - xmin) * (Rout - Rin)
        if xmax > xmin:
            r = self.Rin + (u - xmin) / (xmax - xmin) * (self.Rout - self.Rin)
        else:
            r = u  # Fallback
            
        # Map v -> theta
        # theta = (v - ymin) / (ymax - ymin) * 2pi
        if ymax > ymin:
            theta = (v - ymin) / (ymax - ymin) * 2 * np.pi
        else:
            theta = v # Fallback
            
        # Polar to Cartesian
        X = r * np.cos(theta)
        Y = r * np.sin(theta)
        
        self.nodes[:, 1] = X
        self.nodes[:, 2] = Y

    def _read_data_file(self, filepath, dtype=float):
        """
        Robustly read data file, skipping comments and singleton headers
        """
        valid_lines = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Check if line contains data or just a header count
                # We skip lines with only 1 value (assuming counts/metadata)
                parts = line.split()
                if len(parts) > 1:
                    valid_lines.append(line)
        
        if not valid_lines:
            raise ValueError(f"No valid data lines found in {filepath}")
            
        return np.loadtxt(StringIO('\n'.join(valid_lines)), dtype=dtype)
        
    def _detect_geometry(self):
        """Detect if mesh is annular and compute inner/outer radii"""
        
        # Compute radii from node coordinates
        x = self.nodes[:, 1]
        y = self.nodes[:, 2]
        r = np.sqrt(x**2 + y**2)
        
        self.Rin = r.min()
        self.Rout = r.max()
        
        # Check if it's annular (circular distribution)
        theta = np.arctan2(y, x)
        theta_range = theta.max() - theta.min()
        
        is_annular = theta_range > 5.0  # More than ~286 degrees
        
        if is_annular:
            print(f"  ✓ Detected annular geometry: Rin={self.Rin:.3f}, Rout={self.Rout:.3f}")
        else:
            print(f"  ✓ Detected rectangular geometry: X∈[{x.min():.3f}, {x.max():.3f}], "
                  f"Y∈[{y.min():.3f}, {y.max():.3f}]")
    
    def plot_mesh_overview(self, output_file):
        """Create comprehensive mesh overview plot"""
        
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
        
        # Plot 1: Full mesh with elements
        ax1 = fig.add_subplot(gs[0, 0])
        self._plot_full_mesh(ax1, show_elements=True, show_nodes=False)
        ax1.set_title('Full Mesh with Elements', fontweight='bold')
        
        # Plot 2: Mesh with nodes highlighted
        ax2 = fig.add_subplot(gs[0, 1])
        self._plot_full_mesh(ax2, show_elements=True, show_nodes=True)
        ax2.set_title('Mesh with Node Points', fontweight='bold')
        
        # Plot 3: Boundary nodes only
        ax3 = fig.add_subplot(gs[0, 2])
        self._plot_boundary_nodes(ax3)
        ax3.set_title('Boundary Nodes', fontweight='bold')
        
        # Plot 4: Element quality (aspect ratios)
        ax4 = fig.add_subplot(gs[1, 0])
        self._plot_element_quality(ax4)
        ax4.set_title('Element Quality (Aspect Ratio)', fontweight='bold')
        
        # Plot 5: Radial distribution
        ax5 = fig.add_subplot(gs[1, 1])
        self._plot_radial_distribution(ax5)
        ax5.set_title('Radial Node Distribution', fontweight='bold')
        
        # Plot 6: Mesh statistics
        ax6 = fig.add_subplot(gs[1, 2])
        self._plot_statistics(ax6)
        ax6.axis('off')
        
        plt.suptitle(f'Mesh Visualization - {self.node_elem}-Node Elements', 
                     fontsize=16, fontweight='bold')
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved overview plot: {output_file}")
        plt.close()
    
    def _plot_full_mesh(self, ax, show_elements=True, show_nodes=False):
        """Plot full mesh structure"""
        
        x = self.nodes[:, 1]
        y = self.nodes[:, 2]
        
        if show_elements:
            # Draw element edges
            for elem in self.elements:
                elem_nodes = elem[1:]  # Skip element ID
                
                if self.node_elem == 4:
                    # Q1: Connect corners in order 0-1-2-3-0
                    corners = [0, 1, 2, 3, 0]
                    for i in range(4):
                        n1 = elem_nodes[corners[i]]
                        n2 = elem_nodes[corners[i+1]]
                        ax.plot([x[n1], x[n2]], [y[n1], y[n2]], 
                               'b-', linewidth=0.3, alpha=0.6)
                
                elif self.node_elem == 9:
                    # Q2: Connect corners 0-2-8-6-0
                    corners = [0, 2, 8, 6, 0]
                    for i in range(4):
                        n1 = elem_nodes[corners[i]]
                        n2 = elem_nodes[corners[i+1]]
                        ax.plot([x[n1], x[n2]], [y[n1], y[n2]], 
                               'b-', linewidth=0.3, alpha=0.6)
        
        if show_nodes:
            ax.plot(x, y, 'k.', markersize=1, alpha=0.5, label='Nodes')
        
        # Add circles for inner/outer radius
        if self.Rin > 0:
            circle_inner = Circle((0, 0), self.Rin, fill=False, 
                                 edgecolor='red', linewidth=2, 
                                 linestyle='--', label=f'Rin={self.Rin:.2f}')
            ax.add_patch(circle_inner)
        
        circle_outer = Circle((0, 0), self.Rout, fill=False, 
                             edgecolor='green', linewidth=2, 
                             linestyle='--', label=f'Rout={self.Rout:.2f}')
        ax.add_patch(circle_outer)
        
        ax.set_aspect('equal')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')
        
        # Set axis limits with padding
        padding = 0.1 * self.Rout
        ax.set_xlim(-self.Rout - padding, self.Rout + padding)
        ax.set_ylim(-self.Rout - padding, self.Rout + padding)
    
    def _plot_boundary_nodes(self, ax):
        """Plot boundary nodes highlighted"""
        
        x = self.nodes[:, 1]
        y = self.nodes[:, 2]
        
        # Plot all nodes in light gray
        ax.plot(x, y, 'o', color='lightgray', markersize=2, alpha=0.3, label='Interior')
        
        # Highlight boundary nodes
        bnode_ids = self.boundary_nodes[:, 0].astype(int)
        bx = x[bnode_ids]
        by = y[bnode_ids]
        
        # Color by radius to show inner vs outer boundary
        br = np.sqrt(bx**2 + by**2)
        
        scatter = ax.scatter(bx, by, c=br, s=20, cmap='coolwarm', 
                           edgecolor='black', linewidth=0.5, 
                           label='Boundary', zorder=5)
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Radius', rotation=270, labelpad=15)
        
        # Add reference circles
        if self.Rin > 0:
            circle_inner = Circle((0, 0), self.Rin, fill=False, 
                                 edgecolor='red', linewidth=1.5, linestyle='--')
            ax.add_patch(circle_inner)
        
        circle_outer = Circle((0, 0), self.Rout, fill=False, 
                             edgecolor='green', linewidth=1.5, linestyle='--')
        ax.add_patch(circle_outer)
        
        ax.set_aspect('equal')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')
        
        # Set axis limits
        padding = 0.1 * self.Rout
        ax.set_xlim(-self.Rout - padding, self.Rout + padding)
        ax.set_ylim(-self.Rout - padding, self.Rout + padding)
    
    def _plot_element_quality(self, ax):
        """Plot element quality (aspect ratios)"""
        
        aspect_ratios = []
        centroids_x = []
        centroids_y = []
        
        x = self.nodes[:, 1]
        y = self.nodes[:, 2]
        
        for elem in self.elements:
            elem_nodes = elem[1:]
            
            if self.node_elem == 4:
                corners = [0, 1, 2, 3]
            elif self.node_elem == 9:
                corners = [0, 2, 8, 6]
            else:
                continue
            
            # Get corner coordinates
            ex = x[elem_nodes[corners]]
            ey = y[elem_nodes[corners]]
            
            # Compute edge lengths
            dx1 = ex[1] - ex[0]
            dy1 = ey[1] - ey[0]
            dx2 = ex[3] - ex[0]
            dy2 = ey[3] - ey[0]
            
            len1 = np.sqrt(dx1**2 + dy1**2)
            len2 = np.sqrt(dx2**2 + dy2**2)
            
            # Aspect ratio (max/min)
            if len1 > 0 and len2 > 0:
                aspect = max(len1, len2) / min(len1, len2)
            else:
                aspect = 1.0
            
            aspect_ratios.append(aspect)
            centroids_x.append(ex.mean())
            centroids_y.append(ey.mean())
        
        aspect_ratios = np.array(aspect_ratios)
        
        # Plot elements colored by aspect ratio
        scatter = ax.scatter(centroids_x, centroids_y, c=aspect_ratios, 
                           s=10, cmap='RdYlGn_r', vmin=0.8, vmax=2.0,
                           edgecolor='none', alpha=0.8)
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Aspect Ratio (ideal=1.0)', rotation=270, labelpad=15)
        
        # Add reference circles
        if self.Rin > 0:
            circle_inner = Circle((0, 0), self.Rin, fill=False, 
                                 edgecolor='black', linewidth=1, linestyle='--')
            ax.add_patch(circle_inner)
        
        circle_outer = Circle((0, 0), self.Rout, fill=False, 
                             edgecolor='black', linewidth=1, linestyle='--')
        ax.add_patch(circle_outer)
        
        ax.set_aspect('equal')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.grid(True, alpha=0.3)
        
        # Add quality statistics
        mean_aspect = aspect_ratios.mean()
        max_aspect = aspect_ratios.max()
        good_elements = np.sum(aspect_ratios < 1.5)
        total_elements = len(aspect_ratios)
        
        stats_text = (f"Mean: {mean_aspect:.2f}\n"
                     f"Max: {max_aspect:.2f}\n"
                     f"Good (<1.5): {good_elements}/{total_elements} "
                     f"({100*good_elements/total_elements:.1f}%)")
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        padding = 0.1 * self.Rout
        ax.set_xlim(-self.Rout - padding, self.Rout + padding)
        ax.set_ylim(-self.Rout - padding, self.Rout + padding)
    
    def _plot_radial_distribution(self, ax):
        """Plot node distribution in radial direction"""
        
        x = self.nodes[:, 1]
        y = self.nodes[:, 2]
        r = np.sqrt(x**2 + y**2)
        
        # Histogram of radial positions
        ax.hist(r, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
        ax.set_xlabel('Radius')
        ax.set_ylabel('Number of Nodes')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add vertical lines for Rin and Rout
        if self.Rin > 0:
            ax.axvline(self.Rin, color='red', linestyle='--', 
                      linewidth=2, label=f'Rin={self.Rin:.2f}')
        ax.axvline(self.Rout, color='green', linestyle='--', 
                  linewidth=2, label=f'Rout={self.Rout:.2f}')
        
        ax.legend()
    
    def _plot_statistics(self, ax):
        """Display mesh statistics"""
        
        x = self.nodes[:, 1]
        y = self.nodes[:, 2]
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        
        # Compute statistics
        n_nodes = len(self.nodes)
        n_elements = len(self.elements)
        n_boundary = len(self.boundary_nodes)
        
        # Radial statistics
        r_min = r.min()
        r_max = r.max()
        r_mean = r.mean()
        
        # Angular coverage
        theta_min = theta.min() * 180 / np.pi
        theta_max = theta.max() * 180 / np.pi
        theta_range = theta_max - theta_min
        
        # Create statistics text
        stats_text = f"""
MESH STATISTICS
{'='*40}

Geometry:
  Type: {'Annular' if theta_range > 270 else 'Rectangular'}
  Inner radius: {r_min:.4f}
  Outer radius: {r_max:.4f}
  Mean radius: {r_mean:.4f}

Angular coverage:
  θ range: [{theta_min:.1f}°, {theta_max:.1f}°]
  Coverage: {theta_range:.1f}° ({theta_range/360*100:.1f}%)

Mesh size:
  Total nodes: {n_nodes:,}
  Total elements: {n_elements:,}
  Boundary nodes: {n_boundary:,}
  Element type: Q{1 if self.node_elem==4 else 2} ({self.node_elem}-node)

DOFs per node: 4 (Hermite)
Total DOFs: {n_nodes * 4:,}

Memory estimate: ~{n_nodes * 4 * 8 / 1024 / 1024:.1f} MB
"""
        
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
               fontsize=9, verticalalignment='top', family='monospace',
               bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    def plot_detailed_mesh(self, output_file):
        """Create detailed mesh plot showing node numbering (for small meshes)"""
        
        if len(self.nodes) > 500:
            print("  ⚠ Skipping detailed plot (too many nodes >500)")
            return
        
        fig, ax = plt.subplots(figsize=(12, 12))
        
        x = self.nodes[:, 1]
        y = self.nodes[:, 2]
        
        # Draw elements
        for elem in self.elements:
            elem_id = elem[0]
            elem_nodes = elem[1:]
            
            if self.node_elem == 4:
                corners = [0, 1, 2, 3, 0]
            elif self.node_elem == 9:
                corners = [0, 2, 8, 6, 0]
            else:
                continue
            
            for i in range(4):
                n1 = elem_nodes[corners[i]]
                n2 = elem_nodes[corners[i+1]]
                ax.plot([x[n1], x[n2]], [y[n1], y[n2]], 
                       'b-', linewidth=0.5, alpha=0.6)
            
            # Label element at centroid
            cx = x[elem_nodes].mean()
            cy = y[elem_nodes].mean()
            ax.text(cx, cy, str(elem_id), fontsize=6, ha='center', 
                   color='red', fontweight='bold')
        
        # Plot and label nodes
        ax.plot(x, y, 'ko', markersize=4)
        
        for i, (xi, yi) in enumerate(zip(x, y)):
            ax.text(xi, yi, str(i), fontsize=5, ha='right', va='bottom', 
                   color='blue')
        
        # Highlight boundary nodes
        bnode_ids = self.boundary_nodes[:, 0].astype(int)
        ax.plot(x[bnode_ids], y[bnode_ids], 'go', markersize=6, 
               fillstyle='none', linewidth=2, label='Boundary')
        
        ax.set_aspect('equal')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title(f'Detailed Mesh with Node/Element Numbering', 
                    fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        padding = 0.1 * self.Rout
        ax.set_xlim(-self.Rout - padding, self.Rout + padding)
        ax.set_ylim(-self.Rout - padding, self.Rout + padding)
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved detailed plot: {output_file}")
        plt.close()
    
    def plot_element_connectivity(self, output_file, max_elements=20):
        """Plot first few elements with node connectivity"""
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 14))
        axes = axes.flatten()
        
        x = self.nodes[:, 1]
        y = self.nodes[:, 2]
        
        n_plot = min(max_elements, len(self.elements))
        
        for idx in range(min(4, n_plot)):
            ax = axes[idx]
            elem = self.elements[idx]
            elem_id = elem[0]
            elem_nodes = elem[1:]
            
            # Plot this element's nodes
            ex = x[elem_nodes]
            ey = y[elem_nodes]
            
            # Draw element edges
            if self.node_elem == 4:
                corners = [0, 1, 2, 3, 0]
                ax.plot(ex[corners], ey[corners], 'b-', linewidth=2)
            elif self.node_elem == 9:
                corners = [0, 2, 8, 6, 0]
                ax.plot(ex[corners], ey[corners], 'b-', linewidth=2)
            
            # Plot and label all nodes
            for i, (xi, yi) in enumerate(zip(ex, ey)):
                ax.plot(xi, yi, 'ro', markersize=10)
                ax.text(xi, yi, f'{i}\n({elem_nodes[i]})', 
                       fontsize=8, ha='center', va='bottom',
                       bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
            
            ax.set_aspect('equal')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_title(f'Element {elem_id} Connectivity\n'
                        f'({self.node_elem}-node element)', 
                        fontweight='bold')
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for idx in range(min(4, n_plot), 4):
            axes[idx].axis('off')
        
        plt.suptitle('Element Connectivity Examples', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved connectivity plot: {output_file}")
        plt.close()


def main():
    """Main function to generate all mesh visualization plots"""
    
    print("="*60)
    print("FEM Mesh Visualization Tool")
    print("="*60)
    
    # Parse command-line arguments
    if len(sys.argv) > 1:
        mesh_dir = sys.argv[1]
    else:
        mesh_dir = 'mesh'
    
    if len(sys.argv) > 2:
        node_elem = int(sys.argv[2])
    else:
        # Auto-detect from mesh files
        node_elem = 9  # Default
    
    # Create output directory
    output_dir = Path('../output/figures')
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput directory: {output_dir}")
    
    # Initialize visualizer
    viz = MeshVisualizer(mesh_dir=mesh_dir, node_elem=node_elem)
    
    try:
        # Read mesh files
        viz.read_mesh_files()
        
        print("\nGenerating plots...")
        
        # Generate comprehensive overview
        viz.plot_mesh_overview(output_dir / 'mesh_overview.png')
        
        # Generate detailed plot (only for small meshes)
        viz.plot_detailed_mesh(output_dir / 'mesh_detailed.png')
        
        # Generate connectivity plot
        viz.plot_element_connectivity(output_dir / 'element_connectivity.png')
        
        print("\n" + "="*60)
        print("✓ All plots generated successfully!")
        print("="*60)
        print(f"\nOutput files:")
        print(f"  - {output_dir / 'mesh_overview.png'}")
        print(f"  - {output_dir / 'mesh_detailed.png'} (if mesh < 500 nodes)")
        print(f"  - {output_dir / 'element_connectivity.png'}")
        print()
        
    except FileNotFoundError as e:
        print(f"\n❌ Error: {e}")
        print("Make sure mesh files exist in the specified directory.")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()