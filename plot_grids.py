import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import os

def read_nodes(filename):
    """
    Reads node.dat file.
    Returns a dictionary node_id -> [x, y]
    """
    nodes = {}
    with open(filename, 'r') as f:
        # Skip header comments
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) == 1:
                # Node count line
                continue
            if len(parts) >= 3:
                nid = int(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                nodes[nid] = [x, y]
    return nodes

def read_elements(filename):
    """
    Reads elem.dat file.
    Returns a list of elements, where each element is a list of node IDs.
    We only take the first 4 nodes to define the quad element boundary.
    """
    elements = []
    with open(filename, 'r') as f:
        # Skip header comments
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) == 1:
                # Element count line
                continue
            
            # ElemID node1 node2 ... material
            # We want node1, node2, node3, node4 (indices 1,2,3,4)
            if len(parts) >= 5:
                # 0-based indices in the file?
                # The file has: ElemID n1 n2 n3 n4 ...
                # Let's extract the first 4 node indices
                n1 = int(parts[1])
                n2 = int(parts[2])
                n3 = int(parts[3])
                n4 = int(parts[4])
                elements.append([n1, n2, n3, n4])
    return elements

def plot_grid(ax, folder_name, color='blue'):
    """
    Reads mesh files from folder and plots the grid on the given checks.
    """
    node_file = os.path.join(folder_name, 'node.dat')
    elem_file = os.path.join(folder_name, 'elem.dat')
    
    if not os.path.exists(node_file) or not os.path.exists(elem_file):
        print(f"Error: Files not found in {folder_name}")
        ax.text(0.5, 0.5, f"Files missing in {folder_name}", 
                ha='center', va='center', transform=ax.transAxes)
        return

    print(f"Reading {folder_name}...")
    nodes = read_nodes(node_file)
    elements = read_elements(elem_file)
    print(f"  Nodes: {len(nodes)}, Elements: {len(elements)}")

    if not elements:
        return

    # Create polygon vertices for each element
    verts = []
    for elem in elements:
        # Get coordinates for the 4 corners
        quad = []
        for nid in elem:
            if nid in nodes:
                quad.append(nodes[nid])
        if len(quad) == 4:
            verts.append(quad)
            
    # Create PolyCollection
    pc = PolyCollection(verts, facecolors='none', edgecolors=color, linewidths=0.5, alpha=0.5, label='Elements')
    ax.add_collection(pc)
    
    # Set limits based on nodes
    all_coords = np.array(list(nodes.values()))
    
    # Plot nodes
    # Extract x and y coordinates
    xs = [coord[0] for coord in all_coords]
    ys = [coord[1] for coord in all_coords]
    
    # Plot nodes with small markers
    ax.scatter(xs, ys, s=2, c='red', marker='.', zorder=10, label='Nodes')
    
    if len(all_coords) > 0:
        xmin, ymin = all_coords.min(axis=0)
        xmax, ymax = all_coords.max(axis=0)
        margin_x = (xmax - xmin) * 0.05
        margin_y = (ymax - ymin) * 0.05
        ax.set_xlim(xmin - margin_x, xmax + margin_x)
        ax.set_ylim(ymin - margin_y, ymax + margin_y)
        ax.set_aspect('equal')
        
    ax.set_title(f"Grid: {folder_name}", fontsize=12)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

def main():
    folder1 = "mesh_4"
    folder2 = "mesh"
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    plot_grid(axes[0], folder1, color='blue')
    plot_grid(axes[1], folder2, color='green')
    
    plt.tight_layout()
    output_file = '../output/figures/mesh_comparison.png'
    plt.savefig(output_file, dpi=300)
    print(f"\nSaved plot to {output_file}")
    plt.show()

if __name__ == "__main__":
    main()
