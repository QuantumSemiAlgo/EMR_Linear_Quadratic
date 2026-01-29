/**
 * @file mesh_generator.cpp
 * @brief Automatic Mesh Generation Implementation
 *
 * Generates all required mesh files for FEM analysis from geometric
 * specifications.
 *
 * @author Auto-generated
 * @date 2026-01-23
 */

#include "mesh_generator_1.h"
#include "constants_1.h"
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

// =============================================================================
// MAIN MESH GENERATION FUNCTION
// =============================================================================

PetscErrorCode generate_mesh(data &dat, int node_elem) {
  PetscErrorCode ierr;

  // Determine output directory based on element type
  std::string mesh_dir = (node_elem == 4) ? "../mesh_4/" : "../mesh/";

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n=== Generating mesh for %d-node elements ===\n",
                     node_elem);
  CHKERRQ(ierr);

  ierr =
      PetscPrintf(PETSC_COMM_WORLD, "Target directory: %s\n", mesh_dir.c_str());
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "Mesh parameters: n_x=%d, n_y=%d, annular=%d\n", dat.n_x,
                     dat.n_y, dat.use_annular);
  CHKERRQ(ierr);

  // FORCE DELETE old files
  std::vector<std::string> files;
  files.push_back("node.dat");
  files.push_back("elem.dat");
  files.push_back("bnode.dat");
  files.push_back("belem.dat");
  files.push_back("layer_map.dat");
  files.push_back("layer_def.dat");

  for (size_t i = 0; i < files.size(); ++i) {
    std::string fullpath = mesh_dir + files[i];
    std::remove(fullpath.c_str()); // Ignore errors if file doesn't exist
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Cleaned old mesh files in %s\n",
                     mesh_dir.c_str());
  CHKERRQ(ierr);

  // VALIDATION FOR ANNULAR GEOMETRY
  if (dat.use_annular) {
    // Check domain bounds are correct for periodicity
    if (std::abs(dat.ymin - 0.0) > 1e-10 || std::abs(dat.ymax - 1.0) > 1e-10) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "WARNING: For annular geometry, ymin should be 0.0 "
                         "and ymax should be 1.0\n");
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "         Current: ymin=%.6f, ymax=%.6f\n", dat.ymin,
                         dat.ymax);
      CHKERRQ(ierr);
    }

    if (dat.Rout <= dat.Rin) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
              "Outer radius must be greater than inner radius");
    }

    ierr =
        PetscPrintf(PETSC_COMM_WORLD, "Annular domain: Rin=%.6f, Rout=%.6f\n",
                    dat.Rin, dat.Rout);
    CHKERRQ(ierr);
  }

  // Create directory (portable way)
#ifdef _WIN32
  _mkdir(mesh_dir.c_str());
#else
  mkdir(mesh_dir.c_str(), 0755);
#endif

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "\n=== Generating mesh for %d-node elements in %s ===\n",
                     node_elem, mesh_dir.c_str());
  CHKERRQ(ierr);

  // Calculate mesh parameters
  int n_x = dat.n_x;
  int n_y = dat.n_y;

  int nodes_x, nodes_y, total_nodes, total_elements;

  if (node_elem == 4) {
    // Linear elements: (n_x+1) x (n_y+1) nodes
    nodes_x = n_x + 1;
    nodes_y = n_y + 1;
    total_nodes = nodes_x * nodes_y;
    total_elements = n_x * n_y;
  } else if (node_elem == 9) {
    // Quadratic elements: (2*n_x+1) x (2*n_y+1) nodes
    nodes_x = 2 * n_x + 1;
    nodes_y = 2 * n_y + 1;
    total_nodes = nodes_x * nodes_y;
    total_elements = n_x * n_y;
  } else {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
            "Unsupported element type: node_elem must be 4 or 9");
  }

  double dx = (dat.xmax - dat.xmin) / (nodes_x - 1);
  double dy = (dat.ymax - dat.ymin) / (nodes_y - 1);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Mesh size: %d x %d elements\n", n_x,
                     n_y);
  CHKERRQ(ierr);
  ierr =
      PetscPrintf(PETSC_COMM_WORLD, "  Total nodes: %d, Total elements: %d\n",
                  total_nodes, total_elements);
  CHKERRQ(ierr);

  // Generate all mesh files
  ierr = generate_nodes(mesh_dir, dat, nodes_x, nodes_y, dx, dy);
  CHKERRQ(ierr);

  ierr = generate_elements(mesh_dir, dat, node_elem, n_x, n_y, nodes_x);
  CHKERRQ(ierr);

  ierr = generate_boundary_nodes(mesh_dir, dat, nodes_x, nodes_y, dx, dy);
  CHKERRQ(ierr);

  ierr = generate_boundary_elements(mesh_dir, n_x, n_y, nodes_x, node_elem);
  CHKERRQ(ierr);

  ierr = generate_layer_definitions(mesh_dir, dat);
  CHKERRQ(ierr);

  ierr = generate_layer_mapping(mesh_dir, dat, total_elements);
  CHKERRQ(ierr);

  // Add verification
  ierr = verify_mesh(mesh_dir, node_elem, n_x, n_y);
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "=== Mesh generation complete ===\n\n");
  CHKERRQ(ierr);

  return 0;
}

// =============================================================================
// NODE GENERATION
// =============================================================================

PetscErrorCode generate_nodes(const std::string &mesh_dir, const data &dat,
                              int nodes_x, int nodes_y, double dx, double dy) {
  PetscErrorCode ierr;

  std::string filepath = mesh_dir + "node.dat";
  std::ofstream file(filepath);

  if (!file.is_open()) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Cannot create node.dat");
  }

  int total_nodes = nodes_x * nodes_y;

  // ADD TIMESTAMP HEADER
  std::time_t now = std::time(NULL);
  file << "# Generated: " << std::ctime(&now);
  file << "# Parameters: n_x=" << dat.n_x << ", n_y=" << dat.n_y
       << ", annular=" << dat.use_annular << "\n";
  file << "# Node count: " << total_nodes << "\n";

  file << total_nodes << "\n";

  // Generate nodes row by row (bottom to top)
  for (int j = 0; j < nodes_y; ++j) {
    for (int i = 0; i < nodes_x; ++i) {
      int node_id = j * nodes_x + i;

      double x_fraction = (double)i / (nodes_x - 1);
      if (dat.use_annular && dat.stretch_factor != 1.0) {
        x_fraction = std::pow(x_fraction, dat.stretch_factor);
      }

      double x = dat.xmin + x_fraction * (dat.xmax - dat.xmin);
      double y = dat.ymin + j * dy;

      file << std::setw(8) << node_id << "  " << std::scientific
           << std::setprecision(12) << x << "  " << y << "\n";
    }
  }

  file.close();

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Generated: %s (%d nodes)\n",
                     filepath.c_str(), total_nodes);
  CHKERRQ(ierr);

  return 0;
}

// =============================================================================
// ELEMENT GENERATION
// =============================================================================

PetscErrorCode generate_elements(const std::string &mesh_dir, const data &dat,
                                 int node_elem, int n_x, int n_y, int nodes_x) {
  PetscErrorCode ierr;

  std::string filepath = mesh_dir + "elem.dat";
  std::ofstream file(filepath);

  if (!file.is_open()) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Cannot create elem.dat");
  }

  // ADD TIMESTAMP HEADER
  std::time_t now = std::time(NULL);
  file << "# Generated: " << std::ctime(&now);
  file << "# Parameters: n_x=" << dat.n_x << ", n_y=" << dat.n_y
       << ", annular=" << dat.use_annular << "\n";

  int total_elements = n_x * n_y;
  file << total_elements << "\n";

  if (node_elem == 4) {
    // =========================================================================
    // 4-NODE LINEAR ELEMENTS
    // =========================================================================
    for (int ey = 0; ey < n_y; ++ey) {
      for (int ex = 0; ex < n_x; ++ex) {
        int elem_id = ey * n_x + ex;

        int node0 = ey * nodes_x + ex;           // Bottom-left
        int node1 = ey * nodes_x + ex + 1;       // Bottom-right
        int node2 = (ey + 1) * nodes_x + ex + 1; // Top-right
        int node3 = (ey + 1) * nodes_x + ex;     // Top-left

        int material_id = 1;

        file << std::setw(8) << elem_id << "  " << std::setw(8) << node0 << "  "
             << std::setw(8) << node1 << "  " << std::setw(8) << node2 << "  "
             << std::setw(8) << node3 << "  " << std::setw(4) << material_id
             << "\n";
      }
    }

  } else if (node_elem == 9) {
    // =========================================================================
    // 9-NODE QUADRATIC ELEMENTS
    // =========================================================================
    // CRITICAL: Match the node ordering expected by shapeQ2() in utilrect.cpp
    //
    // Expected ordering (from utilrect.cpp shapeQ2):
    // N[0] = M[0] = corner (-1,-1) = bottom-left
    // N[1] = M[6] = corner (+1,-1) = bottom-right
    // N[2] = M[8] = corner (+1,+1) = top-right
    // N[3] = M[2] = corner (-1,+1) = top-left
    // N[4] = M[3] = midside (0,-1) = bottom-mid
    // N[5] = M[7] = midside (+1,0) = right-mid
    // N[6] = M[5] = midside (0,+1) = top-mid
    // N[7] = M[1] = midside (-1,0) = left-mid
    // N[8] = M[4] = center (0,0)

    for (int ey = 0; ey < n_y; ++ey) {
      for (int ex = 0; ex < n_x; ++ex) {
        int elem_id = ey * n_x + ex;

        // Base node index at bottom-left corner of element
        int i0 = 2 * ex;
        int j0 = 2 * ey;

        // Calculate all 9 node indices
        int node0 = j0 * nodes_x + i0;           // Bottom-left corner
        int node1 = j0 * nodes_x + i0 + 2;       // Bottom-right corner
        int node2 = (j0 + 2) * nodes_x + i0 + 2; // Top-right corner
        int node3 = (j0 + 2) * nodes_x + i0;     // Top-left corner
        int node4 = j0 * nodes_x + i0 + 1;       // Bottom mid-edge
        int node5 = (j0 + 1) * nodes_x + i0 + 2; // Right mid-edge
        int node6 = (j0 + 2) * nodes_x + i0 + 1; // Top mid-edge
        int node7 = (j0 + 1) * nodes_x + i0;     // Left mid-edge
        int node8 = (j0 + 1) * nodes_x + i0 + 1; // Center

        int material_id = 1;

        file << std::setw(8) << elem_id << "  " << std::setw(8) << node0 << "  "
             << std::setw(8) << node1 << "  " << std::setw(8) << node2 << "  "
             << std::setw(8) << node3 << "  " << std::setw(8) << node4 << "  "
             << std::setw(8) << node5 << "  " << std::setw(8) << node6 << "  "
             << std::setw(8) << node7 << "  " << std::setw(8) << node8 << "  "
             << std::setw(4) << material_id << "\n";
      }
    }
  }

  file.close();

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Generated: %s (%d elements)\n",
                     filepath.c_str(), total_elements);
  CHKERRQ(ierr);

  return 0;
}

// =============================================================================
// BOUNDARY NODE GENERATION
// =============================================================================

PetscErrorCode generate_boundary_nodes(const std::string &mesh_dir,
                                       const data &dat, int nodes_x,
                                       int nodes_y, double dx, double dy) {
  PetscErrorCode ierr;

  std::string filepath = mesh_dir + "bnode.dat";
  std::ofstream file(filepath);

  if (!file.is_open()) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Cannot create bnode.dat");
  }

  // ADD TIMESTAMP HEADER
  std::time_t now = std::time(NULL);
  file << "# Generated: " << std::ctime(&now);
  file << "# Parameters: n_x=" << dat.n_x << ", n_y=" << dat.n_y
       << ", annular=" << dat.use_annular << "\n";

  std::vector<int> bnode_ids;
  std::vector<double> bnode_x;
  std::vector<double> bnode_y;
  std::vector<int> bnode_tags; // 1=Left/Inner, 2=Right/Outer, 3=Bottom, 4=Top

  if (dat.use_annular) {
    // ================================================================
    // ANNULAR GEOMETRY: Only radial boundaries (inner + outer circles)
    // Top/bottom are PERIODIC and handled by node overlay
    // ================================================================

    ierr =
        PetscPrintf(PETSC_COMM_WORLD,
                    "  Generating annular boundary nodes (radial only)...\n");
    CHKERRQ(ierr);

    // Inner circle: i = 0 (all j except possibly last if periodic)
    // For n_y elements, we have nodes_y nodes vertically.
    // If standard periodicity is handled later, we list all nodes on the
    // physical boundary. However, usually for periodic BCs, we might exclude
    // the last node if it duplicates the first. Reviewing 'convertnode.cpp', it
    // seems to handle finding duplicates. Ideally, we list all geometric
    // boundary nodes.

    // Inner boundary (Left in rectangular map -> Inner radius)
    for (int j = 0; j < nodes_y; ++j) {
      // Logic: The "bottom" and "top" (theta=0 and theta=2pi) are periodic
      // boundaries. We generally DO apply Dirichlet BCs on the radial
      // boundaries for ALL theta. So we include all j.
      int node_id = j * nodes_x + 0;
      double x = dat.xmin + 0 * dx; // i=0
      double y = dat.ymin + j * dy;

      bnode_ids.push_back(node_id);
      bnode_x.push_back(x);
      bnode_y.push_back(y);
      bnode_tags.push_back(1); // Tag 1 for Inner
    }

    // Outer circle: i = nodes_x-1
    for (int j = 0; j < nodes_y; ++j) {
      int node_id = j * nodes_x + (nodes_x - 1);
      double x = dat.xmin + (nodes_x - 1) * dx; // i=nodes_x-1
      double y = dat.ymin + j * dy;

      bnode_ids.push_back(node_id);
      bnode_x.push_back(x);
      bnode_y.push_back(y);
      bnode_tags.push_back(2); // Tag 2 for Outer
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "  Created %d boundary nodes on radial edges\n",
                       (int)bnode_ids.size());
    CHKERRQ(ierr);

  } else {
    // ================================================================
    // RECTANGULAR GEOMETRY: All four edges
    // ================================================================

    ierr =
        PetscPrintf(PETSC_COMM_WORLD,
                    "  Generating rectangular boundary nodes (all edges)...\n");
    CHKERRQ(ierr);

    std::vector<bool> is_boundary(nodes_x * nodes_y, false);

    // Bottom edge (j=0, all i)
    for (int i = 0; i < nodes_x; ++i) {
      int node_id = i;
      if (!is_boundary[node_id]) {
        bnode_ids.push_back(node_id);
        bnode_x.push_back(dat.xmin + i * dx);
        bnode_y.push_back(dat.ymin);
        bnode_tags.push_back(3); // Tag 3 for Bottom
        is_boundary[node_id] = true;
      }
    }

    // Right edge (i=nodes_x-1)
    for (int j = 0; j < nodes_y; ++j) {
      int node_id = j * nodes_x + (nodes_x - 1);
      if (!is_boundary[node_id]) {
        bnode_ids.push_back(node_id);
        bnode_x.push_back(dat.xmin + (nodes_x - 1) * dx);
        bnode_y.push_back(dat.ymin + j * dy);
        bnode_tags.push_back(2); // Tag 2 for Right
        is_boundary[node_id] = true;
      }
    }

    // Top edge (j=nodes_y-1)
    for (int i = 0; i < nodes_x; ++i) {
      int node_id = (nodes_y - 1) * nodes_x + i;
      if (!is_boundary[node_id]) {
        bnode_ids.push_back(node_id);
        bnode_x.push_back(dat.xmin + i * dx);
        bnode_y.push_back(dat.ymin + (nodes_y - 1) * dy);
        bnode_tags.push_back(4); // Tag 4 for Top
        is_boundary[node_id] = true;
      }
    }

    // Left edge (i=0)
    for (int j = 0; j < nodes_y; ++j) {
      int node_id = j * nodes_x;
      if (!is_boundary[node_id]) {
        bnode_ids.push_back(node_id);
        bnode_x.push_back(dat.xmin);
        bnode_y.push_back(dat.ymin + j * dy);
        bnode_tags.push_back(1); // Tag 1 for Left
        is_boundary[node_id] = true;
      }
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "  Created %d boundary nodes (rectangular)\n",
                       (int)bnode_ids.size());
    CHKERRQ(ierr);
  }

  // Write to file
  // Format: Count \n NumberOfBoundaries \n (NodeID X Y Tag)
  file << bnode_ids.size() << "\n";
  int nboundary = dat.use_annular ? 2 : 4;
  file << nboundary << "\n";

  for (size_t k = 0; k < bnode_ids.size(); ++k) {
    file << std::setw(8) << bnode_ids[k] << "  " << std::scientific
         << std::setprecision(12) << bnode_x[k] << "  " << bnode_y[k] << "  "
         << std::setw(4) << bnode_tags[k] << "\n";
  }

  file.close();

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Generated: %s\n", filepath.c_str());
  CHKERRQ(ierr);

  return 0;
}

// =============================================================================
// BOUNDARY ELEMENT GENERATION
// =============================================================================

PetscErrorCode generate_boundary_elements(const std::string &mesh_dir, int n_x,
                                          int n_y, int nodes_x, int node_elem) {
  PetscErrorCode ierr;

  std::string filepath = mesh_dir + "belem.dat";
  std::ofstream file(filepath);

  if (!file.is_open()) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Cannot create belem.dat");
  }

  // ADD TIMESTAMP HEADER
  std::time_t now = std::time(NULL);
  file << "# Generated: " << std::ctime(&now);
  file << "# Parameters: n_x=" << n_x << ", n_y=" << n_y
       << ", node_elem=" << node_elem << "\n";

  std::vector<std::pair<int, int> > boundary_edges;

  int nodes_y = (node_elem == 4) ? (n_y + 1) : (2 * n_y + 1);
  int stride = (node_elem == 4) ? 1 : 2; // Node spacing along element edges

  // Bottom boundary (j=0)
  for (int i = 0; i < n_x; ++i) {
    int node1 = i * stride;
    int node2 = (i + 1) * stride;
    boundary_edges.push_back(std::make_pair(node1, node2));
  }

  // Right boundary (i=nodes_x-1)
  for (int j = 0; j < n_y; ++j) {
    int node1 = j * stride * nodes_x + (nodes_x - 1);
    int node2 = (j + 1) * stride * nodes_x + (nodes_x - 1);
    boundary_edges.push_back(std::make_pair(node1, node2));
  }

  // Top boundary (j=nodes_y-1) - reverse order
  int top_j = nodes_y - 1;
  for (int i = n_x - 1; i >= 0; --i) {
    int node1 = top_j * nodes_x + (i + 1) * stride;
    int node2 = top_j * nodes_x + i * stride;
    boundary_edges.push_back(std::make_pair(node1, node2));
  }

  // Left boundary (i=0) - reverse order
  for (int j = n_y - 1; j >= 0; --j) {
    int node1 = (j + 1) * stride * nodes_x;
    int node2 = j * stride * nodes_x;
    boundary_edges.push_back(std::make_pair(node1, node2));
  }

  int nbelem = boundary_edges.size();
  file << nbelem << "\n";

  for (size_t k = 0; k < boundary_edges.size(); ++k) {
    file << std::setw(8) << k << "  " << std::setw(8) << boundary_edges[k].first
         << "  " << std::setw(8) << boundary_edges[k].second << "  "
         << std::setw(4) << 1 << "  "  // material left
         << std::setw(4) << 0 << "\n"; // material right
  }

  file.close();

  ierr =
      PetscPrintf(PETSC_COMM_WORLD, "  Generated: %s (%d boundary elements)\n",
                  filepath.c_str(), nbelem);
  CHKERRQ(ierr);

  return 0;
}

// =============================================================================
// LAYER DEFINITION GENERATION
// =============================================================================

PetscErrorCode generate_layer_definitions(const std::string &mesh_dir,
                                          const data &dat) {
  PetscErrorCode ierr;

  std::string filepath = mesh_dir + "layer_def.dat";
  std::ofstream file(filepath);

  if (!file.is_open()) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
            "Cannot create layer_def.dat");
  }

  // ADD TIMESTAMP HEADER
  std::time_t now = std::time(NULL);
  file << "# Generated: " << std::ctime(&now);
  file << "# Parameters: layers=" << dat.n_layers
       << ", annular=" << dat.use_annular << "\n";

  int n_layers = (dat.n_layers > 0) ? dat.n_layers : 1;
  double layer_height = (dat.ymax - dat.ymin) / n_layers;

  for (int L = 0; L < n_layers; ++L) {
    double y_min = dat.ymin + L * layer_height;
    double y_max = dat.ymin + (L + 1) * layer_height;
    int material_id = 1; // Default material

    file << std::setw(4) << L << "  " << std::scientific
         << std::setprecision(12) << y_min << "  " << y_max << "  "
         << std::setw(4) << material_id << "\n";
  }

  file.close();

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Generated: %s (%d layers)\n",
                     filepath.c_str(), n_layers);
  CHKERRQ(ierr);

  return 0;
}

// =============================================================================
// LAYER MAPPING GENERATION
// =============================================================================

PetscErrorCode generate_layer_mapping(const std::string &mesh_dir,
                                      const data &dat, int total_elements) {
  PetscErrorCode ierr;

  std::string filepath = mesh_dir + "layer_map.dat";
  std::ofstream file(filepath);

  if (!file.is_open()) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
            "Cannot create layer_map.dat");
  }

  // ADD TIMESTAMP HEADER
  std::time_t now = std::time(NULL);
  file << "# Generated: " << std::ctime(&now);
  file << "# Parameters: n_x=" << dat.n_x << ", n_y=" << dat.n_y
       << ", layers=" << dat.n_layers << "\n";

  int n_layers = (dat.n_layers > 0) ? dat.n_layers : 1;
  double layer_height = (dat.ymax - dat.ymin) / n_layers;

  int n_x = dat.n_x;
  int n_y = dat.n_y;
  double dy = (dat.ymax - dat.ymin) / n_y;

  for (int ey = 0; ey < n_y; ++ey) {
    for (int ex = 0; ex < n_x; ++ex) {
      int elem_id = ey * n_x + ex;

      // Element centroid y-coordinate
      double elem_y = dat.ymin + (ey + 0.5) * dy;

      // Determine which layer this element belongs to
      int layer_id = (int)((elem_y - dat.ymin) / layer_height);
      if (layer_id >= n_layers)
        layer_id = n_layers - 1;
      if (layer_id < 0)
        layer_id = 0;

      file << std::setw(8) << elem_id << "  " << std::setw(4) << layer_id
           << "\n";
    }
  }

  file.close();

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Generated: %s (%d elements mapped)\n",
                     filepath.c_str(), total_elements);
  CHKERRQ(ierr);

  return 0;
}

// =============================================================================
// DIAGNOSTIC: Verify Generated Mesh
// =============================================================================

PetscErrorCode verify_mesh(const std::string &mesh_dir, int node_elem, int n_x,
                           int n_y) {
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n=== Verifying mesh in %s ===\n",
                     mesh_dir.c_str());
  CHKERRQ(ierr);

  // Read and verify first element
  std::ifstream elem_file(mesh_dir + "elem.dat");
  if (!elem_file.is_open()) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  ERROR: Cannot open elem.dat\n");
    return 1;
  }

  int total_elems;
  elem_file >> total_elems;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Total elements: %d (expected: %d)\n",
                     total_elems, n_x * n_y);
  CHKERRQ(ierr);

  // Read first element
  int elem_id;
  std::vector<int> nodes(node_elem);
  int mat;

  elem_file >> elem_id;
  for (int i = 0; i < node_elem; ++i) {
    elem_file >> nodes[i];
  }
  elem_file >> mat;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  First element (ID=%d):\n", elem_id);
  CHKERRQ(ierr);

  if (node_elem == 4) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "    Corners: BL=%d, BR=%d, TR=%d, TL=%d\n", nodes[0],
                       nodes[1], nodes[2], nodes[3]);
  } else if (node_elem == 9) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "    Corners: BL=%d, BR=%d, TR=%d, TL=%d\n", nodes[0],
                       nodes[1], nodes[2], nodes[3]);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "    Midsides: B=%d, R=%d, T=%d, L=%d, Center=%d\n",
                       nodes[4], nodes[5], nodes[6], nodes[7], nodes[8]);
  }
  CHKERRQ(ierr);

  elem_file.close();

  // Verify expected values for first element
  if (node_elem == 9) {
    int nodes_x = 2 * n_x + 1;
    int expected[9] = {
        0,               // BL corner
        2,               // BR corner
        2 * nodes_x + 2, // TR corner
        2 * nodes_x,     // TL corner
        1,               // Bottom mid
        nodes_x + 2,     // Right mid
        2 * nodes_x + 1, // Top mid
        nodes_x,         // Left mid
        nodes_x + 1      // Center
    };

    bool ok = true;
    for (int i = 0; i < 9; ++i) {
      if (nodes[i] != expected[i]) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                           "  ERROR: Node %d = %d, expected %d\n", i, nodes[i],
                           expected[i]);
        CHKERRQ(ierr);
        ok = false;
      }
    }

    if (ok) {
      ierr =
          PetscPrintf(PETSC_COMM_WORLD, "  ✓ Element connectivity correct\n");
      CHKERRQ(ierr);
    }
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "=== Verification complete ===\n\n");
  CHKERRQ(ierr);

  return 0;
}
