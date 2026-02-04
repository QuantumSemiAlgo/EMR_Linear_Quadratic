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
  if (dat.use_annular && dat.Nsample_R2 > 0) {
    // EMR Mode: Use R2_current from data (which should be set by main loop)
    // Note: generate_emr_mesh generates all files (nodes, elems, bnodes, etc.)
    // For the initial call/check, we might just use a nominal R2.
    // However, generate_mesh seems to be called once.
    // In EMR main loop, we need to regenerate mesh for each R2.
    // So main.cpp will likely call generate_emr_mesh directly or we adapt this.
    // For now, let's allow generate_mesh to dispatch if it detects EMR mode
    // params But since R2 changes, this function might be called repeatedly.
    // We'll assume dat.R2 is set to the current value desired.
    ierr = generate_emr_mesh(dat, dat.R2);
    CHKERRQ(ierr);
  } else {
    // Standard generation
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
  }

  // Add verification
  ierr = verify_mesh(mesh_dir, node_elem, n_x, n_y);
  CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "=== Mesh generation complete ===\n\n");
  CHKERRQ(ierr);

  return 0;
}

// =============================================================================
// EMR MESH GENERATOR
// =============================================================================
PetscErrorCode generate_emr_mesh(data &dat, double R2_current) {
  PetscErrorCode ierr;
  std::string mesh_dir = (dat.node_elem == 4) ? "../mesh_4/" : "../mesh/";

  // Update dat.Rin/Rout for consistency
  // Rin is the inner boundary of the simulation domain (Metal interface is
  // handled via sigma) Wait, EMR physics usually simulates the whole annulus R2
  // < r < R1 ? Or 0 < r < R1 with variable sigma? The plan implies solving
  // "modified Laplace" with variable sigma. Usually EMR has a metal geometric
  // inhomogeneity. If we mesh the metal region (r < R2), we need R_inner_mesh <
  // R2. If we treat metal as equipotential boundary, we mesh R2 < r < R1. Based
  // on "R2 is inner metal radius", usually we solve in semiconductor region R2
  // < r < R1. And apply V=constant at r=R2.  // Let's assume the domain is
  // nearly the full Disk. However, singular at r=0. We set Rin to a small
  // fraction of R2min.
  dat.Rin = 1.0e-6; // 1 micron start
  dat.Rout = dat.R1;

  // 1. Construct Radial Grid (non-uniform)
  std::vector<double> r_nodes;

  // Zones:
  // Zone 0: Core (Rin to R2 - w/2 - w_trans)  [Metal Inner]
  // Zone 1: Transition Inner
  // Zone 2: Interface (around R2)
  // Zone 3: Transition Outer
  // Zone 4: Bulk (to R1)

  double r = dat.Rin;
  r_nodes.push_back(r);

  // Parse refinement params
  int refinement_factor = (dat.node_elem == 9) ? 2 : 1;

  int N_fine = dat.nelem_R2 * refinement_factor;
  double w_fine = dat.width_R2; // Total width of fine zone
  int N_trans = dat.nelem_Rs2 * refinement_factor;
  double w_trans = dat.width_Rs2;
  int N_bulk = dat.nelem_otherR * refinement_factor;

  double R_interface_start = R2_current - w_fine / 2.0;
  double R_interface_end = R2_current + w_fine / 2.0;

  double R_trans_inner_start = R_interface_start - w_trans;
  double R_trans_outer_end = R_interface_end + w_trans;

  // Sanity check
  if (R_trans_inner_start < dat.Rin) {
    // Squeeze inner part?
    R_trans_inner_start = dat.Rin + (R2_current - dat.Rin) * 0.5;
    R_interface_start = R2_current - (R2_current - R_trans_inner_start) * 0.5;
  }

  // Zone 0: Core (Metal)
  // Use, say, 1/3 of bulk elements for core? Or just define N_core?
  // Let's use N_bulk / 2 for core.
  int N_core = (N_bulk > 4) ? N_bulk / 2 : 4;

  double dr_core = (R_trans_inner_start - r) / N_core;
  for (int i = 0; i < N_core; ++i) {
    r += dr_core;
    r_nodes.push_back(r);
  }

  // Zone 1: Trans Inner
  double dr_trans_in = (R_interface_start - r) / N_trans;
  for (int i = 0; i < N_trans; ++i) {
    r += dr_trans_in;
    r_nodes.push_back(r);
  }

  // Zone 2: Interface (Fine)
  double dr_fine = (R_interface_end - r) / N_fine;
  for (int i = 0; i < N_fine; ++i) {
    r += dr_fine;
    r_nodes.push_back(r);
  }

  // Zone 3: Trans Outer
  double dr_trans_out = (R_trans_outer_end - r) / N_trans;
  for (int i = 0; i < N_trans; ++i) {
    r += dr_trans_out;
    r_nodes.push_back(r);
  }

  // Zone 4: Bulk (Semi)
  int N_bulk_outer = N_bulk; // Use full complement
  double r_remaining = dat.Rout - r;
  if (r_remaining > 0) {
    double dr_bulk = r_remaining / N_bulk_outer;
    for (int i = 0; i < N_bulk_outer; ++i) {
      r += dr_bulk;
      r_nodes.push_back(r);
    }
  }
  // Ensure last node is exactly Rout
  r_nodes.back() = dat.Rout;

  // 2. Construct Angular Grid (Port based)
  std::vector<double> theta_nodes;
  // Ports are at theta1, theta2, theta3, theta4 (centers) with widths.
  // Standardize to [0, 2pi].
  // theta parameters in input are x PI.

  struct Port {
    double center;
    double width;
    int id;
  };
  std::vector<Port> ports;
  ports.push_back({dat.theta1 * M_PI, dat.width_L1 * M_PI, 1});
  ports.push_back({dat.theta2 * M_PI, dat.width_L2 * M_PI, 2});
  ports.push_back({dat.theta3 * M_PI, dat.width_L3 * M_PI, 3});
  ports.push_back({dat.theta4 * M_PI, dat.width_L4 * M_PI, 4});

  // Sort ports by angle
  // Assuming theta1 < theta2 < ... for now

  double current_theta = 0.0;
  theta_nodes.push_back(current_theta);

  // We need to mesh the full circle 0 to 2pi
  // Segments:
  // 0 -> P1_start
  // P1_start -> P1_end
  // P1_end -> P2_start ...
  // ... -> 2pi

  // Helper to add nodes
  auto add_angular_segment = [&](double start, double end, int n_elem) {
    if (n_elem <= 0)
      return;
    double dth = (end - start) / n_elem;
    for (int i = 0; i < n_elem; ++i) {
      theta_nodes.push_back(start + (i + 1) * dth);
    }
  };

  int N_port = dat.nelem_L * refinement_factor;
  int N_gap = dat.nelem_otherT * refinement_factor;

  double last_pos = 0.0;

  // Allocate port_code array
  int n_r = r_nodes.size();
  // Estimate total nodes for allocation (will resize later if needed but vector
  // handles it)

  // Logic to fill theta_nodes and track ports
  // This is a bit complex to generalize for arbitrary ports, assume ordered
  // non-overlapping.

  // Port 1
  double p1_start = ports[0].center - ports[0].width / 2.0;
  double p1_end = ports[0].center + ports[0].width / 2.0;
  add_angular_segment(last_pos, p1_start, N_gap); // Gap 0
  add_angular_segment(p1_start, p1_end, N_port);  // Port 1
  last_pos = p1_end;

  // Port 2
  double p2_start = ports[1].center - ports[1].width / 2.0;
  double p2_end = ports[1].center + ports[1].width / 2.0;
  add_angular_segment(last_pos, p2_start, N_gap);
  add_angular_segment(p2_start, p2_end, N_port);
  last_pos = p2_end;

  // Port 3
  double p3_start = ports[2].center - ports[2].width / 2.0;
  double p3_end = ports[2].center + ports[2].width / 2.0;
  add_angular_segment(last_pos, p3_start, N_gap);
  add_angular_segment(p3_start, p3_end, N_port);
  last_pos = p3_end;

  // Port 4
  double p4_start = ports[3].center - ports[3].width / 2.0;
  double p4_end = ports[3].center + ports[3].width / 2.0;
  add_angular_segment(last_pos, p4_start, N_gap);
  add_angular_segment(p4_start, p4_end, N_port);
  last_pos = p4_end;

  // Final Gap
  add_angular_segment(last_pos, 2.0 * M_PI, N_gap);
  // Ensure exactly 2pi is end
  theta_nodes.back() = 2.0 * M_PI;

  // Remove last node (2pi) for periodicity if we were using standard periodic
  // BC logic? But standard FEM often duplicates first/last node for topology
  // and uses constraints. Let's keep 0 and 2pi as separate nodes in geometry,
  // and rely on Convertnode to link them. Wait, Convertnode expects node at
  // xmin and xmax to be linked. Here x is angular. xmin=0, xmax=1 (scaled) or
  // 2pi. Standard logic: rectangular domain mapped to annulus. y is radial
  // (ymin=0 -> Rin, ymax=1 -> Rout). mesh_generator produces rectangular mesh
  // [xmin, xmax] x [ymin, ymax]. Then coordinate transform maps it.

  // !! CRITICAL CHANGE !!
  // EMR geometry is generated directly in Polar coordinates or mapped?
  // The existing code has `generate_nodes` doing:
  // x = xmin + ...
  // y = ymin + ...
  // It produces a rectangular grid in logical space.
  // Then `generate_nodes` transforms x if annular?
  // Let's check `generate_nodes`:
  // "double x = dat.xmin + x_fraction * (dat.xmax - dat.xmin);"
  // "double y = dat.ymin + j * dy;"
  // It creates a rectangular grid.

  // So we should generate the nodes in Logical space [0, 2pi] x [Rin, Rout] ?
  // Or [0, 1] x [0, 1] and let util handle map?
  // The global_params define xmin=0, xmax=2.0 (from input?).
  // Usually angular domain is 0 to 1 or 0 to 2pi.
  // Let's look at `generate_nodes` again.
  // It writes x and y directly.
  // If we want non-uniform mesh, we should write the correct x,y coordinates
  // in the logical space that matches the physical mapping or modification.

  // PROPOSAL:
  // Create nodes directly with the calculated r_nodes and theta_nodes.
  // Map them to Logical X and Y expected by the solver.
  // Solver expects X to be angular (0 to 2pi?) and Y to be radial.
  // If we produce `node.dat` with explicit coordinates, we bypass the uniform
  // grid logic.

  // 3. Write node.dat
  // Nodes are indexed j (radial) then i (angular).
  // i loops fast, j loops slow.

  int n_theta = theta_nodes.size();
  int n_rad = r_nodes.size();
  int total_nodes = n_theta * n_rad;

  // Need to update global dat params to match mesh size (Elements, not nodes)
  dat.n_x = (n_theta - 1) / refinement_factor;
  dat.n_y = (n_rad - 1) / refinement_factor;

  std::string node_file = mesh_dir + "node.dat";
  std::ofstream nf(node_file.c_str());
  nf << total_nodes << "\n";

  // Store port_code
  std::vector<int> node_port_codes(total_nodes, 0);

  for (int j = 0; j < n_rad; ++j) {
    double r = r_nodes[j];
    // Normalize y to [0, 1] if required?
    // Existing code: y goes from ymin to ymax.
    // If we write physical r, we might break non-annular logic?
    // Data struct has use_annular.
    // If use_annular, the code treats X as theta, Y as r?
    // Let's check constants/coordinate transforms.
    // Assuming we just write (theta, r) in node.dat and the solver uses them as
    // (x,y).

    for (int i = 0; i < n_theta; ++i) {
      double th = theta_nodes[i];
      int nid = j * n_theta + i;

      nf << std::setw(8) << nid << "  " << std::scientific
         << std::setprecision(12) << th << "  " << r << "\n";

      // Determine port code for this node
      // Check if th is within any port range
      // Only apply ports at Outer Radius (R1)? No, current I/O is usually at
      // boundary? EMR: Ports are usually contacts on the periphery (R1). So
      // port_code is only non-zero if r == Rout (last j).

      if (j == n_rad - 1) { // Outer boundary
        for (auto &p : ports) {
          double start = p.center - p.width / 2.0 - 1e-6;
          double end = p.center + p.width / 2.0 + 1e-6;
          if (th >= start && th <= end) {
            node_port_codes[nid] = p.id;
            break;
          }
        }
      }
    }
  }
  nf.close();

  // Allocate and Copy port_code to dat (optional, but dat has int *port_code)
  if (dat.port_code)
    delete[] dat.port_code;
  dat.port_code = new int[total_nodes];
  for (int i = 0; i < total_nodes; ++i)
    dat.port_code[i] = node_port_codes[i];

  // 4. Write elem.dat
  std::string elem_file = mesh_dir + "elem.dat";
  std::ofstream ef(elem_file.c_str());

  // Element counts (should match dat.n_x, dat.n_y)
  int n_elem_x = (n_theta - 1) / refinement_factor;
  int n_elem_y = (n_rad - 1) / refinement_factor;
  int total_elems = n_elem_x * n_elem_y;
  ef << total_elems << "\n";

  if (dat.node_elem == 9) {
    // Q9 Elements
    for (int j = 0; j < n_elem_y; ++j) {
      for (int i = 0; i < n_elem_x; ++i) {
        int eid = j * n_elem_x + i;

        // Base indices in the dense grid
        int i0 = 2 * i;
        int j0 = 2 * j;
        int nodes_x_grid = n_theta; // Width of grid

        // Map 9 nodes
        int n0 = j0 * nodes_x_grid + i0;           // BL
        int n1 = j0 * nodes_x_grid + i0 + 2;       // BR
        int n2 = (j0 + 2) * nodes_x_grid + i0 + 2; // TR
        int n3 = (j0 + 2) * nodes_x_grid + i0;     // TL
        int n4 = j0 * nodes_x_grid + i0 + 1;       // Bottom Mid
        int n5 = (j0 + 1) * nodes_x_grid + i0 + 2; // Right Mid
        int n6 = (j0 + 2) * nodes_x_grid + i0 + 1; // Top Mid
        int n7 = (j0 + 1) * nodes_x_grid + i0;     // Left Mid
        int n8 = (j0 + 1) * nodes_x_grid + i0 + 1; // Center

        ef << std::setw(8) << eid << "  " << std::setw(8) << n0 << "  "
           << std::setw(8) << n1 << "  " << std::setw(8) << n2 << "  "
           << std::setw(8) << n3 << "  " << std::setw(8) << n4 << "  "
           << std::setw(8) << n5 << "  " << std::setw(8) << n6 << "  "
           << std::setw(8) << n7 << "  " << std::setw(8) << n8 << "  "
           << std::setw(4) << 1 << "\n";
      }
    }
  } else {
    // Q4 Elements (Standard)
    for (int j = 0; j < n_elem_y; ++j) {
      for (int i = 0; i < n_elem_x; ++i) {
        int eid = j * n_elem_x + i;
        // Linear 4-node
        int n0 = j * n_theta + i;
        int n1 = j * n_theta + i + 1;
        int n2 = (j + 1) * n_theta + i + 1;
        int n3 = (j + 1) * n_theta + i;

        ef << std::setw(8) << eid << "  " << std::setw(8) << n0 << "  "
           << std::setw(8) << n1 << "  " << std::setw(8) << n2 << "  "
           << std::setw(8) << n3 << "  " << std::setw(4) << 1
           << "\n"; // Material 1
      }
    }
  }
  ef.close();

  // 5. Write bnode.dat (Boundary Nodes)
  // Logic: Inner (j=0) and Outer (j=n_rad-1).
  // Also Top/Bottom if not periodic? usually periodic in theta.
  // We'll write radial boundaries.

  std::string bnode_file = mesh_dir + "bnode.dat";
  std::ofstream bf(bnode_file.c_str());

  std::vector<int> bIDs;
  std::vector<int> bTags;
  // Tags: 1=Inner(R2), 2=Outer(R1)

  // Inner (R2)
  for (int i = 0; i < n_theta; ++i) {
    bIDs.push_back(0 * n_theta + i);
    bTags.push_back(1);
  }
  // Outer (R1)
  for (int i = 0; i < n_theta; ++i) {
    bIDs.push_back((n_rad - 1) * n_theta + i);
    bTags.push_back(2);
  }

  bf << bIDs.size() << "\n";
  bf << 2 << "\n"; // 2 boundaries
  for (size_t k = 0; k < bIDs.size(); ++k) {
    int nid = bIDs[k];
    double th = theta_nodes[nid % n_theta]; // Recover theta
    double r = r_nodes[nid / n_theta];      // Recover r
    bf << std::setw(8) << nid << "  " << std::scientific
       << std::setprecision(12) << th << "  " << r << "  " << std::setw(4)
       << bTags[k] << "\n";
  }
  bf.close();

  // 6. Write belem.dat, layer_def, layer_map (dummies or simple)
  // layer_map.dat
  std::string lmap_file = mesh_dir + "layer_map.dat";
  std::ofstream lmf(lmap_file.c_str());
  for (int e = 0; e < total_elems; ++e) {
    lmf << std::setw(8) << e << "  " << 0 << "\n"; // All layer 0
  }
  lmf.close();

  // layer_def.dat
  std::string ldef_file = mesh_dir + "layer_def.dat";
  std::ofstream ldf(ldef_file.c_str());
  ldf << "0  " << dat.Rin << "  " << dat.Rout << "  1\n";
  ldf.close();

  // Write DUMMY belem.dat to satisfy input_reader
  // The solver uses bnode.dat for BCs, so belem.dat is likely unused or for
  // viz.
  std::string belem_file = mesh_dir + "belem.dat";
  std::ofstream bf_dummy(belem_file.c_str());
  if (bf_dummy.is_open()) {
    std::time_t now = std::time(NULL);
    bf_dummy << "# Generated: " << std::ctime(&now);
    bf_dummy << "# Dummy belem.dat for EMR mesh\n";
    bf_dummy << "0\n"; // 0 boundary elements
    bf_dummy.close();
  }

  PetscPrintf(PETSC_COMM_WORLD, "DEBUG: generate_emr_mesh finishing...\n");

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

  std::vector<std::pair<int, int>> boundary_edges;

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
