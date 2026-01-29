/**
 * @file input_reader.cpp
 * @brief Input Parsing and Mesh Loading Module.
 *
 * This file contains functions to read the simulation configuration from an
 * input file (FEMstruct format) and to load the mesh data (nodes, elements,
 * boundary conditions) from external files. It handles parsing, error checking,
 * and populating the global data structure.
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#include "global_params_1.h"
#include "prototypes_1.h"
#include "utilrect_1.h"
#include <cstring>
#include <fstream>
#include <petsc.h>
#include <sstream>
#include <vector>

//------------------------------------------------------------------------------
// Helper Functions for Input Parsing
//------------------------------------------------------------------------------

/**
 * @brief Reads two double values from the next valid line of an input stream.
 *
 * Skips empty lines and comments (lines starting with '#').
 *
 * @param in Input file stream.
 * @param a Reference to store the first double value.
 * @param b Reference to store the second double value.
 * @return true If two doubles were successfully read.
 * @return false If end-of-file was reached or parsing failed.
 */
static bool read_pair(std::ifstream &in, double &a, double &b) {
  std::string line;
  while (std::getline(in, line)) {
    // Skip empty lines and comments
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream ss(line);
    // Attempt to extract two doubles
    return bool(ss >> a >> b);
  }
  return false;
}

/**
 * @brief Reads two integer values from the next valid line.
 *
 * @param in Input file stream.
 * @param i1 Reference to store the first integer.
 * @param i2 Reference to store the second integer.
 * @return true If two integers were successfully read.
 * @return false If EOF or failure.
 */
static bool read_pair_int(std::ifstream &in, int &i1, int &i2) {
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream ss(line);
    return bool(ss >> i1 >> i2);
  }
  return false;
}

/**
 * @brief Reads three integer values from the next valid line.
 *
 * @param in Input file stream.
 * @param i1 Reference to store the first integer.
 * @param i2 Reference to store the second integer.
 * @param i3 Reference to store the third integer.
 * @return true If three integers were successfully read.
 * @return false If EOF or failure.
 */
static bool read_triplet(std::ifstream &in, int &i1, int &i2, int &i3) {
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream ss(line);
    return bool(ss >> i1 >> i2 >> i3);
  }
  return false;
}

/**
 * @brief Reads a single integer count from the next valid line.
 *
 * @param in Input file stream.
 * @param out Reference to store the integer.
 * @return true If successful.
 * @return false If EOF or failure.
 */
static bool read_count(std::ifstream &in, int &out) {
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream ss(line);
    if (ss >> out)
      return true;
  }
  return false;
}

//------------------------------------------------------------------------------
// Main Input Reader Function
//------------------------------------------------------------------------------

/**
 * @brief Reads the main FEM configuration file.
 *
 * Parses the input file to set up domain extents, grid parameters, element
 * types, solver settings, and output paths. Populates the `data` structure.
 *
 * @param fname Path to the input file.
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success, error code otherwise.
 */
PetscErrorCode input_reader(const char *fname, data &dat) {
  PetscErrorCode ierr = 0;
  std::ifstream fs(fname);

  // Check if file opened successfully
  if (!fs.is_open()) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
            "Could not open FEMstruct input file");
  }

  // 1) Debug level
  // Controls the verbosity of the output.
  readin(fs, &dat.ndebug, "Debug level", 'i', "", dat.ndebug);

  // 2) Domain extents: Xmin/Xmax
  {
    double xmin, xmax;
    if (!read_pair(fs, xmin, xmax)) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
              "Failed to read Xmin/Xmax");
    }
    dat.xmin = xmin;
    dat.xmax = xmax;
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  Domain: Xmin=%.6f, Xmax=%.6f\n",
                       dat.xmin, dat.xmax);
    CHKERRQ(ierr);
  }

  // 3) Domain extents: Ymin/Ymax
  {
    double ymin, ymax;
    if (!read_pair(fs, ymin, ymax)) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
              "Failed to read Ymin/Ymax");
    }
    dat.ymin = ymin;
    dat.ymax = ymax;
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  Domain: Ymin=%.6f, Ymax=%.6f\n",
                       dat.ymin, dat.ymax);
    CHKERRQ(ierr);
  }

  // 4) Annular geometry flag
  readin(fs, &dat.use_annular, "Use annular geometry", 'i', "", dat.ndebug);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Annular geometry: %s\n",
                     dat.use_annular ? "YES" : "NO");
  CHKERRQ(ierr);

  // 5) Inner radius (only if annular)
  if (dat.use_annular) {
    readin(fs, &dat.Rin, "Inner radius", 'd', "", dat.ndebug);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  Inner radius: %.6f\n", dat.Rin);
    CHKERRQ(ierr);
  } else {
    dat.Rin = 0.0;
  }

  // 6) Outer radius (only if annular)
  if (dat.use_annular) {
    readin(fs, &dat.Rout, "Outer radius", 'd', "", dat.ndebug);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  Outer radius: %.6f\n", dat.Rout);
    CHKERRQ(ierr);
  } else {
    dat.Rout = 1.0;
  }

  // 6.5) Mesh generation parameters
  {
    int mesh_nx, mesh_ny, mesh_layers, auto_gen;
    if (!read_pair_int(fs, mesh_nx, mesh_ny)) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
              "Failed to read mesh n_x/n_y (expected 2 integers)");
    }
    dat.n_x = mesh_nx;
    dat.n_y = mesh_ny;

    if (!read_pair_int(fs, mesh_layers, auto_gen)) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
              "Failed to read n_layers/auto_generate_mesh");
    }
    dat.n_layers = mesh_layers;
    dat.auto_generate_mesh = auto_gen;

    ierr =
        PetscPrintf(PETSC_COMM_WORLD,
                    "  Mesh generation: n_x=%d, n_y=%d, layers=%d, auto=%d\n",
                    dat.n_x, dat.n_y, dat.n_layers, dat.auto_generate_mesh);
    CHKERRQ(ierr);

    // Read stretch factor (optional, separate line)
    // We assume it follows the layer/auto line if present, or defaults to 1.0
    double sfac = 1.0;
    // Actually, we should enforce strict reading order.
    // Let's modify the input file to include it.
    readin(fs, &dat.stretch_factor, "Stretch factor", 'd', "", dat.ndebug);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  Stretch factor: %.4f\n",
                       dat.stretch_factor);
    CHKERRQ(ierr);
  }

  // 7) Reconstruction grid: ndz_x, ndz_y
  // Defines the grid resolution for post-processing or reconstruction.
  {
    int nx, ny;
    if (!read_pair_int(fs, nx, ny)) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
              "Failed to read ndz_x/ndz_y");
    }
    dat.ndz_x = nx;
    dat.ndz_y = ny;
    dat.dzx = (dat.xmax - dat.xmin) / double(dat.ndz_x);
    dat.dzy = (dat.ymax - dat.ymin) / double(dat.ndz_y);
    ierr = PetscPrintf(
        PETSC_COMM_WORLD,
        "  Reconstruction grid: ndz_x=%d, ndz_y=%d → dzx=%.6f, dzy=%.6f\n",
        dat.ndz_x, dat.ndz_y, dat.dzx, dat.dzy);
    CHKERRQ(ierr);
  }

  // 8) Element topology & DOFs
  // node_elem: Nodes per element (e.g., 4 for Q1, 9 for Q9).
  // dof_per_node: Degrees of freedom per node.
  // Band_DOF: Bandwidth parameter (if applicable).
  {
    int ne, dpn, bd;
    if (!read_triplet(fs, ne, dpn, bd)) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
              "Failed to read node_elem, dof_per_node, Band_DOF");
    }
    dat.node_elem = ne;
    dat.dof_per_node = dpn;
    dat.nband = bd;
    dat.nele_mat = dat.node_elem * dat.dof_per_node; // Total DOFs per element
    dat.ndeg = dat.nele_mat;
    ierr = PetscPrintf(
        PETSC_COMM_WORLD,
        "  Elements: node_elem=%d, dof_per_node=%d, Band_DOF=%d\n"
        "    => nele_mat=%d, ndeg=%d\n",
        dat.node_elem, dat.dof_per_node, dat.nband, dat.nele_mat, dat.ndeg);
    CHKERRQ(ierr);
  }

  // 9) Quadrature Settings
  // Number of Gauss points for numerical integration.
  // Input file specifies points PER SIDE (1D rule).
  int ngaus_per_side;
  readin(fs, &ngaus_per_side, "Gauss pts per side", 'i', "", dat.ndebug);
  dat.ngaus = ngaus_per_side * ngaus_per_side;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Gauss pts: %d per side -> %d total\n",
                     ngaus_per_side, dat.ngaus);
  CHKERRQ(ierr);

  // Pre-calculate Gauss points and weights for the reference rectangle.
  dat.xigaus = new double[dat.ngaus];
  dat.etagaus = new double[dat.ngaus];
  dat.wgaus = new double[dat.ngaus];
  gaussRectangle(dat.ngaus, dat.xigaus, dat.etagaus, dat.wgaus);

  // 10) Solver Tolerances
  readin(fs, &dat.rtol, "Relative tol", 'd', "", dat.ndebug);
  readin(fs, &dat.maxit, "Max iterations", 'i', "", dat.ndebug);
  readin(fs, &dat.atol, "Absolute tol", 'd', "", dat.ndebug);
  readin(fs, &dat.divtol, "Divergence tol", 'd', "", dat.ndebug);

  // 11) Solver Type
  // e.g., "cg", "gmres", etc.
  {
    char name[PETSC_MAX_PATH_LEN];
    readin(fs, name, "Solver type", 'S', "", dat.ndebug);
    std::strncpy(dat.kspType, name, PETSC_MAX_PATH_LEN - 1);
    dat.kspType[PETSC_MAX_PATH_LEN - 1] = '\0';
  }

  // 12) Output Directory
  {
    char path[PETSC_MAX_PATH_LEN];
    readin(fs, path, "Output directory", 'S', "", dat.ndebug);
    std::memset(dat.outpath, 0, sizeof(dat.outpath));
    std::strncpy(dat.outpath, path, sizeof(dat.outpath) - 1);
  }

  fs.close();
  return ierr;
}

//------------------------------------------------------------------------------
// Mesh Input Function
//------------------------------------------------------------------------------

/**
 * @brief Loads mesh data from external files.
 *
 * Reads `node.dat`, `elem.dat`, `bnode.dat`, `belem.dat`, `layer_map.dat`,
 * and `layer_def.dat` from the appropriate mesh directory.
 *
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success, error code otherwise.
 */
PetscErrorCode mesh_input(data &dat) {
  PetscErrorCode ierr;

  // ---------------------------------------------------------------------------
  // Determine Mesh Directory
  // ---------------------------------------------------------------------------
  // Selects the source directory based on the element type (linear vs
  // quadratic).
  std::string base;
  if (dat.node_elem == 4) {
    base = "../mesh_4/"; // Linear (Q1) or Cubic Hermite
  } else if (dat.node_elem == 9) {
    base = "../mesh/"; // Quadratic (Q9) or Quintic Hermite
  } else {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
            "Unsupported element type");
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Reading mesh from: %s\n", base.c_str());
  CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Open Mesh Files
  // ---------------------------------------------------------------------------
  std::ifstream fe(base + "elem.dat"), fn(base + "node.dat"),
      fb(base + "bnode.dat"), fbE(base + "belem.dat");

  if (!fe.is_open() || !fn.is_open() || !fb.is_open() || !fbE.is_open()) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
            "Failed to open one or more mesh files in '../mesh' directory");
  }

  // ---------------------------------------------------------------------------
  // Read Elements (elem.dat)
  // ---------------------------------------------------------------------------
  if (!read_count(fe, dat.nelem)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
            "Could not read number of elements from elem.dat");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Elements: nelem=%d\n", dat.nelem);
  CHKERRQ(ierr);

  // Allocate memory for element connectivity and material IDs.
  dat.elem = new int *[dat.nelem];
  dat.material = new int[dat.nelem];

  int cnt = 0;
  std::string line;
  while (cnt < dat.nelem && std::getline(fe, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream ss(line);

    int id;
    ss >> id; // Read Element ID (unused, assuming sequential)

    dat.elem[cnt] = new int[dat.node_elem];
    // Read node indices for this element
    for (int a = 0; a < dat.node_elem; ++a) {
      ss >> dat.elem[cnt][a];
    }
    // Read material ID
    ss >> dat.material[cnt];
    ++cnt;
  }

  if (cnt < dat.nelem) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
            "elem.dat has fewer entries than declared nelem");
  }
  fe.close();

  // ---------------------------------------------------------------------------
  // Read Nodes (node.dat)
  // ---------------------------------------------------------------------------
  if (!read_count(fn, dat.ngnodes)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
            "Could not read number of nodes from node.dat");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Nodes: ngnodes=%d\n", dat.ngnodes);
  CHKERRQ(ierr);

  // Pass 1: Find the maximum node ID to size the array correctly.
  // This handles cases where node IDs are not contiguous or 0-indexed.
  int max_node_id = -1;
  std::streampos pos = fn.tellg(); // Save current file position
  while (std::getline(fn, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream ss(line);
    int id;
    ss >> id;
    if (id > max_node_id)
      max_node_id = id;
  }
  fn.clear();    // Clear EOF flag
  fn.seekg(pos); // Return to saved position

  ierr = PetscPrintf(PETSC_COMM_WORLD, "  Max node ID: %d\n", max_node_id);
  CHKERRQ(ierr);

  // Allocate node array based on max ID.
  dat.node = new double *[max_node_id + 1];
  for (int i = 0; i <= max_node_id; ++i) {
    dat.node[i] = new double[2];
    dat.node[i][0] = dat.node[i][1] = 0.0;
  }

  // Pass 2: Read coordinates.
  cnt = 0;
  while (cnt < dat.ngnodes && std::getline(fn, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream ss(line);
    int id;
    ss >> id;
    ss >> dat.node[id][0] >> dat.node[id][1];
    ++cnt;
  }

  if (cnt < dat.ngnodes) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
            "node.dat has fewer entries than declared ngnodes");
  }
  fn.close();

  // ---------------------------------------------------------------------------
  // Read Boundary Nodes (bnode.dat)
  // ---------------------------------------------------------------------------
  if (!read_count(fb, dat.nbnode) || !read_count(fb, dat.nboundary)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
            "Could not read boundary counts from bnode.dat");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Boundary: nbnode=%d, nboundary=%d\n",
                     dat.nbnode, dat.nboundary);
  CHKERRQ(ierr);

  dat.bnode = new double *[dat.nbnode];
  dat.bvalue = new double *[dat.nbnode];
  dat.btob = new int[dat.nbnode];
  dat.bnode_id = new int[dat.nbnode];

  cnt = 0;
  while (cnt < dat.nbnode && std::getline(fb, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream ss(line);

    int file_id, btag;
    dat.bnode[cnt] = new double[2];
    dat.bvalue[cnt] = new double[dat.dof_per_node];

    // Format: file_node_id  x  y  boundaryTag
    ss >> file_id >> dat.bnode[cnt][0] >> dat.bnode[cnt][1] >> btag;

    dat.bnode_id[cnt] = file_id;

    // Debug output for boundary nodes
    PetscPrintf(
        PETSC_COMM_WORLD,
        "  read bnode[%d]: file_id=%d  → node_id=%d  coords=(%.3f,%.3f)\n", cnt,
        file_id, dat.bnode_id[cnt], dat.bnode[cnt][0], dat.bnode[cnt][1]);

    dat.btob[cnt] = 0; // Default to Dirichlet (0)
    for (int d = 0; d < dat.dof_per_node; ++d)
      dat.bvalue[cnt][d] = 0.0; // Initialize boundary values to 0

    ++cnt;
  }
  if (cnt < dat.nbnode) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
            "bnode.dat has fewer entries than declared nbnode");
  }
  fb.close();

  // ---------------------------------------------------------------------------
  // Read Boundary Elements (belem.dat)
  // ---------------------------------------------------------------------------
  if (!read_count(fbE, dat.nbelem)) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
            "Could not read number of boundary elements from belem.dat");
  }
  ierr =
      PetscPrintf(PETSC_COMM_WORLD, "Boundary elems: nbelem=%d\n", dat.nbelem);
  CHKERRQ(ierr);

  dat.belem = new int *[dat.nbelem];
  dat.bmaterial = new int *[dat.nbelem];

  cnt = 0;
  while (cnt < dat.nbelem && std::getline(fbE, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::stringstream ss(line);
    int id;
    dat.belem[cnt] = new int[2];
    dat.bmaterial[cnt] = new int[2];

    // Format: ID n1 n2 mat1 mat2
    ss >> id >> dat.belem[cnt][0] >> dat.belem[cnt][1] >>
        dat.bmaterial[cnt][0] >> dat.bmaterial[cnt][1];
    ++cnt;
  }
  if (cnt < dat.nbelem) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
            "belem.dat has fewer entries than declared nbelem");
  }
  fbE.close();

  // ---------------------------------------------------------------------------
  // Read Layer Map (layer_map.dat)
  // ---------------------------------------------------------------------------
  // Maps each element to a specific layer ID.
  {
    std::ifstream flm(base + "layer_map.dat");
    if (!flm.is_open()) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
              "Could not open layer_map.dat");
    }
    dat.layer_id_for_elem = new int[dat.nelem];
    std::string line;
    int count = 0;
    while (std::getline(flm, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      std::istringstream ss(line);
      int eid, lid;
      ss >> eid >> lid;
      if (eid < 0 || eid >= dat.nelem) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_CORRUPT,
                "Invalid ElementID in layer_map.dat");
      }
      dat.layer_id_for_elem[eid] = lid;
      ++count;
    }
    if (count != dat.nelem) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
              "layer_map.dat entries don't match nelem");
    }
    flm.close();
    PetscPrintf(PETSC_COMM_WORLD, "Read layer_map.dat for %d elements\n",
                count);
  }

  // ---------------------------------------------------------------------------
  // Read Layer Definitions (layer_def.dat)
  // ---------------------------------------------------------------------------
  // Defines the physical properties (material, y-range) for each layer.
  {
    std::ifstream fld(base + "layer_def.dat");
    if (!fld.is_open()) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
              "Could not open layer_def.dat");
    }
    std::vector<int> lids, mats;
    std::vector<double> ymins, ymaxs;

    std::string line;
    while (std::getline(fld, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      std::istringstream ss(line);
      int lid, mat;
      double ylo, yhi;
      ss >> lid >> ylo >> yhi >> mat;
      lids.push_back(lid);
      ymins.push_back(ylo);
      ymaxs.push_back(yhi);
      mats.push_back(mat);
    }
    fld.close();

    dat.nlayer = (int)lids.size();
    dat.layer_id = new int[dat.nlayer];
    dat.layer_ymin = new double[dat.nlayer];
    dat.layer_ymax = new double[dat.nlayer];
    dat.layer_material = new int[dat.nlayer];

    for (int i = 0; i < dat.nlayer; ++i) {
      dat.layer_id[i] = lids[i];
      dat.layer_ymin[i] = ymins[i];
      dat.layer_ymax[i] = ymaxs[i];
      dat.layer_material[i] = mats[i];
    }
    PetscPrintf(PETSC_COMM_WORLD, "Read layer_def.dat: %d layers\n",
                dat.nlayer);
  }

  return 0;
}
