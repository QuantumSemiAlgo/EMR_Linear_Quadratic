/**
 * @file mesh_generator.h
 * @brief Automatic Mesh Generation for FEM Analysis
 *
 * This module generates structured rectangular/annular meshes programmatically
 * from geometric specifications, eliminating the need for external mesh files.
 *
 * @author Auto-generated
 * @date 2026-01-23
 */

#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include "global_params_1.h"
#include <petsc.h>

/**
 * @brief Generates complete mesh file set for FEM analysis
 * @param dat Global data structure with geometry parameters
 * @param node_elem Number of nodes per element (4 or 9)
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode generate_mesh(data &dat, int node_elem);

/**
 * @brief Generates EMR-specific mesh with radial refinement and ports
 * @param dat Global data structure
 * @param R2_current Current inner radius R2
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode generate_emr_mesh(data &dat, double R2_current);

/**
 * @brief Generates node coordinate file (node.dat)
 * @param mesh_dir Output directory path
 * @param dat Global data structure
 * @param nodes_x Total nodes in x-direction
 * @param nodes_y Total nodes in y-direction
 * @param dx Node spacing in x
 * @param dy Node spacing in y
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode generate_nodes(const std::string &mesh_dir, const data &dat,
                              int nodes_x, int nodes_y, double dx, double dy);

/**
 * @brief Generates element connectivity file (elem.dat)
 * @param mesh_dir Output directory path
 * @param dat Global data structure
 * @param node_elem Nodes per element (4 or 9)
 * @param n_x Number of elements in x
 * @param n_y Number of elements in y
 * @param nodes_x Total nodes in x-direction
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode generate_elements(const std::string &mesh_dir, const data &dat,
                                 int node_elem, int n_x, int n_y, int nodes_x);

/**
 * @brief Generates boundary node file (bnode.dat)
 * @param mesh_dir Output directory path
 * @param dat Global data structure
 * @param nodes_x Total nodes in x-direction
 * @param nodes_y Total nodes in y-direction
 * @param dx Node spacing in x
 * @param dy Node spacing in y
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode generate_boundary_nodes(const std::string &mesh_dir,
                                       const data &dat, int nodes_x,
                                       int nodes_y, double dx, double dy);

/**
 * @brief Generates boundary element file (belem.dat)
 * @param mesh_dir Output directory path
 * @param n_x Number of elements in x
 * @param n_y Number of elements in y
 * @param nodes_x Total nodes in x-direction
 * @param node_elem Nodes per element (4 or 9)
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode generate_boundary_elements(const std::string &mesh_dir, int n_x,
                                          int n_y, int nodes_x, int node_elem);

/**
 * @brief Generates layer definition file (layer_def.dat)
 * @param mesh_dir Output directory path
 * @param dat Global data structure
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode generate_layer_definitions(const std::string &mesh_dir,
                                          const data &dat);

/**
 * @brief Generates layer mapping file (layer_map.dat)
 * @param mesh_dir Output directory path
 * @param dat Global data structure
 * @param total_elements Total number of elements
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode generate_layer_mapping(const std::string &mesh_dir,
                                      const data &dat, int total_elements);

/**
 * @brief Verifies the generated mesh for correctness (Diagnostic)
 * @param mesh_dir Output directory path
 * @param node_elem Nodes per element (4 or 9)
 * @param n_x Number of elements in x
 * @param n_y Number of elements in y
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode verify_mesh(const std::string &mesh_dir, int node_elem, int n_x,
                           int n_y);

#endif // MESH_GENERATOR_H
