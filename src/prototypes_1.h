/**
 * @file prototypes.h
 * @brief Function Prototypes.
 *
 * This header declares all the major functions used in the FEM application,
 * grouped by their functionality (Main, I/O, Quadrature, Matrix Assembly).
 * It ensures type safety and allows for modular compilation.
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include "global_params_1.h"
#include <fstream>
#include <petsc.h>

// Forward declare core structs to avoid circular includes
struct data;
struct global_matrices;

// =============================================================================
// MAIN APPLICATION ROUTINES
// =============================================================================

/**
 * @brief Reads the main input configuration file (FEMstruct).
 * @param fname Path to the input file.
 * @param dat Reference to the global data structure to populate.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode input_reader(const char *fname, data &dat);

/**
 * @brief Reads mesh data files (elem.dat, node.dat, etc.).
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode mesh_input(data &dat);

/**
 * @brief Assembles the global stiffness matrix (A) and RHS vector (b).
 * @param gmat Reference to the global matrices structure.
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode make_global(global_matrices &gmat, data &dat);

/**
 * @brief Applies boundary conditions (Dirichlet/Neumann) to the system.
 * @param gmat Reference to the global matrices structure.
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode apply_bc(global_matrices &gmat, data &dat);

/**
 * @brief Solves the linear system Ax=b and performs post-processing analysis.
 * @param gmat Reference to the global matrices structure.
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode petsc_diagonalizer(global_matrices &gmat, data &dat);

/**
 * @brief Interpolates the solution onto a regular grid and writes output.
 * @param gmat Reference to the global matrices structure.
 * @param dat Reference to the global data structure.
 * @param psi The solution vector.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode solution(global_matrices &gmat, data &dat, Vec psi);

/**
 * @brief Sets up the node overlay mapping for periodic boundary conditions.
 * @param dat Reference to the global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Convertnode(data &dat);

// =============================================================================
// I/O UTILITY ROUTINES
// =============================================================================

/**
 * @brief Reads a typed value from an input stream, handling comments and
 * errors.
 * @param file Input file stream.
 * @param data Pointer to the variable to store the read value.
 * @param description Description of the parameter (for error messages).
 * @param type Type character ('i'=int, 'd'=double, 's'=string, 'c'=char,
 * 'b'=bool).
 * @param units Units string (for display).
 * @param ndebug Debug level.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode readin(std::ifstream &file, void *data, const char *description,
                      char type, const char *units, int ndebug);

/**
 * @brief Skips characters in the stream until a delimiter is found.
 * @param file Input file stream.
 * @param c Delimiter character.
 */
void ignore_until(std::ifstream &file, char c);

/**
 * @brief Returns the current line number in the file stream (for error
 * reporting).
 * @param file Input file stream.
 * @return Current line number.
 */
int determine_current_line_number(std::ifstream &file);

// =============================================================================
// QUADRATURE & SHAPE ROUTINES
// =============================================================================

/**
 * @brief Generates 1D Gauss-Legendre quadrature points and weights.
 * @param n Number of points.
 * @param xi Output array for points.
 * @param w Output array for weights.
 */
void gaussLegendre1D(int n, double *xi, double *w);

/**
 * @brief Generates 2D tensor-product Gauss quadrature rule.
 * @param ngaus Total number of points.
 * @param xi Output array for xi coordinates.
 * @param eta Output array for eta coordinates.
 * @param w Output array for weights.
 */
void gaussRectangle(int ngaus, double *xi, double *eta, double *w);

// Specific shape function implementations
void shapeQ1(double N[4], double xi, double eta);
void deriv1Q1(double dNdxi[4], double dNdeta[4], double xi, double eta);
void shapeQ2(double N[9], double xi, double eta);
void deriv1Q2(double dNdxi[9], double dNdeta[9], double xi, double eta);
void shapeHermiteCubic(double N[16], double xi, double eta);
void deriv1HermiteCubic(double dNdxi[16], double dNdeta[16], double xi,
                        double eta);
void shapeHermiteQuintic(double N[36], double xi, double eta);
void deriv1HermiteQuintic(double dNdxi[36], double dNdeta[36], double xi,
                          double eta);

/**
 * @brief Dispatches to the correct shape function based on element type.
 * @param N Output array for shape function values.
 * @param xi Local coordinate xi.
 * @param eta Local coordinate eta.
 * @param dat Global data structure (defines element type).
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode shapeRect(double *N, double xi, double eta, const data &dat);

/**
 * @brief Dispatches to the correct shape function derivative based on element
 * type.
 * @param dNdxi Output array for dN/dxi.
 * @param dNdeta Output array for dN/deta.
 * @param xi Local coordinate xi.
 * @param eta Local coordinate eta.
 * @param dat Global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode deriv1Rect(double *dNdxi, double *dNdeta, double xi, double eta,
                          const data &dat);

// =============================================================================
// LOCAL ELEMENT MATRICES
// =============================================================================

/**
 * @brief Assembles the local element stiffness matrix.
 * @param dat Global data structure.
 * @param ie Element index.
 * @param kappa Diffusion coefficient.
 * @param Ae Output stiffness matrix (flattened).
 */
void element_stiffness(const data &dat, int ie, double kappa, PetscScalar *Ae);

/**
 * @brief Assembles the local element mass matrix.
 * @param dat Global data structure.
 * @param ie Element index.
 * @param Me Output mass matrix (flattened).
 */
void element_mass(const data &dat, int ie, PetscScalar *Me);

// =============================================================================
// GEOMETRY UTILITIES
// =============================================================================

/**
 * @brief Locates the element containing a specific point (x,y).
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @param dat Global data structure.
 * @param iel Output element index.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode locelem(double x, double y, const data &dat, int &iel);

#endif // PROTOTYPES_H