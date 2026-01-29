/**
 * @file utilrect.h
 * @brief Rectangular Element Utilities Header.
 *
 * This header declares utility functions for rectangular finite elements,
 * including quadrature rules, shape functions, element matrix assembly,
 * and geometric utilities.
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#pragma once

#include "global_params_1.h" // defines data, PetscScalar

// =============================================================================
// QUADRATURE
// =============================================================================

/**
 * @brief 1D Gauss-Legendre quadrature on [-1,1] for n = 2..6.
 * @param n Number of quadrature points.
 * @param xi Output array for quadrature points (length n).
 * @param w Output array for weights (length n).
 */
void gaussLegendre1D(int n, double *xi, double *w);

/**
 * @brief Legacy wrapper for shape function evaluation.
 * @param nfunc Number of shape functions (4 for Q1).
 * @param xi Local coordinate xi.
 * @param eta Local coordinate eta.
 * @param phi Output array for shape function values.
 */
void shapeRect(int nfunc, double xi, double eta, double *phi);

/**
 * @brief 2D tensor-product Gauss rule on [-1,1]^2.
 * @param ngaus Total number of points (must be perfect square: 4, 9, 16, 36).
 * @param xi Output array for xi coordinates (length ngaus).
 * @param eta Output array for eta coordinates (length ngaus).
 * @param w Output array for weights (length ngaus).
 */
void gaussRectangle(int ngaus, double *xi, double *eta, double *w);

// =============================================================================
// SHAPE FUNCTIONS
// =============================================================================

/**
 * @brief Evaluates all shape functions at (xi, eta).
 * @param N Output array (length = dat.node_elem * dat.dof_per_node).
 * @param xi Local coordinate xi in [-1, 1].
 * @param eta Local coordinate eta in [-1, 1].
 * @param dat Global data structure (determines element type).
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode shapeRect(double *N, double xi, double eta, const data &dat);

/**
 * @brief Evaluates shape function derivatives ∂N/∂ξ and ∂N/∂η.
 * @param dNdxi Output array for ∂N/∂ξ (length = dat.node_elem *
 * dat.dof_per_node).
 * @param dNdeta Output array for ∂N/∂η (length = dat.node_elem *
 * dat.dof_per_node).
 * @param xi Local coordinate xi.
 * @param eta Local coordinate eta.
 * @param dat Global data structure.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode deriv1Rect(double *dNdxi, double *dNdeta, double xi, double eta,
                          const data &dat);

// =============================================================================
// ELEMENT MATRIX ASSEMBLY
// =============================================================================

/**
 * @brief Assembles the local element stiffness matrix.
 * @param dat Global data structure.
 * @param ie Element index.
 * @param kappa Diffusion coefficient.
 * @param Ae Output stiffness matrix (flattened, size nloc*nloc).
 */
void element_stiffness(const data &dat, int ie, double kappa, PetscScalar *Ae);

/**
 * @brief Assembles the local element mass matrix.
 * @param dat Global data structure.
 * @param ie Element index.
 * @param Me Output mass matrix (flattened, size nloc*nloc).
 */
void element_mass(const data &dat, int ie, PetscScalar *Me);

// =============================================================================
// GEOMETRY UTILITIES
// =============================================================================

/**
 * @brief Finds the rectangular element containing point (x, y).
 * @param x X-coordinate of the point.
 * @param y Y-coordinate of the point.
 * @param dat Global data structure.
 * @param iel Output element index.
 * @return PetscErrorCode 0 on success, error if point is outside domain.
 */
PetscErrorCode locelem(double x, double y, const data &dat, int &iel);
