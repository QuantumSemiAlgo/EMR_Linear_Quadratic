/**
 * @file utilrect.cpp
 * @brief Finite Element Utility Functions.
 *
 * This module contains the core mathematical functions for the Finite Element
 * Method (FEM) on 2D rectangular domains. It includes:
 * - Shape functions and derivatives for Q1 (linear), Q2 (quadratic), and
 * Hermite elements.
 * - Gauss-Legendre quadrature rules for numerical integration.
 * - Local element stiffness and mass matrix assembly.
 * - Coordinate transformations for annular domains (Jacobian scaling).
 *
 * @author Arya Retheeshan
 * @date 2025-11-24
 */

#include "utilrect_1.h"
#include "shape_funcs_1.hpp"
#include <cmath>
#include <petsc.h>
#include <vector>

// =============================================================================
// QUADRATURE RULES
// =============================================================================

/**
 * @brief Evaluates shape functions at a specific point (legacy wrapper).
 *
 * @param nfunc Number of shape functions (4 for Q1).
 * @param xi Local coordinate xi [-1, 1].
 * @param eta Local coordinate eta [-1, 1].
 * @param phi Output array for shape function values.
 */
void shapeRect(int nfunc, double xi, double eta, double *phi) {
  // Only support 4-node linear quads for this legacy function
  if (nfunc != 4) {
    for (int i = 0; i < nfunc; ++i)
      phi[i] = 0.0;
    return;
  }
  // Node ordering: BL, BR, TR, TL
  phi[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
  phi[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
  phi[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
  phi[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
}

/**
 * @brief Generates 1D Gauss-Legendre quadrature points and weights.
 *
 * @param n Number of points (order). Supports 2 to 6.
 * @param xi Output vector for quadrature points.
 * @param w Output vector for weights.
 */
void gaussLegendre1D(int n, std::vector<double> &xi, std::vector<double> &w) {
  xi.resize(n);
  w.resize(n);
  switch (n) {
  case 2:
    xi[0] = -0.5773502691896257;
    xi[1] = 0.5773502691896257;
    w[0] = 1.0;
    w[1] = 1.0;
    break;
  case 3:
    xi[0] = -0.7745966692414834;
    xi[1] = 0.0;
    xi[2] = 0.7745966692414834;
    w[0] = 0.5555555555555556;
    w[1] = 0.8888888888888888;
    w[2] = 0.5555555555555556;
    break;
  case 4:
    xi[0] = -0.8611363115940526;
    xi[1] = -0.3399810435848563;
    xi[2] = 0.3399810435848563;
    xi[3] = 0.8611363115940526;
    w[0] = 0.3478548451374539;
    w[1] = 0.6521451548625461;
    w[2] = 0.6521451548625461;
    w[3] = 0.3478548451374539;
    break;
  case 5:
    xi[0] = -0.9061798459386640;
    xi[1] = -0.5384693101056831;
    xi[2] = 0.0;
    xi[3] = 0.5384693101056831;
    xi[4] = 0.9061798459386640;
    w[0] = 0.2369268850561891;
    w[1] = 0.4786286704993665;
    w[2] = 0.5688888888888889;
    w[3] = 0.4786286704993665;
    w[4] = 0.2369268850561891;
    break;
  case 6:
    xi[0] = -0.9324695142031521;
    xi[1] = -0.6612093864662645;
    xi[2] = -0.2386191860831969;
    xi[3] = 0.2386191860831969;
    xi[4] = 0.6612093864662645;
    xi[5] = 0.9324695142031521;
    w[0] = 0.1713244923791704;
    w[1] = 0.3607615730481386;
    w[2] = 0.4679139345726910;
    w[3] = 0.4679139345726910;
    w[4] = 0.3607615730481386;
    w[5] = 0.1713244923791704;
    break;
  default:
    // Fallback to 2-point rule
    xi[0] = -0.5773502691896257;
    xi[1] = 0.5773502691896257;
    w[0] = 1.0;
    w[1] = 1.0;
  }
}

/**
 * @brief Generates 2D tensor-product Gauss quadrature rule on [-1,1]^2.
 *
 * @param ngaus Total number of Gauss points (must be a perfect square, e.g., 4,
 * 9, 16).
 * @param xi Output array for xi coordinates.
 * @param eta Output array for eta coordinates.
 * @param w Output array for weights.
 */
void gaussRectangle(int ngaus, double *xi, double *eta, double *w) {
  // Determine order p such that p*p = ngaus
  int p = int(std::round(std::sqrt(double(ngaus))));
  if (p * p != ngaus)
    p = 2; // Fallback

  std::vector<double> xi1, w1;
  gaussLegendre1D(p, xi1, w1);

  int idx = 0;
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < p; ++i, ++idx) {
      xi[idx] = xi1[i];
      eta[idx] = xi1[j];
      w[idx] = w1[i] * w1[j];
    }
  }
}

// =============================================================================
// ELEMENT MATRIX ASSEMBLY
// =============================================================================

/**
 * @brief Computes the local element stiffness matrix.
 *
 * Calculates Ae_{ab} = ∫ ∇N_a · ∇N_b dΩ.
 * Handles coordinate transformation for both Cartesian and Annular (Polar)
 * domains.
 *
 * @param dat Global data structure.
 * @param ie Element index.
 * @param kappa Diffusion coefficient.
 * @param Ae Output stiffness matrix (flattened).
 */
void element_stiffness(const data &dat, int ie, double kappa, PetscScalar *Ae) {
  int nloc = dat.nele_mat; // Total local DOFs
  int ngaus = dat.ngaus;

  // Zero initialize
  std::fill(Ae, Ae + nloc * nloc, 0.0);

  // Retrieve element node coordinates
  std::vector<double> xcoord(dat.node_elem), ycoord(dat.node_elem);
  for (int a = 0; a < dat.node_elem; ++a) {
    int nid = dat.elem[ie][a];
    xcoord[a] = dat.node[nid][0];
    ycoord[a] = dat.node[nid][1];
  }

  // Precompute scaling factors for annular domain mapping
  // Sr = dr/dx_param, Stheta = dtheta/dy_param
  double Sr = 1.0, Stheta = 1.0;
  if (dat.use_annular) {
    Sr = (dat.Rout - dat.Rin) / (dat.xmax - dat.xmin);
    Stheta = 2.0 * PI / (dat.ymax - dat.ymin);
  }

  // Loop over Gauss points
  for (int gp = 0; gp < ngaus; ++gp) {
    double xi = dat.xigaus[gp];
    double eta = dat.etagaus[gp];
    double w = dat.wgaus[gp];

    // 1. Geometry Interpolation (Jacobian)
    // We use Q1 shape functions for geometry mapping even if the field is
    // higher order.
    std::vector<double> dNdxi_geom(dat.node_elem), dNdeta_geom(dat.node_elem);
    deriv1Q1(dNdxi_geom.data(), dNdeta_geom.data(), xi, eta);

    double N_geom[4];
    shapeQ1(N_geom, xi, eta);

    double dxdxi = 0, dxdeta = 0, dydxi = 0, dydeta = 0;
    double x_gauss = 0.0; // Used as r_param in annular case

    for (int a = 0; a < dat.node_elem; ++a) {
      dxdxi += dNdxi_geom[a] * xcoord[a];
      dxdeta += dNdeta_geom[a] * xcoord[a];
      dydxi += dNdxi_geom[a] * ycoord[a];
      dydeta += dNdeta_geom[a] * ycoord[a];

      if (a < 4)
        x_gauss += N_geom[a] * xcoord[a];
    }

    double detJ = dxdxi * dydeta - dxdeta * dydxi;
    double invJ[2][2] = {{dydeta / detJ, -dxdeta / detJ},
                         {-dydxi / detJ, dxdxi / detJ}};

    // 2. Field Interpolation (Shape Function Derivatives)
    std::vector<double> dNdxi(nloc), dNdeta(nloc);
    deriv1Rect(dNdxi.data(), dNdeta.data(), xi, eta, dat);

    // Map derivatives to physical space (Cartesian parameter space)
    std::vector<double> dNdx(nloc), dNdy(nloc);
    for (int a = 0; a < nloc; ++a) {
      // dNdx = dNdxi * dxi/dx + dNdeta * deta/dx
      // invJ = [[dxi/dx, dxi/dy], [deta/dx, deta/dy]]
      dNdx[a] = invJ[0][0] * dNdxi[a] + invJ[1][0] * dNdeta[a];
      // dNdy = dNdxi * dxi/dy + dNdeta * deta/dy
      dNdy[a] = invJ[0][1] * dNdxi[a] + invJ[1][1] * dNdeta[a];
    }

    // 3. Assemble Stiffness Matrix
    if (dat.use_annular) {
      // --- POLAR COORDINATES ---
      // Transform Laplacian to polar: ∇²u = (1/r) ∂/∂r (r ∂u/∂r) + (1/r²)
      // ∂²u/∂θ² Weak form integrand involves: (∂u/∂r)(∂v/∂r) +
      // (1/r²)(∂u/∂θ)(∂v/∂θ) Integration measure: r dr dθ

      double r_phys = dat.Rin + Sr * (x_gauss - dat.xmin);

      // Weights derived from chain rule and integration measure:
      // Wr     = (Stheta / Sr) * r_phys
      // Wtheta = (Sr / (Stheta * r_phys))

      double Wr = (Stheta / Sr) * r_phys;
      double Wtheta = (Sr / (Stheta * r_phys));

      if (ie == 0 && gp == 0) {
        printf("DEBUG: ie=0 gp=0 r_phys=%.4f Sr=%.4f Stheta=%.4f Wr=%.4f "
               "Wtheta=%.4f detJ=%.4f\n",
               r_phys, Sr, Stheta, Wr, Wtheta, detJ);
      }
      if (std::abs(r_phys - 9.3) < 0.2 && gp == 0 && ie % 100 == 0) {
        printf("DEBUG: near 9.3 ie=%d r_phys=%.4f Wr=%.4f Wtheta=%.4f\n", ie,
               r_phys, Wr, Wtheta);
      }

      for (int a = 0; a < nloc; ++a) {
        for (int b = 0; b < nloc; ++b) {
          Ae[a * nloc + b] +=
              kappa * (Wr * dNdx[a] * dNdx[b] + Wtheta * dNdy[a] * dNdy[b]) *
              detJ * w;
        }
      }
    } else {
      // --- CARTESIAN COORDINATES ---
      // Standard Laplacian: ∇²u = ∂²u/∂x² + ∂²u/∂y²
      for (int a = 0; a < nloc; ++a) {
        for (int b = 0; b < nloc; ++b) {
          Ae[a * nloc + b] +=
              kappa * (dNdx[a] * dNdx[b] + dNdy[a] * dNdy[b]) * detJ * w;
        }
      }
    }
  }
}

/**
 * @brief Computes the local element mass matrix.
 *
 * Calculates Me_{ab} = ∫ N_a N_b dΩ.
 *
 * @param dat Global data structure.
 * @param ie Element index.
 * @param Me Output mass matrix (flattened).
 */
void element_mass(const data &dat, int ie, PetscScalar *Me) {
  int nloc = dat.nele_mat;
  int ngaus = dat.ngaus;
  std::fill(Me, Me + nloc * nloc, 0.0);

  std::vector<double> xcoord(dat.node_elem), ycoord(dat.node_elem);
  for (int a = 0; a < dat.node_elem; ++a) {
    int nid = dat.elem[ie][a];
    xcoord[a] = dat.node[nid][0];
    ycoord[a] = dat.node[nid][1];
  }

  double Sr = 1.0, Stheta = 1.0;
  if (dat.use_annular) {
    Sr = (dat.Rout - dat.Rin) / (dat.xmax - dat.xmin);
    Stheta = 2.0 * PI / (dat.ymax - dat.ymin);
  }

  for (int gp = 0; gp < ngaus; ++gp) {
    double xi = dat.xigaus[gp];
    double eta = dat.etagaus[gp];
    double w = dat.wgaus[gp];

    // Shape functions
    std::vector<double> Nloc(nloc);
    shapeRect(Nloc.data(), xi, eta, dat);

    // Geometry (Q1)
    std::vector<double> dNdxi_geom(dat.node_elem), dNdeta_geom(dat.node_elem);
    deriv1Q1(dNdxi_geom.data(), dNdeta_geom.data(), xi, eta);

    double N_geom[4];
    shapeQ1(N_geom, xi, eta);

    double dxdxi = 0, dxdeta = 0, dydxi = 0, dydeta = 0;
    double x_gauss = 0.0;

    for (int a = 0; a < dat.node_elem; ++a) {
      dxdxi += dNdxi_geom[a] * xcoord[a];
      dxdeta += dNdeta_geom[a] * xcoord[a];
      dydxi += dNdxi_geom[a] * ycoord[a];
      dydeta += dNdeta_geom[a] * ycoord[a];

      if (a < 4)
        x_gauss += N_geom[a] * xcoord[a];
    }
    double detJ = dxdxi * dydeta - dxdeta * dydxi;

    // Integration weight factor
    double weight_factor = 1.0;
    if (dat.use_annular) {
      double r_phys = dat.Rin + Sr * (x_gauss - dat.xmin);
      weight_factor = r_phys * Sr * Stheta; // r dr dtheta
    }

    for (int a = 0; a < nloc; ++a) {
      for (int b = 0; b < nloc; ++b) {
        Me[a * nloc + b] += Nloc[a] * Nloc[b] * weight_factor * detJ * w;
      }
    }
  }
}

// =============================================================================
// SHAPE FUNCTION DISPATCHERS
// =============================================================================

PetscErrorCode shapeRect(double *N, double xi, double eta, const data &dat) {
  int ne = dat.node_elem, nd = dat.dof_per_node;
  if (ne == 4 && nd == 1) {
    shapeQ1(N, xi, eta);
    return 0;
  } else if (ne == 9 && nd == 1) {
    shapeQ2(N, xi, eta);
    return 0;
  } else if (ne == 4 && nd == 4) {
    shapeHermiteCubic(N, xi, eta);
    return 0;
  } else if (ne == 9 && nd == 4) {
    shapeHermiteQuintic(N, xi, eta);
    return 0;
  } else
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
            "Unsupported element type in shapeRect");
  return 0;
}

PetscErrorCode deriv1Rect(double *dNdxi, double *dNdeta, double xi, double eta,
                          const data &dat) {
  int ne = dat.node_elem, nd = dat.dof_per_node;
  if (ne == 4 && nd == 1)
    deriv1Q1(dNdxi, dNdeta, xi, eta);
  else if (ne == 9 && nd == 1)
    deriv1Q2(dNdxi, dNdeta, xi, eta);
  else if (ne == 4 && nd == 4)
    deriv1HermiteCubic(dNdxi, dNdeta, xi, eta);
  else if (ne == 9 && nd == 4)
    deriv1HermiteQuintic(dNdxi, dNdeta, xi, eta);
  else
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
            "Unsupported element type in deriv1Rect");
  return 0; // Should be unreachable
}

// =============================================================================
// Q1 ELEMENTS (Linear)
// =============================================================================
void shapeQ1(double N[4], double xi, double eta) {
  N[0] = 0.25 * (1 - xi) * (1 - eta); // BL
  N[1] = 0.25 * (1 + xi) * (1 - eta); // BR
  N[2] = 0.25 * (1 + xi) * (1 + eta); // TR
  N[3] = 0.25 * (1 - xi) * (1 + eta); // TL
}

void deriv1Q1(double dNdxi[4], double dNdeta[4], double xi, double eta) {
  dNdxi[0] = -0.25 * (1 - eta);
  dNdeta[0] = -0.25 * (1 - xi);
  dNdxi[1] = 0.25 * (1 - eta);
  dNdeta[1] = -0.25 * (1 + xi);
  dNdxi[2] = 0.25 * (1 + eta);
  dNdeta[2] = 0.25 * (1 + xi);
  dNdxi[3] = -0.25 * (1 + eta);
  dNdeta[3] = 0.25 * (1 - xi);
}

// =============================================================================
// Q2 ELEMENTS (Quadratic)
// =============================================================================
void shapeQ2(double N[9], double xi, double eta) {
  // Temporary array M[] in natural grid order
  double M[9];
  M[0] = eta * xi * (eta - 1) * (xi - 1) / 4.0;  // (-1,-1)
  M[1] = -xi * (eta * eta - 1) * (xi - 1) / 2.0; // (-1, 0)
  M[2] = eta * xi * (eta + 1) * (xi - 1) / 4.0;  // (-1,+1)
  M[3] = -eta * (eta - 1) * (xi * xi - 1) / 2.0; // ( 0,-1)
  M[4] = (eta * eta - 1) * (xi * xi - 1);        // ( 0, 0)
  M[5] = -eta * (eta + 1) * (xi * xi - 1) / 2.0; // ( 0,+1)
  M[6] = eta * xi * (eta - 1) * (xi + 1) / 4.0;  // (+1,-1)
  M[7] = -xi * (eta * eta - 1) * (xi + 1) / 2.0; // (+1, 0)
  M[8] = eta * xi * (eta + 1) * (xi + 1) / 4.0;  // (+1,+1)

  // Permute to match element connectivity: Corners, Midsides, Center
  N[0] = M[0]; // BL
  N[1] = M[6]; // BR
  N[2] = M[8]; // TR
  N[3] = M[2]; // TL
  N[4] = M[3]; // Bottom
  N[5] = M[7]; // Right
  N[6] = M[5]; // Top
  N[7] = M[1]; // Left
  N[8] = M[4]; // Center
}

void deriv1Q2(double dNdxi[9], double dNdeta[9], double xi, double eta) {
  double dMx[9], dMy[9];

  // Derivatives of M (natural order)
  // ∂M/∂ξ
  dMx[0] = eta * xi * (eta - 1) / 4.0 + eta * (eta - 1) * (xi - 1) / 4.0;
  dMx[1] = -xi * (eta * eta - 1) / 2.0 + (0.5 - xi / 2.0) * (eta * eta - 1);
  dMx[2] = eta * xi * (eta + 1) / 4.0 + eta * (eta + 1) * (xi - 1) / 4.0;
  dMx[3] = -eta * (eta - 1) * xi;
  dMx[4] = 2.0 * xi * (eta * eta - 1);
  dMx[5] = -eta * (eta + 1) * xi;
  dMx[6] = eta * xi * (eta - 1) / 4.0 + eta * (eta - 1) * (xi + 1) / 4.0;
  dMx[7] = -xi * (eta * eta - 1) / 2.0 + (eta * eta - 1) * (-xi / 2.0 - 0.5);
  dMx[8] = eta * xi * (eta + 1) / 4.0 + eta * (eta + 1) * (xi + 1) / 4.0;

  // ∂M/∂η
  dMy[0] = eta * xi * (xi - 1) / 4.0 + xi * (eta - 1) * (xi - 1) / 4.0;
  dMy[1] = -eta * xi * (xi - 1);
  dMy[2] = eta * xi * (xi - 1) / 4.0 + xi * (eta + 1) * (xi - 1) / 4.0;
  dMy[3] = -eta * (xi * xi - 1) / 2.0 + (0.5 - eta / 2.0) * (xi * xi - 1);
  dMy[4] = 2.0 * eta * (xi * xi - 1);
  dMy[5] = -eta * (xi * xi - 1) / 2.0 + (-eta / 2.0 - 0.5) * (xi * xi - 1);
  dMy[6] = eta * xi * (xi + 1) / 4.0 + xi * (eta - 1) * (xi + 1) / 4.0;
  dMy[7] = -eta * xi * (xi + 1);
  dMy[8] = eta * xi * (xi + 1) / 4.0 + xi * (eta + 1) * (xi + 1) / 4.0;

  // Permute
  dNdxi[0] = dMx[0];
  dNdeta[0] = dMy[0];
  dNdxi[1] = dMx[6];
  dNdeta[1] = dMy[6];
  dNdxi[2] = dMx[8];
  dNdeta[2] = dMy[8];
  dNdxi[3] = dMx[2];
  dNdeta[3] = dMy[2];
  dNdxi[4] = dMx[3];
  dNdeta[4] = dMy[3];
  dNdxi[5] = dMx[7];
  dNdeta[5] = dMy[7];
  dNdxi[6] = dMx[5];
  dNdeta[6] = dMy[5];
  dNdxi[7] = dMx[1];
  dNdeta[7] = dMy[1];
  dNdxi[8] = dMx[4];
  dNdeta[8] = dMy[4];
}

// =============================================================================
// HERMITE ELEMENTS (Cubic)
// =============================================================================
// 1D Hermite basis functions on [0,1]
static inline double h00(double t) { return 2 * t * t * t - 3 * t * t + 1; }
static inline double h10(double t) { return t * t * t - 2 * t * t + t; }
static inline double h01(double t) { return -2 * t * t * t + 3 * t * t; }
static inline double h11(double t) { return t * t * t - t * t; }

static inline double dh00(double t) { return 6 * t * t - 6 * t; }
static inline double dh10(double t) { return 3 * t * t - 4 * t + 1; }
static inline double dh01(double t) { return -6 * t * t + 6 * t; }
static inline double dh11(double t) { return 3 * t * t - 2 * t; }

void shapeHermiteCubic(double N[16], double xi, double eta) {
  // Map [-1,1] -> [0,1]
  double tx = 0.5 * (xi + 1.0);
  double ty = 0.5 * (eta + 1.0);

  double X0 = h00(tx), X1 = h10(tx), X2 = h01(tx), X3 = h11(tx);
  double Y0 = h00(ty), Y1 = h10(ty), Y2 = h01(ty), Y3 = h11(ty);

  // Tensor product construction
  // BL
  N[0] = X0 * Y0;
  N[1] = X1 * Y0;
  N[2] = X0 * Y1;
  N[3] = X1 * Y1;
  // BR
  N[4] = X2 * Y0;
  N[5] = X3 * Y0;
  N[6] = X2 * Y1;
  N[7] = X3 * Y1;
  // TR
  N[8] = X2 * Y2;
  N[9] = X3 * Y2;
  N[10] = X2 * Y3;
  N[11] = X3 * Y3;
  // TL
  N[12] = X0 * Y2;
  N[13] = X1 * Y2;
  N[14] = X0 * Y3;
  N[15] = X1 * Y3;
}

void deriv1HermiteCubic(double dNdxi[16], double dNdeta[16], double xi,
                        double eta) {
  double tx = 0.5 * (xi + 1.0);
  double ty = 0.5 * (eta + 1.0);

  double X0 = h00(tx), X1 = h10(tx), X2 = h01(tx), X3 = h11(tx);
  double dX0 = dh00(tx) * 0.5, dX1 = dh10(tx) * 0.5, dX2 = dh01(tx) * 0.5,
         dX3 = dh11(tx) * 0.5;

  double Y0 = h00(ty), Y1 = h10(ty), Y2 = h01(ty), Y3 = h11(ty);
  double dY0 = dh00(ty) * 0.5, dY1 = dh10(ty) * 0.5, dY2 = dh01(ty) * 0.5,
         dY3 = dh11(ty) * 0.5;

  // BL
  dNdxi[0] = dX0 * Y0;
  dNdeta[0] = X0 * dY0;
  dNdxi[1] = dX1 * Y0;
  dNdeta[1] = X1 * dY0;
  dNdxi[2] = dX0 * Y1;
  dNdeta[2] = X0 * dY1;
  dNdxi[3] = dX1 * Y1;
  dNdeta[3] = X1 * dY1;
  // BR
  dNdxi[4] = dX2 * Y0;
  dNdeta[4] = X2 * dY0;
  dNdxi[5] = dX3 * Y0;
  dNdeta[5] = X3 * dY0;
  dNdxi[6] = dX2 * Y1;
  dNdeta[6] = X2 * dY1;
  dNdxi[7] = dX3 * Y1;
  dNdeta[7] = X3 * dY1;
  // TR
  dNdxi[8] = dX2 * Y2;
  dNdeta[8] = X2 * dY2;
  dNdxi[9] = dX3 * Y2;
  dNdeta[9] = X3 * dY2;
  dNdxi[10] = dX2 * Y3;
  dNdeta[10] = X2 * dY3;
  dNdxi[11] = dX3 * Y3;
  dNdeta[11] = X3 * dY3;
  // TL
  dNdxi[12] = dX0 * Y2;
  dNdeta[12] = X0 * dY2;
  dNdxi[13] = dX1 * Y2;
  dNdeta[13] = X1 * dY2;
  dNdxi[14] = dX0 * Y3;
  dNdeta[14] = X0 * dY3;
  dNdxi[15] = dX1 * Y3;
  dNdeta[15] = X1 * dY3;
}

// =============================================================================
// HERMITE ELEMENTS (Quintic)
// =============================================================================

typedef double (*ShapeFuncPtr)(double xi, double eta);

static ShapeFuncPtr N_funcs[36] = {N0,  N1,  N2,  N3,  N4,  N5,  N6,  N7,  N8,
                                   N9,  N10, N11, N12, N13, N14, N15, N16, N17,
                                   N18, N19, N20, N21, N22, N23, N24, N25, N26,
                                   N27, N28, N29, N30, N31, N32, N33, N34, N35};

static ShapeFuncPtr dNdxi_funcs[36] = {
    dN0_dxi,  dN1_dxi,  dN2_dxi,  dN3_dxi,  dN4_dxi,  dN5_dxi,
    dN6_dxi,  dN7_dxi,  dN8_dxi,  dN9_dxi,  dN10_dxi, dN11_dxi,
    dN12_dxi, dN13_dxi, dN14_dxi, dN15_dxi, dN16_dxi, dN17_dxi,
    dN18_dxi, dN19_dxi, dN20_dxi, dN21_dxi, dN22_dxi, dN23_dxi,
    dN24_dxi, dN25_dxi, dN26_dxi, dN27_dxi, dN28_dxi, dN29_dxi,
    dN30_dxi, dN31_dxi, dN32_dxi, dN33_dxi, dN34_dxi, dN35_dxi};

static ShapeFuncPtr dNdeta_funcs[36] = {
    dN0_deta,  dN1_deta,  dN2_deta,  dN3_deta,  dN4_deta,  dN5_deta,
    dN6_deta,  dN7_deta,  dN8_deta,  dN9_deta,  dN10_deta, dN11_deta,
    dN12_deta, dN13_deta, dN14_deta, dN15_deta, dN16_deta, dN17_deta,
    dN18_deta, dN19_deta, dN20_deta, dN21_deta, dN22_deta, dN23_deta,
    dN24_deta, dN25_deta, dN26_deta, dN27_deta, dN28_deta, dN29_deta,
    dN30_deta, dN31_deta, dN32_deta, dN33_deta, dN34_deta, dN35_deta};

void shapeHermiteQuintic(double N[36], double xi, double eta) {
  for (int k = 0; k < 36; ++k) {
    N[k] = N_funcs[k](xi, eta);
  }
}

void deriv1HermiteQuintic(double dNdxi[36], double dNdeta[36], double xi,
                          double eta) {
  for (int k = 0; k < 36; ++k) {
    dNdxi[k] = dNdxi_funcs[k](xi, eta);
    dNdeta[k] = dNdeta_funcs[k](xi, eta);
  }
}

// =============================================================================
// ELEMENT SEARCH
// =============================================================================

/**
 * @brief Locates the element containing a given point (x,y).
 *
 * Performs a brute-force search over all elements.
 *
 * @param x X-coordinate of the point.
 * @param y Y-coordinate of the point.
 * @param dat Global data structure.
 * @param iel Output element index.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode locelem(double x, double y, const data &dat, int &iel) {
  for (int e = 0; e < dat.nelem; ++e) {
    int n1 = dat.elem[e][0]; // BL
    int n2 = dat.elem[e][1]; // BR
    int n4 = dat.elem[e][3]; // TL

    double xmin = dat.node[n1][0];
    double xmax = dat.node[n2][0];
    double ymin = dat.node[n1][1];
    double ymax = dat.node[n4][1];

    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
      iel = e;
      return 0;
    }
  }

  SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
          "locelem: point (%g, %g) is not inside any element", x, y);
}
