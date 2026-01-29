/**
 * @file constants.h
 * @brief Physical and Mathematical Constants.
 *
 * This file defines various physical constants (electron mass, Planck's
 * constant, etc.) and mathematical constants (PI, square roots) used throughout
 * the simulation.
 *
 * @author Sathwik Bharadwaj (Original)
 * @author Arya Retheeshan (Documentation)
 * @date 2025-11-24
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

// =============================================================================
// MATHEMATICAL CONSTANTS
// =============================================================================

/* Linux has PI defined in /usr/include/math.h, but we define it here for
 * portability */
#ifndef PI
#define PI 3.14159265358979323846
#endif

// Note: The complex constant ei=(0.0,1.0) is typically defined in the dmat
// library.

// =============================================================================
// PHYSICAL CONSTANTS
// =============================================================================

// m_o*c^2 (Rest mass energy of electron in eV)
#define xm0 0.5109990615E+6

// hbar * c (Planck's constant times speed of light in eV*cm)
#define hbarc 1.9732705359E-5

// Length scale: 1 Angstrom in cm
#define xscale 1.0E-8

// Vacuum permittivity in Farad/meter
#define vac_epsilon 8.854187817e-12

// LRR's scaling constants
#define cscale 3.809984038754390
#define del 7.619968077508781
#define xl0 1.0e-8

// Boltzmann's constant (eV/K)
// Value is 1/11604.5
#define Kboltz (1.0e0 / 11604.5e0)

// =============================================================================
// MAGNETIC FIELD PARAMETERS
// =============================================================================
// r^2/B_0: Square of Landau orbit radius per unit magnetic field
// r2b1 is divided by b0 (in Tesla) to get Landau R_o^2 in cm^2
#define r2b1 6.58212224061E-12

// hbar * omega_B per unit magnetic field
// hwb is multiplied by b0 (in Tesla) to get hbar*omega_b in eV
#define hwb 1.1576764E-4

// =============================================================================
// NUMERICAL CONSTANTS (k.P Hamiltonian)
// =============================================================================
// Pre-computed square roots and their inverses for efficiency

#define sqi2 0.70710678118654752440 // 1/sqrt(2)
#define sqi3 0.57735026918962576451 // 1/sqrt(3)
#define sqi6 0.40824829046386301637 // 1/sqrt(6)
#define sq2 1.41421356237309504880  // sqrt(2)
#define sq3 1.73205080756887729353  // sqrt(3)

#define sq6 2.449489742783178098197  // sqrt(6)
#define sq32 1.224744871391589049099 // sqrt(3/2)
#define sq23 0.816496580927726032732 // sqrt(2/3)

#endif // CONSTANTS_H
