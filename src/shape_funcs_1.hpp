/**
 * @file shape_funcs.hpp
 * @brief Auto-Generated Quintic Hermite Shape Functions.
 *
 * This file contains 36 quintic Hermite shape functions (N0-N35) and their
 * first derivatives (dN0_dxi, dN0_deta, etc.) for 9-node Hermite elements.
 * These functions were auto-generated from symbolic mathematics software.
 *
 * DO NOT EDIT MANUALLY - regenerate from symbolic code if changes are needed.
 *
 * @date Auto-generated
 */

inline double N0(double xi, double eta) {
  return (9.0 / 16.0) * pow(eta, 5) * pow(xi, 5) -
         3.0 / 8.0 * pow(eta, 5) * pow(xi, 4) -
         15.0 / 16.0 * pow(eta, 5) * pow(xi, 3) +
         (3.0 / 4.0) * pow(eta, 5) * pow(xi, 2) -
         3.0 / 8.0 * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 4) +
         (5.0 / 8.0) * pow(eta, 4) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 2) -
         15.0 / 16.0 * pow(eta, 3) * pow(xi, 5) +
         (5.0 / 8.0) * pow(eta, 3) * pow(xi, 4) +
         (25.0 / 16.0) * pow(eta, 3) * pow(xi, 3) -
         5.0 / 4.0 * pow(eta, 3) * pow(xi, 2) +
         (3.0 / 4.0) * pow(eta, 2) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 2) * pow(xi, 3) + pow(eta, 2) * pow(xi, 2);
}

inline double N1(double xi, double eta) {
  return (3.0 / 16.0) * pow(eta, 5) * pow(xi, 5) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 4) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 5) * pow(xi, 2) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 3) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 2) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 5) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 3) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 2) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 5) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 2);
}

inline double N2(double xi, double eta) {
  return (3.0 / 16.0) * pow(eta, 5) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 5) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 5) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 5) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 4) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 2) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 2) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 2) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 2);
}

inline double N3(double xi, double eta) {
  return (1.0 / 16.0) * pow(eta, 5) * pow(xi, 5) -
         1.0 / 16.0 * pow(eta, 5) * pow(xi, 4) -
         1.0 / 16.0 * pow(eta, 5) * pow(xi, 3) +
         (1.0 / 16.0) * pow(eta, 5) * pow(xi, 2) -
         1.0 / 16.0 * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 16.0) * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 16.0) * pow(eta, 4) * pow(xi, 3) -
         1.0 / 16.0 * pow(eta, 4) * pow(xi, 2) -
         1.0 / 16.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 16.0) * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 16.0) * pow(eta, 3) * pow(xi, 3) -
         1.0 / 16.0 * pow(eta, 3) * pow(xi, 2) +
         (1.0 / 16.0) * pow(eta, 2) * pow(xi, 5) -
         1.0 / 16.0 * pow(eta, 2) * pow(xi, 4) -
         1.0 / 16.0 * pow(eta, 2) * pow(xi, 3) +
         (1.0 / 16.0) * pow(eta, 2) * pow(xi, 2);
}

inline double N4(double xi, double eta) {
  return -9.0 / 16.0 * pow(eta, 5) * pow(xi, 5) -
         3.0 / 8.0 * pow(eta, 5) * pow(xi, 4) +
         (15.0 / 16.0) * pow(eta, 5) * pow(xi, 3) +
         (3.0 / 4.0) * pow(eta, 5) * pow(xi, 2) +
         (3.0 / 8.0) * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         5.0 / 8.0 * pow(eta, 4) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 2) +
         (15.0 / 16.0) * pow(eta, 3) * pow(xi, 5) +
         (5.0 / 8.0) * pow(eta, 3) * pow(xi, 4) -
         25.0 / 16.0 * pow(eta, 3) * pow(xi, 3) -
         5.0 / 4.0 * pow(eta, 3) * pow(xi, 2) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 4) +
         (5.0 / 4.0) * pow(eta, 2) * pow(xi, 3) + pow(eta, 2) * pow(xi, 2);
}

inline double N5(double xi, double eta) {
  return (3.0 / 16.0) * pow(eta, 5) * pow(xi, 5) +
         (3.0 / 16.0) * pow(eta, 5) * pow(xi, 4) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 2) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 3) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 2) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 5) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 3) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 2) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 2);
}

inline double N6(double xi, double eta) {
  return -3.0 / 16.0 * pow(eta, 5) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 5) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 5) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 5) * pow(xi, 2) +
         (3.0 / 16.0) * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 2) +
         (3.0 / 16.0) * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 3) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 2) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 2) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 2);
}

inline double N7(double xi, double eta) {
  return (1.0 / 16.0) * pow(eta, 5) * pow(xi, 5) +
         (1.0 / 16.0) * pow(eta, 5) * pow(xi, 4) -
         1.0 / 16.0 * pow(eta, 5) * pow(xi, 3) -
         1.0 / 16.0 * pow(eta, 5) * pow(xi, 2) -
         1.0 / 16.0 * pow(eta, 4) * pow(xi, 5) -
         1.0 / 16.0 * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 16.0) * pow(eta, 4) * pow(xi, 3) +
         (1.0 / 16.0) * pow(eta, 4) * pow(xi, 2) -
         1.0 / 16.0 * pow(eta, 3) * pow(xi, 5) -
         1.0 / 16.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 16.0) * pow(eta, 3) * pow(xi, 3) +
         (1.0 / 16.0) * pow(eta, 3) * pow(xi, 2) +
         (1.0 / 16.0) * pow(eta, 2) * pow(xi, 5) +
         (1.0 / 16.0) * pow(eta, 2) * pow(xi, 4) -
         1.0 / 16.0 * pow(eta, 2) * pow(xi, 3) -
         1.0 / 16.0 * pow(eta, 2) * pow(xi, 2);
}

inline double N8(double xi, double eta) {
  return (9.0 / 16.0) * pow(eta, 5) * pow(xi, 5) +
         (3.0 / 8.0) * pow(eta, 5) * pow(xi, 4) -
         15.0 / 16.0 * pow(eta, 5) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 5) * pow(xi, 2) +
         (3.0 / 8.0) * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         5.0 / 8.0 * pow(eta, 4) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 2) -
         15.0 / 16.0 * pow(eta, 3) * pow(xi, 5) -
         5.0 / 8.0 * pow(eta, 3) * pow(xi, 4) +
         (25.0 / 16.0) * pow(eta, 3) * pow(xi, 3) +
         (5.0 / 4.0) * pow(eta, 3) * pow(xi, 2) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 4) +
         (5.0 / 4.0) * pow(eta, 2) * pow(xi, 3) + pow(eta, 2) * pow(xi, 2);
}

inline double N9(double xi, double eta) {
  return -3.0 / 16.0 * pow(eta, 5) * pow(xi, 5) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 4) +
         (3.0 / 16.0) * pow(eta, 5) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 5) * pow(xi, 2) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 3) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 2) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 5) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 3) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 2) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 2);
}

inline double N10(double xi, double eta) {
  return -3.0 / 16.0 * pow(eta, 5) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 5) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 5) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 5) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 4) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 4) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 2) +
         (3.0 / 16.0) * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 3) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 2) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 2) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 2) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 2);
}

inline double N11(double xi, double eta) {
  return (1.0 / 16.0) * pow(eta, 5) * pow(xi, 5) +
         (1.0 / 16.0) * pow(eta, 5) * pow(xi, 4) -
         1.0 / 16.0 * pow(eta, 5) * pow(xi, 3) -
         1.0 / 16.0 * pow(eta, 5) * pow(xi, 2) +
         (1.0 / 16.0) * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 16.0) * pow(eta, 4) * pow(xi, 4) -
         1.0 / 16.0 * pow(eta, 4) * pow(xi, 3) -
         1.0 / 16.0 * pow(eta, 4) * pow(xi, 2) -
         1.0 / 16.0 * pow(eta, 3) * pow(xi, 5) -
         1.0 / 16.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 16.0) * pow(eta, 3) * pow(xi, 3) +
         (1.0 / 16.0) * pow(eta, 3) * pow(xi, 2) -
         1.0 / 16.0 * pow(eta, 2) * pow(xi, 5) -
         1.0 / 16.0 * pow(eta, 2) * pow(xi, 4) +
         (1.0 / 16.0) * pow(eta, 2) * pow(xi, 3) +
         (1.0 / 16.0) * pow(eta, 2) * pow(xi, 2);
}

inline double N12(double xi, double eta) {
  return -9.0 / 16.0 * pow(eta, 5) * pow(xi, 5) +
         (3.0 / 8.0) * pow(eta, 5) * pow(xi, 4) +
         (15.0 / 16.0) * pow(eta, 5) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 5) * pow(xi, 2) -
         3.0 / 8.0 * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 4) +
         (5.0 / 8.0) * pow(eta, 4) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 2) +
         (15.0 / 16.0) * pow(eta, 3) * pow(xi, 5) -
         5.0 / 8.0 * pow(eta, 3) * pow(xi, 4) -
         25.0 / 16.0 * pow(eta, 3) * pow(xi, 3) +
         (5.0 / 4.0) * pow(eta, 3) * pow(xi, 2) +
         (3.0 / 4.0) * pow(eta, 2) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 2) * pow(xi, 3) + pow(eta, 2) * pow(xi, 2);
}

inline double N13(double xi, double eta) {
  return -3.0 / 16.0 * pow(eta, 5) * pow(xi, 5) +
         (3.0 / 16.0) * pow(eta, 5) * pow(xi, 4) +
         (3.0 / 16.0) * pow(eta, 5) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 2) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 8.0) * pow(eta, 4) * pow(xi, 3) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 2) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 5) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 3) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 2) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 5) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 2);
}

inline double N14(double xi, double eta) {
  return (3.0 / 16.0) * pow(eta, 5) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 5) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 5) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 5) * pow(xi, 2) +
         (3.0 / 16.0) * pow(eta, 4) * pow(xi, 5) -
         1.0 / 8.0 * pow(eta, 4) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 3) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 5) +
         (1.0 / 8.0) * pow(eta, 2) * pow(xi, 4) +
         (5.0 / 16.0) * pow(eta, 2) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 2);
}

inline double N15(double xi, double eta) {
  return (1.0 / 16.0) * pow(eta, 5) * pow(xi, 5) -
         1.0 / 16.0 * pow(eta, 5) * pow(xi, 4) -
         1.0 / 16.0 * pow(eta, 5) * pow(xi, 3) +
         (1.0 / 16.0) * pow(eta, 5) * pow(xi, 2) +
         (1.0 / 16.0) * pow(eta, 4) * pow(xi, 5) -
         1.0 / 16.0 * pow(eta, 4) * pow(xi, 4) -
         1.0 / 16.0 * pow(eta, 4) * pow(xi, 3) +
         (1.0 / 16.0) * pow(eta, 4) * pow(xi, 2) -
         1.0 / 16.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 16.0) * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 16.0) * pow(eta, 3) * pow(xi, 3) -
         1.0 / 16.0 * pow(eta, 3) * pow(xi, 2) -
         1.0 / 16.0 * pow(eta, 2) * pow(xi, 5) +
         (1.0 / 16.0) * pow(eta, 2) * pow(xi, 4) +
         (1.0 / 16.0) * pow(eta, 2) * pow(xi, 3) -
         1.0 / 16.0 * pow(eta, 2) * pow(xi, 2);
}

inline double N16(double xi, double eta) {
  return (3.0 / 4.0) * pow(eta, 5) * pow(xi, 4) -
         3.0 / 2.0 * pow(eta, 5) * pow(xi, 2) + (3.0 / 4.0) * pow(eta, 5) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 4) + pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 4) - 5.0 / 4.0 * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 2.0) * pow(eta, 3) * pow(xi, 2) - 5.0 / 4.0 * pow(eta, 3) +
         pow(eta, 2) * pow(xi, 4) - 2 * pow(eta, 2) * pow(xi, 2) + pow(eta, 2);
}

inline double N17(double xi, double eta) {
  return (3.0 / 4.0) * pow(eta, 5) * pow(xi, 5) -
         3.0 / 2.0 * pow(eta, 5) * pow(xi, 3) + (3.0 / 4.0) * pow(eta, 5) * xi -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 5) + pow(eta, 4) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 4) * xi - 5.0 / 4.0 * pow(eta, 3) * pow(xi, 5) +
         (5.0 / 2.0) * pow(eta, 3) * pow(xi, 3) - 5.0 / 4.0 * pow(eta, 3) * xi +
         pow(eta, 2) * pow(xi, 5) - 2 * pow(eta, 2) * pow(xi, 3) +
         pow(eta, 2) * xi;
}

inline double N18(double xi, double eta) {
  return (1.0 / 4.0) * pow(eta, 5) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 2) + (1.0 / 4.0) * pow(eta, 5) -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 4) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 4) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 3) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 2) + (1.0 / 4.0) * pow(eta, 2);
}

inline double N19(double xi, double eta) {
  return (1.0 / 4.0) * pow(eta, 5) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 3) + (1.0 / 4.0) * pow(eta, 5) * xi -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 4) * pow(xi, 3) - 1.0 / 4.0 * pow(eta, 4) * xi -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) - 1.0 / 4.0 * pow(eta, 3) * xi +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 3) + (1.0 / 4.0) * pow(eta, 2) * xi;
}

inline double N20(double xi, double eta) {
  return -3.0 / 4.0 * pow(eta, 4) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 4) +
         (5.0 / 4.0) * pow(eta, 4) * pow(xi, 3) + pow(eta, 4) * pow(xi, 2) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 5) + pow(eta, 2) * pow(xi, 4) -
         5.0 / 2.0 * pow(eta, 2) * pow(xi, 3) - 2 * pow(eta, 2) * pow(xi, 2) -
         3.0 / 4.0 * pow(xi, 5) - 1.0 / 2.0 * pow(xi, 4) +
         (5.0 / 4.0) * pow(xi, 3) + pow(xi, 2);
}

inline double N21(double xi, double eta) {
  return (1.0 / 4.0) * pow(eta, 4) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 2) * pow(xi, 3) +
         (1.0 / 2.0) * pow(eta, 2) * pow(xi, 2) + (1.0 / 4.0) * pow(xi, 5) +
         (1.0 / 4.0) * pow(xi, 4) - 1.0 / 4.0 * pow(xi, 3) -
         1.0 / 4.0 * pow(xi, 2);
}

inline double N22(double xi, double eta) {
  return -3.0 / 4.0 * pow(eta, 5) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 4) +
         (5.0 / 4.0) * pow(eta, 5) * pow(xi, 3) + pow(eta, 5) * pow(xi, 2) +
         (3.0 / 2.0) * pow(eta, 3) * pow(xi, 5) + pow(eta, 3) * pow(xi, 4) -
         5.0 / 2.0 * pow(eta, 3) * pow(xi, 3) - 2 * pow(eta, 3) * pow(xi, 2) -
         3.0 / 4.0 * eta * pow(xi, 5) - 1.0 / 2.0 * eta * pow(xi, 4) +
         (5.0 / 4.0) * eta * pow(xi, 3) + eta * pow(xi, 2);
}

inline double N23(double xi, double eta) {
  return (1.0 / 4.0) * pow(eta, 5) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 5) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 5) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 5) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 2) +
         (1.0 / 4.0) * eta * pow(xi, 5) + (1.0 / 4.0) * eta * pow(xi, 4) -
         1.0 / 4.0 * eta * pow(xi, 3) - 1.0 / 4.0 * eta * pow(xi, 2);
}

inline double N24(double xi, double eta) {
  return -3.0 / 4.0 * pow(eta, 5) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 5) * pow(xi, 2) - 3.0 / 4.0 * pow(eta, 5) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 4) + pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 4) + (5.0 / 4.0) * pow(eta, 3) * pow(xi, 4) -
         5.0 / 2.0 * pow(eta, 3) * pow(xi, 2) + (5.0 / 4.0) * pow(eta, 3) +
         pow(eta, 2) * pow(xi, 4) - 2 * pow(eta, 2) * pow(xi, 2) + pow(eta, 2);
}

inline double N25(double xi, double eta) {
  return -3.0 / 4.0 * pow(eta, 5) * pow(xi, 5) +
         (3.0 / 2.0) * pow(eta, 5) * pow(xi, 3) - 3.0 / 4.0 * pow(eta, 5) * xi -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 5) + pow(eta, 4) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 4) * xi + (5.0 / 4.0) * pow(eta, 3) * pow(xi, 5) -
         5.0 / 2.0 * pow(eta, 3) * pow(xi, 3) + (5.0 / 4.0) * pow(eta, 3) * xi +
         pow(eta, 2) * pow(xi, 5) - 2 * pow(eta, 2) * pow(xi, 3) +
         pow(eta, 2) * xi;
}

inline double N26(double xi, double eta) {
  return (1.0 / 4.0) * pow(eta, 5) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 2) + (1.0 / 4.0) * pow(eta, 5) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 2) + (1.0 / 4.0) * pow(eta, 4) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 3) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 2) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 2);
}

inline double N27(double xi, double eta) {
  return (1.0 / 4.0) * pow(eta, 5) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 3) + (1.0 / 4.0) * pow(eta, 5) * xi +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 3) + (1.0 / 4.0) * pow(eta, 4) * xi -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) - 1.0 / 4.0 * pow(eta, 3) * xi -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 2) * pow(xi, 3) - 1.0 / 4.0 * pow(eta, 2) * xi;
}

inline double N28(double xi, double eta) {
  return (3.0 / 4.0) * pow(eta, 4) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 4) * pow(xi, 3) + pow(eta, 4) * pow(xi, 2) -
         3.0 / 2.0 * pow(eta, 2) * pow(xi, 5) + pow(eta, 2) * pow(xi, 4) +
         (5.0 / 2.0) * pow(eta, 2) * pow(xi, 3) - 2 * pow(eta, 2) * pow(xi, 2) +
         (3.0 / 4.0) * pow(xi, 5) - 1.0 / 2.0 * pow(xi, 4) -
         5.0 / 4.0 * pow(xi, 3) + pow(xi, 2);
}

inline double N29(double xi, double eta) {
  return (1.0 / 4.0) * pow(eta, 4) * pow(xi, 5) -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 2) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 2) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 2) + (1.0 / 4.0) * pow(xi, 5) -
         1.0 / 4.0 * pow(xi, 4) - 1.0 / 4.0 * pow(xi, 3) +
         (1.0 / 4.0) * pow(xi, 2);
}

inline double N30(double xi, double eta) {
  return (3.0 / 4.0) * pow(eta, 5) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 5) * pow(xi, 3) + pow(eta, 5) * pow(xi, 2) -
         3.0 / 2.0 * pow(eta, 3) * pow(xi, 5) + pow(eta, 3) * pow(xi, 4) +
         (5.0 / 2.0) * pow(eta, 3) * pow(xi, 3) - 2 * pow(eta, 3) * pow(xi, 2) +
         (3.0 / 4.0) * eta * pow(xi, 5) - 1.0 / 2.0 * eta * pow(xi, 4) -
         5.0 / 4.0 * eta * pow(xi, 3) + eta * pow(xi, 2);
}

inline double N31(double xi, double eta) {
  return (1.0 / 4.0) * pow(eta, 5) * pow(xi, 5) -
         1.0 / 4.0 * pow(eta, 5) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 5) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 5) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 2) + (1.0 / 4.0) * eta * pow(xi, 5) -
         1.0 / 4.0 * eta * pow(xi, 4) - 1.0 / 4.0 * eta * pow(xi, 3) +
         (1.0 / 4.0) * eta * pow(xi, 2);
}

inline double N32(double xi, double eta) {
  return pow(eta, 4) * pow(xi, 4) - 2 * pow(eta, 4) * pow(xi, 2) + pow(eta, 4) -
         2 * pow(eta, 2) * pow(xi, 4) + 4 * pow(eta, 2) * pow(xi, 2) -
         2 * pow(eta, 2) + pow(xi, 4) - 2 * pow(xi, 2) + 1;
}

inline double N33(double xi, double eta) {
  return pow(eta, 4) * pow(xi, 5) - 2 * pow(eta, 4) * pow(xi, 3) +
         pow(eta, 4) * xi - 2 * pow(eta, 2) * pow(xi, 5) +
         4 * pow(eta, 2) * pow(xi, 3) - 2 * pow(eta, 2) * xi + pow(xi, 5) -
         2 * pow(xi, 3) + xi;
}

inline double N34(double xi, double eta) {
  return pow(eta, 5) * pow(xi, 4) - 2 * pow(eta, 5) * pow(xi, 2) + pow(eta, 5) -
         2 * pow(eta, 3) * pow(xi, 4) + 4 * pow(eta, 3) * pow(xi, 2) -
         2 * pow(eta, 3) + eta * pow(xi, 4) - 2 * eta * pow(xi, 2) + eta;
}

inline double N35(double xi, double eta) {
  return pow(eta, 5) * pow(xi, 5) - 2 * pow(eta, 5) * pow(xi, 3) +
         pow(eta, 5) * xi - 2 * pow(eta, 3) * pow(xi, 5) +
         4 * pow(eta, 3) * pow(xi, 3) - 2 * pow(eta, 3) * xi +
         eta * pow(xi, 5) - 2 * eta * pow(xi, 3) + eta * xi;
}

inline double dN0_dxi(double xi, double eta) {
  return (45.0 / 16.0) * pow(eta, 5) * pow(xi, 4) -
         3.0 / 2.0 * pow(eta, 5) * pow(xi, 3) -
         45.0 / 16.0 * pow(eta, 5) * pow(xi, 2) +
         (3.0 / 2.0) * pow(eta, 5) * xi -
         15.0 / 8.0 * pow(eta, 4) * pow(xi, 4) + pow(eta, 4) * pow(xi, 3) +
         (15.0 / 8.0) * pow(eta, 4) * pow(xi, 2) - pow(eta, 4) * xi -
         75.0 / 16.0 * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 2.0) * pow(eta, 3) * pow(xi, 3) +
         (75.0 / 16.0) * pow(eta, 3) * pow(xi, 2) -
         5.0 / 2.0 * pow(eta, 3) * xi +
         (15.0 / 4.0) * pow(eta, 2) * pow(xi, 4) -
         2 * pow(eta, 2) * pow(xi, 3) - 15.0 / 4.0 * pow(eta, 2) * pow(xi, 2) +
         2 * pow(eta, 2) * xi;
}

inline double dN1_dxi(double xi, double eta) {
  return (15.0 / 16.0) * pow(eta, 5) * pow(xi, 4) -
         3.0 / 4.0 * pow(eta, 5) * pow(xi, 3) -
         9.0 / 16.0 * pow(eta, 5) * pow(xi, 2) +
         (3.0 / 8.0) * pow(eta, 5) * xi - 5.0 / 8.0 * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 4) * pow(xi, 3) +
         (3.0 / 8.0) * pow(eta, 4) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 4) * xi -
         25.0 / 16.0 * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 4.0) * pow(eta, 3) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 3) * pow(xi, 2) -
         5.0 / 8.0 * pow(eta, 3) * xi + (5.0 / 4.0) * pow(eta, 2) * pow(xi, 4) -
         pow(eta, 2) * pow(xi, 3) - 3.0 / 4.0 * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 2.0) * pow(eta, 2) * xi;
}

inline double dN2_dxi(double xi, double eta) {
  return (15.0 / 16.0) * pow(eta, 5) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 5) * pow(xi, 2) +
         (1.0 / 2.0) * pow(eta, 5) * xi -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 4) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 4) * xi - 15.0 / 16.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 3) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 3) * xi +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 2.0) * pow(eta, 2) * xi;
}

inline double dN3_dxi(double xi, double eta) {
  return (5.0 / 16.0) * pow(eta, 5) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 5) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 2) +
         (1.0 / 8.0) * pow(eta, 5) * xi -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 4) * pow(xi, 2) -
         1.0 / 8.0 * pow(eta, 4) * xi - 5.0 / 16.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 3) * pow(xi, 2) -
         1.0 / 8.0 * pow(eta, 3) * xi +
         (5.0 / 16.0) * pow(eta, 2) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 2) + (1.0 / 8.0) * pow(eta, 2) * xi;
}

inline double dN4_dxi(double xi, double eta) {
  return -45.0 / 16.0 * pow(eta, 5) * pow(xi, 4) -
         3.0 / 2.0 * pow(eta, 5) * pow(xi, 3) +
         (45.0 / 16.0) * pow(eta, 5) * pow(xi, 2) +
         (3.0 / 2.0) * pow(eta, 5) * xi +
         (15.0 / 8.0) * pow(eta, 4) * pow(xi, 4) + pow(eta, 4) * pow(xi, 3) -
         15.0 / 8.0 * pow(eta, 4) * pow(xi, 2) - pow(eta, 4) * xi +
         (75.0 / 16.0) * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 2.0) * pow(eta, 3) * pow(xi, 3) -
         75.0 / 16.0 * pow(eta, 3) * pow(xi, 2) - 5.0 / 2.0 * pow(eta, 3) * xi -
         15.0 / 4.0 * pow(eta, 2) * pow(xi, 4) - 2 * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 2) * pow(xi, 2) + 2 * pow(eta, 2) * xi;
}

inline double dN5_dxi(double xi, double eta) {
  return (15.0 / 16.0) * pow(eta, 5) * pow(xi, 4) +
         (3.0 / 4.0) * pow(eta, 5) * pow(xi, 3) -
         9.0 / 16.0 * pow(eta, 5) * pow(xi, 2) - 3.0 / 8.0 * pow(eta, 5) * xi -
         5.0 / 8.0 * pow(eta, 4) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 3) +
         (3.0 / 8.0) * pow(eta, 4) * pow(xi, 2) +
         (1.0 / 4.0) * pow(eta, 4) * xi -
         25.0 / 16.0 * pow(eta, 3) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 3) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 3) * pow(xi, 2) +
         (5.0 / 8.0) * pow(eta, 3) * xi +
         (5.0 / 4.0) * pow(eta, 2) * pow(xi, 4) + pow(eta, 2) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 2) - 1.0 / 2.0 * pow(eta, 2) * xi;
}

inline double dN6_dxi(double xi, double eta) {
  return -15.0 / 16.0 * pow(eta, 5) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 5) * pow(xi, 2) +
         (1.0 / 2.0) * pow(eta, 5) * xi +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 4) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 2) - 1.0 / 2.0 * pow(eta, 4) * xi +
         (15.0 / 16.0) * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 3) * pow(xi, 2) - 1.0 / 2.0 * pow(eta, 3) * xi -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 2.0) * pow(eta, 2) * xi;
}

inline double dN7_dxi(double xi, double eta) {
  return (5.0 / 16.0) * pow(eta, 5) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 5) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 2) - 1.0 / 8.0 * pow(eta, 5) * xi -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 4) * pow(xi, 2) +
         (1.0 / 8.0) * pow(eta, 4) * xi -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 3) * pow(xi, 2) +
         (1.0 / 8.0) * pow(eta, 3) * xi +
         (5.0 / 16.0) * pow(eta, 2) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 2) - 1.0 / 8.0 * pow(eta, 2) * xi;
}

inline double dN8_dxi(double xi, double eta) {
  return (45.0 / 16.0) * pow(eta, 5) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 5) * pow(xi, 3) -
         45.0 / 16.0 * pow(eta, 5) * pow(xi, 2) - 3.0 / 2.0 * pow(eta, 5) * xi +
         (15.0 / 8.0) * pow(eta, 4) * pow(xi, 4) + pow(eta, 4) * pow(xi, 3) -
         15.0 / 8.0 * pow(eta, 4) * pow(xi, 2) - pow(eta, 4) * xi -
         75.0 / 16.0 * pow(eta, 3) * pow(xi, 4) -
         5.0 / 2.0 * pow(eta, 3) * pow(xi, 3) +
         (75.0 / 16.0) * pow(eta, 3) * pow(xi, 2) +
         (5.0 / 2.0) * pow(eta, 3) * xi -
         15.0 / 4.0 * pow(eta, 2) * pow(xi, 4) - 2 * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 2) * pow(xi, 2) + 2 * pow(eta, 2) * xi;
}

inline double dN9_dxi(double xi, double eta) {
  return -15.0 / 16.0 * pow(eta, 5) * pow(xi, 4) -
         3.0 / 4.0 * pow(eta, 5) * pow(xi, 3) +
         (9.0 / 16.0) * pow(eta, 5) * pow(xi, 2) +
         (3.0 / 8.0) * pow(eta, 5) * xi - 5.0 / 8.0 * pow(eta, 4) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 3) +
         (3.0 / 8.0) * pow(eta, 4) * pow(xi, 2) +
         (1.0 / 4.0) * pow(eta, 4) * xi +
         (25.0 / 16.0) * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 4.0) * pow(eta, 3) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 3) * pow(xi, 2) - 5.0 / 8.0 * pow(eta, 3) * xi +
         (5.0 / 4.0) * pow(eta, 2) * pow(xi, 4) + pow(eta, 2) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 2) - 1.0 / 2.0 * pow(eta, 2) * xi;
}

inline double dN10_dxi(double xi, double eta) {
  return -15.0 / 16.0 * pow(eta, 5) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 5) * pow(xi, 2) +
         (1.0 / 2.0) * pow(eta, 5) * xi -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 2) +
         (1.0 / 2.0) * pow(eta, 4) * xi +
         (15.0 / 16.0) * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 3) * pow(xi, 2) - 1.0 / 2.0 * pow(eta, 3) * xi +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 2) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 2) - 1.0 / 2.0 * pow(eta, 2) * xi;
}

inline double dN11_dxi(double xi, double eta) {
  return (5.0 / 16.0) * pow(eta, 5) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 5) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 2) - 1.0 / 8.0 * pow(eta, 5) * xi +
         (5.0 / 16.0) * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 4) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 4) * pow(xi, 2) - 1.0 / 8.0 * pow(eta, 4) * xi -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 3) * pow(xi, 2) +
         (1.0 / 8.0) * pow(eta, 3) * xi -
         5.0 / 16.0 * pow(eta, 2) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 2) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 8.0) * pow(eta, 2) * xi;
}

inline double dN12_dxi(double xi, double eta) {
  return -45.0 / 16.0 * pow(eta, 5) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 5) * pow(xi, 3) +
         (45.0 / 16.0) * pow(eta, 5) * pow(xi, 2) -
         3.0 / 2.0 * pow(eta, 5) * xi - 15.0 / 8.0 * pow(eta, 4) * pow(xi, 4) +
         pow(eta, 4) * pow(xi, 3) + (15.0 / 8.0) * pow(eta, 4) * pow(xi, 2) -
         pow(eta, 4) * xi + (75.0 / 16.0) * pow(eta, 3) * pow(xi, 4) -
         5.0 / 2.0 * pow(eta, 3) * pow(xi, 3) -
         75.0 / 16.0 * pow(eta, 3) * pow(xi, 2) +
         (5.0 / 2.0) * pow(eta, 3) * xi +
         (15.0 / 4.0) * pow(eta, 2) * pow(xi, 4) -
         2 * pow(eta, 2) * pow(xi, 3) - 15.0 / 4.0 * pow(eta, 2) * pow(xi, 2) +
         2 * pow(eta, 2) * xi;
}

inline double dN13_dxi(double xi, double eta) {
  return -15.0 / 16.0 * pow(eta, 5) * pow(xi, 4) +
         (3.0 / 4.0) * pow(eta, 5) * pow(xi, 3) +
         (9.0 / 16.0) * pow(eta, 5) * pow(xi, 2) -
         3.0 / 8.0 * pow(eta, 5) * xi - 5.0 / 8.0 * pow(eta, 4) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 4) * pow(xi, 3) +
         (3.0 / 8.0) * pow(eta, 4) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 4) * xi +
         (25.0 / 16.0) * pow(eta, 3) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 3) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 3) * pow(xi, 2) +
         (5.0 / 8.0) * pow(eta, 3) * xi +
         (5.0 / 4.0) * pow(eta, 2) * pow(xi, 4) - pow(eta, 2) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 2) + (1.0 / 2.0) * pow(eta, 2) * xi;
}

inline double dN14_dxi(double xi, double eta) {
  return (15.0 / 16.0) * pow(eta, 5) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 5) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 5) * pow(xi, 2) +
         (1.0 / 2.0) * pow(eta, 5) * xi +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 4) -
         1.0 / 2.0 * pow(eta, 4) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 2) +
         (1.0 / 2.0) * pow(eta, 4) * xi -
         15.0 / 16.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 3) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 3) * xi - 15.0 / 16.0 * pow(eta, 2) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 2) * xi;
}

inline double dN15_dxi(double xi, double eta) {
  return (5.0 / 16.0) * pow(eta, 5) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 5) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 5) * pow(xi, 2) +
         (1.0 / 8.0) * pow(eta, 5) * xi +
         (5.0 / 16.0) * pow(eta, 4) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 4) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 4) * pow(xi, 2) +
         (1.0 / 8.0) * pow(eta, 4) * xi -
         5.0 / 16.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 3) * pow(xi, 2) -
         1.0 / 8.0 * pow(eta, 3) * xi - 5.0 / 16.0 * pow(eta, 2) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 2) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 2) - 1.0 / 8.0 * pow(eta, 2) * xi;
}

inline double dN16_dxi(double xi, double eta) {
  return 3 * pow(eta, 5) * pow(xi, 3) - 3 * pow(eta, 5) * xi -
         2 * pow(eta, 4) * pow(xi, 3) + 2 * pow(eta, 4) * xi -
         5 * pow(eta, 3) * pow(xi, 3) + 5 * pow(eta, 3) * xi +
         4 * pow(eta, 2) * pow(xi, 3) - 4 * pow(eta, 2) * xi;
}

inline double dN17_dxi(double xi, double eta) {
  return (15.0 / 4.0) * pow(eta, 5) * pow(xi, 4) -
         9.0 / 2.0 * pow(eta, 5) * pow(xi, 2) + (3.0 / 4.0) * pow(eta, 5) -
         5.0 / 2.0 * pow(eta, 4) * pow(xi, 4) + 3 * pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 4) - 25.0 / 4.0 * pow(eta, 3) * pow(xi, 4) +
         (15.0 / 2.0) * pow(eta, 3) * pow(xi, 2) - 5.0 / 4.0 * pow(eta, 3) +
         5 * pow(eta, 2) * pow(xi, 4) - 6 * pow(eta, 2) * pow(xi, 2) +
         pow(eta, 2);
}

inline double dN18_dxi(double xi, double eta) {
  return pow(eta, 5) * pow(xi, 3) - pow(eta, 5) * xi -
         pow(eta, 4) * pow(xi, 3) + pow(eta, 4) * xi -
         pow(eta, 3) * pow(xi, 3) + pow(eta, 3) * xi +
         pow(eta, 2) * pow(xi, 3) - pow(eta, 2) * xi;
}

inline double dN19_dxi(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 5) * pow(xi, 4) -
         3.0 / 2.0 * pow(eta, 5) * pow(xi, 2) + (1.0 / 4.0) * pow(eta, 5) -
         5.0 / 4.0 * pow(eta, 4) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 4) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 4) -
         5.0 / 4.0 * pow(eta, 3) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 3) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 3) +
         (5.0 / 4.0) * pow(eta, 2) * pow(xi, 4) -
         3.0 / 2.0 * pow(eta, 2) * pow(xi, 2) + (1.0 / 4.0) * pow(eta, 2);
}

inline double dN20_dxi(double xi, double eta) {
  return -15.0 / 4.0 * pow(eta, 4) * pow(xi, 4) - 2 * pow(eta, 4) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 4) * pow(xi, 2) + 2 * pow(eta, 4) * xi +
         (15.0 / 2.0) * pow(eta, 2) * pow(xi, 4) +
         4 * pow(eta, 2) * pow(xi, 3) - 15.0 / 2.0 * pow(eta, 2) * pow(xi, 2) -
         4 * pow(eta, 2) * xi - 15.0 / 4.0 * pow(xi, 4) - 2 * pow(xi, 3) +
         (15.0 / 4.0) * pow(xi, 2) + 2 * xi;
}

inline double dN21_dxi(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 4) * pow(xi, 4) + pow(eta, 4) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 4) * pow(xi, 2) - 1.0 / 2.0 * pow(eta, 4) * xi -
         5.0 / 2.0 * pow(eta, 2) * pow(xi, 4) - 2 * pow(eta, 2) * pow(xi, 3) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 2) + pow(eta, 2) * xi +
         (5.0 / 4.0) * pow(xi, 4) + pow(xi, 3) - 3.0 / 4.0 * pow(xi, 2) -
         1.0 / 2.0 * xi;
}

inline double dN22_dxi(double xi, double eta) {
  return -15.0 / 4.0 * pow(eta, 5) * pow(xi, 4) - 2 * pow(eta, 5) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 5) * pow(xi, 2) + 2 * pow(eta, 5) * xi +
         (15.0 / 2.0) * pow(eta, 3) * pow(xi, 4) +
         4 * pow(eta, 3) * pow(xi, 3) - 15.0 / 2.0 * pow(eta, 3) * pow(xi, 2) -
         4 * pow(eta, 3) * xi - 15.0 / 4.0 * eta * pow(xi, 4) -
         2 * eta * pow(xi, 3) + (15.0 / 4.0) * eta * pow(xi, 2) + 2 * eta * xi;
}

inline double dN23_dxi(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 5) * pow(xi, 4) + pow(eta, 5) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 5) * pow(xi, 2) - 1.0 / 2.0 * pow(eta, 5) * xi -
         5.0 / 2.0 * pow(eta, 3) * pow(xi, 4) - 2 * pow(eta, 3) * pow(xi, 3) +
         (3.0 / 2.0) * pow(eta, 3) * pow(xi, 2) + pow(eta, 3) * xi +
         (5.0 / 4.0) * eta * pow(xi, 4) + eta * pow(xi, 3) -
         3.0 / 4.0 * eta * pow(xi, 2) - 1.0 / 2.0 * eta * xi;
}

inline double dN24_dxi(double xi, double eta) {
  return -3 * pow(eta, 5) * pow(xi, 3) + 3 * pow(eta, 5) * xi -
         2 * pow(eta, 4) * pow(xi, 3) + 2 * pow(eta, 4) * xi +
         5 * pow(eta, 3) * pow(xi, 3) - 5 * pow(eta, 3) * xi +
         4 * pow(eta, 2) * pow(xi, 3) - 4 * pow(eta, 2) * xi;
}

inline double dN25_dxi(double xi, double eta) {
  return -15.0 / 4.0 * pow(eta, 5) * pow(xi, 4) +
         (9.0 / 2.0) * pow(eta, 5) * pow(xi, 2) - 3.0 / 4.0 * pow(eta, 5) -
         5.0 / 2.0 * pow(eta, 4) * pow(xi, 4) + 3 * pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 4) + (25.0 / 4.0) * pow(eta, 3) * pow(xi, 4) -
         15.0 / 2.0 * pow(eta, 3) * pow(xi, 2) + (5.0 / 4.0) * pow(eta, 3) +
         5 * pow(eta, 2) * pow(xi, 4) - 6 * pow(eta, 2) * pow(xi, 2) +
         pow(eta, 2);
}

inline double dN26_dxi(double xi, double eta) {
  return pow(eta, 5) * pow(xi, 3) - pow(eta, 5) * xi +
         pow(eta, 4) * pow(xi, 3) - pow(eta, 4) * xi -
         pow(eta, 3) * pow(xi, 3) + pow(eta, 3) * xi -
         pow(eta, 2) * pow(xi, 3) + pow(eta, 2) * xi;
}

inline double dN27_dxi(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 5) * pow(xi, 4) -
         3.0 / 2.0 * pow(eta, 5) * pow(xi, 2) + (1.0 / 4.0) * pow(eta, 5) +
         (5.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         3.0 / 2.0 * pow(eta, 4) * pow(xi, 2) + (1.0 / 4.0) * pow(eta, 4) -
         5.0 / 4.0 * pow(eta, 3) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 3) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 3) -
         5.0 / 4.0 * pow(eta, 2) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 2) - 1.0 / 4.0 * pow(eta, 2);
}

inline double dN28_dxi(double xi, double eta) {
  return (15.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         2 * pow(eta, 4) * pow(xi, 3) - 15.0 / 4.0 * pow(eta, 4) * pow(xi, 2) +
         2 * pow(eta, 4) * xi - 15.0 / 2.0 * pow(eta, 2) * pow(xi, 4) +
         4 * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 2.0) * pow(eta, 2) * pow(xi, 2) - 4 * pow(eta, 2) * xi +
         (15.0 / 4.0) * pow(xi, 4) - 2 * pow(xi, 3) - 15.0 / 4.0 * pow(xi, 2) +
         2 * xi;
}

inline double dN29_dxi(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 4) * pow(xi, 4) - pow(eta, 4) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 4) * pow(xi, 2) + (1.0 / 2.0) * pow(eta, 4) * xi -
         5.0 / 2.0 * pow(eta, 2) * pow(xi, 4) + 2 * pow(eta, 2) * pow(xi, 3) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 2) - pow(eta, 2) * xi +
         (5.0 / 4.0) * pow(xi, 4) - pow(xi, 3) - 3.0 / 4.0 * pow(xi, 2) +
         (1.0 / 2.0) * xi;
}

inline double dN30_dxi(double xi, double eta) {
  return (15.0 / 4.0) * pow(eta, 5) * pow(xi, 4) -
         2 * pow(eta, 5) * pow(xi, 3) - 15.0 / 4.0 * pow(eta, 5) * pow(xi, 2) +
         2 * pow(eta, 5) * xi - 15.0 / 2.0 * pow(eta, 3) * pow(xi, 4) +
         4 * pow(eta, 3) * pow(xi, 3) +
         (15.0 / 2.0) * pow(eta, 3) * pow(xi, 2) - 4 * pow(eta, 3) * xi +
         (15.0 / 4.0) * eta * pow(xi, 4) - 2 * eta * pow(xi, 3) -
         15.0 / 4.0 * eta * pow(xi, 2) + 2 * eta * xi;
}

inline double dN31_dxi(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 5) * pow(xi, 4) - pow(eta, 5) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 5) * pow(xi, 2) + (1.0 / 2.0) * pow(eta, 5) * xi -
         5.0 / 2.0 * pow(eta, 3) * pow(xi, 4) + 2 * pow(eta, 3) * pow(xi, 3) +
         (3.0 / 2.0) * pow(eta, 3) * pow(xi, 2) - pow(eta, 3) * xi +
         (5.0 / 4.0) * eta * pow(xi, 4) - eta * pow(xi, 3) -
         3.0 / 4.0 * eta * pow(xi, 2) + (1.0 / 2.0) * eta * xi;
}

inline double dN32_dxi(double xi, double eta) {
  return 4 * pow(eta, 4) * pow(xi, 3) - 4 * pow(eta, 4) * xi -
         8 * pow(eta, 2) * pow(xi, 3) + 8 * pow(eta, 2) * xi + 4 * pow(xi, 3) -
         4 * xi;
}

inline double dN33_dxi(double xi, double eta) {
  return 5 * pow(eta, 4) * pow(xi, 4) - 6 * pow(eta, 4) * pow(xi, 2) +
         pow(eta, 4) - 10 * pow(eta, 2) * pow(xi, 4) +
         12 * pow(eta, 2) * pow(xi, 2) - 2 * pow(eta, 2) + 5 * pow(xi, 4) -
         6 * pow(xi, 2) + 1;
}

inline double dN34_dxi(double xi, double eta) {
  return 4 * pow(eta, 5) * pow(xi, 3) - 4 * pow(eta, 5) * xi -
         8 * pow(eta, 3) * pow(xi, 3) + 8 * pow(eta, 3) * xi +
         4 * eta * pow(xi, 3) - 4 * eta * xi;
}

inline double dN35_dxi(double xi, double eta) {
  return 5 * pow(eta, 5) * pow(xi, 4) - 6 * pow(eta, 5) * pow(xi, 2) +
         pow(eta, 5) - 10 * pow(eta, 3) * pow(xi, 4) +
         12 * pow(eta, 3) * pow(xi, 2) - 2 * pow(eta, 3) +
         5 * eta * pow(xi, 4) - 6 * eta * pow(xi, 2) + eta;
}

inline double dN0_deta(double xi, double eta) {
  return (45.0 / 16.0) * pow(eta, 4) * pow(xi, 5) -
         15.0 / 8.0 * pow(eta, 4) * pow(xi, 4) -
         75.0 / 16.0 * pow(eta, 4) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 4) * pow(xi, 2) -
         3.0 / 2.0 * pow(eta, 3) * pow(xi, 5) + pow(eta, 3) * pow(xi, 4) +
         (5.0 / 2.0) * pow(eta, 3) * pow(xi, 3) - 2 * pow(eta, 3) * pow(xi, 2) -
         45.0 / 16.0 * pow(eta, 2) * pow(xi, 5) +
         (15.0 / 8.0) * pow(eta, 2) * pow(xi, 4) +
         (75.0 / 16.0) * pow(eta, 2) * pow(xi, 3) -
         15.0 / 4.0 * pow(eta, 2) * pow(xi, 2) +
         (3.0 / 2.0) * eta * pow(xi, 5) - eta * pow(xi, 4) -
         5.0 / 2.0 * eta * pow(xi, 3) + 2 * eta * pow(xi, 2);
}

inline double dN1_deta(double xi, double eta) {
  return (15.0 / 16.0) * pow(eta, 4) * pow(xi, 5) -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 4) -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 2) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 5) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 4) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 2.0) * eta * pow(xi, 5) - 1.0 / 2.0 * eta * pow(xi, 4) -
         1.0 / 2.0 * eta * pow(xi, 3) + (1.0 / 2.0) * eta * pow(xi, 2);
}

inline double dN2_deta(double xi, double eta) {
  return (15.0 / 16.0) * pow(eta, 4) * pow(xi, 5) -
         5.0 / 8.0 * pow(eta, 4) * pow(xi, 4) -
         25.0 / 16.0 * pow(eta, 4) * pow(xi, 3) +
         (5.0 / 4.0) * pow(eta, 4) * pow(xi, 2) -
         3.0 / 4.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 4.0) * pow(eta, 3) * pow(xi, 3) - pow(eta, 3) * pow(xi, 2) -
         9.0 / 16.0 * pow(eta, 2) * pow(xi, 5) +
         (3.0 / 8.0) * pow(eta, 2) * pow(xi, 4) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 2) + (3.0 / 8.0) * eta * pow(xi, 5) -
         1.0 / 4.0 * eta * pow(xi, 4) - 5.0 / 8.0 * eta * pow(xi, 3) +
         (1.0 / 2.0) * eta * pow(xi, 2);
}

inline double dN3_deta(double xi, double eta) {
  return (5.0 / 16.0) * pow(eta, 4) * pow(xi, 5) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 3) +
         (5.0 / 16.0) * pow(eta, 4) * pow(xi, 2) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 5) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 4) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 8.0) * eta * pow(xi, 5) - 1.0 / 8.0 * eta * pow(xi, 4) -
         1.0 / 8.0 * eta * pow(xi, 3) + (1.0 / 8.0) * eta * pow(xi, 2);
}

inline double dN4_deta(double xi, double eta) {
  return -45.0 / 16.0 * pow(eta, 4) * pow(xi, 5) -
         15.0 / 8.0 * pow(eta, 4) * pow(xi, 4) +
         (75.0 / 16.0) * pow(eta, 4) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 4) * pow(xi, 2) +
         (3.0 / 2.0) * pow(eta, 3) * pow(xi, 5) + pow(eta, 3) * pow(xi, 4) -
         5.0 / 2.0 * pow(eta, 3) * pow(xi, 3) - 2 * pow(eta, 3) * pow(xi, 2) +
         (45.0 / 16.0) * pow(eta, 2) * pow(xi, 5) +
         (15.0 / 8.0) * pow(eta, 2) * pow(xi, 4) -
         75.0 / 16.0 * pow(eta, 2) * pow(xi, 3) -
         15.0 / 4.0 * pow(eta, 2) * pow(xi, 2) - 3.0 / 2.0 * eta * pow(xi, 5) -
         eta * pow(xi, 4) + (5.0 / 2.0) * eta * pow(xi, 3) +
         2 * eta * pow(xi, 2);
}

inline double dN5_deta(double xi, double eta) {
  return (15.0 / 16.0) * pow(eta, 4) * pow(xi, 5) +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 4) -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 2) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 5) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 4) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 2.0) * eta * pow(xi, 5) + (1.0 / 2.0) * eta * pow(xi, 4) -
         1.0 / 2.0 * eta * pow(xi, 3) - 1.0 / 2.0 * eta * pow(xi, 2);
}

inline double dN6_deta(double xi, double eta) {
  return -15.0 / 16.0 * pow(eta, 4) * pow(xi, 5) -
         5.0 / 8.0 * pow(eta, 4) * pow(xi, 4) +
         (25.0 / 16.0) * pow(eta, 4) * pow(xi, 3) +
         (5.0 / 4.0) * pow(eta, 4) * pow(xi, 2) +
         (3.0 / 4.0) * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 3) * pow(xi, 3) - pow(eta, 3) * pow(xi, 2) +
         (9.0 / 16.0) * pow(eta, 2) * pow(xi, 5) +
         (3.0 / 8.0) * pow(eta, 2) * pow(xi, 4) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 2) - 3.0 / 8.0 * eta * pow(xi, 5) -
         1.0 / 4.0 * eta * pow(xi, 4) + (5.0 / 8.0) * eta * pow(xi, 3) +
         (1.0 / 2.0) * eta * pow(xi, 2);
}

inline double dN7_deta(double xi, double eta) {
  return (5.0 / 16.0) * pow(eta, 4) * pow(xi, 5) +
         (5.0 / 16.0) * pow(eta, 4) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 3) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 2) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 5) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 5) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 4) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 8.0) * eta * pow(xi, 5) + (1.0 / 8.0) * eta * pow(xi, 4) -
         1.0 / 8.0 * eta * pow(xi, 3) - 1.0 / 8.0 * eta * pow(xi, 2);
}

inline double dN8_deta(double xi, double eta) {
  return (45.0 / 16.0) * pow(eta, 4) * pow(xi, 5) +
         (15.0 / 8.0) * pow(eta, 4) * pow(xi, 4) -
         75.0 / 16.0 * pow(eta, 4) * pow(xi, 3) -
         15.0 / 4.0 * pow(eta, 4) * pow(xi, 2) +
         (3.0 / 2.0) * pow(eta, 3) * pow(xi, 5) + pow(eta, 3) * pow(xi, 4) -
         5.0 / 2.0 * pow(eta, 3) * pow(xi, 3) - 2 * pow(eta, 3) * pow(xi, 2) -
         45.0 / 16.0 * pow(eta, 2) * pow(xi, 5) -
         15.0 / 8.0 * pow(eta, 2) * pow(xi, 4) +
         (75.0 / 16.0) * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 2) * pow(xi, 2) -
         3.0 / 2.0 * eta * pow(xi, 5) - eta * pow(xi, 4) +
         (5.0 / 2.0) * eta * pow(xi, 3) + 2 * eta * pow(xi, 2);
}

inline double dN9_deta(double xi, double eta) {
  return -15.0 / 16.0 * pow(eta, 4) * pow(xi, 5) -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 4) +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 2) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 5) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 4) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 2.0) * eta * pow(xi, 5) + (1.0 / 2.0) * eta * pow(xi, 4) -
         1.0 / 2.0 * eta * pow(xi, 3) - 1.0 / 2.0 * eta * pow(xi, 2);
}

inline double dN10_deta(double xi, double eta) {
  return -15.0 / 16.0 * pow(eta, 4) * pow(xi, 5) -
         5.0 / 8.0 * pow(eta, 4) * pow(xi, 4) +
         (25.0 / 16.0) * pow(eta, 4) * pow(xi, 3) +
         (5.0 / 4.0) * pow(eta, 4) * pow(xi, 2) -
         3.0 / 4.0 * pow(eta, 3) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 4) +
         (5.0 / 4.0) * pow(eta, 3) * pow(xi, 3) + pow(eta, 3) * pow(xi, 2) +
         (9.0 / 16.0) * pow(eta, 2) * pow(xi, 5) +
         (3.0 / 8.0) * pow(eta, 2) * pow(xi, 4) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 2) + (3.0 / 8.0) * eta * pow(xi, 5) +
         (1.0 / 4.0) * eta * pow(xi, 4) - 5.0 / 8.0 * eta * pow(xi, 3) -
         1.0 / 2.0 * eta * pow(xi, 2);
}

inline double dN11_deta(double xi, double eta) {
  return (5.0 / 16.0) * pow(eta, 4) * pow(xi, 5) +
         (5.0 / 16.0) * pow(eta, 4) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 3) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 2) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 3) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 5) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 4) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 3) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 2) -
         1.0 / 8.0 * eta * pow(xi, 5) - 1.0 / 8.0 * eta * pow(xi, 4) +
         (1.0 / 8.0) * eta * pow(xi, 3) + (1.0 / 8.0) * eta * pow(xi, 2);
}

inline double dN12_deta(double xi, double eta) {
  return -45.0 / 16.0 * pow(eta, 4) * pow(xi, 5) +
         (15.0 / 8.0) * pow(eta, 4) * pow(xi, 4) +
         (75.0 / 16.0) * pow(eta, 4) * pow(xi, 3) -
         15.0 / 4.0 * pow(eta, 4) * pow(xi, 2) -
         3.0 / 2.0 * pow(eta, 3) * pow(xi, 5) + pow(eta, 3) * pow(xi, 4) +
         (5.0 / 2.0) * pow(eta, 3) * pow(xi, 3) - 2 * pow(eta, 3) * pow(xi, 2) +
         (45.0 / 16.0) * pow(eta, 2) * pow(xi, 5) -
         15.0 / 8.0 * pow(eta, 2) * pow(xi, 4) -
         75.0 / 16.0 * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 2) * pow(xi, 2) +
         (3.0 / 2.0) * eta * pow(xi, 5) - eta * pow(xi, 4) -
         5.0 / 2.0 * eta * pow(xi, 3) + 2 * eta * pow(xi, 2);
}

inline double dN13_deta(double xi, double eta) {
  return -15.0 / 16.0 * pow(eta, 4) * pow(xi, 5) +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 4) +
         (15.0 / 16.0) * pow(eta, 4) * pow(xi, 3) -
         15.0 / 16.0 * pow(eta, 4) * pow(xi, 2) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 5) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 4) +
         (1.0 / 2.0) * pow(eta, 3) * pow(xi, 3) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 2) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 5) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 4) -
         15.0 / 16.0 * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 2) +
         (1.0 / 2.0) * eta * pow(xi, 5) - 1.0 / 2.0 * eta * pow(xi, 4) -
         1.0 / 2.0 * eta * pow(xi, 3) + (1.0 / 2.0) * eta * pow(xi, 2);
}

inline double dN14_deta(double xi, double eta) {
  return (15.0 / 16.0) * pow(eta, 4) * pow(xi, 5) -
         5.0 / 8.0 * pow(eta, 4) * pow(xi, 4) -
         25.0 / 16.0 * pow(eta, 4) * pow(xi, 3) +
         (5.0 / 4.0) * pow(eta, 4) * pow(xi, 2) +
         (3.0 / 4.0) * pow(eta, 3) * pow(xi, 5) -
         1.0 / 2.0 * pow(eta, 3) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 3) * pow(xi, 3) + pow(eta, 3) * pow(xi, 2) -
         9.0 / 16.0 * pow(eta, 2) * pow(xi, 5) +
         (3.0 / 8.0) * pow(eta, 2) * pow(xi, 4) +
         (15.0 / 16.0) * pow(eta, 2) * pow(xi, 3) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 2) - 3.0 / 8.0 * eta * pow(xi, 5) +
         (1.0 / 4.0) * eta * pow(xi, 4) + (5.0 / 8.0) * eta * pow(xi, 3) -
         1.0 / 2.0 * eta * pow(xi, 2);
}

inline double dN15_deta(double xi, double eta) {
  return (5.0 / 16.0) * pow(eta, 4) * pow(xi, 5) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 4) -
         5.0 / 16.0 * pow(eta, 4) * pow(xi, 3) +
         (5.0 / 16.0) * pow(eta, 4) * pow(xi, 2) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 5) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 4) -
         1.0 / 4.0 * pow(eta, 3) * pow(xi, 3) +
         (1.0 / 4.0) * pow(eta, 3) * pow(xi, 2) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 5) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 4) +
         (3.0 / 16.0) * pow(eta, 2) * pow(xi, 3) -
         3.0 / 16.0 * pow(eta, 2) * pow(xi, 2) - 1.0 / 8.0 * eta * pow(xi, 5) +
         (1.0 / 8.0) * eta * pow(xi, 4) + (1.0 / 8.0) * eta * pow(xi, 3) -
         1.0 / 8.0 * eta * pow(xi, 2);
}

inline double dN16_deta(double xi, double eta) {
  return (15.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         15.0 / 2.0 * pow(eta, 4) * pow(xi, 2) + (15.0 / 4.0) * pow(eta, 4) -
         2 * pow(eta, 3) * pow(xi, 4) + 4 * pow(eta, 3) * pow(xi, 2) -
         2 * pow(eta, 3) - 15.0 / 4.0 * pow(eta, 2) * pow(xi, 4) +
         (15.0 / 2.0) * pow(eta, 2) * pow(xi, 2) - 15.0 / 4.0 * pow(eta, 2) +
         2 * eta * pow(xi, 4) - 4 * eta * pow(xi, 2) + 2 * eta;
}

inline double dN17_deta(double xi, double eta) {
  return (15.0 / 4.0) * pow(eta, 4) * pow(xi, 5) -
         15.0 / 2.0 * pow(eta, 4) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 4) * xi - 2 * pow(eta, 3) * pow(xi, 5) +
         4 * pow(eta, 3) * pow(xi, 3) - 2 * pow(eta, 3) * xi -
         15.0 / 4.0 * pow(eta, 2) * pow(xi, 5) +
         (15.0 / 2.0) * pow(eta, 2) * pow(xi, 3) -
         15.0 / 4.0 * pow(eta, 2) * xi + 2 * eta * pow(xi, 5) -
         4 * eta * pow(xi, 3) + 2 * eta * xi;
}

inline double dN18_deta(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         5.0 / 2.0 * pow(eta, 4) * pow(xi, 2) + (5.0 / 4.0) * pow(eta, 4) -
         pow(eta, 3) * pow(xi, 4) + 2 * pow(eta, 3) * pow(xi, 2) - pow(eta, 3) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 2) - 3.0 / 4.0 * pow(eta, 2) +
         (1.0 / 2.0) * eta * pow(xi, 4) - eta * pow(xi, 2) + (1.0 / 2.0) * eta;
}

inline double dN19_deta(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 4) * pow(xi, 5) -
         5.0 / 2.0 * pow(eta, 4) * pow(xi, 3) + (5.0 / 4.0) * pow(eta, 4) * xi -
         pow(eta, 3) * pow(xi, 5) + 2 * pow(eta, 3) * pow(xi, 3) -
         pow(eta, 3) * xi - 3.0 / 4.0 * pow(eta, 2) * pow(xi, 5) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 3) - 3.0 / 4.0 * pow(eta, 2) * xi +
         (1.0 / 2.0) * eta * pow(xi, 5) - eta * pow(xi, 3) +
         (1.0 / 2.0) * eta * xi;
}

inline double dN20_deta(double xi, double eta) {
  return -3 * pow(eta, 3) * pow(xi, 5) - 2 * pow(eta, 3) * pow(xi, 4) +
         5 * pow(eta, 3) * pow(xi, 3) + 4 * pow(eta, 3) * pow(xi, 2) +
         3 * eta * pow(xi, 5) + 2 * eta * pow(xi, 4) - 5 * eta * pow(xi, 3) -
         4 * eta * pow(xi, 2);
}

inline double dN21_deta(double xi, double eta) {
  return pow(eta, 3) * pow(xi, 5) + pow(eta, 3) * pow(xi, 4) -
         pow(eta, 3) * pow(xi, 3) - pow(eta, 3) * pow(xi, 2) -
         eta * pow(xi, 5) - eta * pow(xi, 4) + eta * pow(xi, 3) +
         eta * pow(xi, 2);
}

inline double dN22_deta(double xi, double eta) {
  return -15.0 / 4.0 * pow(eta, 4) * pow(xi, 5) -
         5.0 / 2.0 * pow(eta, 4) * pow(xi, 4) +
         (25.0 / 4.0) * pow(eta, 4) * pow(xi, 3) +
         5 * pow(eta, 4) * pow(xi, 2) + (9.0 / 2.0) * pow(eta, 2) * pow(xi, 5) +
         3 * pow(eta, 2) * pow(xi, 4) - 15.0 / 2.0 * pow(eta, 2) * pow(xi, 3) -
         6 * pow(eta, 2) * pow(xi, 2) - 3.0 / 4.0 * pow(xi, 5) -
         1.0 / 2.0 * pow(xi, 4) + (5.0 / 4.0) * pow(xi, 3) + pow(xi, 2);
}

inline double dN23_deta(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 4) * pow(xi, 5) +
         (5.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 4) * pow(xi, 3) -
         5.0 / 4.0 * pow(eta, 4) * pow(xi, 2) -
         3.0 / 2.0 * pow(eta, 2) * pow(xi, 5) -
         3.0 / 2.0 * pow(eta, 2) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 3) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 2) + (1.0 / 4.0) * pow(xi, 5) +
         (1.0 / 4.0) * pow(xi, 4) - 1.0 / 4.0 * pow(xi, 3) -
         1.0 / 4.0 * pow(xi, 2);
}

inline double dN24_deta(double xi, double eta) {
  return -15.0 / 4.0 * pow(eta, 4) * pow(xi, 4) +
         (15.0 / 2.0) * pow(eta, 4) * pow(xi, 2) - 15.0 / 4.0 * pow(eta, 4) -
         2 * pow(eta, 3) * pow(xi, 4) + 4 * pow(eta, 3) * pow(xi, 2) -
         2 * pow(eta, 3) + (15.0 / 4.0) * pow(eta, 2) * pow(xi, 4) -
         15.0 / 2.0 * pow(eta, 2) * pow(xi, 2) + (15.0 / 4.0) * pow(eta, 2) +
         2 * eta * pow(xi, 4) - 4 * eta * pow(xi, 2) + 2 * eta;
}

inline double dN25_deta(double xi, double eta) {
  return -15.0 / 4.0 * pow(eta, 4) * pow(xi, 5) +
         (15.0 / 2.0) * pow(eta, 4) * pow(xi, 3) -
         15.0 / 4.0 * pow(eta, 4) * xi - 2 * pow(eta, 3) * pow(xi, 5) +
         4 * pow(eta, 3) * pow(xi, 3) - 2 * pow(eta, 3) * xi +
         (15.0 / 4.0) * pow(eta, 2) * pow(xi, 5) -
         15.0 / 2.0 * pow(eta, 2) * pow(xi, 3) +
         (15.0 / 4.0) * pow(eta, 2) * xi + 2 * eta * pow(xi, 5) -
         4 * eta * pow(xi, 3) + 2 * eta * xi;
}

inline double dN26_deta(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 4) * pow(xi, 4) -
         5.0 / 2.0 * pow(eta, 4) * pow(xi, 2) + (5.0 / 4.0) * pow(eta, 4) +
         pow(eta, 3) * pow(xi, 4) - 2 * pow(eta, 3) * pow(xi, 2) + pow(eta, 3) -
         3.0 / 4.0 * pow(eta, 2) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 2) - 3.0 / 4.0 * pow(eta, 2) -
         1.0 / 2.0 * eta * pow(xi, 4) + eta * pow(xi, 2) - 1.0 / 2.0 * eta;
}

inline double dN27_deta(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 4) * pow(xi, 5) -
         5.0 / 2.0 * pow(eta, 4) * pow(xi, 3) + (5.0 / 4.0) * pow(eta, 4) * xi +
         pow(eta, 3) * pow(xi, 5) - 2 * pow(eta, 3) * pow(xi, 3) +
         pow(eta, 3) * xi - 3.0 / 4.0 * pow(eta, 2) * pow(xi, 5) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 3) - 3.0 / 4.0 * pow(eta, 2) * xi -
         1.0 / 2.0 * eta * pow(xi, 5) + eta * pow(xi, 3) - 1.0 / 2.0 * eta * xi;
}

inline double dN28_deta(double xi, double eta) {
  return 3 * pow(eta, 3) * pow(xi, 5) - 2 * pow(eta, 3) * pow(xi, 4) -
         5 * pow(eta, 3) * pow(xi, 3) + 4 * pow(eta, 3) * pow(xi, 2) -
         3 * eta * pow(xi, 5) + 2 * eta * pow(xi, 4) + 5 * eta * pow(xi, 3) -
         4 * eta * pow(xi, 2);
}

inline double dN29_deta(double xi, double eta) {
  return pow(eta, 3) * pow(xi, 5) - pow(eta, 3) * pow(xi, 4) -
         pow(eta, 3) * pow(xi, 3) + pow(eta, 3) * pow(xi, 2) -
         eta * pow(xi, 5) + eta * pow(xi, 4) + eta * pow(xi, 3) -
         eta * pow(xi, 2);
}

inline double dN30_deta(double xi, double eta) {
  return (15.0 / 4.0) * pow(eta, 4) * pow(xi, 5) -
         5.0 / 2.0 * pow(eta, 4) * pow(xi, 4) -
         25.0 / 4.0 * pow(eta, 4) * pow(xi, 3) + 5 * pow(eta, 4) * pow(xi, 2) -
         9.0 / 2.0 * pow(eta, 2) * pow(xi, 5) + 3 * pow(eta, 2) * pow(xi, 4) +
         (15.0 / 2.0) * pow(eta, 2) * pow(xi, 3) -
         6 * pow(eta, 2) * pow(xi, 2) + (3.0 / 4.0) * pow(xi, 5) -
         1.0 / 2.0 * pow(xi, 4) - 5.0 / 4.0 * pow(xi, 3) + pow(xi, 2);
}

inline double dN31_deta(double xi, double eta) {
  return (5.0 / 4.0) * pow(eta, 4) * pow(xi, 5) -
         5.0 / 4.0 * pow(eta, 4) * pow(xi, 4) -
         5.0 / 4.0 * pow(eta, 4) * pow(xi, 3) +
         (5.0 / 4.0) * pow(eta, 4) * pow(xi, 2) -
         3.0 / 2.0 * pow(eta, 2) * pow(xi, 5) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 4) +
         (3.0 / 2.0) * pow(eta, 2) * pow(xi, 3) -
         3.0 / 2.0 * pow(eta, 2) * pow(xi, 2) + (1.0 / 4.0) * pow(xi, 5) -
         1.0 / 4.0 * pow(xi, 4) - 1.0 / 4.0 * pow(xi, 3) +
         (1.0 / 4.0) * pow(xi, 2);
}

inline double dN32_deta(double xi, double eta) {
  return 4 * pow(eta, 3) * pow(xi, 4) - 8 * pow(eta, 3) * pow(xi, 2) +
         4 * pow(eta, 3) - 4 * eta * pow(xi, 4) + 8 * eta * pow(xi, 2) -
         4 * eta;
}

inline double dN33_deta(double xi, double eta) {
  return 4 * pow(eta, 3) * pow(xi, 5) - 8 * pow(eta, 3) * pow(xi, 3) +
         4 * pow(eta, 3) * xi - 4 * eta * pow(xi, 5) + 8 * eta * pow(xi, 3) -
         4 * eta * xi;
}

inline double dN34_deta(double xi, double eta) {
  return 5 * pow(eta, 4) * pow(xi, 4) - 10 * pow(eta, 4) * pow(xi, 2) +
         5 * pow(eta, 4) - 6 * pow(eta, 2) * pow(xi, 4) +
         12 * pow(eta, 2) * pow(xi, 2) - 6 * pow(eta, 2) + pow(xi, 4) -
         2 * pow(xi, 2) + 1;
}

inline double dN35_deta(double xi, double eta) {
  return 5 * pow(eta, 4) * pow(xi, 5) - 10 * pow(eta, 4) * pow(xi, 3) +
         5 * pow(eta, 4) * xi - 6 * pow(eta, 2) * pow(xi, 5) +
         12 * pow(eta, 2) * pow(xi, 3) - 6 * pow(eta, 2) * xi + pow(xi, 5) -
         2 * pow(xi, 3) + xi;
}

// =============================================================================
// CUBIC HERMITE SHAPE FUNCTIONS (16 DOFs)
// =============================================================================
// Reference: Standard Bicubic Hermite on [-1,1]x[-1,1]
// H0..H3 associated with Node 0, H4..H7 with Node 1, etc.
// Order: u, du/dxi, du/deta, d2u/dxideta at each node.

// 1D Hermite Polynomials
inline double h00(double t) { return (1+2*t)*(1-t)*(1-t); }
inline double h10(double t) { return t*(1-t)*(1-t); }
inline double h01(double t) { return t*t*(3-2*t); }
inline double h11(double t) { return t*t*(t-1); }

inline double dh00(double t) { return 6*t*(t-1); }
inline double dh10(double t) { return (1-t)*(1-3*t); }
inline double dh01(double t) { return 6*t*(1-t); }
inline double dh11(double t) { return t*(3*t-2); }

// Mapping to [-1, 1]: t = (x+1)/2, dt/dx = 1/2
// Let's implement directly for xi, eta in [-1, 1]
// Node 0: (-1, -1) -> Basis functions for Val, dxi, deta, dxideta

// Helper for Bicubic Hermite
// Phi_val = h00(x)*h00(y)
// Phi_dxi = h10(x)*h00(y) * (dx/dxi_map)? 
// Note: Standard definition usually assumes local coords.
// If we use standard 0..1 Hermite:
// H_cubic(u) for u in [0,1]. Map xi in [-1,1] to u = (xi+1)/2.
// Correction factor for derivatives: du/dxi = 0.5.

inline double H_cubic(int type, double xi) {
    double u = (xi + 1.0) * 0.5;
    if (type == 0) return h00(u);       // Value
    if (type == 1) return h10(u) * 0.5; // d/dxi (scaled by du/dxi) -> No wait, h10 is for derivative DOF.
                                        // The shape function for derivative is N_1(u). 
                                        // Effective N for DOF u_xi is: N = h10(u) * (dxi/du)? 
                                        // Let's stick to standard Reference Element definitions.
    if (type == 2) return h01(u);       // Value at u=1
    if (type == 3) return h11(u) * 0.5; // d/dxi at u=1
    return 0.0;
}
// Actually, easier to write explicit formulas for [-1,1]
// N0 = 1/8 * (1-xi)**2 * (2+xi) * (1-eta)**2 * (2+eta) -- Node 0 Val ?
// Let's use the explicit tensor product forms.

// Node 0 (-1, -1)
inline double H0(double xi, double eta) { // Val
    return 0.0625 * (xi-1)*(xi-1)*(2+xi) * (eta-1)*(eta-1)*(2+eta);
}
inline double H1(double xi, double eta) { // dxi
    return 0.0625 * (xi-1)*(xi-1)*(xi+1) * (eta-1)*(eta-1)*(2+eta);
}
inline double H2(double xi, double eta) { // deta
    return 0.0625 * (xi-1)*(xi-1)*(2+xi) * (eta-1)*(eta-1)*(eta+1);
}
inline double H3(double xi, double eta) { // dxideta
    return 0.0625 * (xi-1)*(xi-1)*(xi+1) * (eta-1)*(eta-1)*(eta+1);
}

// Node 1 (1, -1)
inline double H4(double xi, double eta) { // Val
    return 0.0625 * (xi+1)*(xi+1)*(2-xi) * (eta-1)*(eta-1)*(2+eta);
}
inline double H5(double xi, double eta) { // dxi
    return 0.0625 * (xi+1)*(xi+1)*(xi-1) * (eta-1)*(eta-1)*(2+eta);
}
inline double H6(double xi, double eta) { // deta
    return 0.0625 * (xi+1)*(xi+1)*(2-xi) * (eta-1)*(eta-1)*(eta+1);
}
inline double H7(double xi, double eta) { // dxideta
    return 0.0625 * (xi+1)*(xi+1)*(xi-1) * (eta-1)*(eta-1)*(eta+1);
}

// Node 2 (1, 1)
inline double H8(double xi, double eta) { // Val
    return 0.0625 * (xi+1)*(xi+1)*(2-xi) * (eta+1)*(eta+1)*(2-eta);
}
inline double H9(double xi, double eta) { // dxi
    return 0.0625 * (xi+1)*(xi+1)*(xi-1) * (eta+1)*(eta+1)*(2-eta);
}
inline double H10(double xi, double eta) { // deta
    return 0.0625 * (xi+1)*(xi+1)*(2-xi) * (eta+1)*(eta+1)*(eta-1);
}
inline double H11(double xi, double eta) { // dxideta
    return 0.0625 * (xi+1)*(xi+1)*(xi-1) * (eta+1)*(eta+1)*(eta-1);
}

// Node 3 (-1, 1)
inline double H12(double xi, double eta) { // Val
    return 0.0625 * (xi-1)*(xi-1)*(2+xi) * (eta+1)*(eta+1)*(2-eta);
}
inline double H13(double xi, double eta) { // dxi
    return 0.0625 * (xi-1)*(xi-1)*(xi+1) * (eta+1)*(eta+1)*(2-eta);
}
inline double H14(double xi, double eta) { // deta
    return 0.0625 * (xi-1)*(xi-1)*(2+xi) * (eta+1)*(eta+1)*(eta-1);
}
inline double H15(double xi, double eta) { // dxideta
    return 0.0625 * (xi-1)*(xi-1)*(xi+1) * (eta+1)*(eta+1)*(eta-1);
}

// --- DERIVATIVES d/dxi ---
// Node 0
inline double dH0_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi - 3) * (eta-1)*(eta-1)*(2+eta); // Checked: d/dxi ((xi-1)^2(2+xi)) = 3xi^2-3
}
inline double dH1_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi - 2*xi - 1) * (eta-1)*(eta-1)*(2+eta);
}
inline double dH2_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi - 3) * (eta-1)*(eta-1)*(eta+1);
}
inline double dH3_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi - 2*xi - 1) * (eta-1)*(eta-1)*(eta+1);
}

// Node 1
inline double dH4_dxi(double xi, double eta) {
    return 0.0625 * (-3*xi*xi + 3) * (eta-1)*(eta-1)*(2+eta);
}
inline double dH5_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi + 2*xi - 1) * (eta-1)*(eta-1)*(2+eta);
}
inline double dH6_dxi(double xi, double eta) {
    return 0.0625 * (-3*xi*xi + 3) * (eta-1)*(eta-1)*(eta+1);
}
inline double dH7_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi + 2*xi - 1) * (eta-1)*(eta-1)*(eta+1);
}

// Node 2
inline double dH8_dxi(double xi, double eta) {
    return 0.0625 * (-3*xi*xi + 3) * (eta+1)*(eta+1)*(2-eta);
}
inline double dH9_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi + 2*xi - 1) * (eta+1)*(eta+1)*(2-eta);
}
inline double dH10_dxi(double xi, double eta) {
    return 0.0625 * (-3*xi*xi + 3) * (eta+1)*(eta+1)*(eta-1);
}
inline double dH11_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi + 2*xi - 1) * (eta+1)*(eta+1)*(eta-1);
}

// Node 3
inline double dH12_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi - 3) * (eta+1)*(eta+1)*(2-eta);
}
inline double dH13_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi - 2*xi - 1) * (eta+1)*(eta+1)*(2-eta);
}
inline double dH14_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi - 3) * (eta+1)*(eta+1)*(eta-1);
}
inline double dH15_dxi(double xi, double eta) {
    return 0.0625 * (3*xi*xi - 2*xi - 1) * (eta+1)*(eta+1)*(eta-1);
}

// --- DERIVATIVES d/deta ---
// Node 0
inline double dH0_deta(double xi, double eta) {
    return 0.0625 * (xi-1)*(xi-1)*(2+xi) * (3*eta*eta - 3);
}
inline double dH1_deta(double xi, double eta) {
    return 0.0625 * (xi-1)*(xi-1)*(xi+1) * (3*eta*eta - 3);
}
inline double dH2_deta(double xi, double eta) {
    return 0.0625 * (xi-1)*(xi-1)*(2+xi) * (3*eta*eta - 2*eta - 1);
}
inline double dH3_deta(double xi, double eta) {
    return 0.0625 * (xi-1)*(xi-1)*(xi+1) * (3*eta*eta - 2*eta - 1);
}

// Node 1
inline double dH4_deta(double xi, double eta) {
    return 0.0625 * (xi+1)*(xi+1)*(2-xi) * (3*eta*eta - 3);
}
inline double dH5_deta(double xi, double eta) {
    return 0.0625 * (xi+1)*(xi+1)*(xi-1) * (3*eta*eta - 3);
}
inline double dH6_deta(double xi, double eta) {
    return 0.0625 * (xi+1)*(xi+1)*(2-xi) * (3*eta*eta - 2*eta - 1);
}
inline double dH7_deta(double xi, double eta) {
    return 0.0625 * (xi+1)*(xi+1)*(xi-1) * (3*eta*eta - 2*eta - 1);
}

// Node 2
inline double dH8_deta(double xi, double eta) {
    return 0.0625 * (xi+1)*(xi+1)*(2-xi) * (-3*eta*eta + 3);
}
inline double dH9_deta(double xi, double eta) {
    return 0.0625 * (xi+1)*(xi+1)*(xi-1) * (-3*eta*eta + 3);
}
inline double dH10_deta(double xi, double eta) {
    return 0.0625 * (xi+1)*(xi+1)*(2-xi) * (3*eta*eta + 2*eta - 1);
}
inline double dH11_deta(double xi, double eta) {
    return 0.0625 * (xi+1)*(xi+1)*(xi-1) * (3*eta*eta + 2*eta - 1);
}

// Node 3
inline double dH12_deta(double xi, double eta) {
    return 0.0625 * (xi-1)*(xi-1)*(2+xi) * (-3*eta*eta + 3);
}
inline double dH13_deta(double xi, double eta) {
    return 0.0625 * (xi-1)*(xi-1)*(xi+1) * (-3*eta*eta + 3);
}
inline double dH14_deta(double xi, double eta) {
    return 0.0625 * (xi-1)*(xi-1)*(2+xi) * (3*eta*eta + 2*eta - 1);
}
inline double dH15_deta(double xi, double eta) {
    return 0.0625 * (xi-1)*(xi-1)*(xi+1) * (3*eta*eta + 2*eta - 1);
}
