#include <cmath>

/**
 * @brief Fermi function for smooth step transitions
 *
 * Returns a value that transitions from 1 to 0 around r = R2
 *
 * @param r Current radial position
 * @param R2 Transition radius
 * @param delta Width of the transition
 * @return double Value between 0 and 1
 *         ~1 when r < R2 (Inner material)
 *         ~0 when r > R2 (Outer material)
 *         0.5 when r = R2
 */
double fermi(double r, double R2, double delta) {
  // To avoid overflow in exp, clamp the argument
  double arg = (r - R2) / delta;

  if (arg > 100.0)
    return 0.0;
  if (arg < -100.0)
    return 1.0;

  return 1.0 / (1.0 + std::exp(arg));
}
