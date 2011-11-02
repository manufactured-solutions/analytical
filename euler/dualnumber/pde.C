#include "dualnumberarray.h"
#include <iostream>
#include <masa.h>

typedef double RawScalar;

template <std::size_t NDIM, typename Scalar>
void evaluate_q (const NumberArray<NDIM, Scalar>& xyz);

using namespace MASA;

int main(void)
{
  int err = 0;

  const unsigned int NDIM = 2;

  const RawScalar xvecinit[] = {1., 0.};
  const RawScalar yvecinit[] = {0., 1.};

  const NumberArray<NDIM, RawScalar> xvec(xvecinit);
  const NumberArray<NDIM, RawScalar> yvec(yvecinit);

  typedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;

  typedef SecondDerivType ADType;
  // typedef FirstDerivType ADType;

  // initialize the problem in MASA
  err += masa_init("euler-maple","euler_2d");
  
  // call the sanity check routine
  // (tests that all variables have been initialized)
  err += masa_sanity_check();
  err += masa_printid<double>();

  // set up array
  NumberArray<NDIM, ADType> xy;
  xy[0] = ADType(1., xvec);
  xy[1] = ADType(1., yvec);

  // evaluate source terms
  evaluate_q(xy);

}

// Note: ADScalar needs to be a FirstDerivType or better since we have
// first derivatives here.  Adding diffusion will require a
// SecondDerivType or better

template <std::size_t NDIM, typename ADScalar>
void evaluate_q (const NumberArray<NDIM, ADScalar>& xyz)
{
  typedef typename RawType<ADScalar>::value_type Scalar;

  const Scalar PI = std::acos(Scalar(-1));

  const Scalar k = 1.38;
  const Scalar u_0 = 200.23;
  const Scalar u_x = 1.1;
  const Scalar u_y = 1.08;
  const Scalar v_0 = 1.2;
  const Scalar v_x = 1.6;
  const Scalar v_y = .47;
  const Scalar rho_0 = 100.02;
  const Scalar rho_x = 2.22;
  const Scalar rho_y = 0.8;
  const Scalar p_0 = 150.2;
  const Scalar p_x = .91;
  const Scalar p_y = .623;
  const Scalar a_px = .165;
  const Scalar a_py = .612;
  const Scalar a_rhox = 1.0;
  const Scalar a_rhoy = 1.0;
  const Scalar a_ux = .1987;
  const Scalar a_uy = 1.189;
  const Scalar a_vx = 1.91;
  const Scalar a_vy = 1.0;
  const Scalar Gamma = 1.01;
  const Scalar mu = .918;
  const Scalar L = 3.02;

  const ADScalar& x = xyz[0];
  const ADScalar& y = xyz[1];

  // Treat velocity as a vector
  NumberArray<NDIM, ADScalar> U;

  // Arbitrary manufactured solution
  U[0] = u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L);
  U[1] = v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L);
  ADScalar RHO = rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L);
  ADScalar P = p_0 + p_x * std::cos(a_px * PI * x / L) + p_y * std::sin(a_py * PI * y / L);

  // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // Euler equation residuals
  Scalar Q_rho = raw_value(divergence(U));
  NumberArray<NDIM, Scalar> Q_rho_u = raw_value(divergence(RHO*U.outerproduct(U)) + P.derivatives());
  Scalar Q_rho_e = raw_value(divergence((RHO*ET+P)*U));
  

//  std::cout << f;
}
