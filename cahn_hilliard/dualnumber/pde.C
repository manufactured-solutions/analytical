#include "dualnumberarray.h"
#include <iostream>

typedef double RawScalar;

template <std::size_t NDIM, typename Scalar>
double evaluate_q (const NumberArray<NDIM, Scalar>& xyz);

int main(void)
{
  int N   = 10; // mesh pts. in x and y

  const unsigned int NDIM = 2;

  const RawScalar xvecinit[] = {1., 0.};
  const RawScalar yvecinit[] = {0., 1.};

  const NumberArray<NDIM, RawScalar> xvec(xvecinit);
  const NumberArray<NDIM, RawScalar> yvec(yvecinit);

  typedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberArray<NDIM, SecondDerivType> > ThirdDerivType;
  typedef DualNumber<ThirdDerivType, NumberArray<NDIM, ThirdDerivType> > FourthDerivType;

  typedef FourthDerivType ADType;

  // we first set up the DualNumbers that correspond to independent
  // variables, spatial coordinates x and y.

  NumberArray<NDIM, ADType> xy;

  // When main() says "xy[0] = ADType(1., xvec);", that's saying "x = 1, and 
  // the gradient of f(x,y)=x is the constant vector xvec={1,0}"  
  // Likewise "xy[1] = ADType(1., yvec);" means "y = 1, and the gradient of f(x,y)=y 
  xy[0] = ADType(1., xvec);
  xy[1] = ADType(1., yvec);

  // the input argument xyz is another NumberArray 
  // a vector just like Q_rho_u, a spatial location rather 
  // than a vector-valued forcing function.
  double h = 1.0/N;
  for (int i=0; i != N+1; ++i)
    {
      //
      xy[0] = ADType(i*h, xvec);

      for (int j=0; j != N+1; ++j)
	{
          xy[1] = ADType(j*h, yvec);

	  // AD source terms
	  evaluate_q(xy);
	}
    }
 
  // steady as she goes...
  return 0;

}

// Note: ADScalar needs to be a FirstDerivType or better since we have
// first derivatives here.  Adding diffusion will require a
// SecondDerivType or better

template <std::size_t NDIM, typename ADScalar>
double evaluate_q (const NumberArray<NDIM, ADScalar>& xyz)
{
  typedef typename RawType<ADScalar>::value_type Scalar;

  const ADScalar& x = xyz[0];
  const ADScalar& y = xyz[1];

  Scalar PI = std::acos(Scalar(-1.));

  // Solution parameters
  Scalar c_0 = 0.5;
  Scalar c_x = 0.1;
  Scalar c_y = 0.1;
  Scalar a_cx = 2.0;
  Scalar a_cy = 4.0;
  Scalar L = 1.0;

  // PDE parameters
  Scalar D = 1.0;
  Scalar gamma = 0.0001;

  // Concentration
  ADScalar c;

  // Arbitrary manufactured solution
  c = c_0 + c_x * std::sin(a_cx * PI * x / L) + c_y * std::cos(a_cy * PI * y / L);

  // Cahn-Hilliard equation residuals
  Scalar Q_c =
    raw_value(D*divergence((c*c*c-c-gamma*divergence(c.derivatives())).derivatives()));

  return Q_c;
}
