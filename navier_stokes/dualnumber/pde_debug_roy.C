#include "ad_masa.h"
#include <iostream>

// typedef double RawScalar;
typedef ShadowNumber<double, long double> RawScalar;

const unsigned int NDIM = 2;

typedef DualNumber<RawScalar, NumberArray<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberArray<NDIM, FirstDerivType> > SecondDerivType;

typedef SecondDerivType ADType;
// typedef FirstDerivType ADType;

int main(void)
{
  const RawScalar xvecinit[] = {1., 0.};

  const NumberArray<NDIM, RawScalar> xvec(xvecinit);

  ADType x = ADType(1., xvec);

  // typedef typename RawType<ADType>::value_type Scalar;
  typedef double Scalar;

  // The identity tensor I
  NumberArray<NDIM, NumberArray<NDIM, Scalar> > Identity = 
    NumberArrayIdentity<NDIM, Scalar>::identity();

  // Euler equation residuals
  NumberArray<NDIM, Scalar> Q_rho_u = 
    raw_value(divergence(x*Identity));

  return 0;

}
