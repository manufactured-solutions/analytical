
#ifndef __dualnumber_h__
#define __dualnumber_h__

#include <limits> // for numeric_limits
#include <ostream>

#include "compare_types.h"
#include "raw_type.h"

template <typename T, typename D=T>
class DualNumber
{
public:
  typedef T value_type;

  typedef D derivatives_type;

  DualNumber();

  DualNumber(const T& val);

  template <typename T2, typename D2>
  DualNumber(const T2& val, const D2& deriv);

  template <typename T2, typename D2>
  DualNumber(const DualNumber<T2,D2>& src);

  T& value() { return _val; }

  const T& value() const { return _val; }

  D& derivatives() { return _deriv; }

  const D& derivatives() const { return _deriv; }

  DualNumber<T,D> operator- () { return DualNumber<T,D>(-_val, -_deriv); }

  template <typename T2, typename D2>
  DualNumber<T,D>& operator+= (const DualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T,D>& operator+= (const T2& a);

  template <typename T2, typename D2>
  DualNumber<T,D>& operator-= (const DualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T,D>& operator-= (const T2& a);

  template <typename T2, typename D2>
  DualNumber<T,D>& operator*= (const DualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T,D>& operator*= (const T2& a);

  template <typename T2, typename D2>
  DualNumber<T,D>& operator/= (const DualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T,D>& operator/= (const T2& a);


private:
  T _val;
  D _deriv;
};



// Helper class to handle partial specialization for DualNumber
// constructors

template <typename T, typename D>
struct DualNumberConstructor
{
  template <typename T2>
  static T value(const T2& v) { return v; }

  template <typename T2, typename D2>
  static T value(const T2& v, const D2&) { return v; }

  template <typename T2>
  static D deriv(const T2&) { return 0.; }

  template <typename T2, typename D2>
  static D deriv(const DualNumber<T2,D2>& v) { return v.derivatives(); }

  template <typename T2, typename D2>
  static D deriv(const T2&, const D2& d) { return d; }
};

template <typename T, typename D, typename DD>
struct DualNumberConstructor<DualNumber<T,D>, DD>
{
  template <typename T2>
  static DualNumber<T,D> value(const T2& v) { return v; }

  template <typename T2, typename D2>
  static DualNumber<T,D> value(const T2& v, const D2& d) { return DualNumber<T,D>(v,d); }

  template <typename D2>
  static DualNumber<T,D> value(const DualNumber<T,D>& v, const D2&) { return v; }

  template <typename T2>
  static DD deriv(const T2&) { return 0.; }

  template <typename T2, typename D2>
  static DD deriv(const DualNumber<T2,D2>& v) { return v.derivatives(); }

  template <typename T2, typename D2>
  static DD deriv(const T2&, const D2& d) { return d; }
};


//
// Member function definitions
//

template <typename T, typename D>
inline
DualNumber<T,D>::DualNumber() :
  _val(0.), _deriv(0.) {}

template <typename T, typename D>
inline
DualNumber<T,D>::DualNumber(const T& val) :
  _val  (DualNumberConstructor<T,D>::value(val)),
  _deriv(DualNumberConstructor<T,D>::deriv(val)) {}

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D>::DualNumber(const T2& val,
                            const D2& deriv) :
  _val  (DualNumberConstructor<T,D>::value(val,deriv)),
  _deriv(DualNumberConstructor<T,D>::deriv(val,deriv)) {}



template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D>&
DualNumber<T,D>::operator+= (const DualNumber<T2,D2>& a)
{
  this->value() += a.value();
  this->derivatives() += a.derivatives();
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D>&
DualNumber<T,D>::operator-= (const DualNumber<T2,D2>& a)
{
  this->value() -= a.value();
  this->derivatives() -= a.derivatives();
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D>&
DualNumber<T,D>::operator*= (const DualNumber<T2,D2>& a)
{
  this->derivatives() *= a.value();
  this->derivatives() += this->value() * a.derivatives();
  this->value() *= a.value();
  return *this;
}

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D>&
DualNumber<T,D>::operator/= (const DualNumber<T2,D2>& a)
{
  T inva_val = 1. / a.value();
  T this_over_a = this->value() * inva_val;

  this->derivatives() *= inva_val;
  this->derivatives() -= this_over_a * a.derivatives() * inva_val;
  this->value() = this_over_a;
  return *this;
}

template <typename T, typename D>
template <typename T2>
inline
DualNumber<T,D>&
DualNumber<T,D>::operator+= (const T2& a)
{
  this->value() += a;
  return *this;
}

template <typename T, typename D>
template <typename T2>
inline
DualNumber<T,D>&
DualNumber<T,D>::operator-= (const T2& a)
{
  this->value() -= a;
  return *this;
}

template <typename T, typename D>
template <typename T2>
inline
DualNumber<T,D>&
DualNumber<T,D>::operator*= (const T2& a)
{
  this->derivatives() *= a;
  this->value() *= a;
  return *this;
}

template <typename T, typename D>
template <typename T2>
inline
DualNumber<T,D>&
DualNumber<T,D>::operator/= (const T2& a)
{
  T inva_val = 1. / a;

  this->derivatives() *= inva_val;
  this->value() *= inva_val;
  return *this;
}

//
// Non-member functions
//

template <typename T, typename D, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, 
           typename CompareTypes<D,D2>::supertype> 
operator+ (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  typedef typename CompareTypes<D,D2>::supertype DS;
  return DualNumber<TS,DS> 
    (a.value() + b.value(), a.derivatives() + b.derivatives());
}

template <typename T, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, D2>
operator+ (const T& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  return DualNumber<TS,D2>(a + b.value(), b.derivatives());
}

template <typename T, typename D, typename T2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, D>
operator+ (const DualNumber<T,D>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  return DualNumber<TS,D>(a.value() + b, a.derivatives());
}

template <typename T, typename D, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, 
           typename CompareTypes<D,D2>::supertype> 
operator- (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  typedef typename CompareTypes<D,D2>::supertype DS;
  return DualNumber<TS,DS> 
    (a.value() - b.value(), a.derivatives() - b.derivatives());
}

template <typename T, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, D2>
operator- (const T& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  return DualNumber<TS,D2>(a - b.value(), -b.derivatives());
}

template <typename T, typename D, typename T2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, D>
operator- (const DualNumber<T,D>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  return DualNumber<TS,D>(a.value() - b, a.derivatives());
}

template <typename T, typename D, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, 
           typename CompareTypes<D,D2>::supertype> 
operator* (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  typedef typename CompareTypes<D,D2>::supertype DS;
  return DualNumber<TS,DS> 
    (a.value() * b.value(), a.value() * b.derivatives() + b.value() * a.derivatives());
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D2>
>::type
operator* (const T& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  return DualNumber<TS,D2>
    (a * b.value(), a * b.derivatives());
}

template <typename T, typename D, typename T2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D>
>::type
operator* (const DualNumber<T,D>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  return DualNumber<TS,D>
    (a.value() * b, b * a.derivatives());
}

template <typename T, typename D, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, 
           typename CompareTypes<D,D2>::supertype> 
operator/ (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  typedef typename CompareTypes<D,D2>::supertype DS;

  T2 invb_val = 1. / b.value();
  TS a_over_b = a.value() * invb_val;

  return DualNumber<TS,DS> 
    (a_over_b, a.derivatives() * invb_val - a_over_b * b.derivatives() * invb_val);
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D2>
>::type
operator/ (const T& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;

  T2 invb_val = 1. / b.value();
  TS a_over_b = a * invb_val;

  return DualNumber<TS,D2> 
    (a_over_b, - a_over_b * b.derivatives() * invb_val);
}

template <typename T, typename D, typename T2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, D>
operator/ (const DualNumber<T,D>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;

  T2 invb_val = 1. / b;
  TS a_over_b = a.value() * invb_val;

  return DualNumber<TS,D>
    (a_over_b, a.derivatives() * invb_val);
}


namespace std {

// Some forward declarations necessary for recursive DualNumbers

template <typename T, typename D>
DualNumber<T,D> log   (const DualNumber<T,D>& a);

template <typename T, typename D>
DualNumber<T,D> cos   (const DualNumber<T,D>& a);

template <typename T, typename D>
DualNumber<T,D> sqrt  (const DualNumber<T,D>& a);

template <typename T, typename D>
DualNumber<T,D> cosh  (const DualNumber<T,D>& a);

// Now just combined declaration/definitions

template <typename T, typename D, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, 
           typename CompareTypes<D,D2>::supertype> 
pow (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  typedef typename CompareTypes<D,D2>::supertype DS;

  TS a_to_the_b = std::pow(a.value(), b.value());
  return DualNumber<TS,DS>(a_to_the_b,
                           a_to_the_b * (b.value() * a.derivatives() / a.value() +
                                         b.derivatives * log(a.value())));
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D2>
>::type
pow (const T& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;

  TS a_to_the_b = std::pow(a, b.value());

  return DualNumber<TS,D2>(a_to_the_b,
                           a_to_the_b * b.derivatives * log(a));
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D>
>::type
pow (const DualNumber<T,D>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;

  TS a_to_the_b = std::pow(a.value(), b);

  return DualNumber<TS,D>(a_to_the_b,
                          a_to_the_b * (b * a.derivatives() / a.value()));
}

template <typename T, typename D>
inline
DualNumber<T,D> exp   (const DualNumber<T,D>& a)
{
  T expa = std::exp(a.value());
  return DualNumber<T,D>(expa, expa * a.derivatives());
}

template <typename T, typename D>
inline
DualNumber<T,D> log   (const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::log(a.value()), a.derivatives() / a.value());
}

template <typename T, typename D>
inline
DualNumber<T,D> log10 (const DualNumber<T,D>& a)
{
  const static T inv_log_of_10 = 1./std::log(T(10.)); // Require T precision

  return DualNumber<T,D>(std::log10(a.value()), a.derivatives() * inv_log_of_10 / a.value());
}

template <typename T, typename D>
inline
DualNumber<T,D> sin   (const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::sin(a.value()), std::cos(a.value()) * a.derivatives());
}

template <typename T, typename D>
inline
DualNumber<T,D> cos   (const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::cos(a.value()), -std::sin(a.value()) * a.derivatives());
}

template <typename T, typename D>
inline
DualNumber<T,D> tan   (const DualNumber<T,D>& a)
{
  T sec_a = 1. / std::cos(a);
  return DualNumber<T,D>(std::tan(a.value()), sec_a * sec_a * a.derivatives());
}

template <typename T, typename D>
inline
DualNumber<T,D> asin  (const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::asin(a.value()),
                         a.derivatives() / std::sqrt(1 - a.value()*a.value()));
}

template <typename T, typename D>
inline
DualNumber<T,D> acos  (const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::acos(a.value()),
                         -a.derivatives() / std::sqrt(1 - a.value()*a.value()));
}

template <typename T, typename D>
inline
DualNumber<T,D> atan  (const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::atan(a.value()),
                         a.derivatives() /(1 + a.value()*a.value()));
}

template <typename T, typename D, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, 
           typename CompareTypes<D,D2>::supertype> 
atan2 (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  typedef typename CompareTypes<D,D2>::supertype DS;

  return DualNumber<TS,DS>(std::atan2(a.value(), b.value()),
                           (b.value() * a.derivatives() - a.value() * b.derivatives()) /
                           (b.value() * b.value() + a.value() * a.value()));
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D2>
>::type
atan2 (const T& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;

  return DualNumber<TS,D2>(std::atan2(a, b.value()),
                           (-a * b.derivatives()) /
                           (b.value() * b.value() + a * a));
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D>
>::type
atan2 (const DualNumber<T,D>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;

  return DualNumber<TS,D>(std::atan2(a.value(), b),
                          b * a.derivatives() /
                          (b * b + a.value() * a.value()));
}

template <typename T, typename D>
inline
DualNumber<T,D> sinh  (const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::sinh(a.value()), std::cosh(a.value()) * a.derivatives());
}

template <typename T, typename D>
inline
DualNumber<T,D> cosh  (const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::cosh(a.value()), std::sinh(a.value()) * a.derivatives());
}

template <typename T, typename D>
inline
DualNumber<T,D> tanh  (const DualNumber<T,D>& a)
{
  T sech_a = 1./std::cosh(a.value());
  return DualNumber<T,D>(std::tanh(a.value()), sech_a * sech_a * a.derivatives());
}

template <typename T, typename D>
inline
DualNumber<T,D> sqrt  (const DualNumber<T,D>& a)
{
  T sqrt_a = std::sqrt(a.value());

  return DualNumber<T,D>(sqrt_a,
                         a.derivatives() / (2 * sqrt_a));
}

template <typename T, typename D>
inline
DualNumber<T,D> abs  (const DualNumber<T,D>& a)
{
  T sign_a = (a.value() > 0.) - (a.value() < 0.);  // Std < and > return 0 or 1

  return DualNumber<T,D>(std::abs(a.value()),
                         a.derivatives() * sign_a);
}

template <typename T, typename D, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, 
           typename CompareTypes<D,D2>::supertype> 
max (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  bool a_greater = (a.value() > b.value());
  return a_greater ? a : b;
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D2>
>::type
max (const T& a, const DualNumber<T2,D2>& b)
{
  bool a_greater = (a > b.value());
  return a_greater ? a : b;
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D>
>::type
max (const DualNumber<T,D>& a, const T2& b)
{
  bool a_greater = (a.value() > b);
  return a_greater ? a : b;
}

template <typename T, typename D, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, 
           typename CompareTypes<D,D2>::supertype> 
min (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  bool a_greater = (a.value() > b.value());
  return a_greater ? b : a;
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D2>
>::type
min (const T& a, const DualNumber<T2,D2>& b)
{
  bool a_greater = (a > b.value());
  return a_greater ? b : a;
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D>
>::type
min (const DualNumber<T,D>& a, const T2& b)
{
  bool a_greater = (a.value() > b);
  return a_greater ? b : a;
}

template <typename T, typename D>
inline
DualNumber<T,D> ceil (const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::ceil(a.value()), 0.);
}

template <typename T, typename D>
inline
DualNumber<T,D> floor(const DualNumber<T,D>& a)
{
  return DualNumber<T,D>(std::floor(a.value()), 0.);
}

template <typename T, typename D, typename T2, typename D2>
inline
DualNumber<typename CompareTypes<T,T2>::supertype, D>
fmod (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;

  return DualNumber<TS,D>(std::fmod(a.value(), b.value()), a.derivatives());
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
typename CompareTypes<T,T2>::supertype
>::type
fmod (const T& a, const DualNumber<T2,D2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;

  return std::fmod(a, b.value());
}

template <typename T, typename D, typename T2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
DualNumber<typename CompareTypes<T,T2>::supertype, D>
>::type
fmod (const DualNumber<T,D>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;

  return DualNumber<TS,D>(std::fmod(a.value(), b.value()), a.derivatives());
}



template <typename T, typename D>
class numeric_limits<DualNumber<T, D> >
{
public:
  static const bool is_specialized = true;
  static T min() throw() { return DualNumber<T,D>(numeric_limits<T>::min()); }
  static T max() throw() { return DualNumber<T,D>(numeric_limits<T>::max()); }
  static const int  digits = numeric_limits<T>::digits;
  static const int  digits10 = numeric_limits<T>::digits10;
  static const bool is_signed = numeric_limits<T>::is_signed;
  static const bool is_integer = numeric_limits<T>::is_integer;
  static const bool is_exact = numeric_limits<T>::is_exact;
  static const int radix = numeric_limits<T>::radix;
  static T epsilon() throw() {return DualNumber<T,D>(numeric_limits<T>::epsilon()); }
  static T round_error() throw() {return DualNumber<T,D>(numeric_limits<T>::round_error()); }

  static const int  min_exponent = numeric_limits<T>::min_exponent;
  static const int  min_exponent10 = numeric_limits<T>::min_exponent10;
  static const int  max_exponent = numeric_limits<T>::max_exponent;
  static const int  max_exponent10 = numeric_limits<T>::max_exponent10;

  static const bool has_infinity = numeric_limits<T>::has_infinity;
  static const bool has_quiet_NaN = numeric_limits<T>::has_quiet_NaN;
  static const bool has_signaling_NaN = numeric_limits<T>::has_signaling_NaN;
  static const float_denorm_style has_denorm = numeric_limits<T>::has_denorm;
  static const bool has_denorm_loss = numeric_limits<T>::has_denorm_loss;
  static T infinity() throw() {return DualNumber<T,D>(numeric_limits<T>::infinity()); }
  static T quiet_NaN() throw() {return DualNumber<T,D>(numeric_limits<T>::quiet_NaN()); }
  static T signaling_NaN() throw() {return DualNumber<T,D>(numeric_limits<T>::signaling_NaN()); }
  static T denorm_min() throw() {return DualNumber<T,D>(numeric_limits<T>::denorm_min()); }

  static const bool is_iec559 = numeric_limits<T>::is_iec559;
  static const bool is_bounded = numeric_limits<T>::is_bounded;
  static const bool is_modulo = numeric_limits<T>::is_modulo;

  static const bool traps = numeric_limits<T>::traps;
  static const bool tinyness_before = numeric_limits<T>::tinyness_before;
  static const float_round_style round_style = numeric_limits<T>::round_style;
};

} // namespace std

template <typename T, typename D, typename T2, typename D2>
inline
bool
operator<  (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  return (a.value() < b.value());
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
bool
>::type
operator<  (const T& a, const DualNumber<T2,D2>& b)
{
  return (a < b.value());
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
bool
>::type
operator<  (const DualNumber<T,D>& a, const T2& b)
{
  return (a.value() < b.value());
}

template <typename T, typename D, typename T2, typename D2>
inline
bool
operator<=  (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  return (a.value() <= b.value());
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
bool
>::type
operator<=  (const T& a, const DualNumber<T2,D2>& b)
{
  return (a <= b.value());
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
bool
>::type
operator<=  (const DualNumber<T,D>& a, const T2& b)
{
  return (a.value() <= b.value());
}

template <typename T, typename D, typename T2, typename D2>
inline
bool
operator>  (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  return (a.value() > b.value());
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
bool
>::type
operator>  (const T& a, const DualNumber<T2,D2>& b)
{
  return (a > b.value());
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
bool
>::type
operator>  (const DualNumber<T,D>& a, const T2& b)
{
  return (a.value() > b.value());
}

template <typename T, typename D, typename T2, typename D2>
inline
bool
operator>=  (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  return (a.value() >= b.value());
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
bool
>::type
operator>=  (const T& a, const DualNumber<T2,D2>& b)
{
  return (a >= b.value());
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
bool
>::type
operator>=  (const DualNumber<T,D>& a, const T2& b)
{
  return (a.value() >= b.value());
}

template <typename T, typename D, typename T2, typename D2>
inline
bool
operator==  (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  return (a.value() == b.value());
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
bool
>::type
operator==  (const T& a, const DualNumber<T2,D2>& b)
{
  return (a == b.value());
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
bool
>::type
operator==  (const DualNumber<T,D>& a, const T2& b)
{
  return (a.value() == b.value());
}

template <typename T, typename D, typename T2, typename D2>
inline
bool
operator!=  (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b)
{
  return (a.value() != b.value());
}

template <typename T, typename T2, typename D2>
inline
typename boostcopy::enable_if_c<ScalarTraits<T>::value,
bool
>::type
operator!=  (const T& a, const DualNumber<T2,D2>& b)
{
  return (a != b.value());
}

template <typename T, typename T2, typename D>
inline
typename boostcopy::enable_if_c<ScalarTraits<T2>::value,
bool
>::type
operator!=  (const DualNumber<T,D>& a, const T2& b)
{
  return (a.value() != b.value());
}

template <typename T, typename D>
inline
std::ostream&      
operator<< (std::ostream& output, const DualNumber<T,D>& a)
{
  return output << '(' << a.value() << ',' << a.derivatives() << ')';
}


// ScalarTraits, RawType, CompareTypes specializations

template <typename T, typename D>
struct ScalarTraits<DualNumber<T, D> >
{
  static const bool value = true;
};

template <typename T, typename D>
struct RawType<DualNumber<T, D> >
{
  typedef typename RawType<T>::value_type value_type;

  static value_type value(const DualNumber<T, D>& a) { return raw_value(a.value()); }
};

template<typename T, typename T2, typename D>
struct CompareTypes<T, DualNumber<T2, D> > {
  typedef DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<D, T2>::supertype> supertype;
};

template<typename T, typename T2, typename D>
struct CompareTypes<DualNumber<T, D>, T2> {
  typedef DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<D, T2>::supertype> supertype;
};

template<typename T, typename D, typename T2, typename D2>
struct CompareTypes<DualNumber<T, D>, DualNumber<T2, D2> > {
  typedef DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<D, T2>::supertype> supertype;
};

template<typename T, typename D>
struct CompareTypes<DualNumber<T, D>, DualNumber<T, D> > {
  typedef DualNumber<T, D> supertype;
};

#endif // __dualnumber_h__
