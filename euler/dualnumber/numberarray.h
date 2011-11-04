
#ifndef __numberarray_h__
#define __numberarray_h__

#include <algorithm>

#include "compare_types.h"
#include "raw_type.h"

template <std::size_t size, typename T>
class NumberArray
{
public:
  typedef T value_type;

  NumberArray() {}

  NumberArray(const T& val)
    { std::fill(_data, _data+size, val); }

  NumberArray(const T* vals)
    { std::copy(vals, vals+size, _data); }

  template <typename T2>
  NumberArray(NumberArray<size, T2> src)
    { if (size) std::copy(&src[0], &src[0]+size, _data); }

  T& operator[](unsigned int i)
    { return _data[i]; }

  const T& operator[](unsigned int i) const
    { return _data[i]; }

  NumberArray<size,T> operator- () {
    NumberArray<size,T> returnval;
    for (unsigned int i=0; i != size; ++i) returnval[i] = -_data[i];
    return returnval;
  }

  template <typename T2>
  NumberArray<size,T>& operator+= (const NumberArray<size,T2>& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] += a[i]; }

  template <typename T2>
  NumberArray<size,T>& operator+= (const T2& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] += a; }

  template <typename T2>
  NumberArray<size,T>& operator-= (const NumberArray<size,T2>& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] -= a[i]; }

  template <typename T2>
  NumberArray<size,T>& operator-= (const T2& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] -= a; }

  template <typename T2>
  NumberArray<size,T>& operator*= (const NumberArray<size,T2>& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] *= a[i]; }

  template <typename T2>
  NumberArray<size,T>& operator*= (const T2& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] *= a; }

  template <typename T2>
  NumberArray<size,T>& operator/= (const NumberArray<size,T2>& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] /= a[i]; }

  template <typename T2>
  NumberArray<size,T>& operator/= (const T2& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] /= a; }

  template <typename T2>
  typename CompareTypes<T,T2>::supertype
  dot (const NumberArray<size,T2>& a)
  {
    typename CompareTypes<T,T2>::supertype returnval;
    for (unsigned int i=0; i != size; ++i)
      returnval += _data[i] * a[i];
    return returnval;
  }

  template <typename T2>
  NumberArray<size, NumberArray<size, typename CompareTypes<T,T2>::supertype> >
  outerproduct (const NumberArray<size,T2>& a)
  {
    NumberArray<size, NumberArray<size, typename CompareTypes<T,T2>::supertype> > returnval;

    for (unsigned int i=0; i != size; ++i)
      for (unsigned int j=0; j != size; ++j)
        returnval[i][j] = _data[i] * a[j];

    return returnval;
  }

private:
  T _data[size];
};

//
// Non-member functions
//

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator+ (const NumberArray<size,T>& a, const NumberArray<size,T2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval += b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator+ (const T& a, const NumberArray<size,T2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval += b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator+ (const NumberArray<size,T>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval += b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator- (const NumberArray<size,T>& a, const NumberArray<size,T2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval -= b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator- (const T& a, const NumberArray<size,T2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval -= b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator- (const NumberArray<size,T>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval -= b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator* (const NumberArray<size,T>& a, const NumberArray<size,T2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval *= b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator* (const T& a, const NumberArray<size,T2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval *= b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator* (const NumberArray<size,T>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval *= b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator/ (const NumberArray<size,T>& a, const NumberArray<size,T2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval /= b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator/ (const T& a, const NumberArray<size,T2>& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval /= b;
  return returnval;
}

template <std::size_t size, typename T, typename T2>
NumberArray<size, typename CompareTypes<T,T2>::supertype>
operator/ (const NumberArray<size,T>& a, const T2& b)
{
  typedef typename CompareTypes<T,T2>::supertype TS;
  NumberArray<size, TS> returnval(a);
  returnval /= b;
  return returnval;
}

namespace std {

#define MacroComma ,

#define NumberArray_std_unary(funcname) \
template <std::size_t size, typename T> \
NumberArray<size, T> \
funcname (const NumberArray<size, T>& a) \
{ \
  NumberArray<size, T> returnval; \
 \
  for (unsigned int i=0; i != size; ++i) \
    returnval[i] = std::funcname(a[i]); \
 \
  return returnval; \
}


#define NumberArray_std_binary_abab(funcname, atype, btype, aarg, barg) \
template <std::size_t size, typename T, typename T2> \
NumberArray<size, typename CompareTypes<T,T2>::supertype> \
funcname (const atype& a, const btype& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  NumberArray<size, TS> returnval; \
 \
  for (unsigned int i=0; i != size; ++i) \
    returnval[i] = std::funcname(aarg, barg); \
 \
  return returnval; \
}

#define NumberArray_std_binary(funcname) \
NumberArray_std_binary_abab(funcname, NumberArray<size MacroComma T>, NumberArray<size MacroComma T2>, a[i], b[i]) \
NumberArray_std_binary_abab(funcname,                             T , NumberArray<size MacroComma T2>, a,    b[i]) \
NumberArray_std_binary_abab(funcname, NumberArray<size MacroComma T>,                             T2 , a[i], b)

NumberArray_std_binary(pow)
NumberArray_std_unary(exp)
NumberArray_std_unary(log)
NumberArray_std_unary(log10)
NumberArray_std_unary(sin)
NumberArray_std_unary(cos)
NumberArray_std_unary(tan)
NumberArray_std_unary(asin)
NumberArray_std_unary(acos)
NumberArray_std_unary(atan)
NumberArray_std_binary(atan2)
NumberArray_std_unary(sinh)
NumberArray_std_unary(cosh)
NumberArray_std_unary(tanh)
NumberArray_std_unary(sqrt)
NumberArray_std_unary(abs)
NumberArray_std_binary(max)
NumberArray_std_binary(min)
NumberArray_std_unary(ceil)
NumberArray_std_unary(floor)
NumberArray_std_binary(fmod)



template <std::size_t size, typename T>
class numeric_limits<NumberArray<size, T> >
{
public:
  static const bool is_specialized = true;
  static T min() throw() { return NumberArray<size,T>(numeric_limits<T>::min()); }
  static T max() throw() { return NumberArray<size,T>(numeric_limits<T>::max()); }
  static const int  digits = numeric_limits<T>::digits;
  static const int  digits10 = numeric_limits<T>::digits10;
  static const bool is_signed = numeric_limits<T>::is_signed;
  static const bool is_integer = numeric_limits<T>::is_integer;
  static const bool is_exact = numeric_limits<T>::is_exact;
  static const int radix = numeric_limits<T>::radix;
  static T epsilon() throw() {return NumberArray<size,T>(numeric_limits<T>::epsilon()); }
  static T round_error() throw() {return NumberArray<size,T>(numeric_limits<T>::round_error()); }

  static const int  min_exponent = numeric_limits<T>::min_exponent;
  static const int  min_exponent10 = numeric_limits<T>::min_exponent10;
  static const int  max_exponent = numeric_limits<T>::max_exponent;
  static const int  max_exponent10 = numeric_limits<T>::max_exponent10;

  static const bool has_infinity = numeric_limits<T>::has_infinity;
  static const bool has_quiet_NaN = numeric_limits<T>::has_quiet_NaN;
  static const bool has_signaling_NaN = numeric_limits<T>::has_signaling_NaN;
  static const float_denorm_style has_denorm = numeric_limits<T>::has_denorm;
  static const bool has_denorm_loss = numeric_limits<T>::has_denorm_loss;
  static T infinity() throw() {return NumberArray<size,T>(numeric_limits<T>::infinity()); }
  static T quiet_NaN() throw() {return NumberArray<size,T>(numeric_limits<T>::quiet_NaN()); }
  static T signaling_NaN() throw() {return NumberArray<size,T>(numeric_limits<T>::signaling_NaN()); }
  static T denorm_min() throw() {return NumberArray<size,T>(numeric_limits<T>::denorm_min()); }

  static const bool is_iec559 = numeric_limits<T>::is_iec559;
  static const bool is_bounded = numeric_limits<T>::is_bounded;
  static const bool is_modulo = numeric_limits<T>::is_modulo;

  static const bool traps = numeric_limits<T>::traps;
  static const bool tinyness_before = numeric_limits<T>::tinyness_before;
  static const float_round_style round_style = numeric_limits<T>::round_style;
};

} // namespace std

#define NumberArray_operator_binary_abab(opname, atype, btype, aarg, barg) \
template <std::size_t size, typename T, typename T2> \
NumberArray<size, bool> \
operator opname (const atype& a, const btype& b) \
{ \
  NumberArray<size, bool> returnval; \
 \
  for (unsigned int i=0; i != size; ++i) \
    returnval[i] = (aarg opname barg); \
 \
  return returnval; \
}

#define NumberArray_operator_binary(opname) \
NumberArray_operator_binary_abab(opname, NumberArray<size MacroComma T>, NumberArray<size MacroComma T2>, a[i], b[i]) \
NumberArray_operator_binary_abab(opname,                             T , NumberArray<size MacroComma T2>, a,    b[i]) \
NumberArray_operator_binary_abab(opname, NumberArray<size MacroComma T>,                             T2 , a[i], b)

NumberArray_operator_binary(<)
NumberArray_operator_binary(<=)
NumberArray_operator_binary(>)
NumberArray_operator_binary(>=)
NumberArray_operator_binary(==)
NumberArray_operator_binary(!=)

template <std::size_t size, typename T>
std::ostream&      
operator<< (std::ostream& output, const NumberArray<size,T>& a)
{
  output << '{';
  if (size)
    output << a[0];
  for (unsigned int i=1; i<size; ++i)
    output << ',' << a[i];
  output << '}';
  return output;
}


// CompareTypes, RawType specializations

template<std::size_t size, typename T>
struct CompareTypes<NumberArray<size,T>, NumberArray<size,T> > {
  typedef NumberArray<size, T> supertype;
};


template<std::size_t size, typename T, typename T2>
struct CompareTypes<NumberArray<size,T>, NumberArray<size,T2> > {
  typedef NumberArray<size, typename CompareTypes<T, T2>::supertype> supertype;
};


template <>
template <std::size_t size, typename T>
struct RawType<NumberArray<size, T> >
{
  typedef NumberArray<size, typename RawType<T>::value_type> value_type;

  static value_type value(const NumberArray<size, T>& a)
    {
      value_type returnval;
      for (unsigned int i=0; i != size; ++i)
        returnval[i] = RawType<T>::value(a[i]);
      return returnval;
    }
};

#endif // __numberarray_h__
