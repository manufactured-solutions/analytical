#ifndef __raw_type_h__
#define __raw_type_h__


template <typename T>
struct RawType
{
  typedef T value_type;

  static value_type value(const T& a) { return a; }
};

// Make the user syntax slightly nicer
template <typename T>
inline
typename RawType<T>::value_type
raw_value(const T& a) { return RawType<T>::value(a); }

#endif // __raw_type_h__
