#ifndef _macros_h
#define _macros_h

#include <stdio.h>
#include <string>


//// Template function returns x squared.
template<class T>
inline T sqr(const T x)
{
  return x*x;
}

//// Template function returns absolute value of referenced x.
template<class T>
inline T tabs(const T& x) {
  return (x < 0) ? -x : x;
}

//// Template function returns a with b's sign
template<class T, class S>
inline T sign(const T a, const S b)
{
  return ((b) >= 0.0 ? tabs(a) : -tabs(a));
}

//// Template function returns minimum of two numbers
template<class T>
inline T min2(const T x, const T y){
  return ((x < y)? x : y);
}

////Template function returns maximum of two numbers
template<class T>
inline T max2(const T x, const T y){
  return ((x > y)? x : y);
}

////Template function returning maximum of three numbers
template<class T>
inline T max3(const T x, const T y, const T z){
  if ( x > y )
    return ((x > z) ? x : z);
  else
    return ((y > z) ? y : z);
}

//// Template function returnin minimun of three numbers
template<class T>
inline T min3(const T x, const T y, const T z)
{
  if (x < y)
    return ((x < z) ? x : z);
  else
    return ((y < z) ? y : z);
}

template<class T>
inline T max4(const T x, const T y, const T z,const T t){
  if ( x > y ){
    if(x > z)
      return ((x > t) ? x : t);
    else
      return ((z > t) ? z : t);
  }else{
    if(y > z)
      return ((y > t) ? y : t);
    else
      return ((z > t) ? z : t);
  }
}


#endif
