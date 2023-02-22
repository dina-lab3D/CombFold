#ifndef _Rotation3_h

#define _Rotation3_h

#include <iostream>
#include "numerics.h"
#include "Matrix3.h"

/*
CLASS
  Rotation3

  Defines a rotation of 3-dimensional space.

KEYWORDS
  linear, algebra, matrix, real, vector, rotation, inverse, angle, quaternion

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI>
</UL>

GOALS
  Rotation3 hold a 3-dimensional rotation compactly. Using quaternions a
  rotation may be held compactly without sacrificing much on the performance
  side. Rotation3 should support most of the methods implemented in Matrix3.

USAGE
  Rotation3 hosts a set of constructors, operators and inspection methods.
  Below are just several examples of how one may use the class.
  EXAMPLE
    Rotation3 ra(1.0, 1.1, -1.3);  // rotation angles.
    Rotation3 rb(0.9, -2.3, 0.1);  // rotation angles.

    Matrix3 a= ra;
    Matrix3 b= rb;

    Rotation3 rab = ra*rb;
    Matrix3 ab = a*b;  // rab should be = to ab;

    Vector3 v(1.0, 2.0, 3.0);
    cout !rab*v;
  END
*/
class Rotation3 {

public:
  // GROUP: define real type

  //// real is defined using the real defined in Vector3 class. By default
  // this is of type float. For example , changing Vector3::real's defintion
  // to double will turn the gamb library's vector operations to double
  // precision.
  typedef Vector3::real real;

  // GROUP: Constructors

  //// Default constructor: 0 matrix.
  Rotation3();

  //// Instantiate a rotation using 4 quaternion parameters.
  Rotation3(const real a, const real b, const real c, const real d);

  //// Instantiate a Rotation3 using a rotation matrix. If the matrix is not
  // a rotation matrix the parameters will be normalized.
  explicit Rotation3(const Matrix3& m);

  //// Creates a Rotation3 object. given rotational angles around the
  // x axis (xr), around the y axis (yr) and around the z axis (zr) a
  // linear transformation that rotates space in the order x-y-z is created.
  // Rotation in this order is the default representation for the gamb++
  // library rotational transformations.
  Rotation3(const real &xr, const real &yr, const real &zr);

  //// Use a 3x3 real 2D array matrix to initialize a Rotation3. If matrix is
  // not a rotation it is "normalized".
  explicit Rotation3(const real m[3][3]);

  //// Use a 3x3 double prec. 2D array matrix to initialize a Rotation3. If
  // matrix is not a rotation it is "normalized".
  explicit Rotation3(const double m[3][3]);

  //// Constructs a Rotation3 using the 3 given matrix row vectors. If matrix
  // is not a rotation it is "normalized".
  Rotation3(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3);

  // GROUP: Inspection methods and properties.

  //// Returns the ith element of quaternion representing the rotation. index
  // may take values from 0 to 3. Otherwise results are undefined.
  const real operator[](const unsigned short index) const;

  //// Returns a vector matching the r'th row of the rotation matrix
  // representing the rotation object . r must be in the range 1..3 or run
  // time errors may occur.
  Vector3 row(const unsigned short r) const;

  //// Returns a vector matching the col'th column of the rotation matrix
  // representing the rotation object . r must be in the range 1..3 or run
  // time errors may occur.
  Vector3 column(const unsigned short col) const;

  //// Returns the numerical value at position row, col in the 3X3 matrix
  // representing the rotation.
  // row, col should be within the range of 1..3 or run time errors may occur.
  /// Use the matrix(ij) familiy of methods for faster access.
  real element(const unsigned short r, const unsigned short col) const;

  real matrix11() const;
  real matrix12() const;
  real matrix13() const;
  real matrix21() const;
  real matrix22() const;
  real matrix23() const;
  real matrix31() const;
  real matrix32() const;
  real matrix33() const;

  //// Returns the rotation around the x-axis for a rotational matrix assuming
  // rotations are performed in the order of x-y-z.
  real rotX() const;

  //// Returns the rotation around the y-axis for a rotational matrix assuming
  // rotations are performed in the order of x-y-z.
  real rotY() const;

  //// Returns the rotation around the z-axis for a rotational matrix assuming
  // rotations are performed in the order of x-y-z.
  real rotZ() const;

  //// Returns 3 rotation angles assuming
  // rotations are performed in the order of x-y-z.
  Vector3 rotationAngles() const;

  //// Returns the trace of the rotation matrix representing the rotation.
  real trace() const;

  //// Returns squared distance between 2 rotations. Currently this is
  // implemented by taking the euclidean distance between the 4 parameters of
  // the quternion.
  real dist2(const Rotation3& r) const;

  //// Returns distance between 2 rotations. Currently this is
  // implemented by taking the euclidean distance between the 4 parameters of
  // the quternion.
  real dist(const Rotation3& r) const;

  std::pair<Vector3, real> getAxisAndAngle() const;


  // GROUP: Operators

  //// Convert a Rotation3 object to Matrix3 type.
  operator Matrix3() const;

  //// In place multiplication with a given rotation object.
  Rotation3& operator*=(const Rotation3 &);

  //// Multiplication of two rotations.
  friend Rotation3 operator*(const Rotation3&, const Rotation3&);

  //// Inverse rotation.
  friend Rotation3 operator!(const Rotation3&);

  //// Rotation by vector multiplication.
  friend Vector3 operator*(const Rotation3&, const Vector3&);

  //// Output rotation in quaternion format.
  friend std::ostream& operator<<(std::ostream& s, const Rotation3 &t);

protected:

  //// Instantiate a rotation using 4 quaternion parameters. This version
  // differs from the public version. It dows not normalize the quaternion
  // and does not make sure the 4 parameters actually represent a rotation.
  Rotation3(const real a, const real b, const real c, const real d, int);

private:
  real quat[4];
};

#endif
