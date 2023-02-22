#ifndef _Matrix3_h

#define _Matrix3_h

#include <cmath>
#include <iostream>

#include "macros.h"
#include "numerics.h"
#include "Vector3.h"

/*
CLASS
  Matrix3

  Defines a 3x3 matrix with a load of constructors, operators, and methods.

KEYWORDS
  linear, algebra, matrix, real, vector, rotation, determinant, adjoint,
  inverse, angle, vector

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI> Added support for ZXZ rotation: rotZXZ() function that returns rotational
angles for this type of rotation and constructor that can reconstruct matrix
given rotational angles for ZXZ. 21.7.2003 Dina
<\LI>
</UL>

GOALS
  Matrix3 was written in order to define a 3x3 matrix. The class was written
  so that vector-matrix operations in 3D will be performed with highest
  efficiency. The Matrix3 class can be implemented without any loops since
  the dimension is small and known.

USAGE
  Matrix3 hosts a set of constructors, operators and inspection methods. Below
  are just several examples of how one may use the class.
  EXAMPLE
    Matrix3 A(Vector3(1, 2, 3),   // construct matric out of 3 vector3
              Vector3(1, 0, 1),
	      Vector3(4, 5, 6));
    Matrix3 B=A;
    B.transpose();
    C=A+B;                        // + operator. C is symmetric.
    Vector3 v = C.eigenvalues();  // eigenvalues.
    Vector3 u = C*v;              // multiplication by vector.
    Matrix3 D(v);                 // D has v in the diagonal. 0's elsewhere.
    Matrix3 E(v,u);               // E=v*u'
  END
  A special note should be made regarding the rotational matrix constructors
  and methods. No special rotational matrix class was defined for reasons of
  conciseness and efficiency. Matrix3 hosts constructors for rotational
  matrices and inspectors for inspecting angular parameters. These inspection
  methods are probably meaningless when invoked on non-rotational matrices, but
  will not cause run-time errors and their likes.
  EXAMPLE
    Matrix3 A(pi, pi, pi);   // rotational parameters for x-y-z. We get I.
    Matrix3 B(0,  0,  0 );   // again we get I.
    if (A.rotX() == 0) cout << "ok\n";
  END
*/
class Matrix3 {

public:
  // GROUP: define real type

  //// real is defined using the real defined in Vector3 class. By default
  // this is of type float. For example , changing Vector3::real's defintion
  // to double will turn the gamb library's vector operations to double
  // precision.
  typedef Vector3::real real;

  // GROUP: Constructors

  //// Default constructor: 0 matrix.
  Matrix3() {
    vrow[0] = vrow[1] = vrow[2] = Vector3(0, 0, 0);
  }

  //// Constructs a Matrix3 with d's in the diagonal and 0's else where.
  explicit Matrix3(const real &d);

  //// Creates a rotational Matrix3. given rotational values around the
  // x axis (xr), around the y axis (yr) and around the z axis (zr) a
  // linear transformation that rotates space in the order x-y-z is created.
  // Rotation in this order is the default representation for the gamb++
  // library rotational transformations.
  Matrix3(const real &xr, const real &yr, const real &zr);

  //// Creates a rotational Matrix3. given rotational values around the
  // z axis (psi), around the x axis (theta) and around the z axis (phi), a
  // linear transformation that rotates space in the order z-x-z is created.
  // Rotation in this order is NOT the default representation for the gamb++
  // library rotational transformations. The int zxz specifies that the
  // rotation order is ZXZ and not default XYZ.
  Matrix3(const real &psi, const real &theta, const real &phi, const int &zxz);


  //// Initializes a Matrix3 object by stating all the matrix elements
  // explicitly
  Matrix3(const real a11, const real a12, const real a13,
          const real a21, const real a22, const real a23,
          const real a31, const real a32, const real a33);

  //// Creates a Matrix3 out of a real type array with first 3 elements
  // for the 1st row, 2nd 3 elements for second row and last 3 elements for
  // the 3rd row.
  explicit Matrix3(const real trn[9]);

  //// Creates a Matrix3 out of a double type array with first 3 elements
  // for the 1st row, 2nd 3 elements for second row and last 3 elements for
  // the 3rd row.
  explicit Matrix3(const double trn[9]);

  //// Use a 3x3 real 2D array to initialize matrix.
  explicit Matrix3(const real trn[3][3]);

  //// Use a 3x3 double 2D array to initialize matrix.
  explicit Matrix3(const double trn[3][3]);

  //// Constructs a Matrix3 using the 3 given row vectors. For construction
  // with column vectors use this constructor and then the transpose method.
  Matrix3(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3) {
    vrow[0] = Vector3(v1);
    vrow[1] = Vector3(v2);
    vrow[2] = Vector3(v3);
  }

  //// Constructs a Matrix3 with Vector3 diag's values in the matrix diagonal
  // and 0's elsewhere.
  explicit Matrix3(const Vector3 &diag);

  //// Creates the matrix v1*v2' from the two given Vector3s
  Matrix3(const Vector3 &v1, const Vector3 &v2);

  // GROUP: Inspection methods and properties.

  //// Returns the requested row of the matrix. r must be within 0..2. If
  // not run-time errors may be caused. Notice that a const reference is
  // returned so that this method will perform better than the row method
  // and should be used when const rows are good enough.
  const Vector3& operator[](const unsigned short r) const {
    return vrow[r];
  }

  //// Returns a vector matching the r'th row of the matrix. r must be in
  // the range 1..3 or run time errors may occur.
  Vector3 row(const unsigned short r) const {
    return Vector3(vrow[r-1]);
  }

  //// Returns a vector matching the col'th column of the matrix. col must
  // be in the range of 1..3 or run time errors may occur.
  Vector3 column(const unsigned short col) const {
    return Vector3(vrow[0].x(col), vrow[1].x(col), vrow[2].x(col));
  }

  //// Returns a Vector3 with eigenvalues of Matrix3. If eigenvalues are not
  // real (i.e. complex) the function will return the 0 Vector3. So if the
  // Matrix3 is not a 0-matrix and the function returned a 0-vector the
  // Matrix's eigenvalues are not real.
  Vector3 eigenvalues() const;

  //// Returns the numerical value at position row, col in the 3X3 matrix
  // row, col should be within the range of 1..3 or run time errors may occur.
  real element(const unsigned short r, const unsigned short col) const {
    return vrow[r-1][col-1];
  }

  //// Returns the rotation around the x-axis for a rotational matrix assuming
  // rotations are performed in the order of x-y-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  real rotX() const {
    return atan2(element(3,2), element(3,3));
  }

  //// Returns the rotation around the y-axis for a rotational matrix assuming
  // rotations are performed in the order of x-y-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  real rotY() const {
    return atan2(element(3,1), sqrt(sqr(element(2,1))+sqr(element(1,1))));
  }

  //// Returns the rotation around the z-axis for a rotational matrix assuming
  // rotations are performed in the order of x-y-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  real rotZ() const {
    return atan2(element(2,1), element(1,1));
  }

  //// Returns 3 rotation angles assuming
  // rotations are performed in the order of x-y-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  Vector3 rot() const {
    return Vector3(rotX(), rotY(), rotZ());
  }

  //// Returns 3 rotation angles assuming
  // rotations are performed in the order of z-x-z.
  // Note that results are meaningless if Matrix3 is not a rotational matrix.
  Vector3 rotZXZ() const {
    real psi = atan2(element(1,3), element(2,3));
    real theta = atan2(sqrt(sqr(element(1,3))+sqr(element(2,3))), element(3,3));
    real phi = atan2(element(3,1), -element(3,2));
    return Vector3(psi, theta, phi);
  }

  //// Returns the determinant of the matrix
  real determinant() const;

  //// Returns the trace of the matrix
  real trace() const {
    return element(1,1) + element(2,2) + element(3,3);
  }

  //// Returns the adjoint matrix.
  Matrix3 adjoint() const;

  // GROUP: Operators

  //// Transposes the given matrix.
  void  transpose(){
    real temp12, temp13, temp23;

    temp12 = element(1,2);
    temp13 = element(1,3);
    temp23 = element(2,3);
    vrow[0] = Vector3(element(1,1), element(2,1), element(3,1));
    vrow[1] = Vector3(temp12,       element(2,2), element(3,2));
    vrow[2] = Vector3(temp13,       temp23,       element(3,3));
  }

  //// Returns the matrix transpose.
  Matrix3 transposeMatrix() const{
    return Matrix3(Vector3(vrow[0][0], vrow[1][0], vrow[2][0]),
		   Vector3(vrow[0][1], vrow[1][1], vrow[2][1]),
		   Vector3(vrow[0][2], vrow[1][2], vrow[2][2]));
  }

  //// In place addition of given matrice's members.
  Matrix3& operator+=(const Matrix3&);

  //// In-place subtraction of give matrice's members.
  Matrix3& operator-=(const Matrix3&);

  //// Addition of two matrices.
  friend Matrix3 operator+(const Matrix3& t1, const Matrix3& t2) {
    return Matrix3(t1.vrow[0]+t2.vrow[0], t1.vrow[1]+t2.vrow[1],
		   t1.vrow[2]+t2.vrow[2]);
  }

  //// Substraction of two matrices
  friend Matrix3 operator-(const Matrix3& t1, const Matrix3& t2) {
    return Matrix3(t1.vrow[0]-t2.vrow[0], t1.vrow[1]-t2.vrow[1],
		   t1.vrow[2]-t2.vrow[2]);
  }

  //// In place multiplication with given matrix.
  Matrix3& operator*=(const Matrix3 &);

  //// Matrix multiplication of two matrices.
  friend Matrix3 operator*(const Matrix3& t1, const Matrix3& t2) {
    Vector3 c1 = t2.column(1);
    Vector3 c2 = t2.column(2);
    Vector3 c3 = t2.column(3);
    return Matrix3(Vector3(t1.vrow[0]*c1, t1.vrow[0]*c2, t1.vrow[0]*c3),
                 Vector3(t1.vrow[1]*c1, t1.vrow[1]*c2, t1.vrow[1]*c3),
		   Vector3(t1.vrow[2]*c1, t1.vrow[2]*c2, t1.vrow[2]*c3));
  }

  //// In-place matrix mutiplication by real scalar.
  Matrix3& operator*=(const real &m);

  //// In-place matrix division by real scalar.
  Matrix3& operator/=(const real &m);

  //// Matrix mutiplication by real scalar.
  friend Matrix3 operator*(const Matrix3& t, const real &m) {
    return Matrix3(m*t.vrow[0], m*t.vrow[1], m*t.vrow[2]);
  }

  //// Matrix inverse. If matrix inverse is un-defined a 0-matrix is returned
  friend Matrix3 operator!(const Matrix3&);

  //// Matrix by vector multiplication.
  friend Vector3 operator*(const Matrix3& t, const Vector3& v) {
    return Vector3(t.vrow[0]*v, t.vrow[1]*v, t.vrow[2]*v);
  }

  //// Default Matrix output. 3 row vectors are output using Vector3's <<
  // operator with new-lines between rows.
  friend std::ostream& operator<<(std::ostream& s, const Matrix3 &t) {
    return (s << t.vrow[0] << '\n' << t.vrow[1] << '\n' << t.vrow[2] << '\n');
  }


private:
  Vector3  vrow[3];     /* 3 columns of the matrix are held in 3 vectors */
};

#endif
