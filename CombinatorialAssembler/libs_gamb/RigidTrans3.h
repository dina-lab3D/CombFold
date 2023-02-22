#ifndef _RigidTrans3_h
#define _RigidTrans3_h

#include "Matrix3.h"
#include "Rotation3.h"

#include <iostream>

/*
CLASS
  RigidTrans3

  Defines 3-dimensional rigid transformations class. Rigid transformations
  are comprised of a rotation followed by a translation of the vector space.

KEYWORDS
  rigid, linear, rotation, algebra, translation, vector, matrix

AUTHORS
  Meir Fuchs (meirfux@math.tau.ac.il).

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI>
</UL>

GOALS
  The RigidTrans3 class was designed to represent a rigid transformation of
  3 dimensional space. A rigid transformation is a rotation and then a
  translation of vector space. The class is based upon the Matrix3 class and
  the Vector3 class. Although 3-dimensional rotations can be represented
  more compactly using 3 rotational angles instead of 9 rotational matrix
  elements, the latter was chosen for considerations of speed and performance.

USAGE
  RigidTrans3 is comprised of a rotation and a translation. The basic
  constructor requires a rotation specified by either rotational angles
  or by a rotational matrix and a translation specified by a vector.
  EXAMPLE
    RigidTrans3 T(Vector3(pi, pi, pi/2),  // rotation
                  Vector3(0, 0, 1));      // translation
  END
  A second method of construction is used in reference of base construction.
  Given an ordered set of 3 points one can construct an axes-system based
  on these points. A constructor that recieves 3 points returns the
  RigidTrans3 that transforms the normal axes to the point-based axes system.

  Rigid transformations form a group. We can talk of multiplying or
  applying one rigid transformation on top of the other and we can talk
  about the inverse of a rigid transformation. These operators are supplied.
  EXAMPLE
    RigidTrans3 T(Vector3(pi/2, 0, 0), Vector3(1, 0, 0));
    RigidTrans3 S(Vector3(-pi/2, 0, 0), Vector3(-1, 0, 0));
    RigidTrans3 I=T*S  // This gives the identity rigid trans since T=!S
    RigidTrans3 J=!I   // J should equal I (or be very close to it).
  END
*/
class RigidTrans3
{

public:

  // GROUP: real definition.

  //// real is defined to be float by default. By changing real in Vector3 to
  // double one may change the whole library to work with double precision.
  typedef Vector3::real real;

  // GROUP: Contructors.

  //// Instantiates to a an identity rigid transformation. i.e. no rotation and
  // 0 translation.
  RigidTrans3()
    : rot(1), trans()
  {}

  //// Instantiate a rigid transformation using a linear transformation as
  // rotation and a vector for translation. The int is a dummy int argument
  // used to set this constructor apart from the identically defined public
  // constructor. Note that the matrix itself is not used but rather
  // its rotational parameters, so if a non-rotational matrix is
  // given, a different, rotational Matrix3 will be used in the RigidTrans3 object.
  //Warning: This function is deprecated
  RigidTrans3(const Matrix3& r, const Vector3& t, int)
    : rot(Rotation3(r)), trans(t)
  {}

  //// Instantiate a rigid transformation using a linear transformation as
  // rotation and a vector for translation. This constructor is used
  // to instantiate a RigidTrans3 with a given rotational matrix and
  // translational vector without verifying that the matrix is rotational.
  RigidTrans3(const Matrix3& r, const Vector3& t) :
    rot(r), trans(t)
  {}

  //// Instanitiates a rigid transformation using a Rotation3 object and a
  // translational vector.
  RigidTrans3(const Rotation3& r, const Vector3& t)
    : rot(r), trans(t)
  {}

  //// Instanitiates a rigid transformation using the 3 angles given in the
  // rotational vector and the translational vector. The rotational vector is
  // assumes rotations are in the order: around X (Y,Z plane rotation) -
  // around Y (X,Z plane rotation) - around Z (X,Y plane roation).
  RigidTrans3(const Vector3& r, const Vector3& t)
    : rot(r[0], r[1], r[2]) , trans(t)
  {}

  //// Returns the rigid transformations which transforms the natural axes to
  // the axis system created by u,v,w. The x-axis lies on c-u, (where c is the
  // mass center of u,v,w,). The z-axis is perpendicular to the u,v,w
  // plane and the y-axis is perpendicular to both x and z.
  RigidTrans3(const Vector3& u, const Vector3& v, const Vector3& w);


  //// Returns the rigid transformations which transforms the natural axes to
  // the axis system created by u,w and base (where 'base' is the axis center).
  // The x-axis lies on u-base , the z-axis is
  // perpendicular to the u-base , w-base plane and the y-axis is perpendicular
  // to both x and z.(i.e the first variable is set by default as the local x-axis)
  // The int is a dummy int argument used to set this constructor apart from
  // the identically defined public constructor.
  RigidTrans3(const Vector3& u ,const Vector3& w  ,const Vector3& base , int);


  //// Returns the rigid transformations which transforms the new axis system
  // created by u,w and base (where 'base' is the axis center).
  // (this is just inverse of  RigidTrans3(const Vector3&,const Vector3&,const Vector3&, int))
  // The x-axis lies on u-base , the z-axis is perpendicular to the u-trans ,
  // w-trans plane and the y-axis is perpendicular to both x and z
  // (i.e the first variable is set by default as the local x-axis)
  // The int is a dummy int argument used to set this constructor apart from
  // the identically defined public constructors.
  RigidTrans3(int , const Vector3 &u ,const Vector3 &w  ,const Vector3 &base);

  // GROUP: Inspection and propery methods

  //// Returns linear transformation representing rotation of rigid
  // tranformation.
  const Matrix3& rotation() const {
    return rot;
  }

  //// Returns vector with axis rotational angles assuming the order of
  // rotation is x-y-z.
  Vector3 rotationAngles() const {
    return rot.rot();
  }

  //// Returns translation vector of rigid transformation.
  const Vector3& translation() const {
    return trans;
  }

  // GROUP: Operators.

  //// Apply the rigid transformation rn ob vector v. Apply a rotation followed
  // by a translation
  friend Vector3 operator*(const RigidTrans3 &rt, const Vector3 &v) {
    return (rt.rot * v) + rt.trans;
  }

  //// Apply the rigid transformation rn ob vector v. Apply a rotation followed
  // by a translation
  friend Vector3& operator*=(Vector3 &v,const RigidTrans3 &rt) {
    v=(rt.rot * v) + rt.trans;
    return v;
  }



  //// multiply two rigid transformation such that for any vector v
  // (rt1*rt2)*v = rt1*(rt2*v)
  friend RigidTrans3 operator*(const RigidTrans3 &rt1,
			       const RigidTrans3 &rt2) {
    return RigidTrans3(rt1.rot * rt2.rot,  rt1*(rt2.trans));
  }

  //// Gives the inverse of a rigid transformation such that for any vector v
  // (rt*!rt)*v = rt*(!rt*v) = v and rt*!rt = RigidTrans3().
  friend RigidTrans3 operator!(const RigidTrans3 &rt) {
    Matrix3 inv = rt.rot.transposeMatrix();
    return RigidTrans3(inv, -(inv*(rt.trans)));
  }

  //// prints 3 rotations in X-Y-Z order (see above) and translation vector.
  friend std::ostream& operator<<(std::ostream& s, const RigidTrans3 &rt) {
    return s << rt.rotationAngles() << ' ' << rt.trans;
  }

  //// reads 3 rotations in X-Y-Z order (see above) and translation vector.
  friend std::istream& operator>>(std::istream& s, RigidTrans3& rt);

  ////
  // Prints the three rotation angles in X-Y-Z order and the translation vector <br>
  // Author: Oranit Dror (oranit@tau.ac.il)
  void output(FILE* outputFile) const {
    rotationAngles().output(outputFile);
    fprintf(outputFile, " ");
    trans.output(outputFile);
  }

private:
  Matrix3 rot;      // rotation
  Vector3 trans;    // translation
};

#endif
