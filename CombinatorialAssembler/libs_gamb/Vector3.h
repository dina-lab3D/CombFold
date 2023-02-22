#ifndef _Vector3_h
#define _Vector3_h

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <array>

/*
CLASS
  Vector3

  Defines a 3D vector class with operators.

KEYWORDS
  linear, algebra, vector, angle, distance, norm, real

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI> 23/08/2004 - Oranit Dror:<BR>
Adding the '==' operator and the 'getPerpendicularVector' and 'getUnitVector' methods
</LI>
</UL>

GOALS
  Since most of the work on protein structure is done in 3D space, a fast,
  well designed 3D vector class was in place. Knowing the dimension is 3
  allows for a fast, loopless implemetation of Vector3 and Matrix3 to be
  written.

USAGE
  Many 3D vector operations are supplied including vector addition and
  subtraction, dot-product mult., multiplication by scalar, distance and
  angular operators.<P>

  A new Vector3 may be constructed using 3 floats or an array of 3 floats.
  EXAMPLE
    Vector3 v(1.0, 2.0, 3.0);
    float x[3];
    Vector3 v(x);
  END

  Coordinates may be accessed using the x(int) method which returns the
  vector's given coordinate (1..3) or using the [] operator (0..2).
  EXAMPLE
    if (v[0] == v.x(1))  cout << "OK\n";
  END

  Addition, subtraction, multiplication and division by scalar, dot-product
  distance, angle.
  EXAMPLE
    w = u+v;                           // addition
    w-= v;                             // in-place subtraction w==u
    if ((w|u) == 0) cout << "OK\n";    // distance operator |
    if (w^u == 0) cout << "OK\n";      // angle operator ^
    w=2*w;                             // multiplicatio by scalar
    float f=w*v;                       // dot-product;
  END

*/
class Vector3
{

public:
  // GROUP: define real type.

  //// Defines the type with which coordinates are defined. Changing this to
  // double will casue the Vector3 and Matrix3 classes to work with double
  // precision coordinates. float was chosen for considerations of speed.
  // double precision is rarely needed since strctural data is inexact.
  typedef float real;

  // GROUP: Constructors.

  //// Initialize to the 0 vector.
  Vector3() {
    c[0] = c[1] = c[2] = 0;
  }

  //// Initialize using 3 coordinates x,y,z.
  Vector3(const real& nx, const real& ny, const real& nz) {
    c[0]=nx;
    c[1]=ny;
    c[2]=nz;
  }

  //// initialize using 3 coordinates in array.
  Vector3(const real x[3]) {
    c[0] = x[0];
    c[1] = x[1];
    c[2] = x[2];
  }

  //// initialize using 3 double precision coordinates in array.
  Vector3(const double x[3]) {
    c[0] = (real)x[0];
    c[1] = (real)x[1];
    c[2] = (real)x[2];
  }


  // GROUP: Inspection.

  //// Returns position (redundant function)
  const Vector3& position() const { return *this; }

  //// return the x, y or z coordinates.
  real x() const { return c[0]; }

  real y() const { return c[1]; }

  real z() const { return c[2]; }

  //// Returns x=x(1), y=x(2) or z=x(3) coordinate. If other values are given
  // result are unexpected.
  real x(const unsigned short coord) const {
    return c[coord-1];
  }

  ////
  // Returns true if the vector is the zero vector <br>
  // Author: Oranit Dror (oranit@tau.ac.il)
  bool isZero() const {
    return ((c[0] == 0) && (c[1] == 0) && (c[2] == 0));
  }

  // Group Update
  void updateX(const real& x) {
    c[0]=x;
  }

  void updateY(const real& y) {
    c[1]=y;
  }

  void updateZ(const real& z) {
    c[2]=z;
  }

  void update(const Vector3& v) {
    c[0]=v[0];
    c[1]=v[1];
    c[2]=v[2];
  }

  //// Return Vector3's distance from origin.
  real norm() const {
    return sqrt(norm2());
  }

  //// Return Vector3's squared distance ffrom origin.
  real norm2() const {
    return c[0]*c[0]+c[1]*c[1]+c[2]*c[2];
  }

  //// Returns Vector3's distance from Vector3 p.
  real dist(const Vector3& p) const {
    return sqrt(dist2(p));
  }

  //// Returns Vector3's distance squared from Vector3 p.
  real dist2(const Vector3& p) const{
    real d=c[0]-p.c[0];
    real s=d*d;
    d=c[1]-p.c[1];
    s+=d*d;
    d=c[2]-p.c[2];
    return s+d*d;
  }

  //// Treats vector like a real array. Acceptable indices 0, 1, 2. Other
  // indices may cause run-time errors.
  const real& operator[](const unsigned short coord) const{
    return c[coord];
  }

  //// Access coordinates
  real& operator[](const unsigned short coord){
    return c[coord];
  }

  ////
  // Returns a perpendicular vector
  // (A zero vector is returned in case of a zero vector). <br>
  // Author: Oranit Dror (oranit@tau.ac.il)
  Vector3 getPerpendicularVector() const {
    if (x() != 0) {
      return Vector3((-y()-z())/x(),1,1);
    } else if (y() != 0) {
      return Vector3(1,(-x()-z())/y(),1);
    } else if (z() != 0) {
      return Vector3(1,1,(-x()-y())/z());
    } else {
      return Vector3(0,0,0);
    }
  }

  ////
  // Returns a unit vector in the same direction
  // (A zero vector is returned in case of a zero vector). <br>
  // Author: Oranit Dror (oranit@tau.ac.il)
  Vector3 getUnitVector() const {
    return ((isZero()) ? *this : (*this/norm()));
  }

  // GROUP: Operators.

  ////
  // Returns true if the two vectors have the same coordinates
  // Author: Oranit Dror (oranit@tau.ac.il)
  bool operator==(const Vector3& v) const {
    return ((x() == v.x()) && (y() == v.y()) && (z() == v.z()));
  }

  //// Returns Vector3 addition of two Vector3s.
  friend Vector3 operator+(const Vector3& p1, const Vector3& p2) {
    return Vector3(p1.c[0]+p2.c[0], p1.c[1]+p2.c[1], p1.c[2]+p2.c[2]);
  }

  //// Returns Vector3 subtraction of two Vector3s.
  friend Vector3 operator-(const Vector3& p1, const Vector3& p2) {
    return Vector3(p1.c[0]-p2.c[0], p1.c[1]-p2.c[1], p1.c[2]-p2.c[2]);
  }

  //// Returns Vector3 negative.
  friend Vector3 operator-(const Vector3& p) {
    return Vector3(-p.c[0], -p.c[1], -p.c[2]);
  }

  //// Returns dot product of two Vector3s.
  friend real operator*(const Vector3& p1, const Vector3& p2) {
    return p1.c[0]*p2.c[0]+p1.c[1]*p2.c[1]+p1.c[2]*p2.c[2];
  }

  //// Returns vector multiplied by scalar m.
  friend Vector3 operator*(const Vector3& p, const real& m) {
    return Vector3(p.c[0]*m, p.c[1]*m, p.c[2]*m);
  }


  //// Returns vector multiplied by scalar m.
  friend Vector3 operator*(const real& m, const Vector3& p) {
    return Vector3(p.c[0]*m, p.c[1]*m, p.c[2]*m);
  }

  //// Returns vector divided by scalar m.
  friend Vector3 operator/(const Vector3& p, const real& m) {
    return Vector3(p.c[0]/m, p.c[1]/m, p.c[2]/m);
  }

  //// Returns vector that is vertical to both vectors.
  friend Vector3 operator&(const Vector3& p1, const Vector3& p2) {
    return Vector3(p1.c[1]*p2.c[2]-p1.c[2]*p2.c[1],
		   p1.c[2]*p2.c[0]-p1.c[0]*p2.c[2],
		   p1.c[0]*p2.c[1]-p1.c[1]*p2.c[0]);
  }

  //// Returns distance between 2 Vector3s (similar to dist).
  friend real operator|(const Vector3& p1, const Vector3& p2) {
    return p1.dist(p2);
  }

  //// Returns angle between p1 and p2 in radians.
  friend real operator^(const Vector3& p1, const Vector3& p2) {
    float angle=p1*p2/(p1.norm()*p2.norm());
    if(angle<-1)
      angle=-1;
    else
      if(angle>1)
	angle=1;
    return acos(angle);
  }

  //// Vector3 addition of another Vector3's coordinates.
  Vector3& operator+=(const Vector3& p) {
    c[0]+=p.c[0];  c[1]+=p.c[1];  c[2]+=p.c[2];
    return *this;
  }

  //// Vector3 subtraction of another Vector3's coordinates.
  Vector3& operator-=(const Vector3& p) {
    c[0]-=p.c[0];  c[1]-=p.c[1];  c[2]-=p.c[2];
    return *this;
  }

  Vector3& operator+=(const real& m) {
    c[0]+=m; c[1]+=m; c[2]+=m;
    return *this;
  }

  Vector3& operator-=(const real& m) {
    c[0]-=m; c[1]-=m; c[2]-=m;
    return *this;
  }

  //// Multiplies coordinates by scalar m.
  Vector3& operator*=(const real& m) {
    c[0]*=m;  c[1]*=m;  c[2]*=m;
    return *this;
  }

  //// Divides coordinates by scalar m.
  Vector3& operator/=(const real& m) {
    c[0]/=m;  c[1]/=m;  c[2]/=m;
    return *this;
  }

  //// Outputs coordinates delimited by single space.
  friend std::ostream& operator<<(std::ostream& s, const Vector3& v) {
    return s << v.c[0] << ' ' << v.c[1] << ' ' << v.c[2];
  }


  //// Inputs coordinates delimited by single space.
  friend std::istream& operator>>(std::istream& s, Vector3& v) {
    return s >> v.c[0] >> v.c[1] >> v.c[2];
  }

  ////
  // Outputs the vector's coordinates delimited by a single space. <br>
  // Author: Oranit Dror (oranit@tau.ac.il)
  void output(FILE* outputFile) const {
    fprintf(outputFile, "%f %f %f", c[0], c[1], c[2]);
  }

private:
  std::array<real, 3> c;
};

#endif
