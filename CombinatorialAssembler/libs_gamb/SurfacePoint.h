#ifndef _SurfacePoint_h
#define _SurfacePoint_h

#include <iostream>
#include "Vector3.h"
#include "Matrix3.h"
#include "RigidTrans3.h"

/*
CLASS
  SurfacePoint

  A descendant of the Vector3 class. Serves as a base class for surface
  points. Besides their positions holds geometrical information such as nornal
  to the surface and type of surface point (caps, pits and belts).

KEYWORDS
  vector, position, algebra, linear, rigid, matrix, distance, match, vector3,
  normal

AUTHOR
  Zipi Fligelman (zipo@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI>
</UL>

GOALS
  Defines a class to hold information about surface points. Besides the
  position inherited from the Vector3 class this class holds the normal to
  surface at the given point and the type of surface point (caps, pits or
  belts). More surface properties may be added including chemical ones.

  The class also holds a parameter which designates which fragment or part of
  the molecule the surface point belongs to. This is useful for hinge-
  flexible molecules.

USAGE
*/

class SurfacePoint : public Vector3
{
public:

  enum SurfaceType { Undef, Caps , Belts , Pits };

  //// Creates a Surface Point with 0-position and a general empty normal
  SurfacePoint();

  //// copy constructor
  SurfacePoint(const SurfacePoint &t);

  //// Creates a point with a given position newPos and a given normal and
  // a given type
  explicit SurfacePoint(const Vector3& newPos, const Vector3& norm,
			const float surface=0,
                        const int atomIndex1 = -1,
                        const int atomIndex2 = -1,
                        const int atomIndex3 = -1);

  //// Construcs a surface point  using a Shou surface format record.
  explicit SurfacePoint(const char* const shouRec);

  //// Returns normal of surface point.
  Vector3 normal() const;

  //// Returns area of surface represented by surface point and normal
  float surfaceArea() const;

  //// Set the atom indices of the surface point. This method allows the class
  // user to set the indices of the atoms from which the surface point was
  // calculated. These indices are read in through the shou file but user may
  // update them at will.
  void setAtomIndices(const int atomIndex1 = -1,
                      const int atomIndex2 = -1,
                      const int atomIndex3 = -1);

  //// Returns atom's indices that participated in the creation of the Shou
  // surface critical points. Up to 3 atoms may define a surface point
  // depeneding on the surface point type (Caps = 1, Belts = 2, Pits = 3)
  // indx may take value of 0..2. If less then 3 atoms defined the surface
  // point then atomIndex will return a -1.
  int atomIndex(const unsigned int indx) const;

  //// Returns true if all atom indices are greater than or equal to the given
  // atom index. Disregards undefined atom indices (i.e. for caps and belts).
  bool atomIndicesGTE(const unsigned int atomIndex) const;

  //// Returns true if all atom indices are less than the given
  // atom index. Disregards undefined atom indices (i.e. for caps and belts).
  bool atomIndicesLT(const unsigned int atomIndex) const;

  //// Returns surface type (caps, pits or belts)
  SurfaceType surfaceType() const;

  //// Changing a point orientation with a new normal
  void orient(const Vector3 &v);

  //// Defines symmetric adjacency relation. This relation stems from a
  // natural triangulation of the protein surface resulting directly from
  // the defintion of the Shou surface. The triangulation is such that every
  // cap point is adjacent to all belts and pits surrounding it, every belt
  // is adjacent to its two pits and two caps neighbors and hence every pit
  // is adjacent to its three caps and three belts in its immediate
  // neighborhood.
  bool adjacentTo(const SurfacePoint& sp) const;

  //// Shift SurfacePoint position and normal using a linear tranfromation
  SurfacePoint& operator*=(const Matrix3 &);

  //// shift SurfacePoint position and normal using a rigid transformation
  SurfacePoint& operator*=(const RigidTrans3 &);

  //// Shift surface point's position using vector addition. normal doesn't
  // change
  //  V.P.  SurfacePoint& operator+=(const Vector3 &);

  //// Shift surface point's position using vector subtraction. normal
  // doesn't chnage
  SurfacePoint& operator-=(const Vector3 &);

  friend std::ostream& operator<<(std::ostream& s, const SurfacePoint& sp);


protected:
  Vector3 norm;

private:
  // the atoms creating the patch filled according to SurfaceType with 1/2/3
  // atoms
  int atoms[3];

  //The area on which the surface point is based
  float area;
};

#endif
