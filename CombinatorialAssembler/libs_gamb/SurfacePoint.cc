#include "SurfacePoint.h"
#include "Shou.h"

#include <algorithm>

typedef SurfacePoint::SurfaceType SurfaceType;

void sortAtoms(int atoms[3])
{
  if (atoms[1] > atoms[0]) std::swap(atoms[0], atoms[1]);
  if (atoms[2] > atoms[1]) std::swap(atoms[1], atoms[2]);
  if (atoms[1] > atoms[0]) std::swap(atoms[0], atoms[1]);
}

SurfacePoint::SurfacePoint(): Vector3(), norm()
{
  atoms[0] = atoms[1] = atoms[2] = -1;
}

SurfacePoint::SurfacePoint(const SurfacePoint& t)
  : Vector3(t), norm(t.norm), area(t.area)
{
  atoms[0] = t.atoms[0];
  atoms[1] = t.atoms[1];
  atoms[2] = t.atoms[2];
}

SurfacePoint::SurfacePoint(const Vector3& newPos, const Vector3& newNorm,
			   const float surface,
                           const int atomIndex1,
                           const int atomIndex2,
                           const int atomIndex3)
  : Vector3(newPos), norm(newNorm), area(surface)
{
  setAtomIndices(atomIndex1, atomIndex2, atomIndex3);
}

SurfacePoint::SurfacePoint(const char* const shouRec)
  : Vector3(Shou::xCoord(shouRec), Shou::yCoord(shouRec),
		     Shou::zCoord(shouRec)),
    norm(Vector3(Shou::xNorm(shouRec), Shou::yNorm(shouRec),
		 Shou::zNorm(shouRec))),
    area(Shou::surface(shouRec))
{
  for (int i=0; i<3 ; ++i)
    atoms[i] = Shou::atomIndex(shouRec, i)-1;
  sortAtoms(atoms);
}

Vector3 SurfacePoint::normal() const
{
  return norm;
}

float SurfacePoint::surfaceArea() const
{
  return area;
}

void SurfacePoint::setAtomIndices(const int atomIndex1,
                                  const int atomIndex2,
                                  const int atomIndex3)
{
  atoms[0] = atomIndex1;
  atoms[1] = atomIndex2;
  atoms[2] = atomIndex3;
  sortAtoms(atoms);
}

int SurfacePoint::atomIndex(const unsigned int indx) const
{
  return atoms[indx];
}

bool SurfacePoint::atomIndicesGTE(const unsigned int atomIndex) const
{
  bool result = (atoms[0] >= (int)atomIndex);
  if (result && (atoms[1] >= 0)) {
    result = (atoms[1] >= (int)atomIndex);
    if (result && (atoms[2] >= 0))
      result = (atoms[2] >= (int)atomIndex);
  }
  return result;
}

bool SurfacePoint::atomIndicesLT(const unsigned int atomIndex) const
{
  return ((atoms[0] < (int)atomIndex) && (atoms[1] < (int)atomIndex) &&
          (atoms[2] < (int)atomIndex));
}

SurfaceType SurfacePoint::surfaceType() const
{
  if (atoms[0] == -1)
    return Undef;
  if (atoms[1] == -1)
    return Caps;
  if (atoms[2] == -1)
    return Belts;
  return Pits;
}

void SurfacePoint::orient(const Vector3 &v)
{
  norm= v/v.norm();
}

bool SurfacePoint::adjacentTo(const SurfacePoint& sp) const
{
  bool result = false;
  switch (surfaceType()) {
  case Caps:
    result =  (atoms[0] == sp.atoms[0] ||
               atoms[0] == sp.atoms[1] ||
               atoms[0] == sp.atoms[2]);
    break;
  case Belts:
    switch (sp.surfaceType()) {
    case Caps:
      result = (atoms[0] == sp.atoms[0] || atoms[1] == sp.atoms[0]);
      break;
    case Belts:
      result = (atoms[0] == sp.atoms[0] && atoms[1] == sp.atoms[1]);
      break;
    case Pits:
      result = ((atoms[0] == sp.atoms[0] && atoms[1] == sp.atoms[1]) ||
                (atoms[0] == sp.atoms[1] && atoms[1] == sp.atoms[2]) ||
                (atoms[0] == sp.atoms[0] && atoms[1] == sp.atoms[2]));
      break;
    case Undef:
      break;
    }
    break;
  case Pits:
    switch (sp.surfaceType()) {
    case Caps:
      result = (atoms[0] == sp.atoms[0] ||
                atoms[1] == sp.atoms[0] ||
                atoms[2] == sp.atoms[0]);
      break;
    case Belts:
      result = ((atoms[0] == sp.atoms[0] && atoms[1] == sp.atoms[1]) ||
                (atoms[1] == sp.atoms[0] && atoms[2] == sp.atoms[1]) ||
                (atoms[0] == sp.atoms[0] && atoms[2] == sp.atoms[1]));
      break;
    case Pits:
      result = (atoms[0] == sp.atoms[0] &&
                atoms[1] == sp.atoms[1] &&
                atoms[2] == sp.atoms[2]);
      break;
    case Undef:
      break;
    }
  case Undef:
    break;
  }
  return result;
}

SurfacePoint& SurfacePoint::operator*=(const Matrix3& lt) {
  (Vector3&)(*this)=(lt*position());
  orient(lt*normal());
  return *this;
}

SurfacePoint& SurfacePoint::operator*=(const RigidTrans3& rt) {
  (Vector3&)(*this)=(rt*position());
  orient(rt.rotation()*normal());
  return *this;
}

std::ostream& operator<<(std::ostream& s, const SurfacePoint& sp)
{
  s.setf(std::ios::fixed, std::ios::floatfield);
  s.setf(std::ios::right, std::ios::adjustfield);
  s.precision(3);

  s.width(5);
  s << sp.atoms[0]+1;
  s.width(5);
  s << sp.atoms[1]+1;
  s.width(5);
  s << sp.atoms[2]+1;
  s.width(8);
  s << sp.x();
  s.width(8);
  s << sp.y();
  s.width(8);
  s << sp.z();
  s.width(8);
  s << sp.area;
  s.width(7);
  s << sp.norm.x();
  s.width(7);
  s << sp.norm.y();
  s.width(7);
  s << sp.norm.z();
  s << std::endl;
  return s;
}
