#include "Particle.h"
#include "PDB.h"

Particle::Particle(): pos() {};
Particle::Particle(const Vector3& newPos) : pos(newPos) {};
Particle::Particle(const char* const PDBrec) :
  pos(PDB::atomXCoord(PDBrec), 
      PDB::atomYCoord(PDBrec),
      PDB::atomZCoord(PDBrec)) {}  


Vector3 Particle::position() const { 
  return pos;
}
  
void Particle::moveto(const Vector3 &v) {
  pos=v;
}

Particle::operator Vector3() const {
  return pos;
}

Particle::operator const Vector3&() const
{
  return pos;
}
 					  
Particle& Particle::operator*=(const Matrix3& lt) { 
  pos = lt*pos;
  return *this;
}

Particle& Particle::operator*=(const RigidTrans3& rt) { 
  pos = rt*pos; 
  return *this;
}

Particle& Particle::operator+=(const Vector3 &v) {
  pos+=v;
  return *this;
}

Particle& Particle::operator-=(const Vector3 &v) {
  pos-=v;
  return *this;
}

Particle& Particle::operator*=(const float &m) {
  pos*=m;
  return *this;
}

Particle& Particle::operator/=(const float &m) {
  pos/=m;
  return *this;
}
