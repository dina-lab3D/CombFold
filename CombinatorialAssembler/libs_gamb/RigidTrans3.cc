#include <math.h>
#include "RigidTrans3.h"

typedef RigidTrans3::real real;


RigidTrans3::RigidTrans3(const Vector3 &o, const Vector3 &u, const Vector3 &v)
{
  trans = (o+u+v)/3.0;
  Vector3 x = trans - o;
  x/= x.norm();
  Vector3 z = (v-o)&(v-u);
  z/= z.norm();
  Vector3 y = z & x;
  rot = Matrix3(x,y,z);
  rot.transpose();
}

RigidTrans3::RigidTrans3(const Vector3 &u ,const Vector3 &w  ,const Vector3 &base , int)
{
  trans=base;
  Vector3 x = u-trans;
  Vector3 z = x&(w-trans);
  Vector3 y = z&x;
  x/= x.norm();
  y/= y.norm();
  z/= z.norm();
  rot = Matrix3(x,y,z);
  rot.transpose();
}

RigidTrans3::RigidTrans3(int , const Vector3 &u ,const Vector3 &w  ,const Vector3 &base)
{
  trans=base;
  Vector3 x = u-trans;
  Vector3 z = x&(w-trans);
  Vector3 y = z&x;
  x/= x.norm();
  y/= y.norm();
  z/= z.norm();
  rot = Matrix3(x,y,z);
  trans = -1*(rot*trans);
}

std::istream& operator>>(std::istream& s, RigidTrans3& rt)
{
  Vector3 rotvec;
  s >> rotvec >> rt.trans;
  rt.rot=Matrix3(rotvec[0], rotvec[1], rotvec[2]);
  return s;
}
