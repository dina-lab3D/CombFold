#include "Rotation3.h"
#include <cmath>
#include "macros.h"

typedef Rotation3::real real;

real Rotation3::matrix11() const {
  return sqr(quat[0]) + sqr(quat[1]) - sqr(quat[2]) - sqr(quat[3]);
}
real Rotation3::matrix12() const {
  return 2*(quat[1]*quat[2] - quat[0]*quat[3]);
}
real Rotation3::matrix13()  const{
  return 2*(quat[2]*quat[3] + quat[0]*quat[2]);
}
real Rotation3::matrix21() const {
  return 2*(quat[1]*quat[2] + quat[0]*quat[3]);
}
real Rotation3::matrix22() const {
  return sqr(quat[0]) - sqr(quat[1]) + sqr(quat[2]) - sqr(quat[3]);
}
real Rotation3::matrix23() const {
  return 2*(quat[2]*quat[3] - quat[0]*quat[1]);
}
real Rotation3::matrix31() const {
  return 2*(quat[1]*quat[3] - quat[0]*quat[2]);
}
real Rotation3::matrix32() const {
  return 2*(quat[2]*quat[3] + quat[0]*quat[1]);
}
real Rotation3::matrix33() const {
  return sqr(quat[0]) - sqr(quat[1]) - sqr(quat[2]) + sqr(quat[3]);
}

Rotation3::Rotation3()
{
  quat[0] = 1;  // quaternion 1, 0, 0, 0 is the identity rotation.
  quat[1] = quat[2] = quat[3] = 0.0;
}


Rotation3::Rotation3(const Matrix3& m)
{
  quat[0] = fabs(1+m[0][0]+m[1][1]+m[2][2])/4;
  quat[1] = fabs(1+m[0][0]-m[1][1]-m[2][2])/4;
  quat[2] = fabs(1-m[0][0]+m[1][1]-m[2][2])/4;
  quat[3] = fabs(1-m[0][0]-m[1][1]+m[2][2])/4;

  // make sure quat is normalized.
  real sum = quat[0]+quat[1]+quat[2]+quat[3];
  quat[0] = sqrt(quat[0]/sum);
  quat[1] = sqrt(quat[1]/sum);
  quat[2] = sqrt(quat[2]/sum);
  quat[3] = sqrt(quat[3]/sum);

  if (m[2][1]-m[1][2] < 0.0) quat[1] = -quat[1];
  if (m[0][2]-m[2][0] < 0.0) quat[2] = -quat[2];
  if (m[1][0]-m[0][1] < 0.0) quat[3] = -quat[3];
}

Rotation3::Rotation3(const real &xr, const real &yr, const real &zr)
{
  real cx = cos(xr);  real cy = cos(yr);  real cz = cos(zr);
  real sx = sin(xr);  real sy = sin(yr);  real sz = sin(zr);
  real m00 = cz*cy;
  real m11 = -sy*sx*sz + cx*cz;
  real m22 = cy*cx;

  quat[0] = (real)(sqrt(1+m00+m11+m22)/2.0);
  quat[1] = (real)(sqrt(1+m00-m11-m22)/2.0);
  quat[2] = (real)(sqrt(1-m00+m11-m22)/2.0);
  quat[3] = (real)(sqrt(1-m00-m11+m22)/2.0);

  if (cy*sx + sy*cx*sz + sx*cz < 0.0) quat[1] = -quat[1];
  if (sz*sx - sy*cx*cz - sy < 0.0)    quat[2] = -quat[2];
  if (sz*cy + sy*sx*cz + sz*cx < 0.0) quat[3] = -quat[3];

  // vrow[0] = Vector3(cz*cy, -sy*sx*cz - sz*cx, -sy*cx*cz + sz*sx);
  // vrow[1] = Vector3(sz*cy, -sy*sx*sz + cx*cz, -sy*cx*sz - sx*cz);
  // vrow[2] = Vector3(sy,     cy*sx,             cy*cx           );
}

Rotation3::Rotation3(const real a, const real b, const real c, const real d)
{
  real sum = sqrt(a*a + b*b + c*c + d*d);
  if (a<0.0) sum = -sum;
  quat[0] = a/sum;
  quat[1] = b/sum;
  quat[2] = c/sum;
  quat[3] = d/sum;
}

Rotation3::Rotation3(const real a, const real b, const real c, const real d,
                     int)
{
  quat[0] = a;
  quat[1] = b;
  quat[2] = c;
  quat[3] = d;
}

Rotation3::Rotation3(const real m[3][3])
{
  quat[0] = fabs(1+m[0][0]+m[1][1]+m[2][2])/4;
  quat[1] = fabs(1+m[0][0]-m[1][1]-m[2][2])/4;
  quat[2] = fabs(1-m[0][0]+m[1][1]-m[2][2])/4;
  quat[3] = fabs(1-m[0][0]-m[1][1]+m[2][2])/4;

  // make sure quat is normalized.
  real sum = quat[0]+quat[1]+quat[2]+quat[3];
  quat[0] = sqrt(quat[0]/sum);
  quat[1] = sqrt(quat[1]/sum);
  quat[2] = sqrt(quat[2]/sum);
  quat[3] = sqrt(quat[3]/sum);

  if (m[2][1]-m[1][2] < 0.0) quat[1] = -quat[1];
  if (m[0][2]-m[2][0] < 0.0) quat[2] = -quat[2];
  if (m[1][0]-m[0][1] < 0.0) quat[3] = -quat[3];
}

Rotation3::Rotation3(const double m[3][3])
{
  quat[0] = (real)fabs(1+m[0][0]+m[1][1]+m[2][2])/4;
  quat[1] = (real)fabs(1+m[0][0]-m[1][1]-m[2][2])/4;
  quat[2] = (real)fabs(1-m[0][0]+m[1][1]-m[2][2])/4;
  quat[3] = (real)fabs(1-m[0][0]-m[1][1]+m[2][2])/4;

  // make sure quat is normalized.
  real sum = quat[0]+quat[1]+quat[2]+quat[3];
  quat[0] = sqrt(quat[0]/sum);
  quat[1] = sqrt(quat[1]/sum);
  quat[2] = sqrt(quat[2]/sum);
  quat[3] = sqrt(quat[3]/sum);

  if (m[2][1]-m[1][2] < 0.0) quat[1] = -quat[1];
  if (m[0][2]-m[2][0] < 0.0) quat[2] = -quat[2];
  if (m[1][0]-m[0][1] < 0.0) quat[3] = -quat[3];
}

Rotation3::Rotation3(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3)
{
  quat[0] = fabs(1+v1[0]+v2[1]+v3[2])/4;
  quat[1] = fabs(1+v1[0]-v2[1]-v3[2])/4;
  quat[2] = fabs(1-v1[0]+v2[1]-v3[2])/4;
  quat[3] = fabs(1-v1[0]-v2[1]+v3[2])/4;

  // make sure quat is normalized.
  real sum = quat[0]+quat[1]+quat[2]+quat[3];
  quat[0] = sqrt(quat[0]/sum);
  quat[1] = sqrt(quat[1]/sum);
  quat[2] = sqrt(quat[2]/sum);
  quat[3] = sqrt(quat[3]/sum);

  if (v3[1]-v2[2] < 0.0) quat[1] = -quat[1];
  if (v1[2]-v3[0] < 0.0) quat[2] = -quat[2];
  if (v2[0]-v1[1] < 0.0) quat[3] = -quat[3];
}

real Rotation3::trace() const
{
  return sqr(quat[0]) - sqr(quat[1]) - sqr(quat[2]) - sqr(quat[3]);
}

const real Rotation3::operator[](const unsigned short r) const
{
  return quat[r];
}

Vector3 Rotation3::row(const unsigned short r) const
{
  switch (r) {
  case 1: return Vector3(matrix11(), matrix12(), matrix13());
  case 2: return Vector3(matrix21(), matrix22(), matrix23());
  case 3: return Vector3(matrix31(), matrix32(), matrix33());
  }
  return Vector3();
}

Vector3 Rotation3::column(const unsigned short col) const
{
  switch (col) {
  case 1: return Vector3(matrix11(), matrix21(), matrix31());
  case 2: return Vector3(matrix12(), matrix22(), matrix32());
  case 3: return Vector3(matrix13(), matrix23(), matrix33());
  }
  return Vector3();
}

real Rotation3::element(const unsigned short r,
                        const unsigned short col) const
{
  switch (r) {
  case 1:
    switch(col) {
    case 1: return matrix11();
    case 2: return matrix12();
    case 3: return matrix13();
    }
  case 2:
    switch(col) {
    case 1: return matrix21();
    case 2: return matrix22();
    case 3: return matrix23();
    }
  case 3:
    switch(col) {
    case 1: return matrix31();
    case 2: return matrix32();
    case 3: return matrix33();
    }
  }
  return 0.0;
}

real Rotation3::rotX() const
{
  return atan2(matrix32(), matrix33());
}

real Rotation3::rotY() const
{
  return atan2(matrix31(), sqrt(sqr(matrix21())+sqr(matrix11())));
}

real Rotation3::rotZ() const
{
  return atan2(matrix21(), matrix11());
}

Vector3 Rotation3::rotationAngles() const
{
  return Vector3(rotX(), rotY(), rotZ());
}

Rotation3::operator Matrix3() const
{
  const real a = quat[0];
  const real b = quat[1];
  const real c = quat[2];
  const real d = quat[3];
  return Matrix3(a*a+b*b-c*c-d*d, 2*(b*c-a*d)    , 2*(b*d+a*c),
                 2*(b*c+a*d)    , a*a-b*b+c*c-d*d, 2*(c*d-a*b),
                 2*(b*d-a*c)    , 2*(c*d+a*b)    , a*a-b*b-c*c+d*d);
}

real Rotation3::dist2(const Rotation3& r) const
{
  return (sqr(quat[0]-r.quat[0]) + sqr(quat[1]-r.quat[1]) +
          sqr(quat[2]-r.quat[2]) + sqr(quat[3]-r.quat[3]));
}

real Rotation3::dist(const Rotation3& r) const
{
  return sqrt(dist2(r));
}

std::pair<Vector3, real> Rotation3::getAxisAndAngle() const
{
  real a, b, c, d;
  a = quat[0];
  b = quat[1];
  c = quat[2];
  d = quat[3];
  if (fabs(a) > .9999) return std::make_pair(Vector3(1, 0, 0), 0.0);
  double angle = acos(a) * 2;
  double s = sin(angle / 2);
  Vector3 axis(b / s, c / s, d / s);
  return std::make_pair(axis.getUnitVector(), angle);
}

Rotation3 operator*(const Rotation3 &q, const Rotation3 &p)
{
  return Rotation3(q[0]*p[0] - q[1]*p[1] - q[2]*p[2] - q[3]*p[3],
                   q[0]*p[1] + q[1]*p[0] + q[2]*p[3] - q[3]*p[2],
                   q[0]*p[2] - q[1]*p[3] + q[2]*p[0] + q[3]*p[1],
                   q[0]*p[3] + q[1]*p[2] - q[2]*p[1] + q[3]*p[0],
                   0);
}


Vector3 operator*(const Rotation3 &t, const Vector3 &v)
{
  const real a = t.quat[0];
  const real b = t.quat[1];
  const real c = t.quat[2];
  const real d = t.quat[3];
  return Vector3((a*a+b*b-c*c-d*d)*v[0] + 2*(b*c-a*d)*v[1] + 2*(b*d+a*c)*v[2],
                 2*(b*c+a*d)*v[0] + (a*a-b*b+c*c-d*d)*v[1] + 2*(c*d-a*b)*v[2],
                 2*(b*d-a*c)*v[0] + 2*(c*d+a*b)*v[1] + (a*a-b*b-c*c+d*d)*v[2]);
}

Rotation3& Rotation3::operator*=(const Rotation3 &p)
{
  real a = quat[0];
  real b = quat[1];
  real c = quat[2];
  real d = quat[3];

  quat[0] = a*p[0] - b*p[1] - c*p[2] - d*p[3];
  quat[1] = a*p[1] + b*p[0] + c*p[3] - d*p[2];
  quat[2] = a*p[2] - b*p[3] + c*p[0] + d*p[1];
  quat[3] = a*p[3] + b*p[2] - c*p[1] + d*p[0];
  return *this;
}

Rotation3 operator!(const Rotation3 &t)
{
  return Rotation3(t.quat[0], -t.quat[1], -t.quat[2], -t.quat[3], 0);
}

std::ostream& operator<<(std::ostream& s, const Rotation3 &t)
{
  return (s << t.quat[0] << ' ' << t.quat[1] << ' ' << t.quat[2] << ' '
          << t.quat[3]);
}
