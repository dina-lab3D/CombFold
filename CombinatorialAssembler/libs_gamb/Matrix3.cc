#include "Matrix3.h"

typedef Matrix3::real real;

Matrix3::Matrix3(const real &d)
{
  vrow[0] = Vector3(d, 0, 0);
  vrow[1] = Vector3(0, d, 0);
  vrow[2] = Vector3(0, 0, d);
}

Matrix3::Matrix3(const real &xr, const real &yr, const real &zr)
{
  real cx = cos(xr);  real cy = cos(yr);  real cz = cos(zr);
  real sx = sin(xr);  real sy = sin(yr);  real sz = sin(zr);
  vrow[0] = Vector3(cz*cy, -sy*sx*cz - sz*cx, -sy*cx*cz + sz*sx);
  vrow[1] = Vector3(sz*cy, -sy*sx*sz + cx*cz, -sy*cx*sz - sx*cz);
  vrow[2] = Vector3(sy,     cy*sx,             cy*cx           );
}

Matrix3::Matrix3(const real &psi, const real &theta, const real &phi, const int &)
{
  real cpsi = cos(psi);  real ctheta = cos(theta);  real cphi = cos(phi);
  real spsi = sin(psi);  real stheta = sin(theta);  real sphi = sin(phi);
  vrow[0] = Vector3( cpsi*cphi - spsi*sphi*ctheta,  cpsi*sphi + spsi*cphi*ctheta,  stheta*spsi);
  vrow[1] = Vector3(-spsi*cphi - cpsi*sphi*ctheta, -spsi*sphi + cpsi*cphi*ctheta,  stheta*cpsi);
  vrow[2] = Vector3(                  sphi*stheta,                  -cphi*stheta,  ctheta);
}

Matrix3::Matrix3(const real a11, const real a12, const real a13,
                 const real a21, const real a22, const real a23,
                 const real a31, const real a32, const real a33) {
  vrow[0] = Vector3(a11, a12, a13);
  vrow[1] = Vector3(a21, a22, a23);
  vrow[2] = Vector3(a31, a32, a33);
}

Matrix3::Matrix3(const real trn[9]) {
  vrow[0] = Vector3(trn);
  vrow[1] = Vector3(trn+3);
  vrow[2] = Vector3(trn+6);
}

Matrix3::Matrix3(const double trn[9]) {
  vrow[0] = Vector3(trn);
  vrow[1] = Vector3(trn+3);
  vrow[2] = Vector3(trn+6);
}

Matrix3::Matrix3(const real trn[3][3]) {
  vrow[0] = Vector3(trn[0]);
  vrow[1] = Vector3(trn[1]);
  vrow[2] = Vector3(trn[2]);
}

Matrix3::Matrix3(const double trn[3][3]) {
  vrow[0] = Vector3(trn[0]);
  vrow[1] = Vector3(trn[1]);
  vrow[2] = Vector3(trn[2]);
}

Matrix3::Matrix3(const Vector3 &diag)
{
  vrow[0] = Vector3(diag[0], 0, 0);
  vrow[1] = Vector3(0, diag[1], 0);
  vrow[2] = Vector3(0, 0, diag[2]);
}

Matrix3::Matrix3(const Vector3 &v1, const Vector3 &v2)
{
  vrow[0] = v1[0] * v2;
  vrow[1] = v1[1] * v2;
  vrow[2] = v1[2] * v2;
}




real Matrix3::determinant() const
{
  return
    element(1,1) * (element(2,2)*element(3,3) - element(2,3)*element(3,2))
  - element(1,2) * (element(2,1)*element(3,3) - element(2,3)*element(3,1))
  + element(1,3) * (element(2,1)*element(3,2) - element(2,2)*element(3,1));
}

Matrix3& Matrix3::operator*=(const real &m)
{
  vrow[0]*= m;
  vrow[1]*= m;
  vrow[2]*= m;
  return *this;
}

Matrix3& Matrix3::operator/=(const real &m)
{
  vrow[0]/= m;
  vrow[1]/= m;
  vrow[2]/= m;
  return *this;
}

Matrix3& Matrix3::operator+=(const Matrix3 &t)
{
  vrow[0]+= t.vrow[0];
  vrow[1]+= t.vrow[1];
  vrow[2]+= t.vrow[2];
  return *this;
}

Matrix3& Matrix3::operator-=(const Matrix3 &t)
{
  vrow[0]-= t.vrow[0];
  vrow[1]-= t.vrow[1];
  vrow[2]-= t.vrow[2];
  return *this;
}

Matrix3& Matrix3::operator*=(const Matrix3 &t)
{
  Vector3 c1 = t.column(1);
  Vector3 c2 = t.column(2);
  Vector3 c3 = t.column(3);
  vrow[0] = Vector3(vrow[0]*c1, vrow[0]*c2, vrow[0]*c3);
  vrow[1] = Vector3(vrow[1]*c1, vrow[1]*c2, vrow[1]*c3);
  vrow[2] = Vector3(vrow[2]*c1, vrow[2]*c2, vrow[2]*c3);
  return *this;
}

Vector3 Matrix3::eigenvalues() const
{
  const real pi23 = (const real)2.0943951023;  // pi*2/3
  real a0 = -determinant();
  real a1 = (element(1,1)*element(2,2) - element(1,2)*element(2,1) +
	     element(1,1)*element(3,3) - element(1,3)*element(3,1) +
	     element(2,2)*element(3,3) - element(2,3)*element(3,2));
  real a2 = -trace();
  real b1 = (3*a1-sqr(a2))/3;
  if (b1 > 0) return Vector3();  // Error! no eigenvalues.
  real b0 = (2*a2*a2*a2 - 9*a2*a1 + 27*a0)/27;
  real d2 = (b1*b1*b1/27 + sqr(b0)/4);
  if (d2 > 0) return Vector3();  // Error! no eignevalues.
  real phi = acos(-(b0/2)*pow(-3/b1, (real)1.5))/3;
  real sqrtb1 = 2*sqrt(-b1/3);
  a2 /= 3;
  return Vector3(sqrtb1*cos(phi       ) - a2,
		 sqrtb1*cos(phi+2*pi23) - a2,
		 sqrtb1*cos(phi+pi23  ) - a2);
}

Matrix3 Matrix3::adjoint() const
{
  return Matrix3(Vector3(element(2,2)*element(3,3) - element(2,3)*element(3,2),
			 element(1,3)*element(3,2) - element(1,2)*element(3,3),
			 element(1,2)*element(2,3) - element(1,3)*element(2,2)
			 ),
		 Vector3(element(2,3)*element(3,1) - element(2,1)*element(3,3),
			 element(1,1)*element(3,3) - element(1,3)*element(3,1),
			 element(1,3)*element(2,1) - element(1,1)*element(2,3)
			 ),
		 Vector3(element(2,1)*element(3,2) - element(2,2)*element(3,1),
			 element(1,2)*element(3,1) - element(1,1)*element(3,2),
			 element(1,1)*element(2,2) - element(1,2)*element(2,1)
			 )
		 );
}

/* if B is the inverse of A then
   Bij = Aji / |A| where Aji is the determinant of A without row i and col j */
Matrix3 operator!(const Matrix3 &t1)
{

  real det = t1.determinant();
  if (det == 0) return Matrix3();
  Matrix3 inv = t1.adjoint();
  inv *= ((real)1/det);
  return inv;
}
