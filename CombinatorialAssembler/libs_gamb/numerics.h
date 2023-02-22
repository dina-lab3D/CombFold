#ifndef _numerics_h
#define _numerics_h

#define MAX_FLOAT 99999999.0
#define MIN_FLOAT -99999999.0

const float pi     = (float)3.14159265358979;
const float sqrt_2 = (float)1.4142136;
const float sqrt_3 = (float)1.7320508;
const float sqrt_6 = (float)2.4494897;

// calcualtes c where c = sqrt(a^2 + b^2).
float pythag(float a, float b);

// recieves a 3x3 matrix a and returns its DVDecompisions in a,w and v such
// that original = a*w*vt
bool svd3d(double a[3][3], double w[3], double v[3][3]);

#endif
