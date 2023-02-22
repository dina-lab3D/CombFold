#include <math.h>
#include "macros.h"


// Part of the C-code is from the book NUMERICAL RECIPES in C            
// The routine "svd3d" computes the singular value decomposition of a 3x3  
// matrix a. The result is in the form  a = u w v^t.                      
// The maxtix u replaces a on output. The diagonal matrix w              
// of singular values is output as a vector w[0..2]. The matrix v (not   
// the transposed v^t) is output as v[0..2][0..2].
// if an error occurs svd3d will return true. otherwise it'll return false.
// The only error possible is that after 30 iterations the svd calc. have not
// converged.

double pythag(double a, double b)
{
  double absa, absb;

  absa = fabs(a);
  absb = fabs(b);
  if(absa > absb) 
    return (absa * sqrt(1.0 + sqr(absb/absa)));
  else 
    return (absb == 0.0 ? 0.0 : absb * sqrt(1.0+sqr(absa/absb)));
}
  
bool svd3d(double a[3][3], double w[3], double v[3][3])
{
  int flag, i, its, j, jj, k, l, nm=0;
  double anorm, c, f, g, h, s, scale,  x, y, z;
  double rv1[3];
  
  g = scale = anorm = 0.0;
  
  // Householder reduction to bidiagonal form 

  for (i = 0; i < 3; ++i) {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
    if(i < 3) {
      for(k=i; k<3; ++k) scale += fabs(a[k][i]);
      if(scale) {
	for(k=i; k<3; ++k) {
	  a[k][i] /= scale;
	  s += a[k][i] * a[k][i];
	}
	f = a[i][i];
	g = -sign(sqrt(s),f);
	h = f * g - s;
	a[i][i] = f-g;
	for ( j= l; j < 3; ++j) {
	  for(s= 0.0, k= i; k<3; ++k) s+= a[k][i] * a[k][j];
	  f = s/h;
	  for(k=i; k<3;++k) a[k][j] += f * a[k][i];
	}
	for(k=i;k<3;++k) a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if(i<2) {
      for (k = l; k < 3; ++k) scale += fabs(a[i][k]);
      if(scale) {
	for(k= l; k<3; ++k) {
	  a[i][k] /= scale;
	  s += a[i][k] * a[i][k];
	}
	f = a[i][l];
	g = -sign(sqrt(s),f);
	h= f*g-s;
	a[i][l] = f-g;
	for(k = l; k < 3; ++k) rv1[k] = a[i][k]/h;
	for ( j = l; j < 3; ++j) {
	  for (s = 0.0, k= l; k < 3; ++k) s += a[j][k] * a[i][k];
	  for (k = l; k < 3; ++k) a[j][k] += s * rv1[k];
	}
	for(k = l; k < 3; ++k) a[i][k] *= scale;
      }
    }
    anorm = max2(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }

  // Accumulation of right-hand transformations 

  for (i = 2; i >= 0; --i) {
    if (i < 2) {
      if(g) {
	for (j = l; j < 3; ++j)
	  v[j][i] = (a[i][j]/a[i][l])/g;
	for(j = l; j < 3; ++j) {
	  for (s = 0.0, k = l; k < 3; ++k) s += a[i][k] * v[k][j];
	  for(k = l; k < 3; ++k)  v[k][j] += s * v[k][i];
	}
      }
      for(j = l; j < 3; ++j) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g= rv1[i];
    l = i;
  }

  // Accumulation of left-hand transformations 

  for(i = 2; i >= 0; --i) {
    l = i + 1;
    g = w[i];
    for(j = l; j < 3; ++j) a[i][j] = 0.0;
    if(g) {
      g = 1.0/g;
      for ( j = l; j < 3; ++j) {
	for ( s= 0.0, k = l; k < 3; ++k) s+= a[k][i] * a[k][j];
	f = (s/a[i][i]) * g;
	for( k = i; k < 3; ++k) a[k][j] += f * a[k][i];
      }
      for(j = i; j < 3; ++j) a[j][i] *= g;
    }
    else for (j = i; j < 3; ++j) a[j][i] = 0.0;
    ++a[i][i];
  }

  // Diagonalization of bidiagonal form 

  for(k = 2; k >= 0; --k)   {
    for(its= 1; its <= 30; ++its) {
      flag = 1;
      for(l=k; l >= 0; --l) {
	nm = l-1;
	if((float)(fabs(rv1[l]) + anorm) == (float)anorm) {
	  flag = 0;
	  break;
	}
	if((float)(fabs(w[nm]) + anorm) == (float)anorm) break; 
      }
      if(flag) {
	c = 0.0;
	s = 1.0;
	for(i = l; i <= k; ++i) {
	  f = s * rv1[i];
	  rv1[i] = c * rv1[i];
	  if ((float)(fabs(f) + anorm) == (float)anorm) break;
	  g = w[i];
	  h = pythag(f,g);
	  w[i] = h;
	  h = 1.0/h;
	  c = g * h;
	  s = -f * h;
	  for( j = 0; j < 3; ++j) {
	    y = a[j][nm];
	    z = a[j][i];
	    a[j][nm] = y * c + z * s;
	    a[j][i] = z * c - y * s;
	  }
	}
      }
      z = w[k];
      if( l == k) {
	if ( z < 0.0) {
	  w[k] = -z;
	  for ( j = 0; j < 3; ++j) v[j][k] = (-v[j][k]);
	}
	break;
      }
      if(its == 30) return true;
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y-z) * (y + z) + (g - h) * (g + h))/(2.0*h*y);
      g = pythag(f,1.0);
      f= ((x-z) * (x+z) + h * ((y/(f + sign(g,f))) - h))/x;
      c=s = 1.0;
      for(j = l; j <= nm; ++j) {
	i = j + 1;
	g = rv1[i];
	y = w[i];
	h = s * g;
	g = c * g;
	z = pythag(f,h);
	rv1[j] = z;
	c = f/z;
	s = h/z;
	f = x * c + g * s;
	g = g * c - x * s;
	h = y * s;
	y *= c;
	for ( jj= 0; jj < 3; ++jj) {
	  x = v[jj][j];
	  z = v[jj][i];
	  v[jj][j] = x * c + z * s;
	  v[jj][i] = z * c - x * s;
	}
	z = pythag(f,h);
	w[j] = z;
	if(z) {
	  z = 1.0/z;
	  c = f * z;
	  s = h * z;
	}
	f = (c * g) + (s * y);
	x = (c * y) - (s * g);
	for(jj = 0; jj < 3; ++jj) {
	  y = a[jj][j];
	  z = a[jj][i];
	  a[jj][j] = y * c + z * s;
	  a[jj][i] = z * c - y * s;
	}
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  return false;
}



