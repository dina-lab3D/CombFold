#include "DotSphere.h"
#include <numerics.h>

int DotSphere::createSphereDots(float radius, float density) {

  float num_equat = 2*pi*radius*sqrt(density);
  float vert_count = 0.5*num_equat;
  
  for(int i=0; i<vert_count; i++) {
    float phi = (pi*i)/vert_count;
    float z = cos(phi);
    float xy = sin(phi);
    float horz_count = xy*num_equat;
    for(int j=0; j<horz_count-1; j++) {
      float teta = (2*pi*j)/horz_count;
      float x = xy*cos(teta);
      float y = xy*sin(teta);
      push_back(Vector3(radius*x,radius*y,radius*z));
    }
  }
  return size();
}
