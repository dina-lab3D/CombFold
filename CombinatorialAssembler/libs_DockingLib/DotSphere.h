#ifndef DOT_SPHERE_H
#define DOT_SPHERE_H

#include <Vector3.h>

#include <vector>

class DotSphere : public std::vector<Vector3> {
 public:
  // Constructors
  DotSphere(float radius, float density) {
    createSphereDots(radius, density);
  }

  int createSphereDots(float radius, float density);
};

#endif
