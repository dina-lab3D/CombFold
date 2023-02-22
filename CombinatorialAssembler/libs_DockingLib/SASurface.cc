#include "SASurface.h"

SASurface::SASurface(ChemMolecule& molecule, float probeRadius, float d) :
  molecule_(molecule), probeRadius_(probeRadius), density_(d), molIntersection_(&molecule_, probeRadius_)
{}

unsigned int SASurface::buildCapsSurface(Surface& surface) {
  surface.reserve(int(molecule_.size() * density_)); //CHECK!!!
  for(unsigned int i=0; i<molecule_.size(); i++) {
    // create dot sphere for an atom
    const float atomRadius = molecule_[i].getRadius();
    DotSphere dotSphere(atomRadius, density_);
    float pointArea = (4 * pi * atomRadius * atomRadius)/dotSphere.size();
    for(DotSphere::iterator it=dotSphere.begin(); it!=dotSphere.end(); it++) {
      // the sphere dots are around (0,0,0) so we have to move them
      Vector3 dotCoord = *it + molecule_(i);
      Vector3 normal = *it/(*it).norm();
      Vector3 probeCenter = dotCoord + normal*probeRadius_;
      if(!molIntersection_.isIntersecting(probeCenter, i)) {
	surface.push_back(SurfacePoint(dotCoord, normal, pointArea, i , -1 , -1));
      }
    }
  }
  return surface.size();
}

float SASurface::computeASAForAtoms() {
  float ASA = 0;
  for(unsigned int i=0; i<molecule_.size(); i++) {
    const float atomRadius = molecule_[i].getRadius();
    // create dot sphere for an atom
    DotSphere dotSphere(atomRadius, density_);
    float pointArea = (4 * pi * atomRadius * atomRadius)/dotSphere.size();
    float atomASA = 0;
    for(DotSphere::iterator it=dotSphere.begin(); it!=dotSphere.end(); it++) {
      // the sphere dots are around (0,0,0) so we have to move them
      Vector3 dotCoord = *it + molecule_(i);
      Vector3 normal = *it/(*it).norm();
      Vector3 probeCenter = dotCoord + normal*probeRadius_;
      if(!molIntersection_.isIntersecting(probeCenter, i)) {
	atomASA += pointArea;
      }
    }
    molecule_[i].setASA(atomASA);
    ASA += atomASA;
  }
  return ASA;
}

// unsigned int SASurface::buildBeltsSurface(Surface& surface) {
//   surface.reserve(int(molecule.size() * density)); //CHECK!!!
//   for(unsigned int i=0; i<molecule.size(); i++) {
//     for(unsigned int j=i+1; j<molecule.size(); j++) {
//       if(molIntersection.isNeighbours(i, j)) {
// 	Torus torus(molecule(i), molecule(j), molecule[i].getRadius(), molecule[j].getRadius(), probeRadius, density);
// 	const Vector3& center = torus.getCenter();
// 	for(Torus::iterator it=torus.begin(); it!=torus.end(); it++) {
// 	  Vector3 normal = (*it-center)/(*it-center).norm();
// 	  Vector3 probeCenter = *it + normal*probeRadius;
// 	  if(!molIntersection.isIntersecting(probeCenter, i, j)) {
// 	    surface.push_back(SurfacePoint(*it, normal, 1.0 , i , j , -1));
// 	  }
// 	}
//       }
//     }
//   }
//   return surface.size();
// }
