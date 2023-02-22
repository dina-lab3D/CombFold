#ifndef MOL_INTERSECTION_H
#define MOL_INTERSECTION_H

#include <Interface.h>
#include <Vector3.h>

/*
  Class for checking sphere intersections with atoms of the molecule
*/
template<class MoleculeT>
class MolIntersection {
 public:
  // GROUP: Constructors
  MolIntersection();
  MolIntersection(MoleculeT* mol, float probe_radius);

  void setMolecule(MoleculeT* mol) { molecule = mol; }

  // GROUP: neighbours related queries
  bool isNeighbours(unsigned int atom1, unsigned int atom2) const;

  const Interface::ParticleAdjacency& getNeighbours(unsigned int atom) const;

  // GROUP: intersections queries
  bool isIntersectingSpheres(const Vector3& sphereCenter1, const Vector3& sphereCenter2,
			     const float radius1, const float radius2) const;

  bool isIntersecting(const Vector3& probeCenter, unsigned int atomIndex) const;

  bool isIntersecting(const Vector3& probeCenter, unsigned int atomIndex1, unsigned int atomIndex2) const;

  bool isIntersecting(const Vector3& probeCenter,
		      unsigned int atom1Index, unsigned int atom2Index, unsigned int atom3Index) const;
 private:
  void computeNeighbours();

 private:
  MoleculeT* molecule;
  Interface neighbours;
  float probeRadius;
};

template<class MoleculeT>
MolIntersection<MoleculeT>::MolIntersection(MoleculeT* mol, float probe_radius) :
  molecule(mol), probeRadius(probe_radius) {
  computeNeighbours();
}

template<class MoleculeT>
bool MolIntersection<MoleculeT>::isNeighbours(unsigned int atom1, unsigned int atom2) const {
  if(neighbours.isAdjacent(atom1, atom2))
    return true;
  return false;
}

template<class MoleculeT>
const Interface::ParticleAdjacency& MolIntersection<MoleculeT>::getNeighbours(unsigned int atom) const {
  return neighbours.adjacencies(atom);
}

template<class MoleculeT>
bool MolIntersection<MoleculeT>::isIntersectingSpheres(const Vector3& sphereCenter1, const Vector3& sphereCenter2,
						       const float radius1, const float radius2) const {
  float radiusSum2 = (radius1+radius2)*(radius1+radius2);
  float dist2 = sphereCenter1.dist2(sphereCenter2);
  if(fabs(radiusSum2-dist2) < 0.0001)
    return false;
  if(radiusSum2 > dist2)
    return true;
  return false;
}

template<class MoleculeT>
bool MolIntersection<MoleculeT>::isIntersecting(const Vector3& probeCenter, unsigned int atomIndex) const {
  const Interface::ParticleAdjacency& atomNeighbours = neighbours.adjacencies(atomIndex);
  for(Interface::ParticleAdjacency::const_iterator it = atomNeighbours.begin(); it != atomNeighbours.end(); it++)
    if(isIntersectingSpheres(probeCenter, (*molecule)(it->first), probeRadius, (*molecule)[it->first].getRadius()))
      return true;
  return false;
}

template<class MoleculeT>
bool MolIntersection<MoleculeT>::isIntersecting(const Vector3& probeCenter,
						unsigned int atom1Index, unsigned int atom2Index) const {
  std::vector<unsigned int> nUnion;
  neighbours.neighboursUnion(neighbours.adjacencies(atom1Index), neighbours.adjacencies(atom2Index), nUnion);
  for(std::vector<unsigned int>::iterator it=nUnion.begin(); it!=nUnion.end(); it++)
    if(isIntersectingSpheres(probeCenter, (*molecule)(*it), probeRadius, (*molecule)[*it].getRadius()))
      return true;
  return false;
    //  return (isIntersecting(probeCenter, atom1Index) || isIntersecting(probeCenter, atom2Index));
}

template<class MoleculeT>
bool MolIntersection<MoleculeT>::isIntersecting(const Vector3& probeCenter, unsigned int atom1Index,
						unsigned int atom2Index, unsigned int atom3Index) const {
  return (isIntersecting(probeCenter, atom1Index) ||
	  isIntersecting(probeCenter, atom2Index) || isIntersecting(probeCenter, atom3Index));
}

template<class MoleculeT>
void MolIntersection<MoleculeT>::computeNeighbours() {
  for(unsigned int i=0; i<molecule->size(); i++) {
    for(unsigned int j=i+1; j<molecule->size(); j++) {
      float radiusSum = (*molecule)[i].getRadius() + (*molecule)[j].getRadius() + 2*probeRadius;
      float dist2 = (*molecule)(i).dist2((*molecule)(j));
      if(dist2 < radiusSum*radiusSum) {
	float dist = sqrt(dist2);
	neighbours.addAdjacency(i, j, dist);
	neighbours.addAdjacency(j, i, dist);
      }
    }
  }
}

#endif
