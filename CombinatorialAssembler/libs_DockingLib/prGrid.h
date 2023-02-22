#ifndef PR_GRID_H
#define PR_GRID_H

#include <MoleculeGrid.h>
#include <Atom.h>

#include <vector>

class ProbGrid : public MoleculeGrid {
 public:
  ProbGrid(const Surface &surface, const float inDelta, const float maxRadius);

  template<class MoleculeT>
  void markProbs(const MoleculeT &M);

  float getProb(const Vector3 &point) const;
 private:
  std::vector< float > probs;
};

class VolumeGrid : public MoleculeGrid {
 public:
  VolumeGrid(const Surface &surface, const float inDelta, const float maxRadius);
  void computeVolumes(const float radius, const float maxRadius);
  void getVolumeFuncAndNormal(const Vector3 &center, float &volFunc, Vector3& volNormal) const;
  float getVolumeFunc(const Vector3 &center) const;
  const Vector3& getVolumeNormal(const Vector3 &center) const;
 private:
  std::vector< Vector3 > volumeNormals;
  std::vector< float > volumeFunctions;
};

class ResidueGrid : public MoleculeGrid {
 public:
  ResidueGrid(const Surface &surface, const float inDelta, const float maxRadius);

  template<class MoleculeT>
  void markResidues(const MoleculeT &M);
  int getResidueEntry(const Vector3 &point) const;
  void printGrid(std::ofstream& file) const;
 protected:
  std::vector< int > residues;
};

class AtomsGrid : public MoleculeGrid {
 public:
  AtomsGrid(const Surface &surface, const float inDelta, const float maxRadius);

  template<class MoleculeT>
  void markAtoms(const MoleculeT &M, float radius);
  const std::vector<unsigned int>& getAtoms(const Vector3 &point) const;
  unsigned int getAtomsNumber(const Vector3 &point) const;
 private:
  std::vector<std::vector<unsigned int > > atoms;
};

// stores closest atom number for every grid point
class AtomGrid : public MoleculeGrid {
 public:
  AtomGrid(const Surface &surface, const float inDelta, const float maxRadius):
    MoleculeGrid(surface, inDelta, maxRadius), atoms_(maxEntry) {}
  template<class MoleculeT>
  void markAtoms(const MoleculeT &mol, std::vector<float>& radii);

  // get volume corresponding to each atom
  template<class MoleculeT>
  void getVolume(const MoleculeT &mol, std::vector<float>& volume);

  unsigned int getAtom(const Vector3 &point) const;

 private:
  std::vector<unsigned int> atoms_;
};

class AtomTypeGrid : public MoleculeGrid {
 public:
  AtomTypeGrid(const Surface &surface, const float inDelta, const float maxRadius);

  template<class MoleculeT>
  void markAtomTypes(const MoleculeT &M);
  unsigned int getAtomType(const Vector3 &point) const {
     int index = getIndexForPoint(point);
     if(!isValidIndex(index))
       return 0;
     return atomTypes[index];
  }
 private:
  std::vector<unsigned int > atomTypes;
};


class SurfacePointGrid : public MoleculeGrid {
 public:
  SurfacePointGrid(const Surface &surface, const float inDelta, const float maxRadius);
  void computeDistAndPointers(const Surface &surface, bool precise = false);

  const SurfacePoint* getSpPointer(const Vector3& point) const {
    int index = getIndexForPoint(point);
    if (!isValidIndex(index))
      return NULL;
    return spPointers[index];
  }

  bool getDistAndSpPointer(const Vector3& point, float& distance, SurfacePoint* spPoint) const {
    int index = getIndexForPoint(point);
    if (!isValidIndex(index))
      return false;
    distance=grid[index];
    spPoint = (SurfacePoint *)spPointers[index];
    return true;
  }

  float getExactDist2(const Vector3& point) const {
    int index = getIndexForPoint(point);
    if(!isValidIndex(index) || spPointers[index]==NULL)
      return MAX_FLOAT*MAX_FLOAT;
    return point.dist2(spPointers[index]->position());
  }

  float getExactDist(const Vector3& point) const {
    int index = getIndexForPoint(point);
    if(!isValidIndex(index))
      return MAX_FLOAT;
    if(spPointers[index]==NULL)
      return grid[index];
    if(grid[index] < 0)
      return -(point.dist(spPointers[index]->position()));
    else
      return point.dist(spPointers[index]->position());
  }

 protected:
  // sets distance func. falue for index to be the distance from p.
  // and p to be the closest point
  bool setDistAndPointer(int index, const SurfacePoint* p) {
    float d = getPointForIndex(index).dist(p->position());
    if (d >= grid[index]) return false;
    grid[index] = d;
    spPointers[index] = p;
    return true;
  }
 private:
  std::vector<const SurfacePoint*> spPointers;
};

template<class MoleculeT>
void ProbGrid::markProbs(const MoleculeT &M) {
  std::vector<float> weights;
  weights.insert(weights.end(), maxEntry, 0.0);
  for(typename MoleculeT::const_iterator it=M.begin(); it!=M.end(); it++) {
    float atomRadius = it->getAtomRadius();
    float radius = atomRadius*2; //may be +1 is enough

    int centerIndex = getIndexForPoint(it->position());
    if (!isValidIndex(centerIndex)) {
      std::cerr << "Error: Point out of grid" << std::endl;
      exit(1);
    }

    int intRadius = getIntGridRadius(radius);
    int radius2 = intRadius*intRadius;

    int i_bound,j_bound,k_bound;
    i_bound = intRadius;
    for (int i = -i_bound; i <= i_bound; i++ ) {
      j_bound = (int)sqrt(radius2 - i*i);
      for (int j = -j_bound; j <= j_bound; j++) {
	k_bound = (int)sqrt(radius2 - i*i - j*j);
	for (int k = -k_bound; k <= k_bound; k++) {
	  int index = centerIndex + i + xGridNum * j + xyGridNum * k;
	  if(isValidIndex(index) && grid[index] <= 0) {
	    Vector3 point = getPointForIndex(index);
	    float dist = point.dist(it->position());
	    if(dist == 0.0) {
	      probs[index]=it->getSteadyProb();
	      continue; }
	    float weight = atomRadius/dist;
	    if(weight <= weights[index])
	      continue;
	    weights[index]=weight;
	    probs[index]=it->getSteadyProb();
	  }
	}
      }
    }
  }
  weights.clear();
}

template<class MoleculeT>
void ResidueGrid::markResidues(const MoleculeT &M) {
  std::vector<float> weights;
  weights.insert(weights.end(), maxEntry, 0.0);
  for(typename MoleculeT::const_iterator it=M.begin(); it!=M.end(); it++) {
    float atomRadius = it->getRadius();
    float radius = atomRadius*2; //may be +1 is enough

    int centerIndex = getIndexForPoint(it->position());
    if (!isValidIndex(centerIndex)) {
      std::cerr << "Error: Point out of grid" << std::endl;
      exit(1);
    }

    int intRadius = getIntGridRadius(radius);
    int radius2 = intRadius*intRadius;

    int i_bound,j_bound,k_bound;
    i_bound = intRadius;
    for (int i = -i_bound; i <= i_bound; i++ ) {
      j_bound = (int)sqrt(radius2 - i*i);
      for (int j = -j_bound; j <= j_bound; j++) {
	k_bound = (int)sqrt(radius2 - i*i - j*j);
	for (int k = -k_bound; k <= k_bound; k++) {
	  int index = centerIndex + i + xGridNum * j + xyGridNum * k;
	  if(isValidIndex(index) && grid[index] <= 0) {
	    Vector3 point = getPointForIndex(index);
	    float dist = point.dist(it->position());
	    if(dist == 0.0) {
	      if(it->isBackbone())
		residues[index] = M.residueEntry(it->chainId(), it->residueSequenceID()) * -1;
	      else
		residues[index] = M.residueEntry(it->chainId(), it->residueSequenceID());
	      continue;
	    }
	    float weight = atomRadius/dist;
	    if(weight <= weights[index])
	      continue;
	    weights[index] = weight;
	    if(it->isBackbone())
	      residues[index] = M.residueEntry(it->chainId(), it->residueSequenceID()) * -1;
	    else
	      residues[index] = M.residueEntry(it->chainId(), it->residueSequenceID());
	  }
	}
      }
    }
  }
  weights.clear();
}

template<class MoleculeT>
void AtomsGrid::markAtoms(const MoleculeT &M, float radius) {
  int atomIndex=0;
  int intRadius = getIntGridRadius(radius);
  int radius2 = intRadius*intRadius;

  for(typename MoleculeT::const_iterator it=M.begin(); it!=M.end(); it++, atomIndex++) {
    int centerIndex = getIndexForPoint(it->position());
    if (!isValidIndex(centerIndex)) {
      std::cerr << "Error: Point out of grid" << std::endl;
      exit(1);
    }

    int i_bound,j_bound,k_bound;
    i_bound = intRadius;
    for (int i = -i_bound; i <= i_bound; i++ ) {
      j_bound = (int)sqrt(radius2 - i*i);
      for (int j = -j_bound; j <= j_bound; j++) {
	k_bound = (int)sqrt(radius2 - i*i - j*j);
	for (int k = -k_bound; k <= k_bound; k++) {
	  int index = centerIndex + i + xGridNum * j + xyGridNum * k;
	  if(isValidIndex(index)) {
	    atoms[index].push_back(atomIndex);
	  }
	}
      }
    }
  }
}

template<class MoleculeT>
void AtomGrid::markAtoms(const MoleculeT &mol, std::vector<float>& radii) {
  std::vector<float> weights;
  weights.insert(weights.end(), maxEntry, 0.0);
  int atomIndex=0;
  for(typename MoleculeT::const_iterator it=mol.begin(); it!=mol.end(); it++, atomIndex++) {
    int centerIndex = getIndexForPoint(it->position());
    if (!isValidIndex(centerIndex)) {
      std::cerr << "Error: Point out of grid" << std::endl;
      exit(1);
    }
    float atomRadius = radii[atomIndex];
    float radius = atomRadius*2; //may be +1 is enough
    int intRadius = getIntGridRadius(radius);
    int radius2 = intRadius*intRadius;

    int i_bound,j_bound,k_bound;
    i_bound = intRadius;
    for (int i = -i_bound; i <= i_bound; i++ ) {
      j_bound = (int)sqrt(radius2 - i*i);
      for (int j = -j_bound; j <= j_bound; j++) {
	k_bound = (int)sqrt(radius2 - i*i - j*j);
	for (int k = -k_bound; k <= k_bound; k++) {
	  int index = centerIndex + i + xGridNum * j + xyGridNum * k;
	  if(isValidIndex(index)) {
	    Vector3 point = getPointForIndex(index);
	    float dist = point.dist(it->position());
	    if(dist == 0.0) {
	      atoms_[index] = atomIndex;
	      continue;
	    }
	    float weight = radius/dist;
	    if(weight <= weights[index]) continue;
	    weights[index] = weight;
	    atoms_[index] = atomIndex;
	  }
	}
      }
    }
  }
  weights.clear();
}

template<class MoleculeT>
void AtomGrid::getVolume(const MoleculeT &mol, std::vector<float>& volume) {
  std::vector<int> fullCube(mol.size(), 0), halfCube(mol.size(), 0);
  for(int i=0; i< (int)grid.size();i++) {
    if(grid[i] < 0)
      fullCube[atoms_[i]]++;
    if(grid[i] == 0)
      halfCube[atoms_[i]]++;
  }
  float cubeVolume = delta*delta*delta;
  float halfVolume = cubeVolume/2;
  volume.resize(mol.size(), 0.0);
  for(unsigned int i=0; i<mol.size(); i++) {
    volume[i] = cubeVolume * fullCube[i] + halfVolume * halfCube[i];
  }
}

template<class MoleculeT>
void AtomTypeGrid::markAtomTypes(const MoleculeT &M) {
  std::vector<float> weights;
  weights.insert(weights.end(), maxEntry, 0.0);
  int atomIndex=0;
  for(typename MoleculeT::const_iterator it=M.begin(); it!=M.end(); it++, atomIndex++) {
    float atomRadius = it->getRadius();
    float radius = atomRadius*2; //may be +1 is enough
    int centerIndex = getIndexForPoint(it->position());
    if (!isValidIndex(centerIndex)) {
      std::cerr << "Error: Point out of grid" << std::endl;
      exit(1);
    }

    int intRadius = getIntGridRadius(radius);
    int radius2 = intRadius*intRadius;

    int i_bound,j_bound,k_bound;
    i_bound = intRadius;
    for (int i = -i_bound; i <= i_bound; i++ ) {
      j_bound = (int)sqrt(radius2 - i*i);
      for (int j = -j_bound; j <= j_bound; j++) {
	k_bound = (int)sqrt(radius2 - i*i - j*j);
	for (int k = -k_bound; k <= k_bound; k++) {
	  int index = centerIndex + i + xGridNum * j + xyGridNum * k;
	  if(isValidIndex(index)) {
	    Vector3 point = getPointForIndex(index);
	    float dist = point.dist(it->position());
	    if(dist == 0.0) {
	      atomTypes[index] = it->getChemType();
	      continue;
	    }
	    float weight = atomRadius/dist;
	    if(weight <= weights[index])
	      continue;
	    weights[index] = weight;
	    atomTypes[index] = it->getChemType();
	  }
	}
      }
    }
  }
  weights.clear();
}

#endif
