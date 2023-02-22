#include "prGrid.h"

ProbGrid::ProbGrid(const Surface &surface, const float inDelta, const float maxRadius):
  MoleculeGrid(surface, inDelta, maxRadius)
{
  probs.insert(probs.end(), maxEntry, 0.0);
}

float ProbGrid::getProb(const Vector3& point) const
{
  int index = getIndexForPoint(point);
  if (!isValidIndex(index))
    return MAX_FLOAT;
  return probs[index];
}

ResidueGrid::ResidueGrid(const Surface &surface, const float inDelta, const float maxRadius):
  MoleculeGrid(surface, inDelta, maxRadius)
{
  residues.insert(residues.end(), maxEntry, (int)0);
}

int ResidueGrid::getResidueEntry(const Vector3& point) const {
  int index = getIndexForPoint(point);
  if(!isValidIndex(index))
    return 0;
  return residues[index];
}

void ResidueGrid::printGrid(std::ofstream& file) const {
  for(int i=0; i< (int)residues.size();i++) {
    if(residues[i] < 0) {
      Vector3 point = getPointForIndex(i);
      Atom atom(point, 'A', i, residues[i], "PRB", 'X');
      file << atom << std::endl;
    }
  }
}

// const std::vector<unsigned int>& AtomGrid::getAtoms(const Vector3 &point) const {
//   const static std::vector<unsigned int> emptyVec;
//   int index = getIndexForPoint(point);
//   if(isValidIndex(index))
//     return atoms[index];
//   return emptyVec;
// }

// unsigned int AtomGrid::getAtomsNumber(const Vector3 &point) const {
//   int index = getIndexForPoint(point);
//   if(!isValidIndex(index))
//     return 0;
//   return atoms[index].size();
// }

AtomTypeGrid::AtomTypeGrid(const Surface &surface, const float inDelta, const float maxRadius):
  MoleculeGrid(surface, inDelta, maxRadius), atomTypes(maxEntry) {}

SurfacePointGrid::SurfacePointGrid(const Surface &surface, const float inDelta, const float maxRadius):
  MoleculeGrid(surface, inDelta, maxRadius)
{
  spPointers.insert(spPointers.end(), maxEntry, (SurfacePoint*)NULL);
}

void SurfacePointGrid::computeDistAndPointers(const Surface &surface, bool precise) {
  std::vector<long> layer1, layer2, *curr = &layer1, *next = &layer2, *tmp;
  std::vector<const SurfacePoint*> sps1, sps2, *currSps = &sps1, *nextSps = &sps2, *tmpSps;
  std::vector<int> gridLayer;
  gridLayer.insert(gridLayer.end(),maxEntry,-1);

  layer1.reserve(surface.size()*10);
  layer2.reserve(surface.size()*10);
  sps1.reserve(surface.size()*10);
  sps2.reserve(surface.size()*10);

  // mark zero in Grid for each surface point and insert indexes of voxels for first layer
  int index;
  for (Surface::const_iterator it = surface.begin(); it!= surface.end(); ++it) {
    index = getIndexForPoint(it->position());
    if (gridLayer[index] != 0) {
      if (precise) sps1.push_back(&(*it));
      layer1.push_back(index);
      grid[index] = 0.0;
      spPointers[index] = &(*it);
      gridLayer[index] = 0;
    }
  }

  // work on every layer to be computed
  const SurfacePoint* p = NULL;
  for (int layer = 0; layer < maxGridRadius; layer++) {
    std::cerr << "    In Layer " << layer << " out of " << maxGridRadius << " curr->size() " << curr->size() << std::endl;
    // update voxels with current layer distance and insert indices for next layer
    for (unsigned int i = 0; i < curr->size(); i++) {
      int index = (*curr)[i];
      for (int j = 0; j < 26; j++) {
        int nIndex = index + neigbor[j];
        if (isValidIndex(nIndex)) {
          if (precise) {
            p = (*currSps)[i];
            if (setDistAndPointer(nIndex, p)) {
              if (gridLayer[nIndex] < layer + 1) {
                next->push_back(nIndex);
                nextSps->push_back(p);
                gridLayer[nIndex] = layer + 1;
              }
            }
          } else {
            float dist = grid[index] + neigborDist[j];
            if ( dist <= 0 ) {
              std::cerr << " dist = grid[index] + neigborDist[j] grid[nIndex] " << dist << " " <<  grid[index]  << " " << neigborDist[j] << " j " << j << " nIndex " << nIndex << " grid[nIndex] " << grid[nIndex] << std::endl;
            }
            if (grid[nIndex] > dist) {
              grid[nIndex] = dist;
	      spPointers[nIndex] = spPointers[index];
              if (gridLayer[nIndex] < layer + 1) {
                next->push_back(nIndex);
                gridLayer[nIndex] = layer + 1;
              }
            }
          }
        }
      }
    }
    curr->clear();
    tmp = curr; curr = next; next = tmp;
    if (precise) {
      currSps->clear();
      tmpSps = currSps; currSps = nextSps; nextSps = tmpSps;
    }
  }
}

VolumeGrid::VolumeGrid(const Surface &surface, const float inDelta, const float maxRadius):
  MoleculeGrid(surface, inDelta, maxRadius)
{
  volumeNormals.insert(volumeNormals.end(), maxEntry, Vector3());
  volumeFunctions.insert(volumeFunctions.end(), maxEntry, 0);
}

void VolumeGrid::getVolumeFuncAndNormal(const Vector3 &center, float &volFunc, Vector3& volNormal) const {
  int centerIndex = getIndexForPoint(center);
  if (!isValidIndex(centerIndex)) {
    std::cerr << "Error: Point out of grid" << std::endl;
    exit(1);
  }
  volFunc = volumeFunctions[centerIndex];
  volNormal = volumeNormals[centerIndex];
}

float VolumeGrid::getVolumeFunc(const Vector3 &center) const {
  int centerIndex = getIndexForPoint(center);
  if (!isValidIndex(centerIndex)) {
    std::cerr << "Error: Point out of grid" << std::endl;
    exit(1);
  }
  return volumeFunctions[centerIndex];
}

const Vector3& VolumeGrid::getVolumeNormal(const Vector3 &center) const {
  int centerIndex = getIndexForPoint(center);
  if (!isValidIndex(centerIndex)) {
    std::cerr << "Error: Point out of grid" << std::endl;
    exit(1);
  }
  return volumeNormals[centerIndex];
}

void VolumeGrid::computeVolumes(const float radius, const float maxRadius) {
  for(unsigned int gridIndex=0; gridIndex < grid.size(); gridIndex++) {
    if(fabs(grid[gridIndex]) <= maxRadius) {
      float r = (radius > 1+fabs(grid[gridIndex])? radius : 1+fabs(grid[gridIndex]));
      int intRadius = getIntGridRadius(r);
      int radius2 = intRadius*intRadius;

      int outsideNum = 0;
      int insideNum = 0;
      int xOut(0), yOut(0), zOut(0);

      int i_bound,j_bound,k_bound;
      i_bound = intRadius;
      for (int i = -i_bound; i <= i_bound; i++ ) {
	j_bound = (int)sqrt(radius2 - i*i);
	for (int j = -j_bound; j <= j_bound; j++) {
	  k_bound = (int)sqrt(radius2 - i*i - j*j);
	  for (int k = -k_bound; k <= k_bound; k++)  {
	    int index = gridIndex + i + xGridNum * j + xyGridNum * k;
	    if(isValidIndex(index) && grid[index] <= 0.0) {
	      insideNum++;
	    } else {
	      outsideNum++;
	      xOut+=i;
	      yOut+=j;
	      zOut+=k;
	    }
	  }
	}
      }
      volumeFunctions[gridIndex] = ((float)insideNum) / (insideNum + outsideNum);
      volumeNormals[gridIndex] = Vector3((float)xOut, (float)yOut, (float)zOut);
      if(volumeNormals[gridIndex].norm() == 0) {
	std::cerr << "Error in calculateVolumeNormFunc " << std::endl;
      }
      volumeNormals[gridIndex] = volumeNormals[gridIndex]/volumeNormals[gridIndex].norm();
    }
  }
}
