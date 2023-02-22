#include <stdio.h>
#include "MoleculeGrid.h"
#include "Atom.h"

#define MAXIMIZE(a,b) if (a < b) a = b;
#define MINIMIZE(a,b) if (a > b) a = b;

MoleculeGrid::MoleculeGrid(const Surface& surface,
			   const float inDelta, const float maxRadius): delta(inDelta)
{
  maxGridRadius = getIntGridRadius(maxRadius) + 1;

  computeBoundingValues(surface);
  xMin = xMin - maxRadius - delta;
  yMin = yMin - maxRadius - delta;
  zMin = zMin - maxRadius - delta;
  xMax = xMax + maxRadius + delta;
  yMax = yMax + maxRadius + delta;
  zMax = zMax + maxRadius + delta;

  xGridNum = (int) ((xMax - xMin) / delta + 2);
  yGridNum = (int) ((yMax - yMin) / delta + 2);
  zGridNum = (int) ((zMax - zMin) / delta + 2);
  xyGridNum = xGridNum * yGridNum ;

  if ( xGridNum < 1 || yGridNum < 1 || zGridNum < 1 )
  {
    std::cerr << "Invalid parameters for grid generation" << std::endl;
    exit(1);
  }
  maxEntry = xGridNum * yGridNum * zGridNum;
  //cerr << " maxEntry " << maxEntry << endl;

  grid.insert(grid.end(), maxEntry, (float)MAX_FLOAT) ;
  int nc = 0;
  for (int x = -1; x <= 1; x++) {
    for (int y = -1; y <= 1; y++) {
      for (int z = -1; z <= 1; z++) {
        if (x == 0 && y == 0 && z == 0) continue;
        neigbor.push_back(x + xGridNum*y + xyGridNum*z);
        float d = delta * sqrt((float)(x*x + y*y + z*z));
        neigborDist.push_back(d);
        nc++;
      }
    }
  }
}

void MoleculeGrid::computeBoundingValues(const Surface& surface)
{
  xMin=yMin=zMin=(float)MAX_FLOAT;
  xMax=yMax=zMax=(float)MIN_FLOAT;
  for(Surface::const_iterator it = surface.begin(); it!= surface.end(); ++it){
    const Vector3 &v = it->position();
    MINIMIZE(xMin, v[0]);
    MINIMIZE(yMin, v[1]);
    MINIMIZE(zMin, v[2]);
    MAXIMIZE(xMax, v[0]);
    MAXIMIZE(yMax, v[1]);
    MAXIMIZE(zMax, v[2]);
  }
}

void MoleculeGrid::computeDistFromSurface(const Surface &surface, bool precise)
{
  std::vector<Vector3> sps;
  for (Surface::const_iterator it = surface.begin(); it!= surface.end(); ++it) {
    sps.push_back(it->position());
  }

  std::vector<long> layer1, layer2, *curr = &layer1, *next = &layer2, *tmp;
  std::vector<Vector3*> sps1, sps2, *currSps = &sps1, *nextSps = &sps2, *tmpSps;
  std::vector<int> gridLayer;
  gridLayer.insert(gridLayer.end(),maxEntry,-1);

  layer1.reserve(sps.size()*10);
  layer2.reserve(sps.size()*10);
  sps1.reserve(sps.size()*10);
  sps2.reserve(sps.size()*10);

  // mark zero in Grid for each surface point and insert indexes of voxels for first layer
  int index;
  for (std::vector<Vector3>::iterator it = sps.begin(); it!= sps.end(); ++it) {
    index = getIndexForPoint(*it);
    if (gridLayer[index] != 0) {
      if (precise) sps1.push_back( &(*it));
      layer1.push_back(index);
      grid[index] = 0.0;
      gridLayer[index] = 0;
    }
  }

  // work on every layer to be computed
  Vector3* p = NULL;
  for (int layer = 0; layer < maxGridRadius; layer++) {
    std::cerr << "    In Layer " << layer << " out of " << maxGridRadius
              << " curr->size() " << curr->size() << std::endl;
    // update voxels with current layer distance and insert indexes for next layer
    for (unsigned int i = 0; i < curr->size(); i++) {
      int index = (*curr)[i];
      for (int j = 0; j < 26; j++) {
        int nIndex = index + neigbor[j];
        if (isValidIndex(nIndex)) {
          if (precise) {
            p = (*currSps)[i];
            if (setDIST(nIndex, *p)) {
              if (gridLayer[nIndex] < layer + 1) {
                next->push_back(nIndex);
                nextSps->push_back(p);
                gridLayer[nIndex] = layer + 1;
              }
            }
          } else {
            float dist = grid[index] + neigborDist[j];
            if (dist <= 0) {
              std::cerr << " dist = grid[index] + neigborDist[j] grid[nIndex] "
                        << dist << " " << grid[index] << " " << neigborDist[j]
                        << " j " << j << " nIndex " << nIndex << " grid[nIndex] " << grid[nIndex] << std::endl;
            }
            if (grid[nIndex] > dist) {
              grid[nIndex] = dist;
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


float MoleculeGrid::calcVolumeFunc(const Vector3 &center, const float radius) const
{
  int centerIndex = getIndexForPoint(center);
  if (!isValidIndex(centerIndex)) {
    std::cerr << "Error: Point out of grid" << std::endl;
    exit(1);
  }

  int intRadius = getIntGridRadius(radius);
  int radius2 = intRadius*intRadius;
  //maybe compute exact sphere volume??

  int insideNum = 0;
  int outsideNum = 0;
  int i_bound,j_bound,k_bound;
  i_bound = intRadius;
  for (int i = -i_bound; i <= i_bound; i++ )
    {
      j_bound = (int)sqrt((float)radius2 - i*i);
      for (int j = -j_bound; j <= j_bound; j++)
	{
	  k_bound = (int)sqrt((float)radius2 - i*i - j*j);
	  for (int k = -k_bound; k <= k_bound; k++)
	    {
	      int index = centerIndex + i + xGridNum * j + xyGridNum * k;
	      if(isValidIndex(index))
		{
		  if (grid[index] > 0)
		    outsideNum++;
		  else
		    insideNum++;
		}
	      else
		outsideNum++;
	    }
	}
    }

  if(insideNum==0 && outsideNum == 0) {
    std::cerr << "Error in calculateVolumeFunc" << std::endl;
  }
  return ((float)insideNum) / ( insideNum + outsideNum );
}

Vector3 MoleculeGrid::calcVolumeNormal(const Vector3 &center, const float radius) const
{
  int centerIndex = getIndexForPoint(center);
  if (!isValidIndex(centerIndex)) {
    std::cerr << "Error: Point out of grid" << std::endl;
    exit(1);
  }

  int intRadius = getIntGridRadius(radius);
  int radius2 = intRadius*intRadius;

  int xOut = 0;
  int yOut = 0;
  int zOut = 0;

  int i_bound,j_bound,k_bound;
  i_bound = intRadius;
  for (int i = -i_bound; i <= i_bound; i++ )
    {
      j_bound = (int)sqrt((float)radius2 - i*i);
      for (int j = -j_bound; j <= j_bound; j++)
	{
	  k_bound = (int)sqrt((float)radius2 - i*i - j*j);
	  for (int k = -k_bound; k <= k_bound; k++)
	    {
	      int index = centerIndex + i + xGridNum * j + xyGridNum * k;
	      if(isValidIndex(index))
		{
		  if(grid[index] <= 0)
		    continue;
		}
	      xOut+=i;
	      yOut+=j;
	      zOut+=k;
	    }
	}
    }
  Vector3 v((float)xOut, (float)yOut, (float)zOut);

  if(v.norm() == 0) {
    std::cerr << "Error in calculateVolumeNormFunc center" << center << " dist " << grid[centerIndex] << std::endl;
  }
  return v/v.norm();
}

void MoleculeGrid::calcVolumeFuncAndNormal(const Vector3 &center, const float radius,
					   float &volFunc, Vector3& volNormal) const
{
  int centerIndex = getIndexForPoint(center);
  if (!isValidIndex(centerIndex)) {
    std::cerr << "Error: Point out of grid" << std::endl;
    exit(1);
  }

  int intRadius = getIntGridRadius(radius);
  int radius2 = intRadius*intRadius;

  int outsideNum = 0;
  int insideNum = 0;
  int xOut = 0;
  int yOut = 0;
  int zOut = 0;

  int i_bound,j_bound,k_bound;
  i_bound = intRadius;
  for (int i = -i_bound; i <= i_bound; i++ )
    {
      j_bound = (int)sqrt((float)radius2 - i*i);
      for (int j = -j_bound; j <= j_bound; j++)
	{
	  k_bound = (int)sqrt((float)radius2 - i*i - j*j);
	  for (int k = -k_bound; k <= k_bound; k++)
	    {
	      int index = centerIndex + i + xGridNum * j + xyGridNum * k;
	      if(isValidIndex(index) && grid[index] <= 0.0)
		{
		    insideNum++;
		}
	      else
		{
		  outsideNum++;
		  xOut+=i;
		  yOut+=j;
		  zOut+=k;
		}
	    }
	}
    }
  volFunc= ((float)insideNum) / (insideNum + outsideNum);
  volNormal.updateX((float)xOut);
  volNormal.updateY((float)yOut);
  volNormal.updateZ((float)zOut);

  if(volNormal.norm() == 0) {
    std::cerr << "Error in calculateVolumeNormFunc " << center << std::endl;
  }
  volNormal=volNormal/volNormal.norm();
}

void MoleculeGrid::printGrid(std::ofstream& file) const
{
  const float SURFACE_THRESHOLD = 0;//delta * sqrt(3.0);
  for(int i=0; i< (int)grid.size();i++)
    {
      if (grid[i] > SURFACE_THRESHOLD || grid[i] < -SURFACE_THRESHOLD) continue;
      Vector3 point = getPointForIndex(i);
      Atom atom(point, 'A', i, i % 20, "PRB", 'X');
      file << atom << std::endl;
    }
}

float MoleculeGrid::volume() const {
  int fullCube(0), halfCube(0);
  for(int i=0; i< (int)grid.size();i++) {
    if(grid[i] < 0)
      fullCube++;
    if(grid[i] == 0)
      halfCube++;
  }
  float cubeVolume = delta*delta*delta;
  float halfVolume = cubeVolume/2;
  return fullCube*cubeVolume + halfCube*halfVolume;
}
