#ifndef SASURFACE_H
#define SASURFACE_H

#include "DotSphere.h"
//#include "Torus.h"
#include "MolIntersection.h"
#include <Surface.h>
#include "ChemMolecule.h"
#include "Interface.h"
#include <Interface.h>


/*
  Class that computes solvent accessible surface for the molecule
*/
class SASurface {
 public:
  SASurface(ChemMolecule& mol, float probe_radius, float d);

  unsigned int buildCapsSurface(Surface& surface);
  
  //  unsigned int buildBeltsSurface(Surface& surface);

  //// computes ASA for each atom in ChemMolecule and returns ASA of the molecule
  float computeASAForAtoms();

 private:
  ChemMolecule& molecule_;
  float probeRadius_;
  float density_;
  MolIntersection<ChemMolecule> molIntersection_;
};
  
#endif
