#ifndef _Surface_h
#define _Surface_h

#include <iostream>
#include "SurfacePoint.h"
#include "Molecule.h"

/*
CLASS
   Surface

   Defines a collection of SurfaceParticle's that make up the
   surface of a molecule. The class hosts a reader that reads output out of
   a Shou surface file.

KEYWORDS
   surface particle, molecule, atom, vector, linear, rigid, matrix
   PDB, Connolly, Shou, rotation, translation

AUTHORS
  Zipi Fligelman (zipo@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI>
12/01/04 - area method added
</UL>

GOALS
  Holds a list of surface points. Surfcae points are defined using the surface
  point and the normal to the surface at the point. The class is descendant
  of Molecule<SurfacePoint> and hosts all the methods supported by the
  Molecule template.

USAGE

*/

class Surface
  : public Molecule<SurfacePoint>
{
public:

  //// Default constructor - empty Surface.
  Surface();

  //// We read a Surface file/s and we add them to the list of surface point
  // with the appropriate constructor.
  int readShouFile(std::istream& shouStream);

  //// Same as readShouFile
  int readSurfaceFile(std::istream& inFile);

  ////
  int readSurfaceFile(const std::string fileName);

  //// computes surface area by summing the areas of all the points
  float area() const;

  void output2PDB(const std::string outFileName);

  //// Output the molecular surface as a VRML 1.0 model. This function may be
  // used to display the molecular surfac as a 3d model.
  void outputVRML(std::ostream& VRMLstream,
                  const float red = 1.0,
                  const float green = 1.0,
                  const float blue = 1.0) const;

  friend std::ostream& operator<<(std::ostream& s, const Surface& surf);

private:
  int readBinaryFile(std::istream& inFile);

  int readASCIIFile(std::istream& inFile);

  bool isBinaryFile(std::istream& inFile);
};

#endif
