#ifndef _Shou_h

#define _Shou_h
/*
CLASS
  Shou

  Defines a set of methods to read Shou format files. The Shou format files 
  are created using the PDBCP program which outputs a set of critical points 
  describing the molecular surface.

KEYWORDS
  Molecule, Surface

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI>
</UL>

GOALS
  Class Shou is designed to encapsulate a set of methods used for reading the
  Shou file format. Methods are static requiring no object instantiation
  before their use.

USAGE
  The class was written as a takeoff on the PDB class used to read the PDB
  format files. Basically the class encapsulates a set fo methods to extract
  data from a single line found in the Shou surface file. The class's methods
  are all static, they all recieve a char. string as input and they all return
  the requested data as output.
*/
class Shou
{
public:
  //// Checks that the shouRec string is actually a record. This will rule out
  // all REMARKs and # lines.
  static bool isRec(const char* const shouRec);

  static const unsigned short xCoordField  = 15;
  static const unsigned short yCoordField  = 23;
  static const unsigned short zCoordField  = 31;
  static const unsigned short surfaceField = 39;
  static const unsigned short xNormField   = 47;
  static const unsigned short yNormField   = 54;
  static const unsigned short zNormField   = 61;

  //// returns the 1 (for caps), 2 (for belts) or 3 (for pits) atoms 
  // involved in the definition of the surface point. The atom parameter
  // may assume values of 0..2. For each such value the method will return
  // a different atom index involved in the definition of the critical point.
  // For caps and belts which are defined by less then 3 atoms 0s will be
  // returned for the upper indices.
  static unsigned short atomIndex(const char *const shouRec, 
				  unsigned int atom);

  //// Returns x coordinate of the surface point.
  static float xCoord(const char* const shouRec);

  //// Returns y coordinate of the surface point.
  static float yCoord(const char* const shouRec);

  //// Returns z coordinate of the surface point.
  static float zCoord(const char* const shouRec);

  //// Returns surface area of element represented by vector.
  static float surface(const char* const shouRec);

  //// Returns x coordinate of the surface normal.
  static float xNorm(const char* const shouRec);

  //// Returns y coordinate of the surface normal.
  static float yNorm(const char* const shouRec);

  //// Returns z coordinate of the surface normal.
  static float zNorm(const char* const shouRec);

private:

  // Class should not be constructed since all methods are static. The class is
  // used to encapsulate the methods.
  Shou();
};

#endif
