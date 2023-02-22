#ifndef _Particle_h

#define _Particle_h

#include "RigidTrans3.h"
/*
CLASS
  Particle

  A class meant to serve as a basis for particle types such as atoms, molecular
  surface points etc. This class's only property is position and operators to
  alter the Particle's position are supplied. Users of the class are encouraged
  to use the class as an inheritance base in creating their Particle types.

KEYWORDS
  vector, position, algebra, linear, rigid, matrix, distance, match

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI>
</UL>

GOALS
  Class Particle is supposed to serve as a base-class for particle type classes
  such as atoms or surface points. The class may be constructed in the usual
  manner supplying the particle's position or it may be constructed using a
  PDB format ATOM record. The class hosts a set of operators used to
  shift particle positions.

USAGE
  Below is a piece of code used to read in PDB atoms, shift their positions
  using a rigid transformation and print their new positions to the standard
  output.
  EXAMPLE
    RigidTrans3 rt = ....;             // rigid transformation = whatever
    ifstream PDBfile("pdb1tpo.ent");   // open PDB file (e.g. 1TPO)
    char line[85];
    while (!PDBfile.eof()) {
      PDBfile.getline(line, 85);  // Get a PDB file record.
      if (PDB::isATOMrec(line)) {   // check that the PDB rec is an ATOM record
        Particle p(line);           // Construct a Particle using ATOM record.
        p*= rt;                     // Shift particle's position using a rigid
        cout << p.position << '\n'; // transformation and print new position.
      }
    } 
  END
*/
class Particle
{
public:
  //// Creates a Particle with 0-position.
  Particle();
  //// Creates a particle with given position NewPos.
  explicit Particle(const Vector3& newPos);

  //// Constructs a Particle using a PDB ATOM record. The constructor assumes
  // that the given string is an ATOM record. If this is not so the object
  // constructed may hold faulty data.
  explicit Particle(const char* const PDBrec);

  //// Returns particle's position.
  Vector3 position() const;

  //// Moves particle to a new position.
  void moveto(const Vector3 &);

  //// Default Vector3 cast. When vector3 cast is used the particle's position
  // is returned.
  operator Vector3() const;

  //// Default Vector3 cast. When vector3 cast is used the particle's position
  // is returned.
  operator const Vector3&() const;
  
  //// Shift particle's position using a linear transformation.
  Particle& operator*=(const Matrix3 &);

  //// Shift particle's position using a rigid transformation.
  Particle& operator*=(const RigidTrans3 &);

  //// Shift particle's position using vector addition.
  Particle& operator+=(const Vector3 &);

  //// Shift particle's position using vector subtraction.
  Particle& operator-=(const Vector3 &);

  //// Shift particle's position using scalar multiplication.
  Particle& operator*=(const float &);

  //// Shift particle's position using scalar division.
  Particle& operator/=(const float &);

protected:
  // changed to protected by Zipi Fligelman in order to make this information availabnle to descendant
  Vector3 pos;
};

#endif










