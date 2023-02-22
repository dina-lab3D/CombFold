#ifndef _Molecule_h
#define _Molecule_h

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "PDB.h"
#include "CIF.h"
#include "RigidTrans3.h"
// #include "Line.h"

template<class ParticleT>
/*
CLASS
  Molecule

  Template defines a collection of Particles that make up a Molecule.
  Among others class hosts a PDB file readr and linear algebra operators.

KEYWORDS
  particle, molecule, atom, vector, linear, rigid, matrix, centroid, PDB,
  file, rotation, translation

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il) and Zipi Fligelman (zipo@math.tau.ac.il)
  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997-9.

GOALS
  Serves as a container to a set of similarly typed class Particle descendants.
  The Molecule container hosts all methods of the STL vector container. It
  allows linear transfomrations to be applied to all Particles in the molecule.
  A PDB reader, enabling selected atoms to be read off a PDB format file is
  also supplied.

USAGE
  Construction is done using the default constructor. An empty Molecule object
  is returned.

  The Molecule may be built using the add method or using any other method
  supplied by the STL vector. Alternatively Particles may be added from a PDB
  format file using the readPDBfile method. This method may recieve a
  function object which will determine which of the atoms are loaded into
  memory.

  Matrices, rigid transformations, vector additions may be applied to all
  particles simultaneously. See the operators.

  The following is an example of using a Molecule<Particle> object to read in
  all C-alpha atom from a PDB file and print their coordinates after a
  rigid-transformation is applied to all the particles.
  EXAMPLE
  class CAlphaSelector : public PDB::Selector
  {  // This class defines a selector that will pick only C-alpha atoms.
  public:
  bool operator()(const char *PDBrec) const {
  const char *type = PDB::atomType(PDBrec);
  return (*(type+1) == 'C' && *(type+2) == 'A');
  }
  };

  main() {
  Molecule<Particle> M;
  ifstream pdb("pdb1tpo.ent");
  M.readPDBfile(pdb, CAlphaSelector());
  pdb.close();
  RigidTrans3 rt(Whatever);
  M*= rt;
  for (Molecule<Particle>::iterator i= M.begin(); i != M.end(); i++)
  cout << *i << '\n';

  }
  END
  Note that Molecule has 2 base classes. The MoleculeBase class is a virtual
  class that has a sole virtual method - the () operator. This class may be
  used to pass Molecule<T> objects of unknown type T into functions that
  access Molecule<T> vector positions only. This tool is used in the Match
  class.

EXAMPLE
  main() {
  Molecule<Atom> mol;
  Molecule<Atom> antibody;

  ifstream pdb("complex.pdb");
  mol.readPDBfile(pdb);
  pdb.close();

  Molecule<Atom>::ChainSelector selectorAntibody("LH");
  Molecule<Atom>::ChainSelector selectorAntigen("A");

  mol.select(selectorAntibody,antibody);
  mol.select(selectorAntigen);
  }
  END

  Note:
  If the complex contained chains L,H,A than
  after the first select mol contains L,H,A and antibody contains L,H
  after the second select mol contains A.
*/
class Molecule :
  public std::vector<ParticleT> {
public:
  // GROUP: Constructors

  //// Default constructor - empty molecule.
  Molecule() : std::vector<ParticleT>() {}

  virtual ~Molecule() {}

  // GROUP: Adding Particles.

  // Reads all ATOM records from the given PDB file according to the
  // given atom selector and returns the number of read atoms (a
  // negative value is returned in case of an error). <br> The default
  // atom selector is PDB::Selector that selects all ATOM records. To
  // write a different selection function one has to define a
  // descendant of PDB::Selector() and redefine the 'operator()(const
  // char* PDBrec) const' function to operate as desired.
  int readPDBfile(const std::string& pdbFile, const PDB::Selector& selector = PDB::Selector());

  ////
  // Reads all ATOM records from the given PDB file according to the
  // given atom selector and returns the number of read atoms (a
  // negative value is returned in case of an error). <br> The default
  // atom selector is PDB::Selector that selects all records. To write
  // a different selection function one has to define a descendant of
  // PDB::Selector() and redefine the 'operator()(const char* PDBrec)
  // const' function to operate as desired.
  int readPDBfile(std::istream& PDBstream, const PDB::Selector& selector = PDB::Selector());

  // read single MODEL
  int readModelFromPDBfile(std::istream &PDBstream, const PDB::Selector& selector = PDB::Selector());

  //// readHetAtomsFromPDBfile accepts an input stream (file or cin) a
  // PDB file format. Only HETATM records are read.
  int readHetAtomsFromPDBfile(std::istream& PDBstream, const PDB::Selector& selector = PDB::Selector());

  ////
  // Reads all ATOM and HETATM records from the given PDB file
  // according to the given atom selector and returns the number of
  // read atoms (a negative value is returned in case of an
  // error). <br> The default atom selector is PDB::Selector that
  // selects all records. To write a different selection function one
  // has to define a descendant of PDB::Selector() and redefine the
  // 'operator()(const char* PDBrec) const' function to operate as
  // desired. <br> Note that to filter out hydrogen, water molecules
  // or anything else, use the PDB::WaterHydrogenUnSelector
  int readAllPDBfile(const std::string& pdbFile, const PDB::Selector& selector = PDB::Selector());

  ////
  // Reads all ATOM and HETATM records from the given PDB file
  // according to the given atom selector and returns the number of
  // read atoms (a negative value is returned in case of an
  // error). <br> The default atom selector is PDB::Selector that
  // selects all records. To write a different selection function one
  // has to define a descendant of PDB::Selector() and redefine the
  // 'operator()(const char* PDBrec) const' function to operate as
  // desired. <br> Note that to filter out hydrogen, water molecules
  // or anything else, use the PDB::WaterHydrogenUnSelector
  int readAllPDBfile(std::istream& PDBstream, const PDB::Selector& selector = PDB::Selector());

  //// Enables to create selectors that can operate on a Molecule
  //after the Molecule object was created. For example a pdb file containing
  //several chains can be read only once to create a Molecule.
  //Each chain would be selected into a new Molecule object.
  //The selector can be operated on a Molecule using the select function.
  class Selector {
  public:
    virtual bool operator()(const ParticleT*) const {  return true;}
    virtual ~Selector() {}
  };

  ////Selects particles from a Molecule according to their chainId.
  class ChainSelector : public Selector {
  public:
    ChainSelector(const std::string &chainID): chains(chainID) {}
    bool operator()(const ParticleT* particle) const {
      for (int i=0; i< (int)chains.length(); i++)
	if (particle->chainId() == chains[i])
	  return true;
      return false;
    }
  private:
    std::string chains;
  };

  ////Selects particles from a Molecule, which are not water or hydrogens atoms.
  class WaterHydrogenUnSelector : public Selector {
  public:
    bool operator()(const ParticleT* particle) const {
      if (!(particle->isWater() || particle->isH()))
	return true;
      return false;
    }
  };

  ////Selects CA atoms Molecule
  class CaSelector : public Selector {
  public:
    bool operator()(const ParticleT* particle) const {
      if (particle->isCA())
	return true;
      return false;
    }
  };

  ////Selects backbone atoms Molecule
  class BackboneSelector : public Selector {
  public:
    bool operator()(const ParticleT* particle) const {
      if (particle->isBackbone())
	return true;
      return false;
    }
  };

  //// Excludes particles that
  //are not selected by the Selector. The operating Molecule object may
  //be changed as a result of this function.
  int select (const Selector& selector);

  ////Adds particles that
  //are selected by the Selector to the input Molecule.
  //The operating Molecule object is NOT changed as a result of this function.
  int select (const Selector& selector, Molecule<ParticleT>& m) const;

  //// Shorthand for vector push_back. Adds a new Particle at the end of the
  // Particle array.
  void add(const ParticleT& p);

  ////
  // Return a triple of the positions of the most distant atoms
  std::vector<Vector3> getDiameterTriple() const;

  // GROUP: Inspection.

  //// The () operator returns a Particle's position according to index.
  Vector3 operator()(unsigned int particle) const;

  //// Returns centroid of all particle positions.
  Vector3 centroid() const;

  // GROUP: Operators.

  //// Transforms molecule using a linear transformation of type Matrix3.
  void transform(const Matrix3 &lt);

  //// Moves all atom coordinates by given vector3 move.
  inline void translate(const Vector3 &move);

  //// Applies a rigid transformation on all molecule's Particles.
  void rigidTrans(const RigidTrans3 &rt);

  //// Calculates the principal components of the molecules particles.
  // Returns the axis system along which the covariances are 0. Also
  // returns the variance along each axis.
  //void getPrincipalComponents(Matrix3& axes, Vector3& variance);

  //// Concating two molecules together molecule m is added at the
  // end of *this molecule
  void concat(const Molecule& m);

  //// Splitting a molecule on index i all the particles
  // from 0-i (included) will stay in this molecule
  // the other particles : i+1-end-of-molecule are returned
  // as a different molecule
  Molecule splitOnIndex(const int i);

  //// Moves all atom coordinates by given vector3 move.
  Molecule& operator+=(const Vector3 &v);

  //// Applies a linear transformation on all protein coordinates.
  Molecule& operator*=(const Matrix3 &lt);

  //// Applies a rigid transformation on all atom coordinates.
  Molecule& operator*=(const RigidTrans3 &rt);

  //// Equality operator.
  Molecule& operator=(const Molecule& m);

  //// Output operator assuming the ParticleT class has its own
  // overload of the output operator.
  template<class T>
  friend std::ostream& operator<< (std::ostream& s, const Molecule& m);

  void output2cif(std::ostream& s, int modelNum=1) const;

private:
  typedef std::vector<ParticleT> BaseC;
  typedef typename BaseC::iterator base_iterator;
  typedef typename BaseC::const_iterator base_const_iterator;
};


template<class ParticleT>
int Molecule<ParticleT>::readPDBfile(std::istream &PDBstream,
                                     const PDB::Selector& selector)
{
  bool cif = CIF::readAtomSiteTable(PDBstream); // check if CIF
  PDBstream.clear();
  PDBstream.seekg(0);

  while (!PDBstream.eof()) {
    std::string line;
    getline(PDBstream, line);
    // check that line is an ATOM rec and that selector accepts line.
    // if this is the case construct a new Particle using line and add the
    // Particle to the Molecule.
    if(cif) {
      if(CIF::isATOMrec(line)) {
        std::stringstream os;
        os << ParticleT(line, true);
        if(selector(os.str().c_str())) {
          add(ParticleT(line, true));
        }
      }
    } else  {
      if (PDB::isATOMrec(line.c_str()) && selector(line.c_str())) {
        add(ParticleT(line));
      }
    }
  }

  return std::vector<ParticleT>::size();
}

template<class ParticleT>
int Molecule<ParticleT>::readModelFromPDBfile(std::istream &PDBstream,
                                              const PDB::Selector& selector)
{
  while (!PDBstream.eof()) {
    std::string line;
    getline(PDBstream, line);
    // check that line is an ATOM rec and that selector accepts line.
    // if this is the case construct a new Particle using line and add the
    // Particle to the Molecule.
    if (PDB::isATOMrec(line.c_str()) && selector(line.c_str())) {
      add(ParticleT(line.c_str()));
    }
    // finish if ENDMDL
    if (PDB::isENDMDLrec(line)) break;
  }

  return std::vector<ParticleT>::size();
}

template<class ParticleT>
int Molecule<ParticleT>::readPDBfile(const std::string& pdbFile,
                                     const PDB::Selector& selector) {
  std::ifstream pdbStream(pdbFile.c_str());

  if (!pdbStream) {
    return -1;
  }

  unsigned int numOfAtoms = readPDBfile(pdbStream, selector);
  pdbStream.close();
  return numOfAtoms;
}


template<class ParticleT>
int Molecule<ParticleT>::readHetAtomsFromPDBfile(std::istream& PDBstream,
                                                 const PDB::Selector& selector)
{
  bool cif = CIF::readAtomSiteTable(PDBstream); // check if CIF
  PDBstream.clear();
  PDBstream.seekg(0);

  while (!PDBstream.eof()) {
    std::string line;
    getline(PDBstream, line);
    // check that line is an ATOM rec and that selector accepts line.
    // if this is the case construct a new Particle using line and add the
    // Particle to the Molecule.
    if(cif) {
      if(CIF::isHETATMrec(line)) {
        std::stringstream os;
        os << ParticleT(line, true);
        if(selector(os.str().c_str())) {
          add(ParticleT(line, true));
        }
      }
    } else  {
      if (PDB::isHETATMrec(line.c_str()) && selector(line.c_str())) {
        add(ParticleT(line));
      }
    }
  }

  return std::vector<ParticleT>::size();
}

template<class ParticleT>
int Molecule<ParticleT>::readAllPDBfile(std::istream& PDBstream,
                                        const PDB::Selector& selector)
{
  bool cif = CIF::readAtomSiteTable(PDBstream); // check if CIF
  PDBstream.clear();
  PDBstream.seekg(0);

  std::string line;
  while (!PDBstream.eof()) {
    getline(PDBstream, line);
    // check that line is an HETATM or ATOM rec and that selector accepts line.
    // if this is the case construct a new Particle using line and add the
    // Particle to the Molecule.
    if(cif) {
      if((CIF::isATOMrec(line) || CIF::isHETATMrec(line))) {
        std::stringstream os;
        os << ParticleT(line, true);
        if(selector(os.str().c_str())) {
          add(ParticleT(line, true));
        }
      }
    } else  {
      if((PDB::isATOMrec(line) || PDB::isHETATMrec(line)) && selector(line.c_str())) {
        add(ParticleT(line));
      }
    }
  }
  return std::vector<ParticleT>::size();
}


template<class ParticleT>
int Molecule<ParticleT>::readAllPDBfile(const std::string& pdbFile,
                                        const PDB::Selector& selector) {
  std::ifstream pdbStream(pdbFile.c_str());

  if (!pdbStream) {
    return -1;
  }

  unsigned int numOfAtoms = readAllPDBfile(pdbStream, selector);
  pdbStream.close();
  return numOfAtoms;
}

template<class ParticleT>
int Molecule<ParticleT>::select(const Selector& selector){
  Molecule<ParticleT> m;
  for (base_const_iterator it= std::vector<ParticleT>::begin(); it != std::vector<ParticleT>::end(); ++it) {
    if (selector (&(*it)))
      m.add( (*it) );
  }
  *this=m;
  return std::vector<ParticleT>::size();
}

template<class ParticleT>
int Molecule<ParticleT>::select(const Selector& selector, Molecule<ParticleT>& m)const{
  for (base_const_iterator it=std::vector<ParticleT>::begin(); it!=std::vector<ParticleT>::end(); ++it) {
    if (selector (&(*it)))
      m.add( (*it) );
  }
  return m.size();
}

template<class ParticleT>
void Molecule<ParticleT>::add(const ParticleT& p)
{
  this->push_back(p);
}
/*
template<class ParticleT>
std::vector<Vector3> Molecule<ParticleT>::getDiameterTriple() const{
  std::vector<Vector3> triple(3);

  float maxDist2 = 0;
  for (unsigned int i = 0 ; i < this->size() ; i++) {
    for (unsigned int j = 0 ; j < this->size() ;j++) {
      float dist2 = (*this)(i).dist2((*this)(j));

      if (dist2 > maxDist2) {
	maxDist2 = dist2;
	triple[0] = (*this)(i);
	triple[1] = (*this)(j);
      }
    }
  }

  float maxDist = 0;
  for (unsigned int i = 0 ; i < this->size() ; i++) {
    Line line(triple[0], triple[1]);
    float dist = line.distance(operator()(i));

    if (dist > maxDist) {
      maxDist = dist;
      triple[2] = (*this)(i);
    }
  }

  return triple;
}
*/
template<class ParticleT>
Vector3 Molecule<ParticleT>::operator()(unsigned int particle) const
{
  return(*(this->begin()+particle)).position();
}

template<class ParticleT>
Vector3 Molecule<ParticleT>::centroid() const
{
  Vector3 sum;
  for (base_const_iterator i = this->begin(); i != this->end(); ++i)
    sum+= *i;
  sum /= this->size();
  return sum;
}

template<class ParticleT>
void Molecule<ParticleT>::transform(const Matrix3& lt)
{
  for (base_iterator i=this->begin(); i!=this->end(); ++i)   (*i)*=lt;
}

template<class ParticleT>
void Molecule<ParticleT>::translate(const Vector3 &move)
{
  for (base_iterator i=this->begin(); i!=this->end(); ++i)   (*i)+=move;
}

template<class ParticleT>
void Molecule<ParticleT>::rigidTrans(const RigidTrans3& rt)
{
  for (base_iterator i=this->begin(); i!=this->end(); ++i)   (*i)*=rt;
}

// template<class ParticleT>
// void Molecule<ParticleT>::getPrincipalComponents(Matrix3& axes,
//                                                  Vector3& variance)
// {
//   Vector3 mean = centroid();
//   Matrix3 cov;
//   for(base_const_iterator i= begin(); i!= end(); ++i) {
//     Vector3 v = (*i).position() - mean;
//     cov+= Matrix3(v, v);
//   }
//   cov/= size();
//   Matrix3 dummy;
//   cov.svd(dummy, variance, axes);
// }

template<class ParticleT>
void Molecule<ParticleT>::concat(const Molecule<ParticleT>& m)
{
  for (base_const_iterator i = m.begin(); i!= m.end(); ++i)
    this->push_back(*i);
}

template<class ParticleT>
Molecule<ParticleT> Molecule<ParticleT>::splitOnIndex(const int i)
{
  if (i > this->size())
    return Molecule();
  Molecule mol;
  for (base_iterator it = this->begin()+i; it != this->end(); ++it)
    mol.push_back(*it);
  erase(this->begin()+i, this->end());
  return mol;
}

template<class ParticleT>
Molecule<ParticleT>& Molecule<ParticleT>::operator+=(const Vector3 &v)
{
  translate(v);
  return *this;
}

template<class ParticleT>
Molecule<ParticleT>& Molecule<ParticleT>::operator*=(const Matrix3 &lt)
{
  transform(lt);
  return *this;
}

template<class ParticleT>
Molecule<ParticleT>& Molecule<ParticleT>::operator*=(const RigidTrans3 &rt)
{
  rigidTrans(rt);
  return *this;
}

template<class ParticleT>
Molecule<ParticleT>& Molecule<ParticleT>::operator=(const Molecule<ParticleT>& m) {
  return(Molecule<ParticleT>&)std::vector<ParticleT>::operator=(m);
}

template<class ParticleT>
std::ostream& operator<< (std::ostream& s, const Molecule<ParticleT>& m)
{
  for (typename std::vector<ParticleT>::const_iterator it = m.begin();
       it != m.end() ; ++it) {
    s << *it << std::endl;
  }
  return s;
}

template<class ParticleT>
void Molecule<ParticleT>::output2cif(std::ostream& s, int modelNum) const {
  s << "data_1XXX" << std::endl;
  CIF::writeAtomSiteTable(s);
  for (typename std::vector<ParticleT>::const_iterator it = this->begin(); it != this->end() ; ++it) {
    it->output2cif(s, modelNum);
    s << std::endl;
  }
  s << "#" << std::endl;
}

#endif
