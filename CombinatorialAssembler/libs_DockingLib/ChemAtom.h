#ifndef CHEM_ATOM_H
#define CHEM_ATOM_H

#include "Atom.h"

enum HB_TYPE {DONOR = 0, ACCEPTOR = 1, HB_BOTH = 2, NONE = 3};
enum MOL2_TYPE {C2 = 0, C3 =1 , Car = 2, Ccat = 3, N2 = 4, N4 = 5, Nam = 6,
                Nar = 7, Npl3 = 8, O2 = 9, O3 = 10, Oco2 = 11, P3 = 12, S3 = 13, UNK_MOL2_TYPE = 14 };

/*
CLASS
  ChemAtom

  This is a class that addes additional atrributes to Atom

KEYWORD
  ChemAtom, radius, charge, type, ASA, solvation

AUTHORS
  Dina Schneidman (duhovka@tau.ac.il)

  copyright: GAMBA Group , Tel-Aviv Univ. Israel, 2004.

GOALS
  ChemAtom extends Atom and has 4 additional attributes: chemical
  type, radius, charge and probability of the atom of being steady.

CHANGES LOG
<UL>
</UL>

USAGE
*/

class ChemAtom : public Atom { //, public EnergyAtom {
 public:
  // GROUP: Constructors

  //// empty atom
  ChemAtom();

  //// construction with atom data
  ChemAtom(const Vector3& p, const char ch, const unsigned int aid,
	   const unsigned int resid, const char* const t, const char resTp,
	   const int chem_type=0, const float rad=1.5, const float charge=0);

  ChemAtom(const std::string& PDBrec, bool cif = false);


  // GROUP: modifiers

  //// set chemical type
  void setChemType(const int chem_type) {
    if(chem_type <= 0 || chem_type > 18) {
      std::cerr << "Error in ChemType for atom " << *this << std::endl;
    } else {
      chemType_ = chem_type;
    }
  }

  //// set radius
  void setRadius(const float r) { radius_ = r; }

  //// set epsilon
  void setEpsilon(const float e) { epsilon_ = e; }

  //// set charge
  void setCharge(const float c) { charge_ = c; }

  //// set HB data
  void setHBData(const HB_TYPE hbType, const Vector3& hbDirection);

  //// set probability of being steady
  void setSteadyProb(const float p) { steadyProb_ = p; }

  //// set ASA
  void setASA(const float ASA) { ASA_ = ASA; }

  //// set atom position
  void setPosition(const Vector3& v) { update(v); }

  //// set mol2 type
  void setMol2Type(const MOL2_TYPE t) { mol2Type_ = t; }

  // GROUP: Inspectors

  //// get atom chemical type for ACE computation
  int getChemType() const {
    if(chemType_ <= 0 || chemType_ > 18) return 0;
    return chemType_;
  }

  //// get atom radius
  float getRadius() const { return radius_; }

  //// get epsilon value
  float getEpsilon() const { return epsilon_; }

  ////
  bool isHydrogen() const { return isH(); }

  //// get atom charge
  float getCharge() const { return charge_; }

  ////
  bool isDonor() const { return (hbType_ == HB_BOTH || hbType_ == DONOR); }

  ////
  bool isAcceptor() const { return (hbType_ == HB_BOTH || hbType_ == ACCEPTOR); }

  //// returns direction of h-donor/acceptor, in case of undefined direction zero vector is returned
  const Vector3& getHBDirection() const { return hbDirection_; }

  //// get probability of being steady
  float getSteadyProb() const { return steadyProb_; }

  //// get ASA
  float getASA() const { return ASA_; }

  ////
  bool isNonPolar() const { return isNonPolar_; }

  ////
  MOL2_TYPE getMol2Type() const { return mol2Type_; }

  Vector3 position() const {return (Atom(*this)).position(); }

  friend Vector3& operator*=(ChemAtom &v,const RigidTrans3 &rt) {
    //cout << "ChemAtom *= RigidTrans3" << endl;
    v += (rt.rotation() * v) + rt.translation() - v;
    if (!v.hbDirection_.isZero())
      v.hbDirection_ = rt.rotation() * v.hbDirection_;
    return v;
  }


  //// maximal atom radius
  static const float maxRadius;


private:
  void setPolarity();

 private:
  int chemType_; // type required for ACE computation
  float radius_;
  float charge_;
  float epsilon_; // for vdw computation
  HB_TYPE hbType_;
  Vector3 hbDirection_;
  bool isNonPolar_;
  float steadyProb_;
  float ASA_;
  MOL2_TYPE mol2Type_;
};
#endif /* IMP_CHEMATOM_H */
