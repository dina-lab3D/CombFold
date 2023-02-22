#include "ChemAtom.h"

const float ChemAtom::maxRadius = 3.5;
ChemAtom::ChemAtom(){}

ChemAtom::ChemAtom(const Vector3& p, const char ch, const unsigned int aid,
		   const unsigned int resid, const char* const t, const char resTp,
		   const int chemType, const float radius, const float charge):
  Atom(p, ch, aid, resid, t, resTp),
  chemType_(chemType), radius_(radius), charge_(charge), epsilon_(-0.01),
  hbType_(NONE), hbDirection_(Vector3(0,0,0)), ASA_(0) {
  setPolarity();
}

ChemAtom::ChemAtom(const std::string& PDBrec, bool cif) :
  Atom(PDBrec, cif),
  chemType_(0), radius_(1.5), charge_(0), epsilon_(-0.01), hbType_(NONE), hbDirection_(Vector3(0,0,0)), ASA_(0) {
  setPolarity();
}

void ChemAtom::setHBData(const HB_TYPE hbType, const Vector3& hbDirection) {
  // set type
  if(hbType_ == NONE) hbType_ = hbType;
  if((hbType_ == DONOR && hbType == ACCEPTOR) || (hbType_ == ACCEPTOR && hbType == DONOR)) hbType_ = HB_BOTH;

  // set direction
  if(hbDirection_.isZero()) {
    hbDirection_ = hbDirection;
  }
}

void ChemAtom::setPolarity() {
  const char res = residueType();
  const char* atomT = type();

  bool nonPolar = true;

  if (atomT[1] == 'S') nonPolar = true;

  if (atomT[1] != 'C') nonPolar = false;
  if ((atomT[2] == 'A' || atomT[2] == 'B' ) && atomT[3] == ' ') nonPolar = true;
  if ( atomT[2] == 'H') nonPolar = true;
  if ( atomT[2] == 'G') {
    if (res == 'N') nonPolar = false;
    nonPolar = true;
  }
  if ( atomT[2] == 'D') {
    if (res == 'Q') nonPolar = false;
    nonPolar = true;
  }
  if ( atomT[2] == 'E') {
    if (res == 'H') nonPolar = false;
    nonPolar = true;
  }
  if ( atomT[2] == 'Z') {
    if (res == 'R') nonPolar = false;
    nonPolar = true;
  }
  if ( atomT[2] == 'Z') {
    if (res == 'R') nonPolar = false;
    nonPolar = true;
  }


  if(fabs(charge_) > 0.4 && nonPolar) {
    //Logger::warningMessage() << "Non Polar atom " << *this << " charge = " << charge_ << endl;
  }
  isNonPolar_ = nonPolar;
}
