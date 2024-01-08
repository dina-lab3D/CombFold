#ifndef BB_H
#define BB_H

#include "BBGrid.h"
#include "TransformationAndScore.h"
#include "BitId.h"

#include <ChemMolecule.h>
#include <GeomScore.h>
#include <RigidTrans3.h>
#include <Surface.h>
#include <Vector3.h>

#include <vector>

class BBConstructor {
    // idea taken from
    // https://stackoverflow.com/questions/3220009/is-this-key-oriented-access-protection-pattern-a-known-idiom
  private:
    friend class BBContainer;
    BBConstructor() {}
};

class BB {
  public:
    friend class SuperBB;

    BB(int id, const std::string pdbFilename, int groupID, const ChemLib &lib, float gridResolution, float gridMargins,
       float minTempFactor);

    // access
    int getID() const { return id_; }
    BitId bitId() const { return BitId(id_); }
    int groupId() const { return groupId_; }
    std::string getPDBFileName() const { return pdbFileName_; }

    unsigned int getNumOfAtoms() const { return numOfAtoms_; }
    const Vector3 &getCM() const { return cm_; }
    const float getRadius() const { return maxRadius_; }

    const ChemAtom &getChemAtom(int atomIndex) const { return allAtoms_.getChemAtom(atomIndex); }
    const ChemAtom &getChemAtomByIndex(int atomIndex) const { return allAtoms_[atomIndex]; }
    const Atom &getAtomByResId(unsigned int resId) const { return resIndexToCAAtom.at(resId); }

    unsigned int getSurfaceSize() const { return surface_.size(); }
    float getDistFromSurface(const Vector3 &v) const { return grid_->getDist(v); }


    // This uses BBConstructor to make sure that only BBContainer can call this
    void putTransWith(int bbIndex, const std::shared_ptr<TransformationAndScore> &t1, const BBConstructor &) const {
        trans_[bbIndex].push_back(t1);
    }
    void initTrans(unsigned int numberOfBBs, const BBConstructor &) const {
        trans_.insert(trans_.begin(), numberOfBBs, std::vector<std::shared_ptr<TransformationAndScore>>());
    }

    // TODO: do we still need isPenetrating/maxPenetration ?
    bool isPenetrating(const RigidTrans3 &trans, const BB &other, float threshold) const;
    float maxPenetration(const RigidTrans3 &trans, const BB &other) const;

    void getChainConnectivityConstraints(const BB &bb,
                                         std::vector<std::pair<char, std::pair<int, int>>> &) const; // update

    const std::vector<std::shared_ptr<TransformationAndScore>> &getTransformations(int bbIndex) const {
        return trans_[bbIndex];
    }

    bool isIdent(const BB &otherBB) const;
    
  private:
    // after BB is initialized, compute chains and fragment ranges
    void computeFragments(float minTempFactor);

  private:
    // surface points
    // Different methods for computing collision
    Surface surface_;                  // shuo
    Surface msSurface_;                // connolly - dense

    // BB id
    unsigned int id_;

    // transformations to other BBs
    // This is the edge in a graph - The result of a patch dock calculation
    mutable std::vector<std::vector<std::shared_ptr<TransformationAndScore>>> trans_;

    int groupId_;
    std::string pdbFileName_;
    int numOfAtoms_;

    Vector3 cm_;
    float maxRadius_;

    // chain id and residue numbers for each chain fragment in BB
    typedef std::pair<int, int> ResidueRange;
    std::vector<std::pair<char, ResidueRange>> fragmentEndpoints_;

  public: // TODO
    BBGrid *grid_;
    ChemMolecule backBone_;
    ChemMolecule allAtoms_;
    Molecule<Atom> caAtoms_;
    std::map<unsigned int, Atom> resIndexToCAAtom;
};

#endif /* BB_H */
