/**
 * \file BB.h
 * \brief
 *
 * \authors Yuval Inbar, Dina Schneidman
 *
 */
#ifndef BB_H
#define BB_H

#include "BBGrid.h"
#include "ScoreParams.h"
#include "TransformationAndScore.h"

#include <ChemMolecule.h>
#include <GeomScore.h>
#include <RigidTrans3.h>
#include <Surface.h>
#include <Vector3.h>

#include <vector>

#define PROB_RADIUS 1.8

class InterestPoint {
  public:
    InterestPoint(Vector3 position, Vector3 probePosition, bool hydrophobic)
        : position_(position), probePosition_(probePosition), hydrophobic_(hydrophobic) {}
    Vector3 position_;
    Vector3 probePosition_;
    bool hydrophobic_;
};

class SurfaceAtom {
  public:
    Vector3 pos_;
    float radius_;
};

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

    BB(int id, const std::string pdbFilename, int groupID, int groupSize, const ChemLib &lib);

    // access
    int getID() const { return id_; }
    unsigned long bitId() const { return 1 << id_; }
    int groupId() const { return groupId_; }

    std::string getPDBFileName() const { return pdbFileName_; }

    unsigned int getSurfaceSize() const { return surface_.size(); }
    unsigned int getNumOfAtoms() const { return numOfAtoms_; }
    const Vector3 &getCM() const { return cm_; }
    const float getRadius() const { return maxRadius_; }
    const ChemAtom &getChemAtom(int atomIndex) const { return allAtoms_.getChemAtom(atomIndex); }
    const ChemAtom &getChemAtomByIndex(int atomIndex) const { return allAtoms_[atomIndex]; }

    float getDistFromSurface(const Vector3 &v) const { return grid_->getDist(v); }

    const Atom &getAtomByResId(unsigned int resId) const;

    // float scoreForTrans(Score &score, RigidTrans3 trans, BB &other);
    bool isPenetrating(const RigidTrans3 &trans, const BB &other, float threshold) const;
    float maxPenetration(const RigidTrans3 &trans, const BB &other) const;

    void getChainConnectivityConstraints(const BB &bb,
                                         std::vector<std::pair<char, std::pair<int, int>>> &) const; // update

    const std::vector<std::shared_ptr<TransformationAndScore>> &getTransformations(int bbIndex) const {
        return trans_[bbIndex];
    }

    // a range of SuperBB sizes for which the transformations can be used
    const std::pair<int, int> &getTransformsRange(int bbIndex) const { return transformationRanges_[bbIndex]; }

    bool isIdent(const BB &otherBB) const;
    // updates

    void computeSurfaceAtoms(Surface &surface);

    void putTransWith(int bbIndex, const std::shared_ptr<TransformationAndScore> &t1, const BBConstructor &) const {
        trans_[bbIndex].push_back(t1);
    }

    void initTrans(unsigned int numberOfBBs, const BBConstructor &) const {
        trans_.insert(trans_.begin(), numberOfBBs, std::vector<std::shared_ptr<TransformationAndScore>>());
    }

    static float penetrationThreshold_;

    void computeRanges(const std::vector<std::shared_ptr<const BB>> &bbs, int margin) const;

  private:
    // after BB is initialized, compute chains and fragment ranges
    void computeFragments();

  private:
    // surface points
    // Different methods for computing collision
    Surface surface_;                  // shuo
    Surface msSurface_;                // connolly - dense

    // BB id
    int id_;

    // transformations to other BBs
    // This is the edge in a graph - The result of a patch dock calculation
    mutable std::vector<std::vector<std::shared_ptr<TransformationAndScore>>> trans_;

    // group id and size of the BB
    int groupId_;
    int groupSize_;

    // a range of SuperBB sizes for which the transformations can be used
    mutable std::vector<std::pair<int, int>> transformationRanges_;

    Vector3 maxP_, minP_;
    std::vector<SurfaceAtom> sas_;
    std::string pdbFileName_;
    // Triangle t_;
    int numOfAtoms_;

    Vector3 cm_;
    float maxRadius_;

    // chain id and residue numbers for each chain fragment in BB
    typedef std::pair<int, int> ResidueRange;
    std::vector<std::pair<char, ResidueRange>> fragmentEndpoints_;

  public: // TODO
    GeomScore *dockScore_;
    BBGrid *grid_;
    ChemMolecule backBone_;
    ChemMolecule allAtoms_;
    Molecule<Atom> caAtoms_;
    std::map<unsigned int, Atom> resIndexToCAAtom;

  public:
    static float gridResolution_;
    static float gridMargins_;
    static ScoreParams scoreParams_;
};

#endif /* BB_H */
