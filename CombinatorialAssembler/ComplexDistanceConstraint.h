#ifndef DistanceConstraint_H
#define DistanceConstraint_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Molecule.h>
// #include <DistanceConstraint.h>
#include <CrossLink.h>
#include <DistanceRestraint.h>

#include "BB.h"

class ComplexDistanceConstraint {
  public:
    ComplexDistanceConstraint(const std::vector<std::shared_ptr<const BB>> &bbs)
        : bbs_(bbs), noOfSUs_(bbs.size()), constraints_(noOfSUs_ * noOfSUs_), restraints_(noOfSUs_ * noOfSUs_),
          restraintIndsToCrosslinkInds_(noOfSUs_ * noOfSUs_) {}

    int readRestraintsFile(const std::string fileName);
    int addChainConnectivityConstraints();

    int numberOfRestraints(int suInd1, int suInd2) const {
        int ind = suInd1 * noOfSUs_ + suInd2;
        return restraints_[ind].size();
    }

    int numberOfConstraints(int suInd1, int suInd2) const {
        int ind = suInd1 * noOfSUs_ + suInd2;
        return constraints_[ind].size();
    }

    // constraints satisfaction: every constraint needs to be satisfied
    bool areConstraintsSatisfied(int suInd1, int suInd2, RigidTrans3 &trans2) const;

    // restraints satisfaction: a predefined ratio needs to be satisfied
    float getRestraintsRatio(const std::vector<std::shared_ptr<const BB>> &bbs,
                             const std::vector<RigidTrans3> &trans) const;

  private:
    void addConstraint(int suInd1, int suInd2, Vector3 receptorAtom, Vector3 ligandAtom, float maxDistance,
                       float minDistance = 0);

    void addRestraint(int suInd1, int suInd2, Vector3 receptorAtom, Vector3 ligandAtom, float maxDistance,
                      float minDistance = 0);

    // find all SUs with residueSequenceID and chains, maxoffset is +/- few residues
    std::set<int> getSUs(const std::string residueSequenceID, const std::string chains, int maxoffset = 0) const;

  private:
    const std::vector<std::shared_ptr<const BB>> &bbs_;
    int noOfSUs_;
    std::vector<std::vector<DistanceRestraint>> constraints_;
    std::vector<std::vector<DistanceRestraint>> restraints_;
    int numberOfConstraints_ = 0;
    int numberOfRestraints_ = 0;
    std::vector<std::vector<unsigned int>> restraintIndsToCrosslinkInds_;
    std::vector<float> crosslinkIndToWeight_;
};

#endif
