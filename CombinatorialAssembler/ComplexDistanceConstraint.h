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
        : bbs_(bbs), noOfSUs_(bbs.size()), constraints_(noOfSUs_ * noOfSUs_), restraints_(noOfSUs_ * noOfSUs_) {}

    int readRestraintsFile(const std::string fileName);
    int addChainConnectivityConstraints();

    float getDistanceRestraintsRatio() const { return distanceRestraintsRatio_; }
    void setDistanceRestraintsRatio(float r) { distanceRestraintsRatio_ = r; }

    int numberOfRestraints(int suInd1, int suInd2) const {
        int ind = suInd1 * noOfSUs_ + suInd2;
        return restraints_[ind].size();
    }

    int numberOfConstraints(int suInd1, int suInd2) const {
        int ind = suInd1 * noOfSUs_ + suInd2;
        return constraints_[ind].size();
    }

    // constraints satisfaction: every constraint needs to be satisfied
    bool isSatisfied(int suInd1, int suInd2, RigidTrans3 &trans1, RigidTrans3 &trans2) const {
        int ind = suInd1 * noOfSUs_ + suInd2;
        for (unsigned int index = 0; index < constraints_[ind].size(); index++) {
            bool sat = constraints_[ind][index].isSatisfied(trans1, trans2);
            if (!sat) {
                return false;
            }
        }
        return true;
    }

    // restraints satisfaction: a predefined ratio needs to be satisfied
    bool isRatioSatisfied(int suInd1, int suInd2, RigidTrans3 &trans1, RigidTrans3 &trans2) const {
        int ind = suInd1 * noOfSUs_ + suInd2;
        int count = 0;
        if (restraints_[ind].size() > 0) {
            for (unsigned int index = 0; index < restraints_[ind].size(); index++) {
                if (restraints_[ind][index].isSatisfied(trans1, trans2))
                    count++;
            }
            float ratio = (float)count / restraints_[ind].size();
            if (ratio < distanceRestraintsRatio_)
                return false;
            return true;
        }
        return true;
    }

    // constraints satisfaction: every constraint needs to be satisfied
    bool isSatisfied(int suInd1, int suInd2, RigidTrans3 &trans2) const {
        int ind = suInd1 * noOfSUs_ + suInd2;
        for (unsigned int index = 0; index < constraints_[ind].size(); index++) {
            bool sat = constraints_[ind][index].isSatisfied(trans2);
            if (!sat) {
                return false;
            }
        }
        return true;
    }

    // restraints satisfaction: a predefined ratio needs to be satisfied
    bool isRatioSatisfied(int suInd1, int suInd2, RigidTrans3 &trans2) const {
        int ind = suInd1 * noOfSUs_ + suInd2;
        int count = 0;
        if (restraints_[ind].size() > 0) {
            for (unsigned int index = 0; index < restraints_[ind].size(); index++) {
                if (restraints_[ind][index].isSatisfied(trans2))
                    count++;
            }
            float ratio = (float)count / restraints_[ind].size();
            if (ratio < distanceRestraintsRatio_)
                return false;
            return true;
        }
        return true;
    }

    float getRatio(const std::vector<std::shared_ptr<const BB>> &bbs, const std::vector<RigidTrans3> &trans) const {

        unsigned int count = 0;
        unsigned int total = 0;
        for (int suInd1 = 0; suInd1 < (int)bbs.size(); suInd1++) {
            for (int suInd2 = 0; suInd2 < (int)bbs.size(); suInd2++) {
                int index1 = bbs[suInd1]->getID();
                int index2 = bbs[suInd2]->getID();
                if (index1 == index2)
                    continue;
                int ind = index1 * noOfSUs_ + index2;
                total += restraints_[ind].size();
                for (unsigned int index = 0; index < restraints_[ind].size(); index++) {
                    if (restraints_[ind][index].isSatisfied(trans[suInd1], trans[suInd2]))
                        count++;
                }
            }
        }
        float ratio = 1.0;
        if (total > 0)
            ratio = (float)count / total; //(2*numberOfRestraints_);
        return ratio;
    }

    // restraints satisfaction: a predefined ratio needs to be satisfied
    bool isRatioSatisfied(const std::vector<std::shared_ptr<const BB>> &bbs,
                          const std::vector<RigidTrans3> &trans) const {
        float ratio = getRatio(bbs, trans);
        if (ratio < distanceRestraintsRatio_)
            return false;
        return true;
    }

    bool isSatisfied(const std::vector<std::shared_ptr<const BB>> &bbs, const std::vector<RigidTrans3> &trans) const {
        if (((int)trans.size()) != noOfSUs_) {
            std::cerr << " trans.size() != noOfSUs_ " << std::endl;
            return false;
        }
        for (int suInd1 = 0; suInd1 < noOfSUs_; suInd1++) {
            for (int suInd2 = 0; suInd2 < noOfSUs_; suInd2++) {
                int index1 = bbs[suInd1]->getID();
                int index2 = bbs[suInd2]->getID();
                if (index1 == index2)
                    continue;
                int ind = index1 * noOfSUs_ + index2;
                for (unsigned int index = 0; index < constraints_[ind].size(); index++) {
                    bool sat = constraints_[ind][index].isSatisfied(trans[suInd1], trans[suInd2]);
                    if (!sat) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

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
    float distanceRestraintsRatio_ = 1.0; // required ratio of satisfied distances
    int numberOfConstraints_ = 0;
    int numberOfRestraints_ = 0;
    std::vector<CrossLink> crossLinks_;
};

#endif
