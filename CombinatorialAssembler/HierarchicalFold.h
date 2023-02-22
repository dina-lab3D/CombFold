/**
 * \file HierarchicalFold.h
 * \brief
 *
 * \authors Yuval Inbar, Dina Schneidman
 *
 */
#ifndef HIERARCHICALFOLD_H
#define HIERARCHICALFOLD_H

#include "BestK.h"
#include "BestKContainer.h"
#include "ComplexDistanceConstraint.h"
#include <future>
#include <memory>

class HierarchicalFold {
public:
  // N - number of subunits, k - best solutions to save at each step
  HierarchicalFold(const std::vector<std::shared_ptr<const BB> >& bbs, unsigned int k,
                   unsigned int maxResultPerResSet, float minTemperatureToConsiderCollision, float maxBackboneCollisionPercentPerChain) :
    N_(bbs.size()),K_(k),  maxResultPerResSet(maxResultPerResSet), minTemperatureToConsiderCollision(minTemperatureToConsiderCollision),
    maxBackboneCollisionPercentPerChain(maxBackboneCollisionPercentPerChain), finalSizeLimit_(k * N_), clusteredSBBS_(k), complexConst_(bbs)
    {
    unsigned long max = numberOfSubsets(N_); // 2^N
    std::cout << "FinalSizeLimit_ " << finalSizeLimit_ << " max " << max << std::endl;
  }

  void tryToConnect(int id, const SuperBB &sbb1, const SuperBB &sbb2,
                    BestK &results, bool toAdd, std::promise<int> &output,
                    std::vector<std::vector<unsigned int>>& identGroups);

  void createSymmetry(std::vector <std::shared_ptr<SuperBB>> identBBs, BestK &results);

  std::shared_ptr<SuperBB> createJoined(const SuperBB &sbb1, const SuperBB &sbb2, RigidTrans3 &trans,
                                        int numBurried, int bbPen, float newScore, FoldStep &step, int transScore) const;

  bool filterTrans(const SuperBB &sbb1, const SuperBB &sbb2, const RigidTrans3 &trans,
                   float penThreshold, int &bbPenetrations, int maxBBPenetrations, float &newScore,
                   unsigned int sbb1Index, unsigned int sbb2Index) const;

  void fold(const std::string& outFileNamePrefix);

  // 2^N
  static unsigned long numberOfSubsets(unsigned int N) {
    unsigned long ret = 1;
    for (unsigned int i = 0; i < N; i++) ret *= 2;
    return ret;
  }

  void readConstraints(const std::string fileName, float distanceRestraintsRatio) {
    complexConst_.readRestraintsFile(fileName);
    complexConst_.addChainConnectivityConstraints();
    complexConst_.setDistanceRestraintsRatio(distanceRestraintsRatio);
  }

  // add a single BB
  void add(std::shared_ptr<SuperBB> sbb, int index) {
    unsigned long one = 1;
    unsigned long set = one << index;
    BestK* cb = clusteredSBBS_.newBestK(set);
    cb->push(sbb);
  }

  void checkConnectivity() const;

  void outputConnectivityGraph(std::string outFileName = "graph.sif") const;

  std::string findMissingSubunits(unsigned long) const;

  void outputSubsets() const;

  static bool isValid(const SuperBB &sbb1, const SuperBB &sbb2,
                      unsigned int i, unsigned int j) {

    if (sbb1.size() == 1 && sbb2.size() == 1) return true;

    FoldStep foldstep(i, j);
    unsigned int min = std::min(sbb1.size(), sbb2.size());

    if (sbb1.size() > 1) {
      if (!(sbb1.isInChild1(i) || sbb1.isInChild2(i)))
        std::cout << "OOPS not in children" << std::endl;
      unsigned int min1 = sbb1.isInChild1(i) ? std::min(sbb2.size() + sbb1.size1(), sbb1.size2()) :
        std::min(sbb2.size() + sbb1.size2(), sbb1.size1());
      FoldStep foldstep1 = sbb1.lastFoldStep();

      if ((min1 > min || min1 == min) && (foldstep < foldstep1)) return false;
    }

    if (sbb2.size() > 1) {
      if (!(sbb2.isInChild1(j) || sbb2.isInChild2(j)))
        std::cout << "OOPS not in children" << std::endl;
      unsigned int min2 = sbb2.isInChild1(j) ? std::min(sbb1.size() + sbb2.size1(), sbb2.size2()) :
        std::min(sbb1.size() + sbb2.size2(), sbb2.size1());
      if ((min2 > min || min2 == min) && foldstep < sbb2.lastFoldStep()) return false;
    }
    return true;
  }

  // members
  static Timer timer_, timerAll_;
  static unsigned int countResults_;

private:
  const unsigned int N_; // number of subunits
  const unsigned int K_; // number of solutions to save at each stage
  const unsigned int maxResultPerResSet; // number of solutions to save for each resSet
  const float minTemperatureToConsiderCollision;
  const float maxBackboneCollisionPercentPerChain;


  int finalSizeLimit_;
  BestKContainer clusteredSBBS_;
  ComplexDistanceConstraint complexConst_;
};

#endif /* HIERARCHICALFOLD_H */