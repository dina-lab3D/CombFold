#ifndef HIERARCHICALFOLD_H
#define HIERARCHICALFOLD_H

#include "BestK.h"
#include "BestKContainer.h"
#include "ComplexDistanceConstraint.h"
#include "BBContainer.h"
#include <future>
#include <memory>

// TODO: should this class just be a namespace? because each function is called once...
class HierarchicalFold {
  public:
    // N - number of subunits, k - best solutions to save at each step
    HierarchicalFold(BBContainer& bbContainer, unsigned int k, unsigned int maxResultPerResSet,
                     float minTemperatureToConsiderCollision, float maxBackboneCollisionPercentPerChain,
                     float penetrationThreshold, float restraintsRatio)
        : countFilterTras_(0), countFilterTrasSkipped_(0), N_(bbContainer.getBBs().size()), K_(k),
          maxResultPerResSet(maxResultPerResSet), minTemperatureToConsiderCollision(minTemperatureToConsiderCollision),
          maxBackboneCollisionPercentPerChain(maxBackboneCollisionPercentPerChain),
          restraintsRatioThreshold_(restraintsRatio), penetrationThreshold_(penetrationThreshold),
          finalSizeLimit_(k * N_), bestKContainer_(k), complexConst_(bbContainer.getBBs()) {

        // initialize keptResultsByLength and bestKContainer_
        keptResultsByLength[1] = new BestK(N_);
        for (unsigned int i = 0; i < bbContainer.getBBsNumber(); i++) {
            std::shared_ptr<SuperBB> sbb = std::make_shared<SuperBB>(bbContainer.getBB(i));

            keptResultsByLength[1]->push(sbb);
            BitId set = BitId(i);
            BestK *cb = bestKContainer_.newBestK(set);
            cb->push(sbb);
        }
    }
    

    void fold(const std::string &outFileNamePrefix);

    void tryToConnect(int id, const SuperBB &sbb1, const SuperBB &sbb2, BestK &results, bool toAdd,
                      std::promise<int> &output, std::vector<std::vector<unsigned int>> &identGroups);
    
    bool filterTrans(const SuperBB &sbb1, const SuperBB &sbb2, const RigidTrans3 &trans) const;

    void createSymmetry(std::vector<std::shared_ptr<SuperBB>> identBBs, BestK &results);

    // utils
    std::shared_ptr<SuperBB> createJoined(const SuperBB &sbb1, const SuperBB &sbb2, RigidTrans3 &trans, int bbPen,
                                          FoldStep &step, float transScore) const;

    std::shared_ptr<SuperBB> getMatchingSBB(SuperBB sbb1, SuperBB sbb2,
                                            std::vector<std::vector<unsigned int>> &identGroups);

    void readConstraints(const std::string fileName) {
        complexConst_.readRestraintsFile(fileName);
        complexConst_.addChainConnectivityConstraints();
    }

    void checkConnectivity() const;

    void outputConnectivityGraph(std::string outFileName = "graph.sif") const;

    // members
    static Timer timer_, timerAll_;
    static unsigned int countResults_;


    mutable unsigned int countFilterTras_;
    mutable unsigned int countFilterTrasSkipped_;

  private:
    const unsigned int N_;                 // number of subunits
    const unsigned int K_;                 // number of solutions to save at each stage
    const unsigned int maxResultPerResSet; // number of solutions to save for each resSet
    const float minTemperatureToConsiderCollision;
    const float maxBackboneCollisionPercentPerChain;
    const float restraintsRatioThreshold_;

    float penetrationThreshold_; // this is ignored for now

    int finalSizeLimit_;
    BestKContainer bestKContainer_;
    std::map<unsigned int, BestK *> keptResultsByLength;
    ComplexDistanceConstraint complexConst_;
};

#endif /* HIERARCHICALFOLD_H */
