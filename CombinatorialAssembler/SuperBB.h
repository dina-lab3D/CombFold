/**
 * \file SuperBB.h
 * \brief
 *
 * \authors Yuval Inbar, Dina Schneidman
 *
 *This class represents a SuperBB - A subunit of the final protein structure.Can be composed
 *of several smaller subunits.The imporant data members are bbs_ and trans_.
 *The important function is join(the larger one)
 *Comment written by Aharon Pritsker
 */
#ifndef SUPERBB_H
#define SUPERBB_H

#include "BB.h"
#include "FoldStep.h"

class SuperBB {
  public:
    friend class BestK;
    friend class HierarchicalFold;
    friend class TransIterator2;

    SuperBB(std::shared_ptr<const BB> bb, unsigned int numberOfBBs);
    // float calcScoreForTrans(Score &score, RigidTrans3 trans, SuperBB &other);
    // This function takes Two super BBs and joins them
    // trans - This is the transformation that joins the two super BBS together
    // other - This is the second BB to be joined to this one
    void join(const RigidTrans3 &trans, const SuperBB &other, int bbPen, FoldStep &step, float transScore);
    // This function checks for collissions - Receives transformation and a second super BB and decides if
    // they collide
    bool isPenetrating(const RigidTrans3 &trans, const SuperBB &other, float threshold) const;
    void calcFinalScore();
    double calcRmsd(const SuperBB &other) const;
    double calcRmsd(const SuperBB &other, std::vector<std::vector<unsigned int>> &identGroups) const;
    void fullReport(std::ostream &s);
    const FoldStep &lastFoldStep() const { return foldSteps_[size_ - 2]; }
    bool isInChild1(unsigned int i) const { return (child1_ & (1 << i)); }
    bool isInChild2(unsigned int i) const { return (child2_ & (1 << i)); }
    int gap(int id);
    float getRestraintsRatio() const { return restraintsRatio_; }

    unsigned int bitIds() const { return bitIDS_; }
    unsigned int size() const { return size_; }
    unsigned int size1() const { return size1_; }
    unsigned int size2() const { return size2_; }

    void replaceIdentBB(unsigned long oldBBBitId, std::shared_ptr<const BB> bb);

    void setSize(unsigned int s) { size_ = s; }
    void setRestraintsRatio(float r) { restraintsRatio_ = r; }

    friend std::ostream &operator<<(std::ostream &s, const SuperBB &sbb);

  private:
    // members
    unsigned int size_, size1_, size2_;
    unsigned int id_;
    unsigned int numOfAtoms_;
    std::vector<FoldStep> foldSteps_; // FoldStep is a transformation between two SUs
    unsigned int bitIDS_;
    unsigned int child1_, child2_;
    unsigned int numberOfBBs_;
    float restraintsRatio_;

  public: // TODO: Make private
    // BBs that make up the SuperBB
    std::vector<std::shared_ptr<const BB>> bbs_;
    // Holds the transformations of the BBs
    std::vector<RigidTrans3> trans_;

    int backBonePen_;
    float maxPen_;
    float transScore_;
    float weightedTransScore_;
    int multPen_;
    int singlePen_;
};

class TransIterator2 {
  public:
    TransIterator2(const SuperBB &sbb1, const SuperBB &sbb2, int pbb1, int pbb2);

    RigidTrans3 &transformation() { return transformation_; }

    float getScore() { return (*it_)->score_.totalScore_; }

    void generateTransformation() { transformation_ = bb1_.trans_[index1_] * (*it_)->refFrame_ * mediatorTrans_; }

    TransIterator2 &operator++(int) {
        while (!isAtEnd()) {
            it_++;
            if (it_ == trans_->end()) {
                return *this;
            }
            generateTransformation();
            return *this;
        }
        return *this;
    }

    bool isAtEnd() { return it_ == trans_->end(); }

  private:
    const SuperBB &bb1_, &bb2_;
    unsigned int index1_, index2_;
    const std::vector<std::shared_ptr<TransformationAndScore>> *trans_;
    std::vector<std::shared_ptr<TransformationAndScore>>::const_iterator it_;
    RigidTrans3 mediatorTrans_;
    RigidTrans3 transformation_;
};

struct compRestraints {
    bool operator()(const std::shared_ptr<SuperBB> &lhs, const std::shared_ptr<SuperBB> &rhs) const {
        return lhs->getRestraintsRatio() < rhs->getRestraintsRatio();
    }
};

#endif /* SUPERBB_H */
