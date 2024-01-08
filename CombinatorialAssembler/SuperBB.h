/**
 * This class represents a SuperBB - A subcomplex composed by multiple subunits (saved in bbs_) and in certain
 * transformations (saved in trans_)
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

    SuperBB(std::shared_ptr<const BB> bb);

    BitId bitIds() const { return bitIDS_; }
    unsigned int size() const { return size_; }
    float getRestraintsRatio() const { return restraintsRatio_; }
    void setRestraintsRatio(float r) { restraintsRatio_ = r; }

    // This function takes Two super BBs and joins them using a transformation between to BBs
    void join(const RigidTrans3 &trans, const SuperBB &other, int bbPen, FoldStep &step, float transScore);
    void replaceIdentBB(BitId oldBBBitId, std::shared_ptr<const BB> bb);

    // This function checks for collissions - Receives transformation and a second super BB and decides if
    // they collide
    bool isPenetrating(const RigidTrans3 &trans, const SuperBB &other, float threshold) const;
    double calcRmsd(const SuperBB &other, std::vector<std::vector<unsigned int>> &identGroups) const;
    double calcRmsd(const SuperBB &other) const;

    void fullReport(std::ostream &s);
    friend std::ostream &operator<<(std::ostream &s, const SuperBB &sbb);

  private:
    // members
    unsigned int size_;
    std::vector<FoldStep> foldSteps_; // FoldStep is a transformation between two SUs
    BitId bitIDS_;
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

#endif /* SUPERBB_H */
