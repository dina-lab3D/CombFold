#ifndef FOLDSTEP_H
#define FOLDSTEP_H

#include <algorithm>

class FoldStep {
  public:
    // Constructors
    FoldStep(unsigned int i, unsigned int j, RigidTrans3 trans) : i_(i), j_(j), trans_(trans) {}
    FoldStep(unsigned int i, unsigned int j) : i_(i), j_(j) {}
    FoldStep(unsigned int i, unsigned int j, float score) : i_(i), j_(j), tScore_(score) {}

    // bool operator < (const FoldStep& s2) const { return ((i_ < s2.i_) || (i_ == s2.i_ && j_ < s2.j_)); }

    friend bool operator<(const FoldStep &s1, const FoldStep &s2) {
        unsigned int max1 = std::max(s1.i_, s1.j_);
        unsigned int max2 = std::max(s2.i_, s2.j_);
        return (max1 < max2 || (max1 == max2 && std::min(s1.i_, s1.j_) < std::min(s2.i_, s2.j_)));
    }

    static void outputFoldSteps(std::ostream &out, const std::vector<FoldStep> &foldSteps) {
        out << "foldSteps:";
        for (unsigned int i = 0; i < foldSteps.size(); i++) {
            out << " (" << foldSteps[i].i_ << ", " << foldSteps[i].j_ << ")-" << foldSteps[i].tScore_;
        }
    }

  public:
    unsigned int i_, j_;
    float tScore_;
    RigidTrans3 trans_;
};

#endif /* FOLDSTEP_H */
