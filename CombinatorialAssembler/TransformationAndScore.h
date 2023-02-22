#ifndef TRANSFORMATIONANDSCORE_H
#define TRANSFORMATIONANDSCORE_H

#include <RigidTrans3.h>
#include <vector>

#define NO_OF_RANGES 8
typedef struct pointsAndSurface_t {
    int count_;
    float surface_;
} pointsAndSurface;

template <class SCORE_T> class TransformationAndScore_T {
  public:
    const RigidTrans3 &trans() { return refFrame_; }
    float score() { return score_.totalScore_; }

    RigidTrans3 refFrame_;
    SCORE_T score_;
    float dist_; // for debug
    static void outputTrans(std::vector<TransformationAndScore_T *> &transformations);
};

template <class SCORE_T> std::ostream &operator<<(std::ostream &s, const TransformationAndScore_T<SCORE_T> &ts) {
    return s << ts.refFrame_ << " " << ts.score_ << " , " << ts.dist_ << " , " << ((ts.dist_ == 0) ? 1 : 0);
}

template <class SCORE_T> std::istream &operator>>(std::istream &s, TransformationAndScore_T<SCORE_T> &ts) {
    char c;
    return s >> c >> ts.score_ >> c >> c >> ts.refFrame_ >> c >> c >> ts.dist_;
}

class Score {
  public:
    Score() { init(); }

    void init() {
        totalScore_ = 0;
        resCount1_ = resCount2_ = 0;
        interfaceSurface_ = maxPenetrate_ = s_ = c_ = e_ = 0.0;
        for (int i = 0; i < NO_OF_RANGES; i++) {
            ps_[i].count_ = 0;
            ps_[i].surface_ = 0;
        }
    }

    float totalScore_;
    pointsAndSurface ps_[NO_OF_RANGES];
    int resCount1_;
    int resCount2_;
    float interfaceSurface_;
    float maxPenetrate_;
    float s_, c_, e_;
    Score *referenceScore_;

    friend std::ostream &operator<<(std::ostream &s, const Score &score) { return s << score.totalScore_; }

    friend std::istream &operator>>(std::istream &s, Score &ts);
};

// #define CLUSTER_DEBUG
typedef TransformationAndScore_T<Score> TransformationAndScore;

#endif
