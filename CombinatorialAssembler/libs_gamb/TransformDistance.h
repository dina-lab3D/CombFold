#ifndef TRANSFORMDISTANCE_H
#define TRANSFORMDISTANCE_H

#include "RigidTrans3.h"
#include <vector>

/*
CLASS
  TransformDistance

  Class that allows to compute rmsd between two transformations in constant time

KEYWORDS

  Distance, transformations, rmsd

AUTHORS
  Dina, Yuval (duhovka , inbaryuv @tau.ac.il)

GOALS

  Given a set of points P,  and two transformations T1 and T2, the class computes
  the RMSD between the point set after applying the transformations: RMSD (T1(P), T2(P))
  The trivial way to do this, requires O(n) time. However with some preprocessing
  it can be done in O(1). This class preprocess the point set and afterwards each RMSD
  query is answered in constant time.

USAGE

  Molecule<Atom> mol;

  TransformDistance transformDistance(mol);

  RigidTrans3 t1,t2;

  float rmsd = transformDistance.rmsd(t1, t2);
*/
class TransformDistance {
public:

    // GROUP: Constructors.

    //// empty constructor
    TransformDistance() {}

    //// init point set
    template <class ParticleT>
    TransformDistance(const std::vector<ParticleT>& pointSet);

    //// returns rmsd between two transformations
    float rmsd(const RigidTrans3& trans1, const RigidTrans3& trans2) { return sqrt(rmsd2(trans1, trans2)); }

    //// returns squared rmsd between two transformations
    float rmsd2(const RigidTrans3& trans1, const RigidTrans3& trans2);

    //// returns rmsd between point set and point set after applying the transformation
    float rmsd(const RigidTrans3& trans) { return sqrt(rmsd2(trans)); }

    //// returns squared rmsd between point set and point set after applying the transformation
    float rmsd2(const RigidTrans3& trans);

protected:
    // The centroid of the molecule
    Vector3 centroid_;

    // partial multiplications for fast RMSD computation
    float Xij_[3][3];
};

template <class ParticleT>
TransformDistance::TransformDistance(const std::vector<ParticleT>& particles) {
    // init
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            Xij_[i][j]=0;
        }
    }
    centroid_ = Vector3(0,0,0);

    for (unsigned int k = 0 ; k < particles.size() ; k++) {
        // centroid
        centroid_ += particles[k];

        // partial multiplications
        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
                Xij_[i][j] += (particles[k][i]*particles[k][j]);
            }
        }
    }

    // divide by n
    centroid_ /= particles.size();
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            Xij_[i][j] /= particles.size();
        }
    }
}

#endif //TRANSFORMDISTANCE_H