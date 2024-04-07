#include "TransformDistance.h"

float TransformDistance::rmsd2(const RigidTrans3& trans1, const RigidTrans3& trans2) const{
    // compute rotation matrix
    Matrix3 A = trans1.rotation().transposeMatrix() * trans2.rotation();

    float sum = 0;
    for (unsigned short i = 0; i < 3; i++) {
        for (unsigned short j = i; j < 3; j++) {
            if (i==j) sum += (float)((2.0 - 2.0*A[i][i])* Xij_[i][i]);
            else sum += (float)(2.0*( - A[i][j] - A[j][i]) * Xij_[i][j]);
        }
    }

    // (t1 - t2)^2
    Vector3 t1_t2 = trans1.translation()-trans2.translation();
    sum += t1_t2.norm2();

    // 2(t1-t2)(R1-R2)centroid
    sum += 2 * t1_t2 * ((trans1.rotation() - trans2.rotation()) * centroid_);

    return sum;
}

float TransformDistance::rmsd2(const RigidTrans3& trans) {
    // compute rotation matrix
    Matrix3 A = trans.rotation();

    float sum = 0;
    for (unsigned short i = 0; i < 3; i++) {
        for (unsigned short j = i; j < 3; j++) {
            if (i==j) sum += (float)((2.0 -2.0*A[i][i])* Xij_[i][i]);
            else sum += (float)(2*( - A[i][j] - A[j][i]) * Xij_[i][j]);
        }
    }

    // (t1 - t2)^2
    Vector3 t1_t2 = -trans.translation();
    sum += t1_t2.norm2();

    // 2(t1-t2)(R1-R2)centroid
    Matrix3 identity(1);
    sum += 2 * t1_t2 * ((identity - trans.rotation()) * centroid_);

    return sum;
}