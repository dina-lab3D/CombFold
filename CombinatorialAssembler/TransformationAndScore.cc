#include "TransformationAndScore.h"

template <>
void TransformationAndScore_T<Score>::outputTrans(std::vector<TransformationAndScore_T<Score> *> &transformations) {
    int count = 0;
    for (std::vector<TransformationAndScore *>::iterator it = transformations.begin(); it != transformations.end();
         it++) {
        count++;
        std::cout << count << ": " << **it << std::endl;
    }
}
