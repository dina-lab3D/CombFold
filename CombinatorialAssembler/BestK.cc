#include "BestK.h"


bool BestK::push(std::shared_ptr<SuperBB> in) {
    std::lock_guard<std::mutex> locker(_mu);
    if (size() < k_) {
        insert(in);
        return true;
    } else {

        if (internalMinScore() < score(in)) {
            erase(begin());
            insert(in);
            curMinScore = score(*begin());
            return true;
        }
    }

    return false;
}

bool BestK::push_cluster(std::shared_ptr<SuperBB> in, double rmsd, std::vector<std::vector<unsigned int>> &identGroups) {
    std::lock_guard<std::mutex> locker(_mu);

    if (internalMinScore() > score(in))
        return false;

    for (auto it = begin(); it != end(); it++) {
        if (score(in) <= score(*it) && (*it)->calcRmsd(*in, identGroups) < rmsd)
            return false;
    }

    for (auto it = begin(); it != end();) {
        if (score(in) > score(*it) && (*it)->calcRmsd(*in, identGroups) < rmsd) {
            erase(it++);
        } else {
            ++it;
        }
    }

    if (size() >= k_) {
        erase(begin());
    }

    insert(in);
    curMinScore = score(*begin());
    return true;
}

void BestK::cluster(BestK &clusteredBest, double rmsd, std::vector<std::vector<unsigned int>> &identGroups) const {
    if (size() == 0)
        return;
    const std::shared_ptr<SuperBB> firstSBB = *rbegin();
    if (firstSBB->size() == 2) { // don't cluster
        for (auto it = rbegin(); it != rend(); it++) {
            clusteredBest.insert(*it);
        }
    } else {
        std::vector<bool> clustered(size(), false);
        int i = 0;

        for (auto it = rbegin(); it != rend(); it++, i++) {
            const std::shared_ptr<SuperBB> refSBB = *it;
            // std::cout << "clustering " << i << " clsutered:" << clustered[i] << std::endl;
            if (!clustered[i]) {
                clustered[i] = true;
                clusteredBest.insert(refSBB);
            }

            // TODO: maybe shouldn't cluster more if already clustered (can lead to drift)

            // cluster to other SBBs
            int j = 0;
            // for(auto it2 = begin(); it2!= end(); it2++, j++) {
            for (auto it2 = rbegin(); it2 != rend(); it2++, j++) {
                if (!clustered[j] && (*it2)->calcRmsd(*refSBB, identGroups) < rmsd) {
                    clustered[j] = true;
                }
            }
        }
    }
}