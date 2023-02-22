#include "ClusterSBBs.h"
#include "ChemMolecule.h"
#include "Common.h"
#include "Match.h"
#include "SubSet.h"

bool ClusterSBBs::addSBB(SuperBB *sbb) {
    sbb->preClustering();
    unsigned int numOfBBs = sbb->size();
    char *clustStr = new char[numOfBBs * numOfBBs + 1];
    for (unsigned int i = 0; i < sbb->contact_.size(); ++i)
        clustStr[i] = sbb->contact_[i] ? '1' : '0';
    clustStr[numOfBBs * numOfBBs] = '\0';
    std::unordered_map<std::string, std::vector<Cluster>>::iterator it = contactHash_.find(clustStr);
    if (it == contactHash_.end()) { // new contact map, open a new cluster
        contactHash_[clustStr].push_back(sbb);
        size_++;
        return true;
    }
    // else has some with same contact list
    return addSBB(it->second, sbb);
}

bool ClusterSBBs::addSBB(std::vector<Cluster> &vec, SuperBB *sbb) {
    // try to add to existing cluster
    for (std::vector<Cluster>::iterator it = vec.begin(); it != vec.end(); it++) {
        if (withinDist(*it->representive_, *sbb, distance_)) {
            return it->add(sbb);
        }
    }
    // open a new cluster
    vec.push_back(sbb);
    size_++;
    return true;
}

bool ClusterSBBs::Cluster::add(SuperBB *sbb) {
    count_++;
    if (sbb->newScore_ > representive_->newScore_) {
        // if (sbb->transScore_ > representive_->transScore_) {
        SuperBB *tmp = representive_;
        representive_ = sbb;
        max(*representive_, *tmp);
        return true;
    }
    max(*representive_, *sbb);
    return false;
}

void ClusterSBBs::reportAll(std::ostream &s) {
    for (std::vector<Cluster *>::iterator it = finalClusters_.begin(); it != finalClusters_.end(); it++) {
        s << "within_Cluster: " << (*it)->count_ << " ";
        (*it)->representive_->fullReport(s);
    }
    std::cerr << " finalClusters_.size() " << finalClusters_.size() << " size_ " << size_ << " contactHash_.size() "
              << contactHash_.size() << " numOfRmsdCMs_ " << numOfRmsdCMs_ << " numOfRmsdAll_ " << numOfRmsdAll_
              << std::endl;
}

bool ClusterSBBs::withinDist(const SuperBB &sbb1, const SuperBB &sbb2, float distance) const {
    // 1. calculate RMSD over centers of mass
    Match match1;
    Molecule<Atom> A1, B1;
    int j = 0;
    std::vector<int> mapSBBS(sbb1.size());
    for (unsigned int i = 0; i < sbb1.size(); i++) {
        mapSBBS[sbb1.bbs_[i]->getID()] = i;
    }

    for (unsigned int i = 0; i < sbb2.size(); i++) {
        const RigidTrans3 &tr1 = sbb1.trans_[mapSBBS[sbb2.bbs_[i]->getID()]];
        const RigidTrans3 &tr2 = sbb2.trans_[i];
        Atom a(tr1 * sbb2.bbs_[i]->getCM(), j);
        Atom b(tr2 * sbb2.bbs_[i]->getCM(), j);
        A1.add(a);
        B1.add(b);
        match1.add(j, j, 1);
        j++;
    }

    match1.calculateBestFit(A1, B1);
    float rmsd1 = match1.rmsd();
    numOfRmsdCMs_++;
    if (rmsd1 > distance) {
        return false;
    }

    // 2. calculate rmsd over CA
    Match match;
    Molecule<Atom> A, B;
    j = 0;
    for (unsigned int i = 0; i < sbb2.size(); i++) {
        const Molecule<ChemAtom> &backBone = sbb2.bbs_[i]->backBone_;
        const RigidTrans3 &tr1 = sbb1.trans_[mapSBBS[sbb2.bbs_[i]->getID()]];
        const RigidTrans3 &tr2 = sbb2.trans_[i];
        for (Molecule<ChemAtom>::const_iterator it = backBone.begin(); it != backBone.end(); it++) {
            if (it->type()[1] == 'C' && it->type()[2] == 'A') {
                Atom a(*it);
                a *= tr1;
                Atom b(*it);
                b *= tr2;
                A.add(a);
                B.add(b);
                match.add(j, j, 1);
                j++;
            }
        }
    }

    match.calculateBestFit(A, B);
    float rmsd = match.rmsd();
    numOfRmsdAll_++;
    return (rmsd <= distance * 1.5);
}

void ClusterSBBs::merge() {
    for (auto it = contactHash_.begin(); it != contactHash_.end(); it++) {
        for (std::vector<Cluster>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
            bool toPush = true;
            for (std::vector<Cluster *>::iterator it3 = finalClusters_.begin(); it3 != finalClusters_.end(); it3++) {
                if ((*it3)->representive_->contact_ != it2->representive_->contact_ &&
                    withinDist(*(*it3)->representive_, *it2->representive_, distance_)) {
                    (*it3)->add(it2->representive_);
                    (*it3)->count_ += it2->count_ - 1;
                    toPush = false;
                    break;
                }
            }
            // it2->representive_->calcFinalScore();
            if (toPush) {
                finalClusters_.push_back((Cluster *)&(*it2));
            }
        }
    }
    std::cerr << " size_ " << size_ << " contactHash_.size() " << contactHash_.size() << " numOfRmsdCMs_ "
              << numOfRmsdCMs_ << " numOfRmsdAll_ " << numOfRmsdAll_ << std::endl;
}
