#include "BB.h"

#include <Common.h>
#include <connolly_surface.h>

BB::BB(int id, const std::string pdbFileName, int groupID, const ChemLib &lib, float gridResolution, float gridMargins,
       float minTempFactor)
    : id_(id), groupId_(groupID), pdbFileName_(pdbFileName) {
    // read atoms
    Common::readChemMolecule(pdbFileName_, allAtoms_, lib);
    std::cout << "Done reading ChemMolecule " << allAtoms_.size() << std::endl;
    // read backbone atoms
    std::ifstream pdb2(pdbFileName_);
    backBone_.readPDBfile(pdb2, PDB::BBSelector());
    pdb2.close();
    cm_ = backBone_.centroid();
    numOfAtoms_ = allAtoms_.size();

    // read CA atoms
    std::ifstream pdb3(pdbFileName_);
    caAtoms_.readPDBfile(pdb3, PDB::CAlphaSelector());
    pdb3.close();

    // compute ms surface
    msSurface_ = get_connolly_surface(allAtoms_, 10, 1.8);
    std::cout << "Surface size " << msSurface_.size() << std::endl;

    // compute grid
    grid_ = new BBGrid(msSurface_, gridResolution, gridMargins, 1.5);
    grid_->computeDistFromSurface(msSurface_);
    grid_->markTheInside(allAtoms_);
    grid_->markResidues(backBone_);
    std::cout << "Done compute grid " << pdbFileName_ << std::endl;

    std::vector<Atom *> atomsMap;
    atomsMap.push_back(&(*allAtoms_.begin()));
    for (ChemMolecule::iterator i = allAtoms_.begin(); i != allAtoms_.end(); i++) {
        atomsMap.push_back(&(*i));
    }

    computeFragments(minTempFactor); // get the endpoints

    for (Molecule<Atom>::const_iterator it = caAtoms_.begin(); it != caAtoms_.end(); it++) {
        resIndexToCAAtom[it->residueIndex()] = *it;
    }

    maxRadius_ = 0.0;
    for (Molecule<Atom>::const_iterator it = caAtoms_.begin(); it != caAtoms_.end(); it++) {
        float r = (it->position() - cm_).norm();
        if (r > maxRadius_)
            maxRadius_ = r;
    }
    std::cout << "Max radius: " << maxRadius_ << std::endl;

    std::cout << " done reading BB " << pdbFileName_.c_str() << std::endl;
}

void BB::computeFragments(float minTempFactor) {
    // calculate endpoints
    char currChain;
    int firstResIndex, prevResIndex;
    bool currChainSet = false;
    for (auto i = allAtoms_.begin(); i != allAtoms_.end(); i++) {
        // only CA atoms are considered
        if (!i->isCA())
            continue;
        if (i->getTempFactor() < minTempFactor)
            continue;
        char chain = i->chainId();
        int resIndex = i->residueIndex();
        // one more residue of the same chain - advance
        // if over 20 residues diff - new fragment
        if (currChainSet && currChain == chain && resIndex - prevResIndex <= 20) {
            prevResIndex = resIndex;
        } else {                // new chain
            if (currChainSet) { // save currChain
                ResidueRange range(firstResIndex, prevResIndex);
                fragmentEndpoints_.push_back(std::make_pair(currChain, range));
            }
            // update
            currChain = chain;
            firstResIndex = prevResIndex = resIndex;
            currChainSet = true;
        }
    }
    // save last fragment
    if (currChainSet) { // save currChain
        ResidueRange range(firstResIndex, prevResIndex);
        fragmentEndpoints_.push_back(std::make_pair(currChain, range));
    }

    for (int i = 0; i < (int)fragmentEndpoints_.size(); i++) {
        std::cout << "Fragment " << i << " chainId " << fragmentEndpoints_[i].first << " range "
                  << fragmentEndpoints_[i].second.first << ":" << fragmentEndpoints_[i].second.second << std::endl;
    }
}

void BB::getChainConnectivityConstraints(const BB &otherBB,
                                         std::vector<std::pair<char, std::pair<int, int>>> &constraints) const {
    for (int i = 0; i < (int)fragmentEndpoints_.size(); i++) {
        char chainId1 = fragmentEndpoints_[i].first;
        int resIndex1N = fragmentEndpoints_[i].second.first;
        int resIndex1C = fragmentEndpoints_[i].second.second;

        for (int j = 0; j < (int)otherBB.fragmentEndpoints_.size(); j++) {
            char chainId2 = otherBB.fragmentEndpoints_[j].first;
            if (chainId1 != chainId2)
                continue;
            int resIndex2N = otherBB.fragmentEndpoints_[j].second.first;
            int resIndex2C = otherBB.fragmentEndpoints_[j].second.second;
            if (resIndex1N < resIndex2N) { // add constraint on 1C and 2N
                constraints.push_back(std::make_pair(chainId1, std::make_pair(resIndex1C, resIndex2N)));
            } else { // add constraint on 2C and 1N
                constraints.push_back(std::make_pair(chainId1, std::make_pair(resIndex1N, resIndex2C)));
            }
        }
    }
}

bool BB::isPenetrating(const RigidTrans3 &trans, const BB &other, float threshold) const {
    for (Surface::const_iterator it = other.surface_.begin(); it != other.surface_.end(); it++) {
        float penetration = getDistFromSurface(trans * it->position());
        if (threshold > penetration) {
            // cerr << "trans " << trans << " penetrates " << penetration << endl;
            return true;
        }
    }
    return false;
}

float BB::maxPenetration(const RigidTrans3 &trans, const BB &other) const {
    float max = 1000;
    for (Surface::const_iterator it = other.surface_.begin(); it != other.surface_.end(); it++) {
        float penetration = getDistFromSurface(trans * it->position());
        if (max > penetration) {
            max = penetration;
        }
    }
    return max;
}

bool BB::isIdent(const BB &otherBB) const {
    if (getNumOfAtoms() != otherBB.getNumOfAtoms()) {
        std::cout << "different num of atoms" << std::endl;
        return false;
    }
    for (unsigned int i = 0; i < allAtoms_.size(); i++) {
        if (!(allAtoms_[i].position() - otherBB.allAtoms_[i].position()).isZero()) {
            std::cout << "different atom position" << i << std::endl;
            return false;
        }
    }
    if (groupId_ != otherBB.groupId_) {
        std::cout << "different group id" << std::endl;
        return false;
    }
    std::cout << "checking transforms" << std::endl;
    if (trans_.size() != otherBB.trans_.size())
        return false;
    for (unsigned int i = 0; i < trans_.size(); i++) {
        if (i == id_ || i == otherBB.id_)
            continue;

        if (trans_[i].size() != otherBB.trans_[i].size()) {
            std::cout << "different number of trans " << i << " " << trans_[i].size()
                      << " != " << otherBB.trans_[i].size() << std::endl;
            return false;
        }

        for (unsigned int j = 0; j < trans_[i].size(); j++) {
            if (trans_[i][j]->score() != otherBB.trans_[i][j]->score()) {
                std::cout << "different trans score" << i << std::endl;
                return false;
            }
        }
    }
    return true;
}
