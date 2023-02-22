#include "BBGrid.h"

void BBGrid::markResidues(const ChemMolecule &M) {
    std::vector<float> weights;
    weights.insert(weights.end(), maxEntry, 0.0);
    for (Molecule<ChemAtom>::const_iterator it = M.begin(); it != M.end(); it++) {
        float atomRadius = it->getRadius() + radiusAdition_;
        float radius = atomRadius * 2; // may be +1 is enouph
        // cout << " atomRadius " << atomRadius << endl;
        int centerIndex = getIndexForPoint(it->position());
        if (!isValidIndex(centerIndex)) {
            std::cerr << "Error: Point out of grid" << std::endl;
            exit(1);
        }

        int intRadius = getIntGridRadius(radius);
        int radius2 = intRadius * intRadius;

        int i_bound, j_bound, k_bound;
        i_bound = intRadius;
        for (int i = -i_bound; i <= i_bound; i++) {
            j_bound = (int)sqrt(radius2 - i * i);
            for (int j = -j_bound; j <= j_bound; j++) {
                k_bound = (int)sqrt(radius2 - i * i - j * j);
                for (int k = -k_bound; k <= k_bound; k++) {
                    int index = centerIndex + i + xGridNum * j + xyGridNum * k;
                    if (isValidIndex(index) && grid[index] <= 0) {
                        Vector3 point = getPointForIndex(index);
                        float dist = point.dist(it->position());
                        if (dist == 0.0) {
                            if (it->isBackbone())
                                residues[index] = it->residueIndex() * -1;
                            else
                                residues[index] = it->residueIndex();
                            continue;
                        }
                        float weight = atomRadius / dist;
                        if (weight <= weights[index])
                            continue;
                        weights[index] = weight;
                        if (it->isBackbone())
                            residues[index] = it->residueIndex() * -1;
                        else
                            residues[index] = it->residueIndex();
                    }
                }
            }
        }
    }
    weights.clear();
}
