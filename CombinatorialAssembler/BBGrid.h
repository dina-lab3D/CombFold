/**
 * \file BBGrid.h
 * \brief
 *
 * \authors Dina Schneidman
 *
 *
 */
#ifndef BBGRID_H
#define BBGRID_H

#include <ChemMolecule.h>
#include <prGrid.h>

class BBGrid : public ResidueGrid {
  public:
    BBGrid(const Surface &surface, const float inDelta, const float maxRadius, float radiusAdition)
        : ResidueGrid(surface, inDelta, maxRadius), radiusAdition_(radiusAdition){};
    void markResidues(const ChemMolecule &M);

  private:
    float radiusAdition_;
};

#endif /* BBGRID_H */
