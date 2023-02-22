#ifndef ConnollySurface_H
#define ConnollySurface_H


#include "ChemMolecule.h"
#include "Surface.h"

/** The algorithm is taken from Connolly's original MS program, which is
    freely distributable and Copyright 1983, Michael Connolly.

    M.L. Connolly, "Solvent-accessible surfaces of proteins and nucleic acids",
    Science, 221, p709-713 (1983).

    M.L. Connolly, "Analytical molecular surface calculation",
    J. Appl. Cryst. 16, p548-558 (1983).
 */
Surface get_connolly_surface(const ChemMolecule& molecule,
                             float density, float probe_radius);


#endif
