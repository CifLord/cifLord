__author__ = 'richard'

from pymatgen.core.structure import Structure
import numpy as np

# AAAARRRRRGGGGGGHHHHH BUCKING PBC!!!!

def fix_pbc(structure):
    # All sites with c = 1 will be shifted to c = 0 to
    # fix PBC problem and coordinates rounded to five
    # decimal places to avoid problems with boundary

    site_fix = []
    species_fix = []
    for site in structure:
        nsite = site.frac_coords
        nsite[2] = np.around(nsite[2], decimals=5)
        if nsite[2] == 1:
            nsite[2] = 0
        site_fix.append(nsite)
        species_fix.append(site.specie)

    return Structure(structure.lattice.matrix,
                              species_fix, site_fix)