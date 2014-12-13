__author__ = 'richard'

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite

from bucking_pbc import fix_pbc

import math
import numpy as np

def cart_to_polar(site, outer_r):
    site_r = outer_r - site.coords[2]
    ang = site.frac_coords[1]*2*math.pi
    return [site.coords[0], site_r*math.cos(ang), site_r*math.sin(ang)]


def nanotubes(initial_structure, length, radius):
    new_sites = []
    new_species = []

    # b will be the circumference of the nt, a will be
    # its length and c will be the thickness of the nt.
    circumference = 2*math.pi*radius
    n_units = math.ceil(circumference/initial_structure.lattice.b)
    print(n_units)
    super_cell = initial_structure.copy()
    super_cell.make_supercell([length, n_units, 1])
    outer_circumference = super_cell.lattice.b
    outer_r = outer_circumference/(2*math.pi)

    # A sin wave has x (theta) and y components
    # Use asin with domain between 0 and 1,
    # convenient for fractional coordinates.
    min_b = []
    min_c = []

    super_cell = fix_pbc(super_cell)


    for site in super_cell:
        if np.round(site.frac_coords[2], decimals = 5) == 0:
            rep_site = PeriodicSite(site._species,
                                    [site.frac_coords[0],
                                     site.frac_coords[1],
                                     1],
                                    super_cell.lattice)
            print(rep_site.frac_coords)
            # rep_site.frac_coords[2] = 1
            print(rep_site.frac_coords)
            new_site_coord = cart_to_polar(rep_site, outer_r)
            new_sites.append(new_site_coord)
            new_species.append(site.specie)

        # Distance of a site to the center of the nt
        new_site_coord = cart_to_polar(site, outer_r)
        new_sites.append(new_site_coord)

        new_species.append(site.specie)
        min_b.append(new_site_coord[1])
        min_c.append(new_site_coord[2])

    add_c = -1*min(min_c)
    add_b = -1*min(min_b)

    for i, site in enumerate(new_sites):
        new_sites[i][1] = site[1] + add_b
        new_sites[i][2] = site[2] + add_c

    latt = Lattice([super_cell.lattice.matrix[0],
                   [0, outer_r*2, 0],
                   [0, 0, outer_r*2]])

    return Structure(latt, new_species, new_sites, coords_are_cartesian=True)
