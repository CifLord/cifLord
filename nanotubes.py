__author__ = 'richard'

import numpy as np

from pymatgen.transformations.standard_transformations import RotationTransformation
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.operations import SymmOp

import math
import numpy as np

from pymatgen import write_structure
from pymatgen.io.smartio import CifParser

def nanotubes(initial_structure, length, radius):
    new_sites = []
    new_species = []
    lattice = []

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
    coordinates = []

    print(type(super_cell.lattice.matrix[2]))
    for site in super_cell:
        # Distance of a site to the center of the nt
        site_r = outer_r - site.coords[2]
        # print(site.coords[1], site.coords[2])
        # In radians
        ang = site.frac_coords[1]*2*math.pi

        new_sites.append([site.coords[0],
                          site_r*math.cos(ang),
                          site_r*math.sin(ang)])

        new_species.append(site.specie)
        min_b.append(site_r*math.cos(ang))
        min_c.append(site_r*math.sin(ang))

    add_c = -1*min(min_c)
    add_b = -1*min(min_b)

    for i, site in enumerate(new_sites):
        new_sites[i][1] = site[1] + add_b
        new_sites[i][2] = site[2] + add_c

    latt = Lattice([super_cell.lattice.matrix[0],
                   [0, outer_r*2, 0],
                   [0, 0, outer_r*2]])

    print(latt)
    print(super_cell.lattice)

    return Structure(latt, new_species, new_sites, coords_are_cartesian=True)
