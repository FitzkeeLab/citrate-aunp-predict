#
# title:   binding_surface.py
# summary: Calculate binding surface of a PDB with alpha in the bfac column
# date:    April 18, 2023
# author:  Nick Fitzkee (nfitzkee at chemistry.msstate.edu)
#
# history:
# - 20230418 - Added additional chain IDs so that large number of atoms
#              can be displayed without shifting columns
#            - Prints a warning if too many virtual atoms are displayed
#

import sys, math, os
from Bio.PDB import *
from numpy import *

def find_extents(m, extra=5):
    a = m.child_list[0].child_list[0].child_list[0]
    x, y, z = a.get_coord()
    
    xmin, xmax = x, x
    ymin, ymax = y, y
    zmin, zmax = z, z

    for a in m.get_atoms():
        x, y, z = a.get_coord()

        if x > xmax: xmax = x
        if y > ymax: ymax = y
        if z > zmax: zmax = z

        if x < xmin: xmin = x
        if y < ymin: ymin = y
        if z < zmin: zmin = z


    #xmin = int((1.0-math.copysign(extra, xmin))*xmin)
    #ymin = int((1.0-math.copysign(extra, ymin))*ymin)
    #zmin = int((1.0-math.copysign(extra, zmin))*zmin)

    #xmax = int((1.0+math.copysign(extra, xmax))*xmax)
    #ymax = int((1.0+math.copysign(extra, ymax))*ymax)
    #zmax = int((1.0+math.copysign(extra, zmax))*zmax)

    xmin = int(xmin - extra)
    ymin = int(ymin - extra)
    zmin = int(zmin - extra)

    xmax = int(xmax + extra)
    ymax = int(ymax + extra)
    zmax = int(zmax + extra)
    
    #xmax = int((1.0+extra)*xmax)
    #ymax = int((1.0+extra)*ymax)
    #zmax = int((1.0+extra)*zmax)
    
    return xmin, xmax, ymin, ymax, zmin, zmax

def generate_mesh(pdb, out):
    """Generate mesh around a PDB file"""

    s = PDBParser().get_structure('target', pdb)
    m = s[0]

    xmin, xmax, ymin, ymax, zmin, zmax = find_extents(m)

    all_atoms = list(m.get_atoms())
    num_atoms = len(all_atoms)
    
    ns = NeighborSearch(all_atoms)

    dgrid = 1.0

    ngrx = int((xmax - xmin)/dgrid)
    ngry = int((ymax - ymin)/dgrid)
    ngrz = int((zmax - zmin)/dgrid)

    print("x-axis: min=%8.3f max=%8.3f pts=%i" % (xmin, xmax, ngrx))
    print("y-axis: min=%8.3f max=%8.3f pts=%i" % (ymin, ymax, ngry))
    print("z-axis: min=%8.3f max=%8.3f pts=%i" % (zmin, zmax, ngrz))
    
    max_bfac = None
    max_coord = None

    with open(out, 'w') as f:
        w = f.write
        fmt = 'HETATM%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f '\
              '          %1s\n'

        cnams = 'WXYZwxyz'
        aidx = 1
        ridx = 1
        cidx = 0
        cnam = cnams[cidx]
        
        for xi in range(ngrx):
            x = xmin + xi*dgrid
            for yi in range(ngry):
                y = ymin + yi*dgrid
                for zi in range(ngrz):
                    z = zmin + zi*dgrid
                    my_center = array((x, y, z))

                    atoms_tmp = ns.search(my_center, radius=2.0)
                    if len(atoms_tmp) > 0:
                        continue                

                    atoms_tmp = ns.search(my_center, radius=3.0)
                    if len(atoms_tmp) == 0:
                        continue
                    
                    #w(fmt % (aidx, 'O'.center(4), ' ', 'HOH', cnam, ridx,
                    #         ' ', x, y, z, 1.0, 1.0, 'O'))

                    #print(x,y,z)

                    avg_bfac = 0
                    atoms = list(ns.search(array((x, y, z)), radius=10.0))

                    for a in atoms:
                        avg_bfac = avg_bfac + max(a.get_bfactor(), 0.0)

                    if len(atoms):
                        avg_bfac = avg_bfac/len(atoms)

                    if avg_bfac > 30.0:
                        w(fmt % (aidx, 'O'.center(4), ' ', 'HOH', cnam,
                                 ridx, ' ', x, y, z, 1.0, avg_bfac, 'O'))
                        aidx = aidx + 1
                        ridx = ridx + 1

                        if ridx > 9999:
                            ridx = 1
                            cidx = cidx + 1

                            if cidx >= len(cnams):
                                print("Warning: Resetting Chain Counter!")
                                cidx = 0
                                
                            cnam = cnams[cidx]

                    if max_bfac is None or avg_bfac > max_bfac:
                        max_bfac = avg_bfac
                        max_coord = x, y, z

        x, y, z = max_coord
        w(fmt % (aidx, 'OXT'.center(4), ' ', 'HOH', cnam,
                 ridx, ' ', x, y, z, 1.0, avg_bfac, 'O'))

        print(max_bfac, max_coord)
                

                
if __name__ == '__main__':
    try:
        pdb = sys.argv[1]
        out = sys.argv[2]
    except:
        print('usage: %s <pdb> <new_pdb>' % os.path.split(sys.argv[0])[1])
        sys.exit(1)

        
    generate_mesh(pdb, out)
