# citrate-aunp-predict
Tools to predict favorable protein binding regions with a citrate-capped gold nanoparticle

This project consists of three Python 3 scripts:

1. asa-alpha_bfac.py     : A program the populates a PDB B-factor column with alpha*(relative ASA)
2. asa-alpha_all_bfac.py : Similar to above, but operates on all chain IDs in the PDB
3. binding_surface.py    : A program that converts the B-factor column to pseudoatoms at the surface

It requires the following Python libraries:

1. BioPython (Bio)
2. Numeric Python (numpy)
3. Scientific Python (scipy)

These dependencies can be easily installed using the Package Installer for Pythong (pip):

```
pip install --user numpy
pip install --user scipy
pip install --user biopython
```

It also requires NACCESS, which is free for academic users, but less easy to obtain. Instructions
can be found at [the NACCESS website](http://www.bioinf.manchester.ac.uk/naccess/). NACCESSS is 
needed to calculate relative sidechain accessible surface areas (RASA).

## ASA-Alpha Bfac Program

This program takes a PDB file input and writes a new PDB file with a modified B-factor column. The 
B-factor becomes the product of the relative side chain accessible surface area (RASA, 0-100%) with 
the nanoparticle binding alpha value. Higher binding residues are expected to interact with more 
favorably with the nanoparticle surface. Two versions of this program are provided: one where the
chain ID is selectable as a command line argument (asa-alpha_bfac.py) and one where all chain IDs
are selected (asa-alpha_all_bfac.py). Residues with less than 25% RASA are assigned a b-factor of -1.

Both programs leave non-protein atoms unchanged, but if one wishes to remove water atoms and
heteroatoms (which can affect RASA), one can remove these atoms in PyMOL (or manually) first.

The output PDB file contains a suggested command to color the residues by b-factor in PyMOL, 
assuming the complete range is desired (e.g. spectrum_b ...).

## Binding Surface Program

This program takes a PDB file input and writes an output PDB file representing the surface mesh
of pseudoatoms. A 1 angstrom by 1 angstrom mesh is overlaid on the protein. Only atoms on the 
surface are drawn - surface pseudoatoms are selected based on steric clash and their distance 
from actual protein atoms.

Surface pseudoatoms are assigned their own b-factor based on the average alpha * RASA value of 
nearby residues. Surface pseudoatoms that are nearby residues that promote surface binding
have a corresponding higher bfactor themselves. Pseudoatoms below a particular cutoff 
(currently 30) are not drawn. The output is written as a separate PDB file.

In addition to writing the output PDB file, the program writes the highest pseudoatom
b-factor, along with its x, y, z coordinates.

## Usage Example

A sample PDB file is included with this distribution in the example directory (GB3, PDB 2OED). 
If the dependencies are installed, the surface calculations can be performed using the 
following commands:

```
python3 asa-alpha_all_bfac.py 2oed.pdb 2oed_np.pdb
python3 binding-surface.py 2oed_np.pdb 2oed_surf.pdb
```

Example versions of the output files are also included. 

To view the file, we recommend PyMOL. To generate figures such as those shown in the text,
load both the 2oed.pdb and the 2oed_surf.pdb, then type the following commands:

```
hide
show cartoon, 2oed
color grey, 2oed
show spheres, 2oed_surf
spectrum b, white_red, minimum=30, maximum=60, selection=2oed_surf
set sphere_scale, 0.5, 2oed_surf
```

