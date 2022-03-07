# citrate-aunp-predict
Tools to predict favorable protein binding regions with a citrate-capped gold nanoparticle

This program consists of three scripts:

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
