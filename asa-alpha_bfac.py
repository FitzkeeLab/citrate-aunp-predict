#!/usr/bin/python3

#
# title:   asa-alpha_bfac.py
# summary: Alter B-factor column of PDF to the product of RASA and ALPHA
# author:  Nicholas Fitzkee (nfitzkee at chemistry.msstate.edu)
# date:    February 19, 2021
#

from Bio import PDB, SeqUtils
from Bio.PDB import NACCESS as ASA
import os, sys


# Relative side-chain atom residue relative ASA cutoff (in percent)
REL_ASA_MIN = 25.0


# Table of Joanna's Alpha values
ALPHA = {'GLY': 0.980,
         'ALA': 0.660,
         'LEU': 0.400,
         'ILE': 0.402,
         'VAL': 0.455,
         'MET': 0.743,
         'PRO': 0.625,
         'PHE': 0.461,
         'TYR': 0.474,
         'TRP': 0.412,
         'SER': 0.570,
         'CYS': 5.692,
         'THR': 0.502,
         'GLU': 0.352,
         'ASP': 0.466,
         'ASN': 0.541,
         'GLN': 0.401,
         'HIS': 0.791,
         'LYS': 0.812,
         'ARG': 0.837
         }


def residue_surface_area(pdb, out=sys.stdout, chain=None):
    s = PDB.PDBParser().get_structure('target', pdb)
    m = s[0]             # First "model" 

    if not chain is None:
        if chain == '_': chain = ' '
        
        if not chain in m:
            print ('Chain "%s" not found in PDB file.' % chain)
            sys.exit(1)
        
        for c in m.child_list:
            if c.id != chain:
                del m[c.id]

    c = m.child_list[0]  # First "chain"

    # Run NACCESS (dependency) and obtain relative side chain
    # accessible surface area.
    rsa_data, asa_data = ASA.run_naccess(s, pdb)
    rsa_dict = ASA.process_rsa_data(rsa_data)

    avg_prod = 0.0
    max_prod = -1.0
    avg_prdn = 0
    
    # Iterate through all residues in chain
    for res in c.child_list:
        het, id, icode = res.id
        nam = res.resname
        asa_code = (c.id, res.id)

        asa = rsa_dict[asa_code]['side_chain_rel']

        new_bfac = -1.0
        if asa > REL_ASA_MIN:
            new_bfac = asa*ALPHA[nam]
            avg_prod = avg_prod + new_bfac
            avg_prdn = avg_prdn + 1

            if new_bfac > max_prod:
                max_prod = new_bfac
            
        for atom in res.child_list:
            atom.set_bfactor(new_bfac)

    avg_prod = avg_prod/avg_prdn
    
    # Save modified PDB structure; give the option to write it
    # to standard output (screen) or a file
    close_out = False
    if type(out) == type('string'):
        out = open(out, 'w')
        close_out = True

    out.write('REMARK Average bfactor for accessible residues is %8.3f\n' %
              avg_prod)
    out.write('REMARK Maximum bfactor is %8.3f\nREMARK\n' % max_prod)
    out.write('REMARK PyMOL First approx: spectrum b, yellow_red, '
              'minimum=%.3f, maximum=%.3f\nREMARK\n' % (avg_prod, max_prod))

    io = PDB.PDBIO()
    io.set_structure(s)
    io.save(out)

    if close_out: out.close()

if __name__ == '__main__':
    try:
        pdb = sys.argv[1]
        out = sys.argv[2]
    except:
        print('usage: %s <pdb> <new_pdb> [<chain>]' % os.path.split(sys.argv[0])[1])
        print('\n'
              '  <chain>    is the (optional) one character chain ID;\n'
              '             if none specified, the first chain is used.')
        
        sys.exit(1)

    try:
        chain = sys.argv[3]
    except:
        chain = None
    
    residue_surface_area(pdb, out, chain)
