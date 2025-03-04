# This script is used for draw contacts listed in the input excel file
# Usage: 
#   (1) Open pymol; 
#   (2) input "run drawContact.py" in the command line'; 
#   (3) plot_dist(infile, mol, abbr, bondType). Example: plot_dist('/path/to/input.xlsx, '3G10', 3g10A', 'HB')
# Parameters:
#   infile: an excel file generated by contact.py,describing all the contact, hydrogen bonds or saltbridges
#   mol: an object in pymol, in which the contact will be drawn.
#   abbr: an abbreviation for the contacts you want to draw.
#   The contact type you want to draw. You can use "HB" for hydrogen bonds, and 'SB' for salt bridges
    


import os
import pandas as pd
import re


def plot_dist(infile, mol, abbr, bondType):
    df = pd.read_excel(infile)
    df.index = range(len(df))
    for i in range(len(df)):
        cmd.select("resi1", mol + " and chain " + df.Source[i] + " and resi " + str(df.Resi1[i]))
        cmd.select("resi2", mol + " and chain " + df.Target[i] + " and resi " + str(df.Resi2[i])) 
        cmd.show("sticks", "resi1")
        cmd.show("sticks", "resi2")
        cmd.select("atom1", 'resi1' + " and name " + df.Atom1[i])
        cmd.select("atom2", 'resi2' + " and name " + df.Atom2[i])
        cmd.distance(abbr + '_' + bondType, "atom1", "atom2")
        cmd.select("resi1_" + abbr + '_' + bondType, "resi1", merge=1)
        cmd.select("resi2_" + abbr + '_' + bondType, "resi2", merge=1)
    cmd.hide('label')
    cmd.delete('resi1')
    cmd.delete('resi2')
    cmd.delete('atom1')
    cmd.delete('atom2')
    cmd.set('cartoon_transparency', 0.8)
    cmd.set('stick_radius', 0.15)

one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}



#def label_resn(mol):
#    cmd.select(mol + '_CA', mol + '& n. CA')
#    cmd.label(mol + '_CA', "%s%s" % (one_letter[resn],resi))
#cmd.set('label_size', 20)
#cmd.set('label_font_id', 10)

