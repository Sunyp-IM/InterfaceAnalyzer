# This script is used for draw contacts listed in the input excel file

import os
import pandas as pd
import re


def plot_dist(infile, mol, abbr, bondType):
    df = pd.read_excel(infile)
    df.index = range(len(df))
    for i in range(len(df)):
        cmd.select("res1", mol + " and chain " + df.Source[i] + " and resi " + re.findall(r"\d+", df.Res1[i])[0])
        cmd.select("res2", mol + " and chain " + df.Target[i] + " and resi " + re.findall(r"\d+", df.Res2[i])[0])
        cmd.show("sticks", "res1")
        cmd.show("sticks", "res2")
        cmd.select("atom1", 'res1' + " and name " + df.Atom1[i])
        cmd.select("atom2", 'res2' + " and name " + df.Atom2[i])
        cmd.distance(abbr + bondType, "atom1", "atom2")
        cmd.select("res1_" + abbr + bondType, "res1", merge=1)
        cmd.select("res2_" + abbr + bondType, "res2", merge=1)
    cmd.hide('label')
    cmd.delete('res1', 'res2')
    cmd.delete('atom1', 'atom2')


#:
# rootPath = os.path.abspath(os.path.join(os.getcwd(), ".."))
# infile1 = rootPath + '/S309/PT-S309_hbond.xlsx'
# infile2 = rootPath + '/S309/BA.2-S309_hbond.xlsx'
# mol1 = 'PT-S309_7bep'
# abbr1 = "PT"
# mol2 = 'BA.2-S309'
# abbr2 = "BA.2"
#
# mol_dic = {'mol1': [infile1, mol1, abbr1], 'mol2': [infile2, mol2, abbr2]}
# # for key in mol_dic.keys():
# plot_dist(mol_dic['mol2'][0], mol_dic['mol2'][1], mol_dic['mol2'][2])
