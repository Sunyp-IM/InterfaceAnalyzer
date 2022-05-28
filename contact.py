# This script accept a pdb file as the input and generate a excel file listing the atom pairs of contacts in a given
# contact interfaces. Usage:
#       python contact.py -i input -o output -s source -t target
# options:
# -i: Full path to the input pdb file, e.g., '/path/to/directory/????.pdb'
# -o: Full path to the output excel file , e.g., '/path/to/directory/???.xlsx'
# -s: Source residues, including the chain identity and residue index range. There are sever patterns for this option:
#       (1) A   (only one chain is used as the source, no residue index range specified)
#       (2) A 20 100    (only one chain is used as the source, residue index range is specified)
#       (3) A 20 100/B 30 100 (more than one chain is used as the source, residue index range is specified)
# -t: Target residues, having the same patterns as the "-s" option


import os
import sys
import argparse
import pandas as pd
import math
import re
from scipy.spatial import distance
import numpy as np
from biopandas.pdb import PandasPdb


def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("-i", help="Full path to the input pdb file, '/path/to/directory/????.pdb'", type=str)
    parser.add_argument("-o", help="Full path to the output excel file , '/path/to/directory/???.xlsx'", type=str)
    parser.add_argument("-s", help="Source chain(s) and residue ranges", type=str)
    parser.add_argument("-t", help="Target chain(s) and residue ranges", type=str)

    # Parse arguments
    args = parser.parse_args()

    return args


AA_dic = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B', 'CYS': 'C', 'GLU': 'E',
          'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
          'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
          'TYR': 'Y', 'VAL': 'V'}

if __name__ == '__main__':
    opt = parseArguments()
    infile = opt.i
    outfile = opt.o
    sr = opt.s  # "source residues"
    tr = opt.t  # "target residues"
    ppdb = PandasPdb().read_pdb(infile)
    df = pd.concat([ppdb.df['ATOM'], ppdb.df['HETATM']])
    # Define the source residues
    if '/' in opt.s:  # When the source contains more than one chain
        srList = sr.split('/')
        df1 = pd.DataFrame()
        for i in range(len(srList)):
            if len(srList[i]) == 1:  # When the source chain do not contain residue index range
                df10 = df[df['chain_id'] == srList[i]]
            else:
                baseList = srList[i].split(' ')
                df10 = df[(df['chain_id'] == baseList[0]) & (df['residue_number'] > baseList[1]) & (
                            df['residue_number'] < baseList[2])]
            print('df10 = ', df10)
            df1 = df1.append(df10)
            print('df1 = ', df1)
    else:  # When the source contains only one chain
        if len(opt.s) == 1:  # When the source chain do not contain residue index range
            df1 = df[df['chain_id'] == opt.s]
        else:  # When the source contains residue index range
            baseList = opt.s.split(' ')
            df1 = df[(df['chain_id'] == baseList[0]) & (df['residue_number'] > baseList[1]) & (
                        df['residue_number'] < baseList[2])]
    # Define the target residues
    if '/' in opt.t:  # When the target contains more than one chain
        trList = tr.split('/')
        df2 = pd.DataFrame()
        for i in range(len(trList)):
            if len(trList[i]) == 1:  # When the target chain do not contain residue index range
                df20 = df[df['chain_id'] == trList[i]]
            else:  # When the target chain contains residue index range
                baseList = trList[i].split(' ')
                df20 = df[(df['chain_id'] == baseList[0]) & (df['residue_number'] > baseList[1]) & (
                        df['residue_number'] < baseList[2])]
            df2 = df2.append(df20)
    else:  # When the target contains only one chain
        if len(opt.t) == 1:  # When the target chain do not contain residue index range
            df2 = df[df['chain_id'] == opt.t]
        else:  # When the target chain contains residue index range
            baseList = opt.s.split(' ')
            df2 = df[(df['chain_id'] == baseList[0]) & (df['residue_number'] > baseList[1]) & (
                    df['residue_number'] < baseList[2])]

    # Run from Pycharm (one input file)
    # rootPath = os.path.abspath(os.path.join(os.getwdu(), ".."))
    # inPath = rootPath + '/originalStruct/'
    # outPath = rootPath + '/contacts/S309/'
    # infile = inPath + 'BA.1_Omicron-Striker_final.pdb'
    # infile = inPath + 'BA.2-Spike_20220506_zzn_final.pdb'
    # infile = inPath + 'BA2-S309_Local_20220514.pdb'
    # infile = inPath + '7bep_PT-S309.pdb'
    # infile = inPath + 'OmicronRBD-S309_Local_20220201.pdb'
    # f = open(infile, 'r', encoding='UTF-8')
    # outfile = open(outPath + 'simplePDB.pdb', 'w', encoding='UTF-8')
    # try:
    #     for line in f:
    #         g = re.search("^ATOM", line)
    #         if g:
    #             outfile.writelines(line)
    # finally:
    #     outfile.close()
    # df = pd.read_csv(outPath + 'simplePDB.pdb', header=None, sep='\s+')
    # for i in range(len(df)):
    #     if len(df.iloc[i, 9]) > 4:
    #         df.iloc[i, 11] = df.iloc[i, 10]
    #         df.iloc[i, 10] = float(str(df.iloc[i, 9])[3:])
    #         df.iloc[i, 9] = float(str(df.iloc[i, 9])[:3])
    # df.to_excel(outPath + 'simplePDB.xlsx')
    #
    # # Test whether df contains nan
    # if df.isnull().values.any():
    #     print('Input file format error')
    #     sys.exit(1)

    # Calculate Distance between two points

    # df1 = df[(df['chain_id'] == 'C') & (df['residue_number'] > 330) & (df['residue_number'] < 526)]
    # df2 = df[(df['chain_id'] == 'C') & (df['residue_number'] > 332) & (df['residue_number'] < 526)]
    # df1 = df[(df['chain_id'] == 'A') & (df['residue_number'] > 15) & (df['residue_number'] < 306)]
    # df2 = df[(df['chain_id'] == 'C') & (df['residue_number'] > 15) & (df['residue_number'] < 306)]
    # df2 = df[~((df['chain_id'] == 'C') & (df['residue_number'] > 330) & (df['residue_number'] < 526))]

    # Contacts in "BA2-S309_Local_20220514.pdb"
    # df1 = df[(df['chain_id'] == 'B')]
    # df2 = df[(df['chain_id'] == 'A') | (df['chain_id'] == 'C')]

    # Contacts in "BA2-S309_Local_20220514.pdb"
    # df1 = df[(df['chain_id'] == 'B')]
    # df2 = df[(df['chain_id'] == 'A') | (df['chain_id'] == 'C')]

    df1.to_excel('/data/omicron-S/pycharm/test/df1.xlsx')
    df2.to_excel('/data/omicron-S/pycharm/test/df2.xlsx')

    contactList = []
    A1 = df1.iloc[:, 11:14].to_numpy()
    A2 = df2.iloc[:, 11:14].to_numpy()
    A12 = distance.cdist(A1, A2, 'euclidean')  # An array for distance between the two chains)
    T_contact = np.where(A12 <= 4.5)  # A tuple for indices of the contact atoms)
    print("A1 = ", A1)
    print("A2 = ", A2)
    print('T_contact = ', T_contact)
    for i in range(len(T_contact[0])):
        ind1 = T_contact[0][i]
        ind2 = T_contact[1][i]
        source = df1.iloc[ind1, 7]
        if df1.iloc[ind1, 5] in AA_dic.keys():
            res1 = AA_dic[df1.iloc[ind1, 5]] + str(df1.iloc[ind1, 8])
        else:
            res1 = df1.iloc[ind1, 5] + str(df1.iloc[ind1, 8])
        atom1 = df1.iloc[ind1, 3]
        type1 = df1.iloc[ind1, 18]
        target = df2.iloc[ind2, 7]
        if df2.iloc[ind2, 5] in AA_dic.keys():
            res2 = AA_dic[df2.iloc[ind2, 5]] + str(df2.iloc[ind2, 8])
        else:
            res2 = df2.iloc[ind2, 5] + str(df2.iloc[ind2, 8])
        atom2 = df2.iloc[ind2, 3]
        type2 = df2.iloc[ind2, 18]
        dist = A12[ind1, ind2]
        contactList.append([source, res1, atom1, type1, target, res2, atom2, type2, dist])
    col = ['Source', 'Res1', 'Atom1', 'Type1', 'Target', 'Res2', 'Atom2', 'Type2', 'Distance']
    df_contact = pd.DataFrame(contactList, columns=col)
    df_contact.to_excel(opt.o)
