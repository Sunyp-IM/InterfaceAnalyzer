# This script accept a pdb file as the input and generate 5 excel files and a text file listing the atom pairs of contacts in a given
# contact interfaces. Usage:
#       python contact.py -i input -o output -s source -t target -a abbrevition
# arguments:
# -i: Path to the input pdb file, e.g., '/path/to/directory/input.pdb'
# -o: Path to the output directory , e.g., '/path/to/directory/'
# -s: Source residues, including the chain identity and residue index range. There are sever patterns for this argument:
#       (1) A   (only one chain is used as the source, no residue index range specified)
#       (2) A 20 100    (only one chain is used as the source, residue index range is specified)
#       (3) A 20 100/B 30 100 (more than one chain is used as the source, residue index range is specified)
# -t: Target residues, having the same patterns as the "-s" argument
# -a: The abbreviation for the structure studies, a string
# example:  python contact.py -i test/input.pdb -o test/ -s E/J -t H/L -a PT

import argparse
import pandas as pd
from scipy.spatial import distance
import numpy as np
from biopandas.pdb import PandasPdb


def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("-i", help="Path to the input pdb file, '/path/to/directory/input.pdb'", type=str)
    parser.add_argument("-o", help="Path to the output excel file , '/path/to/directory/'", type=str)
    parser.add_argument("-s", help="Source chain(s) and residue ranges", type=str)
    parser.add_argument("-t", help="Target chain(s) and residue ranges", type=str)
    parser.add_argument("-a", help="The abbreviation for the structure studies", type=str)

    # Parse arguments
    args = parser.parse_args()

    return args


def saltBridge(df_in):  # Detect salt bridges (cutoff = 3.5 Å)
    saltBridgeList = []
    for i in range(len(df_in)):
        # Salt bridge case1: Resn1 is an acid reside ('D' or 'E') and  Resn2 is a basic residue ('K' or 'R')
        if ((df_in['Resn1'].iloc[i][0] == 'D' and
             (df_in['Atom1'].iloc[i] == 'OD1' or df_in['Atom1'].iloc[i] == 'OD2')) or
            (df_in['Resn1'].iloc[i][0] == 'E' and
             (df_in['Atom1'].iloc[i] == 'OE1' or df_in['Atom1'].iloc[i] == 'OE2'))) \
                and \
                ((df_in['Resn2'].iloc[i][0] == 'K' and df_in['Atom2'].iloc[i] == 'NZ') or
                 (df_in['Resn2'].iloc[i][0] == 'R' and
                  (df_in['Atom2'].iloc[i] == 'NH1' or df_in['Atom2'].iloc[i] == 'NH2' or df_in['Atom2'].iloc[
                      i] == 'NE'))) \
                and \
                (df_in['Distance'].iloc[i] < 3.5):
            saltBridgeList.append(df_in.iloc[i])

        # Salt bridge case2: Resn1 is a basic residue ('K' or 'R') and  Resn2 is an acid reside ('D' or 'E')
        elif ((df_in['Resn1'].iloc[i][0] == 'K' and df_in['Atom1'].iloc[i] == 'NZ') or
              (df_in['Resn1'].iloc[i][0] == 'R' and
               (df_in['Atom1'].iloc[i] == 'NH1' or df_in['Atom1'].iloc[i] == 'NH2') or df_in['Atom2'].iloc[i] == 'NE')) \
                and \
                ((df_in['Resn2'].iloc[i][0] == 'D' and
                  (df_in['Atom2'].iloc[i] == 'OD1' or df_in['Atom2'].iloc[i] == 'OD2')) or
                 (df_in['Resn2'].iloc[i][0] == 'E' and
                  (df_in['Atom2'].iloc[i] == 'OE1' or df_in['Atom2'].iloc[i] == 'OE2'))) \
                and \
                (df_in['Distance'].iloc[i] <= 3.5):
            saltBridgeList.append(df_in.iloc[i])
    if not saltBridgeList:  # When saltBridgeList is empty
        df_out = 'None'
    else:
        df_out = pd.DataFrame(saltBridgeList)
    return df_out


def hBond(df_in):  # Detect hydrogen bonds (cutoff = 3.5 Å)
    hBondList = []
    for i in range(len(df_in)):
        if (df_in['Type1'].iloc[i] == 'O' or df_in['Type1'].iloc[i] == 'N' or df_in['Type1'].iloc[i] == 'S') \
                and (df_in['Type2'].iloc[i] == 'O' or df_in['Type2'].iloc[i] == 'N' or df_in['Type2'].iloc[i]
                     == 'S') and df_in['Distance'].iloc[i] < 3.5:
            hBondList.append(df_in.iloc[i])
    if not hBondList:  # When hBondList is empty
        df_out = 'None'
    else:
        df_out = pd.DataFrame(hBondList)
        df_out.index = range(len(df_out))
    return df_out


def nonPolar(df_in):  # Detect non-polar interactions
    nonPolarList = []
    for i in range(len(df_in)):
        if df_in['Type1'].iloc[i] == 'C' and df_in['Type2'].iloc[i] == 'C':
            nonPolarList.append(df_in.iloc[i])
    if not nonPolarList:  # When saltBridgeList is not empty
        df_out = 'None'
    else:
        df_out = pd.DataFrame(nonPolarList)
    return df_out


AA_dic = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B', 'CYS': 'C', 'GLU': 'E',
          'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
          'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
          'TYR': 'Y', 'VAL': 'V'}

if __name__ == '__main__':
    # Generate "contact.xlsx", the list of contact atom pairs.
    args = parseArguments()
    infile = args.i
    outDir = args.o
    sr = args.s  # "source residues"
    tr = args.t  # "target residues"
    abbr = args.a
    
    pPdb = PandasPdb().read_pdb(infile)
    df = pd.concat([pPdb.df['ATOM'], pPdb.df['HETATM']])
    # Define the source residues
    if '/' in sr:  # When the source contains more than one chain
        srList = sr.split('/')
        df1 = pd.DataFrame()
        for i in range(len(srList)):
            if len(srList[i]) == 1:  # When the source chain do not contain residue index range
                df10 = df[df['chain_id'] == srList[i]]
            else:
                baseList = srList[i].split(' ')
                df10 = df[(df['chain_id'] == baseList[0]) & (df['residue_number'] > baseList[1]) & (
                        df['residue_number'] < baseList[2])]
            df1 = df1.append(df10)
    else:  # When the source contains only one chain
        srList = [sr]
        if len(sr) == 1:  # When the source chain do not contain residue index range
            df1 = df[df['chain_id'] == sr]
        else:  # When the source contains residue index range
            baseList = sr.split(' ')
            df1 = df[(df['chain_id'] == baseList[0]) & (df['residue_number'] > baseList[1]) & (
                    df['residue_number'] < baseList[2])]
    # Define the target residues
    if '/' in tr:  # When the target contains more than one chain
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
        trList = [tr]
        if len(tr) == 1:  # When the target chain do not contain residue index range
            df2 = df[df['chain_id'] == tr]
        else:  # When the target chain contains residue index range
            baseList = sr.split(' ')
            df2 = df[(df['chain_id'] == baseList[0]) & (df['residue_number'] > baseList[1]) & (
                    df['residue_number'] < baseList[2])]

    contactList = []
    A1 = df1.iloc[:, 11:14].to_numpy()
    A2 = df2.iloc[:, 11:14].to_numpy()
    A12 = distance.cdist(A1, A2, 'euclidean')  # An array for distance between the two chains)
    T_contact = np.where(A12 < 4.5)  # A tuple for indices of the contact atoms)
    for i in range(len(T_contact[0])):
        ind1 = T_contact[0][i]
        ind2 = T_contact[1][i]
        source = df1.iloc[ind1, 7]
        if df1.iloc[ind1, 5] in AA_dic.keys():
            resn1 = AA_dic[df1.iloc[ind1, 5]] + str(df1.iloc[ind1, 8])
        else:
            resn1 = df1.iloc[ind1, 5] + str(df1.iloc[ind1, 8])
        resi1 = df1.iloc[ind1, 8]
        atom1 = df1.iloc[ind1, 3]
        type1 = df1.iloc[ind1, 18]
        target = df2.iloc[ind2, 7]
        if df2.iloc[ind2, 5] in AA_dic.keys():
            resn2 = AA_dic[df2.iloc[ind2, 5]] + str(df2.iloc[ind2, 8])
        else:
            resn2 = df2.iloc[ind2, 5] + str(df2.iloc[ind2, 8])
        resi2 = df2.iloc[ind2, 8]
        atom2 = df2.iloc[ind2, 3]
        type2 = df2.iloc[ind2, 18]
        dist = A12[ind1, ind2]
        contactList.append([source, resn1, resi1, atom1, type1, target, resn2, resi2, atom2, type2, dist])
    col = ['Source', 'Resn1', 'Resi1', 'Atom1', 'Type1', 'Target', 'Resn2', 'Resi2', 'Atom2', 'Type2',
           'Distance']
    df_contact = pd.DataFrame(contactList, columns=col)
    ind = range(len(df_contact))
    df_contact.insert(loc=0, column='Index', value=ind)
    df_contact.to_excel(outDir + abbr + '_contact.xlsx')   
    with open(outDir + abbr + '_contact_residues', 'w') as f:
        for i in srList:
            resi1List = df_contact[df_contact.Source == i].Resi1.drop_duplicates().tolist()
            f.write('Chain' + i + ' (Source): ' + str(resi1List) + '\n')
        for i in trList:
            resi2List = df_contact[df_contact.Target == i].Resi2.drop_duplicates().tolist()
            f.write('Chain' + i + ' (Target): ' + str(resi2List) + '\n') 

    # Generate lists for hydrogen bonds, salt bridges, non-polar contacts and vdw contacts
    df_saltBridge = saltBridge(df_contact)
    if isinstance(df_saltBridge, pd.DataFrame):
        # df_saltBridge.columns = col
        df_saltBridge.to_excel(outDir + abbr + '_saltBridge.xlsx')
        # Remove salt bridges
        df_noSaltBridge = df_contact.drop(df_contact[df_contact.Index.isin(df_saltBridge.Index)].index)
    else:
        print('No salt bridges found')
        df_noSaltBridge = df_contact

    df_hBond = hBond(df_noSaltBridge)
    if isinstance(df_hBond, pd.DataFrame):
        df_hBond.to_excel(outDir + abbr + '_hBond.xlsx')
        # Remove hydrogen bonds
        df_nohBond = df_noSaltBridge.drop(df_noSaltBridge[df_noSaltBridge.Index.isin(df_hBond.Index)].index)
        df_nohBond.to_excel(outDir + abbr + '_vdw.xlsx')
    else:
        print('No hydrogen bond found')

    df_nonPolar = nonPolar(df_nohBond)
    if isinstance(df_nonPolar, pd.DataFrame):
        df_nonPolar.to_excel(outDir + abbr + '_nonPolar.xlsx')
    else:
        print('No non-polar contact found')

    # Generate a summary list for contacts
    group_source = df_contact.groupby(['Source', 'Resn1'], sort=False)
    group_source_list = list(group_source)
    summaryList = []
    for i in range(len(group_source_list)):
        source = group_source_list[i][0][0]
        Resn1 = group_source_list[i][0][1]
        df_source = group_source_list[i][1]
        group_target = df_source.groupby(['Target', 'Resn2'])
        group_target_list = list(group_target)
        res2nList = []  # All contact contact residues in the target chain for a single source chain reside.
        for j in range(len(group_target_list)):
            target = group_target_list[j][0][0]
            resn2 = group_target_list[j][0][1]
            df_target = group_target_list[j][1]
            n_contact = len(df_target)
            if isinstance(df_hBond, pd.DataFrame):
                n_hBond = len(pd.merge(df_hBond, df_target, how='inner', on=['Index']))
            else:
                n_hBond = 0
            if isinstance(df_saltBridge, pd.DataFrame):
                n_saltBridge = len(pd.merge(df_saltBridge, df_target, how='inner', on=['Index']))
            else:
                n_saltBridge = 0
            n_vdw = n_contact - n_hBond - n_saltBridge
            if n_hBond == 0 and n_saltBridge == 0:
                s = target + '.' + resn2 + ' (' + str(n_vdw) + ')'
            elif n_hBond > 0 and n_saltBridge == 0:
                s = target + '.' + resn2 + ' (' + str(n_vdw) + ', ' + str(n_hBond) + ')'
            else:
                s = target + '.' + resn2 + ' (' + str(n_vdw) + ', ' + str(n_hBond) + ', ' + str(n_saltBridge) + ')'
            res2nList.append(s)
        summaryList.append([source + '.' + Resn1, res2nList])
    df_summary = pd.DataFrame(summaryList, columns=['Source', abbr])
    df_summary.to_excel(outDir + abbr + '_summary.xlsx')

