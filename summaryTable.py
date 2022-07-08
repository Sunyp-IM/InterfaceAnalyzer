import os
import sys
import argparse
import pandas as pd

path = '/data/github/InterfaceAnalyzer/test/'
df_contact = pd.read_excel(path + 'PT_contact.xlsx', index_col=[0])
df_hBond = pd.read_excel(path + 'PT_hBond.xlsx',  index_col=[0])
df_saltBridge = pd.read_excel(path + 'PT_saltBridge.xlsx', index_col=[0])
abbr = 'PT'
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
        n_hBond = len(pd.merge(df_hBond, df_target, how='inner', on=['Index']))
        n_saltBridge = len(pd.merge(df_saltBridge, df_target, how='inner', on=['Index']))
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
df_summary.to_excel(path + abbr + '_summary.xlsx')
