# -*- coding: utf-8 -
# 将PDB文件内的氨基酸相关信息读取到dataframe中

import os

import pandas as pd
from pandas import DataFrame

from data import donor, receptor
from epitope_judgment import ab_ag_epi_df
from Wheels import chain_antibody_antigen, dict_to_df, differ_ab_ag, distance,open_pdb_file


# 获取当前.py文件所在的目录
def get_py_path(name_py):
    py_path = os.path.dirname(os.path.realpath(name_py))
    return py_path

"""
将可以形成疏水作用力的氨基酸的坐标的几何中心进行计算,并保存到datafreme中

便于下一步的计算
"""


def Hydrophobic_Geometric_center_coordinates(pdb_file):
    antibody, antigen = differ_ab_ag(pdb_file)
    antibody = dict_to_df(antibody)
    antigen = dict_to_df(antigen)
    hydrophobic_amino_acid = ['MET', 'ALA',
                            'VAL', 'LEU', 'ILE', 'PHE', 'PRO', 'TRP']
    Main_chain_atoms = [' N  ', ' CA ', ' C  ', ' O  ']
    antigen_Hydrophobic = antigen[antigen['RESIDUE_NAME'].isin(
        hydrophobic_amino_acid) &
        ~antigen['ATOM_NAME'].isin(
        Main_chain_atoms)]

    antibody_Hydrophobic = antibody[antibody['RESIDUE_NAME'].isin(
        hydrophobic_amino_acid) &
        ~antibody['ATOM_NAME'].isin(
        Main_chain_atoms)]

    antigen_Hydrophobic_Geometric_center_coordinates = pd.pivot_table(
        antigen_Hydrophobic,
        index=[u'chain', u'RESIDUE_NAME', u're_seq']).reset_index()
    antibody_Hydrophobic_Geometric_center_coordinates = pd.pivot_table(
        antibody_Hydrophobic,
        index=[u'chain', u'RESIDUE_NAME', u're_seq']).reset_index()

    return (antigen_Hydrophobic_Geometric_center_coordinates,
            antibody_Hydrophobic_Geometric_center_coordinates)


"""
将保存在字典中的pdb的原子信息转化为dataframe

将抗原和抗体的信息分别进行保存

对dataframe的每一行都进行判断

如果其中的属性属于donor或者receptor

就将新添加的is_donor或者is_receptor添加 TRUE or FALSE

后续通过loc函数的性质将抗原抗体中,氢键的受体或者供体进行区分
"""


def is_donor_new(antibody, antigen):
    antibody_donor = []
    antigen_donor = []
    for i in range(len(antibody)):
        antibody_donor.append(any([antibody.iloc[i, 3] in donor.keys() and
                                antibody.iloc[i, 2].replace(' ', '') in
                                donor[antibody.iloc[i, 3]]]))
    for i in range(len(antigen)):
        antigen_donor.append(any([antigen.iloc[i, 3] in donor.keys() and
                                antigen.iloc[i, 2].replace(' ', '') in
                                donor[antigen.iloc[i, 3]]]))
    return antibody_donor, antigen_donor


def is_receptor_new(antibody, antigen):
    antibody_receptor = []
    antigen_receptor = []
    for i in range(len(antibody)):
        antibody_receptor.append(any([antibody.iloc[i, 3] in receptor.keys()
                                    and antibody.iloc[i, 2].replace(' ', '')
                                    in receptor[antibody.iloc[i, 3]]]))
    for i in range(len(antigen)):
        antigen_receptor.append(any([antigen.iloc[i, 3] in receptor.keys()
                                    and antigen.iloc[i, 2].replace(' ', '')
                                    in receptor[antigen.iloc[i, 3]]]))
    return antibody_receptor, antigen_receptor

#input pdb_file_df output DataFrame
def H_bond_donor_receptor_coordinates_loc_new(pdb_file):
    antibody, antigen = ab_ag_epi_df(pdb_file)
    antibody_donor, antigen_donor = is_donor_new(antibody, antigen)
    antibody_receptor, antigen_receptor = is_receptor_new(antibody, antigen)
    antibody['is_donor'] = antibody_donor
    ab_H_do_coor = antibody.loc[antibody['is_donor']]
    antigen['is_receptor'] = antigen_receptor
    ag_H_re_coor = antigen.loc[antigen['is_receptor']]
    antibody['is_receptor'] = antibody_receptor
    ab_H_re_coor = antibody.loc[antibody['is_receptor']]
    antigen['is_donor'] = antigen_donor
    ag_H_do_coor = antigen.loc[antigen['is_donor']]
    return ab_H_do_coor, ag_H_re_coor, ab_H_re_coor, ag_H_do_coor


"""
主函数

将一个PDB中的氢键和疏水作用力的信息保存在一个dataframe中
"""


def PDB_force_information(pdb_file):
    (antigen_Hydrophobic,
    antibody_Hydrophobic) = Hydrophobic_Geometric_center_coordinates(pdb_file)
    (ab_H_do_coor,
    ag_H_re_coor,
    ab_H_re_coor,
    ag_H_do_coor) = H_bond_donor_receptor_coordinates_loc_new(pdb_file)
    force_information = []
    for i in range(len(antibody_Hydrophobic)):
        for j in range(len(antigen_Hydrophobic)):
            dis = distance(antibody_Hydrophobic.iloc[i],
                        antigen_Hydrophobic.iloc[j])
            if dis <= 6.5 and dis > 2.5:
                force_information.append(['Hydrophobic',
                                        antibody_Hydrophobic.iloc[i].loc['RESIDUE_NAME'],
                                        antibody_Hydrophobic.iloc[i].loc['chain'],
                                        antibody_Hydrophobic.iloc[i].loc['re_seq'],
                                        antigen_Hydrophobic.iloc[j].loc['RESIDUE_NAME'],
                                        antigen_Hydrophobic.iloc[j].loc['chain'],
                                        antigen_Hydrophobic.iloc[j].loc['re_seq'], dis])
    for i in range(len(ab_H_do_coor)):
        for j in range(len(ag_H_re_coor)):
            dis = distance(ab_H_do_coor.iloc[i], ag_H_re_coor.iloc[j])
            if dis <= 3.5 and dis > 2.5:
                force_information.append(['H_bond',
                                        ab_H_do_coor.iloc[i, 3],
                                        ab_H_do_coor.iloc[i, 4],
                                        ab_H_do_coor.iloc[i, 5],
                                        ag_H_re_coor.iloc[j, 3],
                                        ag_H_re_coor.iloc[j, 4],
                                        ag_H_re_coor.iloc[j, 5],
                                        dis])
    for i in range(len(ab_H_re_coor)):
        for j in range(len(ag_H_do_coor)):
            dis = distance(ab_H_re_coor.iloc[i], ag_H_do_coor.iloc[j])
            if dis <= 3.5 and dis > 2.5:
                force_information.append(['H_bond',
                                        ab_H_re_coor.iloc[i, 3],
                                        ab_H_re_coor.iloc[i, 4],
                                        ab_H_re_coor.iloc[i, 5],
                                        ag_H_do_coor.iloc[j, 3],
                                        ag_H_do_coor.iloc[j, 4],
                                        ag_H_do_coor.iloc[j, 5],
                                        dis])
    force_information = DataFrame(force_information)
    force_information.columns = ['Type_bond',
                                'Amino_antibody', 'chain', 're_seq',
                                'Amino_antigen', 'chain', 're_seq',
                                'distance']
    return force_information


"""
主程序

批量提取PDB中的氢键和疏水作用力

分别保存在各自的txt文件中

并生成一个汇总的txt文件
"""
if __name__ == '__main__':
    absolute_path = get_py_path('Index_list.py')
    pdb_name = os.listdir(absolute_path)
    Index_list = pd.DataFrame()
    Index_list_final = pd.DataFrame()
    for i in pdb_name:
        if i[-3:] == 'ent':
            pdb_file = open_pdb_file(absolute_path + '\\' + i)
            antibody_chain, antigen_chain = chain_antibody_antigen(pdb_file)
            if antibody_chain != [] and antigen_chain != []:
                force_information = PDB_force_information(pdb_file)
                Index_list_final = Index_list.append(force_information)
                force_information.to_csv(absolute_path + '\\Index_list\\'
                                        + i[3:7] + '.txt', sep='\t',
                                        index=False, header=True)

    Index_list_final.to_csv(absolute_path + '\\Index_list\\Index_list_final.txt', sep='\t',
                            index=False, header=True)














