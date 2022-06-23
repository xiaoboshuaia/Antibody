# -*- coding: utf-8 -
from operator import index
import pandas as pd
from data import donor, receptor
from Index_list import differ_ab_ag, dict_to_df
from pandas import DataFrame

# 将Rotamers中的数据保存到dataframe中


def Rotamers(file):

    with open(r'C:\Users\97380\Desktop\第二轮轮转\索引表\索引表\Rotamers.pdb',
              'r') as file:
        file = file.read().split("\n")
    ATOM = []  # 创建一个空的列表
    for i in range(len(file)):  # 对于变量file中的每一个元素进行遍历，即遍历每一行
        atom = []  # 每遍历一个元素，就创建一个空集合atom
        # 如果该字符串是以ATOM开头则将该行进行切片将我需要的信息，添加到atom中，形成一个list
        if file[i].startswith('ATOM'):
            atom = [file[i][0:4],
                    int(file[i][6:11]),
                    file[i][12:16],
                    file[i][17:20],
                    file[i][21],
                    int(file[i][22:26]),
                    float(file[i][30:38]),
                    float(file[i][38:46]),
                    float(file[i][46:54]),
                    file[i][54:60],
                    file[i][60:66],
                    file[i][77:79]]
            ATOM.append(atom)  # 最终将atom添加到ATON中，即一个大列表，其中每个元素是一个小列表，为二维列表
    ATOM = DataFrame(ATOM)  # 将二维列表转换为dataframe格式
    ATOM.columns = ['ATOM', 'ATOM_NUMBER', 'ATOM_NAME',
                    'RESIDUE_NAME', 'chain', 'ro_type',
                    'coordinate_x', 'coordinate_y', 'coordinate_z',
                    'occupancy', 'tempFactor', 'element']
    Rotamers = ATOM
    return Rotamers

# 将rotamers保存到字典中


def dict_ro_atom(rotamers_file):
    amino_name = []
    dict_rotamers = {}
    rotamers_ATOM = []
    for i in rotamers_file:
        if i[:4] == 'ATOM':
            amino_name.append(i[17:20])
            rotamers_ATOM.append(i)
        amino_name_list = set(amino_name)
    for i in amino_name_list:
        dict_rotamers[i] = {}
    for i in rotamers_ATOM:
        dict_rotamers[i[17:20]][int(i[22:26])] = []
    for i in rotamers_ATOM:
        dict_rotamers[i[17:20]][int(i[22:26])].append([i[0:4],
                                                       int(i[6:11]),
                                                       i[12:16],
                                                       i[17:20],
                                                       i[21],
                                                       int(i[22:26]),
                                                       float(i[30:38]),
                                                       float(i[38:46]),
                                                       float(i[46:54]),
                                                       i[54:60],
                                                       i[60:66],
                                                       i[77:79]])
    return dict_rotamers

# 将PDB中ATOM开头的信息储存到命名为‘ATOM’的dataframe中


def atom(file):
    ATOM = []  # 创建一个空的列表
    for i in range(len(file)):  # 对于变量file中的每一个元素进行遍历，即遍历每一行
        atom = []  # 每遍历一个元素，就创建一个空集合atom
        # 如果该字符串是以ATOM开头则将该行进行切片将我需要的信息，添加到atom中，形成一个list
        if file[i].startswith('ATOM'):
            atom = [file[i][0:4],
                    int(file[i][6:11]),
                    file[i][12:16],
                    file[i][17:20],
                    file[i][21],
                    int(file[i][22:26]),
                    float(file[i][30:38]),
                    float(file[i][38:46]),
                    float(file[i][46:54]),
                    file[i][54:60],
                    file[i][60:66],
                    file[i][77:79]]
            ATOM.append(atom)  # 最终将atom添加到ATON中，即一个大列表，其中每个元素是一个小列表，为二维列表
    ATOM = DataFrame(ATOM)  # 将二维列表转换为dataframe格式
    ATOM.columns = ['ATOM', 'ATOM_NUMBER', 'ATOM_NAME',
                    'RESIDUE_NAME', 'chain', 're_seq',
                    'coordinate_x', 'coordinate_y', 'coordinate_z',
                    'occupancy', 'tempFactor', 'element']

    return ATOM

# 找到pdb文件中抗原和抗体链的信息，并将其分别保存为两个list


def find_antibody_antigen_chain(pdb_file):
    COMPND = []
    for i in range(len(pdb_file)):
        B = []
        if pdb_file[i].startswith('COMPND'):
            B = [pdb_file[i]]
            COMPND.append(B)
    antigen_chain_all = []
    antibody_chain_all = []
    for i in range(len(COMPND)):
        if 'SPIKE' in COMPND[i][0]:
            antigen_chain_all.append(COMPND[i+1])
        if 'CHAIN' in COMPND[i][0][11:16] and 'SPIKE' not in COMPND[i-1][0]:
            antibody_chain_all.append(COMPND[i])

    antigen_chain = []
    antibody_chain = []
    for i in antigen_chain_all:
        c = i[0][17:-1].replace(' ', '').replace(';', '').replace(',', '')
        for i in c:
            antigen_chain.append(i)
    for i in antibody_chain_all:
        c = i[0][17:-1].replace(' ', '').replace(';', '').replace(',', '')
        for i in c:
            antibody_chain.append(i)

    return antigen_chain, antibody_chain

# 笨比方法 抗原抗体中的氢键受体与氢键供体区分开


def H_bond_donor_coordinates():
    antibody, antigen = differ_ab_ag()
    antibody_H_bond_donor_coordinates = pd.DataFrame(
        columns=['chain', 'amino', 'ATOM', 're_seq',
                 'coordinate_x', 'coordinate_y', 'coordinate_z'])
    antibody_H_bond_receptor_coordinates = pd.DataFrame(
        columns=['chain', 'amino', 'ATOM', 're_seq',
                 'coordinate_x', 'coordinate_y', 'coordinate_z'])
    antigen_H_bond_donor_coordinates = pd.DataFrame(
        columns=['chain', 'amino', 'ATOM', 're_seq',
                 'coordinate_x', 'coordinate_y', 'coordinate_z'])
    antigen_H_bond_receptor_coordinates = pd.DataFrame(
        columns=['chain', 'amino', 'ATOM', 're_seq',
                 'coordinate_x', 'coordinate_y', 'coordinate_z'])
    for index, antibody in antibody.iterrows():
        if (antibody[1]['amino'] == 'GLY' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'ALA' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'VAL' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'LEU' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'MET' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'ILE' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'SER' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'SER' and antibody[1]['ATOM'] == 'OG'
            or antibody[1]['amino'] == 'THR' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'THR' and antibody[1]['ATOM'] == 'OG1'
            or antibody[1]['amino'] == 'CYS' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'CYS' and antibody[1]['ATOM'] == 'SG'
            or antibody[1]['amino'] == 'ASN' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'ASN' and antibody[1]['ATOM'] == 'ND2'
            or antibody[1]['amino'] == 'GLN' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'GLN' and antibody[1]['ATOM'] == 'NE2'
            or antibody[1]['amino'] == 'PHE' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'TYR' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'TYR' and antibody[1]['ATOM'] == 'OH'
            or antibody[1]['amino'] == 'TRP' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'TRP' and antibody[1]['ATOM'] == 'NE1'
            or antibody[1]['amino'] == 'LYS' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'LYS' and antibody[1]['ATOM'] == 'NZ'
            or antibody[1]['amino'] == 'ARG' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'ARG' and antibody[1]['ATOM'] == 'NE'
            or antibody[1]['amino'] == 'ARG' and antibody[1]['ATOM'] == 'NH1'
            or antibody[1]['amino'] == 'ARG' and antibody[1]['ATOM'] == 'NH2'
            or antibody[1]['amino'] == 'HIS' and antibody[1]['ATOM'] == 'N'
            or antibody[1]['amino'] == 'HIS' and antibody[1]['ATOM'] == 'NE2'
            or antibody[1]['amino'] == 'ASP' and antibody[1]['ATOM'] == 'N'
                or antibody[1]['amino'] == 'GLU' and antibody[1]['ATOM'] == 'N'):
            antibody_H_bond_donor_coordinates.append(
                antibody[1], ignore_index=True)
        elif(antibody[1]['amino'] == 'GLY' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'ALA' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'VAL' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'LEU' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'MET' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'ILE' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'SER' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'THR' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'CYS' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'CYS' and antibody[1]['ATOM'] == 'SG'
             or antibody[1]['amino'] == 'PRO' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'ASN' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'ASN' and antibody[1]['ATOM'] == 'OD1'
             or antibody[1]['amino'] == 'GLN' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'GLN' and antibody[1]['ATOM'] == 'OE1'
             or antibody[1]['amino'] == 'PHE' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'TYR' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'TRP' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'LYS' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'ARG' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'HIS' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'HIS' and antibody[1]['ATOM'] == 'ND1'
             or antibody[1]['amino'] == 'ASP' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'ASP' and antibody[1]['ATOM'] == 'OD1'
             or antibody[1]['amino'] == 'ASP' and antibody[1]['ATOM'] == 'OD2'
             or antibody[1]['amino'] == 'GLU' and antibody[1]['ATOM'] == 'O'
             or antibody[1]['amino'] == 'GLU' and antibody[1]['ATOM'] == 'OE1'
             or antibody[1]['amino'] == 'GLU' and antibody[1]['ATOM'] == 'OE2'):
            antibody_H_bond_receptor_coordinates.append(
                antibody[1], ignore_index=True)
    for antigen in antigen.iterrows():
        if (antigen[1]['amino'] == 'GLY' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'ALA' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'VAL' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'LEU' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'MET' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'ILE' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'SER' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'SER' and antigen[1]['ATOM'] == 'OG'
            or antigen[1]['amino'] == 'THR' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'THR' and antigen[1]['ATOM'] == 'OG1'
            or antigen[1]['amino'] == 'CYS' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'CYS' and antigen[1]['ATOM'] == 'SG'
            or antigen[1]['amino'] == 'ASN' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'ASN' and antigen[1]['ATOM'] == 'ND2'
            or antigen[1]['amino'] == 'GLN' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'GLN' and antigen[1]['ATOM'] == 'NE2'
            or antigen[1]['amino'] == 'PHE' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'TYR' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'TYR' and antigen[1]['ATOM'] == 'OH'
            or antigen[1]['amino'] == 'TRP' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'TRP' and antigen[1]['ATOM'] == 'NE1'
            or antigen[1]['amino'] == 'LYS' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'LYS' and antigen[1]['ATOM'] == 'NZ'
            or antigen[1]['amino'] == 'ARG' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'ARG' and antigen[1]['ATOM'] == 'NE'
            or antigen[1]['amino'] == 'ARG' and antigen[1]['ATOM'] == 'NH1'
            or antigen[1]['amino'] == 'ARG' and antigen[1]['ATOM'] == 'NH2'
            or antigen[1]['amino'] == 'HIS' and antigen[1]['ATOM'] == 'N'
            or antigen[1]['amino'] == 'HIS' and antigen[1]['ATOM'] == 'NE2'
            or antigen[1]['amino'] == 'ASP' and antigen[1]['ATOM'] == 'N'
                or antigen[1]['amino'] == 'GLU' and antigen[1]['ATOM'] == 'N'):
            antigen_H_bond_donor_coordinates.append(
                antigen[1], ignore_index=True)
        elif(antigen[1]['amino'] == 'GLY' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'ALA' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'VAL' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'LEU' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'MET' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'ILE' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'SER' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'THR' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'CYS' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'CYS' and antigen[1]['ATOM'] == 'SG'
             or antigen[1]['amino'] == 'PRO' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'ASN' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'ASN' and antigen[1]['ATOM'] == 'OD1'
             or antigen[1]['amino'] == 'GLN' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'GLN' and antigen[1]['ATOM'] == 'OE1'
             or antigen[1]['amino'] == 'PHE' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'TYR' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'TRP' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'LYS' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'ARG' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'HIS' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'HIS' and antigen[1]['ATOM'] == 'ND1'
             or antigen[1]['amino'] == 'ASP' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'ASP' and antigen[1]['ATOM'] == 'OD1'
             or antigen[1]['amino'] == 'ASP' and antigen[1]['ATOM'] == 'OD2'
             or antigen[1]['amino'] == 'GLU' and antigen[1]['ATOM'] == 'O'
             or antigen[1]['amino'] == 'GLU' and antigen[1]['ATOM'] == 'OE1'
             or antigen[1]['amino'] == 'GLU' and antigen[1]['ATOM'] == 'OE2'):
            antigen_H_bond_receptor_coordinates.append(
                antigen[1], ignore_index=True)

    return (antibody_H_bond_donor_coordinates,
            antigen_H_bond_donor_coordinates,
            antibody_H_bond_receptor_coordinates,
            antigen_H_bond_receptor_coordinates)

    antibody_H_bond_receptor_coordinates = antibody.loc[(antibody['amino'] == 'GLY') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'ALA') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'VAL') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'LEU') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'MET') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'ILE') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'SER') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'THR') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'CYS') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'CYS') & (antibody['ATOM'] == 'SG')
                                                        | (antibody['amino'] == 'PRO') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'ASN') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'ASN') & (antibody['ATOM'] == 'OD1')
                                                        | (antibody['amino'] == 'GLN') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'GLN') & (antibody['ATOM'] == 'OE1')
                                                        | (antibody['amino'] == 'PHE') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'TYR') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'TRP') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'LYS') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'ARG') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'HIS') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'HIS') & (antibody['ATOM'] == 'ND1')
                                                        | (antibody['amino'] == 'ASP') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'ASP') & (antibody['ATOM'] == 'OD1')
                                                        | (antibody['amino'] == 'ASP') & (antibody['ATOM'] == 'OD2')
                                                        | (antibody['amino'] == 'GLU') & (antibody['ATOM'] == 'O')
                                                        | (antibody['amino'] == 'GLU') & (antibody['ATOM'] == 'OE1')
                                                        | (antibody['amino'] == 'GLU') & (antibody['ATOM'] == 'OE2')]
    antigen_H_bond_donor_coordinates = antigen.loc[(antigen['amino'] == 'GLY') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'ALA') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'VAL') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'LEU') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'MET') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'ILE') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'SER') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'SER') & (antigen['ATOM'] == 'OG')
                                                   | (antigen['amino'] == 'THR') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'THR') & (antigen['ATOM'] == 'OG1')
                                                   | (antigen['amino'] == 'CYS') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'CYS') & (antigen['ATOM'] == 'SG')
                                                   | (antigen['amino'] == 'ASN') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'ASN') & (antigen['ATOM'] == 'ND2')
                                                   | (antigen['amino'] == 'GLN') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'GLN') & (antigen['ATOM'] == 'NE2')
                                                   | (antigen['amino'] == 'PHE') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'TYR') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'TYR') & (antigen['ATOM'] == 'OH')
                                                   | (antigen['amino'] == 'TRP') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'TRP') & (antigen['ATOM'] == 'NE1')
                                                   | (antigen['amino'] == 'LYS') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'LYS') & (antigen['ATOM'] == 'NZ')
                                                   | (antigen['amino'] == 'ARG') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'ARG') & (antigen['ATOM'] == 'NE')
                                                   | (antigen['amino'] == 'ARG') & (antigen['ATOM'] == 'NH1')
                                                   | (antigen['amino'] == 'ARG') & (antigen['ATOM'] == 'NH2')
                                                   | (antigen['amino'] == 'HIS') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'HIS') & (antigen['ATOM'] == 'NE2')
                                                   | (antigen['amino'] == 'ASP') & (antigen['ATOM'] == 'N')
                                                   | (antigen['amino'] == 'GLU') & (antigen['ATOM'] == 'N')]
    antigen_H_bond_receptor_coordinates = antigen.loc[(antibody['amino'] == 'GLY') & (antibody['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'ALA') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'VAL') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'LEU') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'MET') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'ILE') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'SER') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'THR') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'CYS') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'CYS') & (antigen['ATOM'] == 'SG')
                                                      | (antigen['amino'] == 'PRO') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'ASN') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'ASN') & (antigen['ATOM'] == 'OD1')
                                                      | (antigen['amino'] == 'GLN') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'GLN') & (antigen['ATOM'] == 'OE1')
                                                      | (antigen['amino'] == 'PHE') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'TYR') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'TRP') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'LYS') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'ARG') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'HIS') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'HIS') & (antigen['ATOM'] == 'ND1')
                                                      | (antigen['amino'] == 'ASP') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'ASP') & (antigen['ATOM'] == 'OD1')
                                                      | (antigen['amino'] == 'ASP') & (antigen['ATOM'] == 'OD2')
                                                      | (antigen['amino'] == 'GLU') & (antigen['ATOM'] == 'O')
                                                      | (antigen['amino'] == 'GLU') & (antigen['ATOM'] == 'OE1')
                                                      | (antigen['amino'] == 'GLU') & (antigen['ATOM'] == 'OE2')]

    return antibody_H_bond_donor_coordinates, antigen_H_bond_donor_coordinates, antibody_H_bond_receptor_coordinates, antigen_H_bond_receptor_coordinates

# 通过便利一遍列表将抗原抗体中的氢键受体与氢键供体区分开，分别保存在不同的dataframe中


def H_bond_donor_receptor_coordinates_for_cycle_old(pdb_file):
    antibody, antigen = differ_ab_ag(pdb_file)
    antibody = dict_to_df(antibody)
    antigen = dict_to_df(antigen)
    antibody_H_bond_donor_coordinates = []
    antibody_H_bond_receptor_coordinates = []
    antigen_H_bond_donor_coordinates = []
    antigen_H_bond_receptor_coordinates = []
    for i in antibody.iterrows():
        if [i[-1]['amino'], i[-1]['ATOM']] in donor:
            antibody_H_bond_donor_coordinates.append(i[-1].tolist())
        elif [i[-1]['amino'], i[-1]['ATOM']] in receptor:
            antibody_H_bond_receptor_coordinates.append(i[-1].tolist())
    for i in antigen.iterrows():
        if [i[-1]['amino'], i[-1]['ATOM']] in donor:
            antigen_H_bond_donor_coordinates.append(i[-1].tolist())
        elif [i[-1]['amino'], i[-1]['ATOM']] in receptor:
            antigen_H_bond_receptor_coordinates.append(i[-1].tolist())
    ab_H_do_coor = DataFrame(antibody_H_bond_donor_coordinates)
    ab_H_re_coor = DataFrame(antibody_H_bond_receptor_coordinates)
    ag_H_do_coor = DataFrame(antigen_H_bond_donor_coordinates)
    ag_H_re_coor = DataFrame(antigen_H_bond_receptor_coordinates)
    return ab_H_do_coor, ab_H_re_coor, ag_H_do_coor, ag_H_re_coor


def H_bond_donor_receptor_coordinates_for_cycle_new(pdb_file):
    antibody, antigen = differ_ab_ag(pdb_file)
    antibody = dict_to_df(antibody)
    antigen = dict_to_df(antigen)
    ab_H_do_coor = []
    ab_H_re_coor = []
    ag_H_do_coor = []
    ag_H_re_coor = []
    for index, row in antibody.iterrows():
        if row['RESIDUE_NAME'] in donor.keys():
            if row['ATOM_NAME'].replace(
                    ' ', '') in donor[row['RESIDUE_NAME']]:
                ab_H_do_coor.append(row)
        if row['RESIDUE_NAME'] in receptor.keys():
            if row['ATOM_NAME'].replace(
                    ' ', '') in receptor[row['RESIDUE_NAME']]:
                ab_H_re_coor.append(row)
    for index, row in antigen.iterrows():
        if row['RESIDUE_NAME'] in donor.keys():
            if row['ATOM_NAME'].replace(
                    ' ', '') in donor[row['RESIDUE_NAME']]:
                ag_H_do_coor.append(row)
        if row['RESIDUE_NAME'] in receptor.keys():
            if row['ATOM_NAME'].replace(
                    ' ', '') in receptor[row['RESIDUE_NAME']]:
                ag_H_re_coor.append(row)
    ab_H_do_coor = DataFrame(ab_H_do_coor)
    ab_H_re_coor = DataFrame(ab_H_re_coor)
    ag_H_do_coor = DataFrame(ag_H_do_coor)
    ag_H_re_coor = DataFrame(ag_H_re_coor)
    return ab_H_do_coor, ab_H_re_coor, ag_H_do_coor, ag_H_re_coor

# 通过dataframe中的loc方法，将抗原抗体中的氢键受体与氢键供体区分开，分别保存在不同的dataframe中


def is_donor_old(antibody, antigen):
    result = []
    for i in range(len(antibody)):

        result.append(
            any([[antibody.iloc[i, 1], antibody.iloc[i, 2]] == j for j in donor]))
    for i in range(len(antigen)):

        result.append(
            any([[antigen.iloc[i, 1], antibody.iloc[i, 2]] == j for j in donor]))

    return result


def is_receptor_old(antibody, antigen):
    result = []
    for i in range(len(antibody)):

        result.append(
            any([[antibody.iloc[i, 1], antibody.iloc[i, 2]] == j for j in receptor]))
    for i in range(len(antigen)):

        result.append(
            any([[antigen.iloc[i, 1], antibody.iloc[i, 2]] == j for j in receptor]))

    return result


def H_bond_donor_receptor_coordinates_loc_old():
    antibody, antigen = differ_ab_ag()
    antibody['is_donor'] = is_donor_old(antibody)
    ab_H_do_coor = antibody.loc[antibody['is_donor']]
    antigen['is_receptor'] = is_receptor_old(antigen)
    ag_H_re_coor = antigen.loc[antigen['is_receptor']]
    antibody['is_receptor'] = is_receptor_old(antibody)
    ab_H_re_coor = antibody.loc[antibody['is_receptor']]
    antigen['is_donor'] = is_donor_old(antigen)
    ag_H_do_coor = antigen.loc[antigen['is_donor']]
    return ab_H_do_coor, ag_H_re_coor, ab_H_re_coor, ag_H_do_coor
















