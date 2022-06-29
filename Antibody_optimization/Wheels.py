'''
Author: xiaobo 973801194@qq.com
Date: 2022-06-23 12:37:47
LastEditors: xiaobo 973801194@qq.com
LastEditTime: 2022-06-28 18:04:42
FilePath: \Project\Antibody\Index\Antibody_optimization\Wheels.py
Description: dict_to_df ,differ_ab_ag ,distance ,
            chain_antibody_antigen ,dict_pdb_atom ,
            get_py_path functions in this file.
'''
# -*- coding: utf-8 -


import os
from pandas import DataFrame
import pandas

# 将字典中的值转化为df的形式


def dict_to_df(dictionary):
    """_summary_

    Args:
        dictionary (_type_): _description_

    Returns:
        _type_: _description_
    """
    Chain = []
    atom_df = []
    for chain, chain_amino in dictionary.items():
        Chain.append(chain)
    for chain in Chain:
        for re_seq, atom_information in dictionary[chain].items():
            for one_atom in atom_information:
                atom_df.append(one_atom)
    ATOM_df = DataFrame(atom_df)
    ATOM_df.columns = ['ATOM', 'ATOM_NUMBER', 'ATOM_NAME',
                    'RESIDUE_NAME', 'chain', 're_seq',
                    'coordinate_x', 'coordinate_y', 'coordinate_z',
                    'occupancy', 'tempFactor', 'element']

    return ATOM_df

# 将ATOM中的信息，通过抗原抗体的链名进行区分，抗原放在一起，抗体放在一起,save format dict


def differ_ab_ag(pdb_file):
    ATOM = dict_pdb_atom(pdb_file)
    antibody = {}
    antigen = {}
    antibody_chain, antigen_chain = chain_antibody_antigen(pdb_file)
    for i in ATOM.keys():
        if i in antibody_chain:
            antibody[i] = ATOM[i]
        else:
            antigen[i] = ATOM[i]
    return antibody, antigen


"""
三位空间两点距离的计算

输入两行数据,提取其中的x,y,z的坐标

进行两点之间距离的计算
"""


def distance(seires_one, series_two):
    dis = ((seires_one['coordinate_x'] - series_two['coordinate_x'])**2 +
        (seires_one['coordinate_y'] - series_two['coordinate_y'])**2 +
        (seires_one['coordinate_z'] - series_two['coordinate_z'])**2)**0.5
    return round(dis, 3)

# 提取pdb中抗原和抗体链的信息


def chain_antibody_antigen(pdb_file):
    antibody_chain = []
    antigen_chain = []
    for i in range(len(pdb_file)):
        if pdb_file[i][11:17] == 'CHAIN:':
            if (('HEAVY' in pdb_file[i-1])
                or ('LIGHT' in pdb_file[i-1])
                    or ('ANTIBODY' in pdb_file[i-1]) or ('HEAVY' in pdb_file[i-2])
                or ('LIGHT' in pdb_file[i-2])
                    or ('ANTIBODY' in pdb_file[i-2])):
                for chain in pdb_file[i][18:].replace(',', '').replace(
                        ';', '').replace(' ', ''):
                    antibody_chain.append(chain)
            else:
                for chain in pdb_file[i][18:].replace(',', '').replace(
                        ';', '').replace(' ', ''):
                    antigen_chain.append(chain)
    return (antibody_chain, antigen_chain)

# 将pdb的原子属性保存在字典中


def dict_pdb_atom(pdb_file):
    chain = []
    dict_pdb = {}
    pdb_file_ATOM = []
    for i in pdb_file:
        if i[:4] == 'ATOM':
            chain.append(i[21])
            pdb_file_ATOM.append(i)
    for i in chain:
        dict_pdb[i] = {}
    for i in pdb_file_ATOM:
        dict_pdb[i[21]][int(i[22:26])] = []
    for i in pdb_file_ATOM:
        dict_pdb[i[21]][int(i[22:26])].append([i[0:4],
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
    return dict_pdb

# 获取当前.py文件所在的目录


def get_py_path(name_py):
    py_path = os.path.dirname(os.path.realpath(name_py))
    return py_path

"""
pen_pdb_file(pdb_path)

pdb_path为绝对路径

打开pdb文件,并将其作为字符串保存到变量中
"""


def open_pdb_file(pdb_path):
    with open(pdb_path, 'r') as pdb_file:
        pdb_file = pdb_file.read().split("\n")
    return pdb_file


# get pdb.file from current directory
def get_pdb_file(path):
    file_list = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.pdb'):
                file_list.append(os.path.join(root, file))
    return file_list



# get txt.file from current directory
def get_txt_file(path):
    file_list = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.txt'):
                file_list.append(os.path.join(root, file))
    return file_list






