# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scipy.linalg as linalg



# 两向量弧度的计算


def vector_radian(revolving, target):
    # 创建两个向量
    x = np.array(revolving)
    y = np.array(target)
    # 计算两个向量的模
    m_x = np.sqrt(x.dot(x))
    m_y = np.sqrt(y.dot(y))
    # 计算两个向量的点乘
    dian = x.dot(y)
    # 计算两个向量之间的夹角的余弦值
    cos_ = dian/(m_x*m_y)
    # 将余弦值转化为弧度
    Rotation_radian = np.arccos(cos_)

    return Rotation_radian

# 平面的法向量及旋转轴


def normal_vector(revolving, target):
    # a,b为两向量，两向量叉乘
    return np.cross(revolving, target)

# 旋转矩阵


def rotatio_matrix(revolving, target):
    axis = normal_vector(revolving, target)
    radian = vector_radian(revolving, target)
    return linalg.expm(np.cross(np.eye(3), axis / linalg.norm(axis)*radian))

# 将原始坐标通过旋转矩阵进行旋转后得到的结果


def new_coor(old_coor, revolving, target):
    new_coor = np.dot(old_coor, rotatio_matrix(revolving, target))
    return new_coor

# 将dataframe的index清零，重新统计


def to_0(dataframe):
    dataframe.index = [0] * len(dataframe)
    return dataframe

# 计算三维坐标中两个点相减的结果，保存为list


def subtraction_list(dataframe1, dataframe2, atom1, atom2):
    list_subtraction = dataframe1.loc[dataframe1['ATOM_NAME']
                                      == atom1].iloc[:, 6:9]
    - dataframe2.loc[dataframe2['ATOM_NAME'] == atom2].iloc[:, 6:9]
    return list_subtraction

# 提取抗体中的某个氨基酸，根据其名字以及序列号


def ab_amino(RESIDUE_NAME, Re_seq, chain):
    return antibody.loc[(antibody['RESIDUE_NAME'] == RESIDUE_NAME) &
                        (antibody['re_seq'] == Re_seq) &
                        (antibody['chain'] == chain)]

# 提取rotamers里的一种氨基酸的所有旋转异构体


def ro_amino(RESIDUE_NAME, ro_type):
    return ro.loc[(ro['RESIDUE_NAME'] == RESIDUE_NAME) &
                  (ro['ro_type'] == ro_type)]

# 两组向量


def vector_pairs(amino_dataframe):
    C = np.array(amino_dataframe.loc[amino_dataframe['ATOM_NAME'] == ' C  '].iloc[:, 6:9])
    CA = np.array(amino_dataframe.loc[amino_dataframe['ATOM_NAME'] == ' CA '].iloc[:, 6:9])
    N = np.array(amino_dataframe.loc[amino_dataframe['ATOM_NAME'] == ' N  '].iloc[:, 6:9])
    CA_C = list(C - CA)
    CA_N = list(N - CA)
    CA_C = [round(i, 3) for i in CA_C[0]]
    CA_N = [round(i, 3) for i in CA_N[0]]
    return CA_C, CA_N

# 新建一个dataframe，对刚体的诶个坐标都进行处理,将刚体每个元素旋转过后的数据添加到新的dataframe中


def new_rotamer(old_rotamer1, vector_1, vector_2):
    rotatio_matrix1 = rotatio_matrix(vector_1, vector_2)
    old_rotamer = old_rotamer1.copy()
    for i in range(len(old_rotamer)):
        new_coor = np.dot(old_rotamer.iloc[i, 6:9].tolist(),
                          rotatio_matrix1)
        old_rotamer.iloc[i, 6:9] = [round(x, 3) for x in
                                    new_coor.tolist()]
    return old_rotamer

# 计算平移的x，y，z的距离，使用中心CA作为计算


def move_x_y_z(ab_data, ro_data):
    ab_data = to_0(ab_data)
    ro_data = to_0(ro_data)
    return subtraction_list(ab_data, ro_data, ' CA ', ' CA ').iloc[0, :].tolist()

# 将经过两次旋转的dataframe进行平移


def coordinate_translation(dataframe1, move_list):
    dataframe = dataframe1.copy()
    coordinate_x = [round(i, 3) for i in [*map(lambda x:x + move_list[0],
                                               dataframe.loc[:, 'coordinate_x']
                                               .tolist())]]
    coordinate_y = [round(i, 3) for i in [*map(lambda x:x + move_list[1],
                                               dataframe.loc[:, 'coordinate_y']
                                               .tolist())]]
    coordinate_z = [round(i, 3) for i in [*map(lambda x:x + move_list[2],
                                               dataframe.loc[:, 'coordinate_z']
                                               .tolist())]]
    dataframe.loc[:, 'coordinate_x'] = coordinate_x
    dataframe.loc[:, 'coordinate_y'] = coordinate_y
    dataframe.loc[:, 'coordinate_z'] = coordinate_z
    return dataframe


"""
旋转思路

将dataframe中一个指定位置的氨基酸进行替换

按照氨基酸的序列号对列表进行切分,分成3份,在将新的dataframe加入进去,后重新生成一列,将原子数进行替换

需要将多个dataframe连接到一起
"""

# 按照re_seq的大小将dataframe进行删除


def del_df(df, re_seq):
    return (df.drop(df[df.re_seq > re_seq - 1].index),
            df.drop(df[df.re_seq < re_seq + 1].index))

# 替换，将列表分为3份，再重新将三个列表进行组合，其中加入需要加入替换的列表


def replace_df(ab_df, ro_df, re_seq):
    ab_df_small_re_seq, ab_df_big_re_seq = del_df(ab_df, re_seq)
    return ab_df_small_re_seq.append(ro_df).append(ab_df_big_re_seq)

# 将抗体按照链 分为需要替换的抗体 和不需要进行替换的抗体两部分


def antibody_divide(antibody_df, chain):
    return (antibody.loc[antibody['chain'] == chain],
            antibody.drop(antibody[antibody.chain == chain].index))

# append之前需要将两个dataframe的列名统一


def Unified_column_names(df1, re_seq, chain):
    df = df1.copy()
    df.rename(columns={'ro_type': 're_seq'}, inplace=True)
    df['re_seq'] = re_seq
    df['chain'] = chain
    return df

# 将新生成的dataframe的 原子数从1开始重新计数


def from_1_to_last(antibodydataframe):
    antibodydataframe['ATOM_NUMBER'] = list(
        range(1, len(antibodydataframe) + 1))
    return antibodydataframe

# 将dataframe中一列的元素，转化为字符串，并将其宽度设置为需要的宽度,左对齐


def formatDataFrame_l(df, column, length):
    df[column] = df[column].astype(str).str.ljust(length)
    return df

# 将dataframe中一列的元素，转化为字符串，并将其宽度设置为需要的宽度,右对齐


def formatDataFrame_r(df, column, length):
    df[column] = df[column].astype(str).str.rjust(length)
    return df

# 将dataframe转化为需要的格式


def all_formatDataFrame(df):
    ATOM = formatDataFrame_l(df, 'ATOM', 6)
    ATOM_NUMBE = formatDataFrame_r(ATOM, 'ATOM_NUMBER', 5)
    ATOM_NAME = formatDataFrame_l(ATOM_NUMBE, 'ATOM_NAME', 6)
    RESIDUE_NAME = formatDataFrame_l(ATOM_NAME, 'RESIDUE_NAME', 4)
    chain = formatDataFrame_l(RESIDUE_NAME, 'chain', 1)
    re_seq = formatDataFrame_r(chain, 're_seq', 4)
    coordinate_x = formatDataFrame_r(re_seq, 'coordinate_x', 12)
    coordinate_y = formatDataFrame_r(coordinate_x, 'coordinate_y', 8)
    coordinate_z = formatDataFrame_r(coordinate_y, 'coordinate_z', 8)
    occupancy = formatDataFrame_r(coordinate_z, 'occupancy', 6)
    tempFactor = formatDataFrame_r(occupancy, 'tempFactor', 6)
    element = formatDataFrame_r(tempFactor, 'element', 12)
    return element


# 获得O的三维坐标
def get_O(antibody_df1):
    antibody_df = antibody_df1.copy()
    return antibody_df.loc[(antibody_df.re_seq == 1) &
                           (antibody_df.RESIDUE_NAME == 'GLN') &
                           (antibody_df.ATOM_NAME == ' O  ')
                           ].iloc[0, 6:9].tolist()

# 将O的三位坐标进行替换


def replace_O(ro_df1, antibody_df, re_seq, RESIDUE_NAME):
    ro_df = ro_df1.copy()
    ro_df.loc[ro_df.ATOM_NAME == ' O  ',
              'coordinate_x':'coordinate_z'] = get_O(antibody_df)
    return ro_df


"""
旋转主函数

将待旋转的rotamer和目标氨基酸的dataframe取出,便于下一步的旋转
    
第一次旋转，将两向量所在的平面的法向量进行旋转，使两两向量共面
    
两平面的法向量
        
将第一次旋转过后的romater的坐标保存到新的dataframe中

进行第二次旋转两两向量已经共面,接下来就只需要将CA_C和CA_N中的一条向量旋转到平行,这里选择CA_C

两次旋转结束后,接下来进行平移,将替换的rotamer移动到目标氨基酸的位置

已经将原来的romater旋转两次,并进行了平移,下一步需要将romater与目标氨基酸进行替换

将目标氨基酸的O原子的坐标,与旋转平移后的坐标进行交换

再合并之前需要将rotamer的dataframe的格式,进行统一

将已经改变的抗体链与剩余的链和抗原进行结合

将重新生成的dataframe保存为pdb文件进行输出
"""


def amino_acid_substitution_main(antibody, antigen, ro, re_seq, ab_RESIDUE_NAME, chain, ro_RESIDUE_NAME):
    for ro_type in set(ro.loc[ro.RESIDUE_NAME == ro_RESIDUE_NAME].ro_type):
        # 将待旋转的rotamer和目标氨基酸的dataframe取出，便于下一步的旋转
        # 第一次旋转，将两向量所在的平面的法向量进行旋转，使两两向量共面
        ro_CA_C, ro_CA_N = vector_pairs(ro_amino(ro_RESIDUE_NAME, ro_type))
        ab_CA_C, ab_CA_N = vector_pairs(
            ab_amino(ab_RESIDUE_NAME, re_seq, chain))
    # 两平面的法向量
        ro_normal_vector = normal_vector(ro_CA_C, ro_CA_N)
        ab_normal_vector = normal_vector(ab_CA_C, ab_CA_N)
    # 将第一次旋转过后的romater的坐标保存到新的dataframe中
        ro_first_rotation = new_rotamer(ro_amino(ro_RESIDUE_NAME, ro_type),
                                        ab_normal_vector, ro_normal_vector)
    # 进行第二次旋转，两两向量已经共面，接下来就只需要将CA_C和CA_N中的一条向量旋转到平行，这里选择CA_C
        ro_CA_C_new, ro_CA_N_new = vector_pairs(ro_first_rotation)
        ro_second_rotation = new_rotamer(
            ro_first_rotation, ab_CA_C, ro_CA_C_new)
    # 两次旋转结束后，接下来进行平移，将替换的rotamer移动到目标氨基酸的位置
        move_list = move_x_y_z(ab_amino(ab_RESIDUE_NAME, re_seq, chain),
                               ro_second_rotation)
        ro_final = coordinate_translation(ro_second_rotation, move_list)
    # 已经将原来的romater旋转两次，并进行了平移，下一步需要将romater与目标氨基酸进行替换
        antibody_goal, antobody_last = antibody_divide(antibody, chain)
    # 将目标氨基酸的O原子的坐标，与旋转平移后的坐标进行交换
        ro_final = replace_O(ro_final, antibody_goal, re_seq, ab_RESIDUE_NAME)
    # 再合并之前需要将rotamer的dataframe的格式，进行统一
        ro_final = Unified_column_names(ro_final, re_seq, chain)
        antibody_replace_chain = replace_df(antibody_goal, ro_final, re_seq)
    # 将已经改变的抗体链与剩余的链和抗原进行结合
        antibody = antobody_last.append(antibody_replace_chain)
        new_antibody = all_formatDataFrame(from_1_to_last(antigen.append
                                                          (antibody)))
        new_antibody.insert(loc=2, column='kong', value=' ')
    # 将重新生成的dataframe保存为pdb文件进行输出

        new_antibody.to_csv(ab_RESIDUE_NAME + str(re_seq) + ro_RESIDUE_NAME
                            + str(ro_type) + '.txt', sep='\t',
                            index=False, header=False)
        with open(ab_RESIDUE_NAME + str(re_seq) + ro_RESIDUE_NAME
                  + str(ro_type) + '.txt', 'r') as want:
            want1 = want.read()
        want2 = want1.replace('\t', '').split('\n')
        huan = '\n'
        f = open(ab_RESIDUE_NAME + str(re_seq) + ro_RESIDUE_NAME
                 + str(ro_type) + '.pdb', 'w', encoding='gbk')
        f.write(huan.join(want2))
        f.close()
