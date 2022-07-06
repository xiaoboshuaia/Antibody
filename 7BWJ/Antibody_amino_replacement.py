"""
Author: xiaobo 973801194@qq.com
Date: 2022-06-24 16:22:02
LastEditors: xiaobo 973801194@qq.com
LastEditTime: 2022-07-02 18:37:56
FilePath: \Project\Antibody\Index\Antibody_optimization\Antibody_amino_replacement.py
Description:
"""
# TODO:   通过字典的方式,将目标氨基酸进行替换
import os
import sys

# 将索引表添加到系统路径中
sys.path.append('d:\Project\Antibody\Index\Antibody_optimization')
import numpy
import pandas
from Data_procress import atom
from amino_acid_substitution import Unified_column_names, all_formatDataFrame, coordinate_translation, move_x_y_z, \
    new_rotamer, normal_vector, vector_pairs
from data import rotamers, top_5_amino_antibody, donor, receptor
from epitope_judgment import ab_ag_correspondence, ab_ag_epi_df
from Wheels import dict_pdb_atom, dict_to_df, differ_ab_ag, distance, get_pdb_file, get_txt_file, open_pdb_file

#input a amino's DataFrame,output two DataFrame,one is the amino's donor,the other is the amino's receptor
#将氨基酸的数据框输入，输出两个数据框，一个是氨基酸的donor，另一个是氨基酸的receptor，通过any的方法来构建一个新的由

def get_donor_receptor(amino_df):
    def judgment(amino_df):
        judgment_donor = []
        judgment_receptor = []
        for i in range(len(amino_df)):
            judgment_donor.append(any([amino_df.iloc[i, 3] in donor.keys() and
                                       amino_df.iloc[i, 2].replace(' ', '') in
                                       donor[amino_df.iloc[i, 3]]]))
            judgment_receptor.append(any([amino_df.iloc[i, 3] in receptor.keys() and
                                          amino_df.iloc[i, 2].replace(' ', '') in
                                          receptor[amino_df.iloc[i, 3]]]))
        return judgment_donor, judgment_receptor

    amino_df = amino_df.copy()  # copy一个新的数据框,因为原数据框是不可变的，因为需要对原始数据来进行操作，但是又不能改变原始数据，所以就需要copy一个新的数据框
    amino_df.loc[:, 'judgment_donor'] = judgment(amino_df)[0]
    amino_df['judgment_receptor'] = judgment(amino_df)[1]
    return amino_df.loc[amino_df['judgment_donor']], amino_df.loc[amino_df['judgment_receptor']]

# get can't replace antibody amino input txt_path output list [['chain','re_seq']]
def get_cannot_replace_amino(txt_file_path):
    cannot_replace_amino = []
    force_info = pandas.read_csv(txt_file_path, delim_whitespace=True)
    for i in range(len(force_info)):
        cannot_replace_amino.append(
            [force_info.iloc[i, 2], force_info.iloc[i, 3], force_info.iloc[i, 1]])  # [antibody_chain,re_seq]
    return cannot_replace_amino

# get amino_dataFrame from big dataFrame by a list which store the amino info like ['chain','re_seq','RESIDUE_NAME']
def get_amino_dataFrame(big_dataFrame, amino_info_list):
    amino_dataFrame = big_dataFrame.loc[
        (big_dataFrame['chain'] == amino_info_list[0]) & (big_dataFrame['re_seq'] == amino_info_list[1])]
    return amino_dataFrame

# get top5 amino list from top_5_amino_antibody
def top5_amino(ag_amino_name):
    top5_amino = []
    for i in top_5_amino_antibody[ag_amino_name]:
        top5_amino.append(i[0])
    return top5_amino

# del cannot replace amino from antibody_amino_list
def del_cannot_replace_amino(antibody_amino_list, cannot_replace_amino):
    for i in cannot_replace_amino:
        if i in antibody_amino_list:
            antibody_amino_list.remove(i)
    return antibody_amino_list

# get after rotation and translation rotamer dataFrame ,input is the rotamer,antibody dataFrame,output is the rotamer dataFrame
def final_rotamer_dataFrame(rotamer_amino_df, antibody_amino_df):
    ro_normal_vector = normal_vector(vector_pairs(rotamer_amino_df)[0], vector_pairs(rotamer_amino_df)[1])
    ab_normal_vector = normal_vector(vector_pairs(antibody_amino_df)[0], vector_pairs(antibody_amino_df)[1])
    ro_first_rotation = new_rotamer(rotamer_amino_df, ab_normal_vector, ro_normal_vector)
    ro_second_rotation = new_rotamer(ro_first_rotation, vector_pairs(antibody_amino_df)[0],
                                     vector_pairs(ro_first_rotation)[0])
    move_list = move_x_y_z(antibody_amino_df, ro_second_rotation)
    ro_final = coordinate_translation(ro_second_rotation, move_list)
    ro_final.loc[ro_final['ATOM_NAME'] == ' O  '] = antibody_amino_df.loc[antibody_amino_df['ATOM_NAME'] == ' O  ']
    return ro_final

# 统一dataFrame的格式
def Unified_df_format(df_goal, df_exp):
    df = df_goal.copy()
    df['chain'] = df_exp.iloc[1,4]
    df['re_seq'] = df_exp.iloc[1,5]
    df['RESIDUE_NAME'] = df_goal.iloc[1,3]
    return df

# 使用字典的方式，将替换的rotamer与抗体上的氨基酸进行替换，并输出为list
def replace_pdb_list(pdb_dict, rotamer_df, antibody_chain, antibody_re_seq):
    pdb_dict = pdb_dict.copy()
    pdb_dict[antibody_chain][antibody_re_seq] = rotamer_df.values.tolist()
    pdb_df = dict_to_df(pdb_dict)
    pdb_df['ATOM_NUMBER'] = list(range(1, len(pdb_df) + 1))
    pdb_list = all_formatDataFrame(pdb_df).values.tolist()
    return pdb_list

# 将该抗体对应的所有抗原氨基酸对应的替换rotamer前五位的氨基酸，放进一个list中，便于下一步替换
def replace_romater_name(ag_amino_list):
    rotamer_name = []
    for ag_amino in ag_amino_list:
        ag_amino_RESIDUE_NAME = ag_amino[2]
        for i in top5_amino(ag_amino_RESIDUE_NAME):
            if i not in rotamer_name:
                rotamer_name.append(i)
    return rotamer_name    







text_pdb = open_pdb_file(r'D:\Project\Antibody\Index\7BWJ\7bwj.pdb')

antibody, antigen = ab_ag_epi_df(text_pdb)

ab_ag_correspondence = ab_ag_correspondence(antibody, antigen)

ag_key = list(ab_ag_correspondence.keys())[0]
ab_key = ab_ag_correspondence[ag_key]

ep_antigen = antigen.loc[(antigen['chain'] == ag_key[0]) & (antigen['re_seq'] == ag_key[1])]
ep_antibody = antibody.loc[(antibody['chain'] == ab_key[0]) & (antibody['re_seq'] == ab_key[1])]

ep_antigen_name = ep_antigen.iloc[0]['RESIDUE_NAME']
ep_antibody_name = ep_antibody.iloc[0]['RESIDUE_NAME']

ag_info = {ep_antigen_name: ag_key}
ab_info = {ep_antibody_name: ab_key}

ab_SER = antibody.loc[
    (antibody['RESIDUE_NAME'] == 'SER') & (antibody['chain'] == ab_key[0]) & (antibody['re_seq'] == ab_key[1])]

# get the range can create the H-bond from antigen's H-bond donor and receptor
x_range, y_range, z_range = [], [], []
for i in range(len(ep_antigen)):
    if donor[ep_antigen.iloc[i, 3]] == ep_antigen.iloc[i, 2].replace(' ', '') or receptor[ep_antigen.iloc[i, 3]] == \
            ep_antigen.iloc[i, 2].replace(' ', ''):
        x_range.append(ep_antigen.iloc[i, 6])
        y_range.append(ep_antigen.iloc[i, 7])
        z_range.append(ep_antigen.iloc[i, 8])
x_range = numpy.arange(min(x_range) - 3.501, max(x_range) + 3.501, 0.001)
y_range = numpy.arange(min(y_range) - 3.501, max(y_range) + 3.501, 0.001)
z_range = numpy.arange(min(z_range) - 3.501, max(z_range) + 3.501, 0.001)

# find the top 5 replace amino in Index_list
replace_amino_name = []
for i in top_5_amino_antibody[ep_antigen_name]:
    replace_amino_name.append(i[0])

ASN_1 = atom(rotamers['ASN'][1])

# rotation and translation the rotamer
# ro_CA_C = vector_pairs(ASN_1)[0];ro_C_N = vector_pairs(ASN_1)[1]
# ab_CA_C = vector_pairs(ab_SER)[0];ab_C_N = vector_pairs(ab_SER)[1]
ro_normal_vector = normal_vector(vector_pairs(ASN_1)[0], vector_pairs(ASN_1)[1])
ab_normal_vector = normal_vector(vector_pairs(ab_SER)[0], vector_pairs(ab_SER)[1])

ro_first_rotation = new_rotamer(ASN_1, ab_normal_vector, ro_normal_vector)

ro_second_rotation = new_rotamer(ro_first_rotation, vector_pairs(ab_SER)[0], vector_pairs(ro_first_rotation)[0])

move_list = move_x_y_z(ab_SER, ro_second_rotation)

ro_final = coordinate_translation(ro_second_rotation, move_list)

ro_final.loc[ro_final['ATOM_NAME'] == ' O  '] = ab_SER.loc[ab_SER['ATOM_NAME'] == ' O  ']
ro_final = Unified_column_names(ro_final, ab_key[1], ab_key[0])

# judgment the H-bond donor and receptor isin the range (unnecessary)
probability_info_list = []
for i in range(len(ro_final)):
    if ro_final.iloc[i, 6] in x_range and ro_final.iloc[i, 7] in y_range and ro_final.iloc[i, 8] in z_range:
        probability_info_list.append([ab_info, ag_info, ro_final.iloc[i, [2, 3]]])

#计算两个氨基酸之间的氢键受体的配体的距离，输入两个dataFrame，输出一个保存距离符合形成氢键的list
def dis_list(ag_df,ro_df):
    force_information = []
    ag_donor, ag_receptor = get_donor_receptor(ag_df)
    ro_donor, ro_receptor = get_donor_receptor(ro_df)
    for i in range(len(ag_donor)):
        for j in range(len(ro_receptor)):
            dis = distance(ag_donor.iloc[i], ro_receptor.iloc[j])
            if dis <= 3.5 and dis > 2.5:
                force_information.append(dis)
    for i in range(len(ag_receptor)):
        for j in range(len(ro_donor)):
            dis = distance(ag_receptor.iloc[i], ro_donor.iloc[j])
            if dis <= 3.5 and dis > 2.5:
                force_information.append(dis)
    return dis


if __name__ == '__main__':
    """将interface的抗体氨基酸进行替换,并计算替换后的抗体能产生的氢键 
    
    """
    pdb_file_list = open_pdb_file(get_pdb_file(os.getcwd())[0])  # get_pdb_file得到的是list，由于只有一个元素，所以要用[0]
    pdb_file_dict = dict_pdb_atom(pdb_file_list)  # 将pdb文件转换成字典
    pdb_file_df = dict_to_df(pdb_file_dict)  # 将字典转换成dataFrame

    cannot_replace_amino_list = get_cannot_replace_amino(get_txt_file(os.getcwd())[0])  # 获得不能替换的抗体氨基酸的list

    antibody_dict, antigen_dict = differ_ab_ag(pdb_file_list)  # 得到抗体和受体的字典

    interface_antibody_df, interface_antigen_df = ab_ag_epi_df(pdb_file_list)  # 得到interface的抗体和受体的dataFrame
#换成表位抗体对应的抗原
    ab_ag_cor_dict = ab_ag_correspondence(interface_antigen_df,
                                          interface_antibody_df)  # 得到interface的抗体和抗原的对应关系即CA距离小于19A

    interface_antibody_list = list(ab_ag_cor_dict.keys())  # 得到interface的抗体的list
    #从所有表位的抗体中删除不能进行替换的抗体
    interface_antibody_list_new = [i for i in interface_antibody_list if list(i) not in cannot_replace_amino_list]
    for ab_amino in interface_antibody_list_new:

        print(ab_amino)

        ab_amino_df = get_amino_dataFrame(interface_antibody_df, ab_amino)  # 得到一个抗原氨基酸的dataFrame
        ab_amino_chain, ab_amino_re_seq, ab_amino_RESIDUE_NAME = ab_amino[0], ab_amino[1], ab_amino[2]
        # 得到抗原氨基酸对应抗体的list，在抗体的list中删除不能替换的抗体氨基酸
        # TODO 这里可能需要再删除重复元素之前，对字典的值做一个copy以免删去原本的值
        #ab_amino_list = del_cannot_replace_amino(ab_ag_cor_dict[ag_amino], cannot_replace_amino_list)  # 删除不能替换的氨基酸
        #取得抗体对应的所有抗原的index中的前五个rotamer的交集
        rotamer_name = []
        for ag_amino in ab_ag_cor_dict[ab_amino]:
            ag_amino_RESIDUE_NAME = ab_ag_cor_dict[ab_amino][2]
            rotamer_name.append(i for i in top5_amino(ag_amino_RESIDUE_NAME) if i not in rotamer_name)  # 得到替换抗体的top5氨基酸名字的list

        for ab_amino in ab_ag_cor_dict[ag_amino]:

            print(ab_amino)

            ab_amino_df = get_amino_dataFrame(interface_antibody_df, ab_amino)  # 得到一个抗体氨基酸的dataFrame
            ab_amino_chain, ab_amino_re_seq, ab_amino_RESIDUE_NAME = ab_amino[0], ab_amino[1], ab_amino[2]

            for rotamer in rotamer_name:
                for ro_type in list(rotamers[rotamer].keys()):

                    print(rotamer, ro_type)

                    rotamer_df = atom(rotamers[rotamer][ro_type])
                    # 得到旋转平移后的rotamer的dataFrame
                    ro_final_df = final_rotamer_dataFrame(rotamer_df, ab_amino_df)
                    # 将旋转平移后的rotamer的dataFrame格式化
                    ro_final_df = Unified_df_format(ro_final_df, ab_amino_chain, ab_amino_re_seq,
                                                    rotamer)  # ab_amino[0]是chain，ab_amino[1]是re_seq
                    
                    if youxindeqingjian:#抗原和rotamer上都是疏水氨基酸就可以输出，或者能形成氢键在进行输出
                        replace_pdb_list = replace_pdb_list(pdb_file_dict, ro_final_df, ab_amino_chain, ab_amino_re_seq)
                        # 将替换后的蛋白质进行输出
                        f = open(ab_amino_chain + ab_amino_re_seq + rotamer + ro_type + ".pdb", "w")
                        for i in replace_pdb_list:
                            f.writelines(i + ['\n'])
                        f.close()

# TODO def x (ro,antigen):
#            return 氢键信息的字典 
#            字典 = x(ro,antigen)
#            if 字典非空：
#               一个替换过rotamer的df = 函数（ro）    
#               格式化df
#               将df生成pdb文件，名字为（链_re_seq_被替换氨基酸名字_替换的rotamer的名字_ro_type）
#                  字典保存成txt文件
#

'''
将第一个抗体作为测试
text
'''
ab_amino_df = get_amino_dataFrame(interface_antibody_df, interface_antibody_list[0])
#ab_amino_list = del_cannot_replace_amino(ab_ag_cor_dict[interface_antigen_list[0]], cannot_replace_amino_list)
interface_antibody_list_new = [i for i in interface_antibody_list if list(i) not in cannot_replace_amino_list]

ab_amino_THR = interface_antibody_list_new[0]
ag_amino_THR = ab_ag_cor_dict[ab_amino_THR]

rotamer_name = replace_romater_name(ag_amino_THR)

rotamer_LEU_type  = list(rotamers[rotamer_name[0]].keys())
rotamer_LEU = atom(rotamers[rotamer_name[0]][rotamer_LEU_type[0]])

ro_LEU = final_rotamer_dataFrame(rotamer_LEU, ab_amino_df)
ro_LEU = Unified_df_format(ro_LEU, ab_amino_df)

LEU_donor,LEU_receptor = get_donor_receptor(ro_LEU)







ab_amino_df = get_amino_dataFrame(interface_antibody_df, ab_ag_cor_dict[interface_antigen_list[0]][0])
rotamer_df = atom(rotamers['ASN'][1])
ro_final_df = final_rotamer_dataFrame(rotamer_df, ab_amino_df)

ro_final_df['re_seq'] = ab_ag_cor_dict[interface_antigen_list[0]][0][1]  # ab_amino_info =
ro_final_df['chain'] = ab_ag_cor_dict[interface_antigen_list[0]][0][0]
ro_final_df['RESIDUE_NAME'] = rotamer_name[0]

ab_chain = ab_ag_cor_dict[interface_antigen_list[0]][0][0]
ab_re_seq = ab_ag_cor_dict[interface_antigen_list[0]][0][1]





pdb_file_dict_replace = pdb_file_dict.copy()
pdb_file_dict_replace[ab_chain][ab_re_seq] = ro_final_df.values.tolist()
pdb_file_df_replace = dict_to_df(pdb_file_dict_replace)

pdb_file_df_replace['ATOM_NUMBER'] = list(range(1, len(pdb_file_df_replace) + 1))

pdb_file_df_replace.to_csv('text.txt', sep=' ',
                           index=False, header=False)

pdb_file_df_replace = all_formatDataFrame(pdb_file_df_replace)

pdb_file_list_replace = pdb_file_df_replace.values.tolist()

f = open("text.pdb", "w")
for i in pdb_file_list_replace:
    f.writelines(i + ['\n'])
f.close()

# 获得旋转后rotamer的pdb文件
ro_final_df = all_formatDataFrame(ro_final_df)
ro_file_list = ro_final_df.values.tolist()
f = open("ASN.pdb", "w")
for i in ro_file_list:
    f.writelines(i + ['\n'])
f.close()

# 初始rotamer
rotamer_df = atom(rotamers['ASN'][1])
rotamer_df = all_formatDataFrame(rotamer_df)
rotamer_df_list = rotamer_df.values.tolist()
f = open("ASN_begin.pdb", "w")
for i in rotamer_df_list:
    f.writelines(i + ['\n'])
f.close()

# 被替换的抗体氨基酸
ab_amino_df = get_amino_dataFrame(interface_antibody_df, ab_ag_cor_dict[interface_antigen_list[0]][0])
ab_amino_df = ab_amino_df.copy()
ab_amino_df = all_formatDataFrame(ab_amino_df)
ab_amino_df_list = ab_amino_df.values.tolist()
f = open("SER.pdb", "w")
for i in ab_amino_df_list:
    f.writelines(i + ['\n'])
f.close()

# 将CA移动到替换的抗体氨基酸的CA上
move_list = move_x_y_z(ab_amino_df, rotamer_df)
ro_final_df = coordinate_translation(rotamer_df, move_list)
ro_final_df = all_formatDataFrame(ro_final_df)
ro_file_list = ro_final_df.values.tolist()
f = open("ASN_shift.pdb", "w")
for i in ro_file_list:
    f.writelines(i + ['\n'])
f.close()

# 将两个CA进行旋转重叠
# 
ro_first_rotamer = new_rotamer(ro_final_df, vector_pairs(ab_amino_df)[0], vector_pairs(ro_final_df)[0])
ro_first_rotamer = ro_first_rotamer.copy()
ro_first_rotamer = all_formatDataFrame(ro_first_rotamer)
ro_first_rotamer_list = ro_first_rotamer.values.tolist()
f = open("ASN_first_rotation.pdb", "w")
for i in ro_first_rotamer_list:
    f.writelines(i + ['\n'])
f.close()


# get a number from a list
def get_number(list):
    number = []
    for i in list:
        number.append(i[0])
    return number


list()
