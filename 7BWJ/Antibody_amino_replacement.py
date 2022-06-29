'''
Author: xiaobo 973801194@qq.com
Date: 2022-06-24 16:22:02
LastEditors: xiaobo 973801194@qq.com
LastEditTime: 2022-06-29 11:52:50
FilePath: \Project\Antibody\Index\Antibody_optimization\Antibody_amino_replacement.py
Description:
'''
import os
import sys
# 将索引表添加到系统路径中
sys.path.append('d:\Project\Antibody\Index\Antibody_optimization')
import numpy
import pandas
from Data_procress import atom
from amino_acid_substitution import Unified_column_names, coordinate_translation, move_x_y_z, new_rotamer, normal_vector, vector_pairs
from data import rotamers,top_5_amino_antibody,donor,receptor
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
    amino_df = amino_df.copy()#copy一个新的数据框,因为原数据框是不可变的，因为需要对原始数据来进行操作，但是又不能改变原始数据，所以就需要copy一个新的数据框
    amino_df.loc[:,'judgment_donor'] = judgment(amino_df)[0]
    amino_df['judgment_receptor'] = judgment(amino_df)[1]
    return amino_df.loc[amino_df['judgment_donor']],amino_df.loc[amino_df['judgment_receptor']]


#get can't replace antibody amino input txt_path output list [['chain','re_seq']]
def get_cannot_replace_amino(txt_file_path):
    cannot_replace_amino = []
    force_info = pandas.read_csv(txt_file_path,delim_whitespace=True)
    for i in range(len(force_info)):
        cannot_replace_amino.append([force_info.iloc[i,2],force_info.iloc[i,3],force_info.iloc[i,1]])#[antibody_chain,re_seq]
    return cannot_replace_amino

# get amino_dataFrame from big dataFrame by a list which store the amino info like ['chain','re_seq','RESIDUE_NAME']
def get_amino_dataFrame(big_dataFrame,amino_info_list):
    amino_dataFrame = big_dataFrame.loc[(big_dataFrame['chain'] == amino_info_list[0]) & (big_dataFrame['re_seq'] == amino_info_list[1])]
    return amino_dataFrame

# get top5 amino list from top_5_amino_antibody
def top5_amino(ag_amino_name):
    top5_amino = []
    for i in top_5_amino_antibody[ag_amino_name]:
        top5_amino.append(i[0])
    return top5_amino









   
text_pdb = open_pdb_file(r'D:\Project\Antibody\Index\7BWJ\7bwj.pdb')


antibody, antigen = ab_ag_epi_df(text_pdb)


ab_ag_correspondence = ab_ag_correspondence(antibody, antigen)

ag_key = list(ab_ag_correspondence.keys())[0]
ab_key = ab_ag_correspondence[ag_key]

ep_antigen = antigen.loc[(antigen['chain'] == ag_key[0]) & (antigen['re_seq'] == ag_key[1])]
ep_antibody = antibody.loc[(antibody['chain'] == ab_key[0]) & (antibody['re_seq'] == ab_key[1])]


ep_antigen_name = ep_antigen.iloc[0]['RESIDUE_NAME'] 
ep_antibody_name = ep_antibody.iloc[0]['RESIDUE_NAME']

ag_info = {ep_antigen_name:ag_key}
ab_info = {ep_antibody_name:ab_key}

ab_SER = antibody.loc[(antibody['RESIDUE_NAME'] == 'SER') & (antibody['chain'] == ab_key[0]) & (antibody['re_seq'] == ab_key[1])]

#get the range can create the H-bond from antigen's H-bond donor and receptor
x_range, y_range, z_range = [], [], []
for i in range(len(ep_antigen)):
    if donor[ep_antigen.iloc[i, 3]] == ep_antigen.iloc[i, 2].replace(' ', '' ) or receptor[ep_antigen.iloc[i, 3]] == ep_antigen.iloc[i, 2].replace(' ', ''):
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


#rotation and translation the rotamer
#ro_CA_C = vector_pairs(ASN_1)[0];ro_C_N = vector_pairs(ASN_1)[1]
#ab_CA_C = vector_pairs(ab_SER)[0];ab_C_N = vector_pairs(ab_SER)[1]
ro_normal_vector = normal_vector(vector_pairs(ASN_1)[0], vector_pairs(ASN_1)[1])
ab_normal_vector = normal_vector(vector_pairs(ab_SER)[0], vector_pairs(ab_SER)[1])

ro_first_rotation = new_rotamer(ASN_1, ab_normal_vector, ro_normal_vector)

ro_second_rotation = new_rotamer(ro_first_rotation, vector_pairs(ab_SER)[0], vector_pairs(ro_first_rotation)[0])

move_list = move_x_y_z(ab_SER,ro_second_rotation)

ro_final = coordinate_translation(ro_second_rotation,move_list)

ro_final.loc[ro_final['ATOM_NAME'] == ' O  '] = ab_SER.loc[ab_SER['ATOM_NAME'] == ' O  ']
ro_final = Unified_column_names(ro_final,ab_key[1],ab_key[0])


#judgment the H-bond donor and receptor isin the range (unnecessary)
probability_info_list = []
for i in range(len(ro_final)):
    if ro_final.iloc[i,6] in x_range and ro_final.iloc[i,7] in y_range and ro_final.iloc[i,8] in z_range:
        probability_info_list.append([ab_info,ag_info,ro_final.iloc[i,[2,3]]])


force_information = []
ag_donor,ag_receptor = get_donor_receptor(ep_antigen)    
ro_donor,ro_receptor = get_donor_receptor(ro_final)
for i in range(len(ag_donor)):
    for j in range(len(ro_receptor)):
        dis = distance(ag_donor.iloc[i],ro_receptor.iloc[j])
        if dis <= 3.5 and dis > 2.5:
            force_information.append(['H_bond',
                                        ag_donor.iloc[i, 3],
                                        ag_donor.iloc[i, 4],
                                        ag_donor.iloc[i, 5],
                                        ro_receptor.iloc[j, 3],
                                        ro_receptor.iloc[j, 4],
                                        ro_receptor.iloc[j, 5],
                                        dis,ep_antibody_name])
            
        
for i in range(len(ag_receptor)):
    for j in range(len(ro_donor)):  
        dis = distance(ag_receptor.iloc[i],ro_donor.iloc[j])
        if dis <= 3.5 and dis > 2.5: 
            force_information.append(['H_bond',
                                        ag_receptor.iloc[i, 3],
                                        ag_receptor.iloc[i, 4],
                                        ag_receptor.iloc[i, 5],
                                        ro_donor.iloc[j, 3],
                                        ro_donor.iloc[j, 4],
                                        ro_donor.iloc[j, 5],
                                        dis,ep_antibody_name])



if __name__ == '__main__':
    pdb_file_list = open_pdb_file(get_pdb_file(os.getcwd())[0])#get_pdb_file得到的是list，由于只有一个元素，所以要用[0]
    pdb_file_dict = dict_pdb_atom(pdb_file_list)
    pdb_file_df = dict_to_df(pdb_file_dict)
    cannot_replace_amino_dict = get_cannot_replace_amino(get_txt_file(os.getcwd())[0])
    antibody_dict,antigen_dict = differ_ab_ag(pdb_file_list)
    interface_antibody_df,interface_antigen_df = ab_ag_epi_df(pdb_file_list)
    ab_ag_cor_dict = ab_ag_correspondence(interface_antibody_df,interface_antigen_df)
    interface_antigen_list = list(ab_ag_cor_dict.keys())
    for ag_amino in interface_antigen_list:
        ag_amino_df = get_amino_dataFrame(interface_antigen_df,ag_amino)
        ab_amino_list = ab_ag_cor_dict[ag_amino]
        rotamer_name = top5_amino(ag_amino[2])
        












