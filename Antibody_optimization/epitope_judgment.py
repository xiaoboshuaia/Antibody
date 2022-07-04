'''
Author: xiaobo 973801194@qq.com
Date: 2022-01-21 17:01:21
LastEditors: xiaobo 973801194@qq.com
LastEditTime: 2022-07-04 18:10:10
FilePath: \Project\Antibody\Index\Antibody_optimization\epitope_judgment.py
Description: 
'''
# -*- coding: utf-8 -*-



from Wheels import dict_to_df, differ_ab_ag, distance

"""
将抗原和抗体的表位附近的氨基酸找出来

有助于减少计算量
"""

# get CA from antigen/antibody


def get_CA(ag_ab_df):
    CA = ag_ab_df.loc[ag_ab_df['ATOM_NAME'] == ' CA ']
    return CA

# get CA's coordinate range,'x','y','z'


def range_CA(CA, x_y_z):
    range_CA = [min(CA['coordinate_'+x_y_z].tolist()),
                max(CA['coordinate_'+x_y_z].tolist())]
    return range_CA

# H-band CA range 19A


def range_H_CA(CA, x_y_z):
    CA_min = range_CA(CA, x_y_z)[0]
    CA_max = range_CA(CA, x_y_z)[1]
    return CA_min - 19, CA_max + 19

# H-band CA range 15A


#def range_H_CA(CA, x_y_z):
    CA_min = range_CA(CA, x_y_z)[0]
    CA_max = range_CA(CA, x_y_z)[1]
    return CA_min - 15, CA_max + 15

# CA from antibody or antigen probably have H_bond,
# CA could be antigen or antibody,the first CA is used to screening


def ab_ag_CA(first_CA, second_CA):
    x_min, x_max = range_H_CA(second_CA, 'x')
    y_min, y_max = range_H_CA(second_CA, 'y')
    z_min, z_max = range_H_CA(second_CA, 'z')
    first_CA = first_CA.copy()
    x = first_CA[(first_CA.coordinate_x > x_min) &
                 (first_CA.coordinate_x < x_max) & 
                 (first_CA.coordinate_y > y_min) & 
                 (first_CA.coordinate_y < y_max) & 
                 (first_CA.coordinate_z > z_min) &
                 (first_CA.coordinate_z < z_max)]
    return x

# delete the repeat element in dict values


def de_rep_val(dict):
    for keys in dict:
        dict[keys] = set(dict[keys])
    return dict

# get antibody_CA and antigen_CA


def ab_ag_CA_2(antibody, antigen):
    ab_CA = ab_ag_CA(get_CA(antibody), get_CA(antigen))
    ag_CA = ab_ag_CA(get_CA(antigen), ab_CA)
    return ab_CA, ag_CA

# get re_seq from potential H-bond amino,save with dict {'chain':['re_seq']},input df.CA,output dict


def ab_ag_req_dic(antibody, antigen):
    ab_CA, ag_CA = ab_ag_CA_2(antibody, antigen)
    ab_re_seq = {}
    ag_re_seq = {}
    for index_ab, ab in ab_CA.iterrows():
        for index_ag, ag in ag_CA.iterrows():
            if distance(ab, ag) < 19:
                if ab['chain'] not in ab_re_seq.keys():
                    ab_re_seq[ab['chain']] = [ab['re_seq']]
                ab_re_seq[ab['chain']].append(ab['re_seq'])
                if ag['chain'] not in ag_re_seq.keys():
                    ag_re_seq[ag['chain']] = [ag['re_seq']]
                ag_re_seq[ag['chain']].append(ag['re_seq'])
    return de_rep_val(ab_re_seq), de_rep_val(ag_re_seq)

#get antigen and antibody amino correspondence input antibody,antigen_df output dict
def ab_ag_correspondence(antibody,antigen):
    ab_CA, ag_CA = ab_ag_CA_2(antibody, antigen)
    ab_ag_correspondence = {}
    for i in range(len(ag_CA)):
        for j in range(len(ab_CA)):
            if distance(ag_CA.iloc[i], ab_CA.iloc[j]) < 19:
                ag = ag_CA.iloc[i]
                ab = ab_CA.iloc[j]
                if (ag['chain'],ag['re_seq'],ag['RESIDUE_NAME']) not in ab_ag_correspondence.keys():
                    ab_ag_correspondence.update({(ag['chain'],ag['re_seq'],ag['RESIDUE_NAME']):[[ab['chain'],ab['re_seq'],ab['RESIDUE_NAME']]]})
                else:
                    ab_ag_correspondence[(ag['chain'],ag['re_seq'],ag['RESIDUE_NAME'])].append([ab['chain'],ab['re_seq'],ab['RESIDUE_NAME']])
    
    
    return ab_ag_correspondence



# create a dict to save epitope amino,from a big dict ,input antibody or antigen dict and a dict to save amino re_seq information
# format ab_ag_re_seq{'chain':['re_seq'],'chain':[re_seq]}
# format ab_ag_dict{'chain':{re_seq:[[ATOM,etc.]]}}


def dict_epito_amino(ab_ag_dict, ab_ag_re_seq):
    epito_amino_dict = {}
    for chain in ab_ag_re_seq.keys():
        epito_amino_dict[chain] = {}
        for re_seq in ab_ag_re_seq[chain]:
            epito_amino_dict[chain][re_seq] = ab_ag_dict[chain][re_seq]
    return epito_amino_dict

# get epitope amino atom dataframe,input pdb_file_df output antibody_dataframe and antigen_dataframe


def ab_ag_epi_df(file):
    ab_re_seq, ag_re_seq = ab_ag_req_dic(dict_to_df(
        differ_ab_ag(file)[0]), dict_to_df(differ_ab_ag(file)[1]))
    ab_epi_ami_df = dict_to_df(dict_epito_amino(
        differ_ab_ag(file)[0], ab_re_seq))
    ag_epi_ami_df = dict_to_df(dict_epito_amino(
        differ_ab_ag(file)[1], ag_re_seq))
    return ab_epi_ami_df, ag_epi_ami_df
