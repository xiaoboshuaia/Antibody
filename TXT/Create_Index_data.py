'''
Author: xiaobo 973801194@qq.com
Date: 2022-06-24 11:19:39
LastEditors: xiaobo 973801194@qq.com
LastEditTime: 2022-07-04 13:44:14
FilePath: \第二轮轮转\索引表\索引表\Antibody_optimization\Create_Index_data.py
Description: 
'''


from collections import Counter
import numpy
from pandas import  read_csv
import os
# from Parent directory import module

numpy.array()

# get txt.file from current directory
def get_txt_file(path):
    file_list = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.txt'):
                file_list.append(os.path.join(root, file))
    return file_list

#get Index_data_dictionary from input txt.file list
def get_Index_data_dictionary(txt_file_list):
    Index_data_dictionary = {}
    for i in txt_file_list:
        data = read_csv(i, delim_whitespace=True)
        for index, row in data.iterrows():
            if row['Amino_antigen'] not in Index_data_dictionary.keys():
                Index_data_dictionary[row['Amino_antigen']] = [row['Amino_antibody']]
            else:
                Index_data_dictionary[row['Amino_antigen']].append(row['Amino_antibody'])
    return Index_data_dictionary

# get top 5 amino_antibody from Index_data_dictionary
def get_top_5_amino_antibody(Index_data_dictionary):
    top_5_amino_antibody = {}
    for key, value in Index_data_dictionary.items():
        top_5_amino_antibody[key] = Counter(value).most_common(5)
    return top_5_amino_antibody

# main function
if __name__ == '__main__':# 判断是否为主程序运行
    pdb_Index_file = get_txt_file(os.getcwd())# 获取当前目录下的txt文件
    Index_data_dictionary = get_Index_data_dictionary(pdb_Index_file)# 创建索引表
    top_5_amino_antibody = get_top_5_amino_antibody(Index_data_dictionary)# 获取索引表中的前5个结果