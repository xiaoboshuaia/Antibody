'''
Author: xiaobo 973801194@qq.com
Date: 2022-06-24 11:19:39
LastEditors: xiaobo 973801194@qq.com
LastEditTime: 2022-06-24 13:56:13
FilePath: \第二轮轮转\索引表\索引表\Antibody_optimization\Create_Index_data.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''



from pandas import DataFrame, read_csv
import os
# from Parent directory import module
import sys
sys.path.append('..')


# get txt.file from current directory
def get_txt_file(path):
    file_list = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.txt'):
                file_list.append(os.path.join(root, file))
    return file_list


# main function
if __name__ == '__main__':# 判断是否为主程序运行
    pdb_Index_file = get_txt_file(os.getcwd())# 获取当前目录下的txt文件
    Index_date = {}# 创建索引表
    for i in pdb_Index_file:# 循环遍历txt文件
        data = read_csv(i, delim_whitespace=True)# 读取txt文件
        for index, row in data.iterrows():# 循环遍历txt文件中的每一行
            if row['Amino_antibody'] not in Index_date:# 判断索引表中是否有该行的索引
                Index_date[row['Amino_antigen']] = [row['Amino_antibody']]# 如果没有则添加该行的索引
            else:
                Index_date[row['Amino_antigen']].append(row['Amino_antibody'])# 如果有则添加该行的值
