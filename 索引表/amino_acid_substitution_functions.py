# -*- coding: utf-8 -
####氨基酸取代
###将一种氨基酸的所有 Rotamer 与抗体上的任意氨基酸进行替换
###将 Rotamer库中的氨基酸的三维坐标，进行旋转以及平移，即可完成
import numpy as np
import scipy.linalg as linalg




###两向量弧度的计算
def vector_radian(revolving,target):
    ##创建两个向量
    x = np.array(revolving)
    y = np.array(target)
    ###计算两个向量的模
    m_x = np.sqrt(x.dot(x))
    m_y = np.sqrt(y.dot(y))
    ###计算两个向量的点乘
    dian = x.dot(y)
    ###计算两个向量之间的夹角的余弦值
    cos_ = dian/(m_x*m_y)
    ###将余弦值转化为弧度
    Rotation_radian = np.arccos(cos_)
    
    return Rotation_radian



###平面的法向量及旋转轴
def normal_vector(revolving,target):#a,b为两向量，两向量叉乘
    return np.cross(revolving, target)
    

###旋转矩阵
def rotatio_matrix(revolving,target):
    axis = normal_vector(revolving,target)
    radian = vector_radian(revolving,target)
    return linalg.expm(np.cross(np.eye(3), axis / linalg.norm(axis)*radian ))


###将原始坐标通过旋转矩阵进行旋转后得到的结果
def new_coor(old_coor, revolving, target):
    new_coor = np.dot(old_coor, rotatio_matrix(revolving,target))
    return new_coor


###平移坐标
###平移矩阵



###两点转换中的dx，dy，dz
#def dxdydz():
    
    #return dx,dy,dz
#def translation_matrix(revolving,target):
    
###平移变换
#def translation_coor():
    #return 







###将dataframe中x，y，z的数据转换成一个list
#def dataframe_to_list(dataframe_row):
   # return 



















