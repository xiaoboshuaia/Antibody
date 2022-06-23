




##将带有ATOM的行从整个PDB文件中提取出来
import pandas as pd
###   PDBread函数，将PDB文件中所有的信息，通过所在字符串中的位置进行提取


a = read






with open("pdb3bgf.txt", "r") as file:
    file = file.read().split("\n")
    
i = 0
while i < len(file):
    
    if file[i].startswith("ATOM") == False:
        del file[i]
        i -= 1
    
    i += 1
###将得到的ATOM通过属于不同的个链来进行区分
###L,S,H为一组 B,A,C为一组
   
antigen=[]
antibodyL=[]
antibodyH=[]
for i in file:
    
    if  i[21] == 'S' :
        antigen.append(i)
    
    elif i[21] == 'L':
        antibodyL.append(i) 
    elif i[21] == 'H':
        antibodyH.append(i)



###将抗原抗体的结构分别存为txt文件
str = '\n'        
f1= open('antigen.txt','w')
f1.write(str.join(antigen))
f1.close
f2 = open('antibodyL.txt','w')
f2.write(str.join(antibodyL))
f2.close
f3 = open('antibodyH.txt','w')
f3.write(str.join(antibodyH))
f3.close


    

        

#将列表中相同的元素形成一个新的列表，生成一个二维列表
a=[]
rt = []
n = 0
for i in range(len(a)-1):
        if a[i] != a[i+1]:
                rt.append(a[n:i+1])
                n = i+1
rt.append(a[n:])


rt = []
n = 0
i = 0
while i < len(ori_list)-1:
    if ori_list[i] != ori_list[i+1]:
        rt.append(ori_list[n:i+1])
        n = i+1
    i += 1

rt.append(ori_list[n:])









#####将antibodyL中的每个氨基酸残基的数据，分别储存在一个列表中，最终组成为一个大的列表
anL = []
n = 0
for i in range(len(antibodyL)-1):
    if antibodyL[i][22:26] != antibodyL[i+1][22:26]:
        anL.append(antibodyL[n:i+1])
        n = i + 1
anL.append(antibodyL[n:])

anH = []
n = 0
for i in range(len(antibodyH)-1):
    if antibodyH[i][22:26] != antibodyH[i+1][22:26]:
        anH.append(antibodyH[n:i+1])
        n = i + 1
anH.append(antibodyH[n:])

ag = []
n = 0
for i in range(len(antigen)-1):
    if antigen[i][22:26] != antigen[i+1][22:26]:
        ag.append(antigen[n:i+1])
        n = i + 1
ag.append(antigen[n:])

    


#####将抗体的每个氨基酸的侧链的中心位置的三位坐标，以及氨基酸的名称，储存在列表中
####首先对于氨基酸是否是GLY进行判断，如果是则使用Cα的三维坐标,如果不是则将除了主链上的4个原子除去，
##计算侧链上的原子的坐标的和再除以原子的数量。
L = []
for i in anL:
    x = []
    y = []
    z = []
    b = []
    if i[0][17:20] == 'GLY':
            b =[i[0][21],i[0][17:20],i[0][22:26].replace(' ',''),'ATOM',float(i[1][30:38].replace(' ','')),
                
                float(i[1][38:46].replace(' ','')),
                float(i[1][46:54].replace(' ',''))]
            L.append(b)
    else: 
        for i in i[4:]:
            x.append(float(i[30:38].replace(' ','')))
            y.append(float(i[38:46].replace(' ','')))
            z.append(float(i[46:54].replace(' ','')))
        b =[i[21],i[17:20],i[22:26].replace(' ',''),'ATOM',sum(x)/len(x),
                        sum(y)/len(y),
                        sum(z)/len(z),]
        L.append(b)

H = []
for i in anH:
    x = []
    y = []
    z = []
    b = []
    if i[0][17:20] == 'GLY':
            b =[i[0][21],i[0][17:20],i[0][22:26].replace(' ',''),'ATOM',
                float(i[1][30:38].replace(' ','')),
                float(i[1][38:46].replace(' ','')),
                float(i[1][46:54].replace(' ',''))]
            H.append(b)
    else: 
        for i in i[4:]:
            x.append(float(i[30:38].replace(' ','')))
            y.append(float(i[38:46].replace(' ','')))
            z.append(float(i[46:54].replace(' ','')))
        b =[i[21],i[17:20],i[22:26].replace(' ',''),'ATOM',sum(x)/len(x),
                        sum(y)/len(y),
                        sum(z)/len(z),]
        H.append(b)
S = []
for i in ag:
    x = []
    y = []
    z = []
    b = []
    if i[0][17:20] == 'GLY':
            b =[i[0][21],i[0][17:20],i[0][22:26].replace(' ',''),'ATOM',
                float(i[1][30:38].replace(' ','')),
                float(i[1][38:46].replace(' ','')),
                float(i[1][46:54].replace(' ',''))]
            S.append(b)
    else: 
        for i in i[4:]:
            x.append(float(i[30:38].replace(' ','')))
            y.append(float(i[38:46].replace(' ','')))
            z.append(float(i[46:54].replace(' ','')))
        b =[i[21],i[17:20],i[22:26].replace(' ',''),'ATOM',sum(x)/len(x),
                        sum(y)/len(y),
                        sum(z)/len(z),]
        S.append(b)


        
            
#####将残基数据导出为txt格式
c=[]

output = open('HH.txt','w',encoding='gbk')
output.write('name,X,Y,Z\n')
for row in c:
	rowtxt = '{},{},{},{}'.format(row[0],row[1],row[2],row[3])
	output.write(rowtxt)
	output.write('\n')
output.close()   




####计算抗体上每个氨基酸到抗原上的距离，
####将三条链的变量分别储存在S,H,L中 
###通过距离进行筛选，将距离小于等于6.5A，且满足成键氨基酸都为疏水氨基酸
hydrophobic_bond_pair=[]
hydrophobic_amino_acid = ['MET','ALA','VAL','LEU','ILE','PHE','PRO','TRP']
body = H + L
for i in body:
    for u in S:
        dis = ((float(i[4]) - float(u[4]))**2 + 
                    (float(i[5]) - float(u[5]))**2 +
                    (float(i[6]) - float(u[6]))**2)**0.5
        b = []
        if dis <= 6.5 and dis > 2.5 and i[1] in hydrophobic_amino_acid and u[1] in hydrophobic_amino_acid :
        
            b = [i[0],i[1],i[2],'NA',u[0],u[1],u[2],'NA',dis]
            hydrophobic_bond_pair.append(b)
        
            
        
#####将每个氨基酸的重原子都分为氢键供体和氢键受体，计算抗原中氢键供体的重原子，与抗原中氢键受体的重原子的距离,
#######计算抗体中氢键受体与抗原中氢键供体的距离
#####将距离小于3.5A，但大于2.5A的氢键对进行保留
###找出每个氨基酸中氢键的供体和受体

        ### GLY
        ##氢键供体
        if donor[17:20] == 'GLY' and donor[13:15] == 'N ':
        ##氢键受体
        if receptor[17:20] == 'GLY' and receptor[13:15] == 'O ':
        ### ALA
        ##氢键供体
        if donor[17:20] == 'ALA' and donor[13:15] == 'N ':
        ##氢键受体
        if receptor[17:20] == 'ALA' and receptor[13:15] == 'O ':    
        ###VAL
        ##氢键供体
        if donor[17:20] == 'VAL' and donor[13:15] == 'N ':
        ##氢键受体
        if receptor[17:20] == 'VAL' and receptor[13:15] == 'O ':
        ###LEU
        ##氢键供体
        if donor[17:20] == 'LEU' and donor[13:15] == 'N ':
        ##氢键受体
        if receptor[17:20] == 'LEU' and receptor[13:15] == 'O ':
        ###MET
        ##氢键供体
        if donor[17:20] == 'MET' and donor[13:15] == 'N ':
        ##氢键受体
        if receptor[17:20] == 'MET' and receptor[13:15] == 'O ':  
        ###ILE
        ##氢键供体
        if donor[17:20] == 'ILE' and donor[13:15] == 'N ':
        ##氢键受体
        if receptor[17:20] == 'ILE' and receptor[13:15] == 'O ':
        ###SER
        ##氢键供体
        if donor[17:20] == 'SER' and donor[13:15] == 'N ':
        if donor[17:20] == 'SER' and donor[13:15] == 'OG':
        ##氢键受体
        if receptor[17:20] == 'SER' and receptor[13:15] == 'O ':
        ###THR
        ##氢键供体
        if donor[17:20] == 'THR' and donor[13:15] == 'N ':
        if donor[17:20] == 'THR' and donor[13:16] == 'OG1':
        ##氢键受体
        if receptor[17:20] == 'THR' and receptor[13:15] == 'O ': 
        ###CYS
        ##氢键供体
        if donor[17:20] == 'CYS' and donor[13:15] == 'N ':
        if donor[17:20] == 'CYS' and donor[13:15] == 'SG':
        ##氢键受体
        if receptor[17:20] == 'CYS' and receptor[13:15] == 'O ':
        if receptor[17:20] == 'CYS' and receptor[13:15] == 'SG':
        ### PRO
        ###没有氢键的供体
        ##氢键受体
        if receptor[17:20] == 'PRO' and receptor[13:15] == 'O ':
        ###ASN
        ##氢键供体
        if donor[17:20] == 'ASN' and donor[13:15] == 'N ':
        if donor[17:20] == 'ASN' and donor[13:16] == 'ND2':
        ##氢键受体
        if receptor[17:20] == 'ASN' and receptor[13:15] == 'O ':
        if receptor[17:20] == 'ASN' and receptor[13:16] == 'OD1':
        ###GLN
        ##氢键供体
        if donor[17:20] == 'GLN' and donor[13:15] == 'N ':
        if donor[17:20] == 'GLN' and donor[13:16] == 'NE2':
        ##氢键受体
        if receptor[17:20] == 'GLN' and receptor[13:15] == 'O ':
        if receptor[17:20] == 'GLN' and receptor[13:16] == 'OE1':
        ### PHE
        ##氢键供体
        if donor[17:20] == 'PHE' and donor[13:15] == 'N ':
        ##氢键受体
        if receptor[17:20] == 'PHE' and receptor[13:15] == 'O ':
        ### TYR
        ##氢键供体
        if donor[17:20] == 'TYR' and donor[13:15] == 'N ':
        if donor[17:20] == 'TYR' and donor[13:15] == 'OH':
        ##氢键受体
        if receptor[17:20] == 'TYR' and receptor[13:15] == 'O ':
        ### TRP
        ##氢键供体
        if donor[17:20] == 'TRP' and donor[13:15] == 'N ':
        if donor[17:20] == 'TRP' and donor[13:16] == 'NE1':
        ##氢键受体
        if receptor[17:20] == 'TRP' and receptor[13:15] == 'O ':
        ### LYS
        ##氢键供体
        if donor[17:20] == 'LYS' and donor[13:15] == 'N ':
        if donor[17:20] == 'LYS' and donor[13:15] == 'NZ':
        ##氢键受体
        if receptor[17:20] == 'LYS' and receptor[13:15] == 'O ':
        ### ARG
        ##氢键供体
        if donor[17:20] == 'ARG' and donor[13:15] == 'N ':
        if donor[17:20] == 'ARG' and donor[13:15] == 'NE':
        if donor[17:20] == 'ARG' and donor[13:16] == 'NH1':
        if donor[17:20] == 'ARG' and donor[13:16] == 'NH2':
        ##氢键受体
        if receptor[17:20] == 'ARG' and receptor[13:15] == 'O ':
        ### HIS
        ##氢键供体
        if donor[17:20] == 'HIS' and donor[13:15] == 'N ':
        if donor[17:20] == 'HIS' and donor[13:16] == 'NE2':
        ##氢键受体
        if receptor[17:20] == 'HIS' and receptor[13:15] == 'O ':
        if receptor[17:20] == 'HIS' and receptor[13:16] == 'ND1':
        ### ASP
        ##氢键供体
        if donor[17:20] == 'ASP' and donor[13:15] == 'N ':
        ##氢键受体
        if receptor[17:20] == 'ASP' and receptor[13:15] == 'O ':
        if receptor[17:20] == 'ASP' and receptor[13:16] == 'OD1':
        if receptor[17:20] == 'ASP' and receptor[13:16] == 'OD2':
        ### GLU
        ##氢键供体
        if donor[17:20] == 'GLU' and donor[13:15] == 'N ':
        ##氢键受体
        if receptor[17:20] == 'GLU' and receptor[13:15] == 'O ':
        if receptor[17:20] == 'GLU' and receptor[13:16] == 'OE1':
        if receptor[17:20] == 'GLU' and receptor[13:16] == 'OE2':

    dis = ((float(i[4]) - float(u[4]))**2 + 
                    (float(i[5]) - float(u[5]))**2 +
                    (float(i[6]) - float(u[6]))**2)**0.5
        b = []
        if dis <= 6.5 and dis > 2.5 and i[1] in hydrophobic_amino_acid and u[1] in hydrophobic_amino_acid :
        
            b = [i[0],i[1],i[2],u[0],u[1],u[2],dis]        



    
#####抗体上的氢键供体和抗原上的氢键受体   
antizong = anH + anL
pair_big = []
pair_big2 = []
for i in antizong:
    for donor in i:
        if (donor[17:20] == 'GLY' and donor[13:15] == 'N '
        or donor[17:20] == 'ALA' and donor[13:15] == 'N '
        or donor[17:20] == 'VAL' and donor[13:15] == 'N '
        or donor[17:20] == 'LEU' and donor[13:15] == 'N '
        or donor[17:20] == 'MET' and donor[13:15] == 'N '
        or donor[17:20] == 'ILE' and donor[13:15] == 'N '
        or donor[17:20] == 'SER' and donor[13:15] == 'N '
        or donor[17:20] == 'SER' and donor[13:15] == 'OG'
        or donor[17:20] == 'THR' and donor[13:15] == 'N '
        or donor[17:20] == 'THR' and donor[13:16] == 'OG1'
        or donor[17:20] == 'CYS' and donor[13:15] == 'N '
        or donor[17:20] == 'CYS' and donor[13:15] == 'SG' 
        or donor[17:20] == 'ASN' and donor[13:15] == 'N '
        or donor[17:20] == 'ASN' and donor[13:16] == 'ND2'
        or donor[17:20] == 'GLN' and donor[13:15] == 'N '
        or donor[17:20] == 'GLN' and donor[13:16] == 'NE2'
        or donor[17:20] == 'PHE' and donor[13:15] == 'N '
        or donor[17:20] == 'TYR' and donor[13:15] == 'N '
        or donor[17:20] == 'TYR' and donor[13:15] == 'OH'
        or donor[17:20] == 'TRP' and donor[13:15] == 'N '
        or donor[17:20] == 'TRP' and donor[13:16] == 'NE1'
        or donor[17:20] == 'LYS' and donor[13:15] == 'N '
        or donor[17:20] == 'LYS' and donor[13:15] == 'NZ'
        or donor[17:20] == 'ARG' and donor[13:15] == 'N '
        or donor[17:20] == 'ARG' and donor[13:15] == 'NE'
        or donor[17:20] == 'ARG' and donor[13:16] == 'NH1'
        or donor[17:20] == 'ARG' and donor[13:16] == 'NH2'
        or donor[17:20] == 'HIS' and donor[13:15] == 'N '
        or donor[17:20] == 'HIS' and donor[13:16] == 'NE2'
        or donor[17:20] == 'ASP' and donor[13:15] == 'N '
        or donor[17:20] == 'GLU' and donor[13:15] == 'N '):
            for i in ag:
                for receptor in i :
                    dis = ((float(donor[30:38].replace(' ','')) - float(receptor[30:38].replace(' ','')))**2 + 
                    (float(donor[38:46].replace(' ','')) - float(receptor[38:46].replace(' ','')))**2 +
                    (float(donor[46:54].replace(' ','')) - float(receptor[46:54].replace(' ','')))**2)**0.5
                    pair_xiao = []
                    if (receptor[17:20] == 'GLY' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5 
                    or receptor[17:20] == 'ALA' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'VAL' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'LEU' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'MET' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ILE' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'SER' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'THR' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'CYS' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5                                                     
                    or receptor[17:20] == 'CYS' and receptor[13:15] == 'SG' and dis <= 3.5 and dis > 2.5                             
                    or receptor[17:20] == 'PRO' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASN' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASN' and receptor[13:16] == 'OD1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLN' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLN' and receptor[13:16] == 'OE1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'PHE' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'TYR' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'TRP' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'LYS' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ARG' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'HIS' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'HIS' and receptor[13:16] == 'ND1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASP' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASP' and receptor[13:16] == 'OD1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASP' and receptor[13:16] == 'OD2'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLU' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLU' and receptor[13:16] == 'OE1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLU' and receptor[13:16] == 'OE2'and dis <= 3.5 and dis > 2.5):
                        pair_xiao = [donor[21],donor[17:20],donor[22:26].replace(' ', ''),donor[13:16].replace(' ',''),
                        receptor[21],receptor[17:20],receptor[22:26].replace(' ', ''),receptor[13:16].replace(' ',''),dis]
                        pair_big.append(pair_xiao)





####抗原上的氢键供体和抗体上的氢键受体                    
for i in ag:
    for donor in i:
        if (donor[17:20] == 'GLY' and donor[13:15] == 'N '
        or donor[17:20] == 'ALA' and donor[13:15] == 'N '
        or donor[17:20] == 'VAL' and donor[13:15] == 'N '
        or donor[17:20] == 'LEU' and donor[13:15] == 'N '
        or donor[17:20] == 'MET' and donor[13:15] == 'N '
        or donor[17:20] == 'ILE' and donor[13:15] == 'N '
        or donor[17:20] == 'SER' and donor[13:15] == 'N '
        or donor[17:20] == 'SER' and donor[13:15] == 'OG'
        or donor[17:20] == 'THR' and donor[13:15] == 'N '
        or donor[17:20] == 'THR' and donor[13:16] == 'OG1'
        or donor[17:20] == 'CYS' and donor[13:15] == 'N '
        or donor[17:20] == 'CYS' and donor[13:15] == 'SG' 
        or donor[17:20] == 'ASN' and donor[13:15] == 'N '
        or donor[17:20] == 'ASN' and donor[13:16] == 'ND2'
        or donor[17:20] == 'GLN' and donor[13:15] == 'N '
        or donor[17:20] == 'GLN' and donor[13:16] == 'NE2'
        or donor[17:20] == 'PHE' and donor[13:15] == 'N '
        or donor[17:20] == 'TYR' and donor[13:15] == 'N '
        or donor[17:20] == 'TYR' and donor[13:15] == 'OH'
        or donor[17:20] == 'TRP' and donor[13:15] == 'N '
        or donor[17:20] == 'TRP' and donor[13:16] == 'NE1'
        or donor[17:20] == 'LYS' and donor[13:15] == 'N '
        or donor[17:20] == 'LYS' and donor[13:15] == 'NZ'
        or donor[17:20] == 'ARG' and donor[13:15] == 'N '
        or donor[17:20] == 'ARG' and donor[13:15] == 'NE'
        or donor[17:20] == 'ARG' and donor[13:16] == 'NH1'
        or donor[17:20] == 'ARG' and donor[13:16] == 'NH2'
        or donor[17:20] == 'HIS' and donor[13:15] == 'N '
        or donor[17:20] == 'HIS' and donor[13:16] == 'NE2'
        or donor[17:20] == 'ASP' and donor[13:15] == 'N '
        or donor[17:20] == 'GLU' and donor[13:15] == 'N '):
            for i in antizong:
                for receptor in i :
                    dis = ((float(donor[30:38].replace(' ','')) - float(receptor[30:38].replace(' ','')))**2 + 
                    (float(donor[38:46].replace(' ','')) - float(receptor[38:46].replace(' ','')))**2 +
                    (float(donor[46:54].replace(' ','')) - float(receptor[46:54].replace(' ','')))**2)**0.5
                    pair_xiao = []
                    if (receptor[17:20] == 'GLY' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5 
                    or receptor[17:20] == 'ALA' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'VAL' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'LEU' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'MET' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ILE' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'SER' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'THR' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'CYS' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5                                                     
                    or receptor[17:20] == 'CYS' and receptor[13:15] == 'SG' and dis <= 3.5 and dis > 2.5                             
                    or receptor[17:20] == 'PRO' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASN' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASN' and receptor[13:16] == 'OD1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLN' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLN' and receptor[13:16] == 'OE1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'PHE' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'TYR' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'TRP' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'LYS' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ARG' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'HIS' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'HIS' and receptor[13:16] == 'ND1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASP' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASP' and receptor[13:16] == 'OD1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'ASP' and receptor[13:16] == 'OD2'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLU' and receptor[13:15] == 'O ' and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLU' and receptor[13:16] == 'OE1'and dis <= 3.5 and dis > 2.5
                    or receptor[17:20] == 'GLU' and receptor[13:16] == 'OE2'and dis <= 3.5 and dis > 2.5):
                        pair_xiao = [receptor[21],receptor[17:20],receptor[22:26].replace(' ', ''),receptor[13:16].replace(' ',''),
                                     donor[21],donor[17:20],donor[22:26].replace(' ', ''),donor[13:16].replace(' ',''),dis]
                        pair_big2.append(pair_xiao) 
                        
                        
H_bond_pair = pair_big + pair_big2                        
complexed =  H_bond_pair + hydrophobic_bond_pair
                       
output = open('3bgf.txt','w',encoding='gbk')
output.write('CLAIN,Antibody_amino_acid,resSeq,ATOM,CLAIN,Antigen_amino_acid,resSeq,ATOM,distance\n')
for row in complexed:
	rowtxt = '{},{},{},{},{},{},{},{},{}'.format(row[0],row[1],row[2],row[3],row[4],
                                              row[5],row[6],row[7],round(row[8],2))
	output.write(rowtxt)
	output.write('\n')
output.close()                        
                        
                        

                            
                        
        
    
 
    
 
    
chain = ['A','B','C']
A = ('A','B','C')
for i in A:
    if i in chain:
        print(i)
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
        
        
        



    

            
















