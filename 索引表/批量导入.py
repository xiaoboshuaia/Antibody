# -*- coding: utf-8 -
import os


###定义函数pair
def pair(file):
    COMPND = []
    for i in range(len(file)):
        B = []
        if file[i].startswith('COMPND'):
            B = [file[i]]
            COMPND.append(B)



###将保存抗原和抗体的信息分别储存在两个列表
    antigen_chain=[]
    antibody_chain=[]         
    
    for i in range(len(COMPND)):
        if 'SPIKE' in COMPND[i][0]:
            antigen_chain.append(COMPND[i+1])
        if 'CHAIN' in COMPND[i][0][11:16] and 'SPIKE'  not in COMPND[i-1][0]:
            antibody_chain.append(COMPND[i])
               
###对列表进行修改，只保留链的信息比如 antibody_chain=[A,B,C] 
###将CHAIN：后面的字符串取出 再将非字母的字符串去除，将剩下的字母（每条链的编号）保存在一个list中，方便调用     
    antigen_chain_1 = []
    antibody_chain_1 = []   
    for i in antigen_chain:
        c = i[0][17:-1].replace(' ','').replace(';','').replace(',','')
        for i in c:
            antigen_chain_1.append(i)
    for i in antibody_chain:
        c = i[0][17:-1].replace(' ','').replace(';','').replace(',','')
        for i in c:
            antibody_chain_1.append(i)    

#####将带有ATOM的行提取出来
#####将抗原和抗体的氨基酸分别存在两个list中
    ATOM = []
    for i in range(len(file)):
        if file[i].startswith('ATOM'):
                ATOM.append(file[i])
    
    antigen=[]
    antibody=[]
    for i in ATOM:
        if  i[21] in antigen_chain_1 :
            antigen.append(i)
        elif i[21]  in  antibody_chain_1:
            antibody.append(i) 
##计算疏水键，将每个氨基酸的所有属性保存到一个list中，便于下一步计算   
    ag = []
    n = 0
    for i in range(len(antigen)-1):
        if antigen[i][22:26] != antigen[i+1][22:26]:
            ag.append(antigen[n:i+1])
            n = i + 1
    ag.append(antigen[n:])        
    
    anti = []
    n = 0
    for i in range(len(antibody)-1):
        if antibody[i][22:26] != antibody[i+1][22:26]:
            anti.append(antibody[n:i+1])
            n = i + 1
    anti.append(antibody[n:])
    
    S = []
    for i in ag:
        x = []
        y = []
        z = []
        b = []
        if i[0][17:20] == 'GLY' :
                b =[i[0][21],i[0][17:20],i[0][22:26].replace(' ',''),'ATOM',
                    float(i[1][30:38].replace(' ','')),
                    float(i[1][38:46].replace(' ','')),
                    float(i[1][46:54].replace(' ',''))]
                S.append(b)
        elif i[0][17:20] != 'GLY' and len(i[4:]) > 1: 
            for i in i[4:]:
                x.append(float(i[30:38].replace(' ','')))
                y.append(float(i[38:46].replace(' ','')))
                z.append(float(i[46:54].replace(' ','')))
            b =[i[21],i[17:20],i[22:26].replace(' ',''),'ATOM',sum(x)/len(x),
                            sum(y)/len(y),
                            sum(z)/len(z),]
            S.append(b)
    
    H = []
    for i in anti:
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
        elif i[0][17:20] != 'GLY' and len(i[4:]) > 1: #将不是GLY的氨基酸取出来并且除去掉因为实验原因而没有测出来的侧链
            for i in i[4:] :            
                x.append(float(i[30:38].replace(' ','')))
                y.append(float(i[38:46].replace(' ','')))
                z.append(float(i[46:54].replace(' ','')))
            b =[i[21],i[17:20],i[22:26].replace(' ',''),'ATOM',sum(x)/len(x),
                            sum(y)/len(y),
                            sum(z)/len(z),]
            H.append(b)
    
    hydrophobic_amino_acid = ['MET','ALA','VAL','LEU','ILE','PHE','PRO','TRP']
    hydrophobic_bond_pair=[]
    for i in H:
        for u in S:
            dis = ((float(i[4]) - float(u[4]))**2 + 
                        (float(i[5]) - float(u[5]))**2 +
                        (float(i[6]) - float(u[6]))**2)**0.5
            b = []
            if dis <= 6.5 and dis > 2.5 and i[1] in hydrophobic_amino_acid and u[1] in hydrophobic_amino_acid :
            
                b = [i[0],i[1],i[2],'NA',u[0],u[1],u[2],'NA',dis]
                hydrophobic_bond_pair.append(b)



###氢键
    pair_big = []
    pair_big2 = []
    for i in anti:
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
                for i in anti:
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
    return complexed





























with open("pdb3bgf.txt", "r") as file:
    file = file.read().split("\n")
    

COMPND = []
for i in range(len(file)):
    B = []
    if file[i].startswith('COMPND'):
        B = [file[i]]
        COMPND.append(B)



###将保存抗原和抗体的信息分别储存在两个列表
antigen_chain=[]
antibody_chain=[]         

for i in range(len(COMPND)):
    if 'SPIKE' in COMPND[i][0]:
        antigen_chain.append(COMPND[i+1])
    if 'CHAIN' in COMPND[i][0][11:16] and 'SPIKE'  not in COMPND[i-1][0]:
        antibody_chain.append(COMPND[i])
           
###对列表进行修改，只保留链的信息比如 antibody_chain=[A,B,C] 
###将CHAIN：后面的字符串取出 再将非字母的字符串去除，将剩下的字母（每条链的编号）保存在一个list中，方便调用     
antigen_chain_1 = []
antibody_chain_1 = []   
for i in antigen_chain:
    c = i[0][17:-1].replace(' ','').replace(';','').replace(',','')
    for i in c:
        antigen_chain_1.append(i)
for i in antibody_chain:
    c = i[0][17:-1].replace(' ','').replace(';','').replace(',','')
    for i in c:
        antibody_chain_1.append(i)    

#####将带有ATOM的行提取出来
#####将抗原和抗体的氨基酸分别存在两个list中
ATOM = []
for i in range(len(file)):
    if file[i].startswith('ATOM'):
            ATOM.append(file[i])

antigen=[]
antibody=[]
for i in ATOM:
    if  i[21] in antigen_chain_1 :
        antigen.append(i)
    elif i[21]  in  antibody_chain_1:
        antibody.append(i) 
##计算疏水键，将每个氨基酸的所有属性保存到一个list中，便于下一步计算   
ag = []
n = 0
for i in range(len(antigen)):
    if antigen[i][22:26] != antigen[i+1][22:26]:
        ag.append(antigen[n:i+1])
        n = i + 1
ag.append(antigen[n:])        

anti = []
n = 0
for i in range(len(antibody)):
    if antibody[i][22:26] != antibody[i+1][22:26]:
        anti.append(antibody[n:i+1])
        n = i + 1
anti.append(antibody[n:])

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

H = []
for i in anti:
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

hydrophobic_amino_acid = ['MET','ALA','VAL','LEU','ILE','PHE','PRO','TRP']
hydrophobic_bond_pair=[]
for i in H:
    for u in S:
        dis = ((float(i[4]) - float(u[4]))**2 + 
                    (float(i[5]) - float(u[5]))**2 +
                    (float(i[6]) - float(u[6]))**2)**0.5
        b = []
        if dis <= 6.5 and dis > 2.5 and i[1] in hydrophobic_amino_acid and u[1] in hydrophobic_amino_acid :
        
            b = [i[0],i[1],i[2],'NA',u[0],u[1],u[2],'NA',dis]
            hydrophobic_bond_pair.append(b)



###氢键
pair_big = []
pair_big2 = []
for i in anti:
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
            for i in anti:
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


#####批量导入
#import glob as gb
#def search():
    #pdb = gb.glob(r'C:\Users\97380\Desktop\第二轮轮转\索引表\索引表\PDB\*.pdb')
     
pdb = os.listdir(r'C:\Users\97380\Desktop\第二轮轮转\索引表\索引表\PDB')
for i in pdb:
    with open(i, "r") as file:

        file = file.read().split("\n")
        complexed = pair(file)
    output = open(i+'.txt','w',encoding='gbk')
    output.write('CLAIN,Antibody_amino_acid,resSeq,ATOM,CLAIN,Antigen_amino_acid,resSeq,ATOM,distance\n')
    for row in complexed:
    	rowtxt = '{},{},{},{},{},{},{},{},{}'.format(row[0],row[1],row[2],row[3],row[4],
                                                  row[5],row[6],row[7],round(row[8],2))
    	output.write(rowtxt)
    	output.write('\n')
    output.close() 





















####将计算相互作用对定义为函数pair

    










