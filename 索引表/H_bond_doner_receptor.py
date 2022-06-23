######将氢键供体与氢键受体的重原子通过列表的形式进行保存
def H_bond_do_re():
    donor = [],receptor = []
    do = ["antibody[1]['amino']=='GLY' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='ALA' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='VAL' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='LEU' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='MET' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='ILE' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='SER' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='SER' and antibody[1]['ATOM']=='OG'",
 "            or antibody[1]['amino']=='THR' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='THR' and antibody[1]['ATOM']=='OG1'",
 "            or antibody[1]['amino']=='CYS' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='CYS' and antibody[1]['ATOM']=='SG'",
 "            or antibody[1]['amino']=='ASN' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='ASN' and antibody[1]['ATOM']=='ND2'",
 "            or antibody[1]['amino']=='GLN' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='GLN' and antibody[1]['ATOM']=='NE2'",
 "            or antibody[1]['amino']=='PHE' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='TYR' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='TYR' and antibody[1]['ATOM']=='OH'",
 "            or antibody[1]['amino']=='TRP' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='TRP' and antibody[1]['ATOM']=='NE1'",
 "            or antibody[1]['amino']=='LYS' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='LYS' and antibody[1]['ATOM']=='NZ'",
 "            or antibody[1]['amino']=='ARG' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='ARG' and antibody[1]['ATOM']=='NE'",
 "            or antibody[1]['amino']=='ARG' and antibody[1]['ATOM']=='NH1'",
 "            or antibody[1]['amino']=='ARG' and antibody[1]['ATOM']=='NH2'",
 "            or antibody[1]['amino']=='HIS' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='HIS' and antibody[1]['ATOM']=='NE2'",
 "            or antibody[1]['amino']=='ASP' and antibody[1]['ATOM']=='N'",
 "            or antibody[1]['amino']=='GLU' and antibody[1]['ATOM']=='N'"]
    
    re = ["  antibody[1]['amino']=='GLY' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='ALA' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='VAL' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='LEU' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='MET' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='ILE' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='SER' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='THR' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='CYS' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='CYS' and antibody[1]['ATOM']=='SG'",
 "            or antibody[1]['amino']=='PRO' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='ASN' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='ASN' and antibody[1]['ATOM']=='OD1'",
 "            or antibody[1]['amino']=='GLN' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='GLN' and antibody[1]['ATOM']=='OE1'",
 "            or antibody[1]['amino']=='PHE' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='TYR' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='TRP' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='LYS' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='ARG' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='HIS' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='HIS' and antibody[1]['ATOM']=='ND1'",
 "            or antibody[1]['amino']=='ASP' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='ASP' and antibody[1]['ATOM']=='OD1'",
 "            or antibody[1]['amino']=='ASP' and antibody[1]['ATOM']=='OD2'",
 "            or antibody[1]['amino']=='GLU' and antibody[1]['ATOM']=='O'",
 "            or antibody[1]['amino']=='GLU' and antibody[1]['ATOM']=='OE1'",
 "            or antibody[1]['amino']=='GLU' and antibody[1]['ATOM']=='OE2'"]
     
    for i in re:
        pair = []
        pair = [i.split('\'')[3],i.split('\'')[-2]]
        receptor.append(pair)        
    for i in do:
        inner_pair = []
        inner_pair = [i.split('\'')[3],i.split('\'')[-2]]
        donor.append(inner_pair)   
    return donor,receptor


donor,receptor = H_bond_do_re()








































        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    