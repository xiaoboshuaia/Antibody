# -*- coding: utf-8 -*-
import pandas as pd
from pandas import DataFrame
with open('7bwj.pdb', "r") as file:
    file = file.read().split("\n")


ATOM = []
for i in range(len(file)):
    anjisuan = []
    if file[i].startswith('ATOM'):
        anjisuan = [file[i][0:4],file[i][6:11],file[i][12:16],file[i][17:20],file[i][21],
                    file[i][22:26],file[i][30:38],file[i][38:46],file[i][46:54]
                   ]
        ATOM.append(anjisuan)

df = DataFrame(ATOM)




 




































