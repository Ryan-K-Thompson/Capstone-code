# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 10:19:14 2022

@author: patag
"""


import csv
from pandas import *
import matplotlib.pyplot as plt
import math
import numpy as np


def extract_values(File, value):
    temp = read_csv(File)
    print(temp)
    # temp = temp.dropna
    # temp=list(temp[temp!=np.isnan]) 
    NaN_present = True in (ele == 'nan' for ele in temp)
    if NaN_present == True:
        while NaN in temp:
            temp.remove(NaN)
    values_array = temp[value].tolist()
    return values_array

def XY_to_pathlength(X,Y):
    S = []
    for count,x in enumerate(X):
        if count == 0:
            S.append(0)
        else:
            s = ((X[count]-X[count-1])**2+(Y[count]-Y[count-1])**2)**0.5 + S[count - 1]
            S.append(s)
    return S

def scale_by(List, scale_factor):
    new_List = [i*scale_factor for i in List]
    return new_List

data = read_csv("laminar.csv")
# data=data.dropna()
x_laminar = data['x'].tolist()
y_laminar = data['y'].tolist()
s_laminar = XY_to_pathlength(x_laminar, y_laminar)
s_div_r_laminar = scale_by(s_laminar, 1/0.0254)
qdot_laminar = data['qdot'].tolist()
qdot_laminar = scale_by(qdot_laminar, 1/(100*100))

data = read_csv("turbulent.csv")
data=data.dropna()
x_turbulent = data['x'].tolist()
y_turbulent = data['y'].tolist()
s_turbulent = XY_to_pathlength(x_turbulent, y_turbulent)
s_div_r_turbulent = scale_by(s_turbulent, 1/0.0254)
qdot_turbulent = data['qdot'].tolist()
qdot_turbulent = scale_by(qdot_turbulent, 1/(100*100))







fig6 = plt.figure("Figure 6")
plt.plot(s_div_r_laminar,qdot_laminar, color = 'r')
plt.plot(s_div_r_turbulent,qdot_turbulent, color = 'k')
plt.grid()
plt.xlabel('S/R')
plt.ylabel('Heat transfer ($\mathregular{W/cm^{2}}$)')
#plt.title('Flight heat transfer')
plt.xlim(0,1.2)
plt.legend(["Laminar","Turbulent"])
# plt.xlim([0,1.2])
# plt.ylim(0,100)
plt.savefig('plt_comparison_qdot.png', dpi=300, bbox_inches="tight")
plt.show()


# # 

