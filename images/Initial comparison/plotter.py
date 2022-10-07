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

x_flight = extract_values("wall_flight.csv", 'X')
y_flight = extract_values("wall_flight.csv", 'Y')
s_flight = XY_to_pathlength(x_flight, y_flight)
s_div_r_flight = scale_by(s_flight, 1/1.35)
qdot_flight = scale_by(extract_values("wall_flight.csv", 'Qdot'),1/(10**4))

x_subsonic = extract_values("wall_subsonic.csv", 'X')
y_subsonic = extract_values("wall_subsonic.csv", 'Y')
s_subsonic = XY_to_pathlength(x_subsonic, y_subsonic)
s_div_r_subsonic = scale_by(s_subsonic, 1/0.0254)
qdot_subsonic = scale_by(extract_values("wall_subsonic.csv", 'Qdot'),1/(10**4))



fig1 = plt.figure("Figure 1")
plt.plot(s_div_r_flight,qdot_flight, color = 'k')
plt.grid()
plt.xlabel('S/R')
plt.ylabel('Heat transfer ($\mathregular{W/cm^{2}}$)')
plt.title('Flight heat transfer')

plt.savefig('plt_Flight.png', dpi=300)
plt.show()

fig2 = plt.figure("Figure 2")
plt.plot(s_div_r_subsonic,qdot_subsonic, color = 'k')
plt.grid()
plt.xlabel('S/R')
plt.ylabel('Heat transfer ($\mathregular{W/cm^{2}}$)')
plt.title('Plasma wind tunnel heat transfer')

plt.savefig('plt_subsonic.png', dpi=300)
plt.show()

fig3 = plt.figure("Figure 3")
plt.plot(s_div_r_flight,qdot_flight, color = 'r')
plt.plot(s_div_r_subsonic,qdot_subsonic, color = 'k')
plt.grid()
plt.xlabel('S/R')
plt.ylabel('Heat transfer ($\mathregular{W/cm^{2}}$)')
# plt.title('Heat transfer comparison')
plt.xlim([0, 1.5])
plt.legend(["Flight ($\mathregular{M_{\infty}}$=16.2, 6.1 km)", "Plasma wind tunnel"])

plt.savefig('plt_subsonic_flight_comparison.png', dpi=300)
plt.show()
