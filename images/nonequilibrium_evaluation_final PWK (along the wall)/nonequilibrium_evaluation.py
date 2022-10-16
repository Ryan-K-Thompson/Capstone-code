# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 12:48:22 2022

@author: patag
"""

import csv
from pandas import *
import cantera as ct
import matplotlib.pyplot as plt
import math

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


# LOAD DATA FROM SPREADSHEET
data = read_csv("species_data_PWK_wall.csv")
data=data.dropna()
x_PWK = data['x'].tolist()
y_PWK = data['y'].tolist()

s_PWK = XY_to_pathlength(x_PWK, y_PWK)

s_div_r_PWK = scale_by(s_PWK, 1/0.0254)


P_PWK = data['P'].tolist()
T_PWK = data['T'].tolist()
CO2_PWK = data['CO2'].tolist()
CO_PWK = data['CO'].tolist()
C2_PWK = data['C2'].tolist()
O2_PWK = data['O2'].tolist()
C_PWK = data['C'].tolist()
O_PWK = data['O'].tolist()

data = read_csv("species_data_flight_wall.csv")
data=data.dropna()
x_flight = data['x'].tolist()
y_flight = data['y'].tolist()

s_flight = XY_to_pathlength(x_flight, y_flight)

s_div_r_flight = scale_by(s_flight, 1/1.325)


P_flight = data['P'].tolist()
T_flight = data['T'].tolist()
CO2_flight = data['CO2'].tolist()
CO_flight = data['CO'].tolist()
C2_flight = data['C2'].tolist()
O2_flight = data['O2'].tolist()
C_flight = data['C'].tolist()
O_flight = data['O'].tolist()


# gas = ct.ThermoPhase("airNASA9noions.cti")  # Import the thermodynamic data as NASA9 format. This data contains the chemical species CO2 CO C O O2 C2

# CO2_eq = []
# CO_eq = []
# C2_eq = []
# O2_eq = []
# C_eq = []
# O_eq = []

# for count in range(572):
#     mass_fraction = "CO2:1"
#     # print(count)
#     # print(P_cfd[count])
#     P = int(P_cfd[count])
#     T = int(T_cfd[count])

#     gas.TPY = T, P, mass_fraction
#     gas.equilibrate('TP')
#     # print(gas.Y)
#     CO2_eq.append(gas.Y[gas.species_index('CO2')])
#     CO_eq.append(gas.Y[gas.species_index('CO')])
#     C2_eq.append(gas.Y[gas.species_index('C2')])
#     O2_eq.append(gas.Y[gas.species_index('O2')])
#     C_eq.append(gas.Y[gas.species_index('C')])
#     O_eq.append(gas.Y[gas.species_index('O')])
    



plt.plot(s_div_r_flight,CO2_flight, color = 'b', linestyle = 'solid')
plt.plot(s_div_r_flight,CO_flight, color = 'r', linestyle = 'solid')
plt.plot(s_div_r_flight,O2_flight, color = 'm', linestyle = 'solid')
plt.plot(s_div_r_flight,O_flight, color = 'c', linestyle = 'solid')
plt.plot(s_div_r_flight,C2_flight, color = 'g', linestyle = 'solid')
plt.plot(s_div_r_flight,C_flight, color = 'y', linestyle = 'solid')

plt.plot(s_div_r_PWK,CO2_PWK, color = 'b', linestyle = 'dashed')
plt.plot(s_div_r_PWK,CO_PWK, color = 'r', linestyle = 'dashed')
plt.plot(s_div_r_PWK,O2_PWK, color = 'm', linestyle = 'dashed')
plt.plot(s_div_r_PWK,O_PWK, color = 'c', linestyle = 'dashed')
plt.plot(s_div_r_PWK,C2_PWK, color = 'g', linestyle = 'dashed')
plt.plot(s_div_r_PWK,C_PWK, color = 'y', linestyle = 'dashed')

# plt.yscale('log')

plt.ylim(10**-20,1)


plt.legend(['$\mathregular{CO_{2}}$ (CFD)','CO (CFD)','$\mathregular{O_{2}}$ (CFD)','O (CFD)','O (Equilibrium)','C (CFD)',],bbox_to_anchor=(1,1))
plt.grid()
plt.xlim(0,1.2)
# plt.xticks([-0.03, -0.02, -0.01,0])
plt.xlabel('S/R')
plt.ylabel('Mass fraction')
# plt.title('Stagnation line chemical non-equilibrium')
plt.savefig('filename.png', dpi=300, bbox_inches='tight')
