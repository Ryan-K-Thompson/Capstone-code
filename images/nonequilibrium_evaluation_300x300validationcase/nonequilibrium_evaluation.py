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

# LOAD DATA FROM SPREADSHEET
data = read_csv("species_data_300x300_viscous.csv")
data=data.dropna()
x_cfd = data['x'].tolist()
x_cfd = [i/(0.0254) for i in x_cfd]
P_cfd = data['P'].tolist()
T_cfd = data['T'].tolist()
CO2_cfd = data['CO2'].tolist()
CO_cfd = data['CO'].tolist()
C2_cfd = data['C2'].tolist()
O2_cfd = data['O2'].tolist()
C_cfd = data['C'].tolist()
O_cfd = data['O'].tolist()
data.dropna()

gas = ct.ThermoPhase("airNASA9noions.cti")  # Import the thermodynamic data as NASA9 format. This data contains the chemical species CO2 CO C O O2 C2

CO2_eq = []
CO_eq = []
C2_eq = []
O2_eq = []
C_eq = []
O_eq = []

for count in range(300):
    mass_fraction = "CO2:1"
    # print(count)
    # print(P_cfd[count])
    P = int(P_cfd[count])
    T = int(T_cfd[count])

    gas.TPY = T, P, mass_fraction
    gas.equilibrate('TP')
    # print(gas.Y)
    CO2_eq.append(gas.Y[gas.species_index('CO2')])
    CO_eq.append(gas.Y[gas.species_index('CO')])
    C2_eq.append(gas.Y[gas.species_index('C2')])
    O2_eq.append(gas.Y[gas.species_index('O2')])
    C_eq.append(gas.Y[gas.species_index('C')])
    O_eq.append(gas.Y[gas.species_index('O')])
    

plt.plot(x_cfd,CO2_cfd, color = 'b', linestyle = 'solid')
plt.plot(x_cfd,CO2_eq, color = 'b', linestyle = 'dashed')
plt.plot(x_cfd,CO_cfd, color = 'r', linestyle = 'solid')
plt.plot(x_cfd,CO_eq, color = 'r', linestyle = 'dashed')
plt.plot(x_cfd,O2_cfd, color = 'm', linestyle = 'solid')
plt.plot(x_cfd,O2_eq, color = 'm', linestyle = 'dashed')
plt.plot(x_cfd,O_cfd, color = 'c', linestyle = 'solid')
plt.plot(x_cfd,O_eq, color = 'c', linestyle = 'dashed')
plt.plot(x_cfd,C2_cfd, color = 'g', linestyle = 'solid')
plt.plot(x_cfd,C2_eq, color = 'g', linestyle = 'dashed')
plt.plot(x_cfd,C_cfd, color = 'y', linestyle = 'solid')
plt.plot(x_cfd,C_eq, color = 'y', linestyle = 'dashed')

plt.yscale('log')

plt.ylim(10**-5,1)


plt.legend(['$\mathregular{CO_{2}}$ (CFD)','$\mathregular{CO_{2}}$ (Equilibrium)','CO (CFD)','CO (Equilibrium)','$\mathregular{O_{2}}$ (CFD)','$\mathregular{O_{2}}$ (Equilibrium)','O (CFD)','O (Equilibrium)','$\mathregular{C_{2}}$ (CFD)','$\mathregular{C_{2}}$ (Equilibrium)','C (CFD)','C (Equilibrium)'],bbox_to_anchor=(1,1))
plt.grid()
plt.xlim(-0.035,0)
plt.xticks([-0.03, -0.02, -0.01,0])
plt.xlabel('X/R')
plt.ylabel('Mass fraction')
# plt.title('Stagnation line chemical non-equilibrium')
plt.savefig('filename.png', dpi=300, bbox_inches='tight')
