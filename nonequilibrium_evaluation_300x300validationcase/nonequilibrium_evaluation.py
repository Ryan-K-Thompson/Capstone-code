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
x_cfd = [i/(-0.0254) for i in x_cfd]
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
    

plt.plot(x_cfd,CO2_cfd, color = 'blue')
plt.plot(x_cfd,CO2_eq, color = 'blue', linestyle = 'dashed')
plt.plot(x_cfd,CO_cfd, color = 'green')
plt.plot(x_cfd,CO_eq, color = 'green', linestyle = 'dashed')
plt.plot(x_cfd,O2_cfd, color = 'pink')
plt.plot(x_cfd,O2_eq, color = 'pink', linestyle = 'dashed')
plt.plot(x_cfd,O_cfd, color = 'yellow')
plt.plot(x_cfd,O_eq, color = 'yellow', linestyle = 'dashed')
plt.plot(x_cfd,C2_cfd, color = 'grey')
plt.plot(x_cfd,C2_eq, color = 'grey', linestyle = 'dashed')
plt.plot(x_cfd,C_cfd, color = 'red')
plt.plot(x_cfd,C_eq, color = 'red', linestyle = 'dashed')

plt.yscale('log')
plt.xscale('log')
plt.ylim(10**-5,1)
plt.xlim(-0.05,0)

plt.legend(['CO_2 (CFD)','CO2 (equilibrium)','CO (CFD)','CO (equilibrium)','O2 (CFD)','O2 (equilibrium)','O (CFD)','O (equilibrium)','C2 (CFD)','C2 (equilibrium)','C (CFD)','C (equilibrium)',])
plt.grid()
plt.xlabel('abs(X/R)')
plt.ylabel('Mass fraction')
plt.title('Stagnation line chemical non-equilibrium')
plt.savefig('filename.png', dpi=300)
