# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:35:24 2022

@author: patag
"""

import sys
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Cantera\SDToolbox\Python3\sdtoolbox')
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Capstone-code')
from real_gas_shock_and_stagnation import radical_oxygen_mass_fraction_PWKnozzleflow
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import contour, contourf

gas = ct.ThermoPhase("airNASA9noions.cti")  # Import the thermodynamic data as NASA9 format. This data contains the chemical species CO2 CO C O O2 C2

enthalpy_linspace = np.linspace(0,25*10**6,250)
O2_MF = 0*enthalpy_linspace
O_MF = 0*enthalpy_linspace
N2_MF = 0*enthalpy_linspace
N_MF = 0*enthalpy_linspace
NO_MF = 0*enthalpy_linspace
C_MF = 0*enthalpy_linspace
C2_MF = 0*enthalpy_linspace
CO_MF = 0*enthalpy_linspace
CO2_MF = 0*enthalpy_linspace
CN_MF = 0*enthalpy_linspace



count = 0
for enthalpy in enthalpy_linspace:
    MF = radical_oxygen_mass_fraction_PWKnozzleflow(gas = gas, enthalpy_addition = enthalpy, composition = "CO2:0.5, N2:0.5", P1 = 90*1000) 
    O2_MF[count] = MF[0]
    O_MF[count] = MF[1]
    N2_MF[count] = MF[2]
    N_MF[count] = MF[3]
    NO_MF[count] = MF[4]
    C_MF[count] = MF[5]
    C2_MF[count] = MF[6]
    CO_MF[count] = MF[7]
    CO2_MF[count] = MF[8]
    CN_MF[count] = MF[9]
    count = count + 1

plt.figure()
plt.title("Mixing chamber gas composition" + "\n" + "Initial compostion (CO$_2$:90%, N$_2$:10%)")
plt.xlabel("Arc chamber enthalpy addition (MJ)")
plt.ylabel("Species mass fraction")
plt.ylim(bottom = 0, top = 1)
plt.xlim(left = 0, right = 25)
enthalpy_linspace = enthalpy_linspace/(10**6)
plt.plot(enthalpy_linspace, O2_MF)
plt.plot(enthalpy_linspace, O_MF)
plt.plot(enthalpy_linspace, N2_MF)
plt.plot(enthalpy_linspace, N_MF)
plt.plot(enthalpy_linspace, NO_MF)
plt.plot(enthalpy_linspace, C_MF)
plt.plot(enthalpy_linspace, C2_MF)
plt.plot(enthalpy_linspace, CO_MF)
plt.plot(enthalpy_linspace, CO2_MF)
plt.plot(enthalpy_linspace, CN_MF)
plt.legend(["$O_2$","$O$", "$N_2$", "$N$", "$NO$", "$C$", "$C_2$", "$CO$", "$CO_2$", "$CN$"], loc = "upper left")

# enthalpy_linspace = np.linspace(0,25*10**6,250)
# O2_MF = 0*enthalpy_linspace
# O_MF = 0*enthalpy_linspace
# N2_MF = 0*enthalpy_linspace
# N_MF = 0*enthalpy_linspace
# NO_MF = 0*enthalpy_linspace
# C_MF = 0*enthalpy_linspace
# C2_MF = 0*enthalpy_linspace
# CO_MF = 0*enthalpy_linspace
# CO2_MF = 0*enthalpy_linspace
# CN_MF = 0*enthalpy_linspace

# count = 0
# for enthalpy in enthalpy_linspace:
#     MF = radical_oxygen_mass_fraction_PWKnozzleflow(gas = gas, enthalpy_addition = enthalpy, composition = "CO2:0.95, N2:0.05", P1 = 90*1000) 
#     O2_MF[count] = MF[0]
#     O_MF[count] = MF[1]
#     N2_MF[count] = MF[2]
#     N_MF[count] = MF[3]
#     NO_MF[count] = MF[4]
#     C_MF[count] = MF[5]
#     C2_MF[count] = MF[6]
#     CO_MF[count] = MF[7]
#     CO2_MF[count] = MF[8]
#     CN_MF[count] = MF[9]
#     count = count + 1

# plt.figure()

# plt.title("Mixing chamber gas composition"+"\n"+"Initial compostion (CO$_2$:95%, N$_2$:5%)", loc = 'center')
# plt.xlabel("Arc chamber enthalpy addition (MJ)")
# plt.ylabel("Species mass fraction")
# plt.ylim(bottom = 0, top = 1)
# plt.xlim(left = 0, right = 25)
# enthalpy_linspace = enthalpy_linspace/(10**6)
# plt.plot(enthalpy_linspace, O2_MF)
# plt.plot(enthalpy_linspace, O_MF)
# plt.plot(enthalpy_linspace, N2_MF)
# plt.plot(enthalpy_linspace, N_MF)
# plt.plot(enthalpy_linspace, NO_MF)
# plt.plot(enthalpy_linspace, C_MF)
# plt.plot(enthalpy_linspace, C2_MF)
# plt.plot(enthalpy_linspace, CO_MF)
# plt.plot(enthalpy_linspace, CO2_MF)
# plt.plot(enthalpy_linspace, CN_MF)
# plt.legend(["$O_2$","$O$", "$N_2$", "$N$", "$NO$", "$C$", "$C_2$", "$CO$", "$CO_2$", "$CN$"], loc = "upper left")
