# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:35:24 2022

@author: patag
"""

import sys
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Cantera\SDToolbox\Python3\sdtoolbox')
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Capstone-code\real_gas_and_stagnation.py')
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Capstone-code')
from real_gas_shock_and_stagnation import shock_perfect_gas, shock_calorically_imperfect, stagnation_conditions_isentropic, freestream_to_shock_to_stagnation_condition_calorically_imperfect
import cantera as ct
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.pyplot import contour, contourf

gas = ct.ThermoPhase("airNASA9noions.cti")  # Import the thermodynamic data as NASA9 format. This data contains the chemical species CO2 CO C O O2 C2


# def radical_oxygen_mass_fraction_PWKnozzleflow(gas, enthalpy_addition, composition, P1):
#     gas = gas
#     enthalpy_0 = gas.h
#     print('enthalpy_0 = ' + str(enthalpy_0))
#     print("initial gas enthalpy = " + str(enthalpy_0))
#     enthalpy_1 = -8.9415*10**6 + enthalpy_addition#7.0585*10**6#enthalpy_0 + enthalpy_addition
#     # arc chamber, enthalpy addition, constant pressure
#     gas.HPX = enthalpy_1, P1, composition
#     # nozzle, equilibriate, constant enthalpy and pressure
#     gas.equilibrate('HP')
    
#     gas()
#     # shock
#     # rho2, u2, p2, h2, e2, gas = shock_calorically_imperfect(gas=gas, u1=2850)
#     # stagnation
#     # stagnation_conditions_isentropic(gas=gas, u2=u2)
#     print(gas.species_names)
#     index = gas.species_names.index("O")
#     print(index)
#     print(gas.Y)
#     O_mass_fraction = gas.Y[index]
#     print(O_mass_fraction)

radical_oxygen_mass_fraction_PWKnozzleflow(gas = gas, enthalpy_addition = 16*10**6, composition = "CO2:1", P1 = 90*1000) 
  


#gas.TPX = T1=, p1, 'CO2:1'