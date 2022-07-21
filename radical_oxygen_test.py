# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:35:24 2022

@author: patag
"""

import sys
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Cantera\SDToolbox\Python3\sdtoolbox')
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Capstone-code')
from real_gas_shock_and_stagnation import shock_perfect_gas, shock_calorically_imperfect, stagnation_conditions_isentropic, freestream_to_shock_to_stagnation_condition_calorically_imperfect
import cantera as ct
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.pyplot import contour, contourf





gas = ct.ThermoPhase("airNASA9noions.cti")  # Import the thermodynamic data as NASA9 format. This data contains the chemical species CO2 CO C O O2 C2

# Define constants
gasconstant = 8.3144598 

# arc chamber:
gas.HPX = 7.0585*10**6,90*10**3,'CO2:1'
# nozzle:
gas.equilibrate('HP')
# shock
rho2, u2, p2, h2, e2, gas = shock_calorically_imperfect(gas=gas, u1=2850)
# stagnation
stagnation_conditions_isentropic(gas=gas, u2=u2)

gas()

#gas.TPX = T1=, p1, 'CO2:1'