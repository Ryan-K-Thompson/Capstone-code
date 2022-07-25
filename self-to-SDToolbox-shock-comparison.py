# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 08:44:38 2022

@author: patag
"""

"""
Purpose of this script is to compare a homebuilt implimentation of Cantera with
the popular python module SDToolbox, which provides many functions to evaluate
the post shock thermodynamic and chemical state of gasses.

We will first take a shock at 1000 m/s, 100 kPa, 25 C
"""

import sys
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Cantera\SDToolbox\Python3\sdtoolbox')
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Capstone-code')
from real_gas_shock_and_stagnation import shock_perfect_gas, shock_calorically_imperfect, stagnation_conditions_isentropic, SDToolbox_postshock_velocity
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import sdtoolbox
from sdtoolbox import postshock
from postshock import PostShock_eq, PostShock_fr

T1 = 182
P1 = 72
u1 = 2273

composition = "CO2:1"

# own script
gas_base = ct.ThermoPhase("airNASA9noions.cti")
gas_base.TPY = T1, P1, composition
gas_base.equilibrate("TP")

gas1 = gas_base
gas2 = gas_base



return_self = shock_calorically_imperfect(gas = gas1, u1 = u1)

u2_self = return_self[1]

# sdtoolbox
gas2_postshock = PostShock_fr(U1 = u1, P1 = P1, T1 = T1, q = composition, mech = "airNASA9noions.cti")
gas_base = ct.ThermoPhase("airNASA9noions.cti")
gas_base.TPY = T1, P1, composition
gas_base.equilibrate("TP")
u2_SDT = SDToolbox_postshock_velocity(gas1=gas_base, gas2=gas2_postshock, u1=u1)


print("self u2: "+str(u2_self))
print("sdt u2: "+str(u2_SDT))
