# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:37:00 2022

@author: patag
"""
import cantera as ct
import numpy as np

# CODE PLAN
# 1. pressure chamber conditions
#     user inputs    
#         pressure (stagnation pressure)
#         temperature
#         Test gas mass flow rate through the mass flow control valve
    
#     calculate
#         gas specific heat capacity at P, T (cantera)
    
#     output
#         Test gas mass flow rate
#         C_p test gas

# USER INPUT
P1_0 = 20000  # Pa
T1 = 300    # K
massflow_inp =  5   # g/s
gas_mass_fractions = "CO2:1"
mechanism = "airNASA9noions.cti"

# CALCULATION
massflow = massflow_inp/1000    # kg/s
gas = ct.ThermoPhase(mechanism)
h_formation = gas.h

gas.TPY = T1, P1_0, gas_mass_fractions
cp_gas1 = gas.cp_mass

# 2. arc chamber conditions
#     input from pressure chamber conditions
#         gas specific heat capacity
#         mass flow rate through the mass flow control valve
    
#     user inputs
#         arc power (VxI) max = 250 kW
#
#         get the following from a convection/conduction heating calculation 
#         alternatively assume a thermal efficiency !!! do this initially to get script running
#         C_p of water   
#             temperature difference between inlet and outlet temperature of the cooling water in ARC CHAMBER
#             temperature difference between inlet and outlet temperature of the cooling water in NOZZLE AND MIXING CHAMBER
#             mass flow of cooling water ARC CHAMBER
#             mass flow of cooling water NOZZLE AND MIXING CHAMBER
    
#     calculate 
#         test gas specific stagnation enthalpy at nozzle exit 
    
#     output
#         test gas specific stagnation enthalpy
#         gas()

# USER INPUT
efficiency = 0.6
power_in_inp = 100  # kW

# CALCULATIONS
power_in = power_in_inp*1000    # W
h0 = efficiency*power_in/massflow
# print('stagnation h (MJ/kg)= '+str(h0/(10**6)))
gas.HP = h0 - h_formation, P1_0
gas.equilibrate("SP")
gas()

# 3. settling chamber
#     Assumption
#         ignore in model
    
# 4. nozzle 
#     Assumption
#         adiabatic !!!! this is a bad assumption, heat transfer is occuring the nozzle
    
#     input 
#         gas()
        
        

    