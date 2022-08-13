# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:37:00 2022

@author: patag
"""
import cantera as ct
import numpy as np
import matplotlib
from matplotlib import pyplot as plt




def PWK(P1_0, power_in_inp, massflow_inp):
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
    # P1_0 = 20000  # Pa
    T1 = 300    # K
    # massflow_inp =  5   # g/s
    gas_mass_fractions = "CO2:1"
    mechanism = "airNASA9noions.cti"
    
    # CALCULATION
    massflow = massflow_inp/1000    # kg/s
    gas = ct.ThermoPhase(mechanism)
    
    gas.TPY = T1, P1_0, gas_mass_fractions
    cp_gas1 = gas.cp_mass
    h_initial = gas.h
    gas()
    
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
    efficiency = 0.5
    # power_in_inp = 18  # kW
    r_chamber = 0.025   # m, arc chamber radius
    A_chamber = np.pi*(r_chamber**2)    # m2, area of arc chamber
    A_ratio = 4
    R_model = 0.05
    
    
    # CALCULATIONS
    
    A_exit = A_chamber*4
    R_jet = (A_exit/np.pi)**0.5
    
    rho1 = gas.density_mass 
    U1 = massflow/(rho1*A_chamber)
    
    
    
    power_in = power_in_inp*1000    # W
    h0 = efficiency*power_in/massflow
    print('stagnation h (MJ/kg)= '+str(h0/(10**6)))
    gas.HP = h0 + h_initial, P1_0
    # gas.equilibrate("HP")
    # gas()
    
    # gas.HP = gas.h,10000
    gas.equilibrate("SP")
    gas()
    rho2 = gas.density_mass
    
    U2 = (rho1/rho2)*U1*(1/A_ratio)
    print(U2)
    
    print('density'+str(rho2))
    
    # 3. settling chamber
    #     Assumption
    #         ignore in model
        
    # 4. nozzle 
    #     Assumption
    #         adiabatic !!!! this is a bad assumption, heat transfer is occuring the nozzle
        
    #     input 
    #         gas()
    L = R_model/R_jet
    
    beta = (U2/R_model)*(1/(2-L-1.68*(L-1)**2-1.28*(L-1)**3))
    
    P_free = P1_0 - 0.5*rho2*U2**2
    print(P1_0/P_free)

    return P1_0, h0, beta
 

P1_0 = range(20000,80000,5000)
massflow = 50
power = range(50,250,25)

P0 = np.zeros([len(P1_0), len(power)])
H0 = np.zeros([len(P1_0), len(power)])
beta = np.zeros([len(P1_0), len(power)])

for count, P in enumerate(P1_0):
    for count2, p in enumerate(power):
        P0[count][count2], H0[count][count2], beta[count][count2] = PWK(P,p,massflow)

fig = plt.figure(figsize =(14, 9))
ax = plt.axes(projection ='3d')
ax.plot_surface(P0,H0/(10**6),beta)
plt.xlabel("P0")
plt.ylabel("H0")
plt.zlabel("beta")