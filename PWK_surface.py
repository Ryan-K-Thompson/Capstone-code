# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:37:00 2022

@author: patag
"""
import cantera as ct
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import sympy 
from sympy import symbols, solve
import scipy
from scipy.optimize import fsolve
import mpmath as mp
from mpmath import findroot

def function():
    # USER INPUT
    P1_0 = 20000  # Pa
    P1_P4 = 1.1
    P4 = P1_0/P1_P4

    
    
    T1 = 300    # K
    massflow_inp =  5   # g/s
    gas_mass_fractions = "CO2:1"
    mechanism = "airNASA9noions.cti"
    
    # CALCULATION
    massflow = massflow_inp/1000    # kg/s
    gas = ct.ThermoPhase(mechanism)
    
    gas.TPY = T1, P1_0, gas_mass_fractions
    cp_gas1 = gas.cp_mass
    h_initial = gas.h
    gas()
    

    efficiency = 0.5
    power_in_inp = 180  # kW
    r_chamber = 0.025   # m, arc chamber radius
    A_chamber = np.pi*(r_chamber**2)    # m2, area of arc chamber
    
    R_model = 0.05
    
    A_exit = A_chamber*4
    R_jet = (A_exit/np.pi)**0.5
    
    rho1 = gas.density_mass 
    U1 = massflow/(rho1*A_chamber)
    
    
    
    power_in = power_in_inp*1000    # W
    h0 = efficiency*power_in/massflow
    print('stagnation h (MJ/kg)= '+str(h0/(10**6)))
    gas.HP = h0 + h_initial, P1_0
    
    gas.equilibrate("SP")
    
    
    gamma = gas.cp_mass/gas.cv_mass
    print("gamma "+str(gamma))
    
    
    
    Me = ((P1_P4**((gamma-1)/gamma)-1)/((gamma-1)/2))**0.5
    print("M = " +str(Me))
    Ue = Me*(gamma*gas.T*8.314)**0.5
    print('Ue = ' +str(Ue))
    Pe = P1_0 - 0.5*gas.density_mass*Ue**2
    print("Pe = "+str(Pe))
    Te = gas.T
    print("Te = "+str(Te))
    
    L = R_model/R_jet
    
    beta = (Ue/R_model)*(1/(2-L-1.68*(L-1)**2-1.28*(L-1)**3))
    print("beta = "+str(beta))

# return P1_0, h0, beta, U2, Temp, Pres

 

# P1_0 = range(20000,80000,5000)
# massflow = 25
# power = range(50,250,25)

# P0 = np.zeros([len(P1_0), len(power)])
# H0 = np.zeros([len(P1_0), len(power)])
# beta = np.zeros([len(P1_0), len(power)])
# T_ = np.zeros([len(P1_0), len(power)])
# P_ = np.zeros([len(P1_0), len(power)])
# U2_ = np.zeros([len(P1_0), len(power)])

# for count, P in enumerate(P1_0):
#     for count2, p in enumerate(power):
#         print(2)
#         P0[count][count2], H0[count][count2], beta[count][count2], U2_[count][count2], T_[count][count2], P_[count][count2] = PWK(P,p,massflow)



# plt.figure()
# fig2 = plt.contourf(P0/1000,H0/10**6,beta)
# plt.xlabel("P_0 (kPa)")
# plt.ylabel("H_0 (MJ)")
# plt.colorbar()
# plt.title("beta")

# plt.figure()
# fig2 = plt.contourf(P0/1000,H0/10**6,U2_)
# plt.xlabel("P_0 (kPa)")
# plt.ylabel("H_0 (MJ)")
# plt.colorbar()
# plt.title("U2")

# plt.figure()
# fig2 = plt.contourf(P0/1000,H0/10**6,T_)
# plt.xlabel("P_0 (kPa)")
# plt.ylabel("H_0 (MJ)")
# plt.colorbar()
# plt.title("T")

# plt.figure()
# fig2 = plt.contourf(P0/1000,H0/10**6,P_)
# plt.xlabel("P_0 (kPa)")
# plt.ylabel("H_0 (MJ)")
# plt.colorbar()
# plt.title("Ps")

# plt.figure()
# fig = plt.figure(figsize =(7,7))
# ax = plt.axes(projection ='3d')
# ax.plot_surface(P0/1000,H0/(10**6),beta)
# plt.xlabel("P_0 (kPa)")
# plt.ylabel("H_0 (MJ)")
# plt.title("Plasma wind tunnel operational surface")
# ax.set_zlabel("Beta (1/s)")
# # plt.zlabel("beta")
