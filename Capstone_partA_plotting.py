# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:50:08 2022

@author: patag
"""

# =============================================================================
# below is the loop for the CAPSTONE A1 graphing
# =============================================================================
import sys
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Cantera\SDToolbox\Python3\sdtoolbox')
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Capstone-code')
from real_gas_shock_and_stagnation import shock_perfect_gas, shock_calorically_imperfect, stagnation_conditions_isentropic, SDToolbox_postshock_velocity
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import contourf, contour
import sdtoolbox
from sdtoolbox import postshock
from postshock import PostShock_eq, PostShock_fr


gas = ct.ThermoPhase("airNASA9noions.cti")  # Import the thermodynamic data as NASA9 format. This data contains the chemical species CO2 CO C O O2 C2

global H_total, P_total

Ts = 180


pdim = 41
Udim = 20

Ps = np.linspace(0,20,num=pdim)
print(Ps)
Us = np.linspace(1000,6000, num = Udim)

tick_dim_x = 9
tick_label_x = np.linspace(min(Us),max(Us),tick_dim_x)




H_stagnation = np.zeros((pdim,Udim))#[ [0]*3 for i in range(3)]
P_stagnation = np.zeros((pdim,Udim))



for counti, i in enumerate(Ps):
    for countj, j in enumerate(Us):
        try:
            p1 = i
            T1 = Ts
            u1 = j
            gas.TPX = T1, p1, 'CO2:1'
            rho1 = gas.density
            h1 = gas.h
            rho2, u2, p2, h2, e2, gas = shock_calorically_imperfect(gas, u1)
            stagnation_conditions_isentropic(gas, u2)
            
            
            H_stagnation[counti][countj] = H_total
            P_stagnation[counti][countj] = P_total
        except:
            print("error at condition: p1 =" +str(p1) + ", u1 = " + str(u1))
            pass
        
#print(H_stagnation)
#print(P_stagnation)

np.reshape(H_stagnation, (pdim,Udim))
np.reshape(P_stagnation, (pdim,Udim))
#print(H_stagnation)
Ps_axis, Us_axis = np.meshgrid(Ps, Us)        

fig1 = plt.figure("Figure 1")
contourf(Us, Ps, P_stagnation/1000,100)       
plt.colorbar()
# plot1a = contour(Us, Ps,P_stagnation/1000, [17.500], colors = 'k')  
# plt.clabel(plot1a,[17.500], inline = True, rightside_up = False, colors = 'k')

plt.title("FLight Stagnation Pressure (kPa)")
plt.yticks(np.arange(0,101,step = 10))
plt.xticks(np.arange(1000,6001,step=1000))
plt.ylabel("Static Pressure (Pa)")
plt.xlabel("Velocity (m/s)")
#plt.xaxis.set_major_locator(plt.MaxNLocator(3))


fig2 = plt.figure("Figure 2")
plot2 = contourf(Us, Ps, H_stagnation/1000000, 100)
plt.colorbar()
# plot2a = contour(Us, Ps,H_stagnation/1000000, [7.058000], colors = 'k')  
# plt.clabel(plot2a,[7.058000], inline = True, rightside_up = False, colors = 'k')

plt.title("Flight Stagnation Enthalpy (MJ/kg)")
plt.yticks(np.arange(0,101,step = 10))
plt.xticks(np.arange(1000,6001,step=1000))
plt.ylabel("Static Pressure (Pa)")
plt.xlabel("Velocity (m/s)")





fig3 = plt.figure("Figure 3")
plot31 = contour(Us, Ps, P_stagnation,[17500], colors = 'k')
plot3b = contour(Us, Ps, H_stagnation,[7058000], colors = 'k')
plt.ylabel("Static Pressure (Pa)")
plt.xlabel("Velocity (m/s)")
plt.title("Matching point")
plt.text(3800,38,"(3725 m/s,36.7 Pa)")

# contourf(H_stagnation)
# contourf(P_stagnation)
# plt.legend
# #print(P1)

 
# rho2, u2, p2, h2, e2, gas = shock_calorically_imperfect(gas)

# print(rho2, u2, p2, h2, e2) 
# gas()
# stagnation_conditions_isentropic(gas)




