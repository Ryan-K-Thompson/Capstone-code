# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 15:55:35 2022

@author: Ryan Thompson

This script evaluates the thermodynamic parameters across a shockwave and into 
the boundary layer of an entry body in chemically frozen flow, that is 
calorically imperfect (i.e. variable ratio of specific heats) and thermally 
perfect.

"""

# Import modules and data
import sys
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Cantera\SDToolbox\Python3\sdtoolbox')
import cantera as ct
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.pyplot import contour, contourf
gas = ct.ThermoPhase("airNASA9noions.cti")  # Import the thermodynamic data as NASA9 format. This data contains the chemical species CO2 CO C O O2 C2

# Define constants
gasconstant = 8.3144598 

# Define initial conditions
condition = "PWK"   # PWK or Pathfinder   

def shock_perfect_gas(gas, u1):
    """
    This function evaluates the change in thermodynamic properties across a 
    shockwave for a perfect gas
    
    Its current application is to provide the initial guess values for the 
    evaluation of thermodynamic properties across a shockwave for a calorically 
    imperfect, thermally perfect gas 
    (i.e. function shock_calorically_imperfect())
    (If other uses are required, this may need to be ammended as it is built
    for purpose)
    
    Need to estimate rho2, u2, p2, h2, e2 post shock using the ideal gas 
    assumptions
    
    """
    gamma = gas.cp/gas.cv
    R = gasconstant/gas.mean_molecular_weight*1000
    speed_of_sound = (gamma*R*gas.T)**0.5
    M1 = u1/speed_of_sound
    M2 = ((1+((gamma-1)/2)*M1**2)/(gamma*M1**2-((gamma-1/2))))**0.5
    P2divP1 = 1+((2*gamma)/(gamma+1))*(M1**2-1)
    rho2divrho1 = ((gamma+1)*M1**2)/(2+(gamma-1)*M1**2)
    V1divV2 = rho2divrho1
    T2divT1 = (1+((2*gamma)/(gamma+1))*(M1**2-1))*((2+(gamma-1)*M1**2)/((gamma+1)*M1**2))
    h2divh1 = T2divT1
    
    
    rho2 = rho2divrho1*gas.density
    u2 = u1/V1divV2
    p2 = P2divP1*gas.P
    T_estimate = p2/(rho2*R)
    print("estimate"+str(T_estimate)+" at condition p1 = " +str(gas.P)+", u1 = " +str(u1))

        
    h2 = h2divh1*gas.h
    
    
    return rho2, u2, p2, h2
    
def shock_calorically_imperfect(gas, u1):
    """
    This function determines the thermodynamic properties of a calorically 
    imperfect, thermally perfect gas passing through a shockwave.
    
    The method is numerical and requires an initial guess, which is determined 
    using the inviscid shock relations
    """
    
    initial_guess = [(*shock_perfect_gas(gas, u1),gas.int_energy_mass)]
    
    print("guess = "+str(initial_guess))
    #initial_guess = 
    rho1 = gas.density
    p1 = gas.P
    h1 = gas.h
    gas_temp = gas
    def equations(vars):
        #print("gamma2 = " + str(gas.cp/gas.cv))
        rho2, u2, p2, h2, e2 = vars
        gas_temp.HP = h2, p2
        eq1 = rho1*u1-rho2*u2
        eq2 = (p1-p2)+(rho1*(u1**2)-rho2*(u2**2))
        eq3 = 0.5*(u1**2-u2**2)+(h1-h2)
        gas_temp.UV = e2, 1/rho2
        eq4 = p2-gas_temp.P
        eq5 = e2-h2+p2/rho2
        #print(rho2, u2, p2, h2, e2)
        
        return [eq1, eq2, eq3, eq4, eq5]
    gas = gas_temp
    rho2, u2, p2, h2, e2 = fsolve(equations, initial_guess, factor=2)
    print("post shock conditions")
    print("U2 = " + str(u2))
    print("P2 = "+str(gas.P))
    print("T2 = "+ str(gas.T))
    return rho2, u2, p2, h2, e2, gas

    

def stagnation_conditions_isentropic(gas, u2):
    global H_total, P_total
    E_total = gas.P/(gas.cp/gas.cv-1)       # from metacomp guide, file:///C:/Users/patag/Desktop/Uni/2022/Capstone/Simulations/CFDpp_V19.1_Manual/Default.htm#ICFD++/TheoryGuide/ConservationEquationsAndClosures/CompRGasesEOS.html%3FTocPath%3DICFD%252B%252B%7CTheory%2520Guide%7CConservation%2520Equations%2520and%2520Closures%7CClosures%7CCompressible%2520Real%2520Gases%7C_____1
    #H_total = (E_total+gas.P)/gas.density   # from metacomp guide, file:///C:/Users/patag/Desktop/Uni/2022/Capstone/Simulations/CFDpp_V19.1_Manual/Default.htm#ICFD++/TheoryGuide/ConservationEquationsAndClosures/CompRGasesEOS.html%3FTocPath%3DICFD%252B%252B%7CTheory%2520Guide%7CConservation%2520Equations%2520and%2520Closures%7CClosures%7CCompressible%2520Real%2520Gases%7C_____1
    #H_total = E_total+(gas.P)/gas.density   # from, https://www.aoe.vt.edu/content/dam/aoe_vt_edu/programs/graduate/forms/lectnotes3-09All101812.pdf
    P_total = gas.P+gas.density*(u2**2)     # from momentum equation
    # assume entropy constant
    s = gas.s
    gas.SP = s, P_total
    H_total = gas.h+float(8.9415*10**6)
    gas()
    print("Stagnation enthalpy (J/kg) = " + str(H_total))
    print("Stagnation pressure (Pa) = " + str(P_total))
    return H_total, P_total

def freestream_to_shock_to_stagnation_condition_calorically_imperfect(condition, gas):
    if condition == "PWK":
        p1 = 796.2 # Pa
        T1 = 1165.8 # K
        gas.TPX = T1, p1, "CO2:0.02629, CO:0.61999, C2:0, O2:0.05779, C:0, O:0.29593"   #equilibriate the gas at the initial conditions, from this point forward, the gas will be assumed to be frozen due to the low dakholmer number
        u1 = 2850 # m/s
        
    elif condition == "Pathfinder":
        p1 = 72.3 # Pa
        T1 = 182.1 # K
        gas.TPX = T1, p1, 'CO2:1'   #equilibriate the gas at the initial conditions, from this point forward, the gas will be assumed to be frozen due to the low dakholmer number
        u1 = 2273.5 # m/s
    else:
        print("The condition input is invalid")
    rho1 = gas.density
    h1 = gas.h
    print(""""Valid input conditions...
    Initial conditions:
    U = """ +str(u1)+" m/s")
    gas()
    
    rho2, u2, p2, h2, e2, gas = shock_calorically_imperfect(gas, u1)
    stagnation_conditions_isentropic(gas, u2)
    gas()
  
def radical_oxygen_mass_fraction_PWKnozzleflow(gas, enthalpy_addition, composition, P1):
    gas = gas
    enthalpy_0 = gas.h
    # print('enthalpy_0 = ' + str(enthalpy_0))
    # print("initial gas enthalpy = " + str(enthalpy_0))
    enthalpy_1 = -8.9415*10**6 + enthalpy_addition#7.0585*10**6#enthalpy_0 + enthalpy_addition
    # arc chamber, enthalpy addition, constant pressure
    gas.HPX = enthalpy_1, P1, composition
    # nozzle, equilibriate, constant enthalpy and pressure
    gas.equilibrate('HP')
    
    # gas()
    # shock
    # rho2, u2, p2, h2, e2, gas = shock_calorically_imperfect(gas=gas, u1=2850)
    # stagnation
    # stagnation_conditions_isentropic(gas=gas, u2=u2)
    # print(gas.species_names)
    index = gas.species_names.index("O")
    # print(index)
    # print(gas.Y)
    O_mass_fraction = gas.Y[index]
    mass_fractions = gas.Y
    # print(O_mass_fraction)
    print(mass_fractions)
    return mass_fractions

def SDToolbox_postshock_velocity(gas1, gas2, u1):
    """
    SDToolbox lacks a function to determine post shock velocity. This fills the 
    gap without interacting witht the functional elements of the SDToolbox 
    algorithms
    
    Parameters
    ----------
    gas1 : TYPE
        DESCRIPTION.
    gas2 : TYPE
        DESCRIPTION.
    u1 : TYPE
        DESCRIPTION.

    Returns
    -------
    u2 : TYPE
        DESCRIPTION.

    """
    density1 = gas1.density_mass
    density2 = gas2.density_mass
    u2 = density1*u1/density2
    return u2

radical_oxygen_mass_fraction_PWKnozzleflow(gas = gas, enthalpy_addition = 16*10**6, composition = "CO2:1", P1 = 90*1000) 
  
    
# freestream_to_shock_to_stagnation_condition_calorically_imperfect(condition = "PWK", gas = gas)  


# =============================================================================
# below is the loop for the CAPSTONE A1 graphing
# =============================================================================
# Ts = 180


# pdim = 20
# Udim = 20

# Ps = np.linspace(1,100,num=pdim)
# Us = np.linspace(1000,6000, num = Udim)

# tick_dim_x = 9
# tick_label_x = np.linspace(min(Us),max(Us),tick_dim_x)




# H_stagnation = np.zeros((pdim,Udim))#[ [0]*3 for i in range(3)]
# P_stagnation = np.zeros((pdim,Udim))



# for counti, i in enumerate(Ps):
#     for countj, j in enumerate(Us):
#         try:
#             p1 = i
#             T1 = Ts
#             u1 = j
#             gas.TPX = T1, p1, 'CO2:1'
#             rho1 = gas.density
#             h1 = gas.h
#             rho2, u2, p2, h2, e2, gas = shock_calorically_imperfect(gas, u1)
#             stagnation_conditions_isentropic(gas, u2)
    
            
#             H_stagnation[counti][countj] = H_total
#             P_stagnation[counti][countj] = P_total
#         except:
#             print("error at condition: p1 =" +str(p1) + ", u1 = " + str(u1))
#             pass
        
# #print(H_stagnation)
# #print(P_stagnation)

# np.reshape(H_stagnation, (pdim,Udim))
# np.reshape(P_stagnation, (pdim,Udim))
# #print(H_stagnation)
# Ps_axis, Us_axis = np.meshgrid(Ps, Us)        

# fig1 = plt.figure("Figure 1")
# contourf(Us, Ps, P_stagnation/1000,100)       
# plt.colorbar()
# plot1a = contour(Us, Ps,P_stagnation/1000, [17.500], colors = 'k')  
# plt.clabel(plot1a,[17.500], inline = True, rightside_up = False, colors = 'k')

# plt.title("FLight Stagnation Pressure (kPa)")
# plt.yticks(np.arange(0,101,step = 10))
# plt.xticks(np.arange(1000,6001,step=1000))
# plt.ylabel("Static Pressure (Pa)")
# plt.xlabel("Velocity (m/s)")
# #plt.xaxis.set_major_locator(plt.MaxNLocator(3))


# fig2 = plt.figure("Figure 2")
# plot2 = contourf(Us, Ps, H_stagnation/1000000, 100)
# plt.colorbar()
# plot2a = contour(Us, Ps,H_stagnation/1000000, [7.058000], colors = 'k')  
# plt.clabel(plot2a,[7.058000], inline = True, rightside_up = False, colors = 'k')

# plt.title("Flight Stagnation Enthalpy (MJ/kg)")
# plt.yticks(np.arange(0,101,step = 10))
# plt.xticks(np.arange(1000,6001,step=1000))
# plt.ylabel("Static Pressure (Pa)")
# plt.xlabel("Velocity (m/s)")





# fig3 = plt.figure("Figure 3")
# plot31 = contour(Us, Ps, P_stagnation,[17500], colors = 'k')
# plot3b = contour(Us, Ps, H_stagnation,[7058000], colors = 'k')
# plt.ylabel("Static Pressure (Pa)")
# plt.xlabel("Velocity (m/s)")
# plt.title("Matching point")
# plt.text(3800,38,"(3725 m/s,36.7 Pa)")

# # contourf(H_stagnation)
# # contourf(P_stagnation)
# # plt.legend
# # #print(P1)

 
# # rho2, u2, p2, h2, e2, gas = shock_calorically_imperfect(gas)

# # print(rho2, u2, p2, h2, e2) 
# # gas()
# # stagnation_conditions_isentropic(gas)





























"""
SCRATCHPAD BELOW
"""


#gas.
#print(gas.T)








# print(u2)
# 
# Pstagnation = gas.P+gas.density*(u2**2)
# print("pstagnation = "+ str(Pstagnation))

# print("enthalpy per kg" +str(gas.h))
# enthalpy_permol = gas.h*gas.mean_molecular_weight/1000
# print("enthalpy per mol"+str(enthalpy_permol))
# enthalpy_formation = -393474 #J/mol
# print("enthalpy formation"+str(enthalpy_formation))
# sensible_enthalpy = enthalpy_permol-enthalpy_formation
# print("sensible enthalpy J/mol = "+str(sensible_enthalpy))
# print("sensible enthalpy J/kg = ", str(sensible_enthalpy/(gas.mean_molecular_weight/1000)))

# hstagnation = 0.5*u2**2+sensible_enthalpy/(gas.mean_molecular_weight/1000)
# print("hstagnation = "+ str(hstagnation))
# print("gamma = " +str(gas.cp/gas.cv))







# E_total = gas.P/(gas.cp/gas.cv-1)
# H_ = (E+gas.P)/gas.density
# print("test  " + str(H))
