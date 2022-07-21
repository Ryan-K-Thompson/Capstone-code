# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:35:24 2022

@author: patag
"""

import sys
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Cantera\SDToolbox\Python3\sdtoolbox')
import cantera as ct
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.pyplot import contour, contourf


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