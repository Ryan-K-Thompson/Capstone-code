# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 08:21:43 2022

@author: patag
"""
print('start')
import math
import numpy as np
import cantera as ct
import sympy as sympy
import scipy as scipy
from sympy import symbols, Eq, solve
from scipy.optimize import fsolve, root


def plenum_to_nozzle_inlet(H_add, P0):
        
    # print('imported')
    
    # Pratio = 1.9
    # Aout_Athroat = 1.2
    
    # pressure chamber
    # U = 200 #m/s
    # P0 = 90*1000 #Pa
    # Ptest = P0/Pratio
    # print(Ptest)
    T0 = 300 # K
    
    gas = ct.ThermoPhase('airNASA9noions.cti')
    
    gas.TPX =T0, P0, "CO2:1"
    H0 = gas.enthalpy_mass
    # print("H0=" +str(H0))
    
    #mixing chamber
    # gas()
    gas.HP = H0 + H_add, P0
    gas.equilibrate('HP')
    # print("temp_mix = "+str(gas0.T))
    
    
    # gamma = gas.cp_mass/gas.cv_mass
    # gasconstant = 8.3144598 
    # R = gasconstant/gas0.mean_molecular_weight*1000
    # a = (gamma*R*gas0.T)**0.5
    # # print("gamma = "+str(gamma))
    
    return gas



def soundspeed(gas):
    """The speed of sound. Assumes an ideal gas."""

    gamma = gas.cp / gas.cv
    return math.sqrt(gamma * ct.gas_constant
                     * gas.T / gas.mean_molecular_weight)


def isentropic(gas):
    """
    In this example, the area ratio vs. Mach number curve is computed. 
    """
    
    # get the stagnation state parameters
    s0 = gas.s
    h0 = gas.h
    p0 = gas.P

    mdot = 0.5  # arbitrary
    #arearatio = 1.2
    amin = 1.e14

    data = np.zeros((200, 10))
    
    gas()
    # compute values for a range of pressure ratios
    for r in range(200):

        p = p0*(r+1)/201.0
        # set the state using (p,s0)
        gas.SP = s0, p

        v = (2.0*(h0 - gas.h))**0.5      # h + V^2/2 = h0
        arearatio = mdot/(gas.density*v)    # rho*v*(A/A*) = constant
        A_throat_choked = mdot/(soundspeed(gas)*gas.density)

        amin = min(amin, arearatio)
        data[r, :] = [arearatio, v/soundspeed(gas), gas.T, p/p0, gas.density_mass, gas.P, v, mdot, A_throat_choked, arearatio]

        gamma = gas.cp_mass/gas.cv_mass
        print(gamma)
        print(gas.T)
        
    data[:, 0] /= amin
    # for r in range(200):
    #     actual_area = arearatio*A_throat_choked
        # data[r, 9] = actual_area
    return data


if __name__ == "__main__":
    print(__doc__)
    
    H_add = 2*10**6
    P0 = 90*10**3
   
    
    gas = plenum_to_nozzle_inlet(H_add, P0)
    data = isentropic(gas)
    try:
        import matplotlib.pyplot as plt
        fig1 = plt.figure("Figure 1")
        plt.plot(data[:, 1], data[:, 0])
        plt.ylabel('Area Ratio')
        plt.xlabel('Mach Number')
        plt.title('Isentropic Flow: Area Ratio vs. Mach Number')
        plt.grid()
        plt.xlim(0,1)
        plt.show()
        

        fig6 = plt.figure("Figure 6")
        plt.plot(data[:, 1], 1/data[:, 3])
        plt.ylabel('P0/P')
        plt.xlabel('M')
        plt.title('Isentropic Flow: P0/P vs mdot')
        plt.xlim(0,1)
        plt.ylim(1,2)
        plt.grid()
        plt.show()
        
        fig2 = plt.figure("Figure 2")
        plt.plot(data[:, 1], data[:, 9])
        plt.ylabel('A actual')
        plt.xlabel("M")
        plt.title('')
        plt.grid()

        plt.show()
        
        fig3 = plt.figure("Figure 3")
        plt.plot(data[:, 1], data[:, 8])
        plt.ylabel('A throat')
        plt.xlabel("M")
        plt.title('')
        plt.grid()

        # plt.show()
        
        # fig3 = plt.figure("Figure 3")
        # plt.plot(data[:, 0], 1/data[:, 3])
        # plt.ylabel('P0/P')
        # plt.xlabel('Area Ratio')
        # plt.title('Isentropic Flow: Area Ratio vs. P0/P')
        # plt.grid()
        # plt.show()
        
        # fig4 = plt.figure("Figure 4")
        # plt.plot(data[:, 1], 1/data[:, 3])
        # plt.ylabel('P0/P')
        # plt.xlabel('M')
        # plt.title('Isentropic Flow: M vs. P0/P')
        # plt.grid()
        # plt.show()

        # fig5 = plt.figure("Figure 5")
        # plt.plot(1/data[:, 3], data[:, 7])
        # plt.ylabel('mdot')
        # plt.xlabel('P0/P')
        # plt.title('Isentropic Flow: P0/P vs mdot')
        # plt.grid()
        # plt.show()
        
    except ImportError:
        print('area ratio,   Mach number,   temperature,   pressure ratio')
        print(data)


# func = lambda M : (Aout_Athroat**2)/((1/(M**2))*((2/(gamma+1))*(1+((gamma-1)/2)*(M**2) ))**((gamma+1)/(gamma-1))) 

# M_initial_guess = 0.58
# M_subsonic_solution = fsolve(func, M_initial_guess)
# print("M="+str(M_subsonic_solution))

# def fun(M):
#     return (Aout_Athroat**2)-((1/(M**2))*((2/(gamma+1))*(1+((gamma-1)/2)*(M**2) ))**((gamma+1)/(gamma-1))) 
# M = root(fun,[0.1,10])
# Msub = M.x[0]
# Msuper = M.x[1]

# P0_P_requiredforchoked = (1+((gamma-1)/2)*Msub**2)**(gamma/(gamma-1))

# if P

# result = scipy.optimize.fsolve(lambda M: (1/M)*( (2/(gamma+1)) * (1+((gamma-1)/2)*(M**2)))**((1+gamma)/(2*(gamma-1))), Aout_Athroat )
# print(result)
# M = symbols("M")
# eq1 = Eq( (1/M)*( (2/(gamma+1)) * (1+((gamma-1)/2)*(M**2)))**((1+gamma)/(2*(gamma-1))) - Aout_Athroat, 0)

# M_array = solve(eq1, M)
# print(M_array)

# # M_out = (((Pratio)**((gamma-1)/gamma)-1)*(2/(gamma-1)))**0.5
# # print("Mout = "+str(M_out))

# # Aout_Athroat = (1/M_out)*( (2/(gamma+1)) * (1+((gamma-1)/2)*(M_out**2)))**((1+gamma)/(2*(gamma-1)))
# # print(Aout_Athroat)


# # gas0.SP = gas0.entropy_mass, Ptest
# # gamma = gas0.cp_mass/gas0.cv_mass
# # gasconstant = 8.3144598 
# # R = gasconstant/gas0.mean_molecular_weight*1000
# # a = (gamma*R*gas0.T)**0.5
# # Uout = M_out*a
# # print("Uout = "+str(Uout))
# # print("tempout = " +str(gas0.T))