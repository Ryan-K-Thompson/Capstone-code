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
from scipy import interpolate
from scipy.optimize import fsolve, root


def plenum_to_nozzle_inlet(H_add, P0):
    # Plenum
    T0 = 300 # Plenum temperature, K
    gas = ct.ThermoPhase('airNASA9noions.cti') # Define the species properties
    gas.TPX = T0, P0, "CO2:1" # set the plenum temperature, pressure, and species
    h_plenum = gas.enthalpy_mass # determine the initial enthalpy of the gas

    # Mixing chamber
    h_mixing = h_plenum + H_add # add the enthalpy from the arc to the gas
    gas.HP = h_mixing, P0 # set the gas enthalpy and pressure in the mixing chamber
    gas.equilibrate('HP') # equilibrate the gas under constant enthalpy and pressure

    return gas

def soundspeed(gas):
    gamma = gas.cp / gas.cv
    return math.sqrt(gamma * ct.gas_constant
                     * gas.T / gas.mean_molecular_weight)

def isentropic(gas, mdot):
    # get the stagnation state parameters
    s0 = gas.s
    h0 = gas.h
    p0 = gas.P

    amin = 1.e14
    data = np.zeros((500, 10))

    # compute values for a range of pressure ratios
    for r in range(500):
        p = p0*(r+1)/501.0
        gas.SP = s0, p # reset the gas entropy and pressure

        v = (2.0*(h0 - gas.h))**0.5      # h + V^2/2 = h0
        area = mdot/(gas.density*v)    # rho*v*(A/A*) = constant
 
        amin = min(amin, area) # determine the throat area
        data[r, :] = [area, v/soundspeed(gas), gas.T, p/p0, gas.density_mass, gas.P, v, mdot, area, (area/(3.14))**0.5]

        
    data[:, 0] /= amin # divide the area by the throat area to determine the area ratio

    return data, amin

def Nozzle_selection(data, R_jet):
    data_design = []
    data_subsonic = []
    # filter data for only subsonic mach numbers
    for count, i in enumerate(data[:,1]):
        if data[count,1]<1:
            array = data[count,:]
            data_subsonic.append(array.tolist())
    data_subsonic = np.array(data_subsonic)
    # interpolate data to determine design data for given radius
    try:
        for r in range(10):
            value_r = float(scipy.interpolate.interp1d(data_subsonic[:,9],data_subsonic[:,r])(R_jet))
            data_design.append(value_r)
        # print(data_design)        
        
    except ValueError:
        print('value below interpolation range')
        
    return data_design

def stagnation_velocity_gradient(Ue,R_jet,R_model):
    L = R_model/R_jet
    beta = (Ue/R_model)*(1/(2-L-1.68*(L-1)**2-1.28*(L-1)**3))
    return beta 

def equivalent_flight_conditions(h0,P0,beta):
    v_flight = (2*h0)**0.5
    density_flight = P0/(v_flight**2)
    
    print(
        "Flight velocity = "+str(v_flight)+' m/s \n'
        "Flight density = "+str(density_flight)+' kg/m3 \n'
        )
    

if __name__ == "__main__":
    print(__doc__)
    mdot = 0.05
    P_arc = 240*10**3
    H_add = 0.5*P_arc/mdot
    P0 = 9*10**3
    R_jet = 0.05
    R_model = 0.0254
    gas = plenum_to_nozzle_inlet(H_add, P0)
    data, amin = isentropic(gas, mdot)
    print("H_add = " +str(H_add))
    
    data_subsonic = []
    count2 = 0
    
    data_design = Nozzle_selection(data, R_jet)
    v_design = data_design[6]
    beta = stagnation_velocity_gradient(v_design, R_jet, R_model)
    print(
        "Nozzle parameters \n"
        "Pressure = " + str(data_design[5]) + " Pa \n"
        "Temperature = " + str(data_design[2]) + " K \n"
        "Velocity = " + str(data_design[6]) + " m/s \n \n"
        "A/A* = " + str(data_design[0]) + " Pa \n""Pressure = " + str(data_design[5]) + " Pa \n"
        "P0/P = " + str(data_design[3]) + "  \n"
        "Outlet radius = " + str(0.05) + " m \n \n"
        "h0 = " + str(H_add/(10**6)) + " MJ/kg \n"
        "P0 = " + str(P0) + " MJ/kg \n"
        "beta = "+str(beta) +" 1/s \n" 
        )
    
    equivalent_flight_conditions(H_add,P0,beta)
    
    try:
        import matplotlib.pyplot as plt
        fig1 = plt.figure("Figure 1")
        plt.plot(data[:, 1], data[:, 0])
        plt.ylabel('Area Ratio')
        plt.xlabel('Mach Number')
        plt.title('Isentropic Flow: Area Ratio vs. Mach Number')
        plt.grid()
        # plt.xlim(0,1)
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
        plt.plot(data[:, 1], data[:, 8])
        plt.ylabel('A actual')
        plt.xlabel("M")
        plt.title('')
        plt.grid()

        plt.show()


        plt.show()
        
        fig3 = plt.figure("Figure 3")
        plt.plot(data[:, 1], data[:, 9])
        plt.ylabel('Radius')
        plt.xlabel('M')
        plt.title('')
        plt.grid()
        plt.show()
        
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