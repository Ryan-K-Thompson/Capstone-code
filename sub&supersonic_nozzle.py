# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 08:21:43 2022

@author: patag
"""
print('start')
import cantera as ct

from sympy import symbols, Eq, solve
# from scipy.optimize import fsolve

print('imported')
H_add = 1*10**6
Pratio = 1.9
Aout_Athroat = 1.2

# pressure chamber
#U = 200 #m/s
P0 = 17500 #Pa
Ptest = P0/Pratio
print(Ptest)
T0 = 300 # K

gas0 = ct.ThermoPhase('airNASA9noions.cti')

gas0.TPX =T0, P0, "CO2:1"
H0 = gas0.enthalpy_mass

#mixing chamber
gas0.HP = H0 + H_add, P0
gas0.equilibrate('HP')
print("temp_mix = "+str(gas0.T))


gamma = gas0.cp_mass/gas0.cv_mass
gasconstant = 8.3144598 
R = gasconstant/gas0.mean_molecular_weight*1000
a = (gamma*R*gas0.T)**0.5
print("gamma = "+str(gamma))

# result = fsolve(lambda M: (1/M)*( (2/(gamma+1)) * (1+((gamma-1)/2)*(M**2)))**((1+gamma)/(2*(gamma-1))), Aout_Athroat )
# print(result)
M = symbols("M")
eq1 = Eq( (1/M)*( (2/(gamma+1)) * (1+((gamma-1)/2)*(M**2)))**((1+gamma)/(2*(gamma-1))) - Aout_Athroat, 0)

M_array = solve(eq1, M)
print(M_array)

# # M_out = (((Pratio)**((gamma-1)/gamma)-1)*(2/(gamma-1)))**0.5
# # print("Mout = "+str(M_out))

# # Aout_Athroat = (1/M_out)*( (2/(gamma+1)) * (1+((gamma-1)/2)*(M_out**2)))**((1+gamma)/(2*(gamma-1)))
# # print(Aout_Athroat)


# gas0.SP = gas0.entropy_mass, Ptest
# gamma = gas0.cp_mass/gas0.cv_mass
# gasconstant = 8.3144598 
# R = gasconstant/gas0.mean_molecular_weight*1000
# a = (gamma*R*gas0.T)**0.5
# Uout = M_out*a
# print("Uout = "+str(Uout))
# print("tempout = " +str(gas0.T))