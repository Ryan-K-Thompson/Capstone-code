# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 08:21:43 2022

@author: patag
"""
print('start')
import math
from mpl_toolkits import mplot3d
import numpy as np
import cantera as ct
import sympy as sympy
import scipy as scipy
from sympy import symbols, Eq, solve
from scipy import interpolate
from scipy.optimize import fsolve, root
from sdtoolbox import postshock, stagnation
from sdtoolbox.postshock import PostShock_eq
from sdtoolbox.stagnation import stgsolve
from matplotlib import pyplot as plt
from matplotlib import cm


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
    data = np.zeros((10000, 10))

    # compute values for a range of pressure ratios
    for r in range(10000):
        p = p0*(r+1)/10001.0
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

def altitude_from_density_long2019(density):
    def fun(alt):
        # function = (2*10**(-4))*math.exp((40000-alt)/7500)-density # Model based on Long et al. "Mars atmospheric entry guidance for optimal terminal altitude" 2019
        function = (0.699*math.exp(-0.00009*alt))/(0.1921*((-23.4-0.00222*alt)+273.1)) - density
        return function
    altitude = scipy.optimize.fsolve(fun,(10000))  
    altitude = list(altitude)
    return altitude

def pressure_from_altitude_basicNASA(altitude):
    # model is a bit sus from the NASA website https://www.grc.nasa.gov/WWW/K-12/rocket/atmosmrm.html
    P_alt = 699*math.exp(-0.00009*altitude[0])
    return P_alt

def temperature_from_altitude_basicNASA(altitude):
    # model is a bit sus from the NASA website https://www.grc.nasa.gov/WWW/K-12/rocket/atmosmrm.html
    if altitude[0] > 7000:
        T_alt = -23.4-0.00222*altitude[0]+273.15
    else:
        T_alt = -31-0.000998*altitude[0]+273.15    
    
    
    return T_alt

def equivalent_flight_conditions(h0,P0,beta):
    mech = 'airNASA9noions.cti'
    q = 'CO2:1'    # assume CO2 atmosphere
    
    v_flight = (2*h0)**0.5
    density_flight = P0/(v_flight**2)
    altitude_flight = altitude_from_density_long2019(density_flight)
    pressure_flight = pressure_from_altitude_basicNASA(altitude_flight) 
    temperature_flight = temperature_from_altitude_basicNASA(altitude_flight)
    
    gas_flight = ct.ThermoPhase(mech) 
    gas_flight.TDY =  temperature_flight, density_flight, q
    # temperature_flight = gas_flight.T
    gas_flight.equilibrate('TV')
    pressure_flight = gas_flight.P
    print(pressure_flight)


    
    # gas_flight()
    
    # evaluate post shock conditions in flight 
    gas_flight_postshock = PostShock_eq(U1=v_flight, P1=pressure_flight, T1=temperature_flight, q=q, mech=mech)
    
    # evaluate stagnation density in flight 
    output = stgsolve(gas = gas_flight_postshock, gas1 = gas_flight, U1 = v_flight, Delta = 0.2)
    stagnation_density_flight = output['rho'][-1] 
    
    radius_flight = v_flight/(beta/(math.sqrt((8/3)*(density_flight/stagnation_density_flight))))
    
    
    
    print(
        "Flight velocity = "+str(v_flight)+' m/s \n'
        "Flight density = "+str(density_flight)+' kg/m3 \n'
        "Flight pressure = "+str(pressure_flight)+' Pa \n'
        "Flight temperature = "+str(temperature_flight)+' K \n'
        "Flight altitude = "+str(altitude_flight[0]/1000)+' km \n'
        "Flight radius (target 1.265m) = "+str(stagnation_density_flight)+' m \n'
        )
    return output
    
def all_in_one_matching(mdot,P_arc,P0,R_jet,R_model):
    H_add = 0.5*P_arc/mdot
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
        "P0/P = " + str(1/data_design[3]) + "  \n"
        "Outlet radius = " + str(0.05) + " m \n \n"
        "h0 = " + str(H_add/(10**6)) + " MJ/kg \n"
        "P0 = " + str(P0) + " Pa \n"
        "beta = "+str(beta) +" 1/s \n" 
        )
    
    output = equivalent_flight_conditions(H_add,P0,beta)
    
    try:
        import matplotlib.pyplot as plt
        fig1 = plt.figure("Figure 1")
        plt.plot(data[:, 1], data[:, 0])
        plt.ylabel('Area Ratio')
        plt.xlabel('Mach Number')
        plt.title('Isentropic Flow: Area Ratio vs. Mach Number')
        plt.grid()
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

    except ImportError:
        print('area ratio,   Mach number,   temperature,   pressure ratio')
        print(data)

def nice_wireframe(H_add,P0,beta_square):
    for iter in range(10):
        for i in range(len(H_add)):
            for j in range(len(H_add)):
                if j == 0:
                    if beta_square[i][j] == 0:
                        beta_square[i][j] = beta_square[i][j+1]
                        H_add[i][j] = H_add[i][j+1]
                        P0[i][j] = P0[i][j+1]
                elif j == len(H_add)-1:
                    if beta_square[i][j] == 0:
                        beta_square[i][j] = beta_square[i][j-1]
                        H_add[i][j] = H_add[i][j-1]
                        P0[i][j] = P0[i][j-1]
                else:
                    if beta_square[i][j] == 0:
                        beta_square[i][j] = beta_square[i][j+1]
                        H_add[i][j] = H_add[i][j+1]
                        P0[i][j] = P0[i][j+1]
                        if beta_square[i][j] == 0:
                            beta_square[i][j] = beta_square[i][j-1]
                            H_add[i][j] = H_add[i][j-1]
                            P0[i][j] = P0[i][j-1]
                            
                
   
    print("done")

def beta_curve_PWK(mdot_range,P_arc,P0_range,R_jet,R_model,resolution):
    H_add_linspace = np.linspace(0.5*P_arc/mdot_range[0],0.5*P_arc/mdot_range[1],resolution)
    mdot_linspace = 0.5*P_arc/H_add_linspace
    P0_linspace = np.linspace(P0_range[0],P0_range[1],resolution)
    beta_square = np.zeros((resolution,resolution)    )

    for count1, mdot in enumerate(mdot_linspace):
        for count2, P0 in enumerate(P0_linspace):
            try:
                H_add = 0.5*P_arc/mdot
                gas = plenum_to_nozzle_inlet(H_add, P0)
                data, amin = isentropic(gas, mdot)
                print("H_add = " +str(H_add))
                
                data_subsonic = []

                
                data_design = Nozzle_selection(data, R_jet)
                v_design = data_design[6]
                beta = stagnation_velocity_gradient(v_design, R_jet, R_model)
                print(
                    "Nozzle parameters \n"
                    "Pressure = " + str(data_design[5]) + " Pa \n"
                    "Temperature = " + str(data_design[2]) + " K \n"
                    "Velocity = " + str(data_design[6]) + " m/s \n \n"
                    "A/A* = " + str(data_design[0]) + " Pa \n""Pressure = " + str(data_design[5]) + " Pa \n"
                    "P0/P = " + str(1/data_design[3]) + "  \n"
                    "Outlet radius = " + str(0.05) + " m \n \n"
                    "h0 = " + str(H_add/(10**6)) + " MJ/kg \n"
                    "P0 = " + str(P0) + " Pa \n"
                    "beta = "+str(beta) +" 1/s \n" 
                    )
                beta_square[count1][count2]=beta

            except:
                beta_square[count1][count2]=0
    
    P0, mdot =np.meshgrid(P0_linspace,mdot_linspace)
    H_add = (0.5*P_arc/mdot)
    print(P0)
    
    fig3 = plt.figure("3",figsize=(10,5))
    nice_wireframe(H_add,P0,beta_square)
    ax = plt.axes(projection='3d')
    # ax.plot_wireframe(H_add/(10**6), P0/(10**3), beta_square)
    surf = ax.plot_surface(H_add/(10**6), P0/(10**3), beta_square/(10**3), cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    plt.gcf().set_size_inches(16, 8)
    ax.set_xlabel("H0 (MJ/kg)")
    ax.set_ylabel("P0 (kPa)")
    ax.set_zlabel("beta (1/s x10^3)")
    # plt.xlim((0,25))
    # plt.ylim((0,400))
    fig3.colorbar(surf, shrink=0.5, aspect=5)
    ax.view_init(20,50)
    plt.title("PWK plot (R_outlet="+str(R_jet)+'mm, R_model='+str(R_model)+"mm)")
    plt.savefig('test.png', dpi=300)

    return beta_square, mdot_linspace, P0_linspace, H_add_linspace


def beta_curve_flight(H0_linspace,P0_range,R_flight_vehicle,resolution):
    H_add_linspace = np.linspace(0.5*P_arc/mdot_range[0],0.5*P_arc/mdot_range[1],resolution)
    mdot_linspace = 0.5*P_arc/H_add_linspace
    P0_linspace = np.linspace(P0_range[0],P0_range[1],resolution)
    beta_square = np.zeros((resolution,resolution))
    beta_max=[0]
    for count1, mdot in enumerate(H0_linspace):
        for count2, P0 in enumerate(P0_linspace):
            try:
                mech = 'airNASA9noions.cti'
                q = 'CO2:1'    # assume CO2 atmosphere
                
                v_flight = (2*h0)**0.5
                density_flight = P0/(v_flight**2)
                altitude_flight = altitude_from_density_long2019(density_flight)
                pressure_flight = pressure_from_altitude_basicNASA(altitude_flight) 
                temperature_flight = temperature_from_altitude_basicNASA(altitude_flight)
                
                gas_flight = ct.ThermoPhase(mech) 
                gas_flight.TDY =  temperature_flight, density_flight, q
                # temperature_flight = gas_flight.T
                gas_flight.equilibrate('TV')
                pressure_flight = gas_flight.P
                print(pressure_flight)
    
    
                
                # gas_flight()
                
                # evaluate post shock conditions in flight 
                gas_flight_postshock = PostShock_eq(U1=v_flight, P1=pressure_flight, T1=temperature_flight, q=q, mech=mech)
                
                # evaluate stagnation density in flight 
                output = stgsolve(gas = gas_flight_postshock, gas1 = gas_flight, U1 = v_flight, Delta = 0.2)
                stagnation_density_flight = output['rho'][-1] 
                
                radius_flight = v_flight/(beta/(math.sqrt((8/3)*(density_flight/stagnation_density_flight))))
                
                
    #             H_add = 0.5*P_arc/mdot
    #             gas = plenum_to_nozzle_inlet(H_add, P0)
    #             data, amin = isentropic(gas, mdot)
    #             print("H_add = " +str(H_add))
                
    #             data_subsonic = []

                
    #             data_design = Nozzle_selection(data, R_jet)
    #             v_design = data_design[6]
    #             beta = stagnation_velocity_gradient(v_design, R_jet, R_model)
    #             print(
    #                 "Nozzle parameters \n"
    #                 "Pressure = " + str(data_design[5]) + " Pa \n"
    #                 "Temperature = " + str(data_design[2]) + " K \n"
    #                 "Velocity = " + str(data_design[6]) + " m/s \n \n"
    #                 "A/A* = " + str(data_design[0]) + " Pa \n""Pressure = " + str(data_design[5]) + " Pa \n"
    #                 "P0/P = " + str(1/data_design[3]) + "  \n"
    #                 "Outlet radius = " + str(0.05) + " m \n \n"
    #                 "h0 = " + str(H_add/(10**6)) + " MJ/kg \n"
    #                 "P0 = " + str(P0) + " Pa \n"
    #                 "beta = "+str(beta) +" 1/s \n" 
    #                 )
    #             beta_square[count1][count2]=beta
    #             # beta_max = max(beta_max,beta) 
    #             # print(beta_max)
    #         except:
    #             beta_square[count1][count2]=0
    
    # P0, mdot =np.meshgrid(P0_linspace,mdot_linspace)
    # H_add = (0.5*P_arc/mdot)
    # print(P0)
    
    # fig3 = plt.figure("3",figsize=(10,5))
    # nice_wireframe(H_add,P0,beta_square)
    # ax = plt.axes(projection='3d')
    # # ax.plot_wireframe(H_add/(10**6), P0/(10**3), beta_square)
    # surf = ax.plot_surface(H_add/(10**6), P0/(10**3), beta_square/(10**3), cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    # plt.gcf().set_size_inches(16, 8)
    # ax.set_xlabel("H0 (MJ/kg)")
    # ax.set_ylabel("P0 (kPa)")
    # ax.set_zlabel("beta (1/s x10^3)")
    # # plt.xlim((0,25))
    # # plt.ylim((0,400))
    # fig3.colorbar(surf, shrink=0.5, aspect=5)
    # ax.view_init(20,50)
    # plt.title("PWK plot (R_outlet="+str(R_jet)+'mm, R_model='+str(R_model)+"mm)")
    # plt.savefig('test.png', dpi=300)

    return beta_square, mdot_linspace, P0_linspace

def ding():
    import winsound
    frequency = 3500  # Set Frequency To 2500 Hertz
    duration = 100  # Set Duration To 1000 ms == 1 second
    winsound.Beep(frequency, duration)

def compare_beta(resolution):
    beta_square, mdot_linspace, P0_linspace, H_add_linspace = beta_curve_PWK(mdot_range=[0.005,0.05], P_arc = 240*10**3, P0_range=[600,400*10**3], R_jet = 0.03, R_model = 0.01, resolution = 10)
    beta_curve_flight(H0_linspace,P0_range,R_flight_vehicle,resolution)
    
    

if __name__ == "__main__":
    print(__doc__)
    #all_in_one_matching(mdot = 0.5, P_arc = 240*10**3, P0=1*10**3, R_jet=0.5, R_model=0.0245)
    compare_beta(resolution=10)
    ding()

