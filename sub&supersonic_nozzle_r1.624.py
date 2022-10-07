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
from scipy.interpolate import griddata
# from mpl_toolkits import Axes3D
from scipy.interpolate import RectBivariateSpline
from sympy import symbols, Eq, solve
from scipy import interpolate
from scipy.optimize import fsolve, root
from sdtoolbox import postshock, stagnation
from sdtoolbox.postshock import PostShock_eq, PostShock_fr
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
    for iter in range(len(beta_square)):
        for i in range(len(H_add)):
            for j in range(len(H_add)):
                if j == 0:
                    if beta_square[i][j] == 0:
                        beta_square[i][j] = beta_square[i][j+1]
                        H_add[i][j] = 0.999*H_add[i][j+1]
                        P0[i][j] = 0.999*P0[i][j+1]
                elif j == len(H_add)-1:
                    if beta_square[i][j] == 0:
                        beta_square[i][j] = beta_square[i][j-1]
                        H_add[i][j] = 1.001*H_add[i][j-1]
                        P0[i][j] = 1.001*P0[i][j-1]
                else:
                    if beta_square[i][j] == 0:
                        beta_square[i][j] = beta_square[i][j+1]
                        H_add[i][j] = 0.999*H_add[i][j+1]
                        P0[i][j] = 0.999*P0[i][j+1]
                        if beta_square[i][j] == 0:
                            beta_square[i][j] = beta_square[i][j-1]
                            H_add[i][j] = 1.001*H_add[i][j-1]
                            P0[i][j] = 1.001*P0[i][j-1]
                            
    print("nice wireframe test \\\\\\\\\\\\\\\\\\\\\5")            
    print(H_add)
    print("done")

def beta_curve_PWK(mdot_range,P_arc,P0_range,R_jet,R_model,resolution):
    number_of_points =resolution**2
    number_done = 0
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
                beta_square[count1][count2]=np.nan
            number_done = number_done + 1
            print(str(100*number_done/number_of_points)+"% through ground calcs")
    
    P0, mdot =np.meshgrid(P0_linspace,mdot_linspace)
    H_add = (0.5*P_arc/mdot)

    
    
    
    
    
    
    
    
    
    
    fig1 = plt.figure("1",figsize=(10,5))
    # nice_wireframe(H_add,P0,beta_square)
    ax = plt.axes(projection='3d')
    # ax.plot_wireframe(H_add/(10**6), P0/(10**3), beta_square)
    surf = ax.plot_surface(H_add/(10**6), P0/(10**3), beta_square/(10**3), cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    plt.gcf().set_size_inches(16, 8)
    ax.set_xlabel("H0 (MJ/kg)")
    ax.set_ylabel("P0 (kPa)")
    ax.set_zlabel("beta (1/s x10^3)")
    # plt.xlim((0,25))
    # plt.ylim((0,400))
    fig1.colorbar(surf, shrink=0.5, aspect=5)
    ax.view_init(20,50)
    plt.title("PWK plot (R_outlet="+str(R_jet)+'m, R_model='+str(R_model)+"m)")
    plt.savefig('test.png', dpi=300)
    ax.view_init(20,50)
    plt.savefig('test3.png', dpi=300)
    ax.view_init(0,0)
    plt.savefig('test4.png', dpi=300)
    ax.view_init(0,90)
    plt.savefig('test5.png', dpi=300)
    ax.view_init(90,0)
    plt.savefig('test6.png', dpi=300)

    return beta_square, P0, H_add, mdot_linspace, P0_linspace, H_add_linspace

def PWK_CFD_inputs_from_intersection(H0,P0,P_arc,R_jet,R_model):
    try:
        H_add = H0 
        mdot = 0.5*P_arc/H_add
        gas = plenum_to_nozzle_inlet(H_add, P0)
        data, amin = isentropic(gas, mdot)
        print("H_add = " +str(H_add))
        
        data_subsonic = []

        
        data_design = Nozzle_selection(data, R_jet)
        v_design = data_design[6]
        beta_pwk = stagnation_velocity_gradient(v_design, R_jet, R_model)
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
            "beta = "+str(beta_pwk) +" 1/s \n" 
            )
        gas()
        density_pwk = gas.density_mass
        temperature_pwk = data_design[2]
        pressure_pwk = data_design[3]*P0
        velocity_pwk = data_design[6]
        
        species = gas.species_names()
        mass_fractions = gas.Y
        print("Species: \n" + str(species))
        print("mass fractions: \n" + str(species))
    except:
        beta_pwk=np.nan
    
    
    
    return  [beta_pwk, density_pwk, temperature_pwk, pressure_pwk, velocity_pwk, mdot]


def beta_curve_flight(H0_linspace,P0_linspace,R_flight_vehicle,resolution):
    beta_square_flight = np.zeros((resolution,resolution)    )
    number_of_points =resolution**2
    number_done = 0
    for count1, h0 in enumerate(H0_linspace):
        for count2, P0 in enumerate(P0_linspace-1):
            try:
                mech = 'airNASA9noions.cti'
                q = 'CO2:1'    # assume CO2 atmosphere
                
                v_flight = (2*h0)**0.5
                density_flight = P0/(v_flight**2)
                altitude_flight = altitude_from_density_long2019(density_flight)
                pressure_flight = pressure_from_altitude_basicNASA(altitude_flight) 
                temperature_flight = temperature_from_altitude_basicNASA(altitude_flight)
                print(
                    "temperature = "+str(temperature_flight) +"\n"
                    'density = '+str(density_flight)+"\n"
                    "altitude = "+str(altitude_flight)
                    )
                # check for reasonable mars atmosphere
                
                if temperature_flight>[130]:
                    if [0.00001]<density_flight<[0.01]:
                        if altitude_flight>[3000]:
                    
                            gas_flight = ct.ThermoPhase(mech) 
                            gas_flight.TDY =  temperature_flight, density_flight, q
                            # temperature_flight = gas_flight.T
                            gas_flight.equilibrate('TV')
                            pressure_flight = gas_flight.P
                            print(pressure_flight)
            
                            # evaluate post shock conditions in flight 
                            gas_flight_postshock = PostShock_eq(U1=v_flight, P1=pressure_flight, T1=temperature_flight, q=q, mech=mech)
                            
                            # evaluate stagnation density in flight 
                            try:
                                output = stgsolve(gas = gas_flight_postshock, gas1 = gas_flight, U1 = v_flight, Delta = 0.5)
                                stagnation_density_flight = output['rho'][-1] 
                            
                                # radius_flight = v_flight/(beta/(math.sqrt((8/3)*(density_flight/stagnation_density_flight))))
                                beta_square_flight[count1][count2] = (v_flight/R_flight_vehicle)*math.sqrt((8/3)*(density_flight/stagnation_density_flight))
                            except:
                                beta_square_flight[count1][count2] = np.nan
                else:
                    beta_square_flight[count1][count2] = np.nan
            except RuntimeWarning:
                beta_square_flight[count1][count2] = np.nan
            number_done = number_done + 1
            print(str(100*number_done/number_of_points)+"% through flight calcs")
    
    P0, H0 =np.meshgrid(P0_linspace,H0_linspace)

    
    fig2 = plt.figure("2",figsize=(10,5))
    # nice_wireframe(H0,P0,beta_square_flight)
    ax = plt.axes(projection='3d')
    # ax.plot_wireframe(H_add/(10**6), P0/(10**3), beta_square)
    surf = ax.plot_surface(H0/(10**6), P0/(10**3), beta_square_flight/(10**3), cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    plt.gcf().set_size_inches(16, 8)
    ax.set_xlabel("H0 (MJ/kg)")
    ax.set_ylabel("P0 (kPa)")
    ax.set_zlabel("beta (1/s x10^3)")
    # plt.xlim((0,25))
    # plt.ylim((0,400))
    fig2.colorbar(surf, shrink=0.5, aspect=5)
    ax.view_init(20,50)
    plt.title("flight plot (R_flight_vehicle="+str(R_flight_vehicle)+'m')
    plt.savefig('test2.png', dpi=300)

    return beta_square_flight, P0, H0

def flight_CFD_inputs_from_intersection(h0,P0,R_flight_vehicle):
    mech = 'airNASA9noions.cti'
    q = 'CO2:1'    # assume CO2 atmosphere
    
    v_flight = (2*h0)**0.5
    density_flight = P0/(v_flight**2)
    altitude_flight = altitude_from_density_long2019(density_flight)
    pressure_flight = pressure_from_altitude_basicNASA(altitude_flight) 
    temperature_flight = temperature_from_altitude_basicNASA(altitude_flight)
    print(
        "temperature = "+str(temperature_flight) +"\n"
        'density = '+str(density_flight)+"\n"
        "altitude = "+str(altitude_flight)
        )
    # check for reasonable mars atmosphere
    
    if temperature_flight>[130]:
        if [0.00001]<density_flight<[0.01]:
            if altitude_flight>[3000]:
        
                gas_flight = ct.ThermoPhase(mech) 
                gas_flight.TDY =  temperature_flight, density_flight, q
                # temperature_flight = gas_flight.T
                gas_flight.equilibrate('TV')
                pressure_flight = gas_flight.P
                print(pressure_flight)

                # evaluate post shock conditions in flight 
                gas_flight_postshock = PostShock_eq(U1=v_flight, P1=pressure_flight, T1=temperature_flight, q=q, mech=mech)
                
                # evaluate stagnation density in flight 
                try:
                    output = stgsolve(gas = gas_flight_postshock, gas1 = gas_flight, U1 = v_flight, Delta = 0.5)
                    stagnation_density_flight = output['rho'][-1] 
                
                    # radius_flight = v_flight/(beta/(math.sqrt((8/3)*(density_flight/stagnation_density_flight))))
                    beta_flight = (v_flight/R_flight_vehicle)*math.sqrt((8/3)*(density_flight/stagnation_density_flight))
                except:
                    beta_flight = np.nan
    else:
        beta_flight = np.nan
    
    velocity_flight = math.sqrt(2*h0)
    
    return [beta_flight, density_flight, temperature_flight, pressure_flight, altitude_flight, velocity_flight]

def CFD_inputs_from_line(intercept_coordinates):
    flight_results =[]
    pwk_results=[]
    
    for i in range(len(intercept_coordinates)):
        values = intercept_coordinates[i]
        P0 = values[1]
        H0 = values[0]
    
        flight_results.append( flight_CFD_inputs_from_intersection(H0,P0,R_flight_vehicle=1.624) )
        pwk_results.append( PWK_CFD_inputs_from_intersection(H0,P0,P_arc=240*10**3,R_jet = 0.03, R_model = 0.0254))
        # results of form [beta_pwk, density_pwk, temperature_pwk, pressure_pwk, velocity_pwk]
    return flight_results, pwk_results

def ding(frequency):
    import winsound
    frequency = frequency  # Set Frequency To 2500 Hertz
    duration = 500  # Set Duration To 1000 ms == 1 second
    winsound.Beep(frequency, duration)

def compare_beta(resolution):
    beta_square_PWK, P0_grid_PWK, H0_grid_PWK, mdot_linspace, P0_linspace, H_add_linspace = beta_curve_PWK(mdot_range=[0.005,0.05], P_arc = 240*10**3, P0_range=[50000,500*10**3], R_jet = 0.03, R_model = 0.0254, resolution=resolution)
    P0_interpolated_grid_PWK, H0_interpolated_grid_PWK, beta_interpolated_grid_PWK = griddata_stack_solution(H0_grid_PWK,P0_grid_PWK,beta_square_PWK,"PWK interpolation")
    
    beta_square_flight, P0_grid_flight, H0_grid_flight = beta_curve_flight(H_add_linspace,P0_linspace,1.624,resolution)
    P0_interpolated_grid_flight, H0_interpolated_grid_flight, beta_interpolated_grid_flight = griddata_stack_solution(H0_grid_flight,P0_grid_flight,beta_square_flight,"Flight interpolation")
    
    
    fig3 = plt.figure("3",figsize=(10,5))

    ax = plt.axes(projection='3d')
    # ax.plot_wireframe(H_add/(10**6), P0/(10**3), beta_square)
    # surf = ax.plot_wireframe(H0_grid_PWK/(10**6), P0_grid_PWK/(10**3), beta_square_PWK/(10**3), color = "red")
    # surf = ax.plot_wireframe(H0_grid_flight/(10**6), P0_grid_flight/(10**3), beta_square_flight/(10**3), color = 'blue')
    surf = ax.plot_surface(H0_grid_PWK/(10**6), P0_grid_PWK/(10**3), beta_square_PWK/(10**3), cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    surf = ax.plot_surface(H0_grid_flight/(10**6), P0_grid_flight/(10**3), beta_square_flight/(10**3), cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    # plt.legend(["Plasma wind tunnel","Flight"])
    plt.gcf().set_size_inches(16, 8)
    ax.set_xlabel("H0 (MJ/kg)")
    ax.set_ylabel("P0 (kPa)")
    ax.set_zlabel("beta (1/s x10^3)")
    # plt.xlim((0,25))
    # plt.ylim((0,400))
    # fig2.colorbar(surf, shrink=0.5, aspect=5)
    # plt.title("combined")
    ax.view_init(20,50)
    plt.savefig('test3.png', dpi=300)
    ax.view_init(0,0)
    plt.savefig('test4.png', dpi=300)
    ax.view_init(0,90)
    plt.savefig('test5.png', dpi=300)
    ax.view_init(90,0)
    plt.savefig('test6.png', dpi=300)
    
    return beta_square_PWK, H0_grid_PWK, P0_grid_PWK, beta_square_flight, H0_grid_flight, P0_grid_flight

def interpolation_of_beta_curve(beta_square_PWK, H0_grid_PWK, P0_grid_PWK, beta_square_flight, H0_grid_flight, P0_grid_flight):
    for i in range(len(beta_square_PWK)):
        for j in range(len(beta_square_PWK)-1):
            
            # define x and ys
            x1_PWK = P0_grid_PWK[i][j]; x2_PWK = P0_grid_PWK[i][j+1]; 
            y1_PWK = beta_square_PWK[i][j]; y2_PWK = beta_square_PWK[i][j+1]; 
            x1_flight = P0_grid_flight[i][j]; x2_flight = P0_grid_flight[i][j+1]; 
            y1_flight = beta_square_flight[i][j]; y2_flight = beta_square_flight[i][j+1];  
            
            # check if there is intersection at this j
            if (y1_PWK > y1_flight and y2_PWK < y2_flight) or (y1_PWK < y1_flight and y2_PWK > y2_flight):
                # then there is intersection
                # define the interpolation functions for each line
                y_PWK_interpolate = scipy.interpolate.interp1d([x1_PWK,x2_PWK],[y1_PWK,y2_PWK])
                y_flight_interpolate = scipy.interpolate.interp1d([x1_flight,x2_flight],[y1_flight,y2_flight])
                
                # define intersection function, i.e. where the two lines intersect, the Delta(y) is 0
                def fun(x):
                    return y_PWK_interpolate(x)-y_flight_interpolate(x)
                try:
                    y_intercept = fsolve(fun, [(x2_PWK-x1_PWK)/2])
                    print(y_intercept)
                except:
                    print("intersection failed")
                
        
def bivariate_interpolation(H0_square,P0_square,beta_square,name):
    
    
    #restructure data to work properly, i.e. x strictly increasing
    beta_square_restructured = beta_square*0
    print(beta_square_restructured)
    P0_restructured = beta_square*0
    H0_restructured = beta_square*0
    for i in range(len(beta_square)):
        for ii in range(len(beta_square)):
            H0_restructured[i][ii]=H0_square[len(beta_square)-i-1][ii]
            P0_restructured[i][ii]=P0_square[len(beta_square)-i-1][ii]
            beta_square_restructured[i][ii]=beta_square[len(beta_square)-i-1][ii]
    
    nice_wireframe(H0_restructured, P0_restructured, beta_square_restructured)

    # remove beta = 0 collumns and rows from all datasets
    # for i in range(len(beta_square)):
    
    # for ii in range(len(beta_square)):
            
    
    def Extract(lst):
        return [item[0] for item in lst]
    
    P0_1d = P0_restructured[0]
    H0_1d = Extract(H0_restructured)
    
    print("H1D = "+str(H0_1d))
    print("P1D = "+str(P0_1d))
    print("beta = "+str(beta_square_restructured))
    interp_spline = RectBivariateSpline(P0_1d,H0_1d,beta_square_restructured)
    
    H0_interp = np.arange(H0_1d[0],H0_1d[-1],(H0_1d[-1]-H0_1d[0])/100).tolist()
    P0_interp = np.arange(P0_1d[0],P0_1d[-1],(P0_1d[-1]-P0_1d[0])/100).tolist()
    
    print(H0_interp)
    print(P0_interp)
    
    print(type(H0_interp))
    print(type(interp_spline))
    X2, Y2 = np.meshgrid(H0_interp,P0_interp)
    Z2 = interp_spline(P0_interp,H0_interp)
    
    ax = plt.axes(projection='3d')
    
    surf = ax.plot_surface(X2, Y2, Z2, cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    plt.gcf().set_size_inches(16, 8)
    ax.set_xlabel("H0 (MJ/kg)")
    ax.set_ylabel("P0 (kPa)")
    ax.set_zlabel("beta (1/s x10^3)")
    plt.title(name)
    print("X2 "+str(X2))
    print("Y2 "+str(Y2))
    print("Z2 "+str(Z2))
    
    
    # fig, ax = plt.subplots(nrows=1,ncols=1)
    # surf = ax.plot_wireframe(X2,Y2,Z2)
    plt.show()
    

def bivariate_interpolation2(H0_square,P0_square,beta_square,name):
    
    def Extract(lst):
        return [item[0] for item in lst]
    
    # P0_1d = P0_square[0]
    # H0_1d = Extract(H0_square)
    P0_1d = P0_square[0]
    H0_1d = Extract(H0_square)
    
    func =  scipy.interpolate.interp2d(P0_square, H0_square, beta_square, kind='cubic', copy=True, bounds_error=False, fill_value=1)
    
    P0_1d_new = np.arange(600,500000,(500000-600)/100)
    H0_1d_new = np.arange(2.4*10**6,25*10**6, (25*10**6-2.4*10**6)/100)
    
    P0_2d_new, H0_2d_new = np.meshgrid(P0_1d_new, H0_1d_new)
    
    
    figB = plt.figure("B",figsize=(10,5))
    ax = plt.axes(projection='3d')
    print( func(P0_1d_new,H0_1d_new))
    surf = ax.plot_surface(H0_2d_new,  P0_2d_new, func(P0_1d_new,H0_1d_new), cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    plt.gcf().set_size_inches(16, 8)
    ax.set_xlabel("H0 (MJ/kg)")
    ax.set_ylabel("P0 (kPa)")
    ax.set_zlabel("beta (1/s x10^3)")
    figB.colorbar(surf, shrink=0.5, aspect=5)
    plt.title(name)
    ax.view_init(20,50)
    plt.savefig('test3.png', dpi=300)
    ax.view_init(0,0)
    plt.savefig('test4.png', dpi=300)
    ax.view_init(0,90)
    plt.savefig('test5.png', dpi=300)
    ax.view_init(90,0)
    plt.savefig('test6.png', dpi=300)
    
    print(P0_square)
    print(H0_square)

    # grid_z = griddata(beta_square, (H0_square, P0_square), method = 'cubic')
    # ax = plt.axes(projection='3d')
    # surf = ax.plot_surface(H0_square,  P0_square, grid_interp, cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    # plt.gcf().set_size_inches(16, 8)
    # ax.set_xlabel("H0 (MJ/kg)")
    # ax.set_ylabel("P0 (kPa)")
    # ax.set_zlabel("beta (1/s x10^3)")
    # plt.title(name)
    
def griddata_stack_solution(H0_square,P0_square,beta_square,name):
    # ref https://stackoverflow.com/questions/34408293/2-d-interpolation-ignoring-nan-values
    x =   P0_square
    y =  H0_square
    z = beta_square
    
    for i in range(len(z)):
        for ii in range(len(z)):
            if z[i][ii] == 0:
                z[i][ii] = np.nan
    
    x=x.ravel()              #Flat input into 1d vector
    x=list(x[x!=np.isnan])   #eliminate any NaN
    y=y.ravel()
    y=list(y[y!=np.isnan])
    z=z.ravel()
    z=list(z[z!=np.isnan])

       
    xnew = np.arange(600,500000,(500000-600)/1000)
    ynew = np.arange(2.4*10**6,25*10**6, (25*10**6-2.4*10**6)/1000)
    znew = griddata((x, y), z, (xnew[None,:], ynew[:,None]), method='linear')

    
    figC = plt.figure("C")
    levels = np.linspace(500, 20000, 15)
    plt.ylabel('Y', size=15)
    plt.xlabel('X', size=15)
    cs = plt.contourf(xnew, ynew, znew, levels=levels, cmap=cm.coolwarm)
    cbar = plt.colorbar(cs)
    # cbar.set_label('Z', rotation=90, fontsize=15) # gas fraction
    plt.show()
    
    figD = plt.figure("D")
    ax = plt.axes(projection='3d')

    X_new, Y_new = np.meshgrid(xnew,ynew)

    surf = ax.plot_surface(X_new/1000,  Y_new/(10**6), znew/1000, cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
    ax.set_xlabel("P0 (kPa)")
    ax.set_ylabel("H0 (MJ/kg)")
    ax.set_zlabel("beta (1/s x10^3)")
    figD.colorbar(surf, shrink=0.5, aspect=5)
    plt.title(name)
    ax.view_init(20,50)
    plt.savefig('test3.png', dpi=300)
    ax.view_init(0,0)
    plt.savefig('test4.png', dpi=300)
    ax.view_init(0,90)
    plt.savefig('test5.png', dpi=300)
    ax.view_init(90,0)
    plt.savefig('test6.png', dpi=300)
    
    return X_new, Y_new, znew

def surface_intercept(x1, y1, z1, x2, y2, z2):
    z_intercept = z2-z1
    fig3 = plt.figure("Figure 3")
    cont = plt.contour(x1, y1, z_intercept,[0], colors = 'k')
    p1 = cont.collections[0].get_paths()[0] # grab the 1st path
    intercept_coordinates = p1.vertices
    print(intercept_coordinates)

    return  intercept_coordinates
    
def knudsen(T,P,R):
    Kn = (3.81*10**(-23))*T/( math.sqrt(2)*((3.3*10**(-10))**2)*P*R )
    return Kn

if __name__ == "__main__":
    # try:
        print(__doc__)
        #all_in_one_matching(mdot = 0.5, P_arc = 240*10**3, P0=1*10**3, R_jet=0.5, R_model=0.0245)
        
        # uncomment for matching surfaces
        beta_square_PWK, H0_grid_PWK, P0_grid_PWK, beta_square_flight, H0_grid_flight, P0_grid_flight = compare_beta(resolution=20)
        
        intercept_coordinates = surface_intercept(H0_grid_PWK, P0_grid_PWK, beta_square_PWK, H0_grid_flight, P0_grid_flight, beta_square_flight)
        
        
        flight_results, pwk_results = CFD_inputs_from_line(intercept_coordinates)
        
        intercept_H0 = []
        intercept_P0 = []
        intercept_beta = []
        for i in range(len(intercept_coordinates)):
            coords = intercept_coordinates[i]
            flight_result_array = flight_results[i]
            
            H0 = coords[0]/(10**6)
            P0 = coords[1]/(10**3)
            beta = flight_result_array[0]/(10**3)
            
            
            intercept_H0.append(H0)
            intercept_P0.append(P0)
            
            intercept_beta.append(beta)
        
        fig6 = plt.figure("6",figsize=(10,5))

        ax = plt.axes(projection='3d')
        
        surf = ax.plot_surface(H0_grid_flight/(10**6), P0_grid_flight/(10**3), beta_square_flight/(10**3), linewidth=1, antialiased=False, rcount=200, ccount=200, alpha=0.5)
        surf = ax.plot_surface(H0_grid_PWK/(10**6), P0_grid_PWK/(10**3), beta_square_PWK/(10**3), linewidth=1, antialiased=False, rcount=200, ccount=200, alpha=0.5)
        ax.plot3D(intercept_H0, intercept_P0, intercept_beta, color='k')
        ax.set_zlim(0, 5)
        # ax.legend(["flight", "PWK", "line"])
        # plt.gcf().set_size_inches(16, 8)
        ax.set_xlabel("H0 (MJ/kg)")
        ax.set_ylabel("P0 (kPa)")
        ax.set_zlabel("beta (1/s x10^3)")

        ax.view_init(20,50)
        plt.savefig('1.png', dpi=300)
        ax.view_init(0,0)
        plt.savefig('2.png', dpi=300)
        ax.view_init(0,90)
        plt.savefig('3.png', dpi=300)
        ax.view_init(90,0)
        plt.savefig('4.png', dpi=300)
        
        fig7 = plt.figure("7")
        ax = plt.axes(projection='3d')
        ax.plot3D(intercept_H0, intercept_P0, intercept_beta)
        plt.show()
        
        for i in range(len(flight_results)):
            result = flight_results[i]
            print(knudsen(result[2], result[3], 1.624))
        
        # figA = plt.figure("A",figsize=(10,5))
        # beta_square_PWK, P0_grid_PWK, H0_grid_PWK, mdot_linspace, P0_linspace, H_add_linspace = beta_curve_PWK(mdot_range=[0.005,0.05], P_arc = 240*10**3, P0_range=[10000,500*10**3], R_jet = 0.03, R_model = 0.0254, resolution=20)
        # ax = plt.axes(projection='3d')
        # surf = ax.plot_surface(H0_grid_PWK/(10**6), P0_grid_PWK/(10**3), beta_square_PWK/(10**3), cmap=cm.coolwarm, linewidth=1, antialiased=False, rcount=200, ccount=200)
        # plt.show()
        #interpolation_of_beta_curve(beta_square_PWK, H0_grid_PWK, P0_grid_PWK, beta_square_flight, H0_grid_flight, P0_grid_flight)
        
        # beta_square_PWK, P0_grid_PWK, H0_grid_PWK, mdot_linspace, P0_linspace, H_add_linspace = beta_curve_PWK(mdot_range=[0.005,0.05], P_arc = 240*10**3, P0_range=[600,500*10**3], R_jet = 0.03, R_model = 0.0254, resolution=10)
        # bivariate_interpolation2(H0_grid_PWK,P0_grid_PWK,beta_square_PWK,"PWK interpolation")
        # griddata_stack_solution(H0_grid_PWK,P0_grid_PWK,beta_square_PWK,"PWK interpolation")
        # bivariate_interpolation(H0_grid_PWK,P0_grid_PWK,beta_square_PWK,"PWK")
        # bivariate_interpolation(H0_grid_flight,P0_grid_flight,beta_square_flight,"flight")
        # line_intercept_of_two_surfs(H0_grid_PWK, P0_grid_PWK, beta_square_PWK, H0_grid_flight, P0_grid_flight, beta_square_flight)
        
        ding(2000)
        
    # except:
    #     ding(200)

