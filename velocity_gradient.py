# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 11:49:59 2022

@author: patag
"""

def Newtonian(U_farfield, R_nose):
    du_dy_ats = U_farfield/R_nose
    return du_dy_ats

def oneD_inviscid_momentum_and_newtonian_pressure(R_nose, P_farfield, P_stagnation, density_stagnation):
    du_dy_ats = (1/R_nose)*(2*(P_stagnation-P_farfield)/density_stagnation)**0.5
    return du_dy_ats

def Stokes_potential_flow(U_post_shock, R_nose, shock_standoff_distance):
    non_dim = shock_standoff_distance/R_nose
    du_dy_ats = (3/2)*(U_post_shock/R_nose)*( (1+non_dim)**3/((1+non_dim)**3-1) )
    return du_dy_ats

# def inviscid_conservation_eqn(U_farfield, R_nose, shock_standoff_distance, P_stagnation, P_postshock, P_farfield, U_farfield, density_postshock, density_stagnation):
#     du_dy_ats =
#     return du_dy_ats
    
def Barbante_and_Chazot(density_farfield, density_BLOE, U_farfield, R_nose):
    Beta = ((8/3)*(density_farfield/density_BLOE))*(U_farfield/R_nose)
    return Beta

Flight_nose_radius = 0.66 # m
Test_nose_radius = Flight_nose_radius/265 # m

