# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 12:32:32 2022

@author: patag
"""

import numpy as np
import csv

import sys

sys.path.append(r'D:\Files\Uni\2022\Capstone\Capstone-code')

'''
Flight data inputs
'''
with open("pathfinder_velocity.csv") as velocity_csv:
    fight_velocity_list = list(csv.reader(velocity_csv))
    
with open("pathfinder_density.csv") as density_csv:
    fight_density1_list = list(csv.reader(density_csv,delimiter = ','))
    print(fight_density1_list)





# flight_density_e_list = [1, 2]
# flight_density_ratio_array = [i/j for i, j in zip(flight_density_1_list, (flight_density_e_list))]

RN = 0.66   # flight nose radius, m
RT_star = np.linspace(0.005,0.025,1)