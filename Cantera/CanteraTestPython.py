# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 15:30:29 2022

@author: patag
"""

import cantera as ct
import numpy as np

#import sys
#sys.path.append("C:\Users\patag\Desktop\Uni\2022\Capstone\Cantera")

gas = ct.Solution('CO2_N2.cti') # https://shepherd.caltech.edu/EDL/PublicResources/sdt/cti_mech.html
gas.TPX = None, None, 'CO2:0.95, N2:0.05'

gas.equilibriate('TP')