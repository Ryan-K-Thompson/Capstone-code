# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 10:19:14 2022

@author: patag
"""


import csv
from pandas import *
import matplotlib.pyplot as plt
import math


def extract_values(File, value):
    data = rea_csv(File)
    data = data.dropna
    values_array = data[value].tolist()
    
    return values_array




