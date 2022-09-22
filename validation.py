# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 15:56:05 2022

@author: patag
"""

import csv
from pandas import *
import matplotlib.pyplot as plt
import math

# LOAD DATA FROM WALL SPREADSHEET
wall = read_csv("wall.csv")
wall=wall.dropna()
x_wall = wall['x'].tolist()
xpR_wall = [i/(0.0254) for i in x_wall]
y_wall = wall['y'].tolist()
s_wall = []
for count,x in enumerate(x_wall):
    if count == 0:
        s_wall.append(0)
    else:
        s = ((x_wall[count]-x_wall[count-1])**2+(y_wall[count]-y_wall[count-1])**2)**0.5 + s_wall[count - 1]
        s_wall.append(s)
        
spR_wall = [i/(0.0254) for i in s_wall]
Qdot_Wpm2_wall = wall['Qdot'].tolist()
Qdot_Wpcm2_wall = [i/(100*100) for i in Qdot_Wpm2_wall]
cell_Re_wall = wall['cell_Re'].tolist()
wall.dropna()

# LOAD DATA FROM SYMM SPREADSHEET
symm = read_csv("symm.csv")
symm=symm.dropna()
x_symm = symm['x'].tolist()
xpR_symm = [i/(0.0254) for i in x_symm]
dVdy_symm = symm['dVdy'].tolist()
M_symm = symm['M'].tolist()
Kn_symm = symm['Kn'].tolist()
gamma_symm = symm['gamma'].tolist()
P_symm = symm['P'].tolist()
T_symm = symm['T'].tolist()
CO2_symm = symm['CO2'].tolist()
CO_symm = symm['CO'].tolist()
C2_symm = symm['C2'].tolist()
O2_symm = symm['O2'].tolist()
C_symm = symm['C'].tolist()
O_symm = symm['O'].tolist()
symm.dropna()

# LOAD EXPERIMENTAL DATA
spR_wall_experimental = [0.009708738,0.097087379,0.203883495,0.300970874,0.398058252,0.495145631,0.601941748,0.699029126,0.796116505,0.902912621,1,1.077669903]
Qdot_Wpcm2_wall_experimental = [824.9132,730.4700058,584.6862342,528.3144731,488.4990034,453.6504211,441.9726098,405.4683984,405.388028,392.0545875,398.5967337,29.32713946]
experimental_error = []
for Qdot in Qdot_Wpcm2_wall_experimental:
    experimental_error.append(Qdot*0.11)

# LOAD CFD COMPARISON DATA (TEMPS)
x_cfdcomp_t_vibrationalelectric = [-0.03919218,-0.037892936,-0.036593693,-0.035294449,-0.033995206,-0.032695962,-0.031396719,-0.030097475,-0.029034458,-0.028384836,-0.028030497,-0.027912384,-0.027774585,-0.027558045,-0.027380875,-0.027321819,-0.027262762,-0.02718402,-0.027085593,-0.027026536,-0.026947794,-0.026849366,-0.02679031,-0.026672197,-0.026593455,-0.026435971,-0.026435971,-0.026357229,-0.026239116,-0.026140688,-0.02600289,-0.025806035,-0.025550123,-0.024782388,-0.023955597,-0.023601258,-0.023246918,-0.022833523,-0.022361071,-0.021829562,-0.02117994,-0.020353149,-0.019231075,-0.017931831,-0.016632588,-0.015321533,-0.014034101,-0.012734857,-0.011435614,-0.01013637,-0.008837127,-0.007537883,-0.00623864,-0.004939396,-0.003817322,-0.003049587,-0.002518079,-0.002104683,-0.0018094,-0.001573174,-0.001336948,-0.001100722,-0.000884181,-0.000687326,-0.0004511,-0.000273931,-0.000155818,0]
T_cfdcomp_t_vibrationalelectric =[1062.969434,1048.344726,1047.805956,1047.267186,1046.728417,1046.189647,1045.650877,1054.502732,1120.914228,1349.988504,1534.915138,1715.635693,1973.820743,2292.229651,2472.925717,2662.278835,2868.848099,3069.670485,3376.650942,3591.828279,3792.650664,4082.414975,4280.376166,4529.961306,4842.688641,5175.468828,4986.09122,5364.813782,5663.178003,5993.113321,6277.12259,6500.850858,6793.419221,7010.65034,6781.280573,6568.801167,6356.32176,6112.25493,5883.945077,5672.826881,5480.310531,5287.003373,5097.356103,4924.655871,4792.648349,4683.642961,4592.802577,4503.052868,4428.954201,4351.725326,4290.147492,4227.004554,4159.166304,4052.200449,3894.246538,3695.51209,3486.976316,3294.557924,3107.927227,2918.451661,2737.584169,2522.284384,2243.866893,1984.10839,1665.511729,1377.067812,1144.60086,860.4936333]
x_cfdcomp_t_transrot = [-0.029861249,-0.02932974,-0.029093514,-0.028857288,-0.028798232,-0.028680119,-0.028562006,-0.028562006,-0.028443893,-0.028286409,-0.028207666,-0.028089553,-0.028089553,-0.027912384,-0.02797144,-0.02797144,-0.027853327,-0.027853327,-0.027853327,-0.027735214,-0.027735214,-0.027735214,-0.027735214,-0.027558045,-0.027617101,-0.027439932,-0.027498988,-0.027498988,-0.027380875,-0.027380875,-0.027203706,-0.027262762,-0.027144649,-0.027144649,-0.025727292,-0.025609179,-0.025491066,-0.025372953,-0.02525484,-0.025136727,-0.024959558,-0.024723332,-0.024487105,-0.024250879,-0.024014653,-0.023778427,-0.023483145,-0.023128805,-0.022685881,-0.022006731,-0.02117994,-0.020084114,-0.018782245,-0.017423945,-0.016293013,-0.014949477,-0.01384512,-0.012734857,-0.011064401,-0.009368635,-0.008069392,-0.006356753,-0.000214874,-9.6761E-05,-0.026984353,-0.02688478,-0.026685633,-0.026088193,-0.025889047]
T_cfdcomp_t_transrot = [1145.851576,1389.116666,1587.004388,1793.500183,1974.245228,2249.654587,2697.225408,2507.8478,2964.026694,3420.189261,3730.047239,4108.753474,3919.375867,4952.271166,4694.053464,4418.595126,5554.811791,5348.218038,5141.624284,6343.262305,6140.111781,5933.518027,5744.14042,6957.805254,6604.898748,7560.321389,7345.144052,7138.550299,7964.876334,7758.28258,8438.246884,8171.421108,8825.585682,8618.991929,9005.767468,8764.692443,8523.617418,8342.798905,8153.372318,7946.729585,7705.630071,7378.425336,7102.86904,6835.920817,6594.796813,6379.497028,6172.780826,5966.040136,5776.478857,5552.387327,5347.602738,5136.489644,4962.258,4825.687187,4721.490931,4619.788936,4547.023168,4463.925263,4374.692379,4290.367898,4220.964544,4139.912332,1032.7204,610.8758407,9180.408657,9702.863738,10312.36026,9615.450681,9267.037185]

# LOAD CFD COMPARISON DATA (CHEMISTRY)

fig1 = plt.figure("Figure 1")
plt.scatter(spR_wall,Qdot_Wpcm2_wall, color = 'blue')
plt.errorbar(spR_wall_experimental,Qdot_Wpcm2_wall_experimental, yerr=experimental_error, fmt="s", color = 'black')
plt.grid()
plt.xlabel('S/R')
plt.ylabel('Heat transfer (W/cm2)')
plt.title('Experimental validation')
plt.xlim(0, 1.2)
plt.savefig('plt_Experimental_validation.png', dpi=300)

fig2 = plt.figure("Figure 2")
plt.scatter(spR_wall,cell_Re_wall, color = 'blue')
plt.grid()
plt.xlabel('S/R')
plt.ylabel('Cell Reynolds number')
plt.title('Mesh sizing - Cell Re')
plt.xlim(0, 1.2)
plt.savefig('plt_Cell_Reynolds_number.png', dpi=300)

fig3 = plt.figure("Figure 3")
plt.scatter(xpR_symm,M_symm, color = 'blue')
plt.grid()
plt.xlabel('X/R')
plt.ylabel('M')
plt.title('M at symmetry')
plt.xlim(-0.04,0)
plt.savefig('plt_M.png', dpi=300)

fig4 = plt.figure("Figure 4")
plt.scatter(xpR_symm,gamma_symm, color = 'blue')
plt.grid()
plt.xlabel('X/R')
plt.ylabel('gamma')
plt.title('gamma at symmetry')
plt.xlim(-0.04,0)
plt.savefig('plt_gamma.png', dpi=300)

fig5 = plt.figure("Figure 5")
plt.scatter(xpR_symm,Kn_symm, color = 'blue')
plt.grid()
plt.xlabel('X/R')
plt.ylabel('Kn')
plt.title('Kn at symmetry')
plt.xlim(-0.04,0)
plt.savefig('plt_Kn.png', dpi=300)

fig6 = plt.figure("Figure 6")
plt.plot(xpR_symm,CO2_symm, color = 'blue')
plt.plot(xpR_symm,CO_symm, color = 'green')
plt.plot(xpR_symm,O2_symm, color = 'pink')
plt.plot(xpR_symm,O_symm, color = 'yellow')
plt.plot(xpR_symm,C2_symm, color = 'grey')
plt.plot(xpR_symm,C_symm, color = 'red')
plt.grid()
plt.xlabel('X/R')
plt.ylabel('mass fraction')
plt.title('Chemical species at symmetry')
plt.xlim(-0.04,0)
plt.legend(['CO2','CO','O2','O','C2','C',])
plt.savefig('plt_chemistry.png', dpi=300)

fig7 = plt.figure("Figure 7")
plt.scatter(xpR_symm,T_symm, color = 'blue')
plt.grid()
plt.xlabel('X/R')
plt.ylabel('T')
plt.title('T at symmetry (K)')
plt.xlim(-0.04,0)
plt.savefig('plt_T.png', dpi=300)

fig8 = plt.figure("Figure 8")
plt.scatter(xpR_symm,P_symm, color = 'blue')
plt.grid()
plt.xlabel('X/R')
plt.ylabel('P')
plt.title('P at symmetry (Pa)')
plt.xlim(-0.04,0)
plt.savefig('plt_P.png', dpi=300)

fig9 = plt.figure("Figure 9")
plt.scatter(xpR_symm,dVdy_symm, color = 'blue')
plt.grid()
plt.xlabel('X/R')
plt.ylabel('dVdy (1/s)')
plt.title('dVdy at symmetry')
plt.xlim(-0.04,0)
plt.savefig('plt_dVdy.png', dpi=300)

fig10 = plt.figure("Figure 10")
plt.scatter(x_cfdcomp_t_vibrationalelectric,T_cfdcomp_t_vibrationalelectric, color = 'red')
plt.scatter(x_cfdcomp_t_transrot,T_cfdcomp_t_transrot, color = 'orange')
plt.plot(xpR_symm,T_symm, color = 'blue')
plt.grid()
plt.xlabel('X/R')
plt.ylabel('T (K)')
plt.title('temperature at symm')
plt.legend(['Vibro-electronic (Moreira et al.)','Trans-rotational (Moreira et al.)', 'This model'])

plt.xlim(-0.04,0)
plt.savefig('plt_Tcomp.png', dpi=300)
