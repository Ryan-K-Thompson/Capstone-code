# -*- coding: utf-8 -*-
"""
Shock and Detonation Toolbox Demo Program

Generate plots and output files for a steady reaction zone between a shock and a blunt body using
Hornung model of linear profile of rho u on stagnation streamline.
 
################################################################################
Theory, numerical methods and applications are described in the following report:

SDToolbox - Numerical Tools for Shock and Detonation Wave Modeling,
Explosion Dynamics Laboratory, Contributors: S. Browne, J. Ziegler,
N. Bitter, B. Schmidt, J. Lawson and J. E. Shepherd, GALCIT
Technical Report FM2018.001 Revised January 2021.

Please cite this report and the website if you use these routines. 

Please refer to LICENCE.txt or the above report for copyright and disclaimers.

http://shepherd.caltech.edu/EDL/PublicResources/sdt/

################################################################################ 
Updated September 2018
Tested with: 
    Python 3.5 and 3.6, Cantera 2.3 and 2.4
Under these operating systems:
    Windows 8.1, Windows 10, Linux (Debian 9)
"""
import sys
sys.path.append(r'C:\Users\patag\Desktop\Uni\2022\Capstone\Cantera\SDToolbox\Python3')
import cantera as ct
from sdtoolbox.postshock import PostShock_fr
from sdtoolbox.stagnation import stgsolve
from sdtoolbox.utilities import znd_plot
import matplotlib.pyplot as plt
from numpy import logspace,log10
import pickle

graphing = True

# =============================================================================
#  PWK
#  p1 = 796.2 # Pa
#  T1 = 1165.8 # K
#  gas.TPX = T1, p1, "CO2:0.02629, CO:0.61999, C2:0, O2:0.05779, C:0, O:0.29593"   #equilibriate the gas at the initial conditions, from this point forward, the gas will be assumed to be frozen due to the low dakholmer number
#  u1 = 2850 # m/s

# Mars
# p1 = 72.3 # Pa
# T1 = 182.1 # K
# gas.TPX = T1, p1, 'CO2:1'   #equilibriate the gas at the initial conditions, from this point forward, the gas will be assumed to be frozen due to the low dakholmer number
# u1 = 2273.5 # m/s
# =============================================================================


P1 = 72.3; T1 = 182.1; U1 = 2273.5; Delta = 0.1
PWK = "CO2:0.02629, CO:0.61999, C2:0, O2:0.05779, C:0, O:0.29593, N:0, NO:0"
CO2 = "CO2:0.95, N2:0.05"
q = CO2
mech = 'airNASA9noions.cti'

gas1 = ct.ThermoPhase("airNASA9noions.cti")
gas1.TPX = T1,P1,q
nsp = gas1.n_species

# FIND POST SHOCK STATE FOR GIVEN SPEED
gas = PostShock_fr(U1, P1, T1, q, mech)
print("post shock gas")
print(gas)
t_end = 0.1

# MATLAB's ode solvers automatically interpolate additional values beyond what is
# needed to meet tolerances. Use the optional t_eval keyword from the integrator
# to do something similar (otherwise plotted solution very sparse at short times)
t_dense = logspace(-10,log10(t_end),num=500)

# SOLVE REACTION ZONE STRUCTURE ODES
out = stgsolve(gas,gas1,U1,Delta,t_end=t_end,max_step=1e-4,t_eval=t_dense)

if graphing:
    # Demonstrate using the figure outputs from znd_plot to make further formatting
    # adjustments, then saving
    # Note that znd_plot still works with stgsolve since the output structures
    # contain many of the same fields
    
    maxx = max(out['distance'])
    figT,figP,figM,figTherm,figSpec = znd_plot(out,maxx=maxx,major_species='All',
                                               xscale='linear',show=False)
    
    figSpec.axes[0].set_xlim((1e-5,None))
    figSpec.axes[0].set_ylim((1e-10,1))
    
    plt.show()
    
    figT.savefig('stg-T.eps',bbox_inches='tight')
    figP.savefig('stg-P.eps',bbox_inches='tight')
    figM.savefig('stg-M.eps',bbox_inches='tight')
    figTherm.savefig('stg-therm.eps',bbox_inches='tight')
    figSpec.savefig('stg-Y.eps',bbox_inches='tight')

# Can't pickle Cantera Solution objects due to underlying C++ structure, 
# so remove it from dictionary first
out.pop('gas1') 
pickle.dump(out,open('stg.p','wb'))



#% save output file for later comparison
#stg = out;
#save 'stg.mat' stg ;

