from __future__ import absolute_import, unicode_literals, print_function
import pymultinest
import math, os
from fm import *
import pdb
import numpy as np
import pickle
from matplotlib.pyplot import *
if not os.path.exists("chains"): os.mkdir("chains")
xsects=xsects(909,3333) #lower wavenumber, upper wavenumber...to convert to wl [um] take 1E4/wno (here, 11 - 3 um)


#prior--only use this to transform "unit (0-1)" coordinates to
#physical values...
#for instance to transform something to have a value between -6 and 6
#would require -6+12*cube[0]. When cube[0]=0, then the transformed 
#cube[0] will equal -6+12.*0, or -6. When cube[0]=1, then the transformed cube[0]
#will be -6+12*1, or 6.
def prior(cube, ndim, nparams): 
	#FINDME--prior ranges (uniform)
        #temperature profile prior ranges
	cube[0]=100.+700.*cube[0]  #Surface temperature range
	cube[1]=0.5+1.*cube[1]  #xRp
	cube[2]=-6+8.*cube[2]  #log cloud-top-pressure
        #background molecular weight
	cube[3]=2.0+42.*cube[3]  #backgroud gas mmw (H2=2.0, CO2=44)
        #gas mixing ratios, log mixing ratios go from -12 to 0 (1E-12 to 1)
	cube[4]=-12.+12*cube[4]  #H2O
	cube[5]=-12.12*+cube[5]  #CH4
	cube[6]=-12.+12*cube[6]  #CO2
	cube[7]=-12.+12*cube[7]  #O3
	cube[8]=-12.+12*cube[8]  #N2O
	cube[9]=-12.+12*cube[9]   #CO


def loglike(cube, ndim, nparams):
   Tsfc,xRp,logPc, Bkg_mmw=cube[0],cube[1],cube[2],cube[3]  #TP params
   logH2O, logCH4, logCO2, logO3, logN2O,logCO= cube[4], cube[5], cube[6],cube[7],cube[8],cube[9]  #gas params

   #planet params--if you changed this in make_spec.py, you must change it here to be consistent
   Rp= 0.910# Planet radius in Earth Radii
   Rstar=0.117   #Stellar Radius in Solar Radii
   M = 0.772  #    Mass in Earth Masses

   #if want to switch off certain parameters, fix their values here and remove them from "cube"-----
   #Tsfc=280.  #this is "surface temp" (isothermal below surface pressure at this temperature)
   logPsfc=2.0  #log surface pressure
   gam_trop=0.0#0.19  #troposphere adiabatic index, gamma (dlnT/dlnP=gamma)
   logPtrop=-0.6  #log tropopause pressure
   gam_strat=0.0#-0.05 #strastophsere adiabatic index -- 
   logPstrat=-3.0  #stratopause pressure--isothermal above this
   #log gas mixing ratios (loosley based off of Hu et al. 2012; Robinson et al. 2011)
   #logH2O=-1.5
   #logCH4=-6.3
   #logCO2=-3.45 #-3.4
   #logO3=-6.5
   #logN2O=-6.3
   #logCO=-7.0
   #Bkg_mmw=28.6  #uknown background gas mmw

    #checking various bad things and returning a -infinity log-liklihood (e.g., VMR's can't sum to be greater than 1)
   if (10**logH2O+10**logCH4+10**logCO2+10**logO3+10**logN2O+10**logCO) >= 1.0: 
      return -np.inf
    
   x=np.array([Tsfc,logPsfc,gam_trop,logPtrop,gam_strat,logPstrat, xRp*Rp, Rstar, M,logPc, Bkg_mmw, logH2O, logCH4, logCO2, logO3, logN2O,logCO])
   y=fx(x,xsects) 
   y_mod=y[0]
   loglikelihood=-0.5*np.sum((y_meas-y_mod)**2/err**2)
   return loglikelihood

########################
#reading in the data from text file
########################
#MCMC output file name--saved to pickle
outfile='MCMC.pic'

########################
#reading in the data from text file
########################
wlgrid, y_meas, err=pickle.load(open('data.pic','rb'))


n_params=10  #FINDME


pymultinest.run(loglike, prior, n_params, outputfiles_basename='./chains/template_',resume=False, verbose=True,n_live_points=1000)
a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename='./chains/template_')
s = a.get_stats()


#
output=a.get_equal_weighted_posterior()
pickle.dump(output,open(outfile,"wb"))


