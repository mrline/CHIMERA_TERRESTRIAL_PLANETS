import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.pyplot import *
from fm import *
import pickle
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import time
xsects=xsects(909,3333) #lower wavenumber, upper wavenumber...to convert to wl [um] take 1E4/wno (here, 11 - 3 um)



#TP profile parameters--using a "4 layer" model--an isothermal region (below surface), a troposphere, a stratosphere, and an isothermal "thermosphere"
Tsfc=280.  #this is "surface temp" (isothermal below surface pressure at this temperature)
logPsfc=0.0  #log surface pressure
gam_trop=0.#0.19  #troposphere adiabatic index, gamma (dlnT/dlnP=gamma)
logPtrop=-0.6  #log tropopause pressure
gam_strat=-0.0#-0.05 #strastophsere adiabatic index -- 
logPstrat=-3.0  #stratopause pressure--isothermal above this
#planet params
Rp= 0.910# Planet radius in Earth Radii
Rstar=0.117   #Stellar Radius in Solar Radii
M = 0.772  #    Mass in Earth Masses
#cloud params
logPc=-0.25 #log cloud-top-pressure bar (here set to refractive boundary)
#log gas mixing ratios (loosley based off of Hu et al. 2012; Robinson et al. 2011)
logH2O=-5.5
logCH4=-6.3
logCO2=-3.45 #-3.4
logO3=-6.5
logN2O=-6.3
logCO=-7.0
Bkg_mmw=28.6  #unknown background gas mmw

'''
#refractive boundary
Rp= 0.910
Rstar=0.117   
v0=2.93E-4
T0=200.
a=0.030
mu=28.6
g=10.
pmax=23E-3*(1.23E-4/v0)*(T0/130.)**1.5*(Rstar)*(5.2/a)*(10.973/(Rp))**0.5*(2.2/mu)**0.5*(24.8/g)**0.5
'''

#state vector
            # 0      1       2         3	4	5       6    7    8   9   10	    11	   12	    13	     14    15     16
          #Tsfc,logPsfc,gam_trop,logPtrop,gam_strat,logPstrat, Rp, Rstar, M,logPc,Bkg_mmw, logH2O, logCH4, logCO2, logO3, logN2O,logCO
x=np.array([Tsfc,logPsfc,gam_trop,logPtrop,gam_strat,logPstrat, Rp, Rstar, M,logPc, Bkg_mmw, logH2O, logCH4, logCO2, logO3, logN2O,logCO])


y_mod,wno,atm=fx(x,xsects)


#read in external noise file if available
wlgrid, junk,junk, err0=np.loadtxt('noise_R100.txt').T  #must have same wlgrid as CK coeffs
err=np.interp(1E4/wno[::-1],wlgrid,err0)
ntran = 25.0
noise_floor = 5.0E-6
err=np.sqrt( (err[::-1]*(1.0/np.sqrt(ntran)))**2.0 + (noise_floor)**2.0 )
fname='DryEarth_1-30um_R100_1tran'

#defining data array
y_meas=np.zeros(len(y_mod))

#adding gaussian noise (note I turned this off now--see justification in Feng et al. 2018)
for i in range(y_meas.shape[0]): y_meas[i]=y_mod[i]#+np.random.randn(1)*err[i]

#computing chi-square of random noise instance
print(np.sum((y_meas-y_mod)**2/err**2)/len(y_meas))

#dumping pickles--model then noised up data data
output=[1E4/wno, y_mod]
pickle.dump(output,open("Model.pic","wb"))  #spectral model to be noised up by instrument noise model

output=[1E4/wno, y_meas,err]  #noised up "synthetic" spectrum
pickle.dump(output,open("data.pic","wb"))

#plotting stuff
wlgrid=1E4/wno
ymin=1E6*np.min(y_mod)*0.98
ymax=1E6*np.max(y_mod)*1.02
fig1, ax=subplots()
xlabel('$\lambda$ ($\mu$m)',fontsize=18)
ylabel('(R$_{p}$/R$_{\star}$)$^{2}$ [ppm]',fontsize=18)
minorticks_on()
plot(1E4/wno, y_mod*1E6)
errorbar(wlgrid, y_meas*1E6, yerr=err*1E6, xerr=None, fmt='ok',alpha=0.25)
plot(wlgrid, y_mod,'ob')
ax.set_xscale('log')
ax.set_xticks([1,1.4,2,2.7,3,3.3,4.3,5,7.7,11,15,17,20,30])
ax.axis([1,30,ymin,ymax])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.tick_params(length=10,width=1,labelsize='large',which='major')



savefig(fname+'_spectrum.pdf',fmt='pdf')
show()
close()
pdb.set_trace()










