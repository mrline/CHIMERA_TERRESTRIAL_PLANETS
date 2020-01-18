import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.pyplot import *
import pickle
import pylab
import pdb
import corner
from scipy import interp
rc('font',family='serif')

#TP profile function for plotting
def TP_simple2(Tsfc, Psfc, gam_trop, Ptrop, gam_strat,Pstrat, P):
        T=np.zeros(len(P))

        #troposphere T--adibat
        Ttrop=Tsfc*(P/Psfc)**gam_trop  #P0=sfc p, Trop T
        
	#stratosphere
	Tpause=Ttrop[P <= Ptrop ][-1]  #tropopause Temp
	PPtrop=P[P <= Ptrop ][-1]
	Tstrat=Tpause*(P/PPtrop)**gam_strat

	#merging troposphere and stratosphfere
	T[P > Ptrop]=Ttrop[P > Ptrop]
	T[P <= Ptrop]=Tstrat[P <= Ptrop]

	#isothermal below surface (making surface a blackbody)
	T[P >= Psfc]=T[P >= Psfc][0]

	#isothermal above "stratopause" pressure, Pstrat
	T[P<=Pstrat]=T[P<=Pstrat][-1]
	T[T<=10]=10
	T[T>=1000]=1000
     
	#pdb.set_trace()
        return T



#load data
data=pickle.load(open("data.pic",'rb'))
wlgrid=data[0]
y_meas=data[1]
err=data[2]


#load chains*************************************
fname1='MCMC.pic'
pic=pickle.load(open(fname1,'rb'))
chain=pic[:,:-1]
lnprob=pic[:,-1]

fname='DryEarth_1-30um_R100_1tran'

chi2=-2.*lnprob/wlgrid.shape[0]

Npars=chain.shape[1]

samples=chain
chi2flat=chi2
chi2best=np.min(chi2flat)
locbest=np.array((chi2flat==chi2best).nonzero())[0,0]
xbest=samples[locbest,:]


#plotting stair-pairs plot via traingle.py******************
truth=np.array([280., 1.0, -0.25, 28.6, -1.5, -6.3, -3.45, -6.5, -6.3, -7.0])

titles=np.array(['T$_{iso}$','$\\times$R$_{p}$','log(CTP)','Bkg MMW', 'log(H$_2$O)','log(CH$_4$)', 'log(CO$_2$)','log(O$_3$)', 'log(N$_2$O)','log(CO)'])
priorlow=np.array([100, 0.5 ,  -6 , 2.0, -12, -12, -12, -12, -12, -12])
priorhigh=np.array([800, 1.5,   1.5, 44.,  0,    0,   0,  0 ,   0,   0 ])
Npars=len(priorlow)
ext=np.zeros([2,Npars])
for i in range(Npars):
	med=np.percentile(samples[:,i],50)
	stdev=np.std(samples[:,i])
	ext[0,i]=med-4.*stdev
	if ext[0,i] < priorlow[i]:
		ext[0,i]=priorlow[i]
	ext[1,i]=med+4.*stdev
	if ext[1,i] > priorhigh[i]:
		ext[1,i]=priorhigh[i]
ext=ext.T
ext[:,0]=priorlow
ext[:,1]=priorhigh
#'''
corner.corner(samples,labels=titles, bins=25,plot_datapoints='False',quantiles=[.16,0.5,.84],show_titles='True',plot_contours='True',extents=ext,levels=(1.-np.exp(-(1)**2/2.),1.-np.exp(-(2)**2/2.),1.-np.exp(-(3)**2/2.)),truths=truth)

savefig(fname+"_stair_pairs.pdf",format='pdf')
show()
close()

'''
#plotting a subset
# 0      1       2         3	    4	     5		6     7        8       9      10      11   12
#Tsfc,logPsfc,gam_trop,logPtrop,gam_strat,logPstrat,Bkg_mmw, logH2O, logCH4, logCO2, logO3, logN2O,logCO
whichpars=np.array([6,7,8,9,10,11,12])  #indicies of paremeters we want to plot
corner.corner(samples[:,whichpars],labels=titles[whichpars], bins=25,plot_datapoints='False',quantiles=[.16,0.5,.84],show_titles='True',plot_contours='True',extents=ext[whichpars,:],levels=(1.-np.exp(-(1)**2/2.),1.-np.exp(-(2)**2/2.),1.-np.exp(-(3)**2/2.)),truths=truth[whichpars])
savefig(fname+"_abundances_stair_pairs.pdf",format='pdf')
show()
close()

#plotting a subset
# 0      1       2         3	    4	     5		6     7        8       9      10      11   12
#Tsfc,logPsfc,gam_trop,logPtrop,gam_strat,logPstrat,Bkg_mmw, logH2O, logCH4, logCO2, logO3, logN2O,logCO
whichpars=np.array([0,1,2,3,4,5])  #indicies of paremeters we want to plot
corner.corner(samples[:,whichpars],labels=titles[whichpars], bins=25,plot_datapoints='False',quantiles=[.16,0.5,.84],show_titles='True',plot_contours='True',extents=ext[whichpars,:],levels=(1.-np.exp(-(1)**2/2.),1.-np.exp(-(2)**2/2.),1.-np.exp(-(3)**2/2.)),truths=truth[whichpars])
savefig(fname+"_TP_stair_pairs.pdf",format='pdf')
show()
close()
'''
#'''

#'''


pdb.set_trace()
#Plotting spectra spread**********************
y_arr=pickle.load(open('y_arr.pic','rb'))
y_mod_arr=y_arr[0]*1000
y_hires_arr=y_arr[1]*1000
y_mod_best=y_arr[2]*1000
y_hires_best=y_arr[3]*1000
wlgrid=y_arr[4]
wnocrop=y_arr[5]
NN=y_mod_arr.shape[0]
y_hires_best_conv=tophatfold(1E4/wnocrop[::-1], y_hires_best[::-1],0.035)

y_conv_arr=[]
for i in range(NN):
    y_hires=y_hires_arr[i,:]
    y_conv=tophatfold(1E4/wnocrop[::-1], y_hires[::-1],0.035)
    y_conv_arr=np.concatenate([y_conv_arr,y_conv])

y_conv_arr=y_conv_arr.reshape(NN,wnocrop.shape[0])

y_median=np.zeros(wnocrop.shape[0])
y_high_1sig=np.zeros(wnocrop.shape[0])
y_high_2sig=np.zeros(wnocrop.shape[0])
y_low_1sig=np.zeros(wnocrop.shape[0])
y_low_2sig=np.zeros(wnocrop.shape[0])

for i in range(wnocrop.shape[0]):
    percentiles=np.percentile(y_conv_arr[:,i],[4.55, 15.9, 50, 84.1, 95.45])
    y_low_2sig[i]=percentiles[0]
    y_low_1sig[i]=percentiles[1]
    y_median[i]=percentiles[2]
    y_high_1sig[i]=percentiles[3]
    y_high_2sig[i]=percentiles[4]



fig, ax=subplots()
def ticklab_minor(x,pos):
	if x % 2 == 0:
		return '%2.0f' %(x)
	else:
		return ""

def ticklab_major(x,pos):
	if x % 1 == 0:
		return '%2.0f' %(x)
	else:
		return ""


ymax=1.2*np.max(y_meas)*1E3
axis([1,10,0,ymax])   
minorticks_on()
semilogx()
ax.xaxis.set_minor_formatter(FuncFormatter(ticklab_minor))
ax.xaxis.set_major_formatter(FuncFormatter(ticklab_major))
ax.tick_params(axis='x',length=10,width=1,which='minor')
ax.tick_params(axis='y',length=10,width=1,which='major')
xlabel('$\lambda$ ($\mu$m)',fontsize=18)
ylabel('F$_{p}$/F$_{\star}$ [10$^{-3}$]',fontsize=18)
minorticks_on()
fill_between(1E4/wnocrop,y_low_2sig[::-1],y_high_2sig[::-1],facecolor='r',alpha=0.1,edgecolor='None')  
fill_between(1E4/wnocrop,y_low_1sig[::-1],y_high_1sig[::-1],facecolor='r',alpha=1.,edgecolor='None')  
plot(1E4/wnocrop,y_median[::-1],color='b')
plot(1E4/wnocrop,y_hires_best_conv[::-1],color='g')
errorbar(wlgrid,y_meas*1000,yerr=err*1000,xerr=None, fmt='D',color='black')
plot(wlgrid, y_mod_best,'og')
#errorbar(wlgrid_bin_center,y_meas_bin*1000,yerr=err_bin*1000,xerr=None, fmt='s',color='black')
savefig(fname+'_spectra.pdf',format='pdf')


close()



pdb.set_trace()
