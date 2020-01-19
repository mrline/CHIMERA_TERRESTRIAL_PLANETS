import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.pyplot import *
import pickle
import pylab
import pdb
import corner
from scipy import interp
rc('font',family='serif')



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

fname='DryEarth_2-11um_R100_25tran'

chi2=-2.*lnprob/wlgrid.shape[0]

Npars=chain.shape[1]

samples=chain
chi2flat=chi2
chi2best=np.min(chi2flat)
locbest=np.array((chi2flat==chi2best).nonzero())[0,0]
xbest=samples[locbest,:]


#plotting stair-pairs plot via traingle.py******************
truth=np.array([280., 1.0, -0.25, 28.6, -5.5, -6.3, -3.45, -6.5, -6.3, -7.0])

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



pdb.set_trace()
