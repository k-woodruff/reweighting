import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import likelihood as lk
import reweight_genie as rw
from tqdm import tqdm

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import cmocean
import csv

'''
Plots of surfaces and distributions
'''

def gridPlot2D(numpts,data_mc,data_dirt,a0s):

    # make empty grid
    a1svec  = np.linspace(-8,5,numpts)
    a2svec   = np.linspace(-50,50,numpts)

    a1s,a2s  = np.meshgrid(a1svec,a2svec)
    lnlikes = np.zeros((len(a1s),len(a2s)))

    # make lnlikelihood matrix
    for i in tqdm(range(len(a1s))):
        for j in range(len(a2s)):
            
            theta = a0s,a1s[i][j],a2s[i][j]
            lnlikes[i][j] = lk.lnlike(theta,np.zeros(6),data_mc,data_dirt)

    '''
    # convert to delta s and MA s
    a3s = -20.*a0s - 10.*a1s - 4.*a2s
    a4s = 45*a0s + 20.*a1s + 6.*a2s
    a5s = -36.*a0s - 15.*a1s - 4.*a2s
    a6s = 10.*a0s + 4.*a1s + a2s

    t_cut = 9.*0.13957**2
    t0  = -0.28
    z0  = (np.sqrt(t_cut) - np.sqrt(t_cut-t0))/(np.sqrt(t_cut) + np.sqrt(t_cut - t0))

    GaS = a0s + a1s*z0 + a2s*z0**2 + a3s*z0**3 + a4s*z0**4 + a5s*z0**5 + a6s*z0**6
    MA  = np.sqrt(np.abs(2.*GaS/(a1s + 2.*a2s*z0 + 3.*a3s*z0**2 + 4.*a4s*z0**3 + 5.*a5s*z0**4 + 6.*a6s*z0**5)))
    '''

    plt.figure(figsize=(5,3.5))
    plt.pcolor(a2s,a1s,np.exp(lnlikes),cmap=cmocean.cm.tempo)
    plt.xlabel(r'$a_2^s$',fontsize=16)
    plt.ylabel(r'$a_1^s$',fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    #plt.xlim(MA.min(),MA.max());
    #plt.ylim(GaS.min(),GaS.max());
    plt.colorbar()

def makeSamplePlot(samples,data_mc,data_sysval):
    # get data
    Q2,D     = lk.getdata()
    data,sig = D
    enu = 1.

    y = np.zeros((len(Q2),len(samples[:,0])))
    for i in range(len(samples[:,0])):
        theta  = samples[i,:]
        nc,cc  = rw.rwgt_ratio_mc(theta[0],theta[1],0.,Q2,data_mc,data_sysval,2)
        y[:,i] = nc/cc

    plt.plot(Q2,y,color='steelblue',alpha=0.05)
    plt.errorbar(Q2,data,yerr=sig,fmt='o',color='tomato')
    #plt.yscale('log')

def makePlot(data_mc,data_sysval):

    MA  = 1.06 # e734
    GaS = -.15

    # Data points
    Q2,D  = lk.getdata()
    data,sig = D

    nc,cc  = rw.rwgt_ratio_mc(GaS,MA,0.,Q2,data_mc,data_sysval,2)
    yTruth = nc/cc

    fig, ax   = plt.subplots(nrows=1,ncols=2,sharey=True,figsize=(10,5))
    gs = gridspec.GridSpec(10,4)
    gs.update(wspace=0.025, hspace=0.05)

    # changing DeltaS in top left
    MA  = 1.06
    GaSvec = np.linspace(-.5,.3,50)
    y = np.zeros((len(Q2),len(GaSvec)))
    ncmin,ccmin = rw.rwgt_ratio_mc(GaSvec[0],MA,0.,Q2,data_mc,data_sysval,2)
    ymin = ncmin/ccmin
    for i,GaS in enumerate(GaSvec):
        nc,cc  = rw.rwgt_ratio_mc(GaS,MA,0.,Q2,data_mc,data_sysval,2)
        y[:,i] = nc/cc
    
    ax[0].plot(Q2,y,color='steelblue',alpha=0.7)
    ax[0].plot(Q2,ymin,color='steelblue',linewidth=3)
    ax[0].plot(Q2,yTruth,color='black')
    ax[0].errorbar(Q2,data,yerr=sig,fmt='o',color='tomato')
    ax[0].set_title(r'$\Delta S$ from -0.5 to 0.3')
    ax[0].set_ylabel(r'NCE/CCQE ratio')

    # changing MA in top right
    MAvec  = np.linspace(0.75,1.5,50)
    GaS = -.15
    y = np.zeros((len(Q2),len(MAvec)))
    ncmin,ccmin = rw.rwgt_ratio_mc(GaS,MAvec[0],0.,Q2,data_mc,data_sysval,2)
    ymin = ncmin/ccmin
    for i,MA in enumerate(MAvec):
        nc,cc  = rw.rwgt_ratio_mc(GaS,MA,0.,Q2,data_mc,data_sysval,2)
        y[:,i] = nc/cc

    ax[1].plot(Q2,y,color='green',alpha=0.7)
    ax[1].plot(Q2,ymin,color='green',linewidth=3)
    ax[1].plot(Q2,yTruth,color='black')
    ax[1].errorbar(Q2,data,yerr=sig,fmt='o',color='tomato')
    ax[1].set_title(r'$M_A^s$ from 0.75 to 1.5')
    ax[1].set_xlabel(r'$Q^2$ [GeV$^2$]')

def makeTriangle(ds,ls,sa,wtruth=0):
   
    fig,((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,sharex='col')
    '''
    ax2r = ax2.twinx()
    ax3r = ax3.twinx()
    ax6r = ax6.twinx()
    '''
    ax2.axis('off')
    ax3.axis('off')
    ax6.axis('off')

    dsmin,dsmax = -0.25,0.01
    lsmin,lsmax = 0.51,1.49
    samin,samax = -.28,0.35

    # 1d hists
    ax1.hist(sa,bins=np.linspace(samin,samax,40),histtype='step',color='midnightblue')
    ax5.hist(ls,bins=np.linspace(lsmin,lsmax,40),histtype='step',color='midnightblue')
    ax9.hist(ds,bins=np.linspace(dsmin,dsmax,40),histtype='step',color='midnightblue')

    # 2d scatter plots
    ax4.plot(sa,ls,'o',markeredgecolor='none',markersize=2,color='midnightblue',alpha=0.5)
    ax7.plot(sa,ds,'o',markeredgecolor='none',markersize=2,color='midnightblue',alpha=0.5)
    ax8.plot(ls,ds,'o',markeredgecolor='none',markersize=2,color='midnightblue',alpha=0.5)

    # 2d hists
    '''
    ax4.hist2d(sa,ls,bins=(np.linspace(samin,samax,15),np.linspace(lsmin,lsmax,15)),cmap=cmocean.cm.tempo)
    ax7.hist2d(sa,ds,bins=(np.linspace(samin,samax,15),np.linspace(dsmin,dsmax,15)),cmap=cmocean.cm.tempo)
    ax8.hist2d(ls,ds,bins=(np.linspace(lsmin,lsmax,15),np.linspace(dsmin,dsmax,15)),cmap=cmocean.cm.tempo)
    '''

    # 2d contours
    '''
    counts2,ls2bins,sa2bins=np.histogram2d(ls,sa,bins=5)
    counts3,ds3bins,sa3bins=np.histogram2d(ds,sa,bins=5)
    counts6,ds6bins,ls6bins=np.histogram2d(ds,ls,bins=5)

    from scipy import interpolate
    n = 1000
    LS2,SA2 = np.meshgrid(ls2bins,sa2bins)
    LS3,SA3 = np.meshgrid(ds3bins,sa3bins)
    LS6,SA6 = np.meshgrid(ds6bins,ls6bins)
    Z2 = counts2/counts2.sum()
    Z3 = counts3/counts3.sum()
    Z6 = counts6/counts6.sum()
    t2 = np.linspace(0,Z2.max(),n)
    t3 = np.linspace(0,Z3.max(),n)
    t6 = np.linspace(0,Z6.max(),n)
    integral2 = ((Z2 >= t2[:,None,None])*Z2).sum(axis=(1,2))
    integral3 = ((Z3 >= t3[:,None,None])*Z3).sum(axis=(1,2))
    integral6 = ((Z6 >= t6[:,None,None])*Z6).sum(axis=(1,2))
    f2 = interpolate.interp1d(integral2,t2)
    f3 = interpolate.interp1d(integral3,t3)
    f6 = interpolate.interp1d(integral6,t6)
    ax2r.contour(Z2.transpose(),levels=f2(np.array([.999,0.9,0.68,0.50])),extent=[ls2bins.min(),ls2bins.max(),sa2bins.min(),sa2bins.max()],cmap=cmocean.cm.deep)
    ax3r.contour(Z3.transpose(),levels=f3(np.array([.999,0.9,0.68,0.50])),extent=[ds3bins.min(),ds3bins.max(),sa3bins.min(),sa3bins.max()],cmap=cmocean.cm.deep)
    ax6r.contour(Z6.transpose(),levels=f6(np.array([.999,0.9,0.68,0.50])),extent=[ds6bins.min(),ds6bins.max(),ls6bins.min(),ls6bins.max()],cmap=cmocean.cm.deep)

    ax2r.get_yaxis().set_visible(False)
    ax3r.get_yaxis().set_visible(True)
    ax6r.get_yaxis().set_visible(True)
    ax2r.set_ylim(-3.0,0.5)
    ax3r.set_ylim(-3.0,0.5)
    ax6r.set_ylim(0.25,2)
    '''

    ax1.get_yaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)
    ax5.get_yaxis().set_visible(False)
    ax6.get_yaxis().set_visible(False)
    ax8.get_yaxis().set_visible(False)
    ax9.get_yaxis().set_visible(False)

    if(wtruth > 0):
        ax1.axvline(0.,color='firebrick',linestyle='--')
        ax4.axvline(0.,color='firebrick',linestyle='--')
        ax7.axvline(0.,color='firebrick',linestyle='--')
        #ax2.axvline(0.990,color='firebrick',linestyle='--')
        ax5.axvline(0.990,color='firebrick',linestyle='--')
        ax8.axvline(0.990,color='firebrick',linestyle='--')
        #ax3.axvline(-0.15204,color='firebrick',linestyle='--')
        #ax6.axvline(-0.15204,color='firebrick',linestyle='--')
        ax9.axvline(-0.15204,color='firebrick',linestyle='--')

        #ax2r.axhline(0.,color='firebrick',linestyle='--')
        #ax3r.axhline(0.,color='firebrick',linestyle='--')
        ax4.axhline(0.990,color='firebrick',linestyle='--')
        #ax6r.axhline(0.990,color='firebrick',linestyle='--')
        ax7.axhline(-0.15204,color='firebrick',linestyle='--')
        ax8.axhline(-0.15204,color='firebrick',linestyle='--')

    fig.subplots_adjust(hspace=0.05,wspace=0.05)

    ax1.yaxis.set_major_locator(MultipleLocator(50.0))
    ax4.yaxis.set_major_locator(MultipleLocator(0.2))
    ax7.yaxis.set_major_locator(MultipleLocator(0.1))

    ax1.xaxis.set_major_locator(MultipleLocator(0.2))
    ax4.xaxis.set_major_locator(MultipleLocator(0.2))
    ax7.xaxis.set_major_locator(MultipleLocator(0.2))

    ax2.xaxis.set_major_locator(MultipleLocator(0.2))
    ax5.xaxis.set_major_locator(MultipleLocator(0.2))
    ax8.xaxis.set_major_locator(MultipleLocator(0.2))

    ax3.xaxis.set_major_locator(MultipleLocator(0.1))
    ax6.xaxis.set_major_locator(MultipleLocator(0.1))
    ax9.xaxis.set_major_locator(MultipleLocator(0.1))
    '''
    ax3r.yaxis.set_major_locator(MultipleLocator(0.75))
    ax6r.yaxis.set_major_locator(MultipleLocator(0.25))
    ax3r.set_ylabel(r'$S_A$')
    ax6r.set_ylabel(r'$M_A^s$')
    '''
    ax1.set_xlim(samin,samax)
    ax4.set_xlim(samin,samax)
    ax7.set_xlim(samin,samax)

    ax5.set_xlim(lsmin,lsmax)
    ax8.set_xlim(lsmin,lsmax)

    ax9.set_xlim(dsmin,dsmax)

    ax4.set_ylim(lsmin,lsmax)
    ax7.set_ylim(dsmin,dsmax)
    ax8.set_ylim(dsmin,dsmax)

    # show +/- 1 sigma
    '''
    ax1.axvline(np.percentile(sa,15.9),color='gray')
    ax1.axvline(np.percentile(sa,84.1),color='gray')
    ax5.axvline(np.percentile(ls,15.9),color='gray')
    ax5.axvline(np.percentile(ls,84.1),color='gray')
    ax9.axvline(np.percentile(ds,15.9),color='gray')
    ax9.axvline(np.percentile(ds,84.1),color='gray')
    '''
    # show +/- 95% credible
    ax1.axvline(np.percentile(sa,5),color='gray')
    ax1.axvline(np.percentile(sa,95),color='gray')
    ax5.axvline(np.percentile(ls,5),color='gray')
    ax5.axvline(np.percentile(ls,95),color='gray')
    ax9.axvline(np.percentile(ds,5),color='gray')
    ax9.axvline(np.percentile(ds,95),color='gray')
    # show median
    ax1.axvline(np.percentile(sa,50),color='black')
    ax5.axvline(np.percentile(ls,50),color='black')
    ax9.axvline(np.percentile(ds,50),color='black')

    ax4.set_ylabel(r'$M_A^s$',fontsize=18)
    ax7.set_ylabel(r'$\Delta s$',fontsize=18)
    ax7.set_xlabel(r'$a_2$',fontsize=18)
    ax8.set_xlabel(r'$M_A^s$',fontsize=18)
    ax9.set_xlabel(r'$\Delta s$',fontsize=18)

def make2DTriangle(ds,ls,wtruth=0):
   
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col')

    ax2.axis('off')

    dsmin,dsmax = -0.5,0.1
    lsmin,lsmax = 0.1,4.9

    # 1d hists
    ax1.hist(ls,bins=np.linspace(lsmin,lsmax,30),histtype='step',color='midnightblue')
    ax4.hist(ds,bins=np.linspace(dsmin,dsmax,30),histtype='step',color='midnightblue')

    # 2d scatter plots
    ax3.plot(ls,ds,'o',markeredgecolor='none',markersize=2,color='midnightblue',alpha=0.5)

    # 2d hists
    '''
    ax3.hist2d(ls,ds,bins=(np.linspace(lsmin,lsmax,15),np.linspace(dsmin,dsmax,15)),cmap=cmocean.cm.tempo)
    '''

    # 2d contours
    '''
    counts3,ds3bins,ls3bins=np.histogram2d(ds,ls,bins=5)

    from scipy import interpolate
    n = 1000
    LS3,SA3 = np.meshgrid(ds3bins,ls3bins)
    Z3 = counts3/counts3.sum()
    t3 = np.linspace(0,Z3.max(),n)
    integral3 = ((Z3 >= t3[:,None,None])*Z3).sum(axis=(1,2))
    f3 = interpolate.interp1d(integral3,t3)
    ax3r.contour(Z3.transpose(),levels=f3(np.array([.999,0.9,0.68,0.50])),extent=[ls3bins.min(),ls3bins.max(),ds3bins.min(),ds3bins.max()],cmap=cmocean.cm.deep)

    ax3r.get_yaxis().set_visible(True)
    ax3r.set_ylim(0.1,2.)
    '''

    ax1.get_yaxis().set_visible(False)
    ax4.get_yaxis().set_visible(False)

    # show +/- 95% credible
    ax1.axvline(np.percentile(ls,5),color='gray')
    ax1.axvline(np.percentile(ls,95),color='gray')
    ax4.axvline(np.percentile(ds,5),color='gray')
    ax4.axvline(np.percentile(ds,95),color='gray')
    # show median
    ax1.axvline(np.percentile(ls,50),color='black')
    ax4.axvline(np.percentile(ds,50),color='black')

    if(wtruth > 0):
        ax1.axvline(0.990,color='firebrick',linestyle='--')
        ax3.axvline(0.990,color='firebrick',linestyle='--')

        ax3.axhline(-0.15204,color='firebrick',linestyle='--')
        ax4.axvline(-0.15204,color='firebrick',linestyle='--')

    fig.subplots_adjust(hspace=0.05,wspace=0.05)

    ax3.yaxis.set_major_locator(MultipleLocator(0.1))

    ax1.xaxis.set_major_locator(MultipleLocator(0.5))
    ax3.xaxis.set_major_locator(MultipleLocator(0.5))

    ax4.xaxis.set_major_locator(MultipleLocator(0.1))

    ax1.set_xlim(lsmin,lsmax)
    ax3.set_xlim(lsmin,lsmax)

    ax3.set_ylim(dsmin,dsmax)

    ax3.set_ylabel(r'$\Delta s$',fontsize=18)
    ax3.set_xlabel(r'$M_A^s$',fontsize=18)
    ax4.set_xlabel(r'$\Delta s$',fontsize=18)

def csv2array(csvname):
    
    with open(csvname,'rb') as f:
        lines = f.readlines()

    ds = np.array([])
    ls = np.array([])
    sa = np.array([])
    lk = np.array([])

    for i,line in enumerate(lines):
        if i%4 == 0:
            ds = np.append(ds,np.array((line.split('\r')[0]).split(',')).astype(float))
        elif i%4 == 1:
            ls = np.append(ls,np.array((line.split('\r')[0]).split(',')).astype(float))
        elif i%4 == 2:
            sa = np.append(sa,np.array((line.split('\r')[0]).split(',')).astype(float))
        elif i%4 == 3:
            lk = np.append(lk,np.array((line.split('\r')[0]).split(',')).astype(float))

    return ds,ls,sa,lk






