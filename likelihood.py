import numpy as np
import pandas as pd
import scipy.stats as st
import reweight_nc as rw

# Q2 bin centers
Q2   = np.linspace(0.1,1.0,13)

# NCEp selection in data
''' 5e19 BNB data - Off-beam (LR score > 0.9) '''
data = np.array([ 154.2291, 103.834, 71.4301, 42.5502, 33.9585, 28.7227,
                   16.1572, 11.5655,  8.5917,  5.0524, -1.9476, -0.4869])
sig  = np.array([ 20.38584929, 14.3244555, 11.0829637, 8.54148744,
                   7.70048741, 6.40548338, 4.98446179,  3.8968395,
                   3.69587625, 2.81927055,     0.9738,     0.4869])

def getdata():

    # put measured data into one tuple for plotting and such
    D    = data,sig

    return Q2,D

def lnprob(theta,grand,data_mc,data_drt,priorform='zanalysis'):

    ''' log probability (posterior*prior) '''

    lp = lnprior(theta,grand,form=priorform)
    return lp + lnlike(theta,grand,data_mc,data_drt)

def lnlike(theta,grand,data_mc,data_drt):

    ''' ratio of cross-sections log likelihood '''

    GaS,lambdaA,SA = theta
    ncevts,ncsig   = rw.rwgt_ratio_mc(GaS,lambdaA,SA,Q2,data_mc,grand)
    drtevts,drtsig = rw.rwgt_ratio_mc(GaS,lambdaA,SA,Q2,data_drt,grand)

    model          = 0.1676*ncevts + 0.4326*0.5*drtevts
    tsig           = np.sqrt(sig**2 + (0.1676*ncsig)**2 + (0.4326*0.5*drtsig)**2)
    lnlikevec      = -0.5*np.square((data-model)/tsig) - np.log(tsig)

    return -len(data)/2.*np.log(2.*np.pi) + np.sum(lnlikevec)

def lnprior(theta,grand,form='zanalysis'):

    ''' set the prior probability for each parameter '''

    formopts = ['zanalysis','systematics','informative','uniform','reasonable','recreate']
    assert form in formopts, 'not a valid prior option'
  
    # prior log prob for each param. in theta
    # needs work

    a0s,a1s,a2s = theta

    if form == 'zanalysis':
        #if MaS >= 2./5. and np.abs(a2) <= 5.*np.abs(ds):
        if np.abs(a0s) > 10 or np.abs(a1s) > 10 or np.abs(a2s) > 10:
            return -np.inf

        a3s = -20.*a0s - 10.*a1s - 4.*a2s
        a4s = 45.*a0s + 20.*a1s + 6.*a2s
        if np.abs(a3s) > 10 or np.abs(a4s) > 10:
            return -np.inf

        a5s = -36.*a0s - 15.*a1s - 4.*a2s
        a6s = 10.*a0s + 4.*a1s + 1.*a2s
        if np.abs(a5s) > 10 or np.abs(a6s) > 10:
            return -np.inf

        lps = st.norm.logpdf(grand)

        return sum(lps)

    if form == 'systematics':
        lps = st.norm.logpdf(grand)

        return sum(lps)

    if form == 'pate13':
        # get prior from previous experiments
        lpds = st.norm.logpdf(ds,-0.30,0.42)
        lpma  = st.norm.logpdf(MaS,1.1,1.1)
        lpa2  = st.norm.logpdf(MaS,0.36,.50)

        return lpds+lpma+lpa2

    if form == 'uniform':
        # uniform prior --> likelihood
        if MaS > 0.:
            return 0.0
        return -np.inf

    if form == 'reasonable':
        # use sane, but loose priors on theta
        # if in these ranges p=1, otherwise p=0
        if -1. < ds < 1. and .5 < MaS < 2.:
            return 0.0
        return -np.inf

    if form == 'recreate':
        # use the exact parameters used in e734 paper --- not implemented yet
        # if in these ranges p=1, otherwise p=0
        if -1. < ds < 1. and .5 < MaS < 2.:
            return 0.0
        return -np.inf

