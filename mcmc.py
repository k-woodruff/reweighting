import numpy as np
import pandas as pd
import likelihood as lk
from tqdm import tqdm

import time
import csv

def run_mcmc(fname_mc,fname_drt,fname_out,startpos,startgrand,stepsize,num_it,chunksize):

    Deltas_chain = np.array([])
    Lambdas_chain = np.array([])
    SA_chain = np.array([])
    syst_chain = np.array([])
    lk_chain = np.array([])
    naccept = 0

    data_mc = pd.read_csv('data/'+fname_mc+'.csv',sep=',')
    data_drt = pd.read_csv('data/'+fname_drt+'.csv',sep=',')

    #startgrand = np.random.randn(6)
    #startgrand = np.zeros(6)

    num_chunks = int(np.floor(num_it/chunksize))
    lst_chunk = num_it % chunksize
    print 'Running mcmc'
    print 'Number of chunks: {}'.format(num_chunks + 1)
    #for ichnk in tqdm(range(num_chunks)):
    for ichnk in range(num_chunks):
        start = time.time()
        a0,a1,a2,sy,lk,na = metropolis_systematics(data_mc,data_drt,startpos,startgrand,stepsize,chunksize)
        startpos = a0[-1],a1[-1],a2[-1]
        startgrand = sy[-1]

        end = time.time()
        print 'This chunk: {} s'.format(end - start)

        Deltas_chain = np.append(Deltas_chain,a0)
        Lambdas_chain = np.append(Lambdas_chain,a1)
        SA_chain = np.append(SA_chain,a2)
        syst_chain = np.append(syst_chain,sy)
        lk_chain = np.append(lk_chain,lk)
        naccept += na

        with open(fname_out,'a') as csvfile:
            chainwriter = csv.writer(csvfile)
            chainwriter.writerow(a0)
            chainwriter.writerow(a1)
            chainwriter.writerow(a2)
            chainwriter.writerow(np.array(sy)[:,0])
            chainwriter.writerow(np.array(sy)[:,1])
            chainwriter.writerow(np.array(sy)[:,2])
            chainwriter.writerow(np.array(sy)[:,3])
            chainwriter.writerow(np.array(sy)[:,4])
            chainwriter.writerow(np.array(sy)[:,5])
            chainwriter.writerow(lk)
        print 'Naccept: {}, percent: {}'.format(naccept,naccept/float(chunksize*(ichnk+1))/2.)

    # last chunk
    if lst_chunk != 0:
        start = time.time()
        a0,a1,a2,sy,lk,na = metropolis_systematics(data_mc,data_drt,startpos,startgrand,stepsize,lst_chunk)

        end = time.time()
        print 'This chunk: {} s'.format(end - start)

        Deltas_chain = np.append(Deltas_chain,a0)
        Lambdas_chain = np.append(Lambdas_chain,a1)
        SA_chain = np.append(SA_chain,a2)
        syst_chain = np.append(syst_chain,sy)
        lk_chain = np.append(lk_chain,lk)
        naccept += na

        with open(fname_out,'a') as csvfile:
            chainwriter = csv.writer(csvfile)
            chainwriter.writerow(a0)
            chainwriter.writerow(a1)
            chainwriter.writerow(a2)
            chainwriter.writerow(np.array(sy)[:,0])
            chainwriter.writerow(np.array(sy)[:,1])
            chainwriter.writerow(np.array(sy)[:,2])
            chainwriter.writerow(np.array(sy)[:,3])
            chainwriter.writerow(np.array(sy)[:,4])
            chainwriter.writerow(np.array(sy)[:,5])
            chainwriter.writerow(lk)
        print 'Naccept: {}, percent: {}'.format(naccept,naccept/float(chunksize*(ichnk+1) + lst_chunk))

    return Deltas_chain,Lambdas_chain,SA_chain,syst_chain,lk_chain,naccept

def metropolis_systematics(data_mc,data_drt,startpos,startgrand,stepsize,num_it):
 
    # num_it is the number of steps
    # num_syst is the number of systematic parameters
    num_syst = 6

    # make empty chains
    a0s_chain = np.zeros(num_it)
    a1s_chain = np.zeros(num_it)
    a2s_chain = np.zeros(num_it)
    syst_chain = np.zeros((num_it,num_syst))
    lk_chain = np.zeros(2*num_it)

    # current Gibbs random sample
    #grand = np.random.randn(num_syst)
    #grand = np.zeros(num_syst)

    # initialize the chains
    a0s_chain[0] = startpos[0]
    a1s_chain[0] = startpos[1]
    a2s_chain[0] = startpos[2]
    syst_chain[0] = startgrand
    lk_chain[0] = lk.lnprob(startpos,startgrand,data_mc,data_drt)


    # start walking
    naccept = 0
    for iw in range(1,num_it):

        #if iw % 10 == 0:
        #    print 'step: {}'.format(iw)

        ''' Axial FF step '''
        pos = a0s_chain[iw-1],a1s_chain[iw-1],a2s_chain[iw-1]

        prop = np.random.multivariate_normal(pos,stepsize)
        prop_ds = prop[0]
        prop_ls = prop[1]
        prop_sa = prop[2]
        prop_theta = prop_ds,prop_ls,prop_sa
        # evaluate at previous syst step
        prop_lk = lk.lnprob(prop_theta,syst_chain[iw-1],data_mc,data_drt)

        if prop_lk - lk_chain[2*iw-2] >= 0:
            # proposal is more likely than current position
            # accept step
            a0s_chain[iw] = prop_ds
            a1s_chain[iw] = prop_ls
            a2s_chain[iw] = prop_sa
            lk_chain[2*iw - 1] = prop_lk
            naccept += 1
        else:
            # proposal is less likely than current position
            # draw between zero and one
            draw = np.random.rand()
            # accept if ratio is greater than draw
            if prop_lk - lk_chain[iw-1] >= np.log(draw):
                # accept step
                a0s_chain[iw] = prop_ds
                a1s_chain[iw] = prop_ls
                a2s_chain[iw] = prop_sa
                lk_chain[2*iw - 1] = prop_lk
                naccept += 1
            else:
                # reject
                a0s_chain[iw] = a0s_chain[iw-1]
                a1s_chain[iw] = a1s_chain[iw-1]
                a2s_chain[iw] = a2s_chain[iw-1]
                lk_chain[2*iw - 1] = lk_chain[2*iw-2]


        ''' systematics step '''
        # Keep axial FF parameters the same
        # take step in systematic parameters
        prop_grand = np.random.randn(num_syst)
        #prop_grand = np.zeros(num_syst)
        # evaluate at current axial FF step
        pos = a0s_chain[iw],a1s_chain[iw],a2s_chain[iw]
        prop_lk = lk.lnprob(pos,prop_grand,data_mc,data_drt)
        #if prop_lk - lk_chain[2*iw - 1] >= 0:
        if True:
            # proposal is more likely than current position
            # accept step
            syst_chain[iw] = prop_grand
            lk_chain[2*iw] = prop_lk
            naccept += 1
        else:
            # proposal is less likely than current position
            # draw between zero and one
            draw = np.random.rand()
            # accept if ratio is greater than draw
            if prop_lk - lk_chain[iw-1] >= np.log(draw):
                # accept step
                syst_chain[iw] = prop_grand
                lk_chain[2*iw] = prop_lk
                naccept += 1
            else:
                # reject
                syst_chain[iw] = syst_chain[iw-1]
                lk_chain[2*iw] = lk_chain[2*iw-1]

    return a0s_chain,a1s_chain,a2s_chain,syst_chain,lk_chain,naccept


def metropolis_gibbs(data_mc,data_drt,startpos,startgrand,stepsize,num_it):
 
    # num_it is the number of gibbs (metropolis) steps
    # ==> total steps in chain is 2*num_it

    # make empty chains
    a0s_chain = np.zeros(2*num_it)
    a1s_chain = np.zeros(2*num_it)
    a2s_chain = np.zeros(2*num_it)
    syst_chain = np.zeros((2*num_it,num_syst))
    lk_chain = np.zeros(2*num_it)

    # current Gibbs random sample
    num_syst = 6
    grand = np.random.randn(num_syst)

    # initialize the chains
    a0s_chain[0] = startpos[0]
    a1s_chain[0] = startpos[1]
    a2s_chain[0] = startpos[2]
    syst_chain[0] = startgrand
    lk_chain[0] = lk.lnprob(startpos,startgrand,data_mc,data_drt)


    # start walking
    naccept = 0
    for iw in range(1,num_it+1):

        im = 2*iw-1  # metropolis step index
        ig = 2*iw    # gibbs step index
        if iw % 10 == 0:
            print 'step: {}'.format(iw)

        ''' Metropolis step '''
        pos = [a0s_chain[im-1],a1s_chain[im-1],a2s_chain[im-1]]

        prop = np.random.multivariate_normal(pos,stepsize)
        prop_ds = prop[0]
        prop_ls = prop[1]
        prop_sa = prop[2]
        prop_theta = prop_ds,prop_ls,prop_sa
        # evaluate at iw-1 (previous Gibbs step)
        prop_lk = lk.lnprob(prop_theta,syst_chain[iw-1],data_mc,data_drt)

        if prop_lk - lk_chain[im-1] >= 0:
            # proposal is more likely than current position
            # accept step
            a0s_chain[im] = prop_ds
            a1s_chain[im] = prop_ls
            a2s_chain[im] = prop_sa
            lk_chain[im] = prop_lk
            naccept += 1
        else:
            # proposal is less likely than current position
            # draw between zero and one
            draw = np.random.rand()
            # accept if ratio is greater than draw
            if prop_lk - lk_chain[im-1] >= np.log(draw):
                # accept step
                a0s_chain[im] = prop_ds
                a1s_chain[im] = prop_ls
                a2s_chain[im] = prop_sa
                lk_chain[im] = prop_lk
                naccept += 1
            else:
                # reject
                a0s_chain[im] = a0s_chain[im-1]
                a1s_chain[im] = a1s_chain[im-1]
                a2s_chain[im] = a2s_chain[im-1]
                lk_chain[im] = lk_chain[im-1]

        if iw == num_it:
            break

        ''' Gibbs step '''
        # Keep metropolis parameters the same
        # take step in nuisance parameters
        grand = np.random.randn(num_syst)
        # always accept
        a0s_chain[ig] = a0s_chain[im]
        a1s_chain[ig] = a1s_chain[im]
        a2s_chain[ig] = a2s_chain[im]
        
        prop_theta = a0s_chain[ig],a1s_chain[ig],a2s_chain[ig]
        lk_chain[ig] = lk.lnprob(prop_theta,data_mc,data_drt,grand)


    return a0s_chain,a1s_chain,a2s_chain,lk_chain,naccept


def metropolis(data_mc,data_dirt,data_drt,startpos,stepsize,num_it):
   
    # num_it is the number of metropolis steps

    # make empty chains
    Deltas_chain = np.zeros(num_it)
    Lambdas_chain = np.zeros(num_it)
    SA_chain = np.zeros(num_it)
    lk_chain = np.zeros(num_it)

    # current Gibbs random sample
    grand = np.random.randn(num_syst)

    # initialize the chains
    Deltas_chain[0] = startpos[0]
    Lambdas_chain[0] = startpos[1]
    SA_chain[0] = startpos[2]
    #lk_chain[0] = lk.lnprob(startpos,data_mc,data_drt,grand)
    # for now, no systematics
    lk_chain[0] = lk.lnprob(startpos,data_mc,data_drt,0.)

    # start walking
    naccept = 0
    for i in range(1,num_it):

        ''' Metropolis step '''
        pos = [Deltas_chain[im-1],Lambdas_chain[im-1],SA_chain[im-1]]

        prop = np.random.multivariate_normal(pos,stepsize)
        prop_ds = prop[0]
        prop_ls = prop[1]
        prop_sa = prop[2]
        prop_theta = prop_ds,prop_ls,prop_sa
        # evaluate with systematic sample
        #grand = np.random.randn(6)
        #prop_lk = lk.lnprob(prop_theta,data_mc,data_drt,grand)
        # for now, no systematics
        prop_lk = lk.lnprob(prop_theta,data_mc,data_drt,0.)

        if prop_lk - lk_chain[i-1] >= 0:
            # proposal is more likely than current position
            # accept step
            Deltas_chain[i] = prop_ds
            Lambdas_chain[i] = prop_ls
            SA_chain[i] = prop_sa
            lk_chain[i] = prop_lk
            naccept += 1
        else:
            # proposal is less likely than current position
            # draw between zero and one
            draw = np.random.rand()
            # accept if ratio is greater than draw
            if prop_lk - lk_chain[i-1] >= np.log(draw):
                # accept step
                Deltas_chain[i] = prop_ds
                Lambdas_chain[i] = prop_ls
                SA_chain[i] = prop_sa
                lk_chain[i] = prop_lk
                naccept += 1
            else:
                # reject
                Deltas_chain[i] = Deltas_chain[i-1]
                Lambdas_chain[i] = Lambdas_chain[i-1]
                SA_chain[i] = SA_chain[i-1]
                lk_chain[i] = lk_chain[i-1]

    return Deltas_chain,Lambdas_chain,SA_chain,lk_chain,naccept

