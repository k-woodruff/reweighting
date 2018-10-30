import numpy as np
import pandas as pd
import xsec._ncformfactor as ncff
import xsec._zeformfactor as zeff
import xsec.NCxsec_genie as ncg

from multiprocessing import Pool
from functools import partial

def rwgt_ratio_mc(deltas,lambdaA,sA,Q2vec,data_mc,grand):

    ''' for reweighting and calculating the ratio need 4 groups '''
    # 1: true NCE reco'd as NCE
    # 2: other reco'd as NCE
    ''' need weights for 1,2 '''

    data_1 = data_mc[(data_mc.mc_ccnc == 1)&(data_mc.mc_mode == 0)]
    data_2 = data_mc[(data_mc.mc_ccnc != 1)|(data_mc.mc_mode != 0)]

    nusign_1 = np.sign(data_1['mc_nuPDG'].values)
    nusign_2 = np.sign(data_2['mc_nuPDG'].values)

    nucsign_1 = np.where(data_1['mc_hitnuc'] == 2112, -1,1)

    num_processes = 6
    pool = Pool(num_processes)
    theta_args = (deltas,lambdaA,sA)

    rwt_list_1 = map(list,zip(*[data_1.mc_Enu_nucRF.values.tolist(),data_1.mc_Q2.values.tolist(),nusign_1.tolist(),nucsign_1.tolist()]))

    axial_weight_1 = pool.map(partial(get_ncxsec,theta_args=theta_args), rwt_list_1)

    pool.close()
    pool.join()

    # systematic proposals
    # DIC
    dicrand = np.abs(grand[0])
    dic_1 = (data_1['syst_dic_wt_m']-1)*dicrand + 1.
    dic_2 = (data_2['syst_dic_wt_m']-1)*dicrand + 1.
    # SPE
    sperand = grand[1]
    spe_1 = (data_1['syst_spe_wt']-1)*sperand + 1.
    spe_2 = (data_2['syst_spe_wt']-1)*sperand + 1.
    # flux
    fluxrand = grand[2]
    flux_1 = (data_1['syst_flux_wt']-1)*fluxrand + 1.
    flux_2 = (data_2['syst_flux_wt']-1)*fluxrand + 1.
    # MEC
    mecrand = grand[3]
    mec_1 = (data_1['syst_mec_wt']-1)*mecrand + 1.
    mec_2 = (data_2['syst_mec_wt']-1)*mecrand + 1.
    # Pauli block
    kfrand = grand[4]
    kf_1 = (data_1['syst_kf_wt_m']-1)*kfrand + 1.
    kf_2 = (data_2['syst_kf_wt_m']-1)*kfrand + 1.
    # Dirt
    drtrand = grand[5]
    drt_1 = (data_1['syst_dirt_wt']-1)*drtrand + 1.
    drt_2 = (data_2['syst_dirt_wt']-1)*drtrand + 1.

    ln_weight_gibbs_1 = np.log(dic_1) + np.log(spe_1) + np.log(flux_1) \
                      + np.log(mec_1) \
                      + np.log(drt_1) \
                      + np.log(kf_1)
    ln_weight_gibbs_2 = np.log(dic_2) + np.log(spe_2) + np.log(flux_2) \
                      + np.log(mec_2) + np.log(drt_2) + np.log(kf_2)

    ln_weight_1 = np.log(axial_weight_1) + ln_weight_gibbs_1
    ln_weight_2 = np.log(data_2['ccqe_wt'].values) + np.log(data_2['mec_wt'].values) + ln_weight_gibbs_2

    events_nc_1,_ = np.histogram(data_1['reco_Q2'], bins=Q2vec, weights=np.exp(ln_weight_1))
    events_nc_2,_ = np.histogram(data_2['reco_Q2'], bins=Q2vec, weights=np.exp(ln_weight_2))
    events_nc = events_nc_1 + events_nc_2

    ncsig = np.sqrt(events_nc)

    return events_nc,ncsig

def get_ncxsec(list_args,theta_args):

    df_ds = -0.15204
    df_MA = 0.990

    deltas,lambdaA,sA = theta_args

    wgt_num = zeff.sigplus(list_args[0],list_args[1],deltas,lambdaA,sA,(int)(list_args[2]),(int)(list_args[3]))

    # use genie/ahrens xsec for denominator
    wgt_den = ncg.NCpxsec(list_args[0],list_args[1],df_ds,df_MA,(int)(list_args[2]),(int)(list_args[3]))

    return wgt_num/wgt_den

