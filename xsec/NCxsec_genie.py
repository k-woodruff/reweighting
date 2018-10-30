import numpy as np
import scipy.stats as st

Mp        = 0.938272      # GeV/c^2
Mv        = 0.840         # GeV/c^2
sin2theta = 0.2277
GF        = 1.16639E-5    # cm/GeV
gA        = -1.267
mup       = 2.7930
mun       = -1.913042
kPi       = 3.1415927     # from genie

def NCpxsec(E,Q2,GaS,MA,nusign,nucsign):
    E2   = E**2
    M    = Mp
    M2   = M**2
    fMv2 = Mv**2
    fMa2 = MA**2
    qmv2 = (1+Q2/fMv2)**2
    qma2 = (1+Q2/fMa2)**2
    fMuP = mup
    fMuN = mun
    fFa0 = gA
    fEta = GaS/gA
    GF2  = GF**2

    fkAlpha = 1.-2.*sin2theta
    fkGamma = -0.66666667*sin2theta

    #-- compute isoscalar form factor terms
    Ge0 = 1.5 * fkGamma / qmv2
    Gm0 = 1.5 * fkGamma * (fMuP+fMuN) / qmv2

    #-- compute isovector form factor terms
    Ge1 = 0.5 * fkAlpha / qmv2
    Gm1 = 0.5 * fkAlpha * (fMuP-fMuN) / qmv2
    Ga1 = -0.5 * fFa0 * (1 + (nucsign) * fEta) / qma2
    
    #-- compute form factors
    Ge  = Ge0 + (nucsign) * Ge1
    Gm  = Gm0 + (nucsign) * Gm1
    Ga  = (nucsign) * Ga1
    Ge2 = Ge**2
    Gm2 = Gm**2
    Ga2 = Ga**2
    
    #-- compute the free nucleon cross section
    tau   = 0.25 * Q2/M2
    fa    = 1-M*tau/E
    fa2   = fa**2
    fb    = tau*(tau+1)*M2/E2
    A     = (Ge2/(1+tau))           * (fa2-fb)
    B     = (Ga2 + tau*Gm2/(1+tau)) * (fa2+fb)
    C     = 4*tau*(M/E)*Gm*Ga       * fa
    xsec0 = 0.5*GF2/kPi;
    xsec  = xsec0 * (A + B + (nusign)*C)

    return xsec

