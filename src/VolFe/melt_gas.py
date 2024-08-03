# melt_gas.py

import pandas as pd
import numpy as np
import math
import warnings as w

#import model_dependent_variables as mdv
import VolFe.model_dependent_variables as mdv
import VolFe.equilibrium_equations as eq

### TO SORT ###
# species X solubility
# species X fugacity coefficient

################
### Contents ###
################

# Fugacity, mole fraction, partial pressure
# Speciation 
# Converting gas and melt compositions
# Volume and density

#################################################################################################################################
########################################### FUGACITY, MOLE FRACTION, PARTIAL PRESSURE ###########################################
#################################################################################################################################

def f_H2O(PT,melt_wf,models):
    
    """ 
    Calculates fugacity of water in the vapor from concentration of water in the melt


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions: pressure (bars) as "P" and temperature ('C) as "T". 
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition. WHAT.
        
    species: pandas.DataFrame
        Dataframe of species.csv file.
    
    models: pandas.DataFrame
        Dataframe of models.csv file.

    Returns
    -------
    fH2O in bars as <class 'mpfr'>

    """

    Hspeciation = models.loc["Hspeciation","option"]
    if Hspeciation in ["none",'none+ideal',"none+regular"]: # fH2O = xmH2OT^2/CH2O
        value = ((xm_H2OT_so(melt_wf))**2.0)/mdv.C_H2O(PT,melt_wf,models)
        return value
    elif Hspeciation == "linear": # fH2O = xmH2OT/CH2O
        value = xm_H2OT_so(melt_wf)/mdv.C_H2O(PT,melt_wf,models)
        return value        
    else: # regular or ideal: fH2O = xmH2Omol/CH2O
        value = xm_H2Omol_so(PT,melt_wf,models)/mdv.C_H2O(PT,melt_wf,models)
        return value
    
def f_CO2(PT,melt_wf,models):
    CO3model = models.loc["carbon dioxide","option"]
    wm_CO2 = 100.*melt_wf['CO2']
    if CO3model == "Shishkina14": # wtppm CO2 modified from Shishkina et al. (2014) Chem. Geol. 388:112-129
        f = (wm_CO2*10000.0)/mdv.C_CO3(PT,melt_wf,models)
    else: # fCO2 = xmCO2/C_CO3
        f = xm_CO2_so(melt_wf)/mdv.C_CO3(PT,melt_wf,models)
    return f

def f_S2(PT,melt_wf,models): # wtppm S2- NOT mole fraction due to parameterisation by O'Neill (2020)
    K = mdv.C_S(PT,melt_wf,models)/1000000.
    fS2 = ((melt_wf["S2-"]/K)**2.)*mdv.f_O2(PT,melt_wf,models)
    return fS2
    
def f_H2(PT,melt_wf,models):
    K = mdv.KHOg(PT,models)
    return f_H2O(PT,melt_wf,models)/(K*mdv.f_O2(PT,melt_wf,models)**0.5)

def f_CO(PT,melt_wf,models):
    K = mdv.KCOg(PT,models)
    return f_CO2(PT,melt_wf,models)/(K*mdv.f_O2(PT,melt_wf,models)**0.5)

def f_H2S(PT,melt_wf,models):
    K = mdv.KHOSg(PT,models)
    return (K*f_S2(PT,melt_wf,models)**0.5*f_H2O(PT,melt_wf,models))/mdv.f_O2(PT,melt_wf,models)**0.5

def f_SO2(PT,melt_wf,models):
    K = mdv.KOSg(PT,models)
    return K*mdv.f_O2(PT,melt_wf,models)*f_S2(PT,melt_wf,models)**0.5

def f_SO3(PT,melt_wf,models):
    K = mdv.KOSg2(PT,models)
    return K*(mdv.f_O2(PT,melt_wf,models))**1.5*(f_S2(PT,melt_wf,models))**0.5

def f_CH4(PT,melt_wf,models):
    K = mdv.KCOHg(PT,models)
    return (f_CO2(PT,melt_wf,models)*f_H2O(PT,melt_wf,models)**2.0)/(K*mdv.f_O2(PT,melt_wf,models)**2.0)

def f_OCS(PT,melt_wf,models):
    OCSmodel = models.loc["carbonylsulfide","option"]
    K = mdv.KOCSg(PT,models)
    if OCSmodel == "COHS":
        if f_H2O(PT,melt_wf,models) > 0.0:
            return (f_CO2(PT,melt_wf,models)*f_H2S(PT,melt_wf,models))/(f_H2O(PT,melt_wf,models)*K)
        else:
            return 0.0
    else:
        if f_CO2(PT,melt_wf,models) > 0.0:
            return ((f_CO(PT,melt_wf,models)**3.0)*f_SO2(PT,melt_wf,models))/((f_CO2(PT,melt_wf,models)**2.0)*K)
        else:
            return 0.0
        
def f_X(PT,melt_wf,models):
    K = mdv.C_X(PT,melt_wf,models)/1000000.
    f_X = melt_wf["XT"]/K
    return f_X


#######################
### oxygen fugacity ###
#######################

# buffers
def fO22Dbuffer(PT,fO2,buffer,models):
    if buffer == "NNO":
        return math.log10(fO2) - mdv.NNO(PT,models)
    elif buffer == "FMQ":
        return math.log10(fO2) - mdv.FMQ(PT,models)
def Dbuffer2fO2(PT,D,buffer,model):
    if buffer == "NNO":
        return 10.0**(D + mdv.NNO(PT,model))
    elif buffer == "FMQ":
        return 10.0**(D + mdv.FMQ(PT,model))
        
def S6S2_2_fO2(S62,melt_wf,PT,models):

    def calc_fO2(S62,melt_wf,PT,models):
        CSO4 = mdv.C_SO4(PT,melt_wf,models)/1000000.
        CS = mdv.C_S(PT,melt_wf,models)/1000000.
        if models.loc["H2S_m","option"] == "False":
            fO2 = ((S62*CS)/CSO4)**0.5
        elif models.loc["H2S_m","option"] == "True":
            KHS = mdv.KHOSg(PT,models)
            CH2S = ((mdv.C_H2S(PT,melt_wf,models)/1000000.)/mdv.species.loc['H2S','M'])*mdv.species.loc['S','M']
            CH2OT = mdv.C_H2O(PT,melt_wf,models)
            xmH2O = xm_H2OT_so(melt_wf)
            W = (CSO4/((KHS*CH2S*(xmH2O**2./CH2OT)) + CS))
            fO2 = (S62/W)**0.5
        return fO2, xmH2O

    melt_wf['Fe3FeT'] = 0.1
    fO2, xmH2O = calc_fO2(S62,melt_wf,PT,models)
    Fe3FeT = mdv.fO22Fe3FeT(fO2,PT,melt_wf,models)
    delta1 = abs(melt_wf['Fe3FeT'] - Fe3FeT)
    conc = eq.melt_speciation(PT,melt_wf,models,1.,1.e-9)
    melt_wf["H2OT"] = conc["wm_H2O"]
    melt_wf["CO2"] = conc["wm_CO2"]
    while delta1 > 0.0001 or delta2 > 1.e-9:
        melt_wf['Fe3FeT'] = Fe3FeT
        fO2, xmH2O = calc_fO2(S62,melt_wf,PT,models)
        Fe3FeT = mdv.fO22Fe3FeT(fO2,PT,melt_wf,models)
        delta1 = abs(melt_wf['Fe3FeT'] - Fe3FeT)
        conc = eq.melt_speciation(PT,melt_wf,models,1.,1.e-9)
        melt_wf["H2OT"] = conc["wm_H2O"]
        melt_wf["CO2"] = conc["wm_CO2"]
        delta2 = abs(xmH2O - xm_H2OT_so(melt_wf))
    else:
        fO2, xmH2O = calc_fO2(S62,melt_wf,PT,models)
    
    #S6p = melt_wf["ST"]*(S62/(1.+S62))
    #S2mT = melt_wf["ST"] - S6p
    #S_H2S_S2m = (mdv.KHOSg(PT,models)*(mdv.C_H2S(PT,melt_wf,models)/1000000.)*(xmH2O**2./mdv.C_H2O(PT,melt_wf,models)))/(mdv.C_S(PT,melt_wf,models)/1000000.)
    #S_H2S = S2mT*(S_H2S_S2m/(1.+S_H2S_S2m))
    #H2S = (S_H2S/mdv.species.loc['S','M'])*mdv.species.loc['H2S','M']
    #S2m = S2mT - S_H2S
    return fO2
    
    
########################        
### partial pressure ###
########################

def p_H2(PT,melt_wf,models):
    return f_H2(PT,melt_wf,models)/mdv.y_H2(PT,models)
def p_H2O(PT,melt_wf,models):
    return f_H2O(PT,melt_wf,models)/mdv.y_H2O(PT,models)
def p_O2(PT,melt_wf,models):
    return mdv.f_O2(PT,melt_wf,models)/mdv.y_O2(PT,models)
def p_SO2(PT,melt_wf,models):
    return f_SO2(PT,melt_wf,models)/mdv.y_SO2(PT,models)
def p_S2(PT,melt_wf,models):
    return f_S2(PT,melt_wf,models)/mdv.y_S2(PT,models)
def p_H2S(PT,melt_wf,models):
    return f_H2S(PT,melt_wf,models)/mdv.y_H2S(PT,models)
def p_CO2(PT,melt_wf,models):
    return f_CO2(PT,melt_wf,models)/mdv.y_CO2(PT,models)
def p_CO(PT,melt_wf,models):
    return f_CO(PT,melt_wf,models)/mdv.y_CO(PT,models)
def p_CH4(PT,melt_wf,models):
    return f_CH4(PT,melt_wf,models)/mdv.y_CH4(PT,models)
def p_OCS(PT,melt_wf,models):
    return f_OCS(PT,melt_wf,models)/mdv.y_OCS(PT,models)
def p_X(PT,melt_wf,models):
    return f_X(PT,melt_wf,models)/mdv.y_X(PT,models)

def p_tot(PT,melt_wf,models):
    return p_H2(PT,melt_wf,models) + p_H2O(PT,melt_wf,models) + p_O2(PT,melt_wf,models) + p_SO2(PT,melt_wf,models) + p_S2(PT,melt_wf,models) + p_H2S(PT,melt_wf,models) + p_CO2(PT,melt_wf,models) + p_CO(PT,melt_wf,models) + p_CH4(PT,melt_wf,models) + p_OCS(PT,melt_wf,models) + p_X(PT,melt_wf,models)


############################       
### vapor molar fraction ###
############################

def xg_H2(PT,melt_wf,models):
    return p_H2(PT,melt_wf,models)/PT['P']
def xg_H2O(PT,melt_wf,models):
    return p_H2O(PT,melt_wf,models)/PT['P']
def xg_O2(PT,melt_wf,models):
    return p_O2(PT,melt_wf,models)/PT['P']
def xg_SO2(PT,melt_wf,models):
    return p_SO2(PT,melt_wf,models)/PT['P']
def xg_S2(PT,melt_wf,models):
    return p_S2(PT,melt_wf,models)/PT['P']
def xg_H2S(PT,melt_wf,models):
    return p_H2S(PT,melt_wf,models)/PT['P']
def xg_CO2(PT,melt_wf,models):
    return p_CO2(PT,melt_wf,models)/PT['P']
def xg_CO(PT,melt_wf,models):
    return p_CO(PT,melt_wf,models)/PT['P']
def xg_CH4(PT,melt_wf,models):
    return p_CH4(PT,melt_wf,models)/PT['P']
def xg_OCS(PT,melt_wf,models):
    return p_OCS(PT,melt_wf,models)/PT['P']
def xg_X(PT,melt_wf,models):
    return p_X(PT,melt_wf,models)/PT['P']
def Xg_tot(PT,melt_wf,models):
    species_X = models.loc["species X","option"]
    Xg_t = xg_CO2(PT,melt_wf,models)*mdv.species.loc["CO2","M"] + xg_CO(PT,melt_wf,models)*mdv.species.loc["CO","M"] + xg_O2(PT,melt_wf,models)*mdv.species.loc["O2","M"] + xg_H2O(PT,melt_wf,models)*mdv.species.loc["H2O","M"] + xg_H2(PT,melt_wf,models)*mdv.species.loc["H2","M"] + xg_CH4(PT,melt_wf,models)*mdv.species.loc["CH4","M"] + xg_SO2(PT,melt_wf,models)*mdv.species.loc["SO2","M"] + xg_S2(PT,melt_wf,models)*mdv.species.loc["S2","M"] + xg_H2S(PT,melt_wf,models)*mdv.species.loc["H2S","M"] + xg_OCS(PT,melt_wf,models)*mdv.species.loc["OCS","M"] + xg_X(PT,melt_wf,models)*mdv.species.loc[species_X,"M"]
    return Xg_t


##########################
### melt mole fraction ###
##########################

# totals
def wm_vol(melt_wf): # wt% total volatiles in the melt
    wm_H2OT = 100.*melt_wf["H2OT"]
    wm_CO2 = 100.*melt_wf["CO2"]
    return wm_H2OT + wm_CO2 #+ wm_S(wm_ST) + wm_SO3(wm_ST)
def wm_nvol(melt_wf): # wt% total of non-volatiles in the melt
    return 100. - wm_vol(melt_wf)

# molecular mass on a singular oxygen basis
def M_m_SO(melt_wf):
    # single oxide, no volatiles, all Fe as FeOT
    melt_comp = melt_single_O(melt_wf,"no","no")
    M_m = 1./melt_comp["Xmtot"]
    return M_m

# molecular mass on a oxide basis
def M_m_ox(melt_wf,models): # no volatiles, all Fe as FeOT
    melt_comp = melt_mole_fraction(melt_wf,models,"no","no") 
    M_m = 1./melt_comp["mol_tot"]
    return M_m
    
# Number of moles in the melt
def Xm_H2OT(melt_wf):
    wm_H2OT = 100.*melt_wf['H2OT']
    return wm_H2OT/mdv.species.loc["H2O","M"]
def Xm_CO2(melt_wf):
    wm_CO2 = 100.*melt_wf['CO2']
    return wm_CO2/mdv.species.loc["CO2","M"]

# Mole fraction in the melt based on mixing between volatile-free melt on a singular oxygen basis and volatiles
def Xm_m_so(melt_wf): # singular oxygen basis
    return wm_nvol(melt_wf)/M_m_SO(melt_wf)    
def Xm_tot_so(melt_wf):
    return Xm_H2OT(melt_wf) + Xm_CO2(melt_wf) + Xm_m_so(melt_wf)
def xm_H2OT_so(melt_wf):
    return Xm_H2OT(melt_wf)/Xm_tot_so(melt_wf)
def xm_CO2_so(melt_wf):
    return Xm_CO2(melt_wf)/Xm_tot_so(melt_wf)
def xm_melt_so(melt_wf):
    return Xm_m_so(melt_wf)/Xm_tot_so(melt_wf)
def Xm_t_so(melt_wf):
    return xm_H2OT_so(melt_wf)*mdv.species.loc["H2O","M"] + xm_CO2_so(melt_wf)*mdv.species.loc["CO2","M"] + xm_melt_so(melt_wf)*M_m_SO(melt_wf)

################################################################################################################################
########################################################## SPECIATION ##########################################################
################################################################################################################################

def ratio2overtotal(x):
    return x/x+1.

def overtotal2ratio(x):
    return x/(1.-x)

###########################
### hydrogen speciation ###
###########################

def xm_OH_so(PT,melt_wf,models):
    Hspeciation = models.loc["Hspeciation","option"]
    Hspeccomp = models.loc["Hspeccomp","option"]
    T_K = PT['T']+273.15
    wm_H2OT = 100.*melt_wf['H2OT']
    
    Z = xm_H2OT_so(melt_wf)
    
    def f(A, B, C, Z, x): # regular solution model rearranged to equal 0 to solve for xm_OH
        return (A + B*x + C*(Z - 0.5*x) + log((x**2.0)/((Z - 0.5*x)*(1.0 - Z - 0.5*x))))

    def df(B, C, Z, x): # derivative of above
        return (B - 0.5*C - (2.0/x) + (0.5/(Z-0.5*x) + (0.5/(1.0-Z-0.5*x))))

    def dx(A,B,C,x):
        return (abs(0-f(A, B, C, Z, x)))
 
    def nr(A,B,C,x0,e1):
        delta1 = dx(A,B,C,x0)
        while delta1 > e1:
            x0 = x0 - f(A, B, C, Z, x0)/df(B, C, Z, x0)
            delta1 = dx(A,B,C,x0)
        return x0
    
    if Z == 0.0: 
        return 0.0
    elif Hspeciation in ["ideal",'none+ideal']: # ideal mixture from Silver & Stolper (1985) J. Geol. 93(2) 161-177
        K = mdv.KHOm(PT,melt_wf,models)
        return (0.5 - (0.25 - ((K - 4.0)/K)*(Z - Z**2.0))**0.5)/((K - 4.0)/(2.0*K)) 
    elif Hspeciation == "none": # all OH-
        return 0.0
    
    elif Hspeciation in ["regular","none+regular"]: # else use regular solution model from Silver & Stolper (1989) J.Pet. 30(3)667-709
        tol = 0.000001 #tolerance
        K = mdv.KHOm(PT,melt_wf,models)
        x0 = (0.5 - (0.25 - ((K - 4.0)/K)*(Z - Z**2.0))**0.5)/((K - 4.0)/(2.0*K))
        R = 83.15 #cm3.bar/mol.K
        A, B, C = mdv.KregH2O(PT,melt_wf,models)
        return nr(A,B,C,x0,tol)

def xm_H2Omol_so(PT,melt_wf,models):
    Z = xm_H2OT_so(melt_wf)
    return Z - 0.5*xm_OH_so(PT,melt_wf,models)

def wm_H2Omol_OH(PT,melt_wf,models): # wt fraction
    H2Omol = xm_H2Omol_so(PT,melt_wf,models)*mdv.species.loc["H2O","M"]
    OH_H2O = xm_OH_so(PT,melt_wf,models)*0.5*mdv.species.loc["H2O","M"]
    melt = xm_melt_so(melt_wf)*M_m_SO(melt_wf)
    CO2T = xm_CO2_so(melt_wf)*mdv.species.loc["CO2","M"]
    wm_H2Omol = H2Omol/(H2Omol + OH_H2O + melt + CO2T)
    wm_OH_H2O = OH_H2O/(H2Omol + OH_H2O + melt + CO2T)
    wm_OH = ((wm_OH_H2O/mdv.species.loc["H2O","M"])*2.)*mdv.species.loc["OH","M"]
    return wm_H2Omol, wm_OH

#########################
### carbon speciation ###
#########################

def xm_CO32_CO2mol(PT,melt_wf,models): # mole fraction
    xm_CO2T = xm_CO2_so(melt_wf)
    xm_H2OT = xm_H2OT_so(melt_wf)
    if models.loc["Cspeccomp","option"] == "Basalt":
        xm_CO32 = xm_CO2T
        xm_CO2mol = 0.
    elif models.loc["Cspeccomp","option"] == "Rhyolite":
        xm_CO32 = 0.
        xm_CO2mol = xm_CO2T
    else: 
        K = mdv.KCOm(PT,melt_wf,models)
        xm_CO32 = (K*xm_CO2T*(1.-xm_CO2T-xm_H2OT))/(1.+ K*(1.-xm_CO2T-xm_H2OT))
        xm_CO2mol = xm_CO2T - xm_CO32
    return xm_CO32, xm_CO2mol

def wm_CO32_CO2mol(PT,melt_wf,models): # wt fraction
    xm_CO32, xm_CO2mol = xm_CO32_CO2mol(PT,melt_wf,models)
    CO2mol = xm_CO2mol*mdv.species.loc["CO2","M"]
    CO2carb = xm_CO32*mdv.species.loc["CO2","M"]
    melt = xm_melt_so(melt_wf)*M_m_SO(melt_wf)
    H2OT = xm_H2OT_so(melt_wf)*mdv.species.loc["H2O","M"]
    wm_CO2mol = CO2mol/(CO2mol + CO2carb + H2OT + melt)
    wm_CO2carb = CO2carb/(CO2mol + CO2carb + H2OT + melt)
    return wm_CO2carb, wm_CO2mol 

##########################
### sulfur speciation ###
##########################

def S6S2(PT,melt_wf,models):
    T_K = PT['T']+273.15
    model = models.loc["sulfate","option"]
    if model == "Nash19":
        A, B = mdv.S_Nash19_terms(PT)
        result = 10.**(A*math.log10(Fe3Fe2(melt_wf)) + B)
    else:
        CSO4 = mdv.C_SO4(PT,melt_wf,models)
        CS = mdv.C_S(PT,melt_wf,models)
        fO2 = mdv.f_O2(PT,melt_wf,models)
        if models.loc["H2S_m","option"] == "False":
            result = (CSO4/CS)*fO2**2.
        elif models.loc["H2S_m","option"] == "True":
            KHS = mdv.KHOSg(PT,models)
            CH2S = mdv.C_H2S(PT,melt_wf,models)
            CH2OT = mdv.C_H2O(PT,melt_wf,models)
            xmH2O = xm_H2OT_so(melt_wf)
            result = (CSO4/((KHS*CH2S*(xmH2O**2./CH2OT)) + CS))*fO2**2.
    return result
    
def S6ST(PT,melt_wf,models):
    S6S2_ = S6S2(PT,melt_wf,models)
    return S6S2_/(S6S2_+1.)

def wm_S(PT,melt_wf,models):
    wm_ST = 100.*melt_wf['ST']
    S6ST_ = S6ST(PT,melt_wf,models)
    return wm_ST*(1.0-S6ST_)

def wm_SO3(PT,melt_wf,models):
    wm_ST = 100.*melt_wf['ST']
    S6ST_ = S6ST(PT,melt_wf,models)    
    return ((wm_ST*S6ST_)/mdv.species.loc["S","M"])*mdv.species.loc["SO3","M"]


#######################
### iron speciation ###
#######################

def Fe3Fe2(melt_wf):
    Fe3FeT = melt_wf['Fe3FeT']
    return Fe3FeT/(1.0 - Fe3FeT)

def Wm_FeT(melt_wf):
    if melt_wf["FeOT"] > 0.0:
        return (melt_wf["FeOT"]/mdv.species.loc["FeO","M"])*mdv.species.loc["Fe","M"]
    elif melt_wf["Fe2O3T"] > 0.0:
        return (melt_wf["Fe2O3T"]/mdv.species.loc["Fe2O3","M"])*mdv.species.loc["Fe","M"]
    else:
        return ((melt_wf["FeO"]/mdv.species.loc["FeO","M"]) + (melt_wf["Fe2O3"]/mdv.species.loc["Fe2O3","M"]))*mdv.species.loc["Fe","M"]

def Wm_FeO(melt_wf):
    Fe3FeT = melt_wf['Fe3FeT']
    return (Wm_FeT(melt_wf)/mdv.species.loc["Fe","M"])*(1.0-Fe3FeT)*mdv.species.loc["FeO","M"]

def Wm_Fe2O3(melt_wf):
    Fe3FeT = melt_wf['Fe3FeT']
    return (Wm_FeT(melt_wf)/mdv.species.loc["Fe","M"])*Fe3FeT*mdv.species.loc["Fe2O3","M"]

def Wm_FeOT(melt_wf):
    return (Wm_FeT(melt_wf)/mdv.species.loc["Fe","M"])*mdv.species.loc["FeO","M"]

def wm_Fe_nv(melt_wf): # no volatiles
    Wm_tot = melt_wf["SiO2"] + melt_wf["TiO2"] + melt_wf["Al2O3"] + melt_wf["MnO"] + melt_wf["MgO"] + melt_wf["MnO"] + melt_wf["CaO"] + melt_wf["Na2O"] + melt_wf["K2O"] + melt_wf["P2O5"] + Wm_FeO(melt_wf) + Wm_Fe2O3(melt_wf)
    FeT = mdv.species.loc["Fe","M"]*((2.0*Wm_Fe2O3(melt_wf)/mdv.species.loc["Fe2O3","M"]) + (Wm_FeO(melt_wf)/mdv.species.loc["FeO","M"]))
    return 100.0*FeT/Wm_tot

def Fe3FeT_i(PT,melt_wf,models):
    model = models.loc["fO2","option"]
    T_K = PT['T']+273.15
    
    if model == "buffered":
        fO2 = 10**(melt_wf["logfO2_i"])
        return mdv.fO22Fe3FeT(fO2,PT,melt_wf,models)
    else:
        if pd.isnull(melt_wf["Fe3FeT_i"]) == False:
            if pd.isnull(melt_wf["logfO2_i"]) == False or pd.isnull(melt_wf["DNNO"]) == False or pd.isnull(melt_wf["DFMQ"]) == False or pd.isnull(melt_wf["S6ST_i"]) == False or pd.isnull(melt_wf["Fe2O3"]) == False or pd.isnull(melt_wf["FeO"]) == False:
                w.warn('you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered Fe3+/FeT')
            return melt_wf["Fe3FeT_i"]
        elif pd.isnull(melt_wf["logfO2_i"]) == False:
            if pd.isnull(melt_wf["DNNO"]) == False or pd.isnull(melt_wf["DFMQ"]) == False or pd.isnull(melt_wf["S6ST_i"]) == False or pd.isnull(melt_wf["Fe2O3"]) == False or pd.isnull(melt_wf["FeO"]) == False:
                w.warn('you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered log(fO2)')
            fO2 = 10.0**(melt_wf["logfO2_i"])
            return mdv.fO22Fe3FeT(fO2,PT,melt_wf,models)
        elif pd.isnull(melt_wf["DNNO"]) == False:
            if pd.isnull(melt_wf["DFMQ"]) == False or pd.isnull(melt_wf["S6ST_i"]) == False or pd.isnull(melt_wf["Fe2O3"]) == False or pd.isnull(melt_wf["FeO"]) == False:
                w.warn('you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered DNNO')
            D = melt_wf["DNNO"]
            fO2 = Dbuffer2fO2(PT,D,"NNO",models)
            return mdv.fO22Fe3FeT(fO2,PT,melt_wf,models)
        elif pd.isnull(melt_wf["DFMQ"]) == False:
            if pd.isnull(melt_wf["S6ST_i"]) == False or pd.isnull(melt_wf["Fe2O3"]) == False or pd.isnull(melt_wf["FeO"]) == False:
                w.warn('you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered DFMQ')
            D = melt_wf["DFMQ"]
            fO2 = Dbuffer2fO2(PT,D,"FMQ",models)
            return mdv.fO22Fe3FeT(fO2,PT,melt_wf,models)
        elif pd.isnull(melt_wf["S6ST_i"]) == False:
            if pd.isnull(melt_wf["Fe2O3"]) == False or pd.isnull(melt_wf["FeO"]) == False:
                w.warn('you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered S6+/ST')
            S6T = melt_wf["S6ST_i"]
            S62 = overtotal2ratio(S6T)
            fO2 = S6S2_2_fO2(S62,melt_wf,PT,models)
            return mdv.fO22Fe3FeT(fO2,PT,melt_wf,models)
        else:
            return ((2.0*melt_wf["Fe2O3"])/mdv.species.loc["Fe2O3","M"])/(((2.0*melt_wf["Fe2O3"])/mdv.species.loc["Fe2O3","M"]) + (melt_wf["FeO"]/mdv.species.loc["FeO","M"]))
        
    
#################################################################################################################################
############################################## converting gas and melt compositions #############################################
#################################################################################################################################

# C/S ratio of the vapor
def gas_CS(PT,melt_wf,models):
    S_values = ((2.*xg_S2(PT,melt_wf,models))+xg_SO2(PT,melt_wf,models)+xg_OCS(PT,melt_wf,models)+xg_H2S(PT,melt_wf,models))
    if S_values > 0:
        xgCS = (xg_CO(PT,melt_wf,models)+xg_CO2(PT,melt_wf,models)+xg_OCS(PT,melt_wf,models)+xg_CH4(PT,melt_wf,models))/S_values
    else:
        xgCS = ""
    return xgCS

def gas_CS_alt(xg):
    S_values = ((2.*xg["S2"])+xg["SO2"]+xg["OCS"]+xg["H2S"])
    if S_values > 0:
        xgCS = (xg["CO"]+xg["CO2"]+xg["OCS"]+xg["CH4"])/S_values
    else:
        xgCS = ""
    return xgCS

# all carbon as CO2 and all hydrogen as H2O
def melt_H2O_CO2_eq(melt_wf):
    wmCO2eq = melt_wf["CO2"] + mdv.species.loc["CO2","M"]*((melt_wf["CO"]/mdv.species.loc["CO","M"])+(melt_wf["CH4"]/mdv.species.loc["CH4","M"]))
    wmH2Oeq = melt_wf["H2OT"] + mdv.species.loc["H2O","M"]*((melt_wf["H2"]/mdv.species.loc["H2","M"])+(2.*melt_wf["CH4"]/mdv.species.loc["CH4","M"])+(melt_wf["H2S"]/mdv.species.loc["H2S","M"]))
    return wmCO2eq, wmH2Oeq

# calculate weight fraction of species in the gas
def gas_wf(gas_mf,models):
    wg_O2 = (mdv.species.loc["O2","M"]*gas_mf["O2"])/gas_mf["Xg_t"]
    wg_H2 = (mdv.species.loc["H2","M"]*gas_mf["H2"])/gas_mf["Xg_t"]
    wg_H2O = (mdv.species.loc["H2O","M"]*gas_mf["H2O"])/gas_mf["Xg_t"]
    wg_H2S = (mdv.species.loc["H2S","M"]*gas_mf["H2S"])/gas_mf["Xg_t"]
    wg_S2 = (mdv.species.loc["S2","M"]*gas_mf["S2"])/gas_mf["Xg_t"]
    wg_SO2 = (mdv.species.loc["SO2","M"]*gas_mf["SO2"])/gas_mf["Xg_t"]
    wg_CO2 = (mdv.species.loc["CO2","M"]*gas_mf["CO2"])/gas_mf["Xg_t"]
    wg_CO = (mdv.species.loc["CO","M"]*gas_mf["CO"])/gas_mf["Xg_t"]
    wg_CH4 = (mdv.species.loc["CH4","M"]*gas_mf["CH4"])/gas_mf["Xg_t"]
    wg_OCS = (mdv.species.loc["OCS","M"]*gas_mf["OCS"])/gas_mf["Xg_t"]
    species_X = models.loc["species X","option"]
    wg_X = (mdv.species.loc[species_X,"M"]*gas_mf["X"])/gas_mf["Xg_t"]
    result = {"wg_O2":wg_O2, "wg_H2":wg_H2, "wg_H2O":wg_H2O, "wg_H2S":wg_H2S, "wg_S2":wg_S2, "wg_SO2":wg_SO2, "wg_CO2":wg_CO2, "wg_CO":wg_CO, "wg_CH4":wg_CH4, "wg_OCS":wg_OCS, "wg_X":wg_X}
    return result

# calculate weight fraction of species in the gas relative to total system
def gas_wft(gas_mf):
    gaswf = gas_wf(gas_mf,models)
    wgt_O2 = gaswf["wg_O2"]*gas_mf['wt_g']
    wgt_H2 = gaswf["wg_H2"]*gas_mf['wt_g']
    wgt_H2O = gaswf["wg_H2O"]*gas_mf['wt_g']
    wgt_H2S = gaswf["wg_H2S"]*gas_mf['wt_g']
    wgt_S2 = gaswf["wg_S2"]*gas_mf['wt_g']
    wgt_SO2 = gaswf["wg_SO2"]*gas_mf['wt_g']
    wgt_CO2 = gaswf["wg_CO2"]*gas_mf['wt_g']
    wgt_CO = gaswf["wg_CO"]*gas_mf['wt_g']
    wgt_CH4 = gaswf["wg_CH4"]*gas_mf['wt_g']
    wgt_OCS = gaswf["wg_OCS"]*gas_mf['wt_g']
    wgt_X = gaswf["wg_X"]*gas_mf['wt_g']
    result = {"wgt_O2":wgt_O2, "wgt_H2":wgt_H2, "wgt_H2O":wgt_H2O, "wgt_H2S":wgt_H2S, "wgt_S2":wgt_S2, "wgt_SO2":wgt_SO2, "wgt_CO2":wgt_CO2, "wgt_CO":wgt_CO, "wgt_CH4":wgt_CH4, "wgt_OCS":wgt_OCS, "wgt_X":wgt_X}
    return result

def gas_weight(gas_mf,bulk_wf):
    gaswft = gas_wft(gas_mf)
    Wg_O2 = gaswtf["wgt_O2"]*bulk_wf['Wt']
    Wg_H2 = gaswtf["wgt_H2"]*bulk_wf['Wt']
    Wg_H2O = gaswtf["wgt_H2O"]*bulk_wf['Wt']
    Wg_H2S = gaswtf["wgt_H2S"]*bulk_wf['Wt']
    Wg_S2 = gaswtf["wgt_S2"]*bulk_wf['Wt']
    Wg_SO2 = gaswtf["wgt_SO2"]*bulk_wf['Wt']
    Wg_CO2 = gaswtf["wgt_CO2"]*bulk_wf['Wt']
    Wg_CO = gaswtf["wgt_CO"]*bulk_wf['Wt']
    Wg_CH4 = gaswtf["wgt_CH4"]*bulk_wf['Wt']
    Wg_OCS = gaswtf["wgt_OCS"]*bulk_wf['Wt']
    Wg_X = gaswtf["wgt_X"]*bulk_wf['Wt']
    Wg_t = gas_mf['wt_g']*bulk_wf['Wt']
    result = {"Wg_O2":Wg_O2, "Wg_H2":Wg_H2, "Wg_H2O":Wg_H2O, "Wg_H2S":Wg_H2S, "Wg_S2":Wg_S2, "Wg_SO2":Wg_SO2, "Wg_CO2":Wg_CO2, "Wg_CO":Wg_CO, "Wg_CH4":Wg_CH4, "Wg_OCS":Wg_OCS, "Wg_X":Wg_X,"Wg_t":Wg_t}
    return result

def gas_moles(gas_mf,bulk_wf,models):
    gasw = gas_weight(gas_mf,bulk_wf)
    Xg_O2 = gasw["Wg_O2"]/mdv.species.loc["O2","M"]
    Xg_H2O = gasw["Wg_H2O"]/mdv.species.loc["H2O","M"]
    Xg_H2 = gasw["Wg_H2"]/mdv.species.loc["H2","M"]
    Xg_H2S = gasw["Wg_H2S"]/mdv.species.loc["H2S","M"]
    Xg_S2 = gasw["Wg_S2"]/mdv.species.loc["S2","M"]
    Xg_SO2 = gasw["Wg_SO2"]/mdv.species.loc["SO2","M"]
    Xg_CO2 = gasw["Wg_CO2"]/mdv.species.loc["CO2","M"]
    Xg_CO = gasw["Wg_CO"]/mdv.species.loc["CO","M"]
    Xg_CH4 = gasw["Wg_CH4"]/mdv.species.loc["CH4","M"]
    Xg_OCS = gasw["g_OCS"]/mdv.species.loc["OCS","M"]
    species_X = models.loc["species X","option"]
    Xg_X = gasw["Wg_X"]/mdv.species.loc[species_X,"M"]
    Xt_g = Xg_O2 + Xg_H2 + Xg_H2O + Xg_H2S + Xg_S2 + Xg_SO2 + Xg_CO2 + Xg_CO + Xg_CH4 + Xg_OCS + Xg_X
    result = {"Xg_O2":Xg_O2, "Xg_H2":Xg_H2, "Xg_H2O":Xg_H2O, "Xg_H2S":Xg_H2S, "Xg_S2":Xg_S2, "Xg_SO2":Xg_SO2, "Xg_CO2":Xg_CO2, "Xg_CO":Xg_CO, "Xg_CH4":Xg_CH4, "Xg_OCS":Xg_OCS, "Xg_X":Xg_X,"Xg_t":Xg_t}
    return result
        
# calculate weight fraction of elements in the gas
def gas_elements(gas_mf,models):
    gaswf = gas_wf(gas_mf,models)
    wg_O = gaswf["wg_O2"] + mdv.species.loc["O","M"]*((gaswf["wg_H2O"]/mdv.species.loc["H2O","M"]) + (2.*gaswf["wg_SO2"]/mdv.species.loc["SO2","M"]) + (2.*gaswf["wg_CO2"]/mdv.species.loc["CO2","M"]) + (gaswf["wg_CO"]/mdv.species.loc["CO","M"]) + (gaswf["wg_OCS"]/mdv.species.loc["OCS","M"]))
    wg_H = gaswf["wg_H2"] + mdv.species.loc["H","M"]*((2.*gaswf["wg_H2O"]/mdv.species.loc["H2O","M"]) + (2.*gaswf["wg_H2S"]/mdv.species.loc["H2S","M"]) + (4.*gaswf["wg_CH4"]/mdv.species.loc["CH4","M"]))
    wg_S = gaswf["wg_S2"] + mdv.species.loc["S","M"]*((gaswf["wg_H2S"]/mdv.species.loc["H2S","M"]) + (gaswf["wg_SO2"]/mdv.species.loc["SO2","M"]) + (gaswf["wg_OCS"]/mdv.species.loc["OCS","M"]))
    wg_C = mdv.species.loc["C","M"]*((gaswf["wg_CO2"]/mdv.species.loc["CO2","M"]) + (gaswf["wg_CO"]/mdv.species.loc["CO","M"]) + (gaswf["wg_OCS"]/mdv.species.loc["OCS","M"]) + (gaswf["wg_CH4"]/mdv.species.loc["CH4","M"]))
    wg_X = gaswf['wg_X']
    result = {"wg_O":wg_O,"wg_H":wg_H, "wg_S":wg_S, "wg_C":wg_C, "wg_X":wg_X}
    return result

# calculate weight fraction of elements in the melt
def melt_elements(PT,melt_wf,bulk_wf,gas_comp,models):
    wm_C = mdv.species.loc["C","M"]*((melt_wf["CO2"]/mdv.species.loc["CO2","M"] + 
                                        melt_wf["CO"]/mdv.species.loc["CO","M"] + melt_wf["CH4"]/mdv.species.loc["CH4","M"]))
    wm_H = 2.*mdv.species.loc["H","M"]*((melt_wf["H2OT"]/mdv.species.loc["H2O","M"] + melt_wf["H2"]/mdv.species.loc["H2","M"] + melt_wf["H2S"]/mdv.species.loc["H2S","M"] + (2.*melt_wf["CH4"])/mdv.species.loc["CH4","M"]))
    wm_S = mdv.species.loc["S","M"]*(((melt_wf["S2-"]+melt_wf["S6+"])/mdv.species.loc["S","M"] + (melt_wf["H2S"]/mdv.species.loc["H2S","M"])))
    wm_Fe = bulk_wf["Fe"]/(1.-gas_comp["wt_g"])
    Fe32 = overtotal2ratio(melt_wf["Fe3FeT"])
    wm_O = mdv.species.loc["O","M"]*((melt_wf["H2OT"]/mdv.species.loc["H2O","M"]) + ((2.0*melt_wf["CO2"])/mdv.species.loc["CO2","M"]) + (3.0*melt_wf["S6+"]/mdv.species.loc["S","M"]) + (melt_wf["CO"]/mdv.species.loc["CO","M"]) + (wm_Fe/mdv.species.loc["Fe","M"])*((1.5*Fe32+1.0)/(Fe32+1.0)))
    wm_X = melt_wf["XT"]
    result = {"wm_C":wm_C, "wm_H":wm_H, "wm_S":wm_S, "wm_Fe":wm_Fe, "wm_O":wm_O, "wm_X":wm_X}
    return result



#################################################################################################################################
###################################################### VOLUME AND DENSITY #######################################################
#################################################################################################################################

# molar volume of individual gas species (J/mol/bar - *10 for cm3/bar)
def gas_molar_volume(gas_PT,models):
    P = PT['P']
    T = PT['T']+273.15
    Pr = P/mdv.species.loc[gas_species,"Pcr"]
    Tr = T/mdv.species.loc[gas_species,"Tcr"]
    PTr = {"P":Pr,"T":Tr}
    R = 8.3145 # J/mol/K
    if gas_species == "O2" or "CO2" or "OCS" or "CH4" or "CO":
        y = mdv.y_SS(PT,models)
    elif gas_species == "SO2":
        y = mdv.y_SO2(PT,models)
    elif gas_species == "H2S":
        y = mdv.y_H2S(PT,models)
    elif gas_species == "H2":
        y = mdv.y_H2(PT,models)
    elif gas_species == "H2O":
        y = mdv.y_H2O(PT,models)
    Vm = (R*T*y)/P
    PT = {"P":P,"T":T-273.15}
    return Vm

def Vm_O2(PT,models):
    gas_species = "O2"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

def Vm_H2(PT,models):
    gas_species = "H2"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

def Vm_H2O(PT,models):
    gas_species = "H2O"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

def Vm_CO2(PT,models):
    gas_species = "CO2"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

def Vm_CH4(PT,models):
    gas_species = "CH4"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

def Vm_CO(PT,models):
    gas_species = "CO"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

def Vm_S2(PT,models):
    gas_species = "S2"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

def Vm_SO2(PT,models):
    gas_species = "SO2"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

def Vm_H2S(PT,models):
    gas_species = "H2S"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

def Vm_OCS(PT,models):
    gas_species = "OCS"
    Vm = gas_molar_volume(gas_PT,models)
    return Vm

# molar volume of the gas in J/bar/mol
def Vm_gas(gas_mf,PT,models):
    Vm = gas_mf["O2"]*Vm_O2(PT,models) + gas_mf["H2"]*Vm_H2(PT,models) + gas_mf["H2O"]*Vm_H2O(PT,models) + gas_mf["CO2"]*Vm_CO2(PT,models) + gas_mf["CO"]*Vm_CO(PT,models) + gas_mf["CH4"]*Vm_CH4(PT,models) + gas_mf["S2"]*Vm_S2(PT,models) + gas_mf["H2S"]*Vm_H2S(PT,models) + gas_mf["SO2"]*Vm_SO2(PT,models) + gas_mf["OCS"]*Vm_OCS(PT,models)
    return Vm

# volume of the gas in cm3
def gas_volume(PT,gas_mf,bulk_wf,models):
    Xg_O2, Xg_H2, Xg_H2O, Xg_H2S, Xg_S2, Xg_SO2, Xg_CO2, Xg_CO, Xg_CH4, Xg_OCS, Xt_g = gas_moles(gas_mf,bulk_wf)
    volume = Xt_g*Vm_gas(gas_mf,PT,models) # in J/bar
    volume_cm3 = volume*10.
    return volume_cm3

# density of the gas in g/cm3
def gas_density(PT,gas_mf,bulk_wf,models):
    volume = gas_volume(PT,gas_mf,bulk_wf,models) # cm3
    mass = (gas_mf['wt_g']*bulk_wf['Wt']) # g
    density = mass/volume
    return density


# volume of the melt in cm3
def melt_volume(run,PT,melt_wf,bulk_wf,gas_mf,setup):
    density = mdv.melt_density(run,PT,melt_wf,setup)
    mass = bulk_wf['Wt']*(1-gas_mf['wt_g'])
    volume = mass/density
    return volume

# volume of the system (melt + gas)
def system_density(run,PT,melt_wf,gas_mf,bulk_wf,setup,models):
    wt_g_ = gas_mf['wt_g']
    wt_m_ = 1. - wt_g_
    density_m = mdv.melt_density(run,PT,melt_wf,setup)
    if wt_g_ > 0.:
        density_g = gas_density(PT,gas_mf,bulk_wf,models)
    else:
        density_g = 0.
    density_sys = density_m*wt_m_ + density_g*wt_g_
    return density_sys

##################################################################################################################################
######################################################## melt composition ########################################################
##################################################################################################################################

def melt_comp(run,setup):
    oxides = setup.columns.tolist()
    if "SiO2" in oxides:
        SiO2 = setup.loc[run,"SiO2"]
    else:
        SiO2 = 0.
    if "TiO2" in oxides:
        TiO2 = setup.loc[run,"TiO2"]
    else:
        TiO2 = 0.
    if "Al2O3" in oxides:    
        Al2O3 = setup.loc[run,"Al2O3"]
    else:
        Al2O3 = 0.
    if "FeOT" in oxides:
        FeOT = setup.loc[run,"FeOT"]
    else:
        FeOT = float('NaN')
    if "Fe2O3T" in oxides:        
        Fe2O3T = setup.loc[run,"Fe2O3T"]
    else:
        Fe2O3T = float('NaN')
    if "FeO" in oxides:
        FeO = setup.loc[run,"FeO"]
    else:
        FeO = float('NaN')
    if "Fe2O3" in oxides:
        Fe2O3 = setup.loc[run,"Fe2O3"]
    else:
        Fe2O3 = float('NaN')
    if "MgO" in oxides:
        MgO = setup.loc[run,"MgO"]
    else:
        MgO = 0.
    if "MnO" in oxides:
        MnO = setup.loc[run,"MnO"]
    else:
        MnO = 0.
    if "CaO" in oxides:
        CaO = setup.loc[run,"CaO"]
    else:
        CaO = 0.
    if "Na2O" in oxides:
        Na2O = setup.loc[run,"Na2O"]
    else:
        Na2O = 0.
    if "K2O" in oxides:
        K2O = setup.loc[run,"K2O"]
    else:
        K2O = 0.
    if "P2O5" in oxides:
        P2O5 = setup.loc[run,"P2O5"]
    else:
        P2O5 = 0.
    if "logfO2" in oxides:        
        logfO2 = setup.loc[run,"logfO2"]
    else:
        logfO2 = float('NaN')
    if "Fe3FeT" in oxides:        
        Fe3FeT = setup.loc[run,"Fe3FeT"]
    else:
        Fe3FeT = float('NaN')
    if "DNNO" in oxides:    
        DNNO = setup.loc[run,"DNNO"]
    else:
        DNNO = float('NaN')
    if "DFMQ" in oxides:
        DFMQ = setup.loc[run,"DFMQ"]
    else:
        DFMQ = float('NaN')
    if "S6ST" in oxides:
        S6ST = setup.loc[run,"S6ST"]    
    else:
        S6ST = float('NaN')
    melt_wf = {'SiO2':SiO2,'TiO2':TiO2,'Al2O3':Al2O3,'FeOT':FeOT,'Fe2O3T':Fe2O3T,'FeO':FeO,'Fe2O3':Fe2O3,'MgO':MgO,'MnO':MnO,'CaO':CaO,'Na2O':Na2O,'K2O':K2O,'P2O5':P2O5,"logfO2_i":logfO2,"Fe3FeT_i":Fe3FeT,"DNNO":DNNO,"DFMQ":DFMQ,"S6ST_i":S6ST}
    return melt_wf

# normalise melt composition in weight fraction
def melt_normalise_wf(melt_wf,volatiles,Fe_speciation):
    if volatiles == "yes": # assumes all H is H2O, all C is CO2, all S is S, all X is X
        H2O = (melt_wf["HT"]/mdv.species.loc["H2","M"])*mdv.species.loc["H2O","M"]
        CO2 = (melt_wf["CT"]/mdv.species.loc["C","M"])*mdv.species.loc["CO2","M"]
        S = melt_wf["ST"]
        X = melt_wf["XT"]
    elif volatiles == "water": # assumes all H is H2O
        H2O = (melt_wf["HT"]/mdv.species.loc["H2","M"])*mdv.species.loc["H2O","M"]
        CO2, S, X = 0.,0.,0.
    elif volatiles == "no":
        H2O, CO2, S, X = 0.,0.,0.,0.
    volatiles = H2O + CO2 + S + X
    if Fe_speciation == "no":
        Wm_FeOT_ = Wm_FeOT(melt_wf)
        Wm_FeO_ = 0.
        Wm_Fe2O3_ = 0.
    elif Fe_speciation == "yes":
        Wm_FeOT_ = 0.
        Wm_FeO_ = Wm_FeO(melt_wf)
        Wm_Fe2O3_ = Wm_Fe2O3(melt_wf)
    tot = (melt_wf["SiO2"] + melt_wf["TiO2"] + melt_wf["Al2O3"] + Wm_FeOT_ + Wm_FeO_ + Wm_Fe2O3_ + melt_wf["MgO"] + melt_wf["MnO"] + melt_wf["CaO"] + melt_wf["Na2O"] + melt_wf["K2O"] + melt_wf["P2O5"])
    result = {"SiO2":(melt_wf["SiO2"]/tot)*(1.-volatiles)}
    result["TiO2"] = (melt_wf["TiO2"]/tot)*(1.-volatiles)
    result["Al2O3"] = (melt_wf["Al2O3"]/tot)*(1.-volatiles)
    result["FeOT"] = (Wm_FeOT_/tot)*(1.-volatiles)
    result["FeO"] = (Wm_FeO_/tot)*(1.-volatiles)
    result["Fe2O3"] = (Wm_Fe2O3_/tot)*(1.-volatiles)
    result["MgO"] = (melt_wf["MgO"]/tot)*(1.-volatiles)
    result["MnO"] = (melt_wf["MnO"]/tot)*(1.-volatiles)
    result["CaO"] = (melt_wf["CaO"]/tot)*(1.-volatiles)
    result["Na2O"] = (melt_wf["Na2O"]/tot)*(1.-volatiles)
    result["K2O"] = (melt_wf["K2O"]/tot)*(1.-volatiles)
    result["P2O5"] = (melt_wf["P2O5"]/tot)*(1.-volatiles)
    result["H2O"] = H2O
    result["CO2"] = CO2
    result["S"] = S
    result["X"] = X
    return result

# calculate cation proportions
def melt_cation_proportion(melt_wf,volatiles,Fe_speciation):
    melt_comp = melt_normalise_wf(melt_wf,volatiles,Fe_speciation)
    tot = ((mdv.species.loc["SiO2","no_cat"]*melt_comp["SiO2"])/mdv.species.loc["SiO2","M"]) + ((mdv.species.loc["TiO2","no_cat"]*melt_comp["TiO2"])/mdv.species.loc["TiO2","M"]) + ((mdv.species.loc["Al2O3","no_cat"]*melt_comp["Al2O3"])/mdv.species.loc["Al2O3","M"]) + ((mdv.species.loc["FeO","no_cat"]*melt_comp["FeOT"])/mdv.species.loc["FeO","M"]) + ((mdv.species.loc["FeO","no_cat"]*melt_comp["FeO"])/mdv.species.loc["FeO","M"]) + ((mdv.species.loc["Fe2O3","no_cat"]*melt_comp["Fe2O3"])/mdv.species.loc["Fe2O3","M"]) + ((mdv.species.loc["MgO","no_cat"]*melt_comp["MgO"])/mdv.species.loc["MgO","M"]) + ((mdv.species.loc["MnO","no_cat"]*melt_comp["MnO"])/mdv.species.loc["MnO","M"]) + ((mdv.species.loc["CaO","no_cat"]*melt_comp["CaO"])/mdv.species.loc["CaO","M"]) + ((mdv.species.loc["Na2O","no_cat"]*melt_comp["Na2O"])/mdv.species.loc["Na2O","M"]) + ((mdv.species.loc["K2O","no_cat"]*melt_comp["K2O"])/mdv.species.loc["K2O","M"]) + ((mdv.species.loc["P2O5","no_cat"]*melt_comp["P2O5"])/mdv.species.loc["P2O5","M"]) + ((mdv.species.loc["H2O","no_cat"]*melt_comp["H2O"])/mdv.species.loc["H2O","M"]) + ((mdv.species.loc["CO2","no_cat"]*melt_comp["CO2"])/mdv.species.loc["CO2","M"])
    result = {"Si":((mdv.species.loc["SiO2","no_cat"]*melt_comp["SiO2"])/mdv.species.loc["SiO2","M"])/tot}
    result["Ti"] = ((mdv.species.loc["TiO2","no_cat"]*melt_comp["TiO2"])/mdv.species.loc["TiO2","M"])/tot
    result["Al"] = ((mdv.species.loc["Al2O3","no_cat"]*melt_comp["Al2O3"])/mdv.species.loc["Al2O3","M"])/tot
    result["FeT"] = ((mdv.species.loc["FeO","no_cat"]*melt_comp["FeOT"])/mdv.species.loc["FeO","M"])/tot
    result["Fe2"] = ((mdv.species.loc["FeO","no_cat"]*melt_comp["FeO"])/mdv.species.loc["FeO","M"])/tot
    result["Fe3"] = ((mdv.species.loc["Fe2O3","no_cat"]*melt_comp["Fe2O3"])/mdv.species.loc["Fe2O3","M"])/tot    
    result["Mg"] = ((mdv.species.loc["MgO","no_cat"]*melt_comp["MgO"])/mdv.species.loc["MgO","M"])/tot
    result["Mn"] = ((mdv.species.loc["MgO","no_cat"]*melt_comp["MnO"])/mdv.species.loc["MnO","M"])/tot
    result["Ca"] = ((mdv.species.loc["CaO","no_cat"]*melt_comp["CaO"])/mdv.species.loc["CaO","M"])/tot
    result["Na"] = ((mdv.species.loc["Na2O","no_cat"]*melt_comp["Na2O"])/mdv.species.loc["Na2O","M"])/tot
    result["K"] = ((mdv.species.loc["K2O","no_cat"]*melt_comp["K2O"])/mdv.species.loc["K2O","M"])/tot
    result["P"] = ((mdv.species.loc["P2O5","no_cat"]*melt_comp["P2O5"])/mdv.species.loc["P2O5","M"])/tot
    result["H"] = ((mdv.species.loc["H2O","no_cat"]*melt_comp["H2O"])/mdv.species.loc["H2O","M"])/tot
    result["C"] = ((mdv.species.loc["CO2","no_cat"]*melt_comp["CO2"])/mdv.species.loc["CO2","M"])/tot
    return result

# calculate mole fractions in the melt
def melt_mole_fraction(melt_wf,models,volatiles,Fe_speciation):
    melt_comp = melt_normalise_wf(melt_wf,volatiles,Fe_speciation)
    species_X = models.loc["species X","option"]
    mol_tot = (melt_comp["SiO2"]/mdv.species.loc["SiO2","M"]) + (melt_comp["TiO2"]/mdv.species.loc["TiO2","M"]) + (melt_comp["Al2O3"]/mdv.species.loc["Al2O3","M"]) + (melt_comp["FeOT"]/mdv.species.loc["FeO","M"]) + (melt_comp["FeO"]/mdv.species.loc["FeO","M"]) + (melt_comp["Fe2O3"]/mdv.species.loc["Fe2O3","M"]) + (melt_comp["MnO"]/mdv.species.loc["MnO","M"]) + (melt_comp["MgO"]/mdv.species.loc["MgO","M"]) + (melt_comp["CaO"]/mdv.species.loc["CaO","M"]) + (melt_comp["Na2O"]/mdv.species.loc["Na2O","M"]) + (melt_comp["K2O"]/mdv.species.loc["K2O","M"]) + (melt_comp["P2O5"]/mdv.species.loc["P2O5","M"]) + (melt_comp["H2O"]/mdv.species.loc["H2O","M"]) + (melt_comp["CO2"]/mdv.species.loc["CO2","M"]) + (melt_comp["S"]/mdv.species.loc["S","M"]) + (melt_comp["X"]/mdv.species.loc[species_X,"M"])
    result = {"SiO2":(melt_comp["SiO2"]/mdv.species.loc["SiO2","M"])/mol_tot}
    result["TiO2"] = (melt_comp["TiO2"]/mdv.species.loc["TiO2","M"])/mol_tot
    result["Al2O3"] = (melt_comp["Al2O3"]/mdv.species.loc["Al2O3","M"])/mol_tot 
    result["FeOT"] = (melt_comp["FeOT"]/mdv.species.loc["FeO","M"])/mol_tot 
    result["FeO"] = (melt_comp["FeO"]/mdv.species.loc["FeO","M"])/mol_tot
    result["Fe2O3"] = (melt_comp["Fe2O3"]/mdv.species.loc["Fe2O3","M"])/mol_tot
    result["MnO"] = (melt_comp["MnO"]/mdv.species.loc["MnO","M"])/mol_tot
    result["MgO"] = (melt_comp["MgO"]/mdv.species.loc["MgO","M"])/mol_tot
    result["CaO"] = (melt_comp["CaO"]/mdv.species.loc["CaO","M"])/mol_tot
    result["Na2O"] = (melt_comp["Na2O"]/mdv.species.loc["Na2O","M"])/mol_tot
    result["K2O"] = (melt_comp["K2O"]/mdv.species.loc["K2O","M"])/mol_tot
    result["P2O5"] = (melt_comp["P2O5"]/mdv.species.loc["P2O5","M"])/mol_tot
    result["H2O"] = (melt_comp["H2O"]/mdv.species.loc["H2O","M"])/mol_tot
    result["CO2"] = (melt_comp["CO2"]/mdv.species.loc["CO2","M"])/mol_tot
    result["S"] = (melt_comp["S"]/mdv.species.loc["S","M"])/mol_tot
    result["X"] = (melt_comp["X"]/mdv.species.loc[species_X,"M"])/mol_tot
    result["mol_tot"] = mol_tot      
    return result

def melt_single_O(melt_wf,volatiles,Fe_speciation):
    melt_comp = melt_normalise_wf(melt_wf,volatiles,Fe_speciation)
    Xmtot = (melt_comp["SiO2"]/(mdv.species.loc["SiO2","M"]/mdv.species.loc["SiO2","no_O"])) + (melt_comp["TiO2"]/(mdv.species.loc["TiO2","M"]/mdv.species.loc["TiO2","no_O"])) + (melt_comp["Al2O3"]/(mdv.species.loc["Al2O3","M"]/mdv.species.loc["Al2O3","no_O"])) + (melt_comp["MnO"]/(mdv.species.loc["MnO","M"]/mdv.species.loc["MnO","no_O"])) + (melt_comp["MgO"]/(mdv.species.loc["MgO","M"]/mdv.species.loc["MgO","no_O"])) + (melt_comp["CaO"]/(mdv.species.loc["CaO","M"]/mdv.species.loc["CaO","no_O"])) + (melt_comp["Na2O"]/(mdv.species.loc["Na2O","M"]/mdv.species.loc["Na2O","no_O"])) + (melt_comp["K2O"]/(mdv.species.loc["K2O","M"]/mdv.species.loc["K2O","no_O"])) + (melt_comp["P2O5"]/(mdv.species.loc["P2O5","M"]/mdv.species.loc["P2O5","no_O"])) + (melt_comp["FeOT"]/(mdv.species.loc["FeO","M"]/mdv.species.loc["FeO","no_O"])) + (melt_comp["FeO"]/(mdv.species.loc["FeO","M"]/mdv.species.loc["FeO","no_O"])) + (melt_comp["Fe2O3"]/(mdv.species.loc["Fe2O3","M"]/mdv.species.loc["Fe2O3","no_O"])) + (melt_comp["H2O"]/(mdv.species.loc["H2O","M"]/mdv.species.loc["H2O","no_O"])) + (melt_comp["CO2"]/(mdv.species.loc["CO2","M"]/mdv.species.loc["CO2","no_O"]))
    result = {"SiO2": (melt_comp["SiO2"]/(mdv.species.loc["SiO2","M"]/mdv.species.loc["SiO2","no_O"]))/Xmtot}
    result["TiO2"] = (melt_comp["TiO2"]/(mdv.species.loc["TiO2","M"]/mdv.species.loc["TiO2","no_O"]))/Xmtot
    result["Al2O3"] = (melt_comp["Al2O3"]/(mdv.species.loc["Al2O3","M"]/mdv.species.loc["Al2O3","no_O"]))/Xmtot 
    result["FeOT"] = (melt_comp["FeOT"]/(mdv.species.loc["FeO","M"]/mdv.species.loc["FeO","no_O"]))/Xmtot 
    result["FeO"] = (melt_comp["FeO"]/(mdv.species.loc["FeO","M"]/mdv.species.loc["FeO","no_O"]))/Xmtot
    result["Fe2O3"] = (melt_comp["Fe2O3"]/(mdv.species.loc["Fe2O3","M"]/mdv.species.loc["Fe2O3","no_O"]))/Xmtot
    result["MnO"] = (melt_comp["MnO"]/(mdv.species.loc["MnO","M"]/mdv.species.loc["MnO","no_O"]))/Xmtot
    result["MgO"] = (melt_comp["MgO"]/(mdv.species.loc["MgO","M"]/mdv.species.loc["MgO","no_O"]))/Xmtot
    result["CaO"] = (melt_comp["CaO"]/(mdv.species.loc["CaO","M"]/mdv.species.loc["CaO","no_O"]))/Xmtot
    result["P2O5"] = (melt_comp["P2O5"]/(mdv.species.loc["P2O5","M"]/mdv.species.loc["P2O5","no_O"]))/Xmtot
    result["Na2O"] = (melt_comp["Na2O"]/(mdv.species.loc["Na2O","M"]/mdv.species.loc["Na2O","no_O"]))/Xmtot
    result["K2O"] = (melt_comp["K2O"]/(mdv.species.loc["K2O","M"]/mdv.species.loc["K2O","no_O"]))/Xmtot
    result["H2O"] = (melt_comp["H2O"]/mdv.species.loc["H2O","M"])/Xmtot
    result["CO2"] = (melt_comp["CO2"]/mdv.species.loc["CO2","M"])/Xmtot 
    result["Xmtot"] = Xmtot
    return result

def melt_pysulfsat(melt_wf): # output dataframe for pysulfsat
    comp = {'SiO2_Liq': [melt_wf["SiO2"]],
            'TiO2_Liq': [melt_wf["TiO2"]],
            'Al2O3_Liq': [melt_wf["Al2O3"]],
            'FeOt_Liq': [Wm_FeOT(melt_wf)],
            'MnO_Liq': [melt_wf["MnO"]],
            'MgO_Liq': [melt_wf["MgO"]],
            'CaO_Liq': [melt_wf["CaO"]],
            'Na2O_Liq': [melt_wf["Na2O"]],
            'K2O_Liq': [melt_wf["K2O"]],
            'P2O5_Liq': [melt_wf["P2O5"]]}
    comp = pd.DataFrame(comp)    
    return comp
