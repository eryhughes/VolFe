# melt_gas.py

import pandas as pd
import numpy as np
import gmpy2 as gp
import math


import model_dependent_variables as mdv

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

def f_H2O(run,PT,melt_wf,setup,species,models):
    Hspeciation = models.loc["Hspeciation","option"]
    if Hspeciation == "none": # fH2O = xmH2OT^2/CH2O
        value = ((xm_H2OT_so(run,melt_wf,setup,species))**2.0)/mdv.C_H2O(run,PT,melt_wf,setup,species,models)
        return value
    else: # regular or ideal: fH2O = xmH2Omol/CH2O
        value = xm_H2Omol_so(run,PT,melt_wf,setup,species,models)/mdv.C_H2O(run,PT,melt_wf,setup,species,models)
        return value
    
def f_CO2(run,PT,melt_wf,setup,species,models):
    CO3model = models.loc["carbonate","option"]
    wm_CO2 = 100.*melt_wf['CO2']
    if CO3model == "Shishkina14": # wtppm CO2 modified from Shishkina et al. (2014) Chem. Geol. 388:112-129
        f = (wm_CO2*10000.0)/mdv.C_CO3(run,PT,melt_wf,setup,species,models)
    else: # fCO2 = xmCO2/C_CO3
        f = xm_CO2_so(run,melt_wf,setup,species)/mdv.C_CO3(run,PT,melt_wf,setup,species,models)
    return f

def f_S2(run,PT,melt_wf,setup,species,models): # wtppm S2- NOT mole fraction due to parameterisation by O'Neill (2020)
    K = mdv.C_S(run,PT,melt_wf,setup,species,models)/1000000.
    fS2 = ((melt_wf["S2-"]/K)**2.)*mdv.f_O2(run,PT,melt_wf,setup,species,models)
    return fS2
    
def f_H2(run,PT,melt_wf,setup,species,models):
    K = mdv.KHOg(PT,models)
    return f_H2O(run,PT,melt_wf,setup,species,models)/(K*mdv.f_O2(run,PT,melt_wf,setup,species,models)**0.5)

def f_CO(run,PT,melt_wf,setup,species,models):
    K = mdv.KCOg(PT,models)
    return f_CO2(run,PT,melt_wf,setup,species,models)/(K*mdv.f_O2(run,PT,melt_wf,setup,species,models)**0.5)

def f_H2S(run,PT,melt_wf,setup,species,models):
    K = mdv.KHOSg(PT,models)
    return (K*f_S2(run,PT,melt_wf,setup,species,models)**0.5*f_H2O(run,PT,melt_wf,setup,species,models))/mdv.f_O2(run,PT,melt_wf,setup,species,models)**0.5

def f_SO2(run,PT,melt_wf,setup,species,models):
    K = mdv.KOSg(PT,models)
    return K*mdv.f_O2(run,PT,melt_wf,setup,species,models)*f_S2(run,PT,melt_wf,setup,species,models)**0.5

def f_SO3(run,PT,melt_wf,setup,species,models):
    K = mdv.KOSg2(PT,models)
    return K*(mdv.f_O2(run,PT,melt_wf,setup,species,models))**1.5*(f_S2(run,PT,melt_wf,setup,species,models))**0.5

def f_CH4(run,PT,melt_wf,setup,species,models):
    K = mdv.KCOHg(PT,models)
    return (f_CO2(run,PT,melt_wf,setup,species,models)*f_H2O(run,PT,melt_wf,setup,species,models)**2.0)/(K*mdv.f_O2(run,PT,melt_wf,setup,species,models)**2.0)

def f_OCS(run,PT,melt_wf,setup,species,models):
    OCSmodel = models.loc["carbonylsulphide","option"]
    K = mdv.KOCSg(PT,models)
    if OCSmodel == "COHS":
        if f_H2O(run,PT,melt_wf,setup,species,models) > 0.0:
            return (f_CO2(run,PT,melt_wf,setup,species,models)*f_H2S(run,PT,melt_wf,setup,species,models))/(f_H2O(run,PT,melt_wf,setup,species,models)*K)
        else:
            return 0.0
    else:
        if f_CO2(run,PT,melt_wf,setup,species,models) > 0.0:
            return ((f_CO(run,PT,melt_wf,setup,species,models)**3.0)*f_SO2(run,PT,melt_wf,setup,species,models))/((f_CO2(run,PT,melt_wf,setup,species,models)**2.0)*K)
        else:
            return 0.0
        
def f_X(run,PT,melt_wf,setup,species,models):
    K = mdv.C_X(run,PT,melt_wf,setup,species,models)/1000000.
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
        
def S6S2_2_fO2(S62,melt_wf,run,PT,setup,species,models):
    CSO4 = mdv.C_SO4(run,PT,melt_wf,setup,species,models)
    CS = mdv.C_S(run,PT,melt_wf,setup,species,models)
    if models.loc["H2S_m","option"] == "no":
        fO2 = ((S62*CS)/CSO4)**0.5
    elif models.loc["H2S_m","option"] == "yes":
        KHS = mdv.KHOSg(PT,models)
        CH2S = mdv.C_H2S(run,PT,melt_wf,setup,species,models)
        CH2OT = mdv.C_H2O(run,PT,melt_wf,setup,species,models)
        xmH2O = xm_H2OT_so(run,melt_wf,setup,species)
        W = (CSO4/((KHS*CH2S*(xmH2O**2./CH2OT)) + CS))
        fO2 = (S62/W)**0.5
    return fO2
    
    
########################        
### partial pressure ###
########################

def p_H2(run,PT,melt_wf,setup,species,models):
    return f_H2(run,PT,melt_wf,setup,species,models)/mdv.y_H2(PT,species,models)
def p_H2O(run,PT,melt_wf,setup,species,models):
    return f_H2O(run,PT,melt_wf,setup,species,models)/mdv.y_H2O(PT,species,models)
def p_O2(run,PT,melt_wf,setup,species,models):
    return mdv.f_O2(run,PT,melt_wf,setup,species,models)/mdv.y_O2(PT,species,models)
def p_SO2(run,PT,melt_wf,setup,species,models):
    return f_SO2(run,PT,melt_wf,setup,species,models)/mdv.y_SO2(PT,species,models)
def p_S2(run,PT,melt_wf,setup,species,models):
    return f_S2(run,PT,melt_wf,setup,species,models)/mdv.y_S2(PT,species,models)
def p_H2S(run,PT,melt_wf,setup,species,models):
    return f_H2S(run,PT,melt_wf,setup,species,models)/mdv.y_H2S(PT,species,models)
def p_CO2(run,PT,melt_wf,setup,species,models):
    return f_CO2(run,PT,melt_wf,setup,species,models)/mdv.y_CO2(PT,species,models)
def p_CO(run,PT,melt_wf,setup,species,models):
    return f_CO(run,PT,melt_wf,setup,species,models)/mdv.y_CO(PT,species,models)
def p_CH4(run,PT,melt_wf,setup,species,models):
    return f_CH4(run,PT,melt_wf,setup,species,models)/mdv.y_CH4(PT,species,models)
def p_OCS(run,PT,melt_wf,setup,species,models):
    return f_OCS(run,PT,melt_wf,setup,species,models)/mdv.y_OCS(PT,species,models)
def p_X(run,PT,melt_wf,setup,species,models):
    return f_X(run,PT,melt_wf,setup,species,models)/mdv.y_X(PT,species,models)

def p_tot(run,PT,melt_wf,setup,species,models):
    return p_H2(run,PT,melt_wf,setup,species,models) + p_H2O(run,PT,melt_wf,setup,species,models) + p_O2(run,PT,melt_wf,setup,species,models) + p_SO2(run,PT,melt_wf,setup,species,models) + p_S2(run,PT,melt_wf,setup,species,models) + p_H2S(run,PT,melt_wf,setup,species,models) + p_CO2(run,PT,melt_wf,setup,species,models) + p_CO(run,PT,melt_wf,setup,species,models) + p_CH4(run,PT,melt_wf,setup,species,models) + p_OCS(run,PT,melt_wf,setup,species,models) + p_X(run,PT,melt_wf,setup,species,models)


############################       
### vapor molar fraction ###
############################

def xg_H2(run,PT,melt_wf,setup,species,models):
    return p_H2(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_H2O(run,PT,melt_wf,setup,species,models):
    return p_H2O(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_O2(run,PT,melt_wf,setup,species,models):
    return p_O2(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_SO2(run,PT,melt_wf,setup,species,models):
    return p_SO2(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_S2(run,PT,melt_wf,setup,species,models):
    return p_S2(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_H2S(run,PT,melt_wf,setup,species,models):
    return p_H2S(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_CO2(run,PT,melt_wf,setup,species,models):
    return p_CO2(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_CO(run,PT,melt_wf,setup,species,models):
    return p_CO(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_CH4(run,PT,melt_wf,setup,species,models):
    return p_CH4(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_OCS(run,PT,melt_wf,setup,species,models):
    return p_OCS(run,PT,melt_wf,setup,species,models)/PT['P']
def xg_X(run,PT,melt_wf,setup,species,models):
    return p_X(run,PT,melt_wf,setup,species,models)/PT['P']
def Xg_tot(run,PT,melt_wf,setup,species,models):
    species_X = models.loc["species X","option"]
    Xg_t = xg_CO2(run,PT,melt_wf,setup,species,models)*species.loc["CO2","M"] + xg_CO(run,PT,melt_wf,setup,species,models)*species.loc["CO","M"] + xg_O2(run,PT,melt_wf,setup,species,models)*species.loc["O2","M"] + xg_H2O(run,PT,melt_wf,setup,species,models)*species.loc["H2O","M"] + xg_H2(run,PT,melt_wf,setup,species,models)*species.loc["H2","M"] + xg_CH4(run,PT,melt_wf,setup,species,models)*species.loc["CH4","M"] + xg_SO2(run,PT,melt_wf,setup,species,models)*species.loc["SO2","M"] + xg_S2(run,PT,melt_wf,setup,species,models)*species.loc["S2","M"] + xg_H2S(run,PT,melt_wf,setup,species,models)*species.loc["H2S","M"] + xg_OCS(run,PT,melt_wf,setup,species,models)*species.loc["OCS","M"] + xg_X(run,PT,melt_wf,setup,species,models)*species.loc[species_X,"M"]
    return Xg_t


##########################
### melt mole fraction ###
##########################

# totals
def wm_vol(melt_wf): # wt% total volatiles in the melt
    wm_H2OT = 100.*melt_wf["H2OT"]
    wm_CO2 = 100.*melt_wf["CO2"]
    return wm_H2OT + wm_CO2 #+ wm_S(wm_ST) + wm_SO3(wm_ST,species)
def wm_nvol(melt_wf): # wt% total of non-volatiles in the melt
    return 100.0 - wm_vol(melt_wf)

# molecular mass on a singular oxygen basis
def M_m_SO(run,melt_wf,setup,species):
    # single oxide, no volatiles, all Fe as FeOT
    tot, SiO2, TiO2, Al2O3, FeOT, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, H2O, CO2 = melt_single_O(run,melt_wf,setup,species,"no","no")
    M_m = 1./tot
    return M_m

# molecular mass on a oxide basis
def M_m_ox(run,melt_wf,setup,species,models): # no volatiles, all Fe as FeOT
    mol_tot, SiO2, TiO2, Al2O3, FeOT, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, H2O, CO2, S, X   = melt_mole_fraction(run,melt_wf,setup,species,models,"no","no") 
    M_m = 1./mol_tot
    return M_m
    
# Number of moles in the melt
def Xm_H2OT(melt_wf,species):
    wm_H2OT = 100.*melt_wf['H2OT']
    return wm_H2OT/species.loc["H2O","M"]
def Xm_CO2(melt_wf,species):
    wm_CO2 = 100.*melt_wf['CO2']
    return wm_CO2/species.loc["CO2","M"]

# Mole fraction in the melt based on mixing between volatile-free melt on a singular oxygen basis and volatiles
def Xm_m_so(run,melt_wf,setup,species): # singular oxygen basis
    return wm_nvol(melt_wf)/M_m_SO(run,melt_wf,setup,species)    
def Xm_tot_so(run,melt_wf,setup,species):
    return Xm_H2OT(melt_wf,species) + Xm_CO2(melt_wf,species) + Xm_m_so(run,melt_wf,setup,species)
def xm_H2OT_so(run,melt_wf,setup,species):
    return Xm_H2OT(melt_wf,species)/Xm_tot_so(run,melt_wf,setup,species)
def xm_CO2_so(run,melt_wf,setup,species):
    return Xm_CO2(melt_wf,species)/Xm_tot_so(run,melt_wf,setup,species)
def xm_melt_so(run,melt_wf,setup,species):
    return Xm_m_so(run,melt_wf,setup,species)/Xm_tot_so(run,melt_wf,setup,species)
def Xm_t_so(run,melt_wf,setup,species):
    return xm_H2OT_so(run,melt_wf,setup,species)*species.loc["H2O","M"] + xm_CO2_so(run,melt_wf,setup,species)*species.loc["CO2","M"] + xm_melt_so(run,melt_wf,setup,species)*M_m_SO(run,melt_wf,setup,species)

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

def xm_OH_so(run,PT,melt_wf,setup,species,models):
    Hspeciation = models.loc["Hspeciation","option"]
    Hspeccomp = models.loc["Hspeccomp","option"]
    T_K = PT['T']+273.15
    wm_H2OT = 100.*melt_wf['H2OT']
    
    Z = xm_H2OT_so(run,melt_wf,setup,species)
    
    def f(A, B, C, Z, x): # regular solution model rearranged to equal 0 to solve for xm_OH
        return (A + B*x + C*(Z - 0.5*x) + gp.log((x**2.0)/((Z - 0.5*x)*(1.0 - Z - 0.5*x))))

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
    elif Hspeciation == "ideal": # ideal mixture from Silver & Stolper (1985) J. Geol. 93(2) 161-177
        K = mdv.KHOm(run,PT,melt_wf,setup,species,models)
        return (0.5 - (0.25 - ((K - 4.0)/K)*(Z - Z**2.0))**0.5)/((K - 4.0)/(2.0*K)) 
    elif Hspeciation == "none": # all OH-
        return 0.0
    
    elif Hspeciation == "regular": # else use regular solution model from Silver & Stolper (1989) J.Pet. 30(3)667-709
        tol = 0.000001 #tolerance
        K = mdv.KHOm(run,PT,melt_wf,setup,species,models)
        x0 = (0.5 - (0.25 - ((K - 4.0)/K)*(Z - Z**2.0))**0.5)/((K - 4.0)/(2.0*K))
        R = 83.15 #cm3.bar/mol.K
        A, B, C = mdv.KregH2O(run,PT,melt_wf,setup,species,models)
        return nr(A,B,C,x0,tol)

def xm_H2Omol_so(run,PT,melt_wf,setup,species,models):
    Z = xm_H2OT_so(run,melt_wf,setup,species)
    return Z - 0.5*xm_OH_so(run,PT,melt_wf,setup,species,models)

def wm_H2Omol_OH(run,PT,melt_wf,setup,species,models):
    H2Omol = xm_H2Omol_so(run,PT,melt_wf,setup,species,models)*species.loc["H2O","M"]
    OH = xm_OH_so(run,PT,melt_wf,setup,species,models)*species.loc["OH","M"]
    melt = xm_melt_so(run,melt_wf,setup,species)*M_m_SO(run,melt_wf,setup,species)
    wm_H2Omol = 100.*H2Omol/(H2Omol + OH + melt)
    wm_OH = 100.*OH/(H2Omol + OH + melt)
    return wm_H2Omol, wm_OH # wt% 


##########################
### sulphur speciation ###
##########################

def S6S2(run,PT,melt_wf,setup,species,models):
    T_K = PT['T']+273.15
    model = models.loc["sulphate","option"]
    if model == "Nash19":
        A, B = mdv.S_Nash19_terms(PT)
        result = 10.**(A*math.log10(Fe3Fe2(melt_wf)) + B)
    else:
        CSO4 = mdv.C_SO4(run,PT,melt_wf,setup,species,models)
        CS = mdv.C_S(run,PT,melt_wf,setup,species,models)
        fO2 = mdv.f_O2(run,PT,melt_wf,setup,species,models)
        if models.loc["H2S_m","option"] == "no":
            result = (CSO4/CS)*fO2**2.
        elif models.loc["H2S_m","option"] == "yes":
            KHS = mdv.KHOSg(PT,models)
            CH2S = mdv.C_H2S(run,PT,melt_wf,setup,species,models)
            CH2OT = mdv.C_H2O(run,PT,melt_wf,setup,species,models)
            xmH2O = xm_H2OT_so(run,melt_wf,setup,species)
            result = (CSO4/((KHS*CH2S*(xmH2O**2./CH2OT)) + CS))*fO2**2.
    return result
    
def S6ST(run,PT,melt_wf,setup,species,models):
    S6S2_ = S6S2(run,PT,melt_wf,setup,species,models)
    return S6S2_/(S6S2_+1.0)

def wm_S(run,PT,melt_wf,setup,species,models):
    wm_ST = 100.*melt_wf['ST']
    S6ST_ = S6ST(run,PT,melt_wf,setup,species,models)
    return wm_ST*(1.0-S6ST_)

def wm_SO3(run,PT,melt_wf,setup,species,models):
    wm_ST = 100.*melt_wf['ST']
    S6ST_ = S6ST(run,PT,melt_wf,setup,species,models)    
    return ((wm_ST*S6ST_)/species.loc["S","M"])*species.loc["SO3","M"]


#######################
### iron speciation ###
#######################

def Fe3Fe2(melt_wf):
    Fe3FeT = melt_wf['Fe3FeT']
    return Fe3FeT/(1.0 - Fe3FeT)

def Wm_FeT(run,setup,species):
    if setup.loc[run,"FeOT"] > 0.0:
        return (setup.loc[run,"FeOT"]/species.loc["FeO","M"])*species.loc["Fe","M"]
    elif setup.loc[run,"Fe2O3T"] > 0.0:
        return (setup.loc[run,"Fe2O3T"]/species.loc["Fe2O3","M"])*species.loc["Fe","M"]
    else:
        return ((setup.loc[run,"FeO"]/species.loc["FeO","M"]) + (setup.loc[run,"Fe2O3"]/species.loc["Fe2O3","M"]))*species.loc["Fe","M"]

def Wm_FeO(run,melt_wf,setup,species):
    Fe3FeT = melt_wf['Fe3FeT']
    return (Wm_FeT(run,setup,species)/species.loc["Fe","M"])*(1.0-Fe3FeT)*species.loc["FeO","M"]

def Wm_Fe2O3(run,melt_wf,setup,species):
    Fe3FeT = melt_wf['Fe3FeT']
    return (Wm_FeT(run,setup,species)/species.loc["Fe","M"])*Fe3FeT*species.loc["Fe2O3","M"]

def Wm_FeOT(run,setup,species):
    return (Wm_FeT(run,setup,species)/species.loc["Fe","M"])*species.loc["FeO","M"]

def wm_Fe_nv(run,melt_wf,setup,species): # no volatiles
    Wm_tot = setup.loc[run,"SiO2"] + setup.loc[run,"TiO2"] + setup.loc[run,"Al2O3"] + setup.loc[run,"MnO"] + setup.loc[run,"MgO"] + setup.loc[run,"MnO"] + setup.loc[run,"CaO"] + setup.loc[run,"Na2O"] + setup.loc[run,"K2O"] + setup.loc[run,"P2O5"] + Wm_FeO(run,melt_wf,setup,species) + Wm_Fe2O3(run,melt_wf,setup,species)
    FeT = species.loc["Fe","M"]*((2.0*Wm_Fe2O3(run,melt_wf,setup,species)/species.loc["Fe2O3","M"]) + (Wm_FeO(run,melt_wf,setup,species)/species.loc["FeO","M"]))
    return 100.0*FeT/Wm_tot

def Fe3FeT_i(run,PT,melt_wf,setup,species,models):
    model = models.loc["fO2","option"]
    T_K = PT['T']+273.15
    
    if model == "buffered":
        fO2 = 10**(setup.loc[run,"logfO2"])
        return mdv.fO22Fe3FeT(fO2,run,PT,melt_wf,setup,species,models)
    else:
        if pd.isnull(setup.loc[run,"Fe3FeT"]) == False:
            return setup.loc[run,"Fe3FeT"]
        elif pd.isnull(setup.loc[run,"logfO2"]) == False:
            fO2 = 10.0**(setup.loc[run,"logfO2"])
            return mdv.fO22Fe3FeT(fO2,run,PT,melt_wf,setup,species,models)
        elif pd.isnull(setup.loc[run,"DNNO"]) == False:
            D = setup.loc[run,"DNNO"]
            fO2 = Dbuffer2fO2(PT,D,"NNO",models)
            return mdv.fO22Fe3FeT(fO2,run,PT,melt_wf,setup,species,models)
        elif pd.isnull(setup.loc[run,"DFMQ"]) == False:
            D = setup.loc[run,"DFMQ"]
            fO2 = Dbuffer2fO2(PT,D,"FMQ",models)
            return mdv.fO22Fe3FeT(fO2,run,PT,melt_wf,setup,species,models)
        else:
            return ((2.0*setup.loc[run,"Fe2O3"])/species.loc["Fe2O3","M"])/(((2.0*setup.loc[run,"Fe2O3"])/species.loc["Fe2O3","M"]) + (setup.loc[run,"FeO"]/species.loc["FeO","M"]))
        
        

    
    
#################################################################################################################################
############################################## converting gas and melt compositions #############################################
#################################################################################################################################

# calculate weight fraction of species in the gas
def gas_wf(gas_mf,species):
    wg_O2 = (species.loc["O2","M"]*gas_mf["O2"])/gas_mf["Xg_t"]
    wg_H2 = (species.loc["H2","M"]*gas_mf["H2"])/gas_mf["Xg_t"]
    wg_H2O = (species.loc["H2O","M"]*gas_mf["H2O"])/gas_mf["Xg_t"]
    wg_H2S = (species.loc["H2S","M"]*gas_mf["H2S"])/gas_mf["Xg_t"]
    wg_S2 = (species.loc["S2","M"]*gas_mf["S2"])/gas_mf["Xg_t"]
    wg_SO2 = (species.loc["SO2","M"]*gas_mf["SO2"])/gas_mf["Xg_t"]
    wg_CO2 = (species.loc["CO2","M"]*gas_mf["CO2"])/gas_mf["Xg_t"]
    wg_CO = (species.loc["CO","M"]*gas_mf["CO"])/gas_mf["Xg_t"]
    wg_CH4 = (species.loc["CH4","M"]*gas_mf["CH4"])/gas_mf["Xg_t"]
    wg_OCS = (species.loc["OCS","M"]*gas_mf["OCS"])/gas_mf["Xg_t"]
    species_X = models.loc["species X","option"]
    wg_X = (species.loc[species_X,"M"]*gas_mf["X"])/gas_mf["Xg_t"]
    return wg_O2, wg_H2, wg_H2O, wg_H2S, wg_S2, wg_SO2, wg_CO2, wg_CO, wg_CH4, wg_OCS, wg_X

# calculate weight fraction of species in the gas relative to total system
def gas_wft(gas_mf,species):
    wg_O2, wg_H2, wg_H2O, wg_H2S, wg_S2, wg_SO2, wg_CO2, wg_CO, wg_CH4, wg_OCS, wg_X = gas_wf(gas_mf,species)
    wt_g = gas_mf['wt_g']
    wgt_O2 = wg_O2*wt_g
    wgt_H2 = wg_H2*wt_g
    wgt_H2O = wg_H2O*wt_g
    wgt_H2S = wg_H2S*wt_g
    wgt_S2 = wg_S2*wt_g
    wgt_SO2 = wg_SO2*wt_g
    wgt_CO2 = wg_CO2*wt_g
    wgt_CO = wg_CO*wt_g
    wgt_CH4 = wg_CH4*wt_g
    wgt_OCS = wg_OCS*wt_g
    wgt_X = wg_X*wt_g
    return wgt_O2, wgt_H2, wgt_H2O, wgt_H2S, wgt_S2, wgt_SO2, wgt_CO2, wgt_CO, wgt_CH4, wgt_OCS, wgt_X

def gas_weight(gas_mf,bulk_wf,species):
    wgt_O2, wgt_H2, wgt_H2O, wgt_H2S, wgt_S2, wgt_SO2, wgt_CO2, wgt_CO, wgt_CH4, wgt_OCS, wgt_X = gas_wft(gas_mf,species)
    Wg_O2 = wgt_O2*bulk_wf['Wt']
    Wg_H2 = wgt_H2*bulk_wf['Wt']
    Wg_H2O = wgt_H2O*bulk_wf['Wt']
    Wg_H2S = wgt_H2S*bulk_wf['Wt']
    Wg_S2 = wgt_S2*bulk_wf['Wt']
    Wg_SO2 = wgt_SO2*bulk_wf['Wt']
    Wg_CO2 = wgt_CO2*bulk_wf['Wt']
    Wg_CO = wgt_CO*bulk_wf['Wt']
    Wg_CH4 = wgt_CH4*bulk_wf['Wt']
    Wg_OCS = wgt_OCS*bulk_wf['Wt']
    Wg_X = wgt_X*bulk_wf['Wt']
    Wg_t = gas_mf['wt_g']*bulk_wf['Wt']
    return Wg_O2, Wg_H2, Wg_H2O, Wg_H2S, Wg_S2, Wg_SO2, Wg_CO2, Wg_CO, Wg_CH4, Wg_OCS, Wg_X, Wg_t 

def gas_moles(gas_mf,bulk_wf,species):
    Wg_O2, Wg_H2, Wg_H2O, Wg_H2S, Wg_S2, Wg_SO2, Wg_CO2, Wg_CO, Wg_CH4, Wg_OCS, Wg_X, Wg_t = gas_weight(gas_mf,bulk_wf,species)
    Xg_O2 = Wg_O2/species.loc["O2","M"]
    Xg_H2O = Wg_H2O/species.loc["H2O","M"]
    Xg_H2 = Wg_H2/species.loc["H2","M"]
    Xg_H2S = Wg_H2S/species.loc["H2S","M"]
    Xg_S2 = Wg_S2/species.loc["S2","M"]
    Xg_SO2 = Wg_SO2/species.loc["SO2","M"]
    Xg_CO2 = Wg_CO2/species.loc["CO2","M"]
    Xg_CO = Wg_CO/species.loc["CO","M"]
    Xg_CH4 = Wg_CH4/species.loc["CH4","M"]
    Xg_OCS = Wg_OCS/species.loc["OCS","M"]
    species_X = models.loc["species X","option"]
    Xg_X = Wg_X/species.loc[species_X,"M"]
    Xt_g = Xg_O2 + Xg_H2 + Xg_H2O + Xg_H2S + Xg_S2 + Xg_SO2 + Xg_CO2 + Xg_CO + Xg_CH4 + Xg_OCS + Xg_X
    return Xg_O2, Xg_H2, Xg_H2O, Xg_H2S, Xg_S2, Xg_SO2, Xg_CO2, Xg_CO, Xg_CH4, Xg_OCS, Xg_X, Xt_g
        
# calculate weight fraction of elements in the gas
def gas_elements(gas_mf,species):
    wg_O2, wg_H2, wg_H2O, wg_H2S, wg_S2, wg_SO2, wg_CO2, wg_CO, wg_CH4, wg_OCS, wg_X = gas_wf(gas_mf,species)
    wg_O = wg_O2 + species.loc["O","M"]*((wg_H2O/species.loc["H2O","M"]) + (2.*wg_SO2/species.loc["SO2","M"]) + (2.*wg_CO2/species.loc["CO2","M"]) + (wg_CO/species.loc["CO","M"]) + (wg_OCS/species.loc["OCS","M"]))
    wg_H = wg_H2 + species.loc["H","M"]*((2.*wg_H2O/species.loc["H2O","M"]) + (2.*wg_H2S/species.loc["H2S","M"]) + (4.*wg_CH4/species.loc["CH4","M"]))
    wg_S = wg_S2 + species.loc["S","M"]*((wg_H2S/species.loc["H2S","M"]) + (wg_SO2/species.loc["SO2","M"]) + (wg_OCS/species.loc["OCS","M"]))
    wg_C = species.loc["S","M"]*((wg_CO2/species.loc["CO2","M"]) + (wg_CO/species.loc["CO","M"]) + (wg_OCS/species.loc["OCS","M"]) + (wg_CH4/species.loc["CH4","M"]))
    return wg_O, wg_H, wg_S, wg_C, wg_X

# calculate weight fraction of elements in the melt
def melt_elements(run,PT,melt_wf,bulk_wf,gas_comp,setup,species,models):
    wm_C = species.loc["C","M"]*((melt_wf["CO2"]/species.loc["CO2","M"] + 
                                        melt_wf["CO"]/species.loc["CO","M"] + melt_wf["CH4"]/species.loc["CH4","M"]))
    wm_H = 2.*species.loc["H","M"]*((melt_wf["H2OT"]/species.loc["H2O","M"] + melt_wf["H2"]/species.loc["H2","M"] + melt_wf["H2S"]/species.loc["H2S","M"] + (2.*melt_wf["CH4"])/species.loc["CH4","M"]))
    wm_S = species.loc["S","M"]*(((melt_wf["S2-"]+melt_wf["S6+"])/species.loc["S","M"] + (melt_wf["H2S"]/species.loc["H2S","M"])))
    wm_Fe = bulk_wf["Fe"]/(1.-gas_comp["wt_g"])
    Fe32 = overtotal2ratio(melt_wf["Fe3FeT"])
    wm_O = species.loc["O","M"]*((melt_wf["H2OT"]/species.loc["H2O","M"]) + ((2.0*melt_wf["CO2"])/species.loc["CO2","M"]) + (3.0*melt_wf["S6+"]/species.loc["S","M"]) + (melt_wf["CO"]/species.loc["CO","M"]) + (wm_Fe/species.loc["Fe","M"])*((1.5*Fe32+1.0)/(Fe32+1.0)))
    wm_X = melt_wf["X"]
    return wm_C, wm_H, wm_S, wm_Fe, wm_O, wm_X



#################################################################################################################################
###################################################### VOLUME AND DENSITY #######################################################
#################################################################################################################################

# molar volume of individual gas species (J/mol/bar - *10 for cm3/bar)
def gas_molar_volume(gas_species,PT,species,models):
    P = PT['P']
    T = PT['T']+273.15
    Pr = P/species.loc[gas_species,"Pcr"]
    Tr = T/species.loc[gas_species,"Tcr"]
    PTr = {"P":Pr,"T":Tr}
    R = 8.3145 # J/mol/K
    if gas_species == "O2" or "CO2" or "OCS" or "CH4" or "CO":
        y = mdv.y_SS(PT,species,models)
    elif gas_species == "SO2":
        y = mdv.y_SO2(PT,species,models)
    elif gas_species == "H2S":
        y = mdv.y_H2S(PT,species,models)
    elif gas_species == "H2":
        y = mdv.y_H2(PT,species,models)
    elif gas_species == "H2O":
        y = mdv.y_H2O(PT,species,models)
    Vm = (R*T*y)/P
    PT = {"P":P,"T":T-273.15}
    return Vm

def Vm_O2(PT,species,models):
    gas_species = "O2"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

def Vm_H2(PT,species,models):
    gas_species = "H2"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

def Vm_H2O(PT,species,models):
    gas_species = "H2O"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

def Vm_CO2(PT,species,models):
    gas_species = "CO2"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

def Vm_CH4(PT,species,models):
    gas_species = "CH4"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

def Vm_CO(PT,species,models):
    gas_species = "CO"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

def Vm_S2(PT,species,models):
    gas_species = "S2"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

def Vm_SO2(PT,species,models):
    gas_species = "SO2"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

def Vm_H2S(PT,species,models):
    gas_species = "H2S"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

def Vm_OCS(PT,species,models):
    gas_species = "OCS"
    Vm = gas_molar_volume(gas_species,PT,species,models)
    return Vm

# molar volume of the gas in J/bar/mol
def Vm_gas(gas_mf,PT,species,models):
    Vm = gas_mf["O2"]*Vm_O2(PT,species,models) + gas_mf["H2"]*Vm_H2(PT,species,models) + gas_mf["H2O"]*Vm_H2O(PT,species,models) + gas_mf["CO2"]*Vm_CO2(PT,species,models) + gas_mf["CO"]*Vm_CO(PT,species,models) + gas_mf["CH4"]*Vm_CH4(PT,species,models) + gas_mf["S2"]*Vm_S2(PT,species,models) + gas_mf["H2S"]*Vm_H2S(PT,species,models) + gas_mf["SO2"]*Vm_SO2(PT,species,models) + gas_mf["OCS"]*Vm_OCS(PT,species,models)
    return Vm

# volume of the gas in cm3
def gas_volume(PT,gas_mf,bulk_wf,species,models):
    Xg_O2, Xg_H2, Xg_H2O, Xg_H2S, Xg_S2, Xg_SO2, Xg_CO2, Xg_CO, Xg_CH4, Xg_OCS, Xt_g = gas_moles(gas_mf,bulk_wf,species)
    volume = Xt_g*Vm_gas(gas_mf,PT,species,models) # in J/bar
    volume_cm3 = volume*10.
    return volume_cm3

# density of the gas in g/cm3
def gas_density(PT,gas_mf,bulk_wf,species,models):
    volume = gas_volume(PT,gas_mf,bulk_wf,species,models) # cm3
    mass = (gas_mf['wt_g']*bulk_wf['Wt']) # g
    density = mass/volume
    return density


# volume of the melt in cm3
def melt_volume(run,PT,melt_wf,bulk_wf,gas_mf,setup,species):
    density = mdv.melt_density(run,PT,melt_wf,setup,species)
    mass = bulk_wf['Wt']*(1-gas_mf['wt_g'])
    volume = mass/density
    return volume

# volume of the system (melt + gas) - could include sulphide and sulphate phases at some point?!
def system_density(run,PT,melt_wf,gas_mf,bulk_wf,setup,species,models):
    wt_g_ = gas_mf['wt_g']
    wt_m_ = 1. - wt_g_
    density_m = mdv.melt_density(run,PT,melt_wf,setup,species)
    if wt_g_ > 0.:
        density_g = gas_density(PT,gas_mf,bulk_wf,species,models)
    else:
        density_g = 0.
    density_sys = density_m*wt_m_ + density_g*wt_g_
    return density_sys

##################################################################################################################################
######################################################## melt composition ########################################################
##################################################################################################################################

# normalise melt composition in weight fraction
def melt_normalise_wf(run,melt_wf,setup,species,volatiles,Fe_speciation):
    if volatiles == "yes": # assumes all H is H2O, all C is CO2, all S is S, all X is X
        
        H2O = (melt_wf["HT"]/species.loc["H2","M"])*species.loc["H2O","M"]
        CO2 = (melt_wf["CT"]/species.loc["C","M"])*species.loc["CO2","M"]
        S = melt_wf["ST"]
        X = melt_wf["XT"]
    elif volatiles == "water": # assumes all H is H2O
        H2O = (melt_wf["HT"]/species.loc["H2","M"])*species.loc["H2O","M"]
        CO2, S, X = 0.,0.,0.
    elif volatiles == "no":
        H2O, CO2, S, X = 0.,0.,0.,0.
    volatiles = H2O + CO2 + S + X
    if Fe_speciation == "no":
        Wm_FeOT_ = Wm_FeOT(run,setup,species)
        Wm_FeO_ = 0.
        Wm_Fe2O3_ = 0.
    elif Fe_speciation == "yes":
        Wm_FeOT_ = 0.
        Wm_FeO_ = Wm_FeO(run,melt_wf,setup,species)
        Wm_Fe2O3_ = Wm_Fe2O3(run,melt_wf,setup,species)
    tot = (setup.loc[run,"SiO2"] + setup.loc[run,"TiO2"] + setup.loc[run,"Al2O3"] + Wm_FeOT_ + Wm_FeO_ + Wm_Fe2O3_ + setup.loc[run,"MgO"] + setup.loc[run,"MnO"] + setup.loc[run,"CaO"] + setup.loc[run,"Na2O"] + setup.loc[run,"K2O"] + setup.loc[run,"P2O5"])
    SiO2 = (setup.loc[run,"SiO2"]/tot)*(1.-volatiles)
    TiO2 = (setup.loc[run,"TiO2"]/tot)*(1.-volatiles)
    Al2O3 = (setup.loc[run,"Al2O3"]/tot)*(1.-volatiles)
    FeOT = (Wm_FeOT_/tot)*(1.-volatiles)
    FeO = (Wm_FeO_/tot)*(1.-volatiles)
    Fe2O3 = (Wm_Fe2O3_/tot)*(1.-volatiles)
    MgO = (setup.loc[run,"MgO"]/tot)*(1.-volatiles)
    MnO = (setup.loc[run,"MnO"]/tot)*(1.-volatiles)
    CaO = (setup.loc[run,"CaO"]/tot)*(1.-volatiles)
    Na2O = (setup.loc[run,"Na2O"]/tot)*(1.-volatiles)
    K2O = (setup.loc[run,"K2O"]/tot)*(1.-volatiles)
    P2O5 = (setup.loc[run,"P2O5"]/tot)*(1.-volatiles)
    return SiO2, TiO2, Al2O3, FeOT, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, H2O, CO2, S, X

# calculate cation proportions
def melt_cation_proportion(run,melt_wf,setup,species,volatiles,Fe_speciation):
    SiO2, TiO2, Al2O3, FeOT, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, H2O, CO2, S, X = melt_normalise_wf(run,melt_wf,setup,species,volatiles,Fe_speciation)
    tot = ((species.loc["SiO2","no_cat"]*SiO2)/species.loc["SiO2","M"]) + ((species.loc["TiO2","no_cat"]*TiO2)/species.loc["TiO2","M"]) + ((species.loc["Al2O3","no_cat"]*Al2O3)/species.loc["Al2O3","M"]) + ((species.loc["FeO","no_cat"]*FeOT)/species.loc["FeO","M"]) + ((species.loc["FeO","no_cat"]*FeO)/species.loc["FeO","M"]) + ((species.loc["Fe2O3","no_cat"]*Fe2O3)/species.loc["Fe2O3","M"]) + ((species.loc["MgO","no_cat"]*MgO)/species.loc["MgO","M"]) + ((species.loc["MnO","no_cat"]*MnO)/species.loc["MnO","M"]) + ((species.loc["CaO","no_cat"]*CaO)/species.loc["CaO","M"]) + ((species.loc["Na2O","no_cat"]*Na2O)/species.loc["Na2O","M"]) + ((species.loc["K2O","no_cat"]*K2O)/species.loc["K2O","M"]) + ((species.loc["P2O5","no_cat"]*P2O5)/species.loc["P2O5","M"]) + ((species.loc["H2O","no_cat"]*H2O)/species.loc["H2O","M"]) + ((species.loc["CO2","no_cat"]*CO2)/species.loc["CO2","M"])
    Si = ((species.loc["SiO2","no_cat"]*SiO2)/species.loc["SiO2","M"])/tot
    Ti = ((species.loc["TiO2","no_cat"]*TiO2)/species.loc["TiO2","M"])/tot
    Al = ((species.loc["Al2O3","no_cat"]*Al2O3)/species.loc["Al2O3","M"])/tot
    FeT = ((species.loc["FeO","no_cat"]*FeOT)/species.loc["FeO","M"])/tot
    Fe2 = ((species.loc["FeO","no_cat"]*FeO)/species.loc["FeO","M"])/tot
    Fe3 = ((species.loc["Fe2O3","no_cat"]*Fe2O3)/species.loc["Fe2O3","M"])/tot    
    Mg = ((species.loc["MgO","no_cat"]*MgO)/species.loc["MgO","M"])/tot
    Mn = ((species.loc["MgO","no_cat"]*MnO)/species.loc["MnO","M"])/tot
    Ca = ((species.loc["CaO","no_cat"]*CaO)/species.loc["CaO","M"])/tot
    Na = ((species.loc["Na2O","no_cat"]*Na2O)/species.loc["Na2O","M"])/tot
    K = ((species.loc["K2O","no_cat"]*K2O)/species.loc["K2O","M"])/tot
    P = ((species.loc["P2O5","no_cat"]*P2O5)/species.loc["P2O5","M"])/tot
    H = ((species.loc["H2O","no_cat"]*H2O)/species.loc["H2O","M"])/tot
    C = ((species.loc["CO2","no_cat"]*CO2)/species.loc["CO2","M"])/tot
    return tot, Si, Ti, Al, FeT, Fe2, Fe3, Mg, Mn, Ca, Na, K, P, H, C

# calculate mole fractions in the melt
def melt_mole_fraction(run,melt_wf,setup,species,models,volatiles,Fe_speciation):
    SiO2, TiO2, Al2O3, FeOT, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, H2O, CO2, S, X = melt_normalise_wf(run,melt_wf,setup,species,volatiles,Fe_speciation)
    species_X = models.loc["species X","option"]
    mol_tot = (SiO2/species.loc["SiO2","M"]) + (TiO2/species.loc["TiO2","M"]) + (Al2O3/species.loc["Al2O3","M"]) + (FeOT/species.loc["FeO","M"]) + (FeO/species.loc["FeO","M"]) + (Fe2O3/species.loc["Fe2O3","M"]) + (MnO/species.loc["MnO","M"]) + (MgO/species.loc["MgO","M"]) + (CaO/species.loc["CaO","M"]) + (Na2O/species.loc["Na2O","M"]) + (K2O/species.loc["K2O","M"]) + (P2O5/species.loc["P2O5","M"]) + (H2O/species.loc["H2O","M"]) + (CO2/species.loc["CO2","M"]) + (S/species.loc["S","M"]) + (X/species.loc[species_X,"M"])
    SiO2 = (SiO2/species.loc["SiO2","M"])/mol_tot
    TiO2 = (TiO2/species.loc["TiO2","M"])/mol_tot
    Al2O3 = (Al2O3/species.loc["Al2O3","M"])/mol_tot 
    FeOT = (FeOT/species.loc["FeO","M"])/mol_tot 
    FeO = (FeO/species.loc["FeO","M"])/mol_tot
    Fe2O3 = (Fe2O3/species.loc["Fe2O3","M"])/mol_tot
    MnO = (MnO/species.loc["MnO","M"])/mol_tot
    MgO = (MgO/species.loc["MgO","M"])/mol_tot
    CaO = (CaO/species.loc["CaO","M"])/mol_tot
    Na2O = (Na2O/species.loc["Na2O","M"])/mol_tot
    K2O = (K2O/species.loc["K2O","M"])/mol_tot
    P2O5 = (P2O5/species.loc["P2O5","M"])/mol_tot
    H2O = (H2O/species.loc["H2O","M"])/mol_tot
    CO2 = (CO2/species.loc["CO2","M"])/mol_tot
    S = (S/species.loc["S","M"])/mol_tot
    X = (X/species.loc[species_X,"M"])/mol_tot
    return mol_tot, SiO2, TiO2, Al2O3, FeOT, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, H2O, CO2, S, X   

def melt_single_O(run,melt_wf,setup,species,volatiles,Fe_speciation):
    SiO2, TiO2, Al2O3, FeOT, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, H2O, CO2, S, X = melt_normalise_wf(run,melt_wf,setup,species,volatiles,Fe_speciation)
    Xmtot = (SiO2/(species.loc["SiO2","M"]/species.loc["SiO2","no_O"])) + (TiO2/(species.loc["TiO2","M"]/species.loc["TiO2","no_O"])) + (Al2O3/(species.loc["Al2O3","M"]/species.loc["Al2O3","no_O"])) + (MnO/(species.loc["MnO","M"]/species.loc["MnO","no_O"])) + (MgO/(species.loc["MgO","M"]/species.loc["MgO","no_O"])) + (CaO/(species.loc["CaO","M"]/species.loc["CaO","no_O"])) + (Na2O/(species.loc["Na2O","M"]/species.loc["Na2O","no_O"])) + (K2O/(species.loc["K2O","M"]/species.loc["K2O","no_O"])) + (P2O5/(species.loc["P2O5","M"]/species.loc["P2O5","no_O"])) + (FeOT/(species.loc["FeO","M"]/species.loc["FeO","no_O"])) + (FeO/(species.loc["FeO","M"]/species.loc["FeO","no_O"])) + (Fe2O3/(species.loc["Fe2O3","M"]/species.loc["Fe2O3","no_O"])) + (H2O/(species.loc["H2O","M"]/species.loc["H2O","no_O"])) + (CO2/(species.loc["CO2","M"]/species.loc["CO2","no_O"]))
    SiO2 = (SiO2/species.loc["SiO2","M"])/Xmtot
    TiO2 = (TiO2/species.loc["TiO2","M"])/Xmtot
    Al2O3 = (Al2O3/species.loc["Al2O3","M"])/Xmtot 
    FeOT = (FeOT/species.loc["FeO","M"])/Xmtot 
    FeO = (FeO/species.loc["FeO","M"])/Xmtot
    Fe2O3 = (Fe2O3/species.loc["Fe2O3","M"])/Xmtot
    MnO = (MnO/species.loc["MnO","M"])/Xmtot
    MgO = (MgO/species.loc["MgO","M"])/Xmtot
    CaO = (CaO/species.loc["CaO","M"])/Xmtot
    Na2O = (Na2O/species.loc["Na2O","M"])/Xmtot
    P2O5 = (P2O5/species.loc["P2O5","M"])/Xmtot
    K2O = (K2O/species.loc["K2O","M"])/Xmtot
    H2O = (H2O/species.loc["H2O","M"])/Xmtot
    CO2 = (CO2/species.loc["CO2","M"])/Xmtot    
    return Xmtot, SiO2, TiO2, Al2O3, FeOT, FeO, Fe2O3, MgO, MnO, CaO, Na2O, K2O, P2O5, H2O, CO2
