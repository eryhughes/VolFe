# calculations.py

import pandas as pd
from datetime import date
import gmpy2 as gp
import numpy as np
import datetime

import melt_gas as mg
import equilibrium_equations as eq
import model_dependent_variables as mdv


###########################
### saturation pressure ###
###########################

# for a given melt composition, calcualte the saturation pressure
def P_sat(run,PT,melt_wf,setup,species,models,Ptol,nr_step,nr_tol):
    ST = melt_wf["ST"]
    H2OT = melt_wf["H2OT"]
    CO2 = melt_wf["CO2"]
    HT = melt_wf["HT"]
    CT = melt_wf["CT"]
    XT = melt_wf["XT"]
    melt_wf1 = {"ST":ST,"CO2":CO2,"H2OT":H2OT,"HT":HT,"CT":CT, "XT":XT} # to work out P_sat
    melt_wf2 = {"ST":ST,"CO2":CO2,"H2OT":H2OT,"HT":HT,"CT":CT, "XT":XT} # to work out sulphur saturation

    def Pdiff(guess,run,melt_wf,setup,species,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(run,PT,melt_wf,setup,species,models))
        return difference

    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    melt_wf1["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
    melt_wf2["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
    #wt_C, wt_O, wt_H, wt_S, wt_Fe, wt_g, Wt = bulk_composition(run,PT,melt_wf1,setup,species,models)
    #bulk_wf = {"H":wt_H,"C":wt_C,"S":wt_S}
    xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf1,setup,species,models,nr_step,nr_tol)
    melt_wf1["H2OT"] = wm_H2O_
    melt_wf2["H2OT"] = wm_H2O_
    melt_wf1["CO2"] = wm_CO2_
    melt_wf2["CO2"] = wm_CO2_
    melt_wf1["S2-"] = wm_S2m_
    melt_wf2["S2-"] = wm_S2m_
    melt_wf1["ST"] = ST
    melt_wf2["ST"] = ST
    if models.loc["sulphur_saturation","option"] == "yes": # must incorporate H2S concentration into S2- for SCSS
        SCSS_,sulphide_sat,SCAS_,sulphate_sat, ss_ST = sulphur_saturation(run,PT,melt_wf2,setup,species,models)
        melt_wf1["ST"] = ss_ST/1000000.
        melt_wf2["ST"] = ST
    delta1 = Pdiff(guess0,run,melt_wf1,setup,species,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,run,melt_wf1,setup,species,models)
        guess0 = mg.p_tot(run,PT,melt_wf1,setup,species,models)
        guess0 = float(guess0)
        PT["P"] = guess0
        melt_wf1["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
        melt_wf2["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
        xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf1,setup,species,models,nr_step,nr_tol)
        melt_wf1["H2OT"] = wm_H2O_
        melt_wf2["H2OT"] = wm_H2O_
        melt_wf1["CO2"] = wm_CO2_
        melt_wf2["CO2"] = wm_CO2_
        melt_wf1["S2-"] = wm_S2m_
        melt_wf2["S2-"] = wm_S2m_
        if models.loc["sulphur_saturation","option"] == "yes":
            SCSS_,sulphide_sat,SCAS_,sulphate_sat,ss_ST = sulphur_saturation(run,PT,melt_wf2,setup,species,models)
            melt_wf1["ST"] = ss_ST/1000000.
            melt_wf2["ST"] = ST
    else:
        P_sat = guess0
        xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf1,setup,species,models,nr_step,nr_tol)
        
    melt_wf["ST"] = ST
    return P_sat, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST

# for a given melt composition, calcualte the saturation pressure
def P_sat_H2O_CO2(run,PT,melt_wf,setup,species,models,Ptol,nr_step,nr_tol): # Pvsat with just H2O and CO2 in vapour
    
    def p_tot_H2O_CO2(run,PT,melt_wf,setup,species,models):
        value = mg.p_H2O(run,PT,melt_wf,setup,species,models) + mg.p_CO2(run,PT,melt_wf,setup,species,models)  
        return value
    
    def Pdiff(guess,run,melt_wf,setup,species,models):
        PT["P"] = guess
        difference = abs(guess - p_tot_H2O_CO2(run,PT,melt_wf,setup,species,models))
        return difference

    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
        guess0 = p_tot_H2O_CO2(run,PT,melt_wf,setup,species,models)
        guess0 = float(guess0)
        PT["P"] = guess0
    else:
        P_sat = guess0
        xg_H2O_ = mg.xg_H2O(run,PT,melt_wf,setup,species,models)
        xg_CO2_ = mg.xg_CO2(run,PT,melt_wf,setup,species,models)
        p_H2O_ = mg.p_H2O(run,PT,melt_wf,setup,species,models)
        p_CO2_ = mg.p_CO2(run,PT,melt_wf,setup,species,models)
        f_H2O_ = mg.f_H2O(run,PT,melt_wf,setup,species,models)
        f_CO2_ = mg.f_CO2(run,PT,melt_wf,setup,species,models)
        xm_H2O_ = (mdv.C_H2O(run,PT,melt_wf,setup,species,models)*f_H2O_)**0.5
        xm_CO2_ = mdv.C_CO3(run,PT,melt_wf,setup,species,models)*f_CO2_
        M_m_ = mg.M_m_SO(run,melt_wf,setup,species)
        Xm_t = xm_CO2_*species.loc["CO2","M"] + xm_H2O_*species.loc["H2O","M"] + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*species.loc["H2O","M"])/Xm_t
        wm_CO2_ = (xm_CO2_*species.loc["CO2","M"])/Xm_t

    return P_sat, xg_H2O_, xg_CO2_, f_H2O_, f_CO2_, p_H2O_, p_CO2_, wm_H2O_, wm_CO2_
        
        
        
###################
### composition ###
###################

# calculate bulk composition, including if a gas phase is present
def bulk_composition(run,PT,melt_wf,setup,species,models):
    bulk_composition = models.loc["bulk_composition","option"]
    eq_Fe = models.loc["eq_Fe","option"]
    wm_CO2 = melt_wf["CO2"]
    wm_H2O = melt_wf["H2OT"]
    wm_H2 = melt_wf["H2"]
    wm_CO = melt_wf["CO"]
    wm_CH4 = melt_wf["CH4"]
    wm_H2S = melt_wf["H2S"]
    wm_S6p = melt_wf["S6+"]
    wm_S2m = melt_wf["S2-"]
    wm_ST = melt_wf["ST"]
    wm_XT = melt_wf["XT"]
    Fe3FeT = melt_wf["Fe3FeT"]
    Fe3Fe2_ = mg.Fe3Fe2(melt_wf)
    S6ST_ = mg.S6ST(run,PT,melt_wf,setup,species,models)
    #SCSS_,sulphide_sat,SCAS_,sulphate_sat = sulphur_saturation(wm_ST/100.0,S6ST_)
    #print(P, S6ST_)

    if bulk_composition == "yes":
        wt_g = 0.
    elif bulk_composition == "wtg":
        wt_g = setup.loc[run,"wt_g"]/100.
    elif bulk_composition == "CO2":
        wt_C_ = ((species.loc['C','M']*(setup.loc[run,"initial_CO2wtpc"]/100.))/species.loc['CO2','M'])
        wt_g = ((wt_C_/species.loc["C","M"]) - (wm_CO2/species.loc["CO2","M"]))/(((mg.xg_CO2(run,PT,melt_wf,setup,species,models)+mg.xg_CO(run,PT,melt_wf,setup,species,models)+mg.xg_CH4(run,PT,melt_wf,setup,species,models)+mg.xg_OCS(run,PT,melt_wf,setup,species,models))/mg.Xg_tot(run,PT,melt_wf,setup,species,models)) - (wm_CO2/species.loc["CO2","M"]))    

    if bulk_composition == "CO2":
        wt_C = wt_C_
    else:
        wt_C = species.loc["C","M"]*((wt_g*(((mg.xg_CO2(run,PT,melt_wf,setup,species,models)+
                                      mg.xg_CO(run,PT,melt_wf,setup,species,models)+
                                      mg.xg_CH4(run,PT,melt_wf,setup,species,models)+
                                      mg.xg_OCS(run,PT,melt_wf,setup,species,models))
                                     /mg.Xg_tot(run,PT,melt_wf,setup,species,models)) 
                                    - (wm_CO2/species.loc["CO2","M"] - wm_CO/species.loc["CO","M"] - 
                                       wm_CH4/species.loc["CH4","M"]))) + (wm_CO2/species.loc["CO2","M"] + 
                                        wm_CO/species.loc["CO","M"] + wm_CH4/species.loc["CH4","M"]))

    wt_H = 2.*species.loc["H","M"]*((wt_g*(((mg.xg_H2O(run,PT,melt_wf,setup,species,models)+mg.xg_H2(run,PT,melt_wf,setup,species,models)+2.*mg.xg_CH4(run,PT,melt_wf,setup,species,models)+mg.xg_H2S(run,PT,melt_wf,setup,species,models))/mg.Xg_tot(run,PT,melt_wf,setup,species,models)) - (wm_H2O/species.loc["H2O","M"] - wm_H2/species.loc["H2","M"] - wm_H2S/species.loc["H2S","M"] - (2.*wm_CH4)/species.loc["CH4","M"]))) + (wm_H2O/species.loc["H2O","M"] + wm_H2/species.loc["H2","M"] + wm_H2S/species.loc["H2S","M"] + (2.*wm_CH4)/species.loc["CH4","M"]))
    
    wt_S = species.loc["S","M"]*((wt_g*(((mg.xg_SO2(run,PT,melt_wf,setup,species,models)
                                      +2.0*mg.xg_S2(run,PT,melt_wf,setup,species,models)
                                      +mg.xg_H2S(run,PT,melt_wf,setup,species,models)
                                      +mg.xg_OCS(run,PT,melt_wf,setup,species,models))
                                     /mg.Xg_tot(run,PT,melt_wf,setup,species,models)) 
                                    - ((wm_S2m+wm_S6p)/species.loc["S","M"] - (wm_H2S/species.loc["H2S","M"])))) 
                              + ((wm_S2m+wm_S6p)/species.loc["S","M"] + (wm_H2S/species.loc["H2S","M"])))
    
    species_X = models.loc["species X","option"]
    wt_X = species.loc[species_X,"M"]*((wt_g*(((mg.xg_X(run,PT,melt_wf,setup,species,models))
                                     /mg.Xg_tot(run,PT,melt_wf,setup,species,models)) 
                                    - wm_XT/species.loc[species_X,"M"])) 
                              + wm_XT/species.loc[species_X,"M"])
    
    if models.loc["mass_volume","option"] == "mass":
        Wt = setup.loc[run, "total_mass_g"]
    elif models.loc["mass_volume","option"] == "volume": ### THIS NEEDS FIXING ###
        Wt = 0.

    if eq_Fe == "no":
        wt_Fe = 0.0
    elif eq_Fe == "yes":
        total_dissolved_volatiles = (wm_CO2 + wm_H2O + wm_ST*(1.-S6ST_) + (species.loc["SO3","M"]*((wm_ST*S6ST_)/species.loc["S","M"])))
        wt_Fe = (1.-wt_g)*(((1.0-total_dissolved_volatiles)*mg.wm_Fe_nv(run,melt_wf,setup,species))/100.0) # wt fraction of Fe

    wt_O = species.loc["O","M"]*((wt_g*(((2.0*mg.xg_CO2(run,PT,melt_wf,setup,species,models) + mg.xg_CO(run,PT,melt_wf,setup,species,models) + 2.0*mg.xg_O2(run,PT,melt_wf,setup,species,models) + mg.xg_H2O(run,PT,melt_wf,setup,species,models) + 2.0*mg.xg_SO2(run,PT,melt_wf,setup,species,models) + mg.xg_OCS(run,PT,melt_wf,setup,species,models))/mg.Xg_tot(run,PT,melt_wf,setup,species,models)) - (wm_H2O/species.loc["H2O","M"]) - ((2.0*wm_CO2)/species.loc["CO2","M"]) - (3.0*wm_S6p/species.loc["S","M"] - (wm_CO/species.loc["CO","M"])))) + (wm_H2O/species.loc["H2O","M"]) + ((2.0*wm_CO2)/species.loc["CO2","M"]) + (3.0*wm_S6p/species.loc["S","M"]) + (wm_CO/species.loc["CO","M"]) + (wt_Fe/species.loc["Fe","M"])*((1.5*Fe3Fe2_+1.0)/(Fe3Fe2_+1.0)))
    return wt_C, wt_O, wt_H, wt_S, wt_X, wt_Fe, wt_g, Wt

# calculate weight fraction of elements in the system when adding gas into a melt
def new_bulk_regas_open(run,PT,melt_wf,bulk_wf,gas_mf,dwtg,setup,species,models):
    wm_C, wm_H, wm_S, wm_Fe, wm_O, wm_X = mg.melt_elements(run,PT,melt_wf,bulk_wf,gas_mf,setup,species,models)
    wg_O, wg_H, wg_S, wg_C, wg_X = mg.gas_elements(gas_mf,species)
    wt_C = (1.-dwtg)*wm_C + dwtg*wg_C
    wt_H = (1.-dwtg)*wm_H + dwtg*wg_H
    wt_O = (1.-dwtg)*wm_O + dwtg*wg_O
    wt_S = (1.-dwtg)*wm_S + dwtg*wg_S
    wt_X = (1.-dwtg)*wm_X + dwtg*wg_X
    wt_Fe = (1.-dwtg)*wm_Fe
    Wt = bulk_wf['Wt']*(1. + dwtg)
    return wt_C, wt_H, wt_S, wt_Fe, wt_O, wt_X, Wt


#########################
### sulphur satuation ###
#########################

# check solid/immiscible liquid sulphur saturation
def sulphur_saturation(run,PT,melt_wf,setup,species,models): # melt weight fraction of ST and S6/ST
    wmST = melt_wf['ST']
    S6T = mg.S6ST(run,PT,melt_wf,setup,species,models)
    wmS2 = wmST*100.0*10000.0*(1.0-S6T)
    wmS6 = wmST*100.0*10000.0*S6T
    SCSS_ = mdv.SCSS(run,PT,melt_wf,setup,species,models)
    SCAS_ = mdv.SCAS(run,PT,melt_wf,setup,species,models)
    StCSS = SCSS_/(1.-S6T)
    StCAS = SCAS_/S6T
    if wmS2 < SCSS_ and wmS6 < SCAS_:
        sulphide_sat = "no"
        sulphate_sat = "no"
        ST = wmST*1000000.
    elif wmS2 >= SCSS_ and wmS6 >= SCAS_:
        sulphide_sat = "yes"
        sulphate_sat = "yes"
        ST = min(StCSS,StCAS)
    elif wmS2 >= SCSS_ and wmS6 < SCAS_:
        sulphide_sat = "yes"
        sulphate_sat = "no"
        ST = StCSS
    elif wmS2 < SCSS_ and wmS6 >= SCAS_:
        sulphide_sat = "no"
        sulphate_sat = "yes"
        ST = StCAS
    else:
        sulphide_sat = "nan"
        sulphate_sat = "nan"
        ST = wmST*1000000.
    return SCSS_, sulphide_sat, SCAS_, sulphate_sat, ST

# fO2 and P of v+sulf+anh saturation
def fO2_P_VSA(run,PT,melt_wf,setup,species,models,nr_step,nr_tol,Ptol):

    def Pdiff(guess,run,melt_wf,setup,species,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(run,PT,melt_wf,setup,species,models))
        return difference

    def fO2_S(run,PT,melt_wf,setup,species,models):
        SCSS_ = mdv.SCSS(run,PT,melt_wf,setup,species,models)/1000000.
        SCAS_ = mdv.SCAS(run,PT,melt_wf,setup,species,models)/1000000.
        CSO4 = mdv.C_SO4(run,PT,melt_wf,setup,species,models)/1000000.
        CS = mdv.C_S(run,PT,melt_wf,setup,species,models)/1000000.
        
        if models.loc["H2S_m","option"] == "no":
            W = CSO4/CS
        elif models.loc["H2S_m","option"] == "yes":
            CH2S = mdv.C_H2S(run,PT,melt_wf,setup,species,models)/1000000.
            KHS = mdv.KHOSg(PT,models)
            CH2OT = mdv.C_H2O(run,PT,melt_wf,setup,species,models)
            xmH2O = mg.xm_H2OT_so(run,melt_wf,setup,species)
            W = (CSO4/((KHS*CH2S*(xmH2O**2./CH2OT)) + CS))
        fO2 = ((1./W)*(SCAS_/SCSS_))**0.5
        ST = SCSS_ + SCAS_
        return fO2, ST
    
    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    fO2_, ST_ = fO2_S(run,PT,melt_wf,setup,species,models)
    melt_wf["ST"] = ST_
    melt_wf["Fe3FeT"] = mdv.fO22Fe3FeT(fO2_,run,PT,melt_wf,setup,species,models)
    xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)
    melt_wf["H2OT"] = wm_H2O_
    melt_wf["CO2"] = wm_CO2_
    melt_wf["S2-"] = wm_S2m_
    delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
        guess0 = mg.p_tot(run,PT,melt_wf,setup,species,models)
        guess0 = float(guess0)
        PT["P"] = guess0
        fO2_, ST_ = fO2_S(run,PT,melt_wf,setup,species,models)
        melt_wf["ST"] = ST_
        melt_wf["Fe3FeT"] = mdv.fO22Fe3FeT(fO2_,run,PT,melt_wf,setup,species,models)
        xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)         
        melt_wf["H2OT"] = wm_H2O_
        melt_wf["CO2"] = wm_CO2_
        melt_wf["S2-"] = wm_S2m_
    else:
        P_sat = guess0
        xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)
        fO2_, ST_ = fO2_S(run,PT,melt_wf,setup,species,models)
        Fe3_FT = mdv.fO22Fe3FeT(fO2_,run,PT,melt_wf,setup,species,models)
    return P_sat, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, Fe3_FT

def P_VSA(run,PT,melt_wf,setup,species,models,nr_step,nr_tol,Ptol):
    P_sat, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, Fe3_FT = fO2_P_VSA(run,PT,melt_wf,setup,species,models,nr_step,nr_tol,Ptol)
    
    def Pdiff(guess,run,melt_wf,setup,species,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(run,PT,melt_wf,setup,species,models))
        return difference
    
    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    melt_wf["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
    if melt_wf["Fe3FeT"] < Fe3_FT:
        SCSS_ = mdv.SCSS(run,PT,melt_wf,setup,species,models)/1000000.
        S6T = mg.S6ST(run,PT,melt_wf,setup,species,models)
        ST_ = SCSS_/(1.-S6T)
    else:
        SCAS_ = mdv.SCAS(run,PT,melt_wf,setup,species,models)/1000000.
        S6T = mg.S6ST(run,PT,melt_wf,setup,species,models)
        ST_ = SCAS_/S6T
    melt_wf["ST"] = ST_
    xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)
    melt_wf["H2OT"] = wm_H2O_
    melt_wf["CO2"] = wm_CO2_
    melt_wf["S2-"] = wm_S2m_
    delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
        guess0 = mg.p_tot(run,PT,melt_wf,setup,species,models)
        guess0 = float(guess0)
        PT["P"] = guess0
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
        if melt_wf["Fe3FeT"] < Fe3_FT:
            SCSS_ = mdv.SCSS(run,PT,melt_wf,setup,species,models)/1000000.
            S6T = mg.S6ST(run,PT,melt_wf,setup,species,models)
            ST_ = SCSS_/(1.-S6T)
        else:
            SCAS_ = mdv.SCAS(run,PT,melt_wf,setup,species,models)/1000000.
            S6T = mg.S6ST(run,PT,melt_wf,setup,species,models)
            ST_ = SCAS_/S6T
        melt_wf["ST"] = ST_
        xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)         
        melt_wf["H2OT"] = wm_H2O_
        melt_wf["CO2"] = wm_CO2_
        melt_wf["S2-"] = wm_S2m_
    else:
        P_sat = guess0
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
        xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)
    return P_sat, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST

# estimated fO2 from sulfur content given you are sulfide saturated and at Pvsat
def P_sat_sulf_anh(run,PT,melt_wf,setup,species,models,Ptol,nr_step,nr_tol):
    ST = melt_wf["ST"]
    H2OT = melt_wf["H2OT"]
    CO2 = melt_wf["CO2"]
    HT = melt_wf["HT"]
    CT = melt_wf["CT"]
 
    def Pdiff(guess,run,melt_wf,setup,species,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(run,PT,melt_wf,setup,species,models))
        return difference
    
    def S62_2_Fe3T(run,PT,melt_wf,setup,species,models,sat):
        ST = melt_wf["ST"]
        if sat == "sulf":
            Ssat = mdv.SCSS(run,PT,melt_wf,setup,species,models)/1000000.
            S2m = Ssat
            S6p = ST - S2m
        elif sat == "anh":
            Ssat = mdv.SCAS(run,PT,melt_wf,setup,species,models)/1000000.
            S6p = Ssat
            S2m = ST - S6p
        S62 = S6p/S2m
        S6T = S6p/ST
        if S62 < 0.:
            return "not possible","","","",Ssat
        else:
            fO2 = mg.S6S2_2_fO2(S62,melt_wf,run,PT,setup,species,models)
            Fe3T = mdv.fO22Fe3FeT(fO2,run,PT,melt_wf,setup,species,models)
            DFMQ = mg.fO22Dbuffer(PT,fO2,"FMQ",models)
            return Fe3T,fO2,S6T,DFMQ,Ssat
        
    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    
    # assume it is sulfide saturated
    melt_wf["Fe3FeT"] = 0.1
    Fe3T_sulf,fO2_sulf,S6T_sulf,DFMQ_sulf,SCSS_ = S62_2_Fe3T(run,PT,melt_wf,setup,species,models,"sulf")
    if Fe3T_sulf == "not possible":
        P_sat_sulf = ""
        Fe3T_sulf = ""
        sulphide_sat = "no"
    else: 
        melt_wf["Fe3FeT"] = Fe3T_sulf
        xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)
        melt_wf["H2OT"] = wm_H2O_
        melt_wf["CO2"] = wm_CO2_
        melt_wf["S2-"] = wm_S2m_
        melt_wf["ST"] = ST
        delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
        while delta1 > Ptol :
            delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
            guess0 = mg.p_tot(run,PT,melt_wf,setup,species,models)
            guess0 = float(guess0)
            PT["P"] = guess0
            Fe3T_sulf,fO2_sulf,S6T_sulf,DFMQ_sulf,SCSS_ = S62_2_Fe3T(run,PT,melt_wf,setup,species,models,"sulf")
            if Fe3T_sulf == "not possible":
                P_sat_sulf = ""
                Fe3T_sulf = ""
                sulphide_sat = "no"
            else: 
                melt_wf["Fe3FeT"] = Fe3T_sulf
                xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)
                melt_wf["H2OT"] = wm_H2O_
                melt_wf["CO2"] = wm_CO2_
                melt_wf["S2-"] = wm_S2m_
        else:
            P_sat_sulf = guess0
            Fe3T_sulf,fO2_sulf,S6T_sulf,DFMQ_sulf,SCSS_ = S62_2_Fe3T(run,PT,melt_wf,setup,species,models,"sulf")
            sulphide_sat = "yes"
           
    # assume it is anhydrite saturated
    Fe3T_anh,fO2_anh,S6T_anh,DFMQ_anh,SCAS_ = S62_2_Fe3T(run,PT,melt_wf,setup,species,models,"anh")
    if Fe3T_anh == "not possible":
        P_sat_anh = ""
        Fe3T_anh = ""
        anhydrite_sat = "no"
    else:  
        melt_wf["Fe3FeT"] = Fe3T_anh
        xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)
        melt_wf["H2OT"] = wm_H2O_
        melt_wf["CO2"] = wm_CO2_
        melt_wf["S2-"] = wm_S2m_
        melt_wf["ST"] = ST
        delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
        while delta1 > Ptol :
            delta1 = Pdiff(guess0,run,melt_wf,setup,species,models)
            guess0 = mg.p_tot(run,PT,melt_wf,setup,species,models)
            guess0 = float(guess0)
            PT["P"] = guess0
            Fe3T_anh,fO2_anh,S6T_anh,DFMQ_anh,SCAS_ = S62_2_Fe3T(run,PT,melt_wf,setup,species,models,"anh")
            if Fe3T_anh == "not possible":
                P_sat_anh = ""
                Fe3T_anh = ""
                anhydrite_sat = "no"
            else:  
                melt_wf["Fe3FeT"] = Fe3T_anh
                xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S2m_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)
                melt_wf["H2OT"] = wm_H2O_
                melt_wf["CO2"] = wm_CO2_
                melt_wf["S2-"] = wm_S2m_
        else:
            P_sat_anh = guess0
            Fe3T_anh,fO2_anh,S6T_anh,DFMQ_anh,SCAS_ = S62_2_Fe3T(run,PT,melt_wf,setup,species,models,"anh")
            anhydrite_sat = "yes"
        

    return P_sat_sulf,P_sat_anh,SCAS_*1000000.,SCSS_*1000000.,sulphide_sat,DFMQ_sulf,fO2_sulf,Fe3T_sulf,S6T_sulf,anhydrite_sat,DFMQ_anh,fO2_anh,Fe3T_anh,S6T_anh

##########################
### graphite satuation ###
##########################

# check graphite saturation
def graphite_saturation(run,PT,melt_wf,setup,species,models): # needs finishing
    K1 = mg.f_CO2(run,PT,melt_wf,setup,species,models)/mdv.f_O2(run,PT,melt_wf,setup,species,models)
    K2 = mdv.KCOs(PT,models) # K for graphite saturation
    if K1 < K2:
        graphite_sat = "no"
    else:
        graphite_sat = "yes"
    fCO2_ = K2*mdv.f_O2(run,PT,melt_wf,setup,species,models)
    xmCO2 = fCO2_*mdv.C_CO3(run,PT,melt_wf,setup,species,models)
    return graphite_sat

                
######################################                
### fO2 range from sulphur content ###
######################################
                
def fO2_range_from_S(run,PT,melt_wf,setup,species,models):
    SCSS_ = mdv.SCSS(run,PT,melt_wf,setup,species,models)
    SCAS_ = mdv.SCAS(run,PT,melt_wf,setup,species,models)

    if setup.loc[run,"STppm"] > SCSS_:
        sulphide_sat = "possible"
        S6ST_1 = (setup.loc[run,"STppm"] - SCSS_)/setup.loc[run,"STppm"]
        S6S2_1 = mg.overtotal2ratio(S6ST_1)
        fO2_1 = mg.S6S2_2_fO2(S6S2_1,melt_wf,run,PT,setup,species,models)
        Fe3FeT_1 = mdv.fO22Fe3FeT(fO2_1,run,PT,melt_wf,setup,species,models)
        Fe3Fe2_1 = mg.overtotal2ratio(Fe3FeT_1)
        DFMQ_1 = mg.fO22Dbuffer(PT,fO2_1,"FMQ",models)
    else:
        sulphide_sat = "no"
        S6ST_1 = ""
        S6S2_1 = ""
        Fe3Fe2_1 = ""
        Fe3FeT_1 = ""
        fO2_1 = ""
        DFMQ_1 = ""
                   
    if setup.loc[run,"STppm"] > SCAS_:
        sulphate_sat = "possible"
        S6ST_2 = (setup.loc[run,"STppm"] - SCAS_)/setup.loc[run,"STppm"]
        S6S2_2 = mg.overtotal2ratio(S6ST_2)
        fO2_2 = mg.S6S2_2_fO2(S6S2_2,melt_wf,run,PT,setup,species,models)
        Fe3FeT_2 = mdv.fO22Fe3FeT(fO2_2,run,PT,melt_wf,setup,species,models)
        Fe3Fe2_2 = mg.overtotal2ratio(Fe3FeT_2)
        DFMQ_2 = mg.fO22Dbuffer(PT,fO2_2,"FMQ",models)
    else:
        sulphate_sat = "no"
        S6ST_2 = ""
        S6S2_2 = ""
        Fe3Fe2_2 = ""
        Fe3FeT_2 = ""
        fO2_2 = ""
        DFMQ_2 = ""

    return SCAS_,SCSS_,sulphide_sat,DFMQ_1,fO2_1,Fe3FeT_1,S6ST_1,sulphate_sat,DFMQ_2,fO2_2,Fe3FeT_2,S6ST_2
    
                             

###################################
### concentration of insolubles ### 
###################################

def conc_insolubles(run,PT,melt_wf,setup,species,models):
    CO2 = melt_wf["CO2"] # weight fraction CO2
    C_CO2_ = (species.loc["C","M"]*CO2)/species.loc["CO2","M"]
    H2O = melt_wf["H2OT"] # weight fraction H2O
    H_H2O = (2.*species.loc["H","M"]*H2O)/species.loc["H2O","M"]
    H2 = (mg.C_H2(run,PT,melt_wf,setup,species,models)*mg.f_H2(run,PT,melt_wf,setup,species,models))/1000000. # weight fraction H2
    H_H2 = (2.*species.loc["H","M"]*H2)/species.loc["H2","M"]
    CH4 = (mg.C_CH4(PT,models)*mg.f_CH4(run,PT,melt_wf,setup,species,models))/1000000. # weight fraction CH4
    H_CH4 = (4.*species.loc["H","M"]*CH4)/species.loc["CH4","M"]
    C_CH4_ = (species.loc["C","M"]*CH4)/species.loc["CH4","M"]
    CO = (mg.C_CO(PT,models)*mg.f_CO(run,PT,melt_wf,setup,species,models))/1000000. # weight fraction CO
    C_CO_ = (species.loc["C","M"]*CO)/species.loc["CO","M"]
    S2m = melt_wf["S2-"] # weight fraction of S2-
    S6p = (mg.C_SO4(run,PT,melt_wf,setup,species,models)*mdv.f_O2(run,PT,melt_wf,setup,species,models)**2*S2m)/mg.C_S(run,PT,melt_wf,setup,species,models) # weight fraction S6+
    H2S = (mg.C_H2S(run,PT,melt_wf,setup,species,models)*mg.f_H2S(run,PT,melt_wf,setup,species,models))/1000000. # weight fraction H2S
    S_H2S = (species.loc["S","M"]*H2S)/species.loc["H2S","M"]
    H_H2S = (2.*species.loc["H","M"]*H2S)/species.loc["H2S","M"]
    C_T = C_CO_ + C_CH4_ + C_CO2_
    H_T = H_H2O + H_H2 + H_CH4 + H_H2S
    S_T = S_H2S + S2m + S6p
    H2O_HT = H_H2O/H_T
    H2_HT = H_H2/H_T
    CH4_HT = H_CH4/H_T
    H2S_HT = H_H2S/H_T
    CO2_CT = C_CO2_/C_T
    CO_CT = C_CO_/C_T
    CH4_CT = C_CH4_/C_T
    S2m_ST = S2m/S_T
    S6p_ST = S6p/S_T
    H2S_ST = S_H2S/S_T
    return H2, CH4, CO, H2S, S6p, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S2m_ST, S6p_ST, H2S_ST, C_T, H_T, S_T

########################################
### measured parameters within error ### 
########################################

def compositions_within_error(run,setup):
    if setup.loc[run,"TC_sd_type"] == "A": # absolute
        TC_sd = setup.loc[run,"TC_sd"]
    else:
        TC_sd = setup.loc[run,"TC_sd"]*setup.loc[run,"T_C"]
    TC = float(np.random.normal(setup.loc[run,"T_C"],TC_sd,1))

    if setup.loc[run,"SiO2_sd_type"] == "A": # absolute
        SiO2_sd = setup.loc[run,"SiO2_sd"]
    else:
        SiO2_sd = setup.loc[run,"SiO2_sd"]*setup.loc[run,"SiO2"]
    SiO2 = float(np.random.normal(setup.loc[run,"SiO2"],SiO2_sd,1))
        
    if setup.loc[run,"TiO2_sd_type"] == "A": # absolute
        TiO2_sd = setup.loc[run,"TiO2_sd"]
    else:
        TiO2_sd = setup.loc[run,"TiO2_sd"]*setup.loc[run,"TiO2"]
    TiO2 = float(np.random.normal(setup.loc[run,"TiO2"],TiO2_sd,1))
        
    if setup.loc[run,"Al2O3_sd_type"] == "A": # absolute
        Al2O3_sd = setup.loc[run,"Al2O3_sd"]
    else:
        Al2O3_sd = setup.loc[run,"Al2O3_sd"]*setup.loc[run,"Al2O3"]
    Al2O3 = float(np.random.normal(setup.loc[run,"Al2O3"],Al2O3_sd,1))

    if setup.loc[run,"FeOT_sd_type"] == "A": # absolute
        FeOT_sd = setup.loc[run,"FeOT_sd"]
    else:
        FeOT_sd = setup.loc[run,"FeOT_sd"]*setup.loc[run,"FeOT"]
    FeOT = float(np.random.normal(setup.loc[run,"FeOT"],FeOT_sd,1))

    if setup.loc[run,"MnO_sd_type"] == "A": # absolute
        MnO_sd = setup.loc[run,"MnO_sd"]
    else:
        MnO_sd = setup.loc[run,"MnO_sd"]*setup.loc[run,"MnO"]
    MnO = float(np.random.normal(setup.loc[run,"MnO"],MnO_sd,1))

    if setup.loc[run,"MgO_sd_type"] == "A": # absolute
        MgO_sd = setup.loc[run,"MgO_sd"]
    else:
        MgO_sd = setup.loc[run,"MgO_sd"]*setup.loc[run,"MgO"]
    MgO = float(np.random.normal(setup.loc[run,"MgO"],MgO_sd,1))
        
    if setup.loc[run,"CaO_sd_type"] == "A": # absolute
        CaO_sd = setup.loc[run,"CaO_sd"]
    else:
        CaO_sd = setup.loc[run,"CaO_sd"]*setup.loc[run,"CaO"]
    CaO = float(np.random.normal(setup.loc[run,"CaO"],CaO_sd,1))

    if setup.loc[run,"Na2O_sd_type"] == "A": # absolute
        Na2O_sd = setup.loc[run,"Na2O_sd"]
    else:
        Na2O_sd = setup.loc[run,"Na2O_sd"]*setup.loc[run,"Na2O"]
    Na2O = float(np.random.normal(setup.loc[run,"Na2O"],Na2O_sd,1))
       
    if setup.loc[run,"K2O_sd_type"] == "A": # absolute
        K2O_sd = setup.loc[run,"K2O_sd"]
    else:
        K2O_sd = setup.loc[run,"K2O_sd"]*setup.loc[run,"K2O"]
    K2O = float(np.random.normal(setup.loc[run,"K2O"],K2O_sd,1))
        
    if setup.loc[run,"P2O5_sd_type"] == "A": # absolute
        P2O5_sd = setup.loc[run,"P2O5_sd"]
    else:
        P2O5_sd = setup.loc[run,"P2O5_sd"]*setup.loc[run,"P2O5"]
    P2O5 = float(np.random.normal(setup.loc[run,"P2O5"],P2O5_sd,1))
        
    if setup.loc[run,"H2O_sd_type"] == "A": # absolute
        H2O_sd = setup.loc[run,"H2O_sd"]
    else:
        H2O_sd = setup.loc[run,"H2O_sd"]*setup.loc[run,"H2O"]
    H2O = float(np.random.normal(setup.loc[run,"H2O"],H2O_sd,1))
        
    if setup.loc[run,"CO2ppm_sd_type"] == "A": # absolute
        CO2ppm_sd = setup.loc[run,"CO2ppm_sd"]
    else:
        CO2ppm_sd = setup.loc[run,"CO2ppm_sd"]*setup.loc[run,"CO2ppm"]
    CO2ppm = float(np.random.normal(setup.loc[run,"CO2ppm"],CO2ppm_sd,1))
        
    if setup.loc[run,"STppm_sd_type"] == "A": # absolute
        STppm_sd = setup.loc[run,"STppm_sd"]
    else:
        STppm_sd = setup.loc[run,"STppm_sd"]*setup.loc[run,"STppm"]
    STppm = float(np.random.normal(setup.loc[run,"STppm"],STppm_sd,1))
        
    if setup.loc[run,"Fe3FeT_sd_type"] == "A": # absolute
        Fe3FeT_sd = setup.loc[run,"Fe3FeT_sd"]
    else:
        Fe3FeT_sd = setup.loc[run,"Fe3FeT_sd"]*setup.loc[run,"Fe3FeT"]
    Fe3FeT = float(np.random.normal(setup.loc[run,"Fe3FeT"],Fe3FeT_sd,1))
        
    return TC,SiO2,TiO2,Al2O3,FeOT,MnO,MgO,CaO,Na2O,K2O,P2O5,H2O,CO2ppm,STppm,Fe3FeT


def calc_isobar_CO2H2O(run,PT,melt_wf,setup,species,models):
    M_H2O = species.loc['H2O','M']
    M_CO2 = species.loc['CO2','M']
    M_m_ = mg.M_m_SO(run,melt_wf,setup,species)

    xm_CO2_ = PT["P"]*mdv.y_CO2(PT,species,models)*mdv.C_CO3(run,PT,melt_wf,setup,species,models) # pure CO2
    Xm_t = xm_CO2_*M_CO2 + (1.-xm_CO2_)*M_m_
    wm_CO2_0 = (xm_CO2_*M_CO2)/Xm_t
            
    if models.loc["Hspeciation","option"] == "none": # fH2O = xmH2OT^2/CH2O
        xm_H2O_ = (PT["P"]*mdv.y_H2O(PT,species,models)*mdv.C_H2O(run,PT,melt_wf,setup,species,models))**0.5 # pure H2O
    else: # regular or ideal: fH2O = xmH2Omol/CH2O
        print("need to sort")
    Xm_t = xm_H2O_*M_H2O + (1.-xm_H2O_)*M_m_
    wm_H2O_0 = (xm_H2O_*M_H2O)/Xm_t
    xm_H2O_step = xm_H2O_/20.
            
    results = pd.DataFrame([[PT["P"],0.,wm_CO2_0*1000000.]])
            
    for m in range(1,20,1):
        xm_H2O_ = xm_H2O_step*m
        Xm_t = xm_H2O_*M_H2O + (1.-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        melt_wf={"H2OT":wm_H2O_}
        melt_wf["CO2"] = 0.
        pH2O = mg.p_H2O(run,PT,melt_wf,setup,species,models)
        pCO2 = PT["P"] - pH2O
        xm_CO2_ = pCO2*mdv.y_CO2(PT,species,models)*mdv.C_CO3(run,PT,melt_wf,setup,species,models)
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wf_CO2 = (xm_CO2_*M_CO2)/Xm_t
        results1 = pd.DataFrame([[PT["P"],melt_wf["H2OT"]*100.,wf_CO2*1000000.]])
        results = results.append(results1, ignore_index=True)
        
    results1 = pd.DataFrame([[PT["P"],wm_H2O_0*100.,0.]])
    results = results.append(results1, ignore_index=True)
    return results

def calc_pure_solubility(run,PT,melt_wf,setup,species,models):
    
    M_H2O = species.loc['H2O','M']
    M_CO2 = species.loc['CO2','M']
    M_m_ = mg.M_m_SO(run,melt_wf,setup,species)
        
    xm_CO2_ = PT["P"]*mdv.y_CO2(PT,species,models)*mdv.C_CO3(run,PT,melt_wf,setup,species,models) # pure CO2
    Xm_t = xm_CO2_*M_CO2 + (1.-xm_CO2_)*M_m_
    wm_CO2 = (xm_CO2_*M_CO2)/Xm_t
        
    xm_H2O_ = (PT["P"]*mdv.y_H2O(PT,species,models)*mdv.C_H2O(run,PT,melt_wf,setup,species,models))**0.5 # pure H2O
    Xm_t = xm_H2O_*M_H2O + (1.-xm_H2O_)*M_m_
    wm_H2O = (xm_H2O_*M_H2O)/Xm_t
       
    results = pd.DataFrame([[PT["P"],wm_H2O*100.,wm_CO2*1000000.]])    
    
    return results
