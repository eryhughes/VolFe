# calculations.py

import pandas as pd
from datetime import date
import numpy as np
import datetime

import VolFe.melt_gas as mg
import VolFe.equilibrium_equations as eq
import VolFe.isotopes as iso
import VolFe.model_dependent_variables as mdv



###########################
### saturation pressure ###
###########################

# for a given melt composition, calcualte the saturation pressure
def P_sat(PT,melt_wf,species,models,Ptol,nr_step,nr_tol):
    ST = melt_wf["ST"]
    #H2OT = melt_wf["H2OT"]
    #CO2 = melt_wf["CO2"]
    #HT = melt_wf["HT"]
    #CT = melt_wf["CT"]
    #XT = melt_wf["XT"]
    #melt_wf1 = {"ST":ST,"CO2":CO2,"H2OT":H2OT,"HT":HT,"CT":CT, "XT":XT} # to work out P_sat
    #melt_wf2 = {"ST":ST,"CO2":CO2,"H2OT":H2OT,"HT":HT,"CT":CT, "XT":XT} # to work out sulfur saturation
    melt_wf1 = {'SiO2': melt_wf["SiO2"], 'TiO2': melt_wf["TiO2"], 'Al2O3': melt_wf["Al2O3"], 'FeOT': melt_wf["FeOT"], 'Fe2O3T': melt_wf["Fe2O3T"], 'FeO': melt_wf["FeO"], 'Fe2O3': melt_wf["Fe2O3"], 'MgO': melt_wf["MgO"], 'MnO': melt_wf["MnO"], 'CaO': melt_wf["CaO"], 'Na2O': melt_wf["Na2O"], 'K2O': melt_wf["K2O"], 'P2O5': melt_wf["P2O5"], 'logfO2_i': melt_wf["logfO2_i"], 'Fe3FeT_i': melt_wf["Fe3FeT_i"], 'DNNO': melt_wf["DNNO"], 'DFMQ': melt_wf["DFMQ"], 'S6ST_i': melt_wf["S6ST_i"], "ST":melt_wf["ST"],"CO2":melt_wf["CO2"],"H2OT":melt_wf["H2OT"],"HT":melt_wf["HT"],"CT":melt_wf["CT"], "XT":melt_wf["XT"]} # to work out P_sat
    melt_wf2 = {'SiO2': melt_wf["SiO2"], 'TiO2': melt_wf["TiO2"], 'Al2O3': melt_wf["Al2O3"], 'FeOT': melt_wf["FeOT"], 'Fe2O3T': melt_wf["Fe2O3T"], 'FeO': melt_wf["FeO"], 'Fe2O3': melt_wf["Fe2O3"], 'MgO': melt_wf["MgO"], 'MnO': melt_wf["MnO"], 'CaO': melt_wf["CaO"], 'Na2O': melt_wf["Na2O"], 'K2O': melt_wf["K2O"], 'P2O5': melt_wf["P2O5"], 'logfO2_i': melt_wf["logfO2_i"], 'Fe3FeT_i': melt_wf["Fe3FeT_i"], 'DNNO': melt_wf["DNNO"], 'DFMQ': melt_wf["DFMQ"], 'S6ST_i': melt_wf["S6ST_i"], "ST":melt_wf["ST"],"CO2":melt_wf["CO2"],"H2OT":melt_wf["H2OT"],"HT":melt_wf["HT"],"CT":melt_wf["CT"], "XT":melt_wf["XT"]} # to work out sulfur saturation
    
    def Pdiff(guess,melt_wf,species,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(PT,melt_wf,species,models))
        return difference

    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    melt_wf1["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
    melt_wf2["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
    #wt_C, wt_O, wt_H, wt_S, wt_Fe, wt_g, Wt = bulk_composition(run,PT,melt_wf1,setup,species,models)
    #bulk_wf = {"H":wt_H,"C":wt_C,"S":wt_S}
    ms_conc = eq.melt_speciation(PT,melt_wf1,species,models,nr_step,nr_tol)
    ms_frac = melt_species_ratios(ms_conc,species)
    melt_wf1["H2OT"] = ms_conc["wm_H2O"]
    melt_wf2["H2OT"] = ms_conc["wm_H2O"]
    melt_wf1["CO2"] = ms_conc["wm_CO2"]
    melt_wf2["CO2"] = ms_conc["wm_CO2"]
    melt_wf1["S2-"] = ms_conc["wm_S2m"]
    melt_wf2["S2-"] = ms_conc["wm_S2m"]
    melt_wf1["ST"] = ST
    melt_wf2["ST"] = ST
    if models.loc["sulfur_saturation","option"] == "yes": # must incorporate H2S concentration into S2- for SCSS
        sulfsat = sulfur_saturation(PT,melt_wf2,species,models)
        melt_wf1["ST"] = sulfsat["ST"]/1000000.
        melt_wf2["ST"] = ST
    delta1 = Pdiff(guess0,melt_wf1,species,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,melt_wf1,species,models)
        guess0 = mg.p_tot(PT,melt_wf1,species,models)
        guess0 = float(guess0)
        PT["P"] = guess0
        melt_wf1["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
        melt_wf2["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
        ms_conc = eq.melt_speciation(PT,melt_wf1,species,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc,species)
        melt_wf1["H2OT"] = ms_conc["wm_H2O"]
        melt_wf2["H2OT"] = ms_conc["wm_H2O"]
        melt_wf1["CO2"] = ms_conc["wm_CO2"]
        melt_wf2["CO2"] = ms_conc["wm_CO2"]
        melt_wf1["S2-"] = ms_conc["wm_S2m"]
        melt_wf2["S2-"] = ms_conc["wm_S2m"]
        if models.loc["sulfur_saturation","option"] == "yes":
            sulfsat = sulfur_saturation(PT,melt_wf2,species,models)
            melt_wf1["ST"] = sulfsat["ST"]/1000000.
            melt_wf2["ST"] = ST
    else:
        P_sat = guess0
        ms_conc = eq.melt_speciation(PT,melt_wf1,species,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc,species)
        
    melt_wf["ST"] = ST
    return P_sat, ms_conc, ms_frac

# for a given melt composition, calcualte the saturation pressure
def P_sat_H2O_CO2(PT,melt_wf,species,models,Ptol,nr_step,nr_tol): # Pvsat with just H2O and CO2 in vapour
    
    def p_tot_H2O_CO2(PT,melt_wf,species,models):
        value = mg.p_H2O(PT,melt_wf,species,models) + mg.p_CO2(PT,melt_wf,species,models)  
        return value
    
    def Pdiff(guess,melt_wf,species,models):
        PT["P"] = guess
        difference = abs(guess - p_tot_H2O_CO2(PT,melt_wf,species,models))
        return difference
    
    if melt_wf['CO2'] == 0 and melt_wf['H2OT'] == 0:
        P_sat = ""
        result = {"xg_H2O":0., "xg_CO2":0., "f_H2O":0., "f_CO2":0., "p_H2O":0., "p_CO2":0., "wm_H2O":0., "wm_CO2":0.}
    
    else:
        guess0 = 40000. # initial guess for pressure
        PT["P"] = guess0
        delta1 = Pdiff(guess0,melt_wf,species,models)
        while delta1 > Ptol :
            delta1 = Pdiff(guess0,melt_wf,species,models)
            guess0 = p_tot_H2O_CO2(PT,melt_wf,species,models)
            guess0 = float(guess0)
            PT["P"] = guess0
        else:
            P_sat = guess0
            xg_H2O_ = mg.xg_H2O(PT,melt_wf,species,models)
            xg_CO2_ = mg.xg_CO2(PT,melt_wf,species,models)
            p_H2O_ = mg.p_H2O(PT,melt_wf,species,models)
            p_CO2_ = mg.p_CO2(PT,melt_wf,species,models)
            f_H2O_ = mg.f_H2O(PT,melt_wf,species,models)
            f_CO2_ = mg.f_CO2(PT,melt_wf,species,models)
            xm_H2O_ = (mdv.C_H2O(PT,melt_wf,species,models)*f_H2O_)**0.5
            xm_CO2_ = mdv.C_CO3(PT,melt_wf,species,models)*f_CO2_
            M_m_ = mg.M_m_SO(melt_wf,species)
            Xm_t = xm_CO2_*species.loc["CO2","M"] + xm_H2O_*species.loc["H2O","M"] + (1.0-xm_CO2_-xm_H2O_)*M_m_
            wm_H2O_ = (xm_H2O_*species.loc["H2O","M"])/Xm_t
            wm_CO2_ = (xm_CO2_*species.loc["CO2","M"])/Xm_t
    
        result = {"xg_H2O":xg_H2O_, "xg_CO2":xg_CO2_, "f_H2O":f_H2O_, "f_CO2":f_CO2_, "p_H2O":p_H2O_, "p_CO2":p_CO2_, "wm_H2O":wm_H2O_, "wm_CO2":wm_CO2_}
    
    return P_sat, result

# for a given fS2 and fO2, calculate Psat
def P_sat_fO2_fS2(PT,melt_wf,species,models,Ptol):
    
    def Pdiff(guess,melt_wf,species,models):
        PT["P"] = guess
        result_fO2fS2 = eq.p_tot_fO2_fS2(PT,melt_wf,species,models)
        difference = abs(guess - result_fO2fS2["P_tot"])
        return difference

    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    delta1 = Pdiff(guess0,melt_wf,species,models)
    
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,melt_wf,species,models)
        result_fO2fS2 = eq.p_tot_fO2_fS2(PT,melt_wf,species,models)
        guess0 = result_fO2fS2["P_tot"]
        guess0 = float(guess0)
        PT["P"] = guess0
    else:
        result = eq.p_tot_fO2_fS2(PT,melt_wf,species,models)
    
    return results
        
        
        
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
    S6ST_ = mg.S6ST(PT,melt_wf,species,models)
    #SCSS_,sulfide_sat,SCAS_,sulfate_sat = sulfur_saturation(wm_ST/100.0,S6ST_)
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
        wt_C = species.loc["C","M"]*((wt_g*(((mg.xg_CO2(PT,melt_wf,species,models)+
                                      mg.xg_CO(PT,melt_wf,species,models)+
                                      mg.xg_CH4(PT,melt_wf,species,models)+
                                      mg.xg_OCS(PT,melt_wf,species,models))
                                     /mg.Xg_tot(PT,melt_wf,species,models)) 
                                    - (wm_CO2/species.loc["CO2","M"] - wm_CO/species.loc["CO","M"] - 
                                       wm_CH4/species.loc["CH4","M"]))) + (wm_CO2/species.loc["CO2","M"] + 
                                        wm_CO/species.loc["CO","M"] + wm_CH4/species.loc["CH4","M"]))

    wt_H = 2.*species.loc["H","M"]*((wt_g*(((mg.xg_H2O(PT,melt_wf,species,models)+mg.xg_H2(PT,melt_wf,species,models)+2.*mg.xg_CH4(PT,melt_wf,species,models)+mg.xg_H2S(PT,melt_wf,species,models))/mg.Xg_tot(PT,melt_wf,species,models)) - (wm_H2O/species.loc["H2O","M"] - wm_H2/species.loc["H2","M"] - wm_H2S/species.loc["H2S","M"] - (2.*wm_CH4)/species.loc["CH4","M"]))) + (wm_H2O/species.loc["H2O","M"] + wm_H2/species.loc["H2","M"] + wm_H2S/species.loc["H2S","M"] + (2.*wm_CH4)/species.loc["CH4","M"]))
    
    wt_S = species.loc["S","M"]*((wt_g*(((mg.xg_SO2(PT,melt_wf,species,models)
                                      +2.0*mg.xg_S2(PT,melt_wf,species,models)
                                      +mg.xg_H2S(PT,melt_wf,species,models)
                                      +mg.xg_OCS(PT,melt_wf,species,models))
                                     /mg.Xg_tot(PT,melt_wf,species,models)) 
                                    - ((wm_S2m+wm_S6p)/species.loc["S","M"] - (wm_H2S/species.loc["H2S","M"])))) 
                              + ((wm_S2m+wm_S6p)/species.loc["S","M"] + (wm_H2S/species.loc["H2S","M"])))
    
    species_X = models.loc["species X","option"]
    wt_X = species.loc[species_X,"M"]*((wt_g*(((mg.xg_X(PT,melt_wf,species,models))
                                     /mg.Xg_tot(PT,melt_wf,species,models)) 
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
        wt_Fe = (1.-wt_g)*(((1.0-total_dissolved_volatiles)*mg.wm_Fe_nv(melt_wf,species))/100.0) # wt fraction of Fe

    wt_O = species.loc["O","M"]*((wt_g*(((2.0*mg.xg_CO2(PT,melt_wf,species,models) + mg.xg_CO(PT,melt_wf,species,models) + 2.0*mg.xg_O2(PT,melt_wf,species,models) + mg.xg_H2O(PT,melt_wf,species,models) + 2.0*mg.xg_SO2(PT,melt_wf,species,models) + mg.xg_OCS(PT,melt_wf,species,models))/mg.Xg_tot(PT,melt_wf,species,models)) - (wm_H2O/species.loc["H2O","M"]) - ((2.0*wm_CO2)/species.loc["CO2","M"]) - (3.0*wm_S6p/species.loc["S","M"] - (wm_CO/species.loc["CO","M"])))) + (wm_H2O/species.loc["H2O","M"]) + ((2.0*wm_CO2)/species.loc["CO2","M"]) + (3.0*wm_S6p/species.loc["S","M"]) + (wm_CO/species.loc["CO","M"]) + (wt_Fe/species.loc["Fe","M"])*((1.5*Fe3Fe2_+1.0)/(Fe3Fe2_+1.0)))
    
    result = {"wt_C":wt_C, "wt_O":wt_O, "wt_H":wt_H, "wt_S":wt_S, "wt_X":wt_X, "wt_Fe":wt_Fe, "wt_g":wt_g, "Wt":Wt}
    
    return result

# calculate weight fraction of elements in the system when adding gas into a melt
def new_bulk_regas_open(PT,melt_wf,bulk_wf,gas_mf,dwtg,species,models):
    me = mg.melt_elements(PT,melt_wf,bulk_wf,gas_mf,species,models)
    mg = mg.gas_elements(gas_mf,species)
    wt_C = (1.-dwtg)*me["wm_C"] + dwtg*ge["wg_C"]
    wt_H = (1.-dwtg)*me["wm_H"] + dwtg*ge["wg_H"]
    wt_O = (1.-dwtg)*me["wm_O"] + dwtg*ge["wg_O"]
    wt_S = (1.-dwtg)*me["wm_S"] + dwtg*ge["wg_S"]
    wt_X = (1.-dwtg)*me["wm_X"] + dwtg*ge["wg_X"]
    wt_Fe = (1.-dwtg)*me["wm_Fe"]
    Wt = bulk_wf['Wt']*(1. + dwtg)
    result = {"wt_C":wt_C, "wt_H":wt_H, "wt_S":wt_S, "wt_Fe":wt_Fe, "wt_O":wt_O, "wt_X":wt_X, "Wt":Wt}
    return result

# calculate mole fraction of species in the cumulated gas from open-system degassing
def gas_comp_all_open(xg,xg_all,models,species):
    species_X = models.loc["species X","option"]
    wt_g_before = xg_all["wt_g"]
    wt_g_inst = xg["wt_g"]*(1.-wt_g_before)
    wt_g_new = wt_g_before + wt_g_inst
    f_before = wt_g_before/wt_g_new
    f_inst = wt_g_inst/wt_g_new
    wf_before = mg.gas_wf(xg_all,species,models)
    wf_inst = mg.gas_wf(xg,species,models)
    gas_wf_all_new = {"O2":wf_inst["wg_O2"]*f_inst+wf_before["wg_O2"]*f_before}
    gas_wf_all_new["CO"] = wf_inst["wg_CO"]*f_inst+wf_before["wg_CO"]*f_before
    gas_wf_all_new["S2"] = wf_inst["wg_S2"]*f_inst+wf_before["wg_S2"]*f_before
    gas_wf_all_new["CO2"] = wf_inst["wg_CO2"]*f_inst+wf_before["wg_CO2"]*f_before
    gas_wf_all_new["H2O"] = wf_inst["wg_H2O"]*f_inst+wf_before["wg_H2O"]*f_before
    gas_wf_all_new["H2"] = wf_inst["wg_H2"]*f_inst+wf_before["wg_H2"]*f_before
    gas_wf_all_new["CH4"] = wf_inst["wg_CH4"]*f_inst+wf_before["wg_CH4"]*f_before
    gas_wf_all_new["SO2"] = wf_inst["wg_SO2"]*f_inst+wf_before["wg_SO2"]*f_before
    gas_wf_all_new["H2S"] = wf_inst["wg_H2S"]*f_inst+wf_before["wg_H2S"]*f_before
    gas_wf_all_new["OCS"] = wf_inst["wg_OCS"]*f_inst+wf_before["wg_OCS"]*f_before
    gas_wf_all_new["X"] = wf_inst["wg_X"]*f_inst+wf_before["wg_X"]*f_before
    tot = (gas_wf_all_new["X"]/species.loc[species_X,"M"]) + (gas_wf_all_new["OCS"]/species.loc["OCS","M"]) + (gas_wf_all_new["H2S"]/species.loc["H2S","M"]) + (gas_wf_all_new["SO2"]/species.loc["SO2","M"]) + (gas_wf_all_new["CH4"]/species.loc["CH4","M"]) + (gas_wf_all_new["H2"]/species.loc["H2","M"]) + (gas_wf_all_new["H2O"]/species.loc["H2O","M"]) + (gas_wf_all_new["CO2"]/species.loc["CO2","M"]) + (gas_wf_all_new["O2"]/species.loc["O2","M"]) + (gas_wf_all_new["CO"]/species.loc["CO","M"]) + (gas_wf_all_new["S2"]/species.loc["S2","M"])
    gas_mf_all_new = {"O2":(gas_wf_all_new["O2"]/species.loc["O2","M"])/tot}
    gas_mf_all_new["CO"] = (gas_wf_all_new["CO"]/species.loc["CO","M"])/tot
    gas_mf_all_new["S2"] = (gas_wf_all_new["S2"]/species.loc["S2","M"])/tot
    gas_mf_all_new["CO2"] = (gas_wf_all_new["CO2"]/species.loc["CO2","M"])/tot
    gas_mf_all_new["H2O"] = (gas_wf_all_new["H2O"]/species.loc["H2O","M"])/tot
    gas_mf_all_new["H2"] = (gas_wf_all_new["H2"]/species.loc["H2","M"])/tot
    gas_mf_all_new["CH4"] = (gas_wf_all_new["CH4"]/species.loc["CH4","M"])/tot
    gas_mf_all_new["SO2"] = (gas_wf_all_new["SO2"]/species.loc["SO2","M"])/tot
    gas_mf_all_new["H2S"] = (gas_wf_all_new["H2S"]/species.loc["H2S","M"])/tot
    gas_mf_all_new["OCS"] = (gas_wf_all_new["OCS"]/species.loc["OCS","M"])/tot
    gas_mf_all_new["X"] = (gas_wf_all_new["X"]/species.loc[species_X,"M"])/tot
    gas_mf_all_new["Xg_t"] = gas_mf_all_new["CO2"]*species.loc["CO2","M"] + gas_mf_all_new["CO"]*species.loc["CO","M"] + gas_mf_all_new["O2"]*species.loc["O2","M"] + gas_mf_all_new["H2O"]*species.loc["H2O","M"] + gas_mf_all_new["H2"]*species.loc["H2","M"] + gas_mf_all_new["CH4"]*species.loc["CH4","M"] + gas_mf_all_new["SO2"]*species.loc["SO2","M"] + gas_mf_all_new["S2"]*species.loc["S2","M"] + gas_mf_all_new["H2S"]*species.loc["H2S","M"] + gas_mf_all_new["OCS"]*species.loc["OCS","M"] + gas_mf_all_new["X"]*species.loc[species_X,"M"]
    gas_mf_all_new["wt_g"] = wt_g_new
    return gas_mf_all_new

# calculate ratios of volatile species in the melt
def melt_species_ratios(conc,species):

    M_H = species.loc['H','M']
    M_S = species.loc['S','M']
    M_C = species.loc['C','M']
    M_CO = species.loc['CO','M']
    M_H2O = species.loc['H2O','M']
    M_H2 = species.loc['H2','M']
    M_CO2 = species.loc['CO2','M']
    M_CH4 = species.loc['CH4','M']
    M_H2S = species.loc['H2S','M']

    wt_H = conc["wm_H2"] + ((2.*(conc["wm_H2O"]/M_H2O))*M_H) + ((4.*(conc["wm_CH4"]/M_CH4))*M_H) + ((2.*(conc["wm_H2S"]/M_H2S))*M_H)
    if wt_H > 0.: # contains H
        H2_HT = conc["wm_H2"]/wt_H
        H2O_HT = ((2.*(conc["wm_H2O"]/M_H2O))*M_H)/wt_H
        CH4_HT = ((4.*(conc["wm_CH4"]/M_CH4))*M_H)/wt_H
        H2S_HT = ((2.*(conc["wm_H2S"]/M_H2S))*M_H)/wt_H
    
    wt_C = ((conc["wm_CO"]/M_CO)*M_C) + ((conc["wm_CO2"]/M_CO2)*M_C) + ((conc["wm_CH4"]/M_CH4)*M_C)
    if wt_C > 0.: # contains C
        CO_CT = ((conc["wm_CO"]/M_CO)*M_C)/wt_C
        CO2_CT = ((conc["wm_CO2"]/M_CO2)*M_C)/wt_C
        CH4_CT = ((conc["wm_CH4"]/M_CH4)*M_C)/wt_C

    wt_S = conc["wm_S2m"] + conc["wm_S6p"] + (M_S*(conc["wm_H2S"]/M_H2S))
    if wt_S > 0.: # contains S
        S2m_ST = conc["wm_S2m"]/wt_S
        S6p_ST = conc["wm_S6p"]/wt_S
        H2S_ST = (M_S*(conc["wm_H2S"]/M_H2S))/wt_S

    frac = {"H2O_HT":H2O_HT, "H2_HT":H2_HT, "CH4_HT":CH4_HT, "CO2_CT":CO2_CT, "CO_CT":CO_CT, "CH4_CT":CH4_CT, "S6p_ST":S6p_ST, "S2m_ST":S2m_ST, "H2S_ST":H2S_ST, "H2S_HT":H2S_HT}
    
    return frac


#########################
### sulfur satuation ###
#########################

# check solid/immiscible liquid sulfur saturation
def sulfur_saturation(PT,melt_wf,species,models): # melt weight fraction of ST and S6/ST
    wmST = melt_wf['ST']
    S6T = mg.S6ST(PT,melt_wf,species,models)
    wmS2 = wmST*100.0*10000.0*(1.0-S6T)
    wmS6 = wmST*100.0*10000.0*S6T
    SCSS_ = mdv.SCSS(PT,melt_wf,species,models)
    SCAS_ = mdv.SCAS(PT,melt_wf,species,models)
    StCSS = SCSS_/(1.-S6T)
    StCAS = SCAS_/S6T
    if wmS2 < SCSS_ and wmS6 < SCAS_:
        sulfide_sat = "no"
        sulfate_sat = "no"
        ST = wmST*1000000.
    elif wmS2 >= SCSS_ and wmS6 >= SCAS_:
        sulfide_sat = "yes"
        sulfate_sat = "yes"
        ST = min(StCSS,StCAS)
    elif wmS2 >= SCSS_ and wmS6 < SCAS_:
        sulfide_sat = "yes"
        sulfate_sat = "no"
        ST = StCSS
    elif wmS2 < SCSS_ and wmS6 >= SCAS_:
        sulfide_sat = "no"
        sulfate_sat = "yes"
        ST = StCAS
    else:
        sulfide_sat = "nan"
        sulfate_sat = "nan"
        ST = wmST*1000000.
    result = {"SCSS":SCSS_, "sulfide_sat":sulfide_sat, "SCAS":SCAS_, "sulfate_sat":sulfate_sat, "ST":ST}
    return result

# fO2 and P of v+sulf+anh saturation
def fO2_P_VSA(PT,melt_wf,species,models,nr_step,nr_tol,Ptol):

    def Pdiff(guess,melt_wf,species,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(PT,melt_wf,species,models))
        return difference

    def fO2_S(run,PT,melt_wf,setup,species,models):
        SCSS_ = mdv.SCSS(PT,melt_wf,species,models)/1000000.
        SCAS_ = mdv.SCAS(PT,melt_wf,species,models)/1000000.
        CSO4 = mdv.C_SO4(PT,melt_wf,species,models)/1000000.
        CS = mdv.C_S(PT,melt_wf,species,models)/1000000.
        
        if models.loc["H2S_m","option"] == "no":
            W = CSO4/CS
        elif models.loc["H2S_m","option"] == "yes":
            CH2S = mdv.C_H2S(PT,melt_wf,species,models)/1000000.
            KHS = mdv.KHOSg(PT,models)
            CH2OT = mdv.C_H2O(PT,melt_wf,species,models)
            xmH2O = mg.xm_H2OT_so(melt_wf,species)
            W = (CSO4/((KHS*CH2S*(xmH2O**2./CH2OT)) + CS))
        fO2 = ((1./W)*(SCAS_/SCSS_))**0.5
        ST = SCSS_ + SCAS_
        return fO2, ST
    
    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    fO2_, ST_ = fO2_S(PT,melt_wf,species,models)
    melt_wf["ST"] = ST_
    melt_wf["Fe3FeT"] = mdv.fO22Fe3FeT(fO2_,PT,melt_wf,species,models)
    ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
    ms_frac = melt_species_ratios(ms_conc,species)
    melt_wf["H2OT"] = ms_cont["wm_H2O"]
    melt_wf["CO2"] = ms_cont["wm_CO2"]
    melt_wf["S2-"] = ms_cont["wm_S2m"]
    delta1 = Pdiff(guess0,melt_wf,species,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,melt_wf,species,models)
        guess0 = mg.p_tot(PT,melt_wf,species,models)
        guess0 = float(guess0)
        PT["P"] = guess0
        fO2_, ST_ = fO2_S(PT,melt_wf,species,models)
        melt_wf["ST"] = ST_
        melt_wf["Fe3FeT"] = mdv.fO22Fe3FeT(fO2_,PT,melt_wf,species,models)
        ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc,species)        
        melt_wf["H2OT"] = ms_conc["wm_H2O"]
        melt_wf["CO2"] = ms_conc["wm_CO2"]
        melt_wf["S2-"] = ms_conc["wm_S2m"]
    else:
        P_sat = guess0
        ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc,species)
        fO2_, ST_ = fO2_S(PT,melt_wf,species,models)
        Fe3_FT = mdv.fO22Fe3FeT(fO2_,PT,melt_wf,species,models)

    ms_frac["Fe3_FT"] = Fe3_FT 
              
    return P_sat, ms_conc, ms_frac

def P_VSA(PT,melt_wf,species,models,nr_step,nr_tol,Ptol):
    
    #P_sat, pVSA_conc, pVSA_frac = fO2_P_VSA(PT,melt_wf,species,models,nr_step,nr_tol,Ptol)
    
    def Pdiff(guess,melt_wf,species,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(PT,melt_wf,species,models))
        return difference
    
    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
    if melt_wf["Fe3FeT"] < Fe3_FT:
        SCSS_ = mdv.SCSS(PT,melt_wf,species,models)/1000000.
        S6T = mg.S6ST(PT,melt_wf,species,models)
        ST_ = SCSS_/(1.-S6T)
    else:
        SCAS_ = mdv.SCAS(PT,melt_wf,species,models)/1000000.
        S6T = mg.S6ST(PT,melt_wf,species,models)
        ST_ = SCAS_/S6T
    melt_wf["ST"] = ST_
    ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
    ms_frac = melt_species_ratios(ms_conc,species)
    melt_wf["H2OT"] = ms_conc["wm_H2O"]
    melt_wf["CO2"] = ms_conc["wm_CO2"]
    melt_wf["S2-"] = ms_conc["wm_S2m"]
    delta1 = Pdiff(guess0,melt_wf,species,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,melt_wf,species,models)
        guess0 = mg.p_tot(PT,melt_wf,species,models)
        guess0 = float(guess0)
        PT["P"] = guess0
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
        if melt_wf["Fe3FeT"] < Fe3_FT:
            SCSS_ = mdv.SCSS(PT,melt_wf,species,models)/1000000.
            S6T = mg.S6ST(PT,melt_wf,species,models)
            ST_ = SCSS_/(1.-S6T)
        else:
            SCAS_ = mdv.SCAS(PT,melt_wf,species,models)/1000000.
            S6T = mg.S6ST(PT,melt_wf,species,models)
            ST_ = SCAS_/S6T
        melt_wf["ST"] = ST_
        ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc,species)    
        melt_wf["H2OT"] = ms_conc["wm_H2O"]
        melt_wf["CO2"] = ms_conc["wm_CO2"]
        melt_wf["S2-"] = ms_conc["wm_S2m"]
    else:
        P_sat = guess0
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
        ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc,species)
    
    return P_sat,ms_conc,ms_frac

# estimated fO2 from sulfur content given you are sulfide saturated and at Pvsat
def P_sat_sulf_anh(PT,melt_wf,species,models,Ptol,nr_step,nr_tol):
    ST = melt_wf["ST"]
    H2OT = melt_wf["H2OT"]
    CO2 = melt_wf["CO2"]
    HT = melt_wf["HT"]
    CT = melt_wf["CT"]
 
    def Pdiff(guess,melt_wf,species,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(PT,melt_wf,species,models))
        return difference
    
    def S62_2_Fe3T(PT,melt_wf,species,models,sat):
        ST = melt_wf["ST"]
        if sat == "sulf":
            Ssat = mdv.SCSS(PT,melt_wf,species,models)/1000000.
            S2m = Ssat
            S6p = ST - S2m
        elif sat == "anh":
            Ssat = mdv.SCAS(PT,melt_wf,species,models)/1000000.
            S6p = Ssat
            S2m = ST - S6p
        S62 = S6p/S2m
        S6T = S6p/ST
        if S62 < 0.:
            return "not possible","","","",Ssat
        else:
            fO2 = mg.S6S2_2_fO2(S62,melt_wf,PT,species,models)
            Fe3T = mdv.fO22Fe3FeT(fO2,PT,melt_wf,species,models)
            DFMQ = mg.fO22Dbuffer(PT,fO2,"FMQ",models)
            return Fe3T,fO2,S6T,DFMQ,Ssat
        
    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    
    # assume it is sulfide saturated
    melt_wf["Fe3FeT"] = 0.1
    Fe3T_sulf,fO2_sulf,S6T_sulf,DFMQ_sulf,SCSS_ = S62_2_Fe3T(PT,melt_wf,species,models,"sulf")
    if Fe3T_sulf == "not possible":
        P_sat_sulf = ""
        Fe3T_sulf = ""
        sulfide_sat = "no"
    else: 
        melt_wf["Fe3FeT"] = Fe3T_sulf
        ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc,species)
        melt_wf["H2OT"] = ms_conc["wm_H2O"]
        melt_wf["CO2"] = ms_conc["wm_CO2"]
        melt_wf["S2-"] = ms_conc["wm_S2m"]
        melt_wf["ST"] = ST
        delta1 = Pdiff(guess0,melt_wf,species,models)
        while delta1 > Ptol :
            delta1 = Pdiff(guess0,melt_wf,species,models)
            guess0 = mg.p_tot(PT,melt_wf,species,models)
            guess0 = float(guess0)
            PT["P"] = guess0
            Fe3T_sulf,fO2_sulf,S6T_sulf,DFMQ_sulf,SCSS_ = S62_2_Fe3T(PT,melt_wf,species,models,"sulf")
            if Fe3T_sulf == "not possible":
                P_sat_sulf = ""
                Fe3T_sulf = ""
                sulfide_sat = "no"
            else: 
                melt_wf["Fe3FeT"] = Fe3T_sulf
                ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
                ms_frac = melt_species_ratios(ms_conc,species)
                melt_wf["H2OT"] = ms_conc["wm_H2O"]
                melt_wf["CO2"] = ms_conc["wm_CO2"]
                melt_wf["S2-"] = ms_conc["wm_S2m"]
        else:
            P_sat_sulf = guess0
            Fe3T_sulf,fO2_sulf,S6T_sulf,DFMQ_sulf,SCSS_ = S62_2_Fe3T(PT,melt_wf,species,models,"sulf")
            sulfide_sat = "yes"
           
    # assume it is anhydrite saturated
    Fe3T_anh,fO2_anh,S6T_anh,DFMQ_anh,SCAS_ = S62_2_Fe3T(PT,melt_wf,species,models,"anh")
    if Fe3T_anh == "not possible":
        P_sat_anh = ""
        Fe3T_anh = ""
        anhydrite_sat = "no"
    else:  
        melt_wf["Fe3FeT"] = Fe3T_anh
        ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc,species)
        melt_wf["H2OT"] = ms_conc["wm_H2O"]
        melt_wf["CO2"] = ms_conc["wm_CO2"]
        melt_wf["S2-"] = ms_conc["wm_S2m"]
        melt_wf["ST"] = ST
        delta1 = Pdiff(guess0,melt_wf,species,models)
        while delta1 > Ptol :
            delta1 = Pdiff(guess0,melt_wf,species,models)
            guess0 = mg.p_tot(PT,melt_wf,species,models)
            guess0 = float(guess0)
            PT["P"] = guess0
            Fe3T_anh,fO2_anh,S6T_anh,DFMQ_anh,SCAS_ = S62_2_Fe3T(PT,melt_wf,species,models,"anh")
            if Fe3T_anh == "not possible":
                P_sat_anh = ""
                Fe3T_anh = ""
                anhydrite_sat = "no"
            else:  
                melt_wf["Fe3FeT"] = Fe3T_anh
                ms_conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
                ms_frac = melt_species_ratios(ms_conc,species)
                melt_wf["H2OT"] = ms_conc["wm_H2O"]
                melt_wf["CO2"] = ms_conc["wm_CO2"]
                melt_wf["S2-"] = ms_conc["wm_S2m"]
        else:
            P_sat_anh = guess0
            Fe3T_anh,fO2_anh,S6T_anh,DFMQ_anh,SCAS_ = S62_2_Fe3T(PT,melt_wf,species,models,"anh")
            anhydrite_sat = "yes"
    
    result = {"P_sat_sulf":P_sat_sulf,"P_sat_anh":P_sat_anh,"SCAS":SCAS_*1000000.,"SCSS":SCSS_*1000000.,"sulf_sat":sulfide_sat,"DFMQ_sulf":DFMQ_sulf,"fO2_sulf":fO2_sulf,"Fe3T_sulf":Fe3T_sulf,"S6T_sulf":S6T_sulf,"anh_sat":anhydrite_sat,"DFMQ_anh":DFMQ_anh,"fO2_anh":fO2_anh,"Fe3T_anh":Fe3T_anh,"S6T_anh":S6T_anh}          

    return result

##########################
### graphite satuation ###
##########################

# check graphite saturation
def graphite_saturation(PT,melt_wf,species,models): # needs finishing
    K1 = mg.f_CO2(PT,melt_wf,species,models)/mdv.f_O2(PT,melt_wf,species,models)
    K2 = mdv.KCOs(PT,models) # K for graphite saturation
    if K1 < K2:
        graphite_sat = "no"
    else:
        graphite_sat = "yes"
    fCO2_ = K2*mdv.f_O2(PT,melt_wf,species,models)
    xmCO2 = fCO2_*mdv.C_CO3(PT,melt_wf,species,models)
    return graphite_sat

                
######################################                
### fO2 range from sulfur content ###
######################################
                
def fO2_range_from_S(PT,melt_wf,species,models):
    SCSS_ = mdv.SCSS(PT,melt_wf,species,models)
    SCAS_ = mdv.SCAS(PT,melt_wf,species,models)

    if melt_wf["ST"]*1000000. > SCSS_:
        sulfide_sat = "possible"
        S6ST_1 = (melt_wf["ST"]*1000000. - SCSS_)/(melt_wf["ST"]*1000000.)
        S6S2_1 = mg.overtotal2ratio(S6ST_1)
        fO2_1 = mg.S6S2_2_fO2(S6S2_1,melt_wf,PT,species,models)
        Fe3FeT_1 = mdv.fO22Fe3FeT(fO2_1,PT,melt_wf,species,models)
        Fe3Fe2_1 = mg.overtotal2ratio(Fe3FeT_1)
        DFMQ_1 = mg.fO22Dbuffer(PT,fO2_1,"FMQ",models)
    else:
        sulfide_sat = "no"
        S6ST_1 = ""
        S6S2_1 = ""
        Fe3Fe2_1 = ""
        Fe3FeT_1 = ""
        fO2_1 = ""
        DFMQ_1 = ""
                   
    if melt_wf["ST"]*1000000. > SCAS_:
        anhydrite_sat = "possible"
        S6ST_2 = SCAS_/(melt_wf["ST"]*1000000.)
        S6S2_2 = mg.overtotal2ratio(S6ST_2)
        fO2_2 = mg.S6S2_2_fO2(S6S2_2,melt_wf,PT,species,models)
        Fe3FeT_2 = mdv.fO22Fe3FeT(fO2_2,PT,melt_wf,species,models)
        Fe3Fe2_2 = mg.overtotal2ratio(Fe3FeT_2)
        DFMQ_2 = mg.fO22Dbuffer(PT,fO2_2,"FMQ",models)
    else:
        anhydrite_sat = "no"
        S6ST_2 = ""
        S6S2_2 = ""
        Fe3Fe2_2 = ""
        Fe3FeT_2 = ""
        fO2_2 = ""
        DFMQ_2 = ""
                                 
    result = {"SCAS":SCAS_, "SCSS":SCSS_, "sulf_sat":sulfide_sat, "DFMQ_sulf":DFMQ_1, "fO2_sulf":fO2_1, "Fe3FeT_sulf":Fe3FeT_1, "S6ST_sulf":S6ST_1, "anh_sat":anhydrite_sat, "DFMQ_anh":DFMQ_2, "fO2_anh":fO2_2, "Fe3FeT_anh":Fe3FeT_2, "S6ST_anh":S6ST_2}

    return result
    
                             
#################################                             
### Mass, volume, and density ###
#################################                             
                             
def mass_vol_rho(PT,melt_wf,gas_mf,bulk_wf,species,models):
    gas_m = gas_mf['wt_g']*bulk_wf['Wt']
    gas_v = mg.gas_volume(PT,gas_mf,bulk_wf,species,models)
    gas_rho = mg.gas_density(PT,gas_mf,bulk_wf,species,models)
    melt_m = (1.-gas_mf['wt_g'])*bulk_wf['Wt']
    melt_v = mg.melt_volume(PT,melt_wf,bulk_wf,gas_mf,species)
    melt_rho = mg.melt_density(PT,melt_wf,species)
    tot_m = bulk_wf['Wt']
    tot_v = gas_v + melt_v
    tot_rho = mg.system_density(PT,melt_wf,gas_mf,bulk_wf,species,models)
    result = {"tot_m":tot_m, "tot_v":tot_v, "tot_rho":tot_rho, "melt_m":melt_m, "melt_v":melt_v, "melt_rho":melt_rho, "gas_m":gas_m, "gas_v":gas_v, "gas_rho":gas_rho}
    return result


######################################################
### mole fraction of elements in different species ###
######################################################

def mf_S_species(melt_wf,gas_mf,species):
    # weight of S in each sulfur-bearing species
    W_S2m = melt_wf["ST"]*(1.-gas_mf["wt_g"])*(1.-melt_wf["S6ST"])
    W_SO4 = melt_wf["ST"]*(1.-gas_mf["wt_g"])*melt_wf["S6ST"]
    W_H2S = ((gas_mf["H2S"]*species.loc["S","M"])/(gas_mf["Xg_t"]))*gas_mf["wt_g"]
    W_SO2 = ((gas_mf["SO2"]*species.loc["S","M"])/(gas_mf["Xg_t"]))*gas_mf["wt_g"]
    W_S2 = ((gas_mf["S2"]*species.loc["S2","M"])/(gas_mf["Xg_t"]))*gas_mf["wt_g"]
    W_OCS = ((gas_mf["OCS"]*species.loc["S","M"])/(gas_mf["Xg_t"]))*gas_mf["wt_g"]
    W_total = W_S2m + W_SO4 + W_H2S + W_SO2 + W_S2 + W_OCS
    # weight and mole fraction of S in each sulfur-bearing species compared to total S
    w_S2m = W_S2m/W_total
    w_SO4 = W_SO4/W_total
    w_SO2 = W_SO2/W_total
    w_H2S = W_H2S/W_total
    w_S2 = W_S2/W_total
    w_OCS = W_OCS/W_total
    mf_S = {"S2-":w_S2m, "SO42-":w_SO4, "SO2":w_SO2, "H2S":w_H2S, "S2": w_S2, "OCS": w_OCS}
    return mf_S
        
##############################################
### fO2 of silm+sulfm+anh at given T and P ###
##############################################

def fO2_silm_sulf_anh(PT,melt_wf,species,models):
    S6 = mdv.SCAS(PT,melt_wf,species)
    S2 = mdv.SCSS(PT,melt_wf,species,models)
    fO2 = ((S6*mdv.C_S(PT,melt_wf,species,models))/(S2*mdv.C_SO4(PT,melt_wf,species,models)))**0.5
    DFMQ = mg.fO22Dbuffer(PT,fO2,"FMQ",models)
    wmST = S6+S2
    S6ST = S6/wmST
    result = {"fO2":fO2, "DFMQ":DFMQ, "wmST":wmST, "S6ST":S6ST, "S6":S6, "S2":S2}
    return result
   

##############################################
### S content at given T, P, fO2, C, and H ###
##############################################

def S_given_T_P_fO2_C_H(PT,melt_wf,species,models,nr_step,nr_tol): # no dissolved H2S  
    P = PT["P"]
    
    conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
    
    melt_wf["H2OT"] = conc["wm_H2O"]
    melt_wf["CO2"] = conc["wm_CO2"]
  
    xg_O2_ = mg.p_O2(PT,melt_wf,species,models)/P
    xg_CO2_ = mg.p_CO2(PT,melt_wf,species,models)/P
    xg_CO_ = mg.p_CO(PT,melt_wf,species,models)/P
    xg_H2_ = mg.p_H2(PT,melt_wf,species,models)/P
    xg_H2O_ = mg.p_H2O(PT,melt_wf,species,models)/P
    xg_CH4_ = mg.p_CH4(PT,melt_wf,species,models)/P
    
    K6_ = mdv.KOSg(PT,models)
    K7_ = mdv.KHOSg(PT,models)
    K8_ = mdv.C_S(PT,melt_wf,species,models)/1000000.0
    K9_ = (mdv.C_SO4(PT,melt_wf,species,models)/1000000.0)
    K10_ = mdv.KOCSg(PT,models)                                           
    y_S2_ = mdv.y_S2(PT,species,models)
    y_SO2_ = mdv.y_SO2(PT,species,models)
    y_H2S_ = mdv.y_H2S(PT,species,models)
    y_O2_ = mdv.y_O2(PT,species,models)
    y_H2O_ = mdv.y_H2O(PT,species,models)
    y_CO2_ = mdv.y_CO2(PT,species,models)
    y_CO_ = mdv.y_CO(PT,species,models)
    y_OCS_ = mdv.y_OCS(PT,species,models)
    M_S = species.loc['S','M']
    M_SO3 = species.loc['SO3','M']
    
    a = 1.
    b = ((K6_*(y_S2_*P)**0.5*(xg_O2_*y_O2_))/y_SO2_) + ((K7_*xg_H2O_*y_H2O_*y_S2_**0.5)/((xg_O2_*y_O2_)**0.5*y_H2S_)) + ((K6_*(y_S2_*P)**0.5*xg_O2_*y_O2_*(xg_CO_*y_CO_)**3.*P)/(K10_*(xg_CO2_*y_CO2_)**2.*y_OCS_))
    c = (xg_O2_ + xg_H2O_ + xg_H2_ + xg_CO_ + xg_CO2_ + xg_CH4_) - 1.
    x = (-b + (b**2 - 4.*a*c)**0.5)/(2.*a)
    xg_S2_ = x**2.
    xg_SO2_ = (K6_*(xg_S2_*P*y_S2_)**0.5*(xg_O2_*P*y_O2_))/(y_SO2_*P)
    wm_S_ = K8_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
    wm_S6p_ = (K9_*(y_S2_*xg_S2_*P)**0.5*(y_O2_*xg_O2_*P)**1.5)
    wm_SO3_ = wm_S6p_*(M_SO3/M_S)
    wm_ST_ = wm_S_ + wm_S6p_
    
    conc["wm_ST"] = wm_ST_
    conc["wm_S2m"] = wm_S_
    conc["wm_S6p"] = wm_S6p_
    conc["wm_SO3"] = wm_SO3_ 
    conc["wm_H2S"] = 0.
    frac = melt_species_ratios(conc,species)                                                     
                                 
    return conc, frac           


###################################
### concentration of insolubles ### 
###################################

def conc_insolubles(PT,melt_wf,species,models):
    CO2 = melt_wf["CO2"] # weight fraction CO2
    C_CO2_ = (species.loc["C","M"]*CO2)/species.loc["CO2","M"]
    H2O = melt_wf["H2OT"] # weight fraction H2O
    H_H2O = (2.*species.loc["H","M"]*H2O)/species.loc["H2O","M"]
    H2 = (mg.C_H2(PT,melt_wf,species,models)*mg.f_H2(PT,melt_wf,species,models))/1000000. # weight fraction H2
    H_H2 = (2.*species.loc["H","M"]*H2)/species.loc["H2","M"]
    CH4 = (mg.C_CH4(PT,models)*mg.f_CH4(PT,melt_wf,species,models))/1000000. # weight fraction CH4
    H_CH4 = (4.*species.loc["H","M"]*CH4)/species.loc["CH4","M"]
    C_CH4_ = (species.loc["C","M"]*CH4)/species.loc["CH4","M"]
    CO = (mg.C_CO(PT,models)*mg.f_CO(PT,melt_wf,species,models))/1000000. # weight fraction CO
    C_CO_ = (species.loc["C","M"]*CO)/species.loc["CO","M"]
    S2m = melt_wf["S2-"] # weight fraction of S2-
    S6p = (mg.C_SO4(PT,melt_wf,species,models)*mdv.f_O2(PT,melt_wf,species,models)**2*S2m)/mg.C_S(PT,melt_wf,species,models) # weight fraction S6+
    H2S = (mg.C_H2S(PT,melt_wf,species,models)*mg.f_H2S(PT,melt_wf,species,models))/1000000. # weight fraction H2S
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
    
    conc = {"H2":H2, "CH4":CH4, "CO":CO, "H2S":H2S, "S6p":S6p, "C_T":C_T, "H_T":H_T, "S_T":S_T}                             
    frac = {"H2O_HT":H2O_HT, "H2_HT":H2_HT, "CH4_HT":CH4_HT, "H2S_HT":H2S_HT, "CO2_CT":CO2_CT, "CO_CT":CO_CT, "CH4_CT":CH4_CT, "S2m_ST":S2m_ST, "S6p_ST":S6p_ST, "H2S_ST":H2S_ST}
                                 
    return conc, frac

########################################
### measured parameters within error ### 
########################################

def compositions_within_error(run,setup):
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
    
    result = {"SiO2":SiO2, "TiO2":TiO2, "Al2O3":Al2O3, "FeOT":FeOT, "MnO":MnO, "MgO":MgO, "CaO":CaO, "Na2O":Na2O, "K2O":K2O, "P2O5":P2O5, "H2O":H2O, "CO2ppm":CO2ppm, "STppm":STppm, "Fe3FeT":Fe3FeT}
    return result


def calc_isobar_CO2H2O(PT,melt_wf,species,models):
    M_H2O = species.loc['H2O','M']
    M_CO2 = species.loc['CO2','M']
    M_m_ = mg.M_m_SO(melt_wf,species)

    xm_CO2_ = PT["P"]*mdv.y_CO2(PT,species,models)*mdv.C_CO3(PT,melt_wf,species,models) # pure CO2
    Xm_t = xm_CO2_*M_CO2 + (1.-xm_CO2_)*M_m_
    wm_CO2_0 = (xm_CO2_*M_CO2)/Xm_t
            
    if models.loc["Hspeciation","option"] == "none": # fH2O = xmH2OT^2/CH2O
        xm_H2O_ = (PT["P"]*mdv.y_H2O(PT,species,models)*mdv.C_H2O(PT,melt_wf,species,models))**0.5 # pure H2O
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
        pH2O = mg.p_H2O(PT,melt_wf,species,models)
        pCO2 = PT["P"] - pH2O
        xm_CO2_ = pCO2*mdv.y_CO2(PT,species,models)*mdv.C_CO3(PT,melt_wf,species,models)
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wf_CO2 = (xm_CO2_*M_CO2)/Xm_t
        results1 = pd.DataFrame([[PT["P"],melt_wf["H2OT"]*100.,wf_CO2*1000000.]])
        results = results.append(results1, ignore_index=True)
        
    results1 = pd.DataFrame([[PT["P"],wm_H2O_0*100.,0.]])
    results = results.append(results1, ignore_index=True)
    return results

def calc_pure_solubility(PT,melt_wf,species,models):
    
    M_H2O = species.loc['H2O','M']
    M_CO2 = species.loc['CO2','M']
    M_m_ = mg.M_m_SO(melt_wf,species)
        
    xm_CO2_ = PT["P"]*mdv.y_CO2(PT,species,models)*mdv.C_CO3(PT,melt_wf,species,models) # pure CO2
    Xm_t = xm_CO2_*M_CO2 + (1.-xm_CO2_)*M_m_
    wm_CO2 = (xm_CO2_*M_CO2)/Xm_t
        
    xm_H2O_ = (PT["P"]*mdv.y_H2O(PT,species,models)*mdv.C_H2O(PT,melt_wf,species,models))**0.5 # pure H2O
    Xm_t = xm_H2O_*M_H2O + (1.-xm_H2O_)*M_m_
    wm_H2O = (xm_H2O_*M_H2O)/Xm_t
       
    results = pd.DataFrame([[PT["P"],wm_H2O*100.,wm_CO2*1000000.]])    
    
    return results

###########################################################################
### P given S content of melt after degassing given conditions of pvsat ### IN PROGRESS
###########################################################################
