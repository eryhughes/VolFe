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
def P_sat(PT,melt_wf,models,Ptol,nr_step,nr_tol):
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
    
    def Pdiff(guess,melt_wf,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(PT,melt_wf,models))
        return difference

    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    melt_wf1["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf1,models)
    melt_wf2["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf2,models)
    #wt_C, wt_O, wt_H, wt_S, wt_Fe, wt_g, Wt = bulk_composition(run,PT,melt_wf1,setup,models)
    #bulk_wf = {"H":wt_H,"C":wt_C,"S":wt_S}
    ms_conc = eq.melt_speciation(PT,melt_wf1,models,nr_step,nr_tol)
    ms_frac = melt_species_ratios(ms_conc)
    melt_wf1["H2OT"] = ms_conc["wm_H2O"]
    melt_wf2["H2OT"] = ms_conc["wm_H2O"]
    melt_wf1["CO2"] = ms_conc["wm_CO2"]
    melt_wf2["CO2"] = ms_conc["wm_CO2"]
    melt_wf1["S2-"] = ms_conc["wm_S2m"]
    melt_wf2["S2-"] = ms_conc["wm_S2m"]
    melt_wf1["ST"] = ST
    melt_wf2["ST"] = ST
    if models.loc["sulfur_saturation","option"] == "True": # must incorporate H2S concentration into S2- for SCSS
        sulfsat = sulfur_saturation(PT,melt_wf2,models)
        melt_wf1["ST"] = sulfsat["ST"]/1000000.
        melt_wf2["ST"] = ST  
    delta1 = Pdiff(guess0,melt_wf1,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,melt_wf1,models)
        guess0 = mg.p_tot(PT,melt_wf1,models)
        guess0 = float(guess0)
        PT["P"] = guess0
        melt_wf1["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf1,models)
        melt_wf2["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf2,models)
        ms_conc = eq.melt_speciation(PT,melt_wf1,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc)
        melt_wf1["H2OT"] = ms_conc["wm_H2O"]
        melt_wf2["H2OT"] = ms_conc["wm_H2O"]
        melt_wf1["CO2"] = ms_conc["wm_CO2"]
        melt_wf2["CO2"] = ms_conc["wm_CO2"]
        melt_wf1["S2-"] = ms_conc["wm_S2m"]
        melt_wf2["S2-"] = ms_conc["wm_S2m"]
        if models.loc["sulfur_saturation","option"] == "True":
            sulfsat = sulfur_saturation(PT,melt_wf2,models)
            melt_wf1["ST"] = sulfsat["ST"]/1000000.
            melt_wf2["ST"] = ST
    else:
        P_sat = guess0
        ms_conc = eq.melt_speciation(PT,melt_wf1,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc)
        
    melt_wf["ST"] = ST
    return P_sat, ms_conc, ms_frac

# for a given melt composition, calcualte the saturation pressure
def P_sat_H2O_CO2(PT,melt_wf,models,Ptol,nr_step,nr_tol): # Pvsat with just H2O and CO2 in vapour
    
    def p_tot_H2O_CO2(PT,melt_wf,models):
        value = mg.p_H2O(PT,melt_wf,models) + mg.p_CO2(PT,melt_wf,models)  
        return value
    
    def Pdiff(guess,melt_wf,models):
        PT["P"] = guess
        difference = abs(guess - p_tot_H2O_CO2(PT,melt_wf,models))
        return difference
    
    if melt_wf['CO2'] == 0 and melt_wf['H2OT'] == 0:
        P_sat = ""
        result = {"xg_H2O":0., "xg_CO2":0., "f_H2O":0., "f_CO2":0., "p_H2O":0., "p_CO2":0., "wm_H2O":0., "wm_CO2":0.}
    
    else:
        guess0 = 40000. # initial guess for pressure
        PT["P"] = guess0
        delta1 = Pdiff(guess0,melt_wf,models)
        while delta1 > Ptol :
            delta1 = Pdiff(guess0,melt_wf,models)
            guess0 = p_tot_H2O_CO2(PT,melt_wf,models)
            guess0 = float(guess0)
            PT["P"] = guess0
        else:
            P_sat = guess0
            xg_H2O_ = mg.xg_H2O(PT,melt_wf,models)
            xg_CO2_ = mg.xg_CO2(PT,melt_wf,models)
            p_H2O_ = mg.p_H2O(PT,melt_wf,models)
            p_CO2_ = mg.p_CO2(PT,melt_wf,models)
            f_H2O_ = mg.f_H2O(PT,melt_wf,models)
            f_CO2_ = mg.f_CO2(PT,melt_wf,models)
            xm_H2O_ = (mdv.C_H2O(PT,melt_wf,models)*f_H2O_)**0.5
            xm_CO2_ = mdv.C_CO3(PT,melt_wf,models)*f_CO2_
            M_m_ = mg.M_m_SO(melt_wf)
            Xm_t = xm_CO2_*mdv.species.loc["CO2","M"] + xm_H2O_*mdv.species.loc["H2O","M"] + (1.0-xm_CO2_-xm_H2O_)*M_m_
            wm_H2O_ = (xm_H2O_*mdv.species.loc["H2O","M"])/Xm_t
            wm_CO2_ = (xm_CO2_*mdv.species.loc["CO2","M"])/Xm_t
    
        result = {"xg_H2O":xg_H2O_, "xg_CO2":xg_CO2_, "f_H2O":f_H2O_, "f_CO2":f_CO2_, "p_H2O":p_H2O_, "p_CO2":p_CO2_, "wm_H2O":wm_H2O_, "wm_CO2":wm_CO2_}
    
    return P_sat, result

# for a given fS2 and fO2, calculate Psat
def P_sat_fO2_fS2(PT,melt_wf,models,Ptol):
    
    def Pdiff(guess,melt_wf,models):
        PT["P"] = guess
        result_fO2fS2 = eq.p_tot_fO2_fS2(PT,melt_wf,models)
        difference = abs(guess - result_fO2fS2["P_tot"])
        return difference

    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    delta1 = Pdiff(guess0,melt_wf,models)
    
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,melt_wf,models)
        result_fO2fS2 = eq.p_tot_fO2_fS2(PT,melt_wf,models)
        guess0 = result_fO2fS2["P_tot"]
        guess0 = float(guess0)
        PT["P"] = guess0
    else:
        result = eq.p_tot_fO2_fS2(PT,melt_wf,models)
    
    return results
        
        
        
###################
### composition ###
###################

# calculate bulk composition, including if a gas phase is present
def bulk_composition(run,PT,melt_wf,setup,models):
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
    S6ST_ = mg.S6ST(PT,melt_wf,models)
    #SCSS_,sulfide_sat,SCAS_,sulfate_sat = sulfur_saturation(wm_ST/100.0,S6ST_)
    #print(P, S6ST_)

    if bulk_composition == "melt-only":
        wt_g = 0.
    elif bulk_composition == "melt+vapor_wtg":
        wt_g = setup.loc[run,"wt_g"]/100.
    elif bulk_composition == "melt+vapor_initialCO2":
        wt_C_ = ((mdv.species.loc['C','M']*(setup.loc[run,"initial_CO2wtpc"]/100.))/mdv.species.loc['CO2','M'])
        wt_g = ((wt_C_/mdv.species.loc["C","M"]) - (wm_CO2/mdv.species.loc["CO2","M"]))/(((mg.xg_CO2(PT,melt_wf,models)+mg.xg_CO(PT,melt_wf,models)+mg.xg_CH4(PT,melt_wf,models)+mg.xg_OCS(PT,melt_wf,models))/mg.Xg_tot(PT,melt_wf,models)) - (wm_CO2/mdv.species.loc["CO2","M"]))    

    if bulk_composition == "melt+vapor_initialCO2":
        wt_C = wt_C_
    else:
        wt_C = mdv.species.loc["C","M"]*((wt_g*(((mg.xg_CO2(PT,melt_wf,models)+
                                      mg.xg_CO(PT,melt_wf,models)+
                                      mg.xg_CH4(PT,melt_wf,models)+
                                      mg.xg_OCS(PT,melt_wf,models))
                                     /mg.Xg_tot(PT,melt_wf,models)) 
                                    - (wm_CO2/mdv.species.loc["CO2","M"] - wm_CO/mdv.species.loc["CO","M"] - 
                                       wm_CH4/mdv.species.loc["CH4","M"]))) + (wm_CO2/mdv.species.loc["CO2","M"] + 
                                        wm_CO/mdv.species.loc["CO","M"] + wm_CH4/mdv.species.loc["CH4","M"]))

    wt_H = 2.*mdv.species.loc["H","M"]*((wt_g*(((mg.xg_H2O(PT,melt_wf,models)+mg.xg_H2(PT,melt_wf,models)+2.*mg.xg_CH4(PT,melt_wf,models)+mg.xg_H2S(PT,melt_wf,models))/mg.Xg_tot(PT,melt_wf,models)) - (wm_H2O/mdv.species.loc["H2O","M"] - wm_H2/mdv.species.loc["H2","M"] - wm_H2S/mdv.species.loc["H2S","M"] - (2.*wm_CH4)/mdv.species.loc["CH4","M"]))) + (wm_H2O/mdv.species.loc["H2O","M"] + wm_H2/mdv.species.loc["H2","M"] + wm_H2S/mdv.species.loc["H2S","M"] + (2.*wm_CH4)/mdv.species.loc["CH4","M"]))
    
    wt_S = mdv.species.loc["S","M"]*((wt_g*(((mg.xg_SO2(PT,melt_wf,models)
                                      +2.0*mg.xg_S2(PT,melt_wf,models)
                                      +mg.xg_H2S(PT,melt_wf,models)
                                      +mg.xg_OCS(PT,melt_wf,models))
                                     /mg.Xg_tot(PT,melt_wf,models)) 
                                    - ((wm_S2m+wm_S6p)/mdv.species.loc["S","M"] - (wm_H2S/mdv.species.loc["H2S","M"])))) 
                              + ((wm_S2m+wm_S6p)/mdv.species.loc["S","M"] + (wm_H2S/mdv.species.loc["H2S","M"])))
    
    species_X = models.loc["species X","option"]
    wt_X = mdv.species.loc[species_X,"M"]*((wt_g*(((mg.xg_X(PT,melt_wf,models))
                                     /mg.Xg_tot(PT,melt_wf,models)) 
                                    - wm_XT/mdv.species.loc[species_X,"M"])) 
                              + wm_XT/mdv.species.loc[species_X,"M"])
    
    if models.loc["mass_volume","option"] == "mass":
        if "total_mass_g" in setup:
            Wt = setup.loc[run, "total_mass_g"]
        else:
            Wt = 1.
    elif models.loc["mass_volume","option"] == "volume": ### THIS NEEDS FIXING ###
        Wt = 0.

    if eq_Fe == "no":
        wt_Fe = 0.0
    elif eq_Fe == "yes":
        total_dissolved_volatiles = (wm_CO2 + wm_H2O + wm_ST*(1.-S6ST_) + (mdv.species.loc["SO3","M"]*((wm_ST*S6ST_)/mdv.species.loc["S","M"])))
        wt_Fe = (1.-wt_g)*(((1.0-total_dissolved_volatiles)*mg.wm_Fe_nv(melt_wf))/100.0) # wt fraction of Fe

    wt_O = mdv.species.loc["O","M"]*((wt_g*(((2.0*mg.xg_CO2(PT,melt_wf,models) + mg.xg_CO(PT,melt_wf,models) + 2.0*mg.xg_O2(PT,melt_wf,models) + mg.xg_H2O(PT,melt_wf,models) + 2.0*mg.xg_SO2(PT,melt_wf,models) + mg.xg_OCS(PT,melt_wf,models))/mg.Xg_tot(PT,melt_wf,models)) - (wm_H2O/mdv.species.loc["H2O","M"]) - ((2.0*wm_CO2)/mdv.species.loc["CO2","M"]) - (3.0*wm_S6p/mdv.species.loc["S","M"] - (wm_CO/mdv.species.loc["CO","M"])))) + (wm_H2O/mdv.species.loc["H2O","M"]) + ((2.0*wm_CO2)/mdv.species.loc["CO2","M"]) + (3.0*wm_S6p/mdv.species.loc["S","M"]) + (wm_CO/mdv.species.loc["CO","M"]) + (wt_Fe/mdv.species.loc["Fe","M"])*((1.5*Fe3Fe2_+1.0)/(Fe3Fe2_+1.0)))
    
    result = {"wt_C":wt_C, "wt_O":wt_O, "wt_H":wt_H, "wt_S":wt_S, "wt_X":wt_X, "wt_Fe":wt_Fe, "wt_g":wt_g, "Wt":Wt}
    
    return result

# calculate weight fraction of elements in the system when adding gas into a melt
def new_bulk_regas_open(PT,melt_wf,bulk_wf,gas_mf,dwtg,models):
    me = mg.melt_elements(PT,melt_wf,bulk_wf,gas_mf,models)
    ge = mg.gas_elements(gas_mf,models)
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
def gas_comp_all_open(xg,xg_all,models):
    species_X = models.loc["species X","option"]
    wt_g_before = xg_all["wt_g"]
    wt_g_inst = xg["wt_g"]*(1.-wt_g_before)
    wt_g_new = wt_g_before + wt_g_inst
    f_before = wt_g_before/wt_g_new
    f_inst = wt_g_inst/wt_g_new
    wf_before = mg.gas_wf(xg_all,models)
    wf_inst = mg.gas_wf(xg,models)
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
    tot = (gas_wf_all_new["X"]/mdv.species.loc[species_X,"M"]) + (gas_wf_all_new["OCS"]/mdv.species.loc["OCS","M"]) + (gas_wf_all_new["H2S"]/mdv.species.loc["H2S","M"]) + (gas_wf_all_new["SO2"]/mdv.species.loc["SO2","M"]) + (gas_wf_all_new["CH4"]/mdv.species.loc["CH4","M"]) + (gas_wf_all_new["H2"]/mdv.species.loc["H2","M"]) + (gas_wf_all_new["H2O"]/mdv.species.loc["H2O","M"]) + (gas_wf_all_new["CO2"]/mdv.species.loc["CO2","M"]) + (gas_wf_all_new["O2"]/mdv.species.loc["O2","M"]) + (gas_wf_all_new["CO"]/mdv.species.loc["CO","M"]) + (gas_wf_all_new["S2"]/mdv.species.loc["S2","M"])
    gas_mf_all_new = {"O2":(gas_wf_all_new["O2"]/mdv.species.loc["O2","M"])/tot}
    gas_mf_all_new["CO"] = (gas_wf_all_new["CO"]/mdv.species.loc["CO","M"])/tot
    gas_mf_all_new["S2"] = (gas_wf_all_new["S2"]/mdv.species.loc["S2","M"])/tot
    gas_mf_all_new["CO2"] = (gas_wf_all_new["CO2"]/mdv.species.loc["CO2","M"])/tot
    gas_mf_all_new["H2O"] = (gas_wf_all_new["H2O"]/mdv.species.loc["H2O","M"])/tot
    gas_mf_all_new["H2"] = (gas_wf_all_new["H2"]/mdv.species.loc["H2","M"])/tot
    gas_mf_all_new["CH4"] = (gas_wf_all_new["CH4"]/mdv.species.loc["CH4","M"])/tot
    gas_mf_all_new["SO2"] = (gas_wf_all_new["SO2"]/mdv.species.loc["SO2","M"])/tot
    gas_mf_all_new["H2S"] = (gas_wf_all_new["H2S"]/mdv.species.loc["H2S","M"])/tot
    gas_mf_all_new["OCS"] = (gas_wf_all_new["OCS"]/mdv.species.loc["OCS","M"])/tot
    gas_mf_all_new["X"] = (gas_wf_all_new["X"]/mdv.species.loc[species_X,"M"])/tot
    gas_mf_all_new["Xg_t"] = gas_mf_all_new["CO2"]*mdv.species.loc["CO2","M"] + gas_mf_all_new["CO"]*mdv.species.loc["CO","M"] + gas_mf_all_new["O2"]*mdv.species.loc["O2","M"] + gas_mf_all_new["H2O"]*mdv.species.loc["H2O","M"] + gas_mf_all_new["H2"]*mdv.species.loc["H2","M"] + gas_mf_all_new["CH4"]*mdv.species.loc["CH4","M"] + gas_mf_all_new["SO2"]*mdv.species.loc["SO2","M"] + gas_mf_all_new["S2"]*mdv.species.loc["S2","M"] + gas_mf_all_new["H2S"]*mdv.species.loc["H2S","M"] + gas_mf_all_new["OCS"]*mdv.species.loc["OCS","M"] + gas_mf_all_new["X"]*mdv.species.loc[species_X,"M"]
    gas_mf_all_new["wt_g"] = wt_g_new
    return gas_mf_all_new

# calculate ratios of volatile species in the melt
def melt_species_ratios(conc):

    M_H = mdv.species.loc['H','M']
    M_S = mdv.species.loc['S','M']
    M_C = mdv.species.loc['C','M']
    M_CO = mdv.species.loc['CO','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_H2S = mdv.species.loc['H2S','M']

    wt_H = conc["wm_H2"] + ((2.*(conc["wm_H2O"]/M_H2O))*M_H) + ((4.*(conc["wm_CH4"]/M_CH4))*M_H) + ((2.*(conc["wm_H2S"]/M_H2S))*M_H)
    if wt_H > 0.: # contains H
        H2_HT = conc["wm_H2"]/wt_H
        H2O_HT = ((2.*(conc["wm_H2O"]/M_H2O))*M_H)/wt_H
        CH4_HT = ((4.*(conc["wm_CH4"]/M_CH4))*M_H)/wt_H
        H2S_HT = ((2.*(conc["wm_H2S"]/M_H2S))*M_H)/wt_H
    else:
        H2_HT, H2O_HT, CH4_HT, H2S_HT = "","","", ""
    
    wt_C = ((conc["wm_CO"]/M_CO)*M_C) + ((conc["wm_CO2"]/M_CO2)*M_C) + ((conc["wm_CH4"]/M_CH4)*M_C)
    if wt_C > 0.: # contains C
        CO_CT = ((conc["wm_CO"]/M_CO)*M_C)/wt_C
        CO2_CT = ((conc["wm_CO2"]/M_CO2)*M_C)/wt_C
        CH4_CT = ((conc["wm_CH4"]/M_CH4)*M_C)/wt_C
    else:
        CO_CT, CO2_CT, CH4_CT = "","",""

    wt_S = conc["wm_S2m"] + conc["wm_S6p"] + (M_S*(conc["wm_H2S"]/M_H2S))
    if wt_S > 0.: # contains S
        S2m_ST = conc["wm_S2m"]/wt_S
        S6p_ST = conc["wm_S6p"]/wt_S
        H2S_ST = (M_S*(conc["wm_H2S"]/M_H2S))/wt_S
    else:
        S2m_ST, S6p_ST, H2S_ST = "","",""

    frac = {"H2O_HT":H2O_HT, "H2_HT":H2_HT, "CH4_HT":CH4_HT, "CO2_CT":CO2_CT, "CO_CT":CO_CT, "CH4_CT":CH4_CT, "S6p_ST":S6p_ST, "S2m_ST":S2m_ST, "H2S_ST":H2S_ST, "H2S_HT":H2S_HT}
    
    return frac


#########################
### sulfur satuation ###
#########################

# check solid/immiscible liquid sulfur saturation
def sulfur_saturation(PT,melt_wf,models): # melt weight fraction of ST and S6/ST
    wmST = melt_wf['ST']
    S6T = mg.S6ST(PT,melt_wf,models)
    wmS2 = wmST*100.0*10000.0*(1.0-S6T)
    wmS6 = wmST*100.0*10000.0*S6T
    SCSS_ = mdv.SCSS(PT,melt_wf,models)
    SCAS_ = mdv.SCAS(PT,melt_wf,models)
    StCSS = SCSS_/(1.-S6T)
    StCAS = SCAS_/S6T
    if wmS2 < SCSS_ and wmS6 < SCAS_:
        sulfide_sat = "False"
        sulfate_sat = "False"
        ST = wmST*1000000.
    elif wmS2 >= SCSS_ and wmS6 >= SCAS_:
        sulfide_sat = "True"
        sulfate_sat = "True"
        ST = min(StCSS,StCAS)
    elif wmS2 >= SCSS_ and wmS6 < SCAS_:
        sulfide_sat = "True"
        sulfate_sat = "False"
        ST = StCSS
    elif wmS2 < SCSS_ and wmS6 >= SCAS_:
        sulfide_sat = "False"
        sulfate_sat = "True"
        ST = StCAS
    else:
        sulfide_sat = "nan"
        sulfate_sat = "nan"
        ST = wmST*1000000.
    result = {"SCSS":SCSS_,"StCSS":StCSS,"sulfide_sat":sulfide_sat, "SCAS":SCAS_, "StCAS":StCAS,"sulfate_sat":sulfate_sat,"ST":ST}
    return result

# fO2 and P of v+sulf+anh saturation
def fO2_P_VSA(PT,melt_wf,models,nr_step,nr_tol,Ptol):

    def Pdiff(guess,melt_wf,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(PT,melt_wf,models))
        return difference

    def fO2_S(PT,melt_wf,models):
        SCSS_ = mdv.SCSS(PT,melt_wf,models)/1000000.
        SCAS_ = mdv.SCAS(PT,melt_wf,models)/1000000.
        CSO4 = mdv.C_SO4(PT,melt_wf,models)/1000000.
        CS = mdv.C_S(PT,melt_wf,models)/1000000.
        
        if models.loc["H2S_m","option"] == "False":
            W = CSO4/CS
        elif models.loc["H2S_m","option"] == "True":
            CH2S = mdv.C_H2S(PT,melt_wf,models)/1000000.
            KHS = mdv.KHOSg(PT,models)
            CH2OT = mdv.C_H2O(PT,melt_wf,models)
            xmH2O = mg.xm_H2OT_so(melt_wf)
            W = (CSO4/((KHS*CH2S*(xmH2O**2./CH2OT)) + CS))
        fO2 = ((1./W)*(SCAS_/SCSS_))**0.5
        ST = SCSS_ + SCAS_
        return fO2, ST
    
    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    melt_wf["Fe3FeT"] = 0.1
    fO2_, ST_ = fO2_S(PT,melt_wf,models)
    melt_wf["ST"] = ST_
    melt_wf["Fe3FeT"] = mdv.fO22Fe3FeT(fO2_,PT,melt_wf,models)
    ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
    ms_frac = melt_species_ratios(ms_conc)
    melt_wf["H2OT"] = ms_conc["wm_H2O"]
    melt_wf["CO2"] = ms_conc["wm_CO2"]
    melt_wf["S2-"] = ms_conc["wm_S2m"]
    delta1 = Pdiff(guess0,melt_wf,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,melt_wf,models)
        guess0 = mg.p_tot(PT,melt_wf,models)
        guess0 = float(guess0)
        PT["P"] = guess0
        fO2_, ST_ = fO2_S(PT,melt_wf,models)
        melt_wf["ST"] = ST_
        melt_wf["Fe3FeT"] = mdv.fO22Fe3FeT(fO2_,PT,melt_wf,models)
        ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc)        
        melt_wf["H2OT"] = ms_conc["wm_H2O"]
        melt_wf["CO2"] = ms_conc["wm_CO2"]
        melt_wf["S2-"] = ms_conc["wm_S2m"]
    else:
        P_sat = guess0
        ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc)
        fO2_, ST_ = fO2_S(PT,melt_wf,models)
        Fe3_FT = mdv.fO22Fe3FeT(fO2_,PT,melt_wf,models)

    ms_frac["Fe3_FT"] = Fe3_FT 
              
    return P_sat, ms_conc, ms_frac

def P_VSA(PT,melt_wf,models,nr_step,nr_tol,Ptol):
    
    P_sat, pVSA_conc, pVSA_frac = fO2_P_VSA(PT,melt_wf,models,nr_step,nr_tol,Ptol)
    Fe3_FT = pVSA_frac['Fe3_FT']
    
    def Pdiff(guess,melt_wf,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(PT,melt_wf,models))
        return difference
    
    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,models)
    
    if melt_wf["Fe3FeT"] < Fe3_FT:
        SCSS_ = mdv.SCSS(PT,melt_wf,models)/1000000.
        S6T = mg.S6ST(PT,melt_wf,models)
        ST_ = SCSS_/(1.-S6T)
    else:
        SCAS_ = mdv.SCAS(PT,melt_wf,models)/1000000.
        S6T = mg.S6ST(PT,melt_wf,models)
        ST_ = SCAS_/S6T
    melt_wf["ST"] = ST_
    ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
    ms_frac = melt_species_ratios(ms_conc)
    melt_wf["H2OT"] = ms_conc["wm_H2O"]
    melt_wf["CO2"] = ms_conc["wm_CO2"]
    melt_wf["S2-"] = ms_conc["wm_S2m"]
    delta1 = Pdiff(guess0,melt_wf,models)
    while delta1 > Ptol :
        delta1 = Pdiff(guess0,melt_wf,models)
        guess0 = mg.p_tot(PT,melt_wf,models)
        guess0 = float(guess0)
        PT["P"] = guess0
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,models)
        if melt_wf["Fe3FeT"] < Fe3_FT:
            SCSS_ = mdv.SCSS(PT,melt_wf,models)/1000000.
            S6T = mg.S6ST(PT,melt_wf,models)
            ST_ = SCSS_/(1.-S6T)
        else:
            SCAS_ = mdv.SCAS(PT,melt_wf,models)/1000000.
            S6T = mg.S6ST(PT,melt_wf,models)
            ST_ = SCAS_/S6T
        melt_wf["ST"] = ST_
        ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc)    
        melt_wf["H2OT"] = ms_conc["wm_H2O"]
        melt_wf["CO2"] = ms_conc["wm_CO2"]
        melt_wf["S2-"] = ms_conc["wm_S2m"]
    else:
        P_sat = guess0
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,models)
        ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc)
    
    return P_sat,ms_conc,ms_frac

# estimated fO2 from sulfur content given you are sulfide saturated and at Pvsat
def P_sat_sulf_anh(PT,melt_wf,models,Ptol,nr_step,nr_tol):
    ST = melt_wf["ST"]
    H2OT = melt_wf["H2OT"]
    CO2 = melt_wf["CO2"]
    HT = melt_wf["HT"]
    CT = melt_wf["CT"]
 
    def Pdiff(guess,melt_wf,models):
        PT["P"] = guess
        difference = abs(guess - mg.p_tot(PT,melt_wf,models))
        return difference
    
    def S62_2_Fe3T(PT,melt_wf,models,sat):
        ST = melt_wf["ST"]
        if sat == "sulf":
            Ssat = mdv.SCSS(PT,melt_wf,models)/1000000.
            S2m = Ssat
            S6p = ST - S2m
        elif sat == "anh":
            Ssat = mdv.SCAS(PT,melt_wf,models)/1000000.
            S6p = Ssat
            S2m = ST - S6p
        S62 = S6p/S2m
        S6T = S6p/ST
        if S62 < 0.:
            return "not possible","","","",Ssat
        else:
            fO2 = mg.S6S2_2_fO2(S62,melt_wf,PT,models)
            Fe3T = mdv.fO22Fe3FeT(fO2,PT,melt_wf,models)
            DFMQ = mg.fO22Dbuffer(PT,fO2,"FMQ",models)
            return Fe3T,fO2,S6T,DFMQ,Ssat
        
    guess0 = 40000. # initial guess for pressure
    PT["P"] = guess0
    
    # assume it is sulfide saturated
    melt_wf["Fe3FeT"] = 0.1
    Fe3T_sulf,fO2_sulf,S6T_sulf,DFMQ_sulf,SCSS_ = S62_2_Fe3T(PT,melt_wf,models,"sulf")
    if Fe3T_sulf == "not possible":
        P_sat_sulf = ""
        Fe3T_sulf = ""
        sulfide_sat = "False"
    else: 
        melt_wf["Fe3FeT"] = Fe3T_sulf
        ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc)
        melt_wf["H2OT"] = ms_conc["wm_H2O"]
        melt_wf["CO2"] = ms_conc["wm_CO2"]
        melt_wf["S2-"] = ms_conc["wm_S2m"]
        melt_wf["ST"] = ST
        delta1 = Pdiff(guess0,melt_wf,models)
        while delta1 > Ptol :
            delta1 = Pdiff(guess0,melt_wf,models)
            guess0 = mg.p_tot(PT,melt_wf,models)
            guess0 = float(guess0)
            PT["P"] = guess0
            Fe3T_sulf,fO2_sulf,S6T_sulf,DFMQ_sulf,SCSS_ = S62_2_Fe3T(PT,melt_wf,models,"sulf")
            if Fe3T_sulf == "not possible":
                P_sat_sulf = ""
                Fe3T_sulf = ""
                sulfide_sat = "False"
            else: 
                melt_wf["Fe3FeT"] = Fe3T_sulf
                ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
                ms_frac = melt_species_ratios(ms_conc)
                melt_wf["H2OT"] = ms_conc["wm_H2O"]
                melt_wf["CO2"] = ms_conc["wm_CO2"]
                melt_wf["S2-"] = ms_conc["wm_S2m"]
        else:
            P_sat_sulf = guess0
            Fe3T_sulf,fO2_sulf,S6T_sulf,DFMQ_sulf,SCSS_ = S62_2_Fe3T(PT,melt_wf,models,"sulf")
            sulfide_sat = "True"
           
    # assume it is anhydrite saturated
    Fe3T_anh,fO2_anh,S6T_anh,DFMQ_anh,SCAS_ = S62_2_Fe3T(PT,melt_wf,models,"anh")
    if Fe3T_anh == "not possible":
        P_sat_anh = ""
        Fe3T_anh = ""
        anhydrite_sat = "False"
    else:  
        melt_wf["Fe3FeT"] = Fe3T_anh
        ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
        ms_frac = melt_species_ratios(ms_conc)
        melt_wf["H2OT"] = ms_conc["wm_H2O"]
        melt_wf["CO2"] = ms_conc["wm_CO2"]
        melt_wf["S2-"] = ms_conc["wm_S2m"]
        melt_wf["ST"] = ST
        delta1 = Pdiff(guess0,melt_wf,models)
        while delta1 > Ptol :
            delta1 = Pdiff(guess0,melt_wf,models)
            guess0 = mg.p_tot(PT,melt_wf,models)
            guess0 = float(guess0)
            PT["P"] = guess0
            Fe3T_anh,fO2_anh,S6T_anh,DFMQ_anh,SCAS_ = S62_2_Fe3T(PT,melt_wf,models,"anh")
            if Fe3T_anh == "not possible":
                P_sat_anh = ""
                Fe3T_anh = ""
                anhydrite_sat = "False"
            else:  
                melt_wf["Fe3FeT"] = Fe3T_anh
                ms_conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
                ms_frac = melt_species_ratios(ms_conc)
                melt_wf["H2OT"] = ms_conc["wm_H2O"]
                melt_wf["CO2"] = ms_conc["wm_CO2"]
                melt_wf["S2-"] = ms_conc["wm_S2m"]
        else:
            P_sat_anh = guess0
            Fe3T_anh,fO2_anh,S6T_anh,DFMQ_anh,SCAS_ = S62_2_Fe3T(PT,melt_wf,models,"anh")
            anhydrite_sat = "True"
    
    result = {"P_sat_sulf":P_sat_sulf,"P_sat_anh":P_sat_anh,"SCAS":SCAS_*1000000.,"SCSS":SCSS_*1000000.,"sulf_sat":sulfide_sat,"DFMQ_sulf":DFMQ_sulf,"fO2_sulf":fO2_sulf,"Fe3T_sulf":Fe3T_sulf,"S6T_sulf":S6T_sulf,"anh_sat":anhydrite_sat,"DFMQ_anh":DFMQ_anh,"fO2_anh":fO2_anh,"Fe3T_anh":Fe3T_anh,"S6T_anh":S6T_anh}          

    return result

##########################
### graphite satuation ###
##########################

# check graphite saturation
def graphite_saturation(PT,melt_wf,models): # needs finishing
    K1 = mg.f_CO2(PT,melt_wf,models)/mdv.f_O2(PT,melt_wf,models)
    K2 = mdv.KCOs(PT,models) # K for graphite saturation
    if K1 < K2:
        graphite_sat = "False"
    else:
        graphite_sat = "True"
    fCO2_ = K2*mdv.f_O2(PT,melt_wf,models)
    xmCO2 = fCO2_*mdv.C_CO3(PT,melt_wf,models)
    return graphite_sat

                
######################################                
### fO2 range from sulfur content ###
######################################
                
def fO2_range_from_S(PT,melt_wf,models):
    SCSS_ = mdv.SCSS(PT,melt_wf,models)
    SCAS_ = mdv.SCAS(PT,melt_wf,models)

    if melt_wf["ST"]*1000000. > SCSS_:
        sulfide_sat = "possible"
        S6ST_1 = (melt_wf["ST"]*1000000. - SCSS_)/(melt_wf["ST"]*1000000.)
        S6S2_1 = mg.overtotal2ratio(S6ST_1)
        fO2_1 = mg.S6S2_2_fO2(S6S2_1,melt_wf,PT,models)
        Fe3FeT_1 = mdv.fO22Fe3FeT(fO2_1,PT,melt_wf,models)
        Fe3Fe2_1 = mg.overtotal2ratio(Fe3FeT_1)
        DFMQ_1 = mg.fO22Dbuffer(PT,fO2_1,"FMQ",models)
    else:
        sulfide_sat = "False"
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
        fO2_2 = mg.S6S2_2_fO2(S6S2_2,melt_wf,PT,models)
        Fe3FeT_2 = mdv.fO22Fe3FeT(fO2_2,PT,melt_wf,models)
        Fe3Fe2_2 = mg.overtotal2ratio(Fe3FeT_2)
        DFMQ_2 = mg.fO22Dbuffer(PT,fO2_2,"FMQ",models)
    else:
        anhydrite_sat = "False"
        S6ST_2 = ""
        S6S2_2 = ""
        Fe3Fe2_2 = ""
        Fe3FeT_2 = ""
        fO2_2 = ""
        DFMQ_2 = ""
                                 
    result = {"SCAS":SCAS_, "SCSS":SCSS_, "sulf_sat":sulfide_sat, "DFMQ_sulf":DFMQ_1, "fO2_sulf":fO2_1, "Fe3T_sulf":Fe3FeT_1, "S6T_sulf":S6ST_1, "anh_sat":anhydrite_sat, "DFMQ_anh":DFMQ_2, "fO2_anh":fO2_2, "Fe3T_anh":Fe3FeT_2, "S6T_anh":S6ST_2}

    return result
    
                             
#################################                             
### Mass, volume, and density ###
#################################                             
                             
def mass_vol_rho(PT,melt_wf,gas_mf,bulk_wf,models):
    gas_m = gas_mf['wt_g']*bulk_wf['Wt']
    gas_v = mg.gas_volume(PT,gas_mf,bulk_wf,models)
    gas_rho = mg.gas_density(PT,gas_mf,bulk_wf,models)
    melt_m = (1.-gas_mf['wt_g'])*bulk_wf['Wt']
    melt_v = mg.melt_volume(PT,melt_wf,bulk_wf,gas_mf)
    melt_rho = mg.melt_density(PT,melt_wf)
    tot_m = bulk_wf['Wt']
    tot_v = gas_v + melt_v
    tot_rho = mg.system_density(PT,melt_wf,gas_mf,bulk_wf,models)
    result = {"tot_m":tot_m, "tot_v":tot_v, "tot_rho":tot_rho, "melt_m":melt_m, "melt_v":melt_v, "melt_rho":melt_rho, "gas_m":gas_m, "gas_v":gas_v, "gas_rho":gas_rho}
    return result


######################################################
### mole fraction of elements in different species ###
######################################################

def mf_S_species_old(melt_wf,gas_mf):
    # weight of S in each sulfur-bearing species
    W_S2m = melt_wf["ST"]*(1.-gas_mf["wt_g"])*(1.-melt_wf["S6ST"])
    W_SO4 = melt_wf["ST"]*(1.-gas_mf["wt_g"])*melt_wf["S6ST"]
    W_H2S = ((gas_mf["H2S"]*mdv.species.loc["S","M"])/(gas_mf["Xg_t"]))*gas_mf["wt_g"]
    W_SO2 = ((gas_mf["SO2"]*mdv.species.loc["S","M"])/(gas_mf["Xg_t"]))*gas_mf["wt_g"]
    W_S2 = ((gas_mf["S2"]*mdv.species.loc["S2","M"])/(gas_mf["Xg_t"]))*gas_mf["wt_g"]
    W_OCS = ((gas_mf["OCS"]*mdv.species.loc["S","M"])/(gas_mf["Xg_t"]))*gas_mf["wt_g"]
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

def mf_S_species(comp):
    # weight of S in each sulfur-bearing species
    wtg = (float(comp["wt_g_wtpc"].iloc[0]))/100.
    wtm = 1. - wtg
    W_S2m = (float(comp["S2-_ppmw"].iloc[0])/1000000.)*wtm
    W_SO4 = (float(comp["S6+_ppmw"].iloc[0])/1000000.)*wtm
    W_H2Smol = (((float(comp["H2S_ppmw"].iloc[0])/1000000.)*mdv.species.loc["S","M"])/mdv.species.loc["H2S","M"])*wtm
    Xg_t = (float(comp["xgO2_mf"].iloc[0])*mdv.species.loc["O2","M"]) + (float(comp["xgH2_mf"].iloc[0])*mdv.species.loc["H2","M"]) + (float(comp["xgH2O_mf"].iloc[0])*mdv.species.loc["H2O","M"]) + (float(comp["xgS2_mf"].iloc[0])*mdv.species.loc["S2","M"]) + (float(comp["xgSO2_mf"].iloc[0])*mdv.species.loc["SO2","M"]) + (float(comp["xgH2S_mf"].iloc[0])*mdv.species.loc["H2S","M"]) + (float(comp["xgCO_mf"].iloc[0])*mdv.species.loc["CO","M"]) + (float(comp["xgCO2_mf"].iloc[0])*mdv.species.loc["CO2","M"]) + (float(comp["xgCH4_mf"].iloc[0])*mdv.species.loc["CH4","M"]) + (float(comp["xgOCS_mf"].iloc[0])*mdv.species.loc["OCS","M"])
    W_H2S = ((float(comp["xgH2S_mf"].iloc[0])*mdv.species.loc["S","M"])/Xg_t)*wtg
    W_SO2 = ((float(comp["xgSO2_mf"].iloc[0])*mdv.species.loc["S","M"])/Xg_t)*wtg
    W_S2 = ((float(comp["xgS2_mf"].iloc[0])*mdv.species.loc["S2","M"])/Xg_t)*wtg
    W_OCS = ((float(comp["xgOCS_mf"].iloc[0])*mdv.species.loc["S","M"])/Xg_t)*wtg
    W_total = W_S2m + W_SO4 + W_H2Smol + W_H2S + W_SO2 + W_S2 + W_OCS
    # weight and mole fraction of S in each sulfur-bearing species compared to total S
    w_S2m = W_S2m/W_total
    w_SO4 = W_SO4/W_total
    w_H2Smol = W_H2Smol/W_total
    w_SO2 = W_SO2/W_total
    w_H2S = W_H2S/W_total
    w_S2 = W_S2/W_total
    w_OCS = W_OCS/W_total
    mf = {"S2-":w_S2m, "SO42-":w_SO4, "H2Smol":w_H2Smol, "SO2":w_SO2, "H2S":w_H2S, "S2": w_S2, "OCS": w_OCS}
    return mf

def mf_C_species(comp):
    # weight of C in each carbon-bearing species
    wtg = (float(comp["wt_g_wtpc"].iloc[0]))/100.
    wtm = 1. - wtg
    W_COmol = (((float(comp["CO_ppmw"].iloc[0])/1000000.)*mdv.species.loc["C","M"])/mdv.species.loc["CO","M"])*wtm
    W_CO2mol = (((float(comp["CO2mol_ppmw"].iloc[0])/1000000.)*mdv.species.loc["C","M"])/mdv.species.loc["CO2","M"])*wtm
    W_carb = (((float(comp["CO2carb_ppmw"].iloc[0])/1000000.)*mdv.species.loc["C","M"])/mdv.species.loc["CO2","M"])*wtm
    W_CH4mol = (((float(comp["CH4_ppmw"].iloc[0])/1000000.)*mdv.species.loc["C","M"])/mdv.species.loc["CH4","M"])*wtm
    Xg_t = (float(comp["xgO2_mf"].iloc[0])*mdv.species.loc["O2","M"]) + (float(comp["xgH2_mf"].iloc[0])*mdv.species.loc["H2","M"]) + (float(comp["xgH2O_mf"].iloc[0])*mdv.species.loc["H2O","M"]) + (float(comp["xgS2_mf"].iloc[0])*mdv.species.loc["S2","M"]) + (float(comp["xgSO2_mf"].iloc[0])*mdv.species.loc["SO2","M"]) + (float(comp["xgH2S_mf"].iloc[0])*mdv.species.loc["H2S","M"]) + (float(comp["xgCO_mf"].iloc[0])*mdv.species.loc["CO","M"]) + (float(comp["xgCO2_mf"].iloc[0])*mdv.species.loc["CO2","M"]) + (float(comp["xgCH4_mf"].iloc[0])*mdv.species.loc["CH4","M"]) + (float(comp["xgOCS_mf"].iloc[0])*mdv.species.loc["OCS","M"])
    W_CO = ((float(comp["xgCO_mf"].iloc[0])*mdv.species.loc["C","M"])/Xg_t)*wtg
    W_CO2 = ((float(comp["xgCO2_mf"].iloc[0])*mdv.species.loc["C","M"])/Xg_t)*wtg
    W_CH4 = ((float(comp["xgCH4_mf"].iloc[0])*mdv.species.loc["C","M"])/Xg_t)*wtg
    W_OCS = ((float(comp["xgOCS_mf"].iloc[0])*mdv.species.loc["C","M"])/Xg_t)*wtg
    W_total = W_COmol + W_CO2mol +  W_carb + W_CH4mol + W_CO + W_CO2 + W_CH4 + W_OCS
    # weight and mole fraction of C in each carbon-bearing species compared to total C
    w_COmol = W_COmol/W_total
    w_CO2mol = W_CO2mol/W_total
    w_carb = W_carb/W_total
    w_CH4mol = W_CH4mol/W_total
    w_CO = W_CO/W_total
    w_CO2 = W_CO2/W_total
    w_CH4 = W_CH4/W_total
    w_OCS = W_OCS/W_total
    mf = {"COmol":w_COmol, "CO2mol":w_CO2mol, "CO32-":w_carb, "CH4mol":w_CH4mol, "CO":w_CO, "CO2":w_CO2, "CH4":w_CH4, "OCS":w_OCS}
    return mf

"xgO2_mf","xgH2_mf","xgH2O_mf","xgS2_mf","xgSO2_mf","xgH2S_mf","xgCO2_mf","xgCO_mf","xgCH4_mf","xgOCS_mf","xgX_mf","xgC_S_mf",
"H2OT_wtpc","OH_wtpc","H2Omol_wtpc","H2_ppmw","CH4_ppmw","CO2T_ppmw","CO2mol_ppmw","CO32-_ppmw","CO_ppmw","S2-_ppmw","S6+_ppmw","H2S_ppmw",

def mf_H_species(comp):
    # weight of H in each hydrogen-bearing species
    wtg = (float(comp["wt_g_wtpc"].iloc[0]))/100.
    wtm = 1. - wtg
    W_H2mol = (float(comp["H2_ppmw"].iloc[0])/1000000.)*wtm
    W_H2Omol = (((float(comp["H2Omol_wtpc"].iloc[0])/100.)*mdv.species.loc["H2","M"])/mdv.species.loc["H2O","M"])*wtm
    W_OH = (((float(comp["OH_wtpc"].iloc[0])/100.)*0.5*mdv.species.loc["H2","M"])/mdv.species.loc["OH","M"])*wtm
    W_CH4mol = (((float(comp["CH4_ppmw"].iloc[0])/1000000.)*2.*mdv.species.loc["H2","M"])/mdv.species.loc["CH4","M"])*wtm
    W_H2Smol = (((float(comp["H2S_ppmw"].iloc[0])/1000000.)*mdv.species.loc["H2","M"])/mdv.species.loc["H2S","M"])*wtm
    Xg_t = (float(comp["xgO2_mf"].iloc[0])*mdv.species.loc["O2","M"]) + (float(comp["xgH2_mf"].iloc[0])*mdv.species.loc["H2","M"]) + (float(comp["xgH2O_mf"].iloc[0])*mdv.species.loc["H2O","M"])  + (float(comp["xgS2_mf"].iloc[0])*mdv.species.loc["S2","M"]) + (float(comp["xgSO2_mf"].iloc[0])*mdv.species.loc["SO2","M"]) + (float(comp["xgH2S_mf"].iloc[0])*mdv.species.loc["H2S","M"]) + (float(comp["xgCO_mf"].iloc[0])*mdv.species.loc["CO","M"]) + (float(comp["xgCO2_mf"].iloc[0])*mdv.species.loc["CO2","M"]) + (float(comp["xgCH4_mf"].iloc[0])*mdv.species.loc["CH4","M"]) + (float(comp["xgOCS_mf"].iloc[0])*mdv.species.loc["OCS","M"])
    W_H2 = ((float(comp["xgH2_mf"].iloc[0])*mdv.species.loc["H2","M"])/Xg_t)*wtg
    W_H2O = ((float(comp["xgH2O_mf"].iloc[0])*mdv.species.loc["H2","M"])/Xg_t)*wtg
    W_CH4 = ((float(comp["xgCH4_mf"].iloc[0])*2.*mdv.species.loc["H2","M"])/Xg_t)*wtg
    W_H2S = ((float(comp["xgH2S_mf"].iloc[0])*mdv.species.loc["H2","M"])/Xg_t)*wtg
    W_total = W_H2mol + W_H2Omol +  W_OH + W_CH4mol + W_H2Smol + W_H2 + W_H2O + W_CH4 + W_H2S
    # weight and mole fraction of H in each hydrogen-bearing species compared to total H
    w_H2mol = W_H2mol/W_total
    w_H2Omol = W_H2Omol/W_total
    w_OH = W_OH/W_total
    w_CH4mol = W_CH4mol/W_total
    w_H2Smol = W_H2Smol/W_total
    w_H2 = W_H2/W_total
    w_H2O = W_H2O/W_total
    w_CH4 = W_CH4/W_total
    w_H2S = W_H2S/W_total
    mf = {"H2mol":w_H2mol, "H2Omol":w_H2Omol, "OH-":w_OH, "CH4mol":w_CH4mol, "H2Smol":w_H2Smol, "H2":w_H2, "H2O":w_H2O, "CH4":w_CH4, "H2S":w_H2S}
    return mf
        
##############################################
### fO2 of silm+sulfm+anh at given T and P ###
##############################################

def fO2_silm_sulf_anh(PT,melt_wf,models):
    S6 = mdv.SCAS(PT,melt_wf)
    S2 = mdv.SCSS(PT,melt_wf,models)
    fO2 = ((S6*mdv.C_S(PT,melt_wf,models))/(S2*mdv.C_SO4(PT,melt_wf,models)))**0.5
    DFMQ = mg.fO22Dbuffer(PT,fO2,"FMQ",models)
    wmST = S6+S2
    S6ST = S6/wmST
    result = {"fO2":fO2, "DFMQ":DFMQ, "wmST":wmST, "S6ST":S6ST, "S6":S6, "S2":S2}
    return result
   

##############################################
### S content at given T, P, fO2, C, and H ###
##############################################

def S_given_T_P_fO2_C_H(PT,melt_wf,models,nr_step,nr_tol): # no dissolved H2S  
    P = PT["P"]
    
    conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
    
    melt_wf["H2OT"] = conc["wm_H2O"]
    melt_wf["CO2"] = conc["wm_CO2"]
  
    xg_O2_ = mg.p_O2(PT,melt_wf,models)/P
    xg_CO2_ = mg.p_CO2(PT,melt_wf,models)/P
    xg_CO_ = mg.p_CO(PT,melt_wf,models)/P
    xg_H2_ = mg.p_H2(PT,melt_wf,models)/P
    xg_H2O_ = mg.p_H2O(PT,melt_wf,models)/P
    xg_CH4_ = mg.p_CH4(PT,melt_wf,models)/P
    
    K6_ = mdv.KOSg(PT,models)
    K7_ = mdv.KHOSg(PT,models)
    K8_ = mdv.C_S(PT,melt_wf,models)/1000000.0
    K9_ = (mdv.C_SO4(PT,melt_wf,models)/1000000.0)
    K10_ = mdv.KOCSg(PT,models)                                           
    y_S2_ = mdv.y_S2(PT,models)
    y_SO2_ = mdv.y_SO2(PT,models)
    y_H2S_ = mdv.y_H2S(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    y_CO2_ = mdv.y_CO2(PT,models)
    y_CO_ = mdv.y_CO(PT,models)
    y_OCS_ = mdv.y_OCS(PT,models)
    M_S = mdv.species.loc['S','M']
    M_SO3 = mdv.species.loc['SO3','M']
    
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
    frac = melt_species_ratios(conc)                                                     
                                 
    return conc, frac           


###################################
### concentration of insolubles ### 
###################################

def conc_insolubles(PT,melt_wf,models):
    CO2 = melt_wf["CO2"] # weight fraction CO2
    C_CO2_ = (mdv.species.loc["C","M"]*CO2)/mdv.species.loc["CO2","M"]
    H2O = melt_wf["H2OT"] # weight fraction H2O
    H_H2O = (2.*mdv.species.loc["H","M"]*H2O)/mdv.species.loc["H2O","M"]
    H2 = (mg.C_H2(PT,melt_wf,models)*mg.f_H2(PT,melt_wf,models))/1000000. # weight fraction H2
    H_H2 = (2.*mdv.species.loc["H","M"]*H2)/mdv.species.loc["H2","M"]
    CH4 = (mg.C_CH4(PT,models)*mg.f_CH4(PT,melt_wf,models))/1000000. # weight fraction CH4
    H_CH4 = (4.*mdv.species.loc["H","M"]*CH4)/mdv.species.loc["CH4","M"]
    C_CH4_ = (mdv.species.loc["C","M"]*CH4)/mdv.species.loc["CH4","M"]
    CO = (mg.C_CO(PT,models)*mg.f_CO(PT,melt_wf,models))/1000000. # weight fraction CO
    C_CO_ = (mdv.species.loc["C","M"]*CO)/mdv.species.loc["CO","M"]
    S2m = melt_wf["S2-"] # weight fraction of S2-
    S6p = (mg.C_SO4(PT,melt_wf,models)*mdv.f_O2(PT,melt_wf,models)**2*S2m)/mg.C_S(PT,melt_wf,models) # weight fraction S6+
    H2S = (mg.C_H2S(PT,melt_wf,models)*mg.f_H2S(PT,melt_wf,models))/1000000. # weight fraction H2S
    S_H2S = (mdv.species.loc["S","M"]*H2S)/mdv.species.loc["H2S","M"]
    H_H2S = (2.*mdv.species.loc["H","M"]*H2S)/mdv.species.loc["H2S","M"]
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

    if setup.loc[run,"Xppm_sd_type"] == "A": # absolute
        Xppm_sd = setup.loc[run,"Xppm_sd"]
    else:
        Xppm_sd = setup.loc[run,"Xppm_sd"]*setup.loc[run,"Xppm"]
    Xppm = float(np.random.normal(setup.loc[run,"Xppm"],Xppm_sd,1))
        
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
    
    result = {"SiO2":SiO2, "TiO2":TiO2, "Al2O3":Al2O3, "FeOT":FeOT, "MnO":MnO, "MgO":MgO, "CaO":CaO, "Na2O":Na2O, "K2O":K2O, "P2O5":P2O5, "H2O":H2O, "CO2ppm":CO2ppm, "Xppm":Xppm,"STppm":STppm, "Fe3FeT":Fe3FeT}
    return result


def calc_isobar_CO2H2O(PT,melt_wf,models):
    M_H2O = mdv.species.loc['H2O','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_m_ = mg.M_m_SO(melt_wf)

    xm_CO2_ = PT["P"]*mdv.y_CO2(PT,models)*mdv.C_CO3(PT,melt_wf,models) # pure CO2
    Xm_t = xm_CO2_*M_CO2 + (1.-xm_CO2_)*M_m_
    wm_CO2_0 = (xm_CO2_*M_CO2)/Xm_t
            
    if models.loc["Hspeciation","option"] == "none": # fH2O = xmH2OT^2/CH2O
        xm_H2O_ = (PT["P"]*mdv.y_H2O(PT,models)*mdv.C_H2O(PT,melt_wf,models))**0.5 # pure H2O
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
        melt_wf["H2OT"]=wm_H2O_
        melt_wf["CO2"] = 0.
        pH2O = mg.p_H2O(PT,melt_wf,models)
        pCO2 = PT["P"] - pH2O
        xm_CO2_ = pCO2*mdv.y_CO2(PT,models)*mdv.C_CO3(PT,melt_wf,models)
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wf_CO2 = (xm_CO2_*M_CO2)/Xm_t
        results1 = pd.DataFrame([[PT["P"],melt_wf["H2OT"]*100.,wf_CO2*1000000.]])
        results = pd.concat([results, results1], ignore_index=True)
        
    results1 = pd.DataFrame([[PT["P"],wm_H2O_0*100.,0.]])
    results = pd.concat([results, results1], ignore_index=True)
    return results

def calc_pure_solubility(PT,melt_wf,models):
    
    M_H2O = mdv.species.loc['H2O','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_m_ = mg.M_m_SO(melt_wf)
        
    xm_CO2_ = PT["P"]*mdv.y_CO2(PT,models)*mdv.C_CO3(PT,melt_wf,models) # pure CO2
    Xm_t = xm_CO2_*M_CO2 + (1.-xm_CO2_)*M_m_
    wm_CO2 = (xm_CO2_*M_CO2)/Xm_t
        
    xm_H2O_ = (PT["P"]*mdv.y_H2O(PT,models)*mdv.C_H2O(PT,melt_wf,models))**0.5 # pure H2O
    Xm_t = xm_H2O_*M_H2O + (1.-xm_H2O_)*M_m_
    wm_H2O = (xm_H2O_*M_H2O)/Xm_t
       
    results = pd.DataFrame([[PT["P"],wm_H2O*100.,wm_CO2*1000000.]])    
    
    return results

def calc_isotopes(PT,comp,R,models,guesses,nr_step,nr_tol,run=0.):
    comp_ = comp[run:run+1]
    R_all_species_C, R_m_g_C = iso.i2s9("C",PT,comp_,R,models,guesses['C'],nr_step,nr_tol)
    R_all_species_S, R_m_g_S = iso.i2s9("S",PT,comp_,R,models,guesses["S"],nr_step,nr_tol)
    R_all_species_H, R_m_g_H = iso.i2s9("H",PT,comp_,R,models,guesses['H'],nr_step,nr_tol)
    return R_all_species_S,R_m_g_S,R_all_species_C,R_m_g_C,R_all_species_H,R_m_g_H

# function to calculate the SCSS and SCAS with varying melt composition
def calc_sulfur_vcomp(setup,models=mdv.default_models):
    for n in range(0,len(setup),1):
        melt_wf=vf.melt_comp(n,setup)
        PT = {"T":setup.loc[n,"T_C"],"P":setup.loc[n,"P_bar"]}
        P = int(PT["P"])
        T = int(PT["T"])
        melt_wf['CO2'] = setup.loc[n,"CO2ppm"]/1000000.
        melt_wf["H2OT"] = setup.loc[n,"H2O"]/100.
        melt_wf['ST'] = setup.loc[n,"STppm"]/1000000.
        melt_wf['CT'] = (melt_wf['CO2']/vf.species.loc['CO2','M'])*vf.species.loc['C','M']
        melt_wf['HT'] = (melt_wf['H2OT']/vf.species.loc['H2O','M'])*(2.*vf.species.loc['H','M'])
        melt_wf['Fe3FeT'] = setup.loc[n,"Fe3+FeT"]
        fO2 = vf.f_O2(PT,melt_wf,models)
        FMQ = vf.fO22Dbuffer(PT,fO2,"FMQ",models)
        sulfsat = vf.sulfur_saturation(PT,melt_wf,models)
        sulfide_capacity = vf.C_S(PT,melt_wf)
        sulfate_capacity = vf.C_SO4(PT,melt_wf)
        # result = {"SCSS":SCSS_,"StCSS":StCSS,"sulfide_sat":sulfide_sat, "SCAS":SCAS_, "StCAS":StCAS,"sulfate_sat":sulfate_sat,"ST":ST}
        results1 = pd.DataFrame([[P,T,FMQ,sulfsat["SCSS"],sulfsat["StCSS"],sulfsat["SCAS"],sulfsat["StCAS"],setup.loc[n,"MgO"],sulfide_capacity,sulfate_capacity]])
        if n == 0.:
            results_headers = pd.DataFrame([["P","T","FMQ","SCSS","StCSS","SCAS","StCAS","MgO","Csulfide","Csulfate"]])
            results = pd.concat([results_headers, results1])
        else:
            results = pd.concat([results, results1])
    results.columns = results.iloc[0]
    results = results[1:]  
    return results

##########################
### check mass balance ###
##########################

def check_mass_balance(xg,melt,melt_and_gas):
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_O = mdv.species.loc['O','M']
    M_S = mdv.species.loc['S','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_OH = mdv.species.loc['OH','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_S2 = mdv.species.loc['S2','M']
    M_SO2 = mdv.species.loc['SO2','M']
    M_SO3 = mdv.species.loc['SO3','M']
    M_H2S = mdv.species.loc['H2S','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_OCS = mdv.species.loc['OCS','M']

    mb_C = melt_and_gas['wt_C'] - (((melt_and_gas['wt_g']*(((xg['xg_CO2']+xg['xg_CO']+xg['xg_CH4']+xg['xg_OCS'])/xg['Xg_t']) - (melt['wm_CO2']/M_CO2) - (melt['wm_CH4']/M_CH4) - (melt['wm_CO']/M_CO))) + (melt['wm_CO2']/M_CO2) + (melt['wm_CH4']/M_CH4) + (melt['wm_CO']/M_CO))*M_C)
    mb_O = melt_and_gas['wt_O'] - (((melt_and_gas['wt_g']*(((2.0*xg['xg_CO2'] + xg['xg_CO'] + 2.0*xg['xg_O2'] + xg['xg_H2O'] + 2.0*xg['xg_SO2'] + xg['xg_OCS'] )/xg['Xg_t']) - (melt['wm_H2O']/M_H2O) - ((2.0*melt['wm_CO2'])/M_CO2) - (3.0*melt['wm_SO3']/M_SO3) - (melt['wm_CO']/M_CO))) + (melt['wm_H2O']/M_H2O) + ((2.0*melt['wm_CO2'])/M_CO2) + (3.0*melt['wm_SO3']/M_SO3) + (melt['wm_CO']/M_CO) + (melt_and_gas['wt_Fe']/M_Fe)*((1.5*melt['Fe32']+1.0)/(melt['Fe32']+1.0)))*M_O)
    mb_H = melt_and_gas['wt_H'] - (((melt_and_gas['wt_g']*(((xg['xg_H2O']+xg['xg_H2']+2.0*xg['xg_CH4']+xg['xg_H2S'])/xg['Xg_t']) - (melt['wm_H2O']/M_H2O) - (melt['wm_H2']/M_H2) - (2.*melt['wm_CH4']/M_CH4) - (melt['wm_H2S']/M_H2S))) + (melt['wm_H2O']/M_H2O) + (melt['wm_H2']/M_H2) + (2.*melt['wm_CH4']/M_CH4) + (melt['wm_H2S']/M_H2S))*(2.0*M_H))
    mb_S = melt_and_gas['wt_S'] - (((melt_and_gas['wt_g']*(((xg['xg_SO2']+2.0*xg['xg_S2']+xg['xg_H2S']+xg['xg_OCS'])/xg['Xg_t']) - (melt['wm_S']/M_S) - (melt['wm_SO3']/M_SO3) - (melt['wm_H2S']/M_H2S))) + (melt['wm_S']/M_S) + (melt['wm_SO3']/M_SO3) + (melt['wm_H2S']/M_H2S))*M_S)

    mb = {'C':mb_C,'O':mb_O,'H':mb_H,'S':mb_S}

    return mb

###########################################################################
### P given S content of melt after degassing given conditions of pvsat ### IN PROGRESS
###########################################################################
