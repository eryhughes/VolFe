# batch_calculations.py

import pandas as pd
from datetime import date
import gmpy2 as gp
import numpy as np
import datetime
import math as math

import melt_gas as mg
import equilibrium_equations as eq
import model_dependent_variables as mdv
import calculations as c

# calculate the saturation pressure for multiple melt compositions in input file
def P_sat_output(first_row,last_row,p_tol,nr_step,nr_tol,setup,species,models):
    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulphide solubility","sulphate solubility","sulphide saturation","ideal gas","carbonylsulphide","mass_volume","insolubles","Saturation calculation","species X","species X solubility","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulphide","option"],models.loc["sulphate","option"],models.loc["sulphur_saturation","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulphide","option"],models.loc["mass_volume","option"],models.loc['insolubles','option'],models.loc['calc_sat','option'],models.loc['species X','option'],models.loc['species X solubility','option'],date.today()]])
    results = results.append(results1, ignore_index=True)
    results1 = ([["Sample","Saturation pressure (bar)","T ('C)","fO2 (DNNO)","fO2 (DFMQ)",
                  "SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)",
                "H2OT (wt%)","CO2T (ppm)","ST (ppm)","X (ppm)","Fe3/FeT","CT input (ppm)","HT input (ppm)","H2OT (wt%)","OH (wt%)","H2Omol (wt%)","H2 (ppm)","CH4 (ppm)","CO2T (ppm)","CO (ppm)","S2- (ppm)","S6+ (ppm)","H2S (ppm)",
                  "H_H2O/HT", "H_H2/HT", "H_CH4/HT", "H_H2S/HT", "C_CO2T/CT", "C_CO/CT", "C_CH4/CT", "S2-/ST", "S6+/ST", "H2S/ST",
                "SCSS (ppm)","sulphide saturated","SCAS (ppm)","anhydrite saturated","S melt (ppm)","graphite saturated",
                "fO2","fH2","fH2O","fS2","fSO2","fH2S","fCO2","fCO","fCH4","fOCS","fX",
                "pO2","pH2","pH2O","pS2","pSO2","pH2S","pCO2","pCO","pCH4","pOCS","pX",
                "xgO2","xgH2","xgH2O","xgS2","xgSO2","xgH2S","xgCO2","xgCO","xgCH4","xgOCS","xgX",
                 "Pvsat (H2O CO2 only)", "xg_H2O (H2O CO2 only)", "xg_CO2 (H2O CO2 only)","f_H2O (H2O CO2 only)", "f_CO2 (H2O CO2 only)","p_H2O (H2O CO2 only)", "p_CO2 (H2O CO2 only)", "Pvsat diff"]])
    results = results.append(results1, ignore_index=True)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        PT={"T":setup.loc[run,"T_C"]}
        melt_wf = {'CO2':setup.loc[run,"CO2ppm"]/1000000.,"H2OT":setup.loc[run,"H2O"]/100.}
        
        # calculate Pvsat assuming only H2O CO2 in vapour and melt
        if melt_wf['CO2'] == 0 and melt_wf['H2OT'] == 0:
            P_sat_H2O_CO2_only, xg_H2O_H2O_CO2_only, xg_CO2_H2O_CO2_only, f_H2O_H2O_CO2_only, f_CO2_H2O_CO2_only, p_H2O_H2O_CO2_only, p_CO2_H2O_CO2_only = 0., 0., 0., 0., 0., 0., 0.
        else:
            if setup.loc[run,"Fe3FeT"] > 0.:
                melt_wf['Fe3FeT'] = setup.loc[run,"Fe3FeT"]
            else:
                melt_wf['Fe3FeT'] = 0.
            P_sat_H2O_CO2_only, xg_H2O_H2O_CO2_only, xg_CO2_H2O_CO2_only, f_H2O_H2O_CO2_only, f_CO2_H2O_CO2_only, p_H2O_H2O_CO2_only, p_CO2_H2O_CO2_only, wm_H2O_H2O_CO2_only, wm_CO2_H2O_CO2_only  = c.P_sat_H2O_CO2(run,PT,melt_wf,setup,species,models,p_tol,nr_step,nr_tol)
        
        if models.loc["calc_sat","option"] == "fO2_fX":
            P_sat_, wm_ST, fSO2, wm_S2m = c.P_sat_fO2_fS2(run,PT,melt_wf,setup,species,models,p_tol)
            PT["P"] = P_sat_
        else:
            wm_ST = setup.loc[run,"STppm"]/1000000.
        melt_wf['ST'] = wm_ST
        melt_wf['CT'] = (melt_wf['CO2']/species.loc['CO2','M'])*species.loc['C','M']
        melt_wf['HT'] = (melt_wf['H2OT']/species.loc['H2O','M'])*(2.*species.loc['H','M'])
        wm_X = setup.loc[run,"Xppm"]/1000000.
        melt_wf['XT'] = wm_X
        if setup.loc[run,"S6ST"] > 0.:
            melt_wf["S6ST"] = setup.loc[run,"S6ST"]
        if models.loc["bulk_composition","option"] == "yes":
            bulk_wf = {"H":(2.*species.loc["H","M"]*melt_wf["H2OT"])/species.loc["H2O","M"],"C":(species.loc["C","M"]*melt_wf["CO2"])/species.loc["CO2","M"],"S":wm_ST, "X":wm_X}
        else:
            print("This needs fixing")
        if models.loc["sulphur_is_sat","option"] == "yes":
            if melt_wf["XT"] > 0.:
                print("this needs fixing")
            P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, Fe3FeT  = c.fO2_P_VSA(run,PT,melt_wf,setup,species,models,nr_step,nr_tol,p_tol)
        elif models.loc["sulphur_saturation","option"] == "no":
            P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST = c.P_sat(run,PT,melt_wf,setup,species,models,p_tol,nr_step,nr_tol)
        elif models.loc["sulphur_saturation","option"] == "yes":
            if melt_wf["XT"] > 0.:
                print("this needs fixing")
            P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST = c.P_VSA(run,PT,melt_wf,setup,species,models,nr_step,nr_tol,p_tol)
        PT["P"] = P_sat_
        melt_wf["H2OT"] = wm_H2O_
        melt_wf["CO2"] = wm_CO2_
        melt_wf["S2-"] = wm_S2m_
        if models.loc["sulphur_is_sat","option"] == "yes":
            melt_wf["Fe3FeT"] = Fe3FeT
        else:
            melt_wf["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
        SCSS_,sulphide_sat,SCAS_,sulphate_sat, ss_ST = c.sulphur_saturation(run,PT,melt_wf,setup,species,models)
        graphite_sat = c.graphite_saturation(run,PT,melt_wf,setup,species,models)
        gas_mf = {"O2":mg.xg_O2(run,PT,melt_wf,setup,species,models),"CO":mg.xg_CO(run,PT,melt_wf,setup,species,models),"CO2":mg.xg_CO2(run,PT,melt_wf,setup,species,models),"H2":mg.xg_H2(run,PT,melt_wf,setup,species,models),"H2O":mg.xg_H2O(run,PT,melt_wf,setup,species,models),"CH4":mg.xg_CH4(run,PT,melt_wf,setup,species,models),"S2":mg.xg_S2(run,PT,melt_wf,setup,species,models),"SO2":mg.xg_SO2(run,PT,melt_wf,setup,species,models),"H2S":mg.xg_H2S(run,PT,melt_wf,setup,species,models),"OCS":mg.xg_OCS(run,PT,melt_wf,setup,species,models),"X":mg.xg_X(run,PT,melt_wf,setup,species,models),"Xg_t":mg.Xg_tot(run,PT,melt_wf,setup,species,models),"wt_g":0.}
        #bulk_wf = {"Wt":setup.loc[run,"total_mass_g"]}

        wm_H2Omol, wm_OH = mg.wm_H2Omol_OH(run,PT,melt_wf,setup,species,models) # wt% - not giving correct answer atm

        # mass, density, volume
        #tot_m, tot_v, tot_rho, melt_m, melt_v, melt_rho, gas_m, gas_v, gas_rho = mass_vol_rho(run,PT,melt_wf,gas_mf,bulk_wf,setup,species,models)
        
        ### store results ###
        results2 = pd.DataFrame([[setup.loc[run,"Sample"],PT["P"],
                setup.loc[run,"T_C"],mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"NNO",models),mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"FMQ",models),setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],
                setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],wm_ST*1000000.,wm_X*1000000.,melt_wf["Fe3FeT"],melt_wf['CT']*1000000.,melt_wf['HT']*1000000.,wm_H2O_*100.,
wm_OH,wm_H2Omol,wm_H2_*1000000.,wm_CH4_*1000000.,wm_CO2_*1000000.,wm_CO_*1000000.,melt_wf['S2-']*1000000.,wm_S6p_*1000000.,wm_H2S_*1000000.,
                H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S2m_ST, S6p_ST, H2S_ST,
                SCSS_,sulphide_sat,SCAS_,sulphate_sat,ss_ST,graphite_sat,
                mdv.f_O2(run,PT,melt_wf,setup,species,models),mg.f_H2(run,PT,melt_wf,setup,species,models),mg.f_H2O(run,PT,melt_wf,setup,species,models),mg.f_S2(run,PT,melt_wf,setup,species,models),mg.f_SO2(run,PT,melt_wf,setup,species,models),mg.f_H2S(run,PT,melt_wf,setup,species,models),mg.f_CO2(run,PT,melt_wf,setup,species,models),mg.f_CO(run,PT,melt_wf,setup,species,models),mg.f_CH4(run,PT,melt_wf,setup,species,models),mg.f_OCS(run,PT,melt_wf,setup,species,models),mg.f_X(run,PT,melt_wf,setup,species,models),
                #mg.y_O2(PT,models),mg.y_H2(PT,species,models),mg.y_H2O(PT,species,models),mg.y_S2(PT,models),mg.y_SO2(PT,models),mdv.y_H2S(PT,models),mg.y_CO2(PT,models),mg.y_CO(PT,models),mg.y_CH4(PT,models),mg.y_OCS(PT,models),
                mg.p_O2(run,PT,melt_wf,setup,species,models),mg.p_H2(run,PT,melt_wf,setup,species,models),mg.p_H2O(run,PT,melt_wf,setup,species,models),mg.p_S2(run,PT,melt_wf,setup,species,models),mg.p_SO2(run,PT,melt_wf,setup,species,models),mg.p_H2S(run,PT,melt_wf,setup,species,models),mg.p_CO2(run,PT,melt_wf,setup,species,models),mg.p_CO(run,PT,melt_wf,setup,species,models),mg.p_CH4(run,PT,melt_wf,setup,species,models),mg.p_OCS(run,PT,melt_wf,setup,species,models),mg.p_X(run,PT,melt_wf,setup,species,models),
                mg.xg_O2(run,PT,melt_wf,setup,species,models),mg.xg_H2(run,PT,melt_wf,setup,species,models),mg.xg_H2O(run,PT,melt_wf,setup,species,models),mg.xg_S2(run,PT,melt_wf,setup,species,models),mg.xg_SO2(run,PT,melt_wf,setup,species,models),mg.xg_H2S(run,PT,melt_wf,setup,species,models),mg.xg_CO2(run,PT,melt_wf,setup,species,models),mg.xg_CO(run,PT,melt_wf,setup,species,models),mg.xg_CH4(run,PT,melt_wf,setup,species,models),mg.xg_OCS(run,PT,melt_wf,setup,species,models),mg.xg_X(run,PT,melt_wf,setup,species,models),
                 P_sat_H2O_CO2_only, xg_H2O_H2O_CO2_only, xg_CO2_H2O_CO2_only,f_H2O_H2O_CO2_only, f_CO2_H2O_CO2_only,p_H2O_H2O_CO2_only, p_CO2_H2O_CO2_only,P_sat_H2O_CO2_only-PT["P"]]])
                             
        results = results.append(results2, ignore_index=True)
        results.to_csv('saturation_pressures.csv', index=False, header=False)
        print(n, setup.loc[run,"Sample"],PT["P"])

        

###############
### gassing ###
###############

def gassing(run,gassing_inputs,setup,species,models):
    
    print(setup.loc[run,"Sample"])
     
    if models.loc["fO2","option"] != "Kress91A":
        print("stop - model requires Kress91A to be used")
        return
    if models.loc["fO2","option"] == "regas" and models.loc["starting_P","option"] == "bulk":
        print("not possible - will calculate negative gas weight fractions")
        return
    
    # set T, volatile composition of the melt, and tolerances
    PT={"T":setup.loc[run,"T_C"]}
    melt_wf = {'CO2':setup.loc[run,"CO2ppm"]/1000000.,"H2OT":setup.loc[run,"H2O"]/100.,"ST":setup.loc[run,"STppm"]/1000000.,"H2":0.,"XT":setup.loc[run,"Xppm"]/1000000.}
    melt_wf["CT"] = (melt_wf["CO2"]/species.loc["CO2","M"])*species.loc["C","M"]
    melt_wf["HT"] = (2.*melt_wf["H2OT"]/species.loc["H2O","M"])*species.loc["H","M"]
    melt_wf["ST"] = melt_wf["ST"]
    nr_step = gassing_inputs["nr_step"]
    nr_tol = gassing_inputs["nr_tol"]
    dp_step = gassing_inputs["dp_step"]
    psat_tol = gassing_inputs["psat_tol"]
    dwtg = gassing_inputs["dwtg"]
    i_nr_step = gassing_inputs["i_nr_step"]
    i_nr_tol = gassing_inputs["i_nr_tol"]

    # Calculate saturation pressure for composition given in inputs file
    if models.loc["insolubles","option"] == "H2O-CO2 only":  
        P_sat_, xg_H2O_, xg_CO2_, f_H2O_, f_CO2_, p_H2O_, p_CO2_, wm_H2O_, wm_CO2_ = c.P_sat_H2O_CO2(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
        wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, wm_SO3_, wm_X = 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.
    else:
        P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST = c.P_sat(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
        wm_SO3_ = (wm_S6p_*species.loc["SO3","M"])/species.loc["S","M"]
    
    PT["P"] = P_sat_
    print(PT["T"],PT["P"],datetime.datetime.now())

    # update melt composition at saturation pressure and check for sulphur saturation
    melt_wf["H2OT"] = wm_H2O_
    melt_wf["CO2"] = wm_CO2_
    melt_wf["CO"] = wm_CO_
    melt_wf["CH4"] = wm_CH4_
    melt_wf["H2"] = wm_H2_
    melt_wf["S2-"] = wm_S2m_
    melt_wf["S6+"] = wm_S6p_
    melt_wf["H2S"] = wm_H2S_
    melt_wf["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
    melt_wf["S6ST"] = mg.S6ST(run,PT,melt_wf,setup,species,models)
    #wm_H2Omol_, wm_OH_ = mg.wm_H2Omol_OH(run,PT,melt_wf,setup,species,models) #wt%
    #melt_wf["H2Omol"] = wm_H2Omol_/100.
    #melt_wf["OH"] = wm_OH_/100.
    SCSS_,sulphide_sat,SCAS_,sulphate_sat,ST_ = c.sulphur_saturation(run,PT,melt_wf,setup,species,models)
    
    xg_CS = (mg.xg_CO(run,PT,melt_wf,setup,species,models)+mg.xg_CO2(run,PT,melt_wf,setup,species,models)+mg.xg_OCS(run,PT,melt_wf,setup,species,models)+mg.xg_CH4(run,PT,melt_wf,setup,species,models))/(mg.xg_S2(run,PT,melt_wf,setup,species,models)+mg.xg_SO2(run,PT,melt_wf,setup,species,models)+mg.xg_OCS(run,PT,melt_wf,setup,species,models)+mg.xg_H2S(run,PT,melt_wf,setup,species,models))
    wm_CO2eq = wm_CO2_ + species.loc["CO2","M"]*((wm_CO_/species.loc["CO","M"])+(wm_CH4_/species.loc["CH4","M"]))
    wm_H2Oeq = wm_H2O_ + species.loc["H2O","M"]*((wm_H2_/species.loc["H2","M"])+(2.*wm_CH4_/species.loc["CH4","M"])+(wm_H2S_/species.loc["H2S","M"]))
    
    # Set bulk composition
    wt_C, wt_O, wt_H, wt_S, wt_X, wt_Fe, wt_g_, Wt_ = c.bulk_composition(run,PT,melt_wf,setup,species,models)
    bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"Fe":wt_Fe,"Wt":Wt_,"X":wt_X}
        
    # set system and initial guesses
    system = eq.set_system(melt_wf,models)
    guessx, guessy, guessz, guessw = eq.initial_guesses(run,PT,melt_wf,setup,species,models,system)
    
    # calculating equilibrium isotopic fractionation?
    if models.loc["isotopes","option"] == "yes":
        if setup.loc[run,"bulk 34S/32S"] > 0.:
            R_i = {"S": setup.loc[run,"bulk 34S/32S"]}
        elif setup.loc[run,"bulk d34S"] >= 0.:
            d = setup.loc[run,"bulk d34S"]
            R_i = {"S": delta2ratio("VCDT",d)}
        elif setup.loc[run,"bulk d34S"] <= 0.:
            d = setup.loc[run,"bulk d34S"]
            R_i = {"S": delta2ratio("VCDT",setup.loc[run,"bulk d34S"])}
        else:
            print("input bulk isotope value for system or change 'isotopes' to 'no'")
            return
        R_S_SO4_, R_S_S2_ = iso.i2s2("S",PT,R_i,melt_wf)
    else:
        R_i = {"S":""}
        
    # create results table
    results_header1 = pd.DataFrame([["System","run","Sample","H2O (mwf)","CO2 (mwf)","ST (mwf)","X (mwf)","Fe3FeT","S6ST","O (twf)","C (twf)","H (wtf)","S (twf)","Fe (twf)","Saturation P (bars)",
"SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)","34S/32S",
"oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulphide solubility","sulphate solubility","ideal gas","carbonylsulphide","bulk composition","equilibrate Fe","starting pressure","gassing direction","gassing style","mass_volume","crystallisation","isotopes","Date"]])
    results_header2 = pd.DataFrame([[system,run,setup.loc[run,"Sample"],melt_wf["H2OT"],melt_wf["CO2"],melt_wf["ST"],melt_wf["XT"],melt_wf["Fe3FeT"],"SORT",bulk_wf["O"],bulk_wf["C"],bulk_wf["H"],bulk_wf["S"],bulk_wf["Fe"],P_sat_,
setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],R_i["S"],
models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulphide","option"],models.loc["sulphate","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulphide","option"],models.loc["bulk_composition","option"],models.loc["eq_Fe","option"],models.loc["starting_P","option"],models.loc["gassing_direction","option"],models.loc["gassing_style","option"],models.loc["mass_volume","option"],models.loc["crystallisation","option"],models.loc["isotopes","option"],date.today()]])
    results_header = results_header1.append(results_header2, ignore_index=True)
    results_chemistry1 = pd.DataFrame([["P","T('C)","xg_O2","xg_CO","xg_CO2","xg_H2","xg_H2O","xg_CH4","xg_S2","xg_SO2","xg_H2S","xg_OCS","xg_X","Xg_t","xg_CS",
               "xm_CO2","xm_H2O","Xm_t_SO","Xm_t_ox",
               "wm_CO2-eq","wm_H2O-eq","wm_CO2T","wm_H2OT","wm_H2","wm_CO","wm_CH4","wm_S","wm_SO3","wm_H2S","wm_ST","wm_X","Fe32","Fe3T","S62","S6T",
               "DFMQ","DNNO","SCSS","sulphide sat?","SCAS","sulphate sat?","wt_g","wt_g_O","wt_g_C","wt_g_H","wt_g_S","wt_g_X","wt_O","wt_C","wt_H","wt_S","wt_X",
               "fO2","fH2","fH2O","fS2","fSO2","fH2S","fCO2","fCO","fCH4","fOCS","fX",
               "yO2","yH2","yH2O","yS2","ySO2","yH2S","yCO2","yCO","yCH4","yOCS","yX",
               "M_m_SO","M_m_ox","C_H2O","C_H2","C_CO3","C_CO","C_CH4","C_S","C_SO4","C_H2S","C_X",
               "KD1","KHOg","KHOm","KHOSg","KCOg","KCOHg","KOCSg","KSOg","KSOg2",
               "Solving species","tot_mass (g)","tot_vol (cm3)","tot_rho (g/cm3)","melt_mass (g)","melt_vol (cm3)","melt_rho (g/cm3)","gas_mass (g)","gas_vol (cm3)","gas_rho (g/cm3)"]])
    results_chemistry = results_header.append(results_chemistry1, ignore_index=True)
    results1 = pd.DataFrame([[PT["P"],PT["T"],mg.xg_O2(run,PT,melt_wf,setup,species,models),mg.xg_CO(run,PT,melt_wf,setup,species,models),mg.xg_CO2(run,PT,melt_wf,setup,species,models),mg.xg_H2(run,PT,melt_wf,setup,species,models),mg.xg_H2O(run,PT,melt_wf,setup,species,models),mg.xg_CH4(run,PT,melt_wf,setup,species,models),mg.xg_S2(run,PT,melt_wf,setup,species,models),mg.xg_SO2(run,PT,melt_wf,setup,species,models),mg.xg_H2S(run,PT,melt_wf,setup,species,models),mg.xg_OCS(run,PT,melt_wf,setup,species,models),mg.xg_X(run,PT,melt_wf,setup,species,models),mg.Xg_tot(run,PT,melt_wf,setup,species,models),xg_CS,
mg.xm_CO2_so(run,melt_wf,setup,species),mg.xm_H2OT_so(run,melt_wf,setup,species),mg.Xm_t_so(run,melt_wf,setup,species),"",
wm_CO2eq,wm_H2Oeq,melt_wf["CO2"],melt_wf["H2OT"],wm_H2_,wm_CO_,wm_CH4_,wm_S2m_,wm_SO3_,wm_H2S_,melt_wf["ST"],melt_wf["XT"],mg.Fe3Fe2(melt_wf),melt_wf["Fe3FeT"],mg.S6S2(run,PT,melt_wf,setup,species,models),mg.S6ST(run,PT,melt_wf,setup,species,models),
mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"FMQ",models),mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"NNO",models),SCSS_,sulphide_sat,SCAS_,sulphate_sat,wt_g_,"","","","","",
bulk_wf["O"],bulk_wf["C"],bulk_wf["H"],bulk_wf["S"],bulk_wf["X"],
mdv.f_O2(run,PT,melt_wf,setup,species,models),mg.f_H2(run,PT,melt_wf,setup,species,models),mg.f_H2O(run,PT,melt_wf,setup,species,models),mg.f_S2(run,PT,melt_wf,setup,species,models),mg.f_SO2(run,PT,melt_wf,setup,species,models),mg.f_H2S(run,PT,melt_wf,setup,species,models),mg.f_CO2(run,PT,melt_wf,setup,species,models),mg.f_CO(run,PT,melt_wf,setup,species,models),mg.f_CH4(run,PT,melt_wf,setup,species,models),mg.f_OCS(run,PT,melt_wf,setup,species,models),mg.f_X(run,PT,melt_wf,setup,species,models),
mdv.y_O2(PT,species,models),mdv.y_H2(PT,species,models),mdv.y_H2O(PT,species,models),mdv.y_S2(PT,species,models),mdv.y_SO2(PT,species,models),mdv.y_H2S(PT,species,models),mdv.y_CO2(PT,species,models),mdv.y_CO(PT,species,models),mdv.y_CH4(PT,species,models),mdv.y_OCS(PT,species,models),mdv.y_X(PT,species,models),
mg.M_m_SO(run,melt_wf,setup,species),mg.M_m_ox(run,melt_wf,setup,species,models),mdv.C_H2O(run,PT,melt_wf,setup,species,models),mdv.C_H2(run,PT,melt_wf,setup,species,models),mdv.C_CO3(run,PT,melt_wf,setup,species,models),mdv.C_CO(run,PT,melt_wf,setup,species,models),mdv.C_CH4(run,PT,melt_wf,setup,species,models),mdv.C_S(run,PT,melt_wf,setup,species,models),mdv.C_SO4(run,PT,melt_wf,setup,species,models),mdv.C_H2S(run,PT,melt_wf,setup,species,models),mdv.C_X(run,PT,melt_wf,setup,species,models),
"",mdv.KHOg(PT,models),mdv.KHOm(run,PT,melt_wf,setup,species,models),mdv.KHOSg(PT,models),mdv.KCOg(PT,models),mdv.KCOHg(PT,models),mdv.KOCSg(PT,models),mdv.KOSg(PT,models),mdv.KOSg2(PT,models),"na","Fix tot_mass (g)","Fix tot_vol (cm3)","Fix tot_rho (g/cm3)","Fix melt_mass (g)","Fix melt_vol (cm3)","Fix melt_rho (g/cm3)",0.,0.,""]])
    results_chemistry = results_chemistry.append(results1, ignore_index=True)
    results_chemistry.to_csv('results_gassing_chemistry.csv', index=False, header=False)
    
    # results for isotope calculations...
    if models.loc["isotopes","option"] == "yes":
        a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_ = iso.i2s6_S_alphas(PT)
        results_isotopes1 = pd.DataFrame([["P","T_C","xg_O2","xg_CO","xg_CO2","xg_H2","xg_H2O","xg_CH4","xg_S2","xg_SO2","xg_SO3","xg_H2S","xg_OCS","wt_g",
                        "wm_CO2","wm_H2O","wm_H2","wm_S","wm_SO3","wm_ST","Fe3T","S6T",
                         "DFMQ","DNNO","SCSS","sulphide sat?","SCAS","sulphate sat?",
                         "RS_S2-","RS_SO42-","RS_H2S","RS_SO2","RS_S2","R_OCS","R_m","R_g",
                         "dS_S2-","dS_SO42-","dS_H2S","dS_SO2","dS_S2","dS_OCS","dS_m","dS_g",
                         "a_H2S_S2-","a_SO42-_S2-","a_S2_S2-","a_SO2_S2-","a_OCS_S2-","a_g_m"]])
        results_isotopes = results_header.append(results_isotopes1, ignore_index=True)
        results1 = pd.DataFrame([[PT["P"],PT["T"],mg.xg_O2(run,PT,melt_wf,setup,species,models),mg.xg_CO(run,PT,melt_wf,setup,species,models),mg.xg_CO2(run,PT,melt_wf,setup,species,models),mg.xg_H2(run,PT,melt_wf,setup,species,models),mg.xg_H2O(run,PT,melt_wf,setup,species,models),mg.xg_CH4(run,PT,melt_wf,setup,species,models),mg.xg_S2(run,PT,melt_wf,setup,species,models),mg.xg_SO2(run,PT,melt_wf,setup,species,models),mg.xg_H2S(run,PT,melt_wf,setup,species,models),mg.xg_OCS(run,PT,melt_wf,setup,species,models),wt_g_,
melt_wf["CO2"],melt_wf["H2OT"],0,(mg.wm_S(run,PT,melt_wf,setup,species,models)/100),(mg.wm_SO3(run,PT,melt_wf,setup,species,models)/100),melt_wf["ST"],melt_wf["Fe3FeT"],mg.S6ST(run,PT,melt_wf,setup,species,models),
mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"FMQ"),mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"NNO"),SCSS_,sulphide_sat,SCAS_,sulphate_sat,
R_S_S2_,R_S_SO4_,"","","","",R_i["S"],"",ratio2delta("VCDT",R_S_S2_),ratio2delta("VCDT",R_S_SO4_),"","","","",ratio2delta("VCDT",R_i["S"]),"",
a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_,""]])
        results_isotopes = results_isotopes.append(results1, ignore_index= True)
        results_isotopes.to_csv('results_gassing_isotopes.csv', index=False, header=False)
    
    if models.loc["P_variation","option"] == "polybaric":
        # pressure ranges and options
        starting_P = models.loc["starting_P","option"]
        if starting_P == "bulk":
            if models.loc["bulk_composition","option"] != "yes":
                print("not possible - starting pressure is the Pvsat for the melt")
        if starting_P == "set":
            initial = int(setup.loc[run,"P_bar"])
        else:
            if models.loc["gassing_direction","option"] == "degas":
                answer = PT["P"]/dp_step
                answer = round(answer)
                initial = round(answer*dp_step)
            elif models.loc["gassing_direction","option"] == "regas":
                initial = round(PT["P"])+1 
        if models.loc["gassing_direction","option"] == "degas":
            step = int(-1*dp_step) # pressure step in bars
            final = 0
        elif models.loc["gassing_direction","option"] == "regas":
            step = int(dp_step)
            final = int(setup.loc[run,"final_P"])
    elif models.loc["T_variation","option"] == "polythermal": # temperature ranges and options
        PT["P"] = setup.loc[run,"P_bar"]
        final = int(setup.loc[run,"final_T"])
        if setup.loc[run,"final_T"] > setup.loc[run,"T_C"]:
            initial = int(round(PT["T"])) 
            step = int(dp_step) # temperature step in 'C
        elif setup.loc[run,"final_T"] < setup.loc[run,"T_C"]:
            initial = int(round(PT["T"]))
            step = int(-1.*dp_step) # temperature step in 'C
    
    # add some gas to the system if doing open-system regassing
    if models.loc["gassing_direction","option"] == "regas" and models.loc["gassing_style","option"] == "open":
        gas_mf = {"O2":mg.xg_O2(run,PT,melt_wf,setup,species,models),"CO":mg.xg_CO(run,PT,melt_wf,setup,species,models),"CO2":mg.xg_CO2(run,PT,melt_wf,setup,species,models),"H2":mg.xg_H2(run,PT,melt_wf,setup,species,models),"H2O":mg.xg_H2O(run,PT,melt_wf,setup,species,models),"CH4":mg.xg_CH4(run,PT,melt_wf,setup,species,models),"S2":mg.xg_S2(run,PT,melt_wf,setup,species,models),"SO2":mg.xg_SO2(run,PT,melt_wf,setup,species,models),"SO3":mg.xg_SO3(run,PT,melt_wf,setup,species,models),"H2S":mg.xg_H2S(run,PT,melt_wf,setup,species,models),"OCS":mg.xg_OCS(run,PT,melt_wf,setup,species,models),"X":mg.xg_X(run,PT,melt_wf,setup,species,models),"Xg_t":mg.Xg_tot(run,PT,melt_wf,setup,species,models),"wt_g":0.}
        wt_C, wt_H, wt_S, wt_X, wt_Fe, wt_O, Wt = new_bulk_regas_open(run,PT,melt_wf,bulk_wf,gas_mf,dwtg,setup,species,models)
        bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"Fe":wt_Fe,"X":wt_X,"Wt":Wt}
    
    # run over different pressures #
    number_of_step = 0.
    for i in range(initial,final,step): # P is pressure in bars or T is temperature in 'C
        number_of_step = number_of_step + 1.
        eq_Fe = models.loc["eq_Fe","option"]
        
        if models.loc["gassing_style","option"] == "open": # check melt is still vapor-saturated
            if models.loc["insolubles","option"] == "H2O-CO2 only":  
                P_sat_, xg_H2O_, xg_CO2_, f_H2O_, f_CO2_, p_H2O_, p_CO2_, wm_H2O_, wm_CO2_ = c.P_sat_H2O_CO2(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
                wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST, wm_SO3_ = 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.
            else:
                P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST = c.P_sat(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
                wm_SO3_ = (wm_S6p_*species.loc["SO3","M"])/species.loc["S","M"]
        
        if models.loc["P_variation","option"] == "polybaric": 
            P = i - dp_step
            if P < dp_step:
                P = 1.
            PT["P"] = P
        elif models.loc["T_variation","option"] == "polythermal":
            T = i - dp_step
            PT["T"] = T
        if P_sat_ > PT["P"]:  
            # work out equilibrium partitioning between melt and gas phase
            xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, xg_X_, Xg_t, xm_H2O_, xm_CO2_, Xm_t, wm_H2O_, wm_H2_, wm_CO2_, wm_CO_, wm_CH4_, wm_S_, wm_SO3_, wm_H2S_, wm_ST_, wm_X_, Fe32, Fe3T, S62, S6T, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g_X, wt_g, wt_O_, wt_C_, wt_H_, wt_S_, wt_X_, guessx, guessy, guessz, guessw = eq.mg_equilibrium(run,PT,melt_wf,bulk_wf,setup,species,models,nr_step,nr_tol,guessx,guessy,guessz,guessw)
            # gas composition
            gas_mf = {"O2":xg_O2_,"CO":xg_CO_,"S2":xg_S2_,"CO2":xg_CO2_,"H2O":xg_H2O_,"H2":xg_H2_,"CH4":xg_CH4_,"SO2":xg_SO2_,"H2S":xg_H2S_,"OCS":xg_OCS_,"X":xg_X_,"Xg_t":Xg_t,"wt_g":wt_g}
            
            if models.loc["insolubles","option"] == "H2O-CO2 only":
                Fe3T = melt_wf["Fe3FeT"]
                Fe32 = mg.overtotal2ratio(Fe3T)
        else:
            xm_H2O_, wm_H2O_, xm_CO2_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_H2S_, wm_S_, wm_S6p_, H2O_HT, H2_HT, CH4_HT, CO2_CT, CO_CT, CH4_CT, S6T, S2m_ST, H2S_ST, H2S_HT = eq.melt_speciation(run,PT,melt_wf,setup,species,models,nr_step,nr_tol)
            wm_ST_ = wm_S_ + wm_S6p_
            S62 = S6T/S2m_ST
            Fe3T = melt_wf["Fe3FeT"]
            Fe32 = mg.overtotal2ratio(Fe3T)
            xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, xg_X, Xg_t, Xm_t, Xm_t_ox, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g = "","","","","","","","","","","","","","",0.,0.,0.,0.,0.
            guessx, guessy, guessz, guessw = eq.initial_guesses(run,PT,melt_wf,setup,species,models,system)
        
        # set melt composition for forward calculation
        melt_wf["CO2"] = wm_CO2_
        melt_wf["H2OT"] = wm_H2O_
        melt_wf["H2"] = wm_H2_
        melt_wf["CO"] = wm_CO_
        melt_wf["CH4"] = wm_CH4_
        melt_wf["H2S"] = wm_H2S_
        melt_wf["S6+"] = wm_S6p_
        melt_wf["S2-"] = wm_S_
        melt_wf["ST"] = wm_ST_
        melt_wf["XT"] = wm_X_
        melt_wf["Fe3FeT"] = Fe3T
        if P_sat_ < PT["P"]:  
            wt_C_, wt_O_, wt_H_, wt_S_, wt_X_, wt_Fe, wt_g_, Wt_ = c.bulk_composition(run,PT,melt_wf,setup,species,models)
    
        # check for sulphur saturation and display warning in outputs
        SCSS_,sulphide_sat,SCAS_,sulphate_sat, ST_ = c.sulphur_saturation(run,PT,melt_wf,setup,species,models)
        if sulphide_sat == "yes":
            warning = "WARNING: sulphide-saturated"
        elif sulphate_sat == "yes":
            warning = "WARNING: sulphate-saturated"
        else:
            warning = ""
        
        # calculate fO2
        if eq_Fe == "yes":
            fO2_ = mdv.f_O2(run,PT,melt_wf,setup,species,models)
        elif eq_Fe == "no":
            fO2_ = (xg_O2_*mdv.y_O2(PT,models)*PT["P"])
        
        tot_m, tot_v, tot_rho, melt_m, melt_v, melt_rho, gas_m, gas_v, gas_rho = "","","","","","","","",""
        if system in ["SCHOFe"]:
            xg_CS = (xg_CO_+xg_CO2_+xg_OCS_+xg_CH4_)/(xg_S2_+xg_SO2_+xg_OCS_+xg_H2S_)
        else:
            xg_CS = ""
        wm_CO2eq = wm_CO2_ + species.loc["CO2","M"]*((wm_CO_/species.loc["CO","M"])+(wm_CH4_/species.loc["CH4","M"]))
        wm_H2Oeq = wm_H2O_ + species.loc["H2O","M"]*((wm_H2_/species.loc["H2","M"])+(2.*wm_CH4_/species.loc["CH4","M"])+(wm_H2S_/species.loc["H2S","M"]))
    
        # store results
        results2 = pd.DataFrame([[PT["P"],PT["T"],xg_O2_,xg_CO_,xg_CO2_,xg_H2_,xg_H2O_,xg_CH4_,xg_S2_,xg_SO2_,xg_H2S_,xg_OCS_,xg_X_,Xg_t,xg_CS,
               xm_CO2_,xm_H2O_,Xm_t,"",
               wm_CO2eq,wm_H2Oeq,wm_CO2_,wm_H2O_,wm_H2_,wm_CO_,wm_CH4_,wm_S_,wm_SO3_,wm_H2S_,wm_ST_,wm_X_,Fe32,Fe3T,S62,S6T,
               mg.fO22Dbuffer(PT,fO2_,"FMQ",models),mg.fO22Dbuffer(PT,fO2_,"NNO",models),SCSS_,sulphide_sat,SCAS_,sulphate_sat,
               wt_g,wt_g_O,wt_g_C,wt_g_H,wt_g_S,wt_g_X,wt_O_,wt_C_,wt_H_,wt_S_,wt_X_,
fO2_,mg.f_H2(run,PT,melt_wf,setup,species,models),mg.f_H2O(run,PT,melt_wf,setup,species,models),mg.f_S2(run,PT,melt_wf,setup,species,models),mg.f_SO2(run,PT,melt_wf,setup,species,models),mg.f_H2S(run,PT,melt_wf,setup,species,models),mg.f_CO2(run,PT,melt_wf,setup,species,models),mg.f_CO(run,PT,melt_wf,setup,species,models),mg.f_CH4(run,PT,melt_wf,setup,species,models),mg.f_OCS(run,PT,melt_wf,setup,species,models),mg.f_X(run,PT,melt_wf,setup,species,models),
mdv.y_O2(PT,species,models),mdv.y_H2(PT,species,models),mdv.y_H2O(PT,species,models),mdv.y_S2(PT,species,models),mdv.y_SO2(PT,species,models),mdv.y_H2S(PT,species,models),mdv.y_CO2(PT,species,models),mdv.y_CO(PT,species,models),mdv.y_CH4(PT,species,models),mdv.y_OCS(PT,species,models),mdv.y_X(PT,species,models),
mg.M_m_SO(run,melt_wf,setup,species),mg.M_m_ox(run,melt_wf,setup,species,models),mdv.C_H2O(run,PT,melt_wf,setup,species,models),mdv.C_H2(run,PT,melt_wf,setup,species,models),mdv.C_CO3(run,PT,melt_wf,setup,species,models),mdv.C_CO(run,PT,melt_wf,setup,species,models),mdv.C_CH4(run,PT,melt_wf,setup,species,models),mdv.C_S(run,PT,melt_wf,setup,species,models),mdv.C_SO4(run,PT,melt_wf,setup,species,models),mdv.C_H2S(run,PT,melt_wf,setup,species,models),mdv.C_X(run,PT,melt_wf,setup,species,models),
"",mdv.KHOg(PT,models),mdv.KHOm(run,PT,melt_wf,setup,species,models),mdv.KHOSg(PT,models),mdv.KCOg(PT,models),mdv.KCOHg(PT,models),mdv.KOCSg(PT,models),mdv.KOSg(PT,models),mdv.KOSg2(PT,models),
models.loc["solve_species","option"],tot_m,tot_v,tot_rho,melt_m,melt_v,melt_rho,gas_m,gas_v,gas_rho]])
        results_chemistry = results_chemistry.append(results2, ignore_index=True)
        results_chemistry.to_csv('results_gassing_chemistry.csv', index=False, header=False)
        
        # equilibrium isotope fractionation
        if models.loc["isotopes", "option"] == "yes":
            if models.loc["H2S","option"] == "yes":
                print("not currently possible")
            A, B = iso.i2s6("S",PT,R_i,melt_wf,gas_mf,species,i_nr_step,i_nr_tol,guessx)
            RS_Sm, RS_H2S, RS_SO4, RS_S2, RS_SO2, RS_OCS = A
            RS_m, RS_g = B
            a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_ = iso.i2s6_S_alphas(PT)
            xg_SO3_ = 0.
            results2 = pd.DataFrame([[PT["P"],PT["T"],xg_O2_,xg_CO_,xg_CO2_,xg_H2_,xg_H2O_,xg_CH4_,xg_S2_,xg_SO2_,xg_SO3_,xg_H2S_,xg_OCS_,wt_g,
                        wm_CO2_,wm_H2O_,wm_H2_,wm_S_,wm_SO3_,wm_ST_,Fe3T,S6T,
                        mg.fO22Dbuffer(PT,fO2_,"FMQ"),mg.fO22Dbuffer(PT,fO2_,"NNO"),SCSS_,sulphide_sat,SCAS_,sulphate_sat,
                        RS_Sm, RS_SO4, RS_H2S, RS_SO2, RS_S2, RS_OCS, RS_m, RS_g, ratio2delta("VCDT",RS_Sm),ratio2delta("VCDT",RS_SO4),ratio2delta("VCDT",RS_H2S),ratio2delta("VCDT",RS_SO2),ratio2delta("VCDT",RS_S2),ratio2delta("VCDT",RS_OCS),ratio2delta("VCDT",RS_m),ratio2delta("VCDT",RS_g),a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_,RS_g/RS_m]])
            results_isotopes = results_isotopes.append(results2, ignore_index=True)
            results_isotopes.to_csv('results_gassing_isotopes.csv', index=False, header=False)
        
        if(number_of_step % 100==0):
            print(PT["T"],PT["P"],mg.fO22Dbuffer(PT,fO2_,"FMQ",models),warning,datetime.datetime.now())

        # recalculate bulk composition if needed
        if models.loc["gassing_style","option"] == "open":
            wm_C_, wm_H_, wm_S_, wm_Fe_, wm_O_, wm_X_ = mg.melt_elements(run,PT,melt_wf,bulk_wf,gas_mf,setup,species,models)
            if models.loc["gassing_direction","option"] == "degas":
                Wt_ = bulk_wf['Wt']
                #if wm_C_ < 1.e-6:
                    #wm_C_ = 0.
                bulk_wf = {"C":wm_C_,"H":wm_H_,"O":wm_O_,"S":wm_S_,"X":wm_X_,"Fe":wm_Fe_,"Wt":(Wt_*(1. - wt_g))}
                melt_wf["CT"] = wm_C_
                melt_wf["HT"] = wm_H_
                melt_wf["ST"] = wm_S_
                melt_wf["XT"] = wm_X_
            elif models.loc["gassing_direction","option"] == "regas":
                wt_C, wt_H, wt_S, wt_X, wt_Fe, wt_O, Wt = c.new_bulk_regas_open(run,PT,melt_wf,bulk_wf,gas_mf,dwtg,setup,species,models)
                bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"X":wt_X,"Fe":wt_Fe,"Wt":Wt}
                melt_wf["CT"] = wm_C_
                melt_wf["HT"] = wm_H_
                melt_wf["ST"] = wm_S_
                melt_wf["XT"] = wm_X_
        if models.loc["crystallisation","option"] == "yes":
            wt_C_ = bulk_wf["C"]
            wt_H_ = bulk_wf["H"]
            wt_O_ = bulk_wf["O"]
            wt_S_ = bulk_wf["S"]
            wt_X_ = bulk_wf["X"]
            wt_Fe_ = bulk_wf["Fe"]
            wt_ = bulk_wf["Wt"]
            Xst = setup.loc[run,"crystallisation_pc"]/100.
            bulk_wf = {"C":wt_C_*(1./(1.-Xst)),"H":wt_H_*(1./(1.-Xst)),"O":wt_O_*(1./(1.-Xst)),"S":wt_S_*(1./(1.-Xst)),"X":wt_X_*(1./(1.-Xst)),"Fe":wt_Fe_*(1./(1.-Xst)),"Wt":wt_*(1.-Xst)}
            
    print("done", datetime.datetime.now())
            
def calc_isobar(run,setup,species,models,initial_P,final_P,step_P):
    
    if models.loc["insolubles","option"] == "H2O-CO2 only":
        PT={"T":setup.loc[run,"T_C"]}
        
        # set up results table
        results = pd.DataFrame([["P (bar)","H2O (wt%)","CO2 (ppm)"]])
    
        initial_P = int(initial_P)
        final_P = int(final_P+1)
        step_P = int(step_P)
        
        for n in range(initial_P,final_P,step_P):
            PT["P"] = n # pressure in bars
            melt_wf={"not needed":0.}
            results1 = c.calc_isobar_CO2H2O(run,PT,melt_wf,setup,species,models)
            results = results.append(results1, ignore_index=True)
            print(setup.loc[run,"Sample"],n)
    else:
        print("insolubles must be H2O-CO2 only")
    results.to_csv('results_isobars.csv', index=False, header=False)
    
def calc_pure_solubility(run,setup,species,models,initial_P):
    print(setup.loc[run,"Sample"],initial_P)
    PT={"T":setup.loc[run,"T_C"]}
        
    # set up results table
    results = pd.DataFrame([["P (bar)","H2O (wt%)","CO2 (ppm)"]])
    
    initial_P = int(initial_P)
        
    for n in range(initial_P,1,-1):
        PT["P"] = n # pressure in bars
        melt_wf={"not needed":0.}
        results1 = c.calc_pure_solubility(run,PT,melt_wf,setup,species,models)
        results = results.append(results1, ignore_index=True)
        
    results.to_csv('results_pure_solubility.csv', index=False, header=False)
    print("done")

    
    
##################
### Capacities ###
##################

# print capacities for multiple melt compositions in input file
def capacities_output(first_row,last_row,setup,species,models):
    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulphide solubility","sulphate solubility","sulphide saturation","ideal gas","carbonylsulphide","mass_volume","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulphide","option"],models.loc["sulphate","option"],models.loc["sulphur_saturation","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulphide","option"],models.loc["mass_volume","option"],date.today()]])
    results = results.append(results1, ignore_index=True)
    results1 = ([["Sample","Pressure (bar)","T ('C)","SiO2","TiO2","Al2O3","FeOT","MnO","MgO","CaO","Na2O","K2O","P2O5",
                "H2O","CO2 (ppm)","ST (ppm)","Fe3+/FeT","fO2 DFMQ","ln[C_CO32-]","ln[C_H2OT]","ln[C_S2-]","ln[C_S6+]","ln[C_H2S]","ln[C_H2]","ln[C_CO]","ln[C_CH4]","ln[C_X]","M_m_SO",]])
    results = results.append(results1, ignore_index=True)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        PT={"T":setup.loc[run,"T_C"]}
        CO2 = setup.loc[run,"CO2ppm"]/1000000.
        CT = (CO2/species.loc["CO2","M"])*species.loc["C","M"]
        H2O = setup.loc[run,"H2O"]/100.
        HT = (H2O/species.loc["H2O","M"])*species.loc["H2","M"]
        XT = setup.loc[run,"Xppm"]/1000000.
        ST = setup.loc[run,"STppm"]/1000000.
        melt_wf = {'CO2':CO2,"H2OT":H2O,"ST":ST,"X":XT,"HT":HT,"CT":CT}
        PT["P"] = setup.loc[run,"P_bar"]
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
        C_CO32 = mdv.C_CO3(run,PT,melt_wf,setup,species,models)
        C_H2OT = mdv.C_H2O(run,PT,melt_wf,setup,species,models)
        C_S2 = mdv.C_S(run,PT,melt_wf,setup,species,models)
        C_S6 = mdv.C_SO4(run,PT,melt_wf,setup,species,models)
        C_H2S = mdv.C_H2S(run,PT,melt_wf,setup,species,models)
        C_H2 = mdv.C_H2(run,PT,melt_wf,setup,species,models)
        C_CO = mdv.C_CO(run,PT,melt_wf,setup,species,models)
        C_CH4 = mdv.C_CH4(run,PT,melt_wf,setup,species,models)
        C_X = mdv.C_X(run,PT,melt_wf,setup,species,models)
        fO2_ = mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"FMQ",models)
        M_m = mg.M_m_SO(run,melt_wf,setup,species)
                
        ### store results ###
        results2 = pd.DataFrame([[setup.loc[run,"Sample"],setup.loc[run,"P_bar"],setup.loc[run,"T_C"],setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"STppm"],melt_wf["Fe3FeT"],fO2_,gp.log(C_CO32),gp.log(C_H2OT),gp.log(C_S2),gp.log(C_S6),gp.log(C_H2S),gp.log(C_H2),gp.log(C_CO),gp.log(C_CH4),gp.log(C_X),M_m]])
                             
        results = results.append(results2, ignore_index=True)
        results.to_csv('capacities.csv', index=False, header=False)
        print(n, setup.loc[run,"Sample"],gp.log(C_CO32),gp.log(C_H2OT),gp.log(C_S2),gp.log(C_S6),gp.log(C_H2S),gp.log(C_H2),gp.log(C_CO),gp.log(C_CH4),M_m)

        
#############################        
### fugacity coefficients ###
#############################

def fugacity_coefficients(first_row,last_row,setup,species,models):
    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulphide solubility","sulphate solubility","sulphide saturation","ideal gas","carbonylsulphide","mass_volume","insolubles","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulphide","option"],models.loc["sulphate","option"],models.loc["sulphur_saturation","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulphide","option"],models.loc["mass_volume","option"],models.loc['insolubles','option'],date.today()]])
    results = results.append(results1, ignore_index=True)
    results1 = ([["Sample","Pressure (bar)","T ('C)","yH2O","yCO2","yH2","yCO","yO2","yCH4","yS2","ySO2","yH2S","yOCS"]])
    results = results.append(results1, ignore_index=True)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        PT={"T":setup.loc[run,"T_C"]}
        PT["P"] = setup.loc[run,"P_bar"]
        
        ### store results ###
        results2 = pd.DataFrame([[setup.loc[run,"Sample"],PT["P"],setup.loc[run,"T_C"],mdv.y_H2O(PT,species,models),mdv.y_CO2(PT,species,models),mdv.y_H2(PT,species,models),mdv.y_CO(PT,species,models),mdv.y_O2(PT,species,models),mdv.y_CH4(PT,species,models),mdv.y_S2(PT,species,models),mdv.y_SO2(PT,species,models),mdv.y_H2S(PT,species,models),mdv.y_OCS(PT,species,models)]])
                             
        results = results.append(results2, ignore_index=True)
        results.to_csv('fugacity_coefficients_outputs.csv', index=False, header=False)
        print(n, setup.loc[run,"Sample"],mdv.y_H2O(PT,species,models))

######################################                
### fO2 range from sulphur content ###
######################################        
        
def fO2_range_from_S_output(first_row,last_row,setup,species,models,p_tol,nr_step,nr_tol):

    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulphide solubility","sulphate solubility","SCSS","SCAS","ideal gas","carbonylsulphide","insolubles","H2S","species X","species X solubility","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulphide","option"],models.loc["sulphate","option"],models.loc["SCSS","option"],models.loc["SCAS","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulphide","option"],models.loc['insolubles','option'],models.loc["H2S_m","option"],models.loc['species X','option'],models.loc['species X solubility','option'],date.today()]])
    results = results.append(results1, ignore_index=True) 
    results1 = pd.DataFrame([["Sample","T ('C)","wm_STppm","wm_H2OT","wm_CO2ppm","P (bar) sulf","S2- SCSS","sulphide saturated?","DFMQ-sulphide","fO2-sulphide","Fe3FeT-sulphide","S6ST-sulphide","P (bar) anh","S6+ SCAS","sulphate saturated?","DFMQ-sulphate","fO2-sulphate","Fe3FeT-sulphate","S6ST-sulphate"]])
    results = results.append(results1, ignore_index=True)

    # run over rows in file
    for n in range(first_row,last_row,1): # number of rows in the table
   
        # Select run conditions
        run = n # row number
        PT = {"T":setup.loc[run,"T_C"]}
        melt_wf = {"H2OT":(setup.loc[run, "H2O"]/100.)}
        melt_wf["CO2"] = setup.loc[run, "CO2ppm"]/1000000.
        melt_wf["HT"] = ((setup.loc[run, "H2O"]/100.)/species.loc["H2O","M"])*species.loc["H2","M"]
        if setup.loc[run,"P_bar"] > 0.:
            PT["P"] = setup.loc[run,"P_bar"]
            P_sat_sulf, P_sat_anh = setup.loc[run,"P_bar"],setup.loc[run,"P_bar"]
            melt_wf["Fe3FeT"] = setup.loc[run, "Fe3FeT"]
            SCAS_,SCSS_,sulphide_sat,DFMQ_1,fO2_1,Fe3FeT_1,S6ST_1,sulphate_sat,DFMQ_2,fO2_2,Fe3FeT_2,S6ST_2 = c.fO2_range_from_S(run,PT,melt_wf,setup,species,models)
        else:
            melt_wf["ST"] = setup.loc[run, "STppm"]/1000000.
            melt_wf["XT"] = setup.loc[run, "Xppm"]/1000000.
            melt_wf["CT"] = ((setup.loc[run, "CO2ppm"]/1000000.)/species.loc["CO2","M"])*species.loc["C","M"]
            P_sat_sulf,P_sat_anh,SCAS_,SCSS_,sulphide_sat,DFMQ_1,fO2_1,Fe3FeT_1,S6ST_1,sulphate_sat,DFMQ_2,fO2_2,Fe3FeT_2,S6ST_2 = c.P_sat_sulf_anh(run,PT,melt_wf,setup,species,models,p_tol,nr_step,nr_tol)

        # store results
        results1 = pd.DataFrame([[setup.loc[run,"Sample"],PT["T"],setup.loc[run,"STppm"],setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],P_sat_sulf,SCSS_,sulphide_sat,DFMQ_1,fO2_1,Fe3FeT_1,S6ST_1,P_sat_anh,SCAS_,sulphate_sat,DFMQ_2,fO2_2,Fe3FeT_2,S6ST_2]])
        results = results.append(results1, ignore_index=True)
        results.to_csv('fO2_range_from_S.csv', index=False, header=False)
        print(setup.loc[run,"Sample"])


        
########################################
### measured parameters within error ### 
########################################
        
        
def compositions_within_error_output(setup,run,iterations):
    
    # set up results table
    results = pd.DataFrame([["Sample","T_C",
                  "SiO2","TiO2","Al2O3","FeOT","MnO","MgO","CaO","Na2O","K2O","P2O5",
                "H2OT","CO2ppm","STppm","Fe3FeT"]])
    results1 = pd.DataFrame([[setup.loc[run,"Sample"],setup.loc[run,"T_C"],setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],setup.loc[run,"FeOT"],setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],
                setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"STppm"], setup.loc[run,"Fe3FeT"]]])
                             
    results = results.append(results1, ignore_index=True)
    results1 = pd.DataFrame([["sds",setup.loc[run,"TC_sd"],setup.loc[run,"SiO2_sd"],setup.loc[run,"TiO2_sd"],setup.loc[run,"Al2O3_sd"],setup.loc[run,"FeOT_sd"],setup.loc[run,"MnO_sd"],setup.loc[run,"MgO_sd"],setup.loc[run,"CaO_sd"],setup.loc[run,"Na2O_sd"],setup.loc[run,"K2O_sd"],setup.loc[run,"P2O5_sd"],
                setup.loc[run,"H2O_sd"],setup.loc[run,"CO2ppm_sd"],setup.loc[run,"STppm_sd"], setup.loc[run,"Fe3FeT_sd"]]])
                             
    results = results.append(results1, ignore_index=True)
    results1 = pd.DataFrame([["sd types",setup.loc[run,"TC_sd_type"],setup.loc[run,"SiO2_sd_type"],setup.loc[run,"TiO2_sd_type"],setup.loc[run,"Al2O3_sd_type"],setup.loc[run,"FeOT_sd_type"],setup.loc[run,"MnO_sd_type"],setup.loc[run,"MgO_sd_type"],setup.loc[run,"CaO_sd_type"],setup.loc[run,"Na2O_sd_type"],setup.loc[run,"K2O_sd_type"],setup.loc[run,"P2O5_sd_type"],
                setup.loc[run,"H2O_sd_type"],setup.loc[run,"CO2ppm_sd_type"],setup.loc[run,"STppm_sd_type"], setup.loc[run,"Fe3FeT_sd_type"]]])
                             
    results = results.append(results1, ignore_index=True)
    for n in range(0,iterations,1): # n is number of rows of data in conditions file
        TC,SiO2,TiO2,Al2O3,FeOT,MnO,MgO,CaO,Na2O,K2O,P2O5,H2O,CO2ppm,STppm,Fe3FeT = c.compositions_within_error(run,setup)
        results1 = pd.DataFrame([[run,TC,SiO2,TiO2,Al2O3,FeOT,MnO,MgO,CaO,Na2O,K2O,P2O5,H2O,CO2ppm,STppm,Fe3FeT]])
        results = results.append(results1, ignore_index=True)
        results.to_csv('random_compositions.csv', index=False, header=False)
        print(n, setup.loc[run,"Sample"],SiO2)
        
