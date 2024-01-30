# batch_calculations.py

import pandas as pd
from datetime import date
import gmpy2 as gp
import numpy as np
import datetime
import math as math

import VolFe.melt_gas as mg
import VolFe.equilibrium_equations as eq
import VolFe.isotopes as iso

import VolFe.model_dependent_variables as mdv
import VolFe.calculations as c



# building results tables
# outputing sample name
def results_table_sample_name(setup,run):
    results_headers = pd.DataFrame([["sample"]])
    results_values = pd.DataFrame([[setup.loc[run,"Sample"]]])
    return results_headers, results_values
# outputting melt composition, T, P
def results_table_melt_comp_etc(PT,melt_comp,conc,frac,melt_wf):
    results_headers = pd.DataFrame([["T_C","P_bar",
        "SiO2_wtpc", "TiO2_wtpc", "Al2O3_wtpc", "FeOT_wtpc", "MnO_wtpc", "MgO_wtpc", "CaO_wtpc", "Na2O_wtpc", "K2O_wtpc", "P2O5_wtpc",
        "H2OT_wtpc","OH_wtpc","H2Omol_wtpc","H2_ppmw","CH4_ppmw","CO2T_ppmw","CO2mol_ppmw","CO32-_ppmw","CO_ppmw","S2-_ppmw","S6+_ppmw","H2S_ppmw",
        "H_H2OT/HT", "H_H2/HT", "H_CH4/HT", "H_H2S/HT", "C_CO2T/CT", "C_CO/CT", "C_CH4/CT", "S2-/ST", "S6+/ST", "H2S/ST", "Fe3+/FeT"]])
    results_values = pd.DataFrame([[PT["T"],PT["P"],
                melt_comp["SiO2"]*100., melt_comp["TiO2"]*100., melt_comp["Al2O3"]*100., melt_comp["FeOT"]*100., melt_comp["MnO"]*100., melt_comp["MgO"]*100., melt_comp["CaO"]*100., melt_comp["Na2O"]*100., melt_comp["K2O"]*100., melt_comp["P2O5"]*100.,
                conc["wm_H2O"]*100.,conc["wm_OH"]*100,conc["wm_H2Omol"]*100.,conc["wm_H2"]*1000000.,conc["wm_CH4"]*1000000.,conc["wm_CO2"]*1000000.,conc["wm_CO2mol"]*1000000,conc["wm_CO2carb"]*1000000,conc["wm_CO"]*1000000.,conc["wm_S2m"]*1000000.,conc["wm_S6p"]*1000000.,conc["wm_H2S"]*1000000.,
                frac["H2O_HT"], frac["H2_HT"], frac["CH4_HT"], frac["H2S_HT"], frac["CO2_CT"], frac["CO_CT"], frac["CH4_CT"], frac["S2m_ST"], frac["S6p_ST"], frac["H2S_ST"],melt_wf["Fe3FeT"]]])
    return results_headers, results_values
def results_table_melt_vol():
    results_headers = pd.DataFrame([["H2OT-eq_wtpc","CO2T-eq_ppmw","ST_ppmw","X_ppmw"]])
    return results_headers
# outputting model options used in the calculation
def results_table_model_options(models): 
    results_headers = pd.DataFrame([["setup opt","insoluble opt","H2S_m opt","species X opt","Hspeciation opt",
                 "fO2 opt","NNObuffer opt","FMQbuffer opt",
                 "carbon dioxide opt","water opt","hydrogen opt","sulfide opt","sulfate opt","hydrogen sulfide opt","methane opt","carbon monoxide opt","species X solubility opt","Cspeccomp opt","Hspeccomp opt",
                 "SCSS opt","SCAS opt","sulfur_saturation opt","sulfur_is_sat opt","graphite_saturation opt","ideal_gas opt",
                 "y_CO2 opt","y_SO2 opt","y_H2S opt","y_H2 opt","y_O2 opt","y_S2 opt","y_CO opt","y_CH4 opt","y_H2O opt","y_OCS opt","y_X opt",
                 "KHOg opt","KHOSg opt","KOSg opt","KOSg2 opt","KCOg opt","KCOHg opt","KOCSg opt","KCOs opt","carbonylsulfide opt",
                 "density opt","Date"]])
    results_values = pd.DataFrame([[models.loc["setup","option"],models.loc["insolubles","option"], models.loc["H2S_m","option"], models.loc["species X","option"],models.loc["Hspeciation","option"], 
                models.loc["fO2","option"], models.loc["NNObuffer","option"], models.loc["FMQbuffer","option"],
                 models.loc["carbon dioxide","option"], models.loc["water","option"], models.loc["hydrogen","option"], models.loc["sulfide","option"], models.loc["sulfate","option"], models.loc["hydrogen sulfide","option"], models.loc["methane","option"], models.loc["carbon monoxide","option"], models.loc["species X solubility","option"], models.loc["Cspeccomp","option"], models.loc["Hspeccomp","option"],
                 models.loc["SCSS","option"], models.loc["SCAS","option"], models.loc["sulfur_saturation","option"], models.loc["sulfur_is_sat","option"], models.loc["graphite_saturation","option"], models.loc["ideal_gas","option"],
                 models.loc["y_CO2","option"], models.loc["y_SO2","option"], models.loc["y_H2S","option"], models.loc["y_H2","option"], models.loc["y_O2","option"], models.loc["y_S2","option"], models.loc["y_CO","option"], models.loc["y_CH4","option"], models.loc["y_H2O","option"],models.loc["y_OCS","option"], models.loc["y_X","option"],
                 models.loc["KHOg","option"], models.loc["KHOSg","option"], models.loc["KOSg","option"], models.loc["KOSg2","option"], models.loc["KCOg","option"], models.loc["KCOHg","option"],models.loc["KOCSg","option"], models.loc["KCOs","option"],models.loc["carbonylsulfide","option"],
                 models.loc["density","option"],datetime.datetime.now()]])
    return results_headers, results_values
# outputting fugacities, partial pressures, gas mole fraction, fugacity coefficients, molecular masses, solubility constants, equilibrium constants, melt density
def results_table_f_p_xg_y_M_C_K_d(PT,melt_wf,species,models): 
    results_headers = pd.DataFrame([["fO2_DNNO","fO2_DFMQ",
                "fO2_bar","fH2_bar","fH2O_bar","fS2_bar","fSO2_bar","fH2S_bar","fCO2_bar","fCO_bar","fCH4_bar","fOCS_bar","fX_bar",
                "pO2_bar","pH2_bar","pH2O_bar","pS2_bar","pSO2_bar","pH2S_bar","pCO2_bar","pCO_bar","pCH4_bar","pOCS_bar","pX_bar",
                "xgO2_mf","xgH2_mf","xgH2O_mf","xgS2_mf","xgSO2_mf","xgH2S_mf","xgCO2_mf","xgCO_mf","xgCH4_mf","xgOCS_mf","xgX_mf","xgC_S_mf",
                "yO2","yH2","yH2O","yS2","ySO2","yH2S","yCO2","yCO","yCH4","yOCS","yX",
                "M_m_SO","M_m_ox",
                "C_H2O_mf_bar","C_H2_ppm_bar","C_CO2T_mf_bar","C_CO_ppm_bar","C_CH4_ppm_bar","C_S_ppm","C_SO4_ppm_bar","C_H2S_ppm_bar","C_X_ppm_bar",
                "KHOg","KHOSg","KCOg","KCOHg","KOCSg","KSOg","KSOg2","KHOm","KCOm","KCOs",
                "density_gcm3"]])
    results_values = pd.DataFrame([[mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"NNO",models),mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"FMQ",models),
                mdv.f_O2(PT,melt_wf,species,models),mg.f_H2(PT,melt_wf,species,models),mg.f_H2O(PT,melt_wf,species,models),mg.f_S2(PT,melt_wf,species,models),mg.f_SO2(PT,melt_wf,species,models),mg.f_H2S(PT,melt_wf,species,models),mg.f_CO2(PT,melt_wf,species,models),mg.f_CO(PT,melt_wf,species,models),mg.f_CH4(PT,melt_wf,species,models),mg.f_OCS(PT,melt_wf,species,models),mg.f_X(PT,melt_wf,species,models),
                mg.p_O2(PT,melt_wf,species,models),mg.p_H2(PT,melt_wf,species,models),mg.p_H2O(PT,melt_wf,species,models),mg.p_S2(PT,melt_wf,species,models),mg.p_SO2(PT,melt_wf,species,models),mg.p_H2S(PT,melt_wf,species,models),mg.p_CO2(PT,melt_wf,species,models),mg.p_CO(PT,melt_wf,species,models),mg.p_CH4(PT,melt_wf,species,models),mg.p_OCS(PT,melt_wf,species,models),mg.p_X(PT,melt_wf,species,models),
                mg.xg_O2(PT,melt_wf,species,models),mg.xg_H2(PT,melt_wf,species,models),mg.xg_H2O(PT,melt_wf,species,models),mg.xg_S2(PT,melt_wf,species,models),mg.xg_SO2(PT,melt_wf,species,models),mg.xg_H2S(PT,melt_wf,species,models),mg.xg_CO2(PT,melt_wf,species,models),mg.xg_CO(PT,melt_wf,species,models),mg.xg_CH4(PT,melt_wf,species,models),mg.xg_OCS(PT,melt_wf,species,models),mg.xg_X(PT,melt_wf,species,models),mg.gas_CS(PT,melt_wf,species,models),
                mdv.y_O2(PT,species,models),mdv.y_H2(PT,species,models),mdv.y_H2O(PT,species,models),mdv.y_S2(PT,species,models),mdv.y_SO2(PT,species,models),mdv.y_H2S(PT,species,models),mdv.y_CO2(PT,species,models),mdv.y_CO(PT,species,models),mdv.y_CH4(PT,species,models),mdv.y_OCS(PT,species,models),mdv.y_X(PT,species,models),
                mg.M_m_SO(melt_wf,species),mg.M_m_ox(melt_wf,species,models),
                mdv.C_H2O(PT,melt_wf,species,models),mdv.C_H2(PT,melt_wf,species,models),mdv.C_CO3(PT,melt_wf,species,models),mdv.C_CO(PT,melt_wf,species,models),mdv.C_CH4(PT,melt_wf,species,models),mdv.C_S(PT,melt_wf,species,models),mdv.C_SO4(PT,melt_wf,species,models),mdv.C_H2S(PT,melt_wf,species,models),mdv.C_X(PT,melt_wf,species,models),
                mdv.KHOg(PT,models),mdv.KHOSg(PT,models),mdv.KCOg(PT,models),mdv.KCOHg(PT,models),mdv.KOCSg(PT,models),mdv.KOSg(PT,models),mdv.KOSg2(PT,models),mdv.KHOm(PT,melt_wf,species,models),mdv.KCOm(PT,melt_wf,species,models),mdv.KCOs(PT,models),
                mdv.melt_density(PT,melt_wf,species,models)]])
    return results_headers, results_values
# headers for open system degassing all gas
def results_table_open_all_gas():
    results_headers = pd.DataFrame([["xgO2_all_mf","xgH2_all_mf","xgH2O_all_mf","xgS2_all_mf","xgSO2_all_mf","xgH2S_all_mf","xgCO2_all_mf","xgCO_all_mf","xgCH4_all_mf","xgOCS_all_mf","xgX_all_mf","xgC_S_all_mf"]])
    return results_headers
# saturation conditions
def results_table_sat(sulf_sat_result,PT,melt_wf,species,models):
    results_headers = pd.DataFrame([["SCSS_ppm","sulfide saturated","SCAS_ppm","anhydrite saturated","ST melt if sat","graphite saturated"]])
    results_values = pd.DataFrame([[sulf_sat_result["SCSS"],sulf_sat_result["sulfide_sat"],sulf_sat_result["SCAS"],sulf_sat_result["sulfate_sat"],sulf_sat_result["ST"],c.graphite_saturation(PT,melt_wf,species,models)]])
    return results_headers, results_values

def options_from_setup(run,models,setup):
    """ 
    Allows model options to be read from the setup file rather than models file.


    Parameters
    ----------
    run: float
        Integer of the row in the setup file to read from (note the first row under the headers is row 0).   
    setup: pandas.DataFrame
        Dataframe with melt compositions to be used, require header using the same labels as row labels from models file if you want to use that option.
    models: pandas.DataFrame
        Dataframe of models.csv file.

    Returns
    -------
    results: pandas.DataFrame

    """
    if models.loc["setup","option"] == "no":
        return models
    elif models.loc["setup","option"] == "yes":
        models = pd.read_csv("models.csv", index_col = [0])
    # species
    if models.loc["insolubles","option"] == "setup":
            models.loc["insolubles","option"] = setup.loc[run,"insolubles"]
    if models.loc["H2S_m","option"] == "setup":
            models.loc["H2S_m","option"] = setup.loc[run,"H2S_m"]
    if models.loc["species X","option"] == "setup":
            models.loc["species X","option"] = setup.loc[run,"species X"]
    if models.loc["Hspeciation","option"] == "setup":
            models.loc["Hspeciation","option"] = setup.loc[run,"Hspeciation"]
    # oxygen fugacity
    if models.loc["fO2","option"] == "setup":
            models.loc["fO2","option"] = setup.loc[run,"fO2"]
    if models.loc["NNObuffer","option"] == "setup":
            models.loc["NNObuffer","option"] = setup.loc[run,"NNObuffer"]
    if models.loc["FMQbuffer","option"] == "setup":
            models.loc["FMQbuffer","option"] = setup.loc[run,"FMQbuffer"]
    # solubility constants
    if models.loc["carbon dioxide","option"] == "setup":
            models.loc["carbon dioxide","option"] = setup.loc[run,"carbon dioxide"]
    if models.loc["water","option"] == "setup":
            models.loc["water","option"] = setup.loc[run,"water"]
    if models.loc["hydrogen","option"] == "setup":
            models.loc["hydrogen","option"] = setup.loc[run,"hydrogen"]
    if models.loc["sulfide","option"] == "setup":
            models.loc["sulfide","option"] = setup.loc[run,"sulfide"]
    if models.loc["sulfate","option"] == "setup":
            models.loc["sulfate","option"] = setup.loc[run,"sulfate"]
    if models.loc["hydrogen sulfide","option"] == "setup":
            models.loc["hydrogen sulfide","option"] = setup.loc[run,"hydrogen sulfide"]
    if models.loc["methane","option"] == "setup":
            models.loc["methane","option"] = setup.loc[run,"methane"]
    if models.loc["carbon monoxide","option"] == "setup":
            models.loc["carbon monoxide","option"] = setup.loc[run,"carbon monoxide"]
    if models.loc["species X solubility","option"] == "setup":
            models.loc["species X solubility","option"] = setup.loc[run,"species X solubility"]
    if models.loc["Cspeccomp","option"] == "setup":
            models.loc["Cspeccomp","option"] = setup.loc[run,"Cspeccomp"]
    if models.loc["Hspeccomp","option"] == "setup":
            models.loc["Hspeccomp","option"] = setup.loc[run,"Hspeccomp"]
    # saturation conditions
    if models.loc["SCSS","option"] == "setup":
            models.loc["SCSS","option"] = setup.loc[run,"SCSS"]
    if models.loc["SCAS","option"] == "setup":
            models.loc["SCAS","option"] = setup.loc[run,"SCAS"]
    if models.loc["sulfur_saturation","option"] == "setup":
            models.loc["sulfur_saturation","option"] = setup.loc[run,"sulfur_saturation"]
    if models.loc["sulfur_is_sat","option"] == "setup":
            models.loc["sulfur_is_sat","option"] = setup.loc[run,"sulfur_is_sat"]
    if models.loc["graphite_saturation","option"] == "setup":
            models.loc["graphite_saturation","option"] = setup.loc[run,"graphite_saturation"]
    # fugacity coefficients
    if models.loc["ideal_gas","option"] == "setup":
            models.loc["ideal_gas","option"] = setup.loc[run,"ideal_gas"]
    if models.loc["y_CO2","option"] == "setup":
            models.loc["y_CO2","option"] = setup.loc[run,"y_CO2","option"]
    if models.loc["y_SO2","option"] == "setup":
            models.loc["y_SO2","option"] = setup.loc[run,"y_SO2","option"]
    if models.loc["y_H2S","option"] == "setup":
            models.loc["y_H2S","option"] = setup.loc[run,"y_H2S","option"]
    if models.loc["y_H2","option"] == "setup":
            models.loc["y_H2","option"] = setup.loc[run,"y_H2","option"]
    if models.loc["y_O2","option"] == "setup":
            models.loc["y_O2","option"] = setup.loc[run,"y_O2","option"]
    if models.loc["y_S2","option"] == "setup":
            models.loc["y_S2","option"] = setup.loc[run,"y_S2","option"]
    if models.loc["y_CO","option"] == "setup":
            models.loc["y_CO","option"] = setup.loc[run,"y_CO","option"]
    if models.loc["y_CH4","option"] == "setup":
            models.loc["y_CH4","option"] = setup.loc[run,"y_CH4","option"]
    if models.loc["y_H2O","option"] == "setup":
            models.loc["y_H2O","option"] = setup.loc[run,"y_H2O","option"]
    if models.loc["y_OCS","option"] == "setup":
            models.loc["y_OCS","option"] = setup.loc[run,"y_OCS","option"]
    if models.loc["y_X","option"] == "setup":
            models.loc["y_X","option"] = setup.loc[run,"y_X","option"]
    # equilibrium constants
    if models.loc["KHOg","option"] == "setup":
            models.loc["KHOg","option"] = setup.loc[run,"KHOg","option"]
    if models.loc["KHOSg","option"] == "setup":
            models.loc["KHOSg","option"] = setup.loc[run,"KHOSg","option"]
    if models.loc["KOSg","option"] == "setup":
            models.loc["KOSg","option"], = setup.loc[run,"KOSg","option"]
    if models.loc["KOSg2","option"] == "setup":
            models.loc["KOSg2","option"] = setup.loc[run,"KOSg2","option"]
    if models.loc["KCOg","option"] == "setup":
            models.loc["KCOg","option"] = setup.loc[run,"KCOg","option"]
    if models.loc["KCOHg","option"] == "setup":
            models.loc["KCOHg","option"] = setup.loc[run,"KCOHg","option"]
    if models.loc["KOCSg","option"] == "setup":
            models.loc["KOCSg","option"] = setup.loc[run,"KOCSg","option"]
    if models.loc["KCOs","option"] == "setup":
            models.loc["KCOs","option"] = setup.loc[run,"KCOs","option"]
    if models.loc["carbonylsulfide","option"] == "setup":
            models.loc["carbonylsulfide","option"] = setup.loc[run,"carbonylsulfide","option"]
    # other
    if models.loc["density","option"] == "setup":
            models.loc["density","option"] = setup.loc[run,"density","option"]
    return models


# calculate the saturation pressure for multiple melt compositions in setup file
def P_sat_output(first_row,last_row,p_tol,nr_step,nr_tol,setup,species,models):
    
    """ 
    Calculates the pressure of vapor saturation for multiple melt compositions given volatile-free melt composition, volatile content, temperature, and an fO2 estimate.


    Parameters
    ----------
    first_row: float
        Integer of the first row in the setup file to run (note the first row under the headers is row 0).   
    last_row: float
        Integer of the last row in the setup file to run (note the first row under the headers is row 0).
    p_tol: float
        Required tolerance for convergence of Pvsat in bars (0.1 bars is normally used).
    nr_step: float
        Step size for Newton-Raphson solver for melt speciation (typically 1 is fine, but this can be made smaller if there are problems with convergence.).
    nr_tol: float
        Tolerance for the Newton-Raphson solver for melt speciation in weight fraction (typically 1e-9 is sufficient but this can be made larger if there are problems with convergence).
    setup: pandas.DataFrame
        Dataframe with melt compositions to be used, requires following headers): 
        Sample, T_C, 
        DNNO or DFMQ or logfO2 or (Fe2O3 and FeO) or Fe3FeT or S6ST
        SiO2, TiO2, Al2O3, (Fe2O3T or FeOT unless Fe2O3 and FeO given), MnO, MgO, CaO, Na2O, K2O, P2O5, 
        H2O and/or CO2ppm and/or STppm and/or Xppm
        Note: concentrations (unless otherwise stated) are in wt%
    species: pandas.DataFrame
        Dataframe of species.csv file.
    models: pandas.DataFrame
        Dataframe of models.csv file.

    Returns
    -------
    results: pandas.DataFrame

    Outputs
    -------
    results_saturation_pressures: csv file (if output csv = yes in models)

    """
    
    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        PT={"T":setup.loc[run,"T_C"]}
        melt_wf=mg.melt_comp(run,setup)
        melt_wf['CO2'] = setup.loc[run,"CO2ppm"]/1000000.
        melt_wf["H2OT"] = setup.loc[run,"H2O"]/100.

        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)
        
        # calculate Pvsat assuming only H2O CO2 in vapour and melt
        if setup.loc[run,"Fe3FeT"] > 0.:
            melt_wf['Fe3FeT'] = setup.loc[run,"Fe3FeT"]
        else:
            melt_wf['Fe3FeT'] = 0.
        P_sat_H2O_CO2_only, P_sat_H2O_CO2_result = c.P_sat_H2O_CO2(PT,melt_wf,species,models,p_tol,nr_step,nr_tol)
        
        if models.loc["calc_sat","option"] == "fO2_fX":
            P_sat_fO2_fS2_result = c.P_sat_fO2_fS2(PT,melt_wf,species,models,p_tol)
            PT["P"] = P_sat_fO2_fS2_result["P_tot"]
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
            raise TypeError('This is not currently possible')
        if models.loc["sulfur_is_sat","option"] == "yes":
            if melt_wf["XT"] > 0.:
                raise TypeError('This is not currently possible')
            P_sat_, conc, frac  = c.fO2_P_VSA(PT,melt_wf,species,models,nr_step,nr_tol,p_tol)
        elif models.loc["sulfur_saturation","option"] == "no":
            P_sat_, conc, frac = c.P_sat(PT,melt_wf,species,models,p_tol,nr_step,nr_tol)
        elif models.loc["sulfur_saturation","option"] == "yes":
            if melt_wf["XT"] > 0.:
                raise TypeError('This is not currently possible')
            P_sat_, conc, frac = c.P_VSA(PT,melt_wf,species,models,nr_step,nr_tol,p_tol)
        PT["P"] = P_sat_
        melt_wf["H2OT"] = conc["wm_H2O"]
        melt_wf["CO2"] = conc["wm_CO2"]
        melt_wf["S2-"] = conc["wm_S2m"]
        if models.loc["sulfur_is_sat","option"] == "yes":
            melt_wf["Fe3FeT"] = frac["Fe3FeT"]
        else:
            melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
        
        sulf_sat_result = c.sulfur_saturation(PT,melt_wf,species,models)
        # gas_mf = {"O2":mg.xg_O2(PT,melt_wf,species,models),"CO":mg.xg_CO(PT,melt_wf,species,models),"CO2":mg.xg_CO2(PT,melt_wf,species,models),"H2":mg.xg_H2(PT,melt_wf,species,models),"H2O":mg.xg_H2O(PT,melt_wf,species,models),"CH4":mg.xg_CH4(PT,melt_wf,species,models),"S2":mg.xg_S2(PT,melt_wf,species,models),"SO2":mg.xg_SO2(PT,melt_wf,species,models),"H2S":mg.xg_H2S(PT,melt_wf,species,models),"OCS":mg.xg_OCS(PT,melt_wf,species,models),"X":mg.xg_X(PT,melt_wf,species,models),"Xg_t":mg.Xg_tot(PT,melt_wf,species,models),"wt_g":0.}     
        melt_comp = mg.melt_normalise_wf(melt_wf,species,"yes","no")  
        
        # create results
        results_headers_table_sample_name, results_values_table_sample_name = results_table_sample_name(setup,run)
        results_headers_table_melt_comp_etc, results_values_table_melt_comp_etc = results_table_melt_comp_etc(PT,melt_comp,conc,frac,melt_wf)
        results_headers_table_model_options, results_values_table_model_options = results_table_model_options(models)    
        results_headers_table_f_p_xg_y_M_C_K_d, results_values_table_f_p_xg_y_M_C_K_d = results_table_f_p_xg_y_M_C_K_d(PT,melt_wf,species,models)
        results_headers_table_sat, results_values_table_sat = results_table_sat(sulf_sat_result,PT,melt_wf,species,models)
        results_headers_table_melt_vol = results_table_melt_vol() # "H2OT-eq_wtpc","CO2T-eq_ppmw","ST_ppmw","X_ppmw"
        results_values_table_melt_vol = pd.DataFrame([[setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"STppm"],setup.loc[run,"Xppm"]]])
        results_headers_table_H2OCO2only = pd.DataFrame([["Pvsat (H2O CO2 only)", "xg_H2O (H2O CO2 only)", "xg_CO2 (H2O CO2 only)","f_H2O (H2O CO2 only)", "f_CO2 (H2O CO2 only)","p_H2O (H2O CO2 only)", "p_CO2 (H2O CO2 only)", "Pvsat_diff_bar"]])
        results_values_table_H2OCO2only = pd.DataFrame([[P_sat_H2O_CO2_only, P_sat_H2O_CO2_result["xg_H2O"], P_sat_H2O_CO2_result["xg_CO2"], P_sat_H2O_CO2_result["f_H2O"], P_sat_H2O_CO2_result["f_CO2"], P_sat_H2O_CO2_result["p_H2O"], P_sat_H2O_CO2_result["p_CO2"], (P_sat_H2O_CO2_only-PT["P"])]])
        results_headers = pd.concat([results_headers_table_sample_name,results_headers_table_melt_comp_etc,results_headers_table_melt_vol,results_headers_table_sat,results_headers_table_H2OCO2only,results_headers_table_f_p_xg_y_M_C_K_d,results_headers_table_model_options],axis=1)
        results1 = pd.concat([results_values_table_sample_name,results_values_table_melt_comp_etc,results_values_table_melt_vol,results_values_table_sat,results_values_table_H2OCO2only,results_values_table_f_p_xg_y_M_C_K_d,results_values_table_model_options],axis=1)
    
        if n == first_row:
            results = pd.concat([results_headers, results1])
        else:                         
            results = pd.concat([results, results1])
        
        if models.loc["print status","option"] == "yes":
            print(n, setup.loc[run,"Sample"],PT["P"])
    
    results.columns = results.iloc[0]
    results = results[1:]  
    if models.loc["output csv","option"] == "yes":
        results.to_csv('results_saturation_pressures.csv', index=False, header=True)
    
    return results

### NEEDS CHECKING ###        
def P_sat_output_fS2(first_row,last_row,p_tol,nr_step,nr_tol,setup,species,models):
    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","sulfide saturation","ideal gas","carbonylsulfide","mass_volume","insolubles","Saturation calculation","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["sulfur_saturation","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["mass_volume","option"],models.loc['insolubles','option'],models.loc['calc_sat','option'],date.today()]])
    results = pd.concat([results, results1], ignore_index=True)
    results1 = ([["Sample","Pressure (bar)","Saturation pressure (bars)","T ('C)","fO2 (DNNO)","fO2 (DFMQ)",
                  "SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)","ST (ppm)","S6/ST","Fe3/FeT",
                  "SCSS (ppm)","sulfide saturated","SCAS (ppm)","anhydrite saturated","S melt (ppm)",
                "f-fO2","f-fS2","f-fSO2","f-pO2","f-pS2","f-pSO2","f-xgO2","f-xgS2","f-xgSO2",
                 "b-fO2","b-fS2","b-fSO2","b-pO2","b-pS2","b-pSO2","b-xgO2","b-xgS2","b-xgSO2","ySO2"]])
    results = pd.concat([results, results1], ignore_index=True)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n

        PT={"T":setup.loc[run,"T_C"]}
        melt_wf=mg.melt_comp(run,setup)
        melt_wf['CO2'] = 0.
        melt_wf["H2OT"] = 0.

        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        P_sat_, wm_ST, fSO2, wm_S2m, wm_S6p, pS2, pO2, pSO2, xgS2, xgO2, xgSO2 = c.P_sat_fO2_fS2(PT,melt_wf,species,models,p_tol)
        if setup.loc[run,"P_bar"] > 0.:
            PT["P"] = setup.loc[run,"P_bar"]
        else:
            PT["P"] = P_sat_
        melt_wf["ST"] = wm_ST
        melt_wf["S2-"] = wm_S2m
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
        SCSS_,sulfide_sat,SCAS_,sulfate_sat, ss_ST = c.sulfur_saturation(PT,melt_wf,species,models)
        gas_mf = {"O2":mg.xg_O2(PT,melt_wf,species,models),"S2":mg.xg_S2(PT,melt_wf,species,models),"SO2":mg.xg_SO2(PT,melt_wf,species,models)}
        
        ### store results ###
        results2 = pd.DataFrame([[setup.loc[run,"Sample"],PT["P"],P_sat_,
                setup.loc[run,"T_C"],mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"NNO",models),mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"FMQ",models),setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(melt_wf,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],wm_ST,mg.S6ST(PT,melt_wf,species,models),melt_wf["Fe3FeT"],SCSS_,sulfide_sat,SCAS_,sulfate_sat,ss_ST,
                mdv.f_O2(PT,melt_wf,species,models),setup.loc[run,"fS2"],fSO2, pO2,pS2,pSO2, xgO2,xgS2,xgSO2,mdv.f_O2(PT,melt_wf,species,models),mg.f_S2(PT,melt_wf,species,models),mg.f_SO2(PT,melt_wf,species,models), mg.p_O2(PT,melt_wf,species,models),mg.p_S2(PT,melt_wf,species,models),mg.p_SO2(PT,melt_wf,species,models), mg.xg_O2(PT,melt_wf,species,models),mg.xg_S2(PT,melt_wf,species,models),mg.xg_SO2(PT,melt_wf,species,models),mg.y_SO2(PT,species,models)]])
                             
        results = pd.concat([results, results1], ignore_index=True)
        if models.loc["output csv","option"] == "yes":
            results.to_csv('saturation_pressures_fS2.csv', index=False, header=False)
        if models.loc["print status","option"] == "yes":
            print(n, setup.loc[run,"Sample"],PT["P"])

###############
### gassing ###
###############

def gassing(run,gassing_inputs,setup,species,models):
     
    """ 
    Calculates the pressure of vapor saturation for multiple melt compositions given volatile-free melt composition, volatile content, temperature, and an fO2 estimate.


    Parameters
    ----------
    run: float
        Integer of the row in the setup file to run (note the first row under the headers is row 0).
    gassing_inputs: pandas.DataFrame
        Dataframe containing the following:
        dp_step: Pressure step size for gassing calculation in bars
        nr_step: Step size for Newton-Raphson solver for melt speciation (typically 1 is fine, but this can be made smaller if there are problems with convergence).
        nr_tol: Tolerance for the Newton-Raphson solver for melt speciation in weight fraction (typically 1e-6 is sufficient but this can be made larger if there are problems with convergence).
        psat_tol: Required tolerance for convergence of Pvsat in bars (0.1 bars is normally used).
        dwtg: Amount of gas to add at each step if regassing in an open-system in wt fraction total system
        i_nr_step = Step-size for newton-raphson convergence for isotopes (typically 1e-1 is fine, but this can be made smaller if there are problems with convergence).
        i_nr_tol = Tolerance for newton-raphson convergence for isotopes (typically 1e-9 is sufficient but this can be made larger if there are problems with convergence).
    setup: pandas.DataFrame
        Dataframe with melt composition to be used, requires following headers (notes in [] are not part of the headers): 
        Sample, T_C, 
        DNNO or DFMQ or logfO2 or (Fe2O3 and FeO) or Fe3FeT or S6ST, [at initial pressure]
        SiO2, TiO2, Al2O3, (Fe2O3T or FeOT unless Fe2O3 and FeO given), MnO, MgO, CaO, Na2O, K2O, P2O5, [concentrations are in wt%]
        (H2O and/or CO2ppm and/or STppm and/or Xppm) [concentration of H2O in wt%]
        P_bar [IF starting from a given pressure]
        final_P [IF regassing, pressure calculation stops at in bars]
        wt_g [IF starting from given pressure and gas is present, can specifiy the gas present in wt%]
        initial_CO2wtpc [IF starting from given pressure and gas is present, can specifiy initial composition using initial CO2 dissolved in the melt in wt%]
        xg_O2, xg_CO, xg_H2, xg_S2, xg_X [IF starting from a given pressure, need initial guesses for mole fraction of gas species]
    species: pandas.DataFrame
        Dataframe of species.csv file.
    models: pandas.DataFrame
        Dataframe of models.csv file.

    Returns
    -------
    results: pandas.DataFrame

    Outputs
    -------
    If output csv = yes in models
    results_gassing_chemistry: csv file

    """

    if models.loc["print status","option"] == "yes":
        print(setup.loc[run,"Sample"])

    # check if any options need to be read from the setup file rather than the models file
    models = options_from_setup(run,models,setup)

    if models.loc["fO2","option"] != "Kress91A":
        raise TypeError("Change 'fO2' option in models to 'Kress91A' (other fO2 options are not currently supported)")
    if models.loc["gassing_direction","option"] == "regas" and models.loc["starting_P","option"] == "bulk":
        raise TypeError("gassing_direction = regas and starting_P = bulk in models are not compatible")

    # set T, volatile composition of the melt, and tolerances
    PT={"T":setup.loc[run,"T_C"]}
    melt_wf = mg.melt_comp(run,setup)
    melt_wf['CO2'] = setup.loc[run,"CO2ppm"]/1000000.
    melt_wf["H2OT"] = setup.loc[run,"H2O"]/100.
    melt_wf["ST"] = setup.loc[run,"STppm"]/1000000.
    melt_wf["H2"] = 0.
    melt_wf["XT"] = setup.loc[run,"Xppm"]/1000000.
    melt_wf["CT"] = (melt_wf["CO2"]/species.loc["CO2","M"])*species.loc["C","M"]
    melt_wf["HT"] = (2.*melt_wf["H2OT"]/species.loc["H2O","M"])*species.loc["H","M"]
    melt_wf["ST"] = melt_wf["ST"]
    if setup.loc[run,"S6ST"] > 0.:
        melt_wf["S6ST"] = setup.loc[run,"S6ST"]
    nr_step = gassing_inputs["nr_step"]
    nr_tol = gassing_inputs["nr_tol"]
    dp_step = gassing_inputs["dp_step"]
    psat_tol = gassing_inputs["psat_tol"]
    if models.loc["gassing_style","option"] == "open":
        dwtg = gassing_inputs["dwtg"]

    # Calculate saturation pressure for composition given in setup file
    if models.loc["insolubles","option"] == "H2O-CO2 only":  
        P_sat_, P_sat_H2O_CO2_result = c.P_sat_H2O_CO2(PT,melt_wf,species,models,psat_tol,nr_step,nr_tol)
        conc = {"wm_H2O":P_sat_H2O_CO2_result["wm_H2O"], "wm_CO2":P_sat_H2O_CO2_result["wm_CO2"], "wm_H2":0., "wm_CO":0., "wm_CH4":0., "wm_H2S":0., "wm_S2m":0., "wm_S6p":0., "ST": 0.}
        frac = c.melt_species_ratios(conc,species)
    else:
        P_sat_, conc, frac = c.P_sat(PT,melt_wf,species,models,psat_tol,nr_step,nr_tol)
    PT["P"] = P_sat_
    if models.loc["print status","option"] == "yes":
        print(PT["T"],PT["P"],datetime.datetime.now())

    # update melt composition at saturation pressure, check for sulfur saturation, and calculate some things
    melt_wf["H2OT"] = conc["wm_H2O"]
    melt_wf["CO2"] = conc["wm_CO2"]
    melt_wf["CO"] = conc["wm_CO"]
    melt_wf["CH4"] = conc["wm_CH4"]
    melt_wf["H2"] = conc["wm_H2"]
    melt_wf["S2-"] = conc["wm_S2m"]
    melt_wf["S6+"] = conc["wm_S6p"]
    melt_wf["H2S"] = conc["wm_H2S"]
    melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
    melt_wf["S6ST"] = mg.S6ST(PT,melt_wf,species,models)
    sulf_sat_result = c.sulfur_saturation(PT,melt_wf,species,models)    
    wm_CO2eq, wm_H2Oeq = mg.melt_H2O_CO2_eq(melt_wf,species)
    melt_comp = mg.melt_normalise_wf(melt_wf,species,"yes","no")
    
    # Set bulk composition
    bulk_comp = c.bulk_composition(run,PT,melt_wf,setup,species,models)
    bulk_wf = {"C":bulk_comp["wt_C"],"H":bulk_comp["wt_H"],"O":bulk_comp["wt_O"],"S":bulk_comp["wt_S"],"Fe":bulk_comp["wt_Fe"],"Wt":bulk_comp["Wt"],"X":bulk_comp["wt_X"]}

    # set system and initial guesses
    system = eq.set_system(melt_wf,models)
    guesses = eq.initial_guesses(run,PT,melt_wf,setup,species,models,system)    

    # create results
    results_headers_table_sample_name, results_values_table_sample_name = results_table_sample_name(setup,run)
    results_headers_table_melt_comp_etc, results_values_table_melt_comp_etc = results_table_melt_comp_etc(PT,melt_comp,conc,frac,melt_wf)
    results_headers_table_model_options, results_values_table_model_options = results_table_model_options(models)    
    results_headers_table_f_p_xg_y_M_C_K_d, results_values_table_f_p_xg_y_M_C_K_d = results_table_f_p_xg_y_M_C_K_d(PT,melt_wf,species,models)
    results_headers_table_sat, results_values_table_sat = results_table_sat(sulf_sat_result,PT,melt_wf,species,models)
    results_headers_table_melt_vol = results_table_melt_vol() # "H2OT-eq_wtpc","CO2T-eq_ppmw","ST_ppmw","X_ppmw"
    results_values_table_melt_vol = pd.DataFrame([[wm_H2Oeq*100.,wm_CO2eq*1000000.,conc["ST"]*1000000.,melt_wf["XT"]*1000000.]])
    results_headers_table_wtg_etc = pd.DataFrame([["wt_g_wtpc","wt_g_O_wtf","wt_g_C_wtf","wt_g_H_wtf","wt_g_S_wtf","wt_g_X_wtf","wt_O_wtpc","wt_C_wtpc","wt_H_wtpc","wt_S_wtpc","wt_X_wtpc","Solving species"]])
    results_values_table_wtg_etc = pd.DataFrame([[bulk_comp["wt_g"]*100.,"","","","","",bulk_wf["O"]*100.,bulk_wf["C"]*100.,bulk_wf["H"]*100.,bulk_wf["S"]*100.,bulk_wf["X"]*100.,""]])
    if models.loc["gassing_style","option"] == "open" and models.loc["gassing_direction","option"] == "degas":
        results_headers_table_open_all_gas = results_table_open_all_gas()
        results_values_table_open_all_gas = pd.DataFrame([[mg.xg_O2(PT,melt_wf,species,models),mg.xg_H2(PT,melt_wf,species,models),mg.xg_H2O(PT,melt_wf,species,models),mg.xg_S2(PT,melt_wf,species,models),mg.xg_SO2(PT,melt_wf,species,models),mg.xg_H2S(PT,melt_wf,species,models),mg.xg_CO2(PT,melt_wf,species,models),mg.xg_CO(PT,melt_wf,species,models),mg.xg_CH4(PT,melt_wf,species,models),mg.xg_OCS(PT,melt_wf,species,models),mg.xg_X(PT,melt_wf,species,models),mg.gas_CS(PT,melt_wf,species,models)]])
        results_headers = pd.concat([results_headers_table_sample_name,results_headers_table_melt_comp_etc,results_headers_table_melt_vol,results_headers_table_sat,results_headers_table_f_p_xg_y_M_C_K_d,results_headers_table_wtg_etc,results_headers_table_open_all_gas,results_headers_table_model_options],axis=1)
        results1 = pd.concat([results_values_table_sample_name,results_values_table_melt_comp_etc,results_values_table_melt_vol,results_values_table_sat,results_values_table_f_p_xg_y_M_C_K_d,results_values_table_wtg_etc,results_values_table_open_all_gas,results_values_table_model_options],axis=1)
    else:    
        results_headers = pd.concat([results_headers_table_sample_name,results_headers_table_melt_comp_etc,results_headers_table_melt_vol,results_headers_table_sat,results_headers_table_f_p_xg_y_M_C_K_d,results_headers_table_wtg_etc,results_headers_table_model_options],axis=1)
        results1 = pd.concat([results_values_table_sample_name,results_values_table_melt_comp_etc,results_values_table_melt_vol,results_values_table_sat,results_values_table_f_p_xg_y_M_C_K_d,results_values_table_wtg_etc,results_values_table_model_options],axis=1)
    results = pd.concat([results_headers, results1])
    
    # results for isotope calculations...
    if models.loc["isotopes","option"] == "yes":
        raise TypeError("This is not currently supported")
        a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_ = iso.i2s6_S_alphas(PT)
        results_isotopes1 = pd.DataFrame([["P","T_C","xg_O2","xg_CO","xg_CO2","xg_H2","xg_H2O","xg_CH4","xg_S2","xg_SO2","xg_SO3","xg_H2S","xg_OCS","wt_g",
                        "wm_CO2","wm_H2O","wm_H2","wm_S","wm_SO3","wm_ST","Fe3T","S6T",
                         "DFMQ","DNNO","SCSS","sulfide sat?","SCAS","sulfate sat?",
                         "RS_S2-","RS_SO42-","RS_H2S","RS_SO2","RS_S2","R_OCS","R_m","R_g",
                         "dS_S2-","dS_SO42-","dS_H2S","dS_SO2","dS_S2","dS_OCS","dS_m","dS_g",
                         "a_H2S_S2-","a_SO42-_S2-","a_S2_S2-","a_SO2_S2-","a_OCS_S2-","a_g_m"]])
        results_isotopes = results_header.append(results_isotopes1, ignore_index=True)
        results1 = pd.DataFrame([[PT["P"],PT["T"],mg.xg_O2(PT,melt_wf,species,models),mg.xg_CO(PT,melt_wf,species,models),mg.xg_CO2(PT,melt_wf,species,models),mg.xg_H2(PT,melt_wf,species,models),mg.xg_H2O(PT,melt_wf,species,models),mg.xg_CH4(PT,melt_wf,species,models),mg.xg_S2(PT,melt_wf,species,models),mg.xg_SO2(PT,melt_wf,species,models),mg.xg_H2S(PT,melt_wf,species,models),mg.xg_OCS(PT,melt_wf,species,models),wt_g_,
melt_wf["CO2"],melt_wf["H2OT"],0,(mg.wm_S(PT,melt_wf,species,models)/100),(mg.wm_SO3(PT,melt_wf,species,models)/100),melt_wf["ST"],melt_wf["Fe3FeT"],mg.S6ST(PT,melt_wf,species,models),
mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"FMQ"),mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"NNO"),SCSS_,sulfide_sat,SCAS_,sulfate_sat,
R_S_S2_,R_S_SO4_,"","","","",R_i["S"],"",ratio2delta("VCDT",R_S_S2_),ratio2delta("VCDT",R_S_SO4_),"","","","",ratio2delta("VCDT",R_i["S"]),"",
a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_,""]])
        results_isotopes = pd.concat([results_isotopes, results1], ignore_index=True) 
        if models.loc["output csv","option"] == "yes":
            results_isotopes.to_csv('results_gassing_isotopes.csv', index=False, header=False)
    
    if models.loc["P_variation","option"] == "polybaric":
        # pressure ranges and options
        starting_P = models.loc["starting_P","option"]
        if starting_P == "bulk":
            if models.loc["bulk_composition","option"] != "yes":
                raise warning("not that if starting_P = bulk and bulk_composition != yes, Pvsat calculated is not for bulk system")
        if starting_P == "set":
            initial = int(setup.loc[run,"P_bar"])
        else:
            if models.loc["gassing_direction","option"] == "degas":
                answer = PT["P"]/dp_step
                answer = round(answer)
                initial = round(answer*dp_step)
            elif models.loc["gassing_direction","option"] == "regas":
                answer = PT["P"]/dp_step
                answer = round(answer)
                answer = round(answer*dp_step)
                initial = round(answer+dp_step)
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
        gas_mf = {"O2":mg.xg_O2(PT,melt_wf,species,models),"CO":mg.xg_CO(PT,melt_wf,species,models),"CO2":mg.xg_CO2(run,PT,melt_wf,setup,species,models),"H2":mg.xg_H2(PT,melt_wf,species,models),"H2O":mg.xg_H2O(PT,melt_wf,species,models),"CH4":mg.xg_CH4(PT,melt_wf,species,models),"S2":mg.xg_S2(PT,melt_wf,species,models),"SO2":mg.xg_SO2(PT,melt_wf,species,models),"SO3":mg.xg_SO3(PT,melt_wf,species,models),"H2S":mg.xg_H2S(run,PT,melt_wf,setup,species,models),"OCS":mg.xg_OCS(PT,melt_wf,species,models),"X":mg.xg_X(PT,melt_wf,species,models),"Xg_t":mg.Xg_tot(PT,melt_wf,species,models),"wt_g":0.}
        wt_C, wt_H, wt_S, wt_X, wt_Fe, wt_O, Wt = new_bulk_regas_open(PT,melt_wf,bulk_wf,gas_mf,dwtg,species,models)
        bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"Fe":wt_Fe,"X":wt_X,"Wt":Wt}
    
    # run over different pressures #
    number_of_step = 0.
    for i in range(initial,final,step): # P is pressure in bars or T is temperature in 'C
        number_of_step = number_of_step + 1.
        eq_Fe = models.loc["eq_Fe","option"]
        
        if models.loc["gassing_style","option"] == "open": # check melt is still vapor-saturated
            if models.loc["insolubles","option"] == "H2O-CO2 only":  
                P_sat_, P_sat_H2O_CO2_result = c.P_sat_H2O_CO2(PT,melt_wf,species,models,psat_tol,nr_step,nr_tol)
                conc = {"wm_H2O":P_sat_H2O_CO2_result["wm_H2O"], "wm_CO2":P_sat_H2O_CO2_result["wm_CO2"], "wm_H2":0., "wm_CO":0., "wm_CH4":0., "wm_H2S":0., "wm_S2m":0., "wm_S6p":0., "ST": 0.}
                frac = c.melt_species_ratios(conc,species)
            else:
                P_sat_, conc, frac = c.P_sat(PT,melt_wf,species,models,psat_tol,nr_step,nr_tol)
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
            xg, conc, melt_and_gas, guesses = eq.mg_equilibrium(PT,melt_wf,bulk_wf,species,models,nr_step,nr_tol,guesses)
            # gas composition
            gas_mf = {"O2":xg["xg_O2"],"CO":xg["xg_CO"],"S2":xg["xg_S2"],"CO2":xg["xg_CO2"],"H2O":xg["xg_H2O"],"H2":xg["xg_H2"],"CH4":xg["xg_CH4"],"SO2":xg["xg_SO2"],"H2S":xg["xg_H2S"],"OCS":xg["xg_OCS"],"X":xg["xg_X"],"Xg_t":xg["Xg_t"],"wt_g":melt_and_gas["wt_g"]}
            if models.loc["gassing_style","option"] == "open" and models.loc["gassing_direction","option"] == "degas": 
                if i == initial:
                    gas_mf_all = gas_mf
                else:
                    gas_mf_all = c.gas_comp_all_open(gas_mf,gas_mf_all,models,species)
            if models.loc["insolubles","option"] == "H2O-CO2 only":
                Fe3T = melt_wf["Fe3FeT"]
                Fe32 = mg.overtotal2ratio(Fe3T)
        else: # NEEDS SORTING ###
            conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
            frac = c.melt_species_ratios(conc,species)
            #wm_ST_ = wm_S_ + wm_S6p_
            #S62 = S6T/S2m_ST
            #Fe3T = melt_wf["Fe3FeT"]
            #Fe32 = mg.overtotal2ratio(Fe3T)
            xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, xg_X, Xg_t, Xm_t, Xm_t_ox, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g = "","","","","","","","","","","","","","",0.,0.,0.,0.,0.
            guesses = eq.initial_guesses(run,PT,melt_wf,setup,species,models,system)

        # set melt composition for forward calculation
        melt_wf["CO2"] = conc["wm_CO2"]
        melt_wf["H2OT"] = conc["wm_H2O"]
        melt_wf["H2"] = conc["wm_H2"]
        melt_wf["CO"] = conc["wm_CO"]
        melt_wf["CH4"] = conc["wm_CH4"]
        melt_wf["H2S"] = conc["wm_H2S"]
        melt_wf["S6+"] = (conc["wm_SO3"]/species.loc["SO3","M"])*species.loc["S","M"]
        melt_wf["S2-"] = conc["wm_S"]
        melt_wf["ST"] = conc["wm_ST"]
        melt_wf["XT"] = conc["wm_X"]
        melt_wf["Fe3FeT"] = conc["Fe3T"]
        if P_sat_ < PT["P"]:  
            bulk_comp = c.bulk_composition(run,PT,melt_wf,setup,species,models)
    
        # check for sulfur saturation and display warning in outputs
        sulf_sat_result = c.sulfur_saturation(PT,melt_wf,species,models)
        if sulf_sat_result["sulfide_sat"] == "yes":
            warning = "WARNING: sulfide-saturated"
        elif sulf_sat_result["sulfate_sat"] == "yes":
            warning = "WARNING: sulfate-saturated"
        else:
            warning = ""
        
        # calculate fO2
        if eq_Fe == "yes":
            fO2_ = mdv.f_O2(PT,melt_wf,species,models)
        elif eq_Fe == "no":
            fO2_ = (xg_O2_*mdv.y_O2(PT,models)*PT["P"])
        
        wm_CO2eq, wm_H2Oeq = mg.melt_H2O_CO2_eq(melt_wf,species)
        melt_comp = mg.melt_normalise_wf(melt_wf,species,"yes","no")
        frac = c.melt_species_ratios(conc,species)

        # store results
        results_headers_table_sample_name, results_values_table_sample_name = results_table_sample_name(setup,run)
        results_headers_table_melt_comp_etc, results_values_table_melt_comp_etc = results_table_melt_comp_etc(PT,melt_comp,conc,frac,melt_wf)
        results_headers_table_model_options, results_values_table_model_options = results_table_model_options(models)    
        results_headers_table_f_p_xg_y_M_C_K_d, results_values_table_f_p_xg_y_M_C_K_d = results_table_f_p_xg_y_M_C_K_d(PT,melt_wf,species,models)
        results_headers_table_sat, results_values_table_sat = results_table_sat(sulf_sat_result,PT,melt_wf,species,models)
        results_values_table_melt_vol = pd.DataFrame([[wm_H2Oeq*100.,wm_CO2eq*1000000.,conc["wm_ST"]*1000000.,melt_wf["XT"]*1000000.]])
        results_values_table_wtg_etc = pd.DataFrame([[melt_and_gas["wt_g"]*100.,melt_and_gas["wt_g_O"],melt_and_gas["wt_g_C"],melt_and_gas["wt_g_H"],melt_and_gas["wt_g_S"],melt_and_gas["wt_g_X"],melt_and_gas["wt_O"]*100.,melt_and_gas["wt_C"]*100.,melt_and_gas["wt_H"]*100.,melt_and_gas["wt_S"]*100.,melt_and_gas["wt_X"]*100.,models.loc["solve_species","option"]]])
        if models.loc["gassing_style","option"] == "open" and models.loc["gassing_direction","option"] == "degas":
            results_values_table_open_all_gas = pd.DataFrame([[gas_mf_all["O2"],gas_mf_all["H2"],gas_mf_all["H2O"],gas_mf_all["S2"],gas_mf_all["SO2"],gas_mf_all["H2S"],gas_mf_all["CO2"],gas_mf_all["CO"],gas_mf_all["CH4"],gas_mf_all["OCS"],gas_mf_all["X"],mg.gas_CS_alt(gas_mf_all)]])
            results1 = pd.concat([results_values_table_sample_name,results_values_table_melt_comp_etc,results_values_table_melt_vol,results_values_table_sat,results_values_table_f_p_xg_y_M_C_K_d,results_values_table_wtg_etc,results_values_table_open_all_gas,results_values_table_model_options],axis=1)
        else:
            results1 = pd.concat([results_values_table_sample_name,results_values_table_melt_comp_etc,results_values_table_melt_vol,results_values_table_sat,results_values_table_f_p_xg_y_M_C_K_d,results_values_table_wtg_etc,results_values_table_model_options],axis=1)
        results = pd.concat([results, results1])
        
        # equilibrium isotope fractionation
        if models.loc["isotopes", "option"] == "yes":
            raise TypeError("This is not currently supported")
            if models.loc["H2S","option"] == "yes":
                print("not currently possible")
            A, B = iso.i2s6("S",PT,R_i,melt_wf,gas_mf,species,i_nr_step,i_nr_tol,guessx)
            RS_Sm, RS_H2S, RS_SO4, RS_S2, RS_SO2, RS_OCS = A
            RS_m, RS_g = B
            a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_ = iso.i2s6_S_alphas(PT)
            xg_SO3_ = 0.
            results2 = pd.DataFrame([[PT["P"],PT["T"],xg_O2_,xg_CO_,xg_CO2_,xg_H2_,xg_H2O_,xg_CH4_,xg_S2_,xg_SO2_,xg_SO3_,xg_H2S_,xg_OCS_,wt_g,
                        wm_CO2_,wm_H2O_,wm_H2_,wm_S_,wm_SO3_,wm_ST_,Fe3T,S6T,
                        mg.fO22Dbuffer(PT,fO2_,"FMQ"),mg.fO22Dbuffer(PT,fO2_,"NNO"),SCSS_,sulfide_sat,SCAS_,sulfate_sat,
                        RS_Sm, RS_SO4, RS_H2S, RS_SO2, RS_S2, RS_OCS, RS_m, RS_g, ratio2delta("VCDT",RS_Sm),ratio2delta("VCDT",RS_SO4),ratio2delta("VCDT",RS_H2S),ratio2delta("VCDT",RS_SO2),ratio2delta("VCDT",RS_S2),ratio2delta("VCDT",RS_OCS),ratio2delta("VCDT",RS_m),ratio2delta("VCDT",RS_g),a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_,RS_g/RS_m]])
            results_isotopes = pd.concat([results_isotopes, results2], ignore_index=True)
            if models.loc["output csv","option"] == "yes":
                results_isotopes.to_csv('results_gassing_isotopes.csv', index=False, header=False)
        
        if models.loc["print status","option"] == "yes":
            if(number_of_step % 100==0):
                print(PT["T"],PT["P"],mg.fO22Dbuffer(PT,fO2_,"FMQ",models),warning,datetime.datetime.now())

        # recalculate bulk composition if needed
        if models.loc["gassing_style","option"] == "open":
            results_me = mg.melt_elements(PT,melt_wf,bulk_wf,gas_mf,species,models)
            if models.loc["gassing_direction","option"] == "degas":
                Wt_ = bulk_wf['Wt']
                if results_me["wm_C"] < 1.e-6:
                    results_me["wm_C"] = 0.
                bulk_wf = {"C":results_me["wm_C"],"H":results_me["wm_H"],"O":results_me["wm_O"],"S":results_me["wm_S"],"X":results_me["wm_X"],"Fe":results_me["wm_Fe"],"Wt":(Wt_*(1. - melt_and_gas["wt_g"]))}
                melt_wf["CT"] = results_me["wm_C"]
                melt_wf["HT"] = results_me["wm_H"]
                melt_wf["ST"] = results_me["wm_S"]
                melt_wf["XT"] = results_me["wm_X"]
            elif models.loc["gassing_direction","option"] == "regas":
                results_nbro = c.new_bulk_regas_open(PT,melt_wf,bulk_wf,gas_mf,dwtg,species,models)
                bulk_wf = {"C":results_nbro["wt_C"],"H":results_nbro["wt_H"],"O":results_nbro["wt_O"],"S":results_nbro["wt_S"],"X":results_nbro["wt_X"],"Fe":results_nbro["wt_Fe"],"Wt":results_nbro["Wt"]}
                melt_wf["CT"] = results_nbro["wm_C"]
                melt_wf["HT"] = results_nbro["wm_H"]
                melt_wf["ST"] = results_nbro["wm_S"]
                melt_wf["XT"] = results_nbro["wm_X"]
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
    
    results.columns = results.iloc[0]
    results = results[1:]  
    if models.loc["output csv","option"] == "yes":
        results.to_csv('results_gassing_chemistry.csv', index=False, header=True)
    
    if models.loc["print status","option"] == "yes":
        print("done", datetime.datetime.now())

    return results


################
### isotopes ###
################

def isotopes(run,isotope_inputs,setup,species,models):

    i_nr_step = gassing_inputs["i_nr_step"]
    i_nr_tol = gassing_inputs["i_nr_tol"]
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
            raise TypeError("input bulk isotope value for system or change 'isotopes' to 'no' in models")
            return
        R_S_SO4_, R_S_S2_ = iso.i2s2("S",PT,R_i,melt_wf)
    else:
        R_i = {"S":""}
   
    
    
def calc_isobar(run,setup,species,models,initial_P,final_P,step_P):
    
    if models.loc["insolubles","option"] == "H2O-CO2 only":
        PT={"T":setup.loc[run,"T_C"]}
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        # set up results table
        results = pd.DataFrame([["P (bar)","H2O (wt%)","CO2 (ppm)"]])
    
        initial_P = int(initial_P)
        final_P = int(final_P+1)
        step_P = int(step_P)
        
        for n in range(initial_P,final_P,step_P):
            PT["P"] = n # pressure in bars
            melt_wf=mg.melt_comp(run,setup)
            results1 = c.calc_isobar_CO2H2O(PT,melt_wf,species,models)
            results = pd.concat([results, results1], ignore_index=True)
            if models.loc["print status","option"] == "yes":
                print(setup.loc[run,"Sample"],n)
    else:
        raise TypeError("insolubles must be H2O-CO2 only")
    if models.loc["output csv","option"] == "yes":
        results.to_csv('results_isobars.csv', index=False, header=False)

    return results
    
def calc_pure_solubility(run,setup,species,models,initial_P):
    if models.loc["print status","option"] == "yes":
        print(setup.loc[run,"Sample"],initial_P)
    PT={"T":setup.loc[run,"T_C"]}

    # check if any options need to be read from the setup file rather than the models file
    models = options_from_setup(run,models,setup)

    # set up results table
    results = pd.DataFrame([["P (bar)","H2O (wt%)","CO2 (ppm)"]])
    
    initial_P = int(initial_P)
        
    for n in range(initial_P,1,-1):
        PT["P"] = n # pressure in bars
        melt_wf=mg.melt_comp(run,setup)
        results1 = c.calc_pure_solubility(PT,melt_wf,species,models)
        results = pd.concat([results, results1], ignore_index=True)
    if models.loc["output csv","option"] == "yes":    
        results.to_csv('results_pure_solubility.csv', index=False, header=False)
    if models.loc["print status","option"] == "yes":
        print("done")

    return results

        
#########################
### Sulfate capacity ###
#########################

# calculate the Csulfate for multiple melt compositions in input file
def Csulfate_output(first_row,last_row,setup,species,models):
    
    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","sulfide saturation","ideal gas","carbonylsulfide","mass_volume","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["sulfur_saturation","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["mass_volume","option"],date.today()]])
    results = pd.concat([results, results1], ignore_index=True)
    results1 = ([["Sample","Pressure (bar)","T ('C)","SiO2","TiO2","Al2O3","FeOT","MnO","MgO","CaO","Na2O","K2O","P2O5",
                "H2O","CO2 (ppm)","ST (ppm)","S6/ST","Fe3/FeT","ln[Csulfide]","ln[Csulfate]","fO2","DNNO","DFMQ"]])
    results = pd.concat([results, results1], ignore_index=True)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        PT={"T":setup.loc[run,"T_C"]}
        melt_wf=mg.melt_comp(run,setup)
        melt_wf['CO2'] = setup.loc[run,"CO2ppm"]/1000000.
        melt_wf["H2OT"] = setup.loc[run,"H2O"]/100.
        melt_wf["ST"] = setup.loc[run,"STppm"]/1000000.
        PT["P"] = setup.loc[run,"P_bar"]
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
        if setup.loc[run,"S6ST"] >= 0.:
            melt_wf["S6ST"] = setup.loc[run,"S6ST"]
        else:
            melt_wf["S6ST"] = ""
        Csulfate_ = mdv.C_SO4(PT,melt_wf,species,models)
                
        ### store results ###
        results2 = pd.DataFrame([[setup.loc[run,"Sample"],setup.loc[run,"P_bar"],setup.loc[run,"T_C"],setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"STppm"],melt_wf["S6ST"],melt_wf["Fe3FeT"],gp.log(mg.C_S(PT,melt_wf,species,models)),gp.log(Csulfate_),mdv.f_O2(PT,melt_wf,species,models),mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"NNO"),mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"FMQ")]])
        results = pd.concat([results, results2], ignore_index=True)                     
        if models.loc["output csv","option"] == "yes":
            results.to_csv('Csulfate.csv', index=False, header=False)
        if models.loc["print status","option"] == "yes":
            print(n, setup.loc[run,"Sample"],gp.log(Csulfate_),gp.log(mg.C_S(PT,melt_wf,species,models)))
        return results

##################
### Capacities ###
##################

# print capacities for multiple melt compositions in input file
def capacities_output(first_row,last_row,setup,species,models):
    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","sulfide saturation","ideal gas","carbonylsulfide","mass_volume","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["sulfur_saturation","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["mass_volume","option"],date.today()]])
    results = pd.concat([results, results1], ignore_index=True)
    results1 = ([["Sample","Pressure (bar)","T ('C)","SiO2","TiO2","Al2O3","FeOT","MnO","MgO","CaO","Na2O","K2O","P2O5",
                "H2O","CO2 (ppm)","ST (ppm)","Fe3+/FeT","fO2 DFMQ","ln[C_CO32-]","ln[C_H2OT]","ln[C_S2-]","ln[C_S6+]","ln[C_H2S]","ln[C_H2]","ln[C_CO]","ln[C_CH4]","ln[C_X]","M_m_SO",]])
    results = pd.concat([results, results1], ignore_index=True)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        PT={"T":setup.loc[run,"T_C"]}
        melt_wf=mg.melt_comp(run,setup)
        CO2 = setup.loc[run,"CO2ppm"]/1000000.
        CT = (CO2/species.loc["CO2","M"])*species.loc["C","M"]
        H2O = setup.loc[run,"H2O"]/100.
        HT = (H2O/species.loc["H2O","M"])*species.loc["H2","M"]
        XT = setup.loc[run,"Xppm"]/1000000.
        ST = setup.loc[run,"STppm"]/1000000.
        melt_wf = {'CO2':CO2,"H2OT":H2O,"ST":ST,"X":XT,"HT":HT,"CT":CT}
        PT["P"] = setup.loc[run,"P_bar"]
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)
        C_CO32 = mdv.C_CO3(PT,melt_wf,species,models)
        C_H2OT = mdv.C_H2O(PT,melt_wf,species,models)
        C_S2 = mdv.C_S(PT,melt_wf,species,models)
        C_S6 = mdv.C_SO4(PT,melt_wf,species,models)
        C_H2S = mdv.C_H2S(PT,melt_wf,species,models)
        C_H2 = mdv.C_H2(PT,melt_wf,species,models)
        C_CO = mdv.C_CO(PT,melt_wf,species,models)
        C_CH4 = mdv.C_CH4(PT,melt_wf,species,models)
        C_X = mdv.C_X(PT,melt_wf,species,models)
        fO2_ = mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"FMQ",models)
        M_m = mg.M_m_SO(melt_wf,species)
                
        ### store results ###
        results2 = pd.DataFrame([[setup.loc[run,"Sample"],setup.loc[run,"P_bar"],setup.loc[run,"T_C"],setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"STppm"],melt_wf["Fe3FeT"],fO2_,gp.log(C_CO32),gp.log(C_H2OT),gp.log(C_S2),gp.log(C_S6),gp.log(C_H2S),gp.log(C_H2),gp.log(C_CO),gp.log(C_CH4),gp.log(C_X),M_m]])
        results = pd.concat([results, results2], ignore_index=True)                    
        if models.loc["output csv","option"] == "yes":
            results.to_csv('capacities.csv', index=False, header=False)
        if models.loc["print status","option"] == "yes":
            print(n, setup.loc[run,"Sample"],gp.log(C_CO32),gp.log(C_H2OT),gp.log(C_S2),gp.log(C_S6),gp.log(C_H2S),gp.log(C_H2),gp.log(C_CO),gp.log(C_CH4),M_m)
        return results

        
##########################
### Fe3+/Fe2+ from fO2 ###        
##########################

def Fe3Fe2_output(first_row,last_row,setup,species,models):
    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","sulfide saturation","ideal gas","carbonylsulfide","mass_volume","insolubles","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["sulfur_saturation","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["mass_volume","option"],models.loc['insolubles','option'],date.today()]])
    results = pd.concat([results, results1], ignore_index=True)
    results1 = ([["Sample","Pressure (bar)","T ('C)","fO2 (DNNO)","fO2 (DFMQ)",
                  "SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)","Fe3/FeT"]])
    results = pd.concat([results, results1], ignore_index=True)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        PT={"T":setup.loc[run,"T_C"]}
        PT["P"] = setup.loc[run,"P_bar"]
        melt_wf=mg.melt_comp(run,setup)
        melt_wf['Fe3FeT'] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
      ### store results ###
        results2 = pd.DataFrame([[setup.loc[run,"Sample"],PT["P"],
                setup.loc[run,"T_C"],mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"NNO"),mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,species,models),"FMQ"),setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],melt_wf['Fe3FeT']]])
        results = pd.concat([results, results2], ignore_index=True)                     
        if models.loc["output csv","option"] == "yes":
            results.to_csv('Fe3FeT_outputs.csv', index=False, header=False)
        if models.loc["print status","option"] == "yes":
            print(n, setup.loc[run,"Sample"],melt_wf['Fe3FeT'])
        return results
        

#############################        
### fugacity coefficients ###
#############################

def fugacity_coefficients(first_row,last_row,setup,species,models):
    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","sulfide saturation","ideal gas","carbonylsulfide","mass_volume","insolubles","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["sulfur_saturation","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["mass_volume","option"],models.loc['insolubles','option'],date.today()]])
    results = pd.concat([results, results1], ignore_index=True)
    results1 = ([["Sample","Pressure (bar)","T ('C)","yH2O","yCO2","yH2","yCO","yO2","yCH4","yS2","ySO2","yH2S","yOCS"]])
    results = pd.concat([results, results1], ignore_index=True)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        PT={"T":setup.loc[run,"T_C"]}
        PT["P"] = setup.loc[run,"P_bar"]
        
        ### store results ###
        results2 = pd.DataFrame([[setup.loc[run,"Sample"],PT["P"],setup.loc[run,"T_C"],mdv.y_H2O(PT,species,models),mdv.y_CO2(PT,species,models),mdv.y_H2(PT,species,models),mdv.y_CO(PT,species,models),mdv.y_O2(PT,species,models),mdv.y_CH4(PT,species,models),mdv.y_S2(PT,species,models),mdv.y_SO2(PT,species,models),mdv.y_H2S(PT,species,models),mdv.y_OCS(PT,species,models)]])
        results = pd.concat([results, results2], ignore_index=True)                     
        if models.loc["output csv","option"] == "yes":
            results.to_csv('fugacity_coefficients_outputs.csv', index=False, header=False)
        if models.loc["print status","option"] == "yes":
            print(n, setup.loc[run,"Sample"],mdv.y_H2O(PT,species,models))
        return results

######################################                
### fO2 range from sulfur content ###
######################################        
        
def fO2_range_from_S_output(first_row,last_row,setup,species,models,p_tol,nr_step,nr_tol):

    # set up results table
    results = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","SCSS","SCAS","ideal gas","carbonylsulfide","insolubles","H2S","species X","species X solubility","Date"]])
    results1 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["SCSS","option"],models.loc["SCAS","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc['insolubles','option'],models.loc["H2S_m","option"],models.loc['species X','option'],models.loc['species X solubility','option'],date.today()]])
    results = pd.concat([results, results1], ignore_index=True)
    results1 = pd.DataFrame([["Sample","T ('C)","wm_STppm","wm_H2OT","wm_CO2ppm","P (bar) sulf","S2- SCSS","sulfide saturated?","DFMQ-sulfide","fO2-sulfide","Fe3FeT-sulfide","S6ST-sulfide","P (bar) anh","S6+ SCAS","sulfate saturated?","DFMQ-sulfate","fO2-sulfate","Fe3FeT-sulfate","S6ST-sulfate"]])
    results = pd.concat([results, results1], ignore_index=True)

    # run over rows in file
    for n in range(first_row,last_row,1): # number of rows in the table
   
        # Select run conditions
        run = n # row number
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        PT = {"T":setup.loc[run,"T_C"]}
        melt_wf=mg.melt_comp(run,setup)
        melt_wf["H2OT"] = (setup.loc[run, "H2O"]/100.)
        melt_wf["CO2"] = setup.loc[run, "CO2ppm"]/1000000.
        melt_wf["HT"] = ((setup.loc[run, "H2O"]/100.)/species.loc["H2O","M"])*species.loc["H2","M"]
        if setup.loc[run,"P_bar"] > 0.:
            PT["P"] = setup.loc[run,"P_bar"]
            P_sat_sulf, P_sat_anh = setup.loc[run,"P_bar"],setup.loc[run,"P_bar"]
            melt_wf["Fe3FeT"] = setup.loc[run, "Fe3FeT"]
            SCAS_,SCSS_,sulfide_sat,DFMQ_1,fO2_1,Fe3FeT_1,S6ST_1,sulfate_sat,DFMQ_2,fO2_2,Fe3FeT_2,S6ST_2 = c.fO2_range_from_S(PT,melt_wf,species,models)
        else:
            melt_wf["ST"] = setup.loc[run, "STppm"]/1000000.
            melt_wf["XT"] = setup.loc[run, "Xppm"]/1000000.
            melt_wf["CT"] = ((setup.loc[run, "CO2ppm"]/1000000.)/species.loc["CO2","M"])*species.loc["C","M"]
            P_sat_sulf,P_sat_anh,SCAS_,SCSS_,sulfide_sat,DFMQ_1,fO2_1,Fe3FeT_1,S6ST_1,sulfate_sat,DFMQ_2,fO2_2,Fe3FeT_2,S6ST_2 = c.P_sat_sulf_anh(PT,melt_wf,species,models,p_tol,nr_step,nr_tol)

        # store results
        results1 = pd.DataFrame([[setup.loc[run,"Sample"],PT["T"],setup.loc[run,"STppm"],setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],P_sat_sulf,SCSS_,sulfide_sat,DFMQ_1,fO2_1,Fe3FeT_1,S6ST_1,P_sat_anh,SCAS_,sulfate_sat,DFMQ_2,fO2_2,Fe3FeT_2,S6ST_2]])
        results = pd.concat([results, results1], ignore_index=True)
        if models.loc["output csv","option"] == "yes":
            results.to_csv('fO2_range_from_S.csv', index=False, header=False)
        if models.loc["print status","option"] == "yes":
            print(setup.loc[run,"Sample"])
        return results

############################################
### melt-vapour equilibrium at given fO2 ###
############################################

def eq_given_fO2(inputs,setup,species,models): # only S atm
    
    option = inputs["option"]
    
    # set T, volatile composition of the melt, and tolerances
    #nr_step = inputs["nr_step"]
    #nr_tol = inputs["nr_tol"]
    
    if option == "loop":
        run = inputs["run"]
        PT={'T':setup.loc[run,"T_C"]}
        melt_wf=mg.melt_comp(run,setup)
        melt_wf['CO2'] = setup.loc[run,"CO2ppm"]/1000000.
        melt_wf["H2OT"]  = setup.loc[run,"H2O"]/100.
        melt_wf["ST"] = setup.loc[run,"STppm"]/1000000.
        melt_wf["H2"] = 0.
        step = inputs["dfO2_step"]
        initial = inputs["fO2_i"]
        final = inputs["fO2_f"]
        start = 0
        end = inputs["no_steps"]
        difference = (final - initial) + 1
        step_size = difference/end
        # Set bulk composition
        wt_C = (melt_wf["CO2"]/species.loc["CO2","M"])*species.loc["C","M"]
        wt_H = (melt_wf["H2OT"]/species.loc["H2O","M"])*(2.*species.loc["H","M"])
        wt_S = melt_wf["ST"]
        bulk_wf = {"C":wt_C,"H":wt_H,"S":wt_S}
        system = eq.set_system(bulk_wf,models)
        
    elif option == "spreadsheet":
        start = inputs["first row"]
        end = inputs["last row"]

    # update Fe3+/FeT and water speciation at saturation pressure and check for sulfur saturation
    #melt_wf["Fe3FeT"] = mg.fO22Fe3FeT(10.**initial,run,PT,setup,species,models)
    #melt_wf["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
    #melt_wf["S6ST"] = mg.S6ST(run,PT,melt_wf,setup,species,models)
    #melt_wf["H2Omol"] = 0. #mg.wm_H2Omol(run,PT,melt_wf,setup,species,models)
    #melt_wf["OH"] = 0., #mg.wm_OH(run,PT,melt_wf,setup,species,models)
    #SCSS_,sulfide_sat,SCAS_,sulfate_sat,ST_ = sulfur_saturation(run,PT,melt_wf,setup,species,models)

        
    # create results table
    results_header1 = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfur speciation","ideal gas","carbonylsulfide","bulk composition","equilibrate Fe","starting pressure","gassing direction","gassing style","mass_volume","crystallisation","isotopes","Date"]])
    results_header2 = pd.DataFrame([[
models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["bulk_composition","option"],models.loc["eq_Fe","option"],models.loc["starting_P","option"],models.loc["gassing_direction","option"],models.loc["gassing_style","option"],models.loc["mass_volume","option"],models.loc["crystallisation","option"],models.loc["isotopes","option"],date.today()]])
    results_header = pd.concat([results_header1, results_header2], ignore_index=True)
    results_chemistry1 = pd.DataFrame([["P","T('C)","System","run","Sample","H2O (mwf)","CO2 (mwf)","ST (mwf)","Fe3FeT","S6ST","O (twf)","C (twf)","H (wtf)","S (twf)","Fe (twf)",
"SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)","xg_O2","xg_CO","xg_CO2","xg_H2","xg_H2O","xg_CH4","xg_S2","xg_SO2","xg_SO3","xg_H2S","xg_OCS","Xg_t",
               "xm_CO2","xm_H2O","xm_H2Omol","xm_OH","xm_H2","Xm_t_SO","Xm_t_ox",
               "wm_CO2","wm_H2O","wm_H2Omol","wm_OH","wm_H2","wm_S","wm_SO3","wm_ST","Fe3T","S6T",
               "DFMQ","DNNO","SCSS","sulfide sat?","SCAS","sulfate sat?","wt_g","wt_O","wt_C","wt_H","wt_S",
               "fO2","fH2","fH2O","fS2","fSO2","fSO3","fH2S","fCO2","fCO","fCH4","fOCS",
               "yO2","yH2","yH2O","yS2","ySO2","ySO3","yH2S","yCO2","yCO","yCH4","yOCS",
               "M_m_SO","M_m_ox","C_H2O","C_H2","C_CO3","C_S","C_SO4",
               "KD1","KHOg","KHOm","KHOSg","KCOg","KCOHg","KOCSg","KSOg","KSOg2"]])
    results_chemistry = pd.concat([results_header, results_chemistry1], ignore_index=True)
    
    
    # run over different fO2 #
    for i in range(start,end,1): 
        if option == "loop":
            i_ = (i*step_size)+initial # fO2 in log units 
            fO2_ = 10.**i_
            PT["P"] = inputs["P"]
            if models.loc["print status","option"] == "yes":
                print(setup.loc[run,"Sample"],fO2_)
        elif option == "spreadsheet":
            run = i
            PT={'T':setup.loc[run,"T_C"]}
            melt_wf=mg.melt_comp(run,setup)
            melt_wf['CO2'] = setup.loc[run,"CO2ppm"]/1000000.
            melt_wf["H2OT"]  = setup.loc[run,"H2O"]/100.
            melt_wf["ST"] = setup.loc[run,"STppm"]/1000000.
            melt_wf["H2"] = 0.
            PT["P"] = setup.loc[run,"P_bar"]
            melt_wf['Fe3FeT'] = mg.Fe3FeT_i(PT,melt_wf,species,models)
            fO2_ = mdv.f_O2(PT,melt_wf,species,models)
            # Set bulk composition
            wt_C = (melt_wf["CO2"]/species.loc["CO2","M"])*species.loc["C","M"]
            wt_H = (melt_wf["H2OT"]/species.loc["H2O","M"])*(2.*species.loc["H","M"])
            wt_S = melt_wf["ST"]
            bulk_wf = {"C":wt_C,"H":wt_H,"S":wt_S,"CT":wt_C,"HT":wt_H,"ST":wt_S}
            system = eq.set_system(bulk_wf,models)
            if models.loc["print status","option"] == "yes":
                print(setup.loc[run,"Sample"],fO2_,PT["P"])
        
        # work out equilibrium partitioning between melt and gas phase
        melt_wf["Fe3FeT"] = mg.fO22Fe3FeT(fO2_,PT,species,models)
        Fe3T = melt_wf["Fe3FeT"]
        #xg_S2_, result1, result2, result3 = eq.eq_SOFe_fO2(run,PT,10.**i_,bulk_wf,melt_wf,species,setup,models,nr_step,nr_tol,guessx)
        #xg_SO2_, xg_O2_, Xg_t, Fe32, Fe3T, wm_S_, wm_SO3_, S62, S6T, wm_ST_ = result1
        #wt_g, wt_O_, wt_S_ = result3
        #guessx = xg_S2_
        #xg_O2_, xg_S2_, xg_SO2_, wm_S_, wm_SO3_, wm_ST_, S6T, Xg_t = eq.S_P_fO2_1(run,PT,fO2_,melt_wf,setup,species,models)
        #wm_CO2_, xm_CO2_, wm_H2O_, xm_H2O_, wm_H2_, xm_H2_, wm_H2Omol_, xm_H2Omol_, wm_OH_, xm_OH_, xg_CO_, xg_CO2_, xg_H2_, xg_H2O_, xg_CH4_, xg_H2S_, xg_OCS_, xg_SO3_, Xm_t, Xm_t_ox, wt_C_, wt_H_, wt_g = 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
        xg_O2_, xg_S2_, xg_SO2_, xg_SO3_, xg_H2_, xg_H2O_, xg_H2S_, xg_CO_, xg_CO2_, xg_CH4_, xg_OCS_, wm_S_, wm_SO3_, wm_ST_, S6T, wm_H2O_, xm_H2O_, wm_H2_, xm_H2_, wm_H2Omol_, xm_H2Omol_, wm_OH_, xm_OH_, wm_CO2_, xm_CO2_, Xg_t, Xm_t, Xm_t_ox, wt_C_, wt_H_, wt_S_, wt_g = eq.S_P_fO2(PT,fO2_,melt_wf,species,models)
                
        # set melt composition for forward calculation
        melt_wf = {"CO2":wm_CO2_,"H2OT":wm_H2O_,"ST":wm_ST_,"S6ST":S6T,"Fe3FeT":Fe3T,"H2":wm_H2_,"H2Omol":wm_H2Omol_,"OH":wm_OH_,"S2-":wm_S_}
    
        # check for sulfur saturation and display warning in outputs
        SCSS_,sulfide_sat,SCAS_,sulfate_sat, ST_ = sulfur_saturation(PT,melt_wf,species,models)
        if sulfide_sat == "yes":
            warning = "WARNING: sulfide-saturated"
        elif sulfate_sat == "yes":
            warning = "WARNING: sulfate-saturated"
        else:
            warning = ""

        # store results
               
        results2 = pd.DataFrame([[PT["P"],PT["T"],system,run,setup.loc[run,"Sample"],melt_wf["H2OT"],melt_wf["CO2"],melt_wf["ST"],melt_wf["Fe3FeT"],"SORT","SORT",bulk_wf["C"],bulk_wf["H"],bulk_wf["S"],"SORT",
setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(melt_wf,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],xg_O2_,xg_CO_,xg_CO2_,xg_H2_,xg_H2O_,xg_CH4_,xg_S2_,xg_SO2_,xg_SO3_,xg_H2S_,xg_OCS_,Xg_t,
               xm_CO2_,xm_H2O_,xm_H2Omol_,xm_OH_,xm_H2_,Xm_t,Xm_t_ox,
               wm_CO2_,wm_H2O_,wm_H2Omol_,wm_OH_,wm_H2_,wm_S_,wm_SO3_,wm_ST_,Fe3T,S6T,
               mg.fO22Dbuffer(PT,fO2_,"FMQ"),mg.fO22Dbuffer(PT,fO2_,"NNO"),SCSS_,sulfide_sat,SCAS_,sulfate_sat,
               wt_g,"SORT",wt_C,wt_H,wt_S,
fO2_,mg.f_H2(PT,melt_wf,species,models),mg.f_H2O(PT,melt_wf,species,models),mg.f_S2(PT,melt_wf,species,models),mg.f_SO2(PT,melt_wf,species,models),mg.f_SO3(PT,melt_wf,species,models),mg.f_H2S(PT,melt_wf,species,models),mg.f_CO2(PT,melt_wf,species,models),mg.f_CO(PT,melt_wf,species,models),mg.f_CH4(PT,melt_wf,species,models),mg.f_OCS(PT,melt_wf,species,models),
mg.y_O2(PT,species,models),mg.y_H2(PT,species,models),mg.y_H2O(PT,species,models),mg.y_S2(PT,species,models),mg.y_SO2(PT,species,models),mg.y_SO3(PT,species,models),mdv.y_H2S(PT,species,models),mg.y_CO2(PT,species,models),mg.y_CO(PT,species,models),mg.y_CH4(PT,species,models),mg.y_OCS(PT,species,models),
mg.M_m_SO(melt_wf,species),mg.M_m_ox(melt_wf,species,models),mg.C_H2O(PT,melt_wf,species,models),mg.C_H2(PT,melt_wf,species,models),mdv.C_CO3(PT,melt_wf,species,models),mg.C_S(PT,melt_wf,species,models),mg.C_SO4(PT,melt_wf,species,models),
mg.KD1(PT,melt_wf,species,models),mg.KHOg(PT,models),mg.KHOm(PT,melt_wf,species,models),mg.KHOSg(PT,models),mg.KCOg(PT,models),mg.KCOHg(PT,models),mg.KOCSg(PT,models),mg.KOSg(PT,models),mg.KOSg2(PT,models)]])
        results_chemistry = pd.concat([results_chemistry, results2], ignore_index=True)
        if models.loc["output csv","option"] == "yes":
            results_chemistry.to_csv('results_fO2_chemistry.csv', index=False, header=False)
        return results_chemistry

##############################################
### fO2 of silm+sulfm+anh at given T and P ###
##############################################

def fO2_SSA_output(first_row,last_row,setup,species,models):
    # set up results table
    results = pd.DataFrame([["Sample","P (bar)","T ('C)","fO2 (DFMQ)",
                  "SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)",
                "H2OT (wt%)","CO2 (ppm)","ST (ppm)","SCSS (ppm)","SCAS (ppm)", "S6+/ST","Fe3+/FeT"]])

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        PT={"T":setup.loc[run,"T_C"], "P":setup.loc[run,"P_bar"]}
        melt_wf = {'CO2':setup.loc[run,"CO2ppm"]/1000000.,"H2OT":setup.loc[run,"H2O"]/100.,"Fe3FeT":0.}
        
        fO2, DFMQ, wmST, S6ST, S6, S2 = c.fO2_silm_sulf_anh(run,PT,melt_wf,setup,species,models)
        
        Fe3FeT_n0 = melt_wf["Fe3FeT"]
        Fe3FeT_n1 = mdv.fO22Fe3FeT(fO2,run,PT,setup,species,models)
                            
        if ((Fe3FeT_n0 - Fe3FeT_n1)**2)**0.5 > 0.01:
            melt_wf["Fe3FeT"] = Fe3FeT_n1
            fO2, DFMQ, wmST, S6ST, S6, S2 = c.fO2_silm_sulf_anh(run,PT,melt_wf,setup,species,models)
            Fe3FeT_n0 = melt_wf["Fe3FeT"]
            Fe3FeT_n1 = mdv.fO22Fe3FeT(fO2,run,PT,setup,species,models)
       
        ### store results ###
        results2 = pd.DataFrame([[setup.loc[run,"Sample"],PT["P"],setup.loc[run,"T_C"],DFMQ,setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],
                setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],wmST, S2, S6, S6ST, melt_wf["Fe3FeT"]]])
        results = pd.concat([results, results2], ignore_index=True)
        if models.loc["output csv","option"] == "yes":
            results.to_csv('fO2_silm_sulf_anh.csv', index=False, header=False)
        if models.loc["print status","option"] == "yes":
            print(n, setup.loc[run,"Sample"],PT["P"])
        return results
        
##############################################
### S content at given T, P, fO2, C, and H ###
##############################################
            
def S_given_T_P_fO2_C_H_output(first_row,last_row,setup,species,models,nr_step,nr_tol):
                
    if models.loc["H2S_m","option"] != "no":
        raise TypeError("This calculation assumes H2S is insoluble in the melt")

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        PT={"T":setup.loc[run,"T_C"], "P":setup.loc[run,"P_bar"]}
        melt_wf=mg.melt_comp(run,setup)
        melt_wf['CO2']= setup.loc[run,"CO2ppm"]/1000000.
        melt_wf["H2OT"] = setup.loc[run,"H2O"]/100.
        melt_wf["XT"] = setup.loc[run,"Xppm"]/1000000.
        melt_wf["ST"] = 0.
        melt_wf["CT"] = (melt_wf["CO2"]*species.loc['C','M'])/species.loc['CO2','M']
        melt_wf["HT"] = (melt_wf["H2OT"]*species.loc['H','M'])/species.loc['H2O','M']
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,species,models)                                              
        
        conc, frac = c.S_given_T_P_fO2_C_H(PT,melt_wf,species,models,nr_step,nr_tol)
        melt_wf["ST"] = conc["wm_ST"]
        melt_wf["S2-"] = conc["wm_S2m"]
        melt_comp = mg.melt_normalise_wf(melt_wf,species,"yes","no")
        sulf_sat_result = c.sulfur_saturation(PT,melt_wf,species,models)

        # create results
        results_headers_table_sample_name, results_values_table_sample_name = results_table_sample_name(setup,run)
        results_headers_table_melt_comp_etc, results_values_table_melt_comp_etc = results_table_melt_comp_etc(PT,melt_comp,conc,frac,melt_wf)
        results_headers_table_model_options, results_values_table_model_options = results_table_model_options(models)    
        results_headers_table_f_p_xg_y_M_C_K_d, results_values_table_f_p_xg_y_M_C_K_d = results_table_f_p_xg_y_M_C_K_d(PT,melt_wf,species,models)
        results_headers_table_sat, results_values_table_sat = results_table_sat(sulf_sat_result,PT,melt_wf,species,models)
        results_headers_table_melt_vol = results_table_melt_vol() # "H2OT-eq_wtpc","CO2T-eq_ppmw","ST_ppmw","X_ppmw"
        results_values_table_melt_vol = pd.DataFrame([[setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],conc["wm_ST"]*1000000.,setup.loc[run,"Xppm"]]])
        results_headers = pd.concat([results_headers_table_sample_name,results_headers_table_melt_comp_etc,results_headers_table_melt_vol,results_headers_table_sat,results_headers_table_f_p_xg_y_M_C_K_d,results_headers_table_model_options],axis=1)
        results1 = pd.concat([results_values_table_sample_name,results_values_table_melt_comp_etc,results_values_table_melt_vol,results_values_table_sat,results_values_table_f_p_xg_y_M_C_K_d,results_values_table_model_options],axis=1)
    
        if n == first_row:
            results = pd.concat([results_headers, results1])
        else:                         
            results = pd.concat([results, results1])
        
        if models.loc["print status","option"] == "yes":
            print(n, setup.loc[run,"Sample"],PT["P"])
    
    results.columns = results.iloc[0]
    results = results[1:]  
    if models.loc["output csv","option"] == "yes":
        results.to_csv('results_S_given_T_P_fO2_C_H.csv', index=False, header=True)

    return results
        
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
                             
    results = pd.concat([results, results1], ignore_index=True)
    results1 = pd.DataFrame([["sds","",setup.loc[run,"SiO2_sd"],setup.loc[run,"TiO2_sd"],setup.loc[run,"Al2O3_sd"],setup.loc[run,"FeOT_sd"],setup.loc[run,"MnO_sd"],setup.loc[run,"MgO_sd"],setup.loc[run,"CaO_sd"],setup.loc[run,"Na2O_sd"],setup.loc[run,"K2O_sd"],setup.loc[run,"P2O5_sd"],
                setup.loc[run,"H2O_sd"],setup.loc[run,"CO2ppm_sd"],setup.loc[run,"STppm_sd"], setup.loc[run,"Fe3FeT_sd"]]])
                             
    results = pd.concat([results, results1], ignore_index=True)
    results1 = pd.DataFrame([["sd types","",setup.loc[run,"SiO2_sd_type"],setup.loc[run,"TiO2_sd_type"],setup.loc[run,"Al2O3_sd_type"],setup.loc[run,"FeOT_sd_type"],setup.loc[run,"MnO_sd_type"],setup.loc[run,"MgO_sd_type"],setup.loc[run,"CaO_sd_type"],setup.loc[run,"Na2O_sd_type"],setup.loc[run,"K2O_sd_type"],setup.loc[run,"P2O5_sd_type"],
                setup.loc[run,"H2O_sd_type"],setup.loc[run,"CO2ppm_sd_type"],setup.loc[run,"STppm_sd_type"], setup.loc[run,"Fe3FeT_sd_type"]]])
                             
    results = pd.concat([results, results1], ignore_index=True)
    for n in range(0,iterations,1): # n is number of rows of data in conditions file
        SiO2,TiO2,Al2O3,FeOT,MnO,MgO,CaO,Na2O,K2O,P2O5,H2O,CO2ppm,STppm,Fe3FeT = c.compositions_within_error(run,setup)
        results1 = pd.DataFrame([[run,setup.loc[run,"T_C"],SiO2,TiO2,Al2O3,FeOT,MnO,MgO,CaO,Na2O,K2O,P2O5,H2O,CO2ppm,STppm,Fe3FeT]])
        results = pd.concat([results, results1], ignore_index=True)
        if models.loc["output csv","option"] == "yes":
            results.to_csv('random_compositions.csv', index=False, header=False)
        if models.loc["print status","option"] == "yes":
            print(n, setup.loc[run,"Sample"],SiO2)
    
    return results
        
##################################################################################
### vapor undersaturated cooling for S2-, S6+, Fe3+, Fe2+, H2OT, CO32- cooling ###
##################################################################################

def cooling(run,cooling_inputs,setup,species,models):
    
    if models.loc["print status","option"] == "yes":
        print(setup.loc[run,"Sample"])
    
    # check if any options need to be read from the setup file rather than the models file
    models = options_from_setup(run,models,setup)

    # set T, volatile composition of the melt, and tolerances
    PT={"P":setup.loc[run,"P_bar"]}
    PT["T"]=setup.loc[run,"T_C"]
    melt_wf = {'CO2':setup.loc[run,"CO2ppm"]/1000000.,"H2OT":setup.loc[run,"H2O"]/100.,"H2":0.}
    melt_wf["CT"] = (melt_wf["CO2"]/species.loc["CO2","M"])*species.loc["C","M"]
    melt_wf["HT"] = (2.*melt_wf["H2OT"]/species.loc["H2O","M"])*species.loc["H","M"]
    nr_step = cooling_inputs["nr_step"]
    nr_tol = cooling_inputs["nr_tol"]
    dt_step = cooling_inputs["dt_step"]
    psat_tol = cooling_inputs["psat_tol"]
        
    # Calculate S and Fe composition at initial T and check for sulfur saturation
    melt_wf["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
    Fe3Fe2_ = mg.Fe3Fe2(melt_wf)
    melt_wf["CO"] = 0.
    melt_wf["CH4"] = 0.
    melt_wf["H2S"] = 0. 
    wm_H2_, wm_CH4_, wm_H2S_, wm_CO_ = 0., 0., 0., 0.
    
    if models.loc["bulk_O","option"] == "inc S":
        melt_wf["ST"] = setup.loc[run,"STppm"]/1000000.
        melt_wf["S6ST"] = mg.S6ST(run,PT,melt_wf,setup,species,models)
        melt_wf["S2-"] = (1.-melt_wf["S6ST"])*melt_wf["ST"]
        melt_wf["S6+"] = melt_wf["S6ST"]*melt_wf["ST"]     
    elif models.loc["bulk_O","option"] == "exc S":
        melt_wf["ST"] = 0.
        melt_wf["S2-"] = 0.
        melt_wf["S6+"] = 0.
    
    SCSS_,sulfide_sat,SCAS_,sulfate_sat,ST_ = c.sulfur_saturation(run,PT,melt_wf,setup,species,models)
    wm_S2m_ = melt_wf["S2-"]
    wm_S6p_ = melt_wf["S6+"]
    wm_SO3_ = (wm_S6p_/(species.loc["S","M"]))*(species.loc["S","M"] + 3.*species.loc["O","M"])
        
    # Set bulk composition
    wt_C, wt_O, wt_H, wt_S, wt_Fe, wt_g_, Wt_ = c.bulk_composition(run,PT,melt_wf,setup,species,models)
    
    if models.loc["bulk_O","option"] == "exc S":
        melt_wf["ST"] = setup.loc[run,"STppm"]/1000000.           
        melt_wf["S6ST"] = mg.S6ST(run,PT,melt_wf,setup,species,models)
        melt_wf["S2-"] = (1.-melt_wf["S6ST"])*melt_wf["ST"]
        melt_wf["S6+"] = melt_wf["S6ST"]*melt_wf["ST"]     
        wt_S = melt_wf["ST"]
    
    bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"Fe":wt_Fe,"Wt":Wt_}
        
    # set system and initial guesses
    system = eq.set_system(melt_wf,models)
    guessx = mdv.f_O2(run,PT,melt_wf,setup,species,models)

    if models.loc["bulk_O","option"] == "exc S":
        fO2_, A, B, C = eq.eq_SOFe_melt(run,PT,bulk_wf,melt_wf,species,setup,models,nr_step/10.,nr_tol,guessx)
        Fe32, Fe3T, S62, S6T = A
        melt_wf["S6ST"] = S6T
        melt_wf["Fe3FeT"] = Fe3T          
        melt_wf["S2-"] = (1.-melt_wf["S6ST"])*melt_wf["ST"]
        melt_wf["S6+"] = melt_wf["S6ST"]*melt_wf["ST"]
        guessx = mdv.f_O2(run,PT,melt_wf,setup,species,models)
    
    # create results table
    results_header1 = pd.DataFrame([["System","run","Sample","H2O (mwf)","CO2 (mwf)","ST (mwf)","Fe3FeT","S6ST","O (twf)","C (twf)","H (wtf)","S (twf)","Fe (twf)",
"SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)",
"oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","ideal gas","carbonylsulfide","bulk composition","equilibrate Fe","starting pressure","gassing direction","gassing style","mass_volume","crystallisation","isotopes","Date"]])
    results_header2 = pd.DataFrame([[system,run,setup.loc[run,"Sample"],melt_wf["H2OT"],melt_wf["CO2"],melt_wf["ST"],melt_wf["Fe3FeT"],"SORT",bulk_wf["O"],bulk_wf["C"],bulk_wf["H"],bulk_wf["S"],bulk_wf["Fe"],
setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],
models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["bulk_composition","option"],models.loc["eq_Fe","option"],models.loc["starting_P","option"],models.loc["gassing_direction","option"],models.loc["gassing_style","option"],models.loc["mass_volume","option"],models.loc["crystallisation","option"],models.loc["isotopes","option"],date.today()]])
    results_header = pd.concat([results_header1, results_header2], ignore_index=True)
    results_chemistry1 = pd.DataFrame([["P","T('C)","xg_O2","xg_CO","xg_CO2","xg_H2","xg_H2O","xg_CH4","xg_S2","xg_SO2","xg_H2S","xg_OCS","Xg_t",
               "xm_CO2","xm_H2O","Xm_t_SO","Xm_t_ox",
               "wm_CO2","wm_H2O","wm_H2","wm_CO","wm_CH4","wm_S","wm_SO3","wm_H2S","wm_ST","Fe32","Fe3T","S62","S6T",
               "DFMQ","DNNO","SCSS","sulfide sat?","SCAS","sulfate sat?","wt_g","wt_g_O","wt_g_C","wt_g_H","wt_g_S","wt_O","wt_C","wt_H","wt_S",
               "fO2","fH2","fH2O","fS2","fSO2","fH2S","fCO2","fCO","fCH4","fOCS",
               "yO2","yH2","yH2O","yS2","ySO2","yH2S","yCO2","yCO","yCH4","yOCS",
               "M_m_SO","M_m_ox","C_H2O","C_H2","C_CO3","C_CO","C_CH4","C_S","C_SO4","C_H2S",
               "KD1","KHOg","KHOm","KHOSg","KCOg","KCOHg","KOCSg","KSOg","KSOg2"]])
    results_chemistry = pd.concat([results_header, results_chemistry1], ignore_index=True)
    results1 = pd.DataFrame([[PT["P"],PT["T"],mg.xg_O2(run,PT,melt_wf,setup,species,models),mg.xg_CO(run,PT,melt_wf,setup,species,models),mg.xg_CO2(run,PT,melt_wf,setup,species,models),mg.xg_H2(run,PT,melt_wf,setup,species,models),mg.xg_H2O(run,PT,melt_wf,setup,species,models),mg.xg_CH4(run,PT,melt_wf,setup,species,models),mg.xg_S2(run,PT,melt_wf,setup,species,models),mg.xg_SO2(run,PT,melt_wf,setup,species,models),mg.xg_H2S(run,PT,melt_wf,setup,species,models),mg.xg_OCS(run,PT,melt_wf,setup,species,models),mg.Xg_tot(run,PT,melt_wf,setup,species,models),
mg.xm_CO2_so(run,melt_wf,setup,species),mg.xm_H2OT_so(run,melt_wf,setup,species),mg.Xm_t_so(run,melt_wf,setup,species),mg.Xm_t_ox(run,melt_wf,setup,species),
melt_wf["CO2"],melt_wf["H2OT"],wm_H2_,wm_CO_,wm_CH4_,wm_S2m_,wm_SO3_,wm_H2S_,melt_wf["ST"],mg.Fe3Fe2(melt_wf),melt_wf["Fe3FeT"],mg.S6S2(run,PT,melt_wf,setup,species,models),mg.S6ST(run,PT,melt_wf,setup,species,models),
mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"FMQ",models),mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"NNO",models),SCSS_,sulfide_sat,SCAS_,sulfate_sat,wt_g_,"","","","",
bulk_wf["O"],bulk_wf["C"],bulk_wf["H"],bulk_wf["S"],
mdv.f_O2(run,PT,melt_wf,setup,species,models),mg.f_H2(run,PT,melt_wf,setup,species,models),mg.f_H2O(run,PT,melt_wf,setup,species,models),mg.f_S2(run,PT,melt_wf,setup,species,models),mg.f_SO2(run,PT,melt_wf,setup,species,models),mg.f_H2S(run,PT,melt_wf,setup,species,models),mg.f_CO2(run,PT,melt_wf,setup,species,models),mg.f_CO(run,PT,melt_wf,setup,species,models),mg.f_CH4(run,PT,melt_wf,setup,species,models),mg.f_OCS(run,PT,melt_wf,setup,species,models),
mdv.y_O2(PT,species,models),mdv.y_H2(PT,species,models),mdv.y_H2O(PT,species,models),mdv.y_S2(PT,species,models),mdv.y_SO2(PT,species,models),mdv.y_SO3(PT,species,models),mdv.y_H2S(PT,species,models),mdv.y_CO2(PT,species,models),mdv.y_CO(PT,species,models),mdv.y_CH4(PT,species,models),mdv.y_OCS(PT,species,models),
mg.M_m_SO(run,melt_wf,setup,species),mg.M_m_ox(run,melt_wf,setup,species,models),mdv.C_H2O(run,PT,melt_wf,setup,species,models),mdv.C_H2(run,PT,melt_wf,setup,species,models),mdv.C_CO3(run,PT,melt_wf,setup,species,models),mdv.C_CO(run,PT,melt_wf,setup,species,models),mdv.C_CH4(run,PT,melt_wf,setup,species,models),mdv.C_S(run,PT,melt_wf,setup,species,models),mdv.C_SO4(run,PT,melt_wf,setup,species,models),mdv.C_H2S(run,PT,melt_wf,setup,species,models),
mg.KD1(run,PT,setup,species,models),mdv.KHOg(PT,models),mdv.KHOm(run,PT,melt_wf,setup,species,models),mdv.KHOSg(PT,models),mdv.KCOg(PT,models),mdv.KCOHg(PT,models),mdv.KOCSg(PT,models),mdv.KOSg(PT,models),mdv.KOSg2(PT,models)]])
    results_chemsitry = pd.concat([results_chemistry, results1], ignore_index=True)
    if models.loc["output csv","option"] == "yes":
        results_chemistry.to_csv('results_cooling_chemistry.csv', index=False, header=False)
    
    if models.loc["bulk_O","option"] == "inc S":
        if setup.loc[run,"final_T"] > setup.loc[run,"T_C"]:
            initial = int(round(PT["T"])+1) 
            step = 1 # temperature step in 'C
        elif setup.loc[run,"final_T"] < setup.loc[run,"T_C"]:
            initial = int(round(PT["T"])-1)
            final = int(round(setup.loc[run,"final_T"]))
            step = -1 # temperature step in 'C
    elif models.loc["bulk_O","option"] == "exc S":
        initial = 10
        final = 5000
        step = 10
    
    # run over different temperatures #
    for i in range(initial,final,step): # P is pressure in bars or T is temperature in 'C
        eq_Fe = models.loc["eq_Fe","option"]
        if models.loc["bulk_O","option"] == "inc S":
            T = i/dt_step
            PT["T"] = T
        elif models.loc["bulk_O","option"] == "exc S":
            PT["T"] = setup.loc[run,"T_C"]
            melt_wf["ST"]=i/1000000.
            bulk_wf["S"]=i/1000000.
        
        # work out equilibrium partitioning between melt and gas phase
        fO2_, A, B, C = eq.eq_SOFe_melt(run,PT,bulk_wf,melt_wf,species,setup,models,nr_step,nr_tol,guessx)
        Fe32, Fe3T, S62, S6T = A
        
        # set melt composition for forward calculation
        wm_CO2_ = melt_wf["CO2"]
        wm_H2O_ = melt_wf["H2OT"]
        wm_H2_, wm_CO_, wm_CH4_, wm_H2S_ = "","","","" 
        melt_wf["Fe3FeT"] = Fe3T
        melt_wf["S6ST"] = S6T
        if models.loc["bulk_O","option"] == "inc S":
            wm_ST_ = wt_S
        elif models.loc["bulk_O","option"] == "exc S":
            wm_ST_ = melt_wf["ST"]
            wm_S = wm_ST_
        melt_wf["S2-"] = wt_S*(1.-S6T)    
        wm_S_ = wt_S*(1.-S6T)
        wm_SO3_ = ((wt_S*S6T)/species.loc["S","M"])*species.loc["SO3","M"]
    
        # check for sulfur saturation and display warning in outputs
        SCSS_,sulfide_sat,SCAS_,sulfate_sat, ST_ = sulfur_saturation(run,PT,melt_wf,setup,species,models)
        if sulfide_sat == "yes":
            warning = "WARNING: sulfide-saturated"
        elif sulfate_sat == "yes":
            warning = "WARNING: sulfate-saturated"
        else:
            warning = ""
        
        # calculate fO2
        if eq_Fe == "yes":
            fO2_ = mdv.f_O2(run,PT,melt_wf,setup,species,models)
        elif eq_Fe == "no":
            fO2_ = (xg_O2_*mg.y_O2(PT,models)*PT["P"])
            
        guessx = fO2_
            
        # volume, density, and mass
        #gas_mf = {"O2":xg_O2_,"CO":xg_CO_,"S2":xg_S2_,"CO2":xg_CO2_,"H2O":xg_H2O_,"H2":xg_H2_,"CH4":xg_CH4_,"SO2":xg_SO2_,"H2S":xg_H2S_,"OCS":xg_OCS_,"Xg_t":Xg_t,"wt_g":wt_g}
        xg_O2_,xg_CO_,xg_CO2_,xg_H2_,xg_H2O_,xg_CH4_,xg_S2_,xg_SO2_,xg_H2S_,xg_OCS_,Xg_t= "","","","","","","","","","",""
        wt_g,wt_g_O,wt_g_C,wt_g_H,wt_g_S  = "","","","",""
        wt_O_,wt_C_,wt_H_,wt_S_ = wt_O, wt_C, wt_H, wt_S
        
        if models.loc["print status","option"] == "yes":
            print(i,fO2_)
        
        # store results
        results2 = pd.DataFrame([[PT["P"],PT["T"],xg_O2_,xg_CO_,xg_CO2_,xg_H2_,xg_H2O_,xg_CH4_,xg_S2_,xg_SO2_,xg_H2S_,xg_OCS_,Xg_t,
               "","","","",
               wm_CO2_,wm_H2O_,wm_H2_,wm_CO_,wm_CH4_,wm_S_,wm_SO3_,wm_H2S_,wm_ST_,Fe32,Fe3T,S62,S6T,
               mg.fO22Dbuffer(PT,fO2_,"FMQ",models),mg.fO22Dbuffer(PT,fO2_,"NNO",models),SCSS_,sulfide_sat,SCAS_,sulfate_sat,
               wt_g,wt_g_O,wt_g_C,wt_g_H,wt_g_S,wt_O_,wt_C_,wt_H_,wt_S_,
fO2_,mg.f_H2(run,PT,melt_wf,setup,species,models),mg.f_H2O(run,PT,melt_wf,setup,species,models),mg.f_S2(run,PT,melt_wf,setup,species,models),mg.f_SO2(run,PT,melt_wf,setup,species,models),mg.f_SO3(run,PT,melt_wf,setup,species,models),mg.f_H2S(run,PT,melt_wf,setup,species,models),mg.f_CO2(run,PT,melt_wf,setup,species,models),mg.f_CO(run,PT,melt_wf,setup,species,models),mg.f_CH4(run,PT,melt_wf,setup,species,models),mg.f_OCS(run,PT,melt_wf,setup,species,models),
mg.y_O2(PT,species,models),mg.y_H2(PT,species,models),mg.y_H2O(PT,species,models),mg.y_S2(PT,species,models),mg.y_SO2(PT,species,models),mg.y_SO3(PT,species,models),mdv.y_H2S(PT,species,models),mg.y_CO2(PT,species,models),mg.y_CO(PT,species,models),mg.y_CH4(PT,species,models),mg.y_OCS(PT,species,models),
mg.M_m_SO(run,melt_wf,setup,species),mg.M_m_ox(run,melt_wf,setup,species,models),mg.C_H2O(run,PT,melt_wf,setup,species,models),mg.C_H2(run,PT,melt_wf,setup,species,models),mdv.C_CO3(run,PT,melt_wf,setup,species,models),mg.C_CO(run,PT,melt_wf,setup,species,models),mg.C_CH4(run,PT,melt_wf,setup,species,models),mg.C_S(run,PT,melt_wf,setup,species,models),mg.C_SO4(run,PT,melt_wf,setup,species,models),mg.C_H2S(run,PT,melt_wf,setup,species,models),
mg.KD1(run,PT,setup,species,models),mg.KHOg(PT,models),mg.KHOm(run,PT,melt_wf,setup,species,models),mg.KHOSg(PT,models),mg.KCOg(PT,models),mg.KCOHg(PT,models),mg.KOCSg(PT,models),mg.KOSg(PT,models),mg.KOSg2(PT,models)]])
        results_chemsitry = pd.concat([results_chemsitry, results2], ignore_index=True)
        if models.loc["output csv","option"] == "yes":
            results_chemistry.to_csv('results_cooling_chemistry.csv', index=False, header=False)

    return results

def titratingS(run,cooling_inputs,setup,species,models):
    if models.loc["print status","option"] == "yes":
        print(setup.loc[run,"Sample"])
    
    # check if any options need to be read from the setup file rather than the models file
    models = options_from_setup(run,models,setup)

    # set T, volatile composition of the melt, and tolerances
    PT={"P":setup.loc[run,"P_bar"]}
    PT["T"]=setup.loc[run,"T_C"]
    melt_wf = {'CO2':setup.loc[run,"CO2ppm"]/1000000.,"H2OT":setup.loc[run,"H2O"]/100.,"H2":0.}
    melt_wf["CT"] = (melt_wf["CO2"]/species.loc["CO2","M"])*species.loc["C","M"]
    melt_wf["HT"] = (2.*melt_wf["H2OT"]/species.loc["H2O","M"])*species.loc["H","M"]
    nr_step = cooling_inputs["nr_step"]
    nr_tol = cooling_inputs["nr_tol"]
    dt_step = cooling_inputs["dt_step"]
    psat_tol = cooling_inputs["psat_tol"]
        
    # Calculate S and Fe composition at initial T and check for sulfur saturation
    melt_wf["Fe3FeT"] = mg.Fe3FeT_i(run,PT,melt_wf,setup,species,models)
    Fe3Fe2_ = mg.Fe3Fe2(melt_wf)
    melt_wf["CO"] = 0.
    melt_wf["CH4"] = 0.
    melt_wf["H2S"] = 0. 
    wm_H2_, wm_CH4_, wm_H2S_, wm_CO_ = 0., 0., 0., 0.
    melt_wf["ST"] = 0.
    melt_wf["S2-"] = 0.
    melt_wf["S6+"] = 0.
    
    SCSS_,sulfide_sat,SCAS_,sulfate_sat,ST_ = sulfur_saturation(run,PT,melt_wf,setup,species,models)
    wm_S2m_ = melt_wf["S2-"]
    wm_S6p_ = melt_wf["S6+"]
    wm_SO3_ = (wm_S6p_/(species.loc["S","M"]))*(species.loc["S","M"] + 3.*species.loc["O","M"])
        
    # Set bulk composition
    wt_C, wt_O, wt_H, wt_S, wt_Fe, wt_g_, Wt_ = bulk_composition(run,PT,melt_wf,setup,species,models)
    
    melt_wf["ST"] = setup.loc[run,"STppm"]/1000000.           
    melt_wf["S6ST"] = mg.S6ST(run,PT,melt_wf,setup,species,models)
    melt_wf["S2-"] = (1.-melt_wf["S6ST"])*melt_wf["ST"]
    melt_wf["S6+"] = melt_wf["S6ST"]*melt_wf["ST"]     
    wt_S = melt_wf["ST"]
    
    bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"Fe":wt_Fe,"Wt":Wt_}
    
    P_sat_, wm_H2O_psat, wm_CO2_psat, wm_H2_psat, wm_CO_psat, wm_CH4_psat, wm_S2m_psat, wm_S6p_psat, wm_H2S_psat, H2O_HTpsat, H2_HTpsat, CH4_HTpsat, H2S_HTpsat, CO2_CTpsat, CO_CTpsat, CH4_CTpsat, S6p_STpsat, S2m_STpsat, H2S_STpsat = P_sat(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
    PT["P"] = setup.loc[run,"P_bar"]
    if models.loc["print status","option"] == "yes":
        print(P_sat_, PT["P"])
        
    # set system and initial guesses
    system = eq.set_system(melt_wf,models)

    if PT["P"] > P_sat_:
        guessx = mdv.f_O2(run,PT,melt_wf,setup,species,models)
        fO2_, A, B, C = eq.eq_SOFe_melt(run,PT,bulk_wf,melt_wf,species,setup,models,nr_step/10.,nr_tol,guessx)
        Fe32, Fe3T, S62, S6T = A
        melt_wf["S6ST"] = S6T
        melt_wf["Fe3FeT"] = Fe3T          
        melt_wf["S2-"] = (1.-melt_wf["S6ST"])*melt_wf["ST"]
        melt_wf["S6+"] = melt_wf["S6ST"]*melt_wf["ST"]
        guessx = mdv.f_O2(run,PT,melt_wf,setup,species,models)
    else:
        #PT["P"] = P_sat_
        #guessx, guessy, guessz = eq.initial_guesses(run,PT,melt_wf,setup,species,models,system)
        guessx, guessy, guessz = 3.28098558298437E-10,0.0485362525980284,0.
        #PT["P"] = setup.loc[run,"P_bar"]
        xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, Xm_t, Xm_t_ox, wm_H2O_, wm_H2_, wm_CO2_, wm_CO_, wm_CH4_, wm_S_, wm_SO3_, wm_H2S_, wm_ST_, Fe32, Fe3T, S62, S6T, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g, wt_O_, wt_C_, wt_H_, wt_S_, guessx, guessy, guessz, guessw = eq.mg_equilibrium(run,PT,melt_wf,bulk_wf,setup,species,models,nr_step,nr_tol,guessx,guessy,guessz,guessw)
        guessz = 0.00139401157503231
        
    # create results table
    results_header1 = pd.DataFrame([["System","run","Sample","H2O (mwf)","CO2 (mwf)","ST (mwf)","Fe3FeT","S6ST","O (twf)","C (twf)","H (wtf)","S (twf)","Fe (twf)",
"SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)",
"oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","ideal gas","carbonylsulfide","bulk composition","equilibrate Fe","starting pressure","gassing direction","gassing style","mass_volume","crystallisation","isotopes","Date"]])
    results_header2 = pd.DataFrame([[system,run,setup.loc[run,"Sample"],melt_wf["H2OT"],melt_wf["CO2"],melt_wf["ST"],melt_wf["Fe3FeT"],"SORT",bulk_wf["O"],bulk_wf["C"],bulk_wf["H"],bulk_wf["S"],bulk_wf["Fe"],
setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],
models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["bulk_composition","option"],models.loc["eq_Fe","option"],models.loc["starting_P","option"],models.loc["gassing_direction","option"],models.loc["gassing_style","option"],models.loc["mass_volume","option"],models.loc["crystallisation","option"],models.loc["isotopes","option"],date.today()]])
    results_header = pd.concat([results_header1, results_header2], ignore_index=True)
    results_chemistry1 = pd.DataFrame([["P","T('C)","xg_O2","xg_CO","xg_CO2","xg_H2","xg_H2O","xg_CH4","xg_S2","xg_SO2","xg_H2S","xg_OCS","Xg_t",
               "xm_CO2","xm_H2O","Xm_t_SO","Xm_t_ox",
               "wm_CO2","wm_H2O","wm_H2","wm_CO","wm_CH4","wm_S","wm_SO3","wm_H2S","wm_ST","Fe32","Fe3T","S62","S6T",
               "DFMQ","DNNO","SCSS","sulfide sat?","SCAS","sulfate sat?","wt_g","wt_g_O","wt_g_C","wt_g_H","wt_g_S","wt_O","wt_C","wt_H","wt_S",
               "fO2","fH2","fH2O","fS2","fSO2","fSO3","fH2S","fCO2","fCO","fCH4","fOCS",
               "yO2","yH2","yH2O","yS2","ySO2","ySO3","yH2S","yCO2","yCO","yCH4","yOCS",
               "M_m_SO","M_m_ox","C_H2O","C_H2","C_CO3","C_CO","C_CH4","C_S","C_SO4","C_H2S",
               "KD1","KHOg","KHOm","KHOSg","KCOg","KCOHg","KOCSg","KSOg","KSOg2"]])
    results_chemistry = pd.concat([results_header, results_chemistry1], ignore_index=True)
    results1 = pd.DataFrame([[PT["P"],PT["T"],mg.xg_O2(run,PT,melt_wf,setup,species,models),mg.xg_CO(run,PT,melt_wf,setup,species,models),mg.xg_CO2(run,PT,melt_wf,setup,species,models),mg.xg_H2(run,PT,melt_wf,setup,species,models),mg.xg_H2O(run,PT,melt_wf,setup,species,models),mg.xg_CH4(run,PT,melt_wf,setup,species,models),mg.xg_S2(run,PT,melt_wf,setup,species,models),mg.xg_SO2(run,PT,melt_wf,setup,species,models),mg.xg_H2S(run,PT,melt_wf,setup,species,models),mg.xg_OCS(run,PT,melt_wf,setup,species,models),mg.Xg_tot(run,PT,melt_wf,setup,species,models),
mg.xm_CO2_so(run,melt_wf,setup,species),mg.xm_H2OT_so(run,melt_wf,setup,species),mg.Xm_t_so(run,melt_wf,setup,species),mg.Xm_t_ox(run,melt_wf,setup,species),
melt_wf["CO2"],melt_wf["H2OT"],wm_H2_,wm_CO_,wm_CH4_,wm_S2m_,wm_SO3_,wm_H2S_,melt_wf["ST"],mg.Fe3Fe2(melt_wf),melt_wf["Fe3FeT"],mg.S6S2(run,PT,melt_wf,setup,species,models),mg.S6ST(run,PT,melt_wf,setup,species,models),
mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"FMQ",models),mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"NNO",models),SCSS_,sulfide_sat,SCAS_,sulfate_sat,wt_g_,"","","","",
bulk_wf["O"],bulk_wf["C"],bulk_wf["H"],bulk_wf["S"],
mdv.f_O2(run,PT,melt_wf,setup,species,models),mg.f_H2(run,PT,melt_wf,setup,species,models),mg.f_H2O(run,PT,melt_wf,setup,species,models),mg.f_S2(run,PT,melt_wf,setup,species,models),mg.f_SO2(run,PT,melt_wf,setup,species,models),mg.f_SO3(run,PT,melt_wf,setup,species,models),mg.f_H2S(run,PT,melt_wf,setup,species,models),mg.f_CO2(run,PT,melt_wf,setup,species,models),mg.f_CO(run,PT,melt_wf,setup,species,models),mg.f_CH4(run,PT,melt_wf,setup,species,models),mg.f_OCS(run,PT,melt_wf,setup,species,models),
mg.y_O2(PT,species,models),mg.y_H2(PT,species,models),mg.y_H2O(PT,species,models),mg.y_S2(PT,species,models),mg.y_SO2(PT,species,models),mg.y_SO3(PT,species,models),mdv.y_H2S(PT,species,models),mg.y_CO2(PT,species,models),mg.y_CO(PT,species,models),mg.y_CH4(PT,species,models),mg.y_OCS(PT,species,models),
mg.M_m_SO(run,melt_wf,setup,species),mg.M_m_ox(run,melt_wf,setup,species,models),mg.C_H2O(run,PT,melt_wf,setup,species,models),mg.C_H2(run,PT,melt_wf,setup,species,models),mdv.C_CO3(run,PT,melt_wf,setup,species,models),mg.C_CO(run,PT,melt_wf,setup,species,models),mg.C_CH4(run,PT,melt_wf,setup,species,models),mg.C_S(run,PT,melt_wf,setup,species,models),mg.C_SO4(run,PT,melt_wf,setup,species,models),mg.C_H2S(run,PT,melt_wf,setup,species,models),
mg.KD1(run,PT,setup,species,models),mg.KHOg(PT,models),mg.KHOm(run,PT,melt_wf,setup,species,models),mg.KHOSg(PT,models),mg.KCOg(PT,models),mg.KCOHg(PT,models),mg.KOCSg(PT,models),mg.KOSg(PT,models),mg.KOSg2(PT,models)]])
    results_chemistry = pd.concat([results_chemistry, results1], ignore_index=True)
    if models.loc["output csv","option"] == "yes":
        results_chemistry.to_csv('results_titrating_chemistry.csv', index=False, header=False)
    
    initial = 100
    final = 5000
    step = 100
    
    # run over different S concentrations #
    for i in range(initial,final,step): # P is pressure in bars or T is temperature in 'C
        eq_Fe = models.loc["eq_Fe","option"]

        melt_wf["ST"]=i/1000000.
        bulk_wf["S"]=i/1000000.
        
        system = eq.set_system(melt_wf,models)
        P_sat_, wm_H2O_psat, wm_CO2_psat, wm_H2_psat, wm_CO_psat, wm_CH4_psat, wm_S2m_psat, wm_S6p_psat, wm_H2S_psat, H2O_HTpsat, H2_HTpsat, CH4_HTpsat, H2S_HTpsat, CO2_CTpsat, CO_CTpsat, CH4_CTpsat, S6p_STpsat, S2m_STpsat, H2S_STpsat = P_sat(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
        PT["P"]=setup.loc[run,"P_bar"]
        
        # work out equilibrium partitioning between melt and gas phase
        if PT["P"] > P_sat_:
            fO2_, A, B, C = eq.eq_SOFe_melt(run,PT,bulk_wf,melt_wf,species,setup,models,nr_step,nr_tol,guessx)
            Fe32, Fe3T, S62, S6T = A
            wm_CO2_ = melt_wf["CO2"]
            wm_H2O_ = melt_wf["H2OT"]
            wm_H2_, wm_CO_, wm_CH4_, wm_H2S_ = "","","","" 
            melt_wf["Fe3FeT"] = Fe3T
            melt_wf["S6ST"] = S6T
            wm_ST_ = melt_wf["ST"]
            wm_S = wm_ST_
            melt_wf["S2-"] = wt_S*(1.-S6T) 
            wm_S_ = wt_S*(1.-S6T)
            wm_SO3_ = ((wt_S*S6T)/species.loc["S","M"])*species.loc["SO3","M"]
            xg_O2_,xg_CO_,xg_CO2_,xg_H2_,xg_H2O_,xg_CH4_,xg_S2_,xg_SO2_,xg_H2S_,xg_OCS_,Xg_t= "","","","","","","","","","",""
            wt_g,wt_g_O,wt_g_C,wt_g_H,wt_g_S  = 0.,"","","",""
            wt_O_,wt_C_,wt_H_,wt_S_ = wt_O, wt_C, wt_H, wt_S
        else:
            xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, Xm_t, Xm_t_ox, wm_H2O_, wm_H2_, wm_CO2_, wm_CO_, wm_CH4_, wm_S_, wm_SO3_, wm_H2S_, wm_ST_, Fe32, Fe3T, S62, S6T, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g, wt_O_, wt_C_, wt_H_, wt_S_, guessx, guessy, guessz, guessw = eq.mg_equilibrium(run,PT,melt_wf,bulk_wf,setup,species,models,nr_step,nr_tol,guessx,guessy,guessz,guessw)
            melt_wf["Fe3FeT"] = Fe3T
        
        # check for sulfur saturation and display warning in outputs
        SCSS_,sulfide_sat,SCAS_,sulfate_sat, ST_ = sulfur_saturation(run,PT,melt_wf,setup,species,models)
        if sulfide_sat == "yes":
            warning = "WARNING: sulfide-saturated"
        elif sulfate_sat == "yes":
            warning = "WARNING: sulfate-saturated"
        else:
            warning = ""
        
        # calculate fO2
        if eq_Fe == "yes":
            fO2_ = mdv.f_O2(run,PT,melt_wf,setup,species,models)
        elif eq_Fe == "no":
            fO2_ = (xg_O2_*mg.y_O2(PT,models)*PT["P"])
            
        if PT["P"] > P_sat_:
            guessx = fO2_
            
        if models.loc["print status","option"] == "yes":
            print(i,fO2_)
        
        # store results
        results2 = pd.DataFrame([[PT["P"],PT["T"],xg_O2_,xg_CO_,xg_CO2_,xg_H2_,xg_H2O_,xg_CH4_,xg_S2_,xg_SO2_,xg_H2S_,xg_OCS_,Xg_t,
               "","","","",
               wm_CO2_,wm_H2O_,wm_H2_,wm_CO_,wm_CH4_,wm_S_,wm_SO3_,wm_H2S_,wm_ST_,Fe32,Fe3T,S62,S6T,
               mg.fO22Dbuffer(PT,fO2_,"FMQ",models),mg.fO22Dbuffer(PT,fO2_,"NNO",models),SCSS_,sulfide_sat,SCAS_,sulfate_sat,
               wt_g,wt_g_O,wt_g_C,wt_g_H,wt_g_S,wt_O_,wt_C_,wt_H_,wt_S_,
fO2_,mg.f_H2(run,PT,melt_wf,setup,species,models),mg.f_H2O(run,PT,melt_wf,setup,species,models),mg.f_S2(run,PT,melt_wf,setup,species,models),mg.f_SO2(run,PT,melt_wf,setup,species,models),mg.f_SO3(run,PT,melt_wf,setup,species,models),mg.f_H2S(run,PT,melt_wf,setup,species,models),mg.f_CO2(run,PT,melt_wf,setup,species,models),mg.f_CO(run,PT,melt_wf,setup,species,models),mg.f_CH4(run,PT,melt_wf,setup,species,models),mg.f_OCS(run,PT,melt_wf,setup,species,models),
mg.y_O2(PT,species,models),mg.y_H2(PT,species,models),mg.y_H2O(PT,species,models),mg.y_S2(PT,species,models),mg.y_SO2(PT,species,models),mg.y_SO3(PT,species,models),mdv.y_H2S(PT,species,models),mg.y_CO2(PT,species,models),mg.y_CO(PT,species,models),mg.y_CH4(PT,species,models),mg.y_OCS(PT,species,models),
mg.M_m_SO(run,melt_wf,setup,species),mg.M_m_ox(run,melt_wf,setup,species,models),mg.C_H2O(run,PT,melt_wf,setup,species,models),mg.C_H2(run,PT,melt_wf,setup,species,models),mdv.C_CO3(run,PT,melt_wf,setup,species,models),mg.C_CO(run,PT,melt_wf,setup,species,models),mg.C_CH4(run,PT,melt_wf,setup,species,models),mg.C_S(run,PT,melt_wf,setup,species,models),mg.C_SO4(run,PT,melt_wf,setup,species,models),mg.C_H2S(run,PT,melt_wf,setup,species,models),
mg.KD1(run,PT,setup,species,models),mg.KHOg(PT,models),mg.KHOm(run,PT,melt_wf,setup,species,models),mg.KHOSg(PT,models),mg.KCOg(PT,models),mg.KCOHg(PT,models),mg.KOCSg(PT,models),mg.KOSg(PT,models),mg.KOSg2(PT,models)]])
        results_chemistry = pd.concat([results_chemistry, results2], ignore_index=True)
        if models.loc["output csv","option"] == "yes":
            results_chemistry.to_csv('results_titrating_chemistry.csv', index=False, header=False)

    return results_chemistry     
        
##########################
### X% sulfur degassed ###
##########################

def Xpc_S_degassed(first_row,last_row,inputs,setup,species,models):
     
    if models.loc["fO2","option"] != "Kress91A":
        raise TypeError("Kress91A must be used for fO2 option to be used")
    if models.loc["P_variation","option"] == "isobaric":
        raise TypeError("P_variation must be polybaric")

    nr_step = inputs["nr_step"]
    nr_tol = inputs["nr_tol"]
    psat_tol = inputs["psat_tol"]
    dp_step = inputs["dp_step"]
    
    # create results table
    results_header1 = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","ideal gas","carbonylsulfide","bulk composition","equilibrate Fe","starting pressure","gassing direction","gassing style","mass_volume","crystallisation","isotopes","Date"]])
    results_header2 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["bulk_composition","option"],models.loc["eq_Fe","option"],models.loc["starting_P","option"],models.loc["gassing_direction","option"],models.loc["gassing_style","option"],models.loc["mass_volume","option"],models.loc["crystallisation","option"],models.loc["isotopes","option"],date.today()]])
    results_header = pd.concat([results_header1, results_header2], ignore_index=True)
    results_chemistry1 = pd.DataFrame([["run","Sample","% S degassed","initial fO2 (DFMQ)","H2O (wt%)","CO2 (ppm)","ST (ppm)","X (ppm)","Fe3FeT","S6ST","O (twf)","C (twf)","H (wtf)","S (twf)","Fe (twf)","Saturation P (bars)",
"SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)","P","T('C)","xg_O2","xg_CO","xg_CO2","xg_H2","xg_H2O","xg_CH4","xg_S2","xg_SO2","xg_H2S","xg_OCS","xg_X","Xg_t","xg_CS",
               "xm_CO2","xm_H2O","Xm_t_SO","Xm_t_ox",
               "wm_CO2-eq","wm_H2O-eq","wm_CO2","wm_H2O","wm_H2","wm_CO","wm_CH4","wm_S","wm_SO3","wm_H2S","wm_ST","wm_X","Fe32","Fe3T","S62","S6T",
               "DFMQ","DNNO","wt_g",
               "fO2"]])
    results_chemistry = pd.concat([results_header, results_chemistry1], ignore_index=True)
    
    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        # set T, volatile composition of the melt, and tolerances
        PT={"T":setup.loc[run,"T_C"]}
        melt_wf = {'CO2':setup.loc[run,"CO2ppm"]/1000000.,"H2OT":setup.loc[run,"H2O"]/100.,"ST":setup.loc[run,"STppm"]/1000000.,"H2":0.,"XT":setup.loc[run,"Xppm"]/1000000.}
        melt_wf["CT"] = (melt_wf["CO2"]/species.loc["CO2","M"])*species.loc["C","M"]
        melt_wf["HT"] = (2.*melt_wf["H2OT"]/species.loc["H2O","M"])*species.loc["H","M"]
        melt_wf["ST"] = melt_wf["ST"]
        target_S = melt_wf["ST"]*((100.-setup.loc[run,"Xpc"])/100.)
        if models.loc["print status","option"] == "yes":
            print(run,setup.loc[run,"Sample"],melt_wf["ST"]*1000000.,target_S*1000000.,datetime.datetime.now())
        
        # Calculate saturation pressure
        P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST = c.P_sat(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
        wm_SO3_ = (wm_S6p_*species.loc["SO3","M"])/species.loc["S","M"]
        PT["P"] = P_sat_

        # update melt composition at saturation pressure and check for sulfur saturation
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
    
        # Set bulk composition
        wt_C, wt_O, wt_H, wt_S, wt_X, wt_Fe, wt_g_, Wt_ = c.bulk_composition(run,PT,melt_wf,setup,species,models)
        bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"Fe":wt_Fe,"Wt":Wt_,"X":wt_X}
    
        # set system and initial guesses
        system = eq.set_system(melt_wf,models)
        guessx, guessy, guessz, guessw = eq.initial_guesses(run,PT,melt_wf,setup,species,models,system)
            
        
        if models.loc["P_variation","option"] == "polybaric":
            # pressure ranges and options
            starting_P = models.loc["starting_P","option"]
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
    
        # add some gas to the system if doing open-system regassing
        if models.loc["gassing_direction","option"] == "regas" and models.loc["gassing_style","option"] == "open":
            gas_mf = {"O2":mg.xg_O2(run,PT,melt_wf,setup,species,models),"CO":mg.xg_CO(run,PT,melt_wf,setup,species,models),"CO2":mg.xg_CO2(run,PT,melt_wf,setup,species,models),"H2":mg.xg_H2(run,PT,melt_wf,setup,species,models),"H2O":mg.xg_H2O(run,PT,melt_wf,setup,species,models),"CH4":mg.xg_CH4(run,PT,melt_wf,setup,species,models),"S2":mg.xg_S2(run,PT,melt_wf,setup,species,models),"SO2":mg.xg_SO2(run,PT,melt_wf,setup,species,models),"SO3":mg.xg_SO3(run,PT,melt_wf,setup,species,models),"H2S":mg.xg_H2S(run,PT,melt_wf,setup,species,models),"OCS":mg.xg_OCS(run,PT,melt_wf,setup,species,models),"X":mg.xg_X(run,PT,melt_wf,setup,species,models),"Xg_t":mg.Xg_tot(run,PT,melt_wf,setup,species,models),"wt_g":0.}
            wt_C, wt_H, wt_S, wt_X, wt_Fe, wt_O, Wt = new_bulk_regas_open(run,PT,melt_wf,bulk_wf,gas_mf,dwtg,setup,species,models)
            bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"Fe":wt_Fe,"X":wt_X,"Wt":Wt}
    
        # run over different pressures #
        number_of_step = 0.
        
        P = initial
        wm_ST_ = melt_wf["ST"]
        while wm_ST_ > target_S: # step size = initial
            eq_Fe = models.loc["eq_Fe","option"]
            if models.loc["gassing_style","option"] == "open": # check melt is still vapor-saturated
                P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST = c.P_sat(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
                wm_SO3_ = (wm_S6p_*species.loc["SO3","M"])/species.loc["S","M"]
            P_ = P
            wm_ST__ = wm_ST_
            P = P - dp_step
            PT["P"] = P
            if P_sat_ > PT["P"]:  
                # work out equilibrium partitioning between melt and gas phase
                xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, xg_X_, Xg_t, xm_H2O_, xm_CO2_, Xm_t, wm_H2O_, wm_H2_, wm_CO2_, wm_CO_, wm_CH4_, wm_S_, wm_SO3_, wm_H2S_, wm_ST_, wm_X_, Fe32, Fe3T, S62, S6T, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g_X, wt_g, wt_O_, wt_C_, wt_H_, wt_S_, wt_X_, guessx, guessy, guessz, guessw = eq.mg_equilibrium(run,PT,melt_wf,bulk_wf,setup,species,models,nr_step,nr_tol,guessx,guessy,guessz,guessw)
                # gas composition
                gas_mf = {"O2":xg_O2_,"CO":xg_CO_,"S2":xg_S2_,"CO2":xg_CO2_,"H2O":xg_H2O_,"H2":xg_H2_,"CH4":xg_CH4_,"SO2":xg_SO2_,"H2S":xg_H2S_,"OCS":xg_OCS_,"X":xg_X_,"Xg_t":Xg_t,"wt_g":wt_g}
            
            else:
                conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
                frac = c.melt_species_ratios(conc,species)
                wm_ST_ = wm_S_ + wm_S6p_
                S62 = S6T/S2m_ST
                Fe3T = melt_wf["Fe3FeT"]
                Fe32 = mg.overtotal2ratio(Fe3T)
                xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, xg_X, Xg_t, Xm_t, Xm_t_ox, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g = "","","","","","","","","","","","","","",0.,0.,0.,0.,0.
                guessx, guessy, guessz, guessw = eq.initial_guesses(run,PT,melt_wf,setup,species,models,system)
        
        P = P_
        wm_ST_ = wm_ST__
        while wm_ST_ > target_S: # step size = 1.
            eq_Fe = models.loc["eq_Fe","option"]
            if models.loc["gassing_style","option"] == "open": # check melt is still vapor-saturated
                P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST = c.P_sat(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
                wm_SO3_ = (wm_S6p_*species.loc["SO3","M"])/species.loc["S","M"]
            P_ = P
            P = P - 1.
            PT["P"] = P
            if P_sat_ > PT["P"]:  
                # work out equilibrium partitioning between melt and gas phase
                xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, xg_X_, Xg_t, xm_H2O_, xm_CO2_, Xm_t, wm_H2O_, wm_H2_, wm_CO2_, wm_CO_, wm_CH4_, wm_S_, wm_SO3_, wm_H2S_, wm_ST_, wm_X_, Fe32, Fe3T, S62, S6T, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g_X, wt_g, wt_O_, wt_C_, wt_H_, wt_S_, wt_X_, guessx, guessy, guessz, guessw = eq.mg_equilibrium(run,PT,melt_wf,bulk_wf,setup,species,models,nr_step,nr_tol,guessx,guessy,guessz,guessw)
                # gas composition
                gas_mf = {"O2":xg_O2_,"CO":xg_CO_,"S2":xg_S2_,"CO2":xg_CO2_,"H2O":xg_H2O_,"H2":xg_H2_,"CH4":xg_CH4_,"SO2":xg_SO2_,"H2S":xg_H2S_,"OCS":xg_OCS_,"X":xg_X_,"Xg_t":Xg_t,"wt_g":wt_g}
            
            else:
                conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
                frac = c.melt_species_ratios(conc,species)
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
    
        # check for sulfur saturation and display warning in outputs
        SCSS_,sulfide_sat,SCAS_,sulfate_sat, ST_ = c.sulfur_saturation(run,PT,melt_wf,setup,species,models)
        if sulfide_sat == "yes":
            warning = "WARNING: sulfide-saturated"
        elif sulfate_sat == "yes":
            warning = "WARNING: sulfate-saturated"
        else:
            warning = ""
        
        # calculate fO2
        if eq_Fe == "yes":
            fO2_ = mdv.f_O2(run,PT,melt_wf,setup,species,models)
        elif eq_Fe == "no":
            fO2_ = (xg_O2_*mdv.y_O2(PT,models)*PT["P"])
        
        xg_CS = mg.gas_CS(run,PT,melt_wf,setup,species,models)
        wm_CO2eq, wm_H2Oeq = mg.melt_H2O_CO2_eq(melt_wf,species)
        
        # store results
        results1 = pd.DataFrame([[run,setup.loc[run,"Sample"],setup.loc[run,"Xpc"],setup.loc[run,"DFMQ"],setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"STppm"],setup.loc[run,"Xppm"],melt_wf["Fe3FeT"],"SORT",bulk_wf["O"],bulk_wf["C"],bulk_wf["H"],bulk_wf["S"],bulk_wf["Fe"],P_sat_,
setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],PT["P"],PT["T"],mg.xg_O2(run,PT,melt_wf,setup,species,models),mg.xg_CO(run,PT,melt_wf,setup,species,models),mg.xg_CO2(run,PT,melt_wf,setup,species,models),mg.xg_H2(run,PT,melt_wf,setup,species,models),mg.xg_H2O(run,PT,melt_wf,setup,species,models),mg.xg_CH4(run,PT,melt_wf,setup,species,models),mg.xg_S2(run,PT,melt_wf,setup,species,models),mg.xg_SO2(run,PT,melt_wf,setup,species,models),mg.xg_H2S(run,PT,melt_wf,setup,species,models),mg.xg_OCS(run,PT,melt_wf,setup,species,models),mg.xg_X(run,PT,melt_wf,setup,species,models),mg.Xg_tot(run,PT,melt_wf,setup,species,models),xg_CS,
mg.xm_CO2_so(run,melt_wf,setup,species),mg.xm_H2OT_so(run,melt_wf,setup,species),mg.Xm_t_so(run,melt_wf,setup,species),"",
wm_CO2eq,wm_H2Oeq,melt_wf["CO2"],melt_wf["H2OT"],wm_H2_,wm_CO_,wm_CH4_,wm_S2m_,wm_SO3_,wm_H2S_,melt_wf["ST"],melt_wf["XT"],mg.Fe3Fe2(melt_wf),melt_wf["Fe3FeT"],mg.S6S2(run,PT,melt_wf,setup,species,models),mg.S6ST(run,PT,melt_wf,setup,species,models),
mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"FMQ",models),mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"NNO",models),wt_g,
mdv.f_O2(run,PT,melt_wf,setup,species,models)]])
        results_chemistry = pd.concat([results_chemistry, results1], ignore_index=True)
        if models.loc["output csv","option"] == "yes":
            results_chemistry.to_csv('results_XpcSdegasses.csv', index=False, header=False)

    if models.loc["print status","option"] == "yes":            
        print("done", datetime.datetime.now())

    return results  
    
####################
### fO2 at 1 bar ###
####################

def fO2_at_1bar(first_row,last_row,inputs,setup,species,models):
     
    if models.loc["fO2","option"] != "Kress91A":
        raise TypeError("Kress91A must be used for fO2 option")
    if models.loc["P_variation","option"] == "isobaric":
        raise TypeError("P_variation must be polybaric")

    nr_step = inputs["nr_step"]
    nr_tol = inputs["nr_tol"]
    psat_tol = inputs["psat_tol"]
    P1 = inputs["P1"]
    P2 = inputs["P2"]
    dp_step = inputs["dp_step"]
    
    # create results table
    results_header1 = pd.DataFrame([["oxygen fugacity","carbon dioxide solubility","C speciation composition","water solubility","water speciation","water speciation composition","sulfide solubility","sulfate solubility","ideal gas","carbonylsulfide","bulk composition","equilibrate Fe","starting pressure","gassing direction","gassing style","mass_volume","crystallisation","isotopes","Date"]])
    results_header2 = pd.DataFrame([[models.loc["fO2","option"],models.loc["carbon dioxide","option"],models.loc["Cspeccomp","option"],models.loc["water","option"],models.loc["Hspeciation","option"],models.loc["Hspeccomp","option"],models.loc["sulfide","option"],models.loc["sulfate","option"],models.loc["ideal_gas","option"],models.loc["carbonylsulfide","option"],models.loc["bulk_composition","option"],models.loc["eq_Fe","option"],models.loc["starting_P","option"],models.loc["gassing_direction","option"],models.loc["gassing_style","option"],models.loc["mass_volume","option"],models.loc["crystallisation","option"],models.loc["isotopes","option"],date.today()]])
    results_header = pd.concat([results_header1, results_header2], ignore_index=True)
    results_chemistry1 = pd.DataFrame([["run","Sample","initial fO2 (DFMQ)","H2O (wt%)","CO2 (ppm)","ST (ppm)","X (ppm)","Fe3FeT","S6ST","O (twf)","C (twf)","H (wtf)","S (twf)","Fe (twf)","Saturation P (bars)",
"SiO2 (wt%)","TiO2 (wt%)","Al2O3 (wt%)","FeOT (wt%)","MnO (wt%)","MgO (wt%)","CaO (wt%)","Na2O (wt%)","K2O (wt%)","P2O5 (wt%)","P","T('C)","xg_O2","xg_CO","xg_CO2","xg_H2","xg_H2O","xg_CH4","xg_S2","xg_SO2","xg_H2S","xg_OCS","xg_X","Xg_t","xg_CS",
               "xm_CO2","xm_H2O","Xm_t_SO","Xm_t_ox",
               "wm_CO2-eq","wm_H2O-eq","wm_CO2","wm_H2O","wm_H2","wm_CO","wm_CH4","wm_S","wm_SO3","wm_H2S","wm_ST","wm_X","Fe32","Fe3T","S62","S6T",
               "DFMQ","DNNO","wt_g",
               "fO2","D-DFMQ"]])
    results_chemistry1 = pd.concat([results_header, results_chemistry1], ignore_index=True)
    results_chemistry_P1 = pd.concat([results_header, results_chemistry1], ignore_index=True) 
    results_chemistry_P2 = pd.concat([results_header, results_chemistry1], ignore_index=True)
    
    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        # set T, volatile composition of the melt, and tolerances
        PT={"T":setup.loc[run,"T_C"]}
        melt_wf = {'CO2':setup.loc[run,"CO2ppm"]/1000000.,"H2OT":setup.loc[run,"H2O"]/100.,"ST":setup.loc[run,"STppm"]/1000000.,"H2":0.,"XT":setup.loc[run,"Xppm"]/1000000.}
        melt_wf["CT"] = (melt_wf["CO2"]/species.loc["CO2","M"])*species.loc["C","M"]
        melt_wf["HT"] = (2.*melt_wf["H2OT"]/species.loc["H2O","M"])*species.loc["H","M"]
        melt_wf["ST"] = melt_wf["ST"]
        target_S = melt_wf["ST"]*((100.-setup.loc[run,"Xpc"])/100.)
        if models.loc["print status","option"] == "yes":
            print(run,setup.loc[run,"Sample"],melt_wf["ST"]*1000000.,target_S*1000000.,datetime.datetime.now())
        
        # Calculate saturation pressure
        P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST = c.P_sat(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
        wm_SO3_ = (wm_S6p_*species.loc["SO3","M"])/species.loc["S","M"]
        PT["P"] = P_sat_

        # update melt composition at saturation pressure and check for sulfur saturation
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
    
        # Set bulk composition
        wt_C, wt_O, wt_H, wt_S, wt_X, wt_Fe, wt_g_, Wt_ = c.bulk_composition(run,PT,melt_wf,setup,species,models)
        bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"Fe":wt_Fe,"Wt":Wt_,"X":wt_X}
    
        # set system and initial guesses
        system = eq.set_system(melt_wf,models)
        guessx, guessy, guessz, guessw = eq.initial_guesses(run,PT,melt_wf,setup,species,models,system)
            
        
        if models.loc["P_variation","option"] == "polybaric":
            # pressure ranges and options
            starting_P = models.loc["starting_P","option"]
            if starting_P == "set":
                initial = int(setup.loc[run,"P_bar"])
            else:
                if models.loc["gassing_direction","option"] == "degas":
                    answer = PT["P"]/dp_step
                    answer = math.floor(answer)
                    initial = round(answer*dp_step)
                elif models.loc["gassing_direction","option"] == "regas":
                    initial = round(PT["P"])+1 
            if models.loc["gassing_direction","option"] == "degas":
                step = int(-1*dp_step) # pressure step in bars
                final = 0
            elif models.loc["gassing_direction","option"] == "regas":
                step = int(dp_step)
                final = int(setup.loc[run,"final_P"])
    
        # add some gas to the system if doing open-system regassing
        if models.loc["gassing_direction","option"] == "regas" and models.loc["gassing_style","option"] == "open":
            gas_mf = {"O2":mg.xg_O2(run,PT,melt_wf,setup,species,models),"CO":mg.xg_CO(run,PT,melt_wf,setup,species,models),"CO2":mg.xg_CO2(run,PT,melt_wf,setup,species,models),"H2":mg.xg_H2(run,PT,melt_wf,setup,species,models),"H2O":mg.xg_H2O(run,PT,melt_wf,setup,species,models),"CH4":mg.xg_CH4(run,PT,melt_wf,setup,species,models),"S2":mg.xg_S2(run,PT,melt_wf,setup,species,models),"SO2":mg.xg_SO2(run,PT,melt_wf,setup,species,models),"SO3":mg.xg_SO3(run,PT,melt_wf,setup,species,models),"H2S":mg.xg_H2S(run,PT,melt_wf,setup,species,models),"OCS":mg.xg_OCS(run,PT,melt_wf,setup,species,models),"X":mg.xg_X(run,PT,melt_wf,setup,species,models),"Xg_t":mg.Xg_tot(run,PT,melt_wf,setup,species,models),"wt_g":0.}
            wt_C, wt_H, wt_S, wt_X, wt_Fe, wt_O, Wt = new_bulk_regas_open(run,PT,melt_wf,bulk_wf,gas_mf,dwtg,setup,species,models)
            bulk_wf = {"C":wt_C,"H":wt_H,"O":wt_O,"S":wt_S,"Fe":wt_Fe,"X":wt_X,"Wt":Wt}
    
        # run over different pressures #
        number_of_step = 0.
        
        P = initial
        wm_ST_ = melt_wf["ST"]
        for i in range(initial,final,step):
            eq_Fe = models.loc["eq_Fe","option"]
            if models.loc["gassing_style","option"] == "open": # check melt is still vapor-saturated
                P_sat_, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_, H2O_HT, H2_HT, CH4_HT, H2S_HT, CO2_CT, CO_CT, CH4_CT, S6p_ST, S2m_ST, H2S_ST = c.P_sat(run,PT,melt_wf,setup,species,models,psat_tol,nr_step,nr_tol)
                wm_SO3_ = (wm_S6p_*species.loc["SO3","M"])/species.loc["S","M"]
            P_ = P
            wm_ST__ = wm_ST_
            P = P - dp_step
            if P <= 0.:
                P = 1
            PT["P"] = P
            if P_sat_ > PT["P"]:  
                # work out equilibrium partitioning between melt and gas phase
                xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, xg_X_, Xg_t, xm_H2O_, xm_CO2_, Xm_t, wm_H2O_, wm_H2_, wm_CO2_, wm_CO_, wm_CH4_, wm_S_, wm_SO3_, wm_H2S_, wm_ST_, wm_X_, Fe32, Fe3T, S62, S6T, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g_X, wt_g, wt_O_, wt_C_, wt_H_, wt_S_, wt_X_, guessx, guessy, guessz, guessw = eq.mg_equilibrium(run,PT,melt_wf,bulk_wf,setup,species,models,nr_step,nr_tol,guessx,guessy,guessz,guessw)
                # gas composition
                gas_mf = {"O2":xg_O2_,"CO":xg_CO_,"S2":xg_S2_,"CO2":xg_CO2_,"H2O":xg_H2O_,"H2":xg_H2_,"CH4":xg_CH4_,"SO2":xg_SO2_,"H2S":xg_H2S_,"OCS":xg_OCS_,"X":xg_X_,"Xg_t":Xg_t,"wt_g":wt_g}
            
            else:
                conc = eq.melt_speciation(PT,melt_wf,species,models,nr_step,nr_tol)
                frac = c.melt_species_ratios(conc,species)
                wm_ST_ = wm_S_ + wm_S6p_
                S62 = S6T/S2m_ST
                Fe3T = melt_wf["Fe3FeT"]
                Fe32 = mg.overtotal2ratio(Fe3T)
                xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, xg_X, Xg_t, Xm_t, Xm_t_ox, wt_g_O, wt_g_C, wt_g_H, wt_g_S, wt_g = "","","","","","","","","","","","","","",0.,0.,0.,0.,0.
                guessx, guessy, guessz, guessw = eq.initial_guesses(run,PT,melt_wf,setup,species,models,system)
        
            if P == 1 or P == P1 or P == P2:
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
    
                # check for sulfur saturation and display warning in outputs
                SCSS_,sulfide_sat,SCAS_,sulfate_sat, ST_ = c.sulfur_saturation(run,PT,melt_wf,setup,species,models)
                if sulfide_sat == "yes":
                    warning = "WARNING: sulfide-saturated"
                elif sulfate_sat == "yes":
                    warning = "WARNING: sulfate-saturated"
                else:
                    warning = ""
        
                # calculate fO2
                if eq_Fe == "yes":
                    fO2_ = mdv.f_O2(run,PT,melt_wf,setup,species,models)
                elif eq_Fe == "no":
                    fO2_ = (xg_O2_*mdv.y_O2(PT,models)*PT["P"])
        
                xg_CS = mg.gas_CS(run,PT,melt_wf,setup,species,models)
                wm_CO2eq, wm_H2Oeq = mg.melt_H2O_CO2_eq(melt_wf,species)                
    
                # store results
                results1 = pd.DataFrame([[run,setup.loc[run,"Sample"],setup.loc[run,"DFMQ"],setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"STppm"],setup.loc[run,"Xppm"],melt_wf["Fe3FeT"],"SORT",bulk_wf["O"],bulk_wf["C"],bulk_wf["H"],bulk_wf["S"],bulk_wf["Fe"],P_sat_,
setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],mg.Wm_FeOT(run,setup,species),setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],PT["P"],PT["T"],mg.xg_O2(run,PT,melt_wf,setup,species,models),mg.xg_CO(run,PT,melt_wf,setup,species,models),mg.xg_CO2(run,PT,melt_wf,setup,species,models),mg.xg_H2(run,PT,melt_wf,setup,species,models),mg.xg_H2O(run,PT,melt_wf,setup,species,models),mg.xg_CH4(run,PT,melt_wf,setup,species,models),mg.xg_S2(run,PT,melt_wf,setup,species,models),mg.xg_SO2(run,PT,melt_wf,setup,species,models),mg.xg_H2S(run,PT,melt_wf,setup,species,models),mg.xg_OCS(run,PT,melt_wf,setup,species,models),mg.xg_X(run,PT,melt_wf,setup,species,models),mg.Xg_tot(run,PT,melt_wf,setup,species,models),xg_CS,
mg.xm_CO2_so(run,melt_wf,setup,species),mg.xm_H2OT_so(run,melt_wf,setup,species),mg.Xm_t_so(run,melt_wf,setup,species),"",
wm_CO2eq,wm_H2Oeq,melt_wf["CO2"],melt_wf["H2OT"],wm_H2_,wm_CO_,wm_CH4_,wm_S2m_,wm_SO3_,wm_H2S_,melt_wf["ST"],melt_wf["XT"],mg.Fe3Fe2(melt_wf),melt_wf["Fe3FeT"],mg.S6S2(run,PT,melt_wf,setup,species,models),mg.S6ST(run,PT,melt_wf,setup,species,models),
mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"FMQ",models),mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"NNO",models),wt_g,
mdv.f_O2(run,PT,melt_wf,setup,species,models),setup.loc[run,"DFMQ"]-mg.fO22Dbuffer(PT,mdv.f_O2(run,PT,melt_wf,setup,species,models),"FMQ",models)]])
                
                if P == P1:
                    results_chemistry_P1 = pd.concat([results_chemistry_P1, results1], ignore_index=True)
                    if models.loc["output csv","option"] == "yes":
                        results_chemistry_P1.to_csv('results_P1.csv', index=False, header=False)
                if P == P2:
                    results_chemistry_P1 = pd.concat([results_chemistry_P2, results1], ignore_index=True)
                    if models.loc["output csv","option"] == "yes":
                        results_chemistry_P2.to_csv('results_P2.csv', index=False, header=False)
                if P == 1:
                    results_chemistry_1 = pd.concat([results_chemistry_1, results1], ignore_index=True)
                    if models.loc["output csv","option"] == "yes":
                        results_chemistry_1.to_csv('results_1bar.csv', index=False, header=False)
    if models.loc["print status","option"] == "yes":            
        print("done", datetime.datetime.now()) 
    return results_chemistry_P1, results_chemistry_P2, results_chemistry_1       