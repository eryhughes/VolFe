# batch_calculations.py

import pandas as pd
from datetime import date
import numpy as np
import datetime
import math as math
import warnings as w

import VolFe.melt_gas as mg
import VolFe.equilibrium_equations as eq
import VolFe.isotopes as iso
import VolFe.model_dependent_variables as mdv
import VolFe.calculations as c

################
### Contents ###
################
# building results tables
# options from setup file
# calculate the pressure of vapor saturation
# calculate de/regassing paths
# calculate isobars
# calculate solubility constants
# calculate fugacity coefficients
# Use melt S oxybarometer
# measured parameters within error
# Below this: in development

###############################
### building results tables ###
###############################
# outputing sample name
def results_table_sample_name(setup,run):
    results_headers = pd.DataFrame([["sample"]])
    results_values = pd.DataFrame([[setup.loc[run,"Sample"]]])
    return results_headers, results_values
# outputting melt composition, T, P
def results_table_melt_comp_etc(PT,melt_comp,conc,frac,melt_wf):
    results_headers = pd.DataFrame([["T_C","P_bar",
        "SiO2_wtpc", "TiO2_wtpc", "Al2O3_wtpc", "FeOT_wtpc", "MnO_wtpc", "MgO_wtpc", "CaO_wtpc", "Na2O_wtpc", "K2O_wtpc", "P2O5_wtpc",
        "H2OT_wtpc","OH_wtpc","H2Omol_wtpc","H2_ppmw","CH4_ppmw","CO2T_ppmw","CO2mol_ppmw","CO2carb_ppmw","CO_ppmw","S2-_ppmw","S6+_ppmw","H2S_ppmw",
        "H_H2OT/HT", "H_H2/HT", "H_CH4/HT", "H_H2S/HT", "C_CO2T/CT", "C_CO/CT", "C_CH4/CT", "S2-/ST", "S6+/ST", "H2S/ST", "Fe3+/FeT","sulf_XFe","sulf_XCu","sulf_XNi"]])
    if "sulf_XFe" in melt_wf:
        melt_wf
    else:
        melt_wf["sulf_XFe"] = 1.
    if "sulf_XCu" in melt_wf:
        melt_wf
    else:
        melt_wf["sulf_XCu"] = 0.
    if "sulf_XNi" in melt_wf:
        melt_wf
    else:
        melt_wf["sulf_XNi"] = 0.
    results_values = pd.DataFrame([[PT["T"],PT["P"],
                melt_comp["SiO2"]*100., melt_comp["TiO2"]*100., melt_comp["Al2O3"]*100., melt_comp["FeOT"]*100., melt_comp["MnO"]*100., melt_comp["MgO"]*100., melt_comp["CaO"]*100., melt_comp["Na2O"]*100., melt_comp["K2O"]*100., melt_comp["P2O5"]*100.,
                conc["wm_H2O"]*100.,conc["wm_OH"]*100,conc["wm_H2Omol"]*100.,conc["wm_H2"]*1000000.,conc["wm_CH4"]*1000000.,conc["wm_CO2"]*1000000.,conc["wm_CO2mol"]*1000000,conc["wm_CO2carb"]*1000000,conc["wm_CO"]*1000000.,conc["wm_S2m"]*1000000.,conc["wm_S6p"]*1000000.,conc["wm_H2S"]*1000000.,
                frac["H2O_HT"], frac["H2_HT"], frac["CH4_HT"], frac["H2S_HT"], frac["CO2_CT"], frac["CO_CT"], frac["CH4_CT"], frac["S2m_ST"], frac["S6p_ST"], frac["H2S_ST"],melt_wf["Fe3FeT"],melt_wf["sulf_XFe"],melt_wf["sulf_XCu"],melt_wf["sulf_XNi"]]])
    return results_headers, results_values
def results_table_melt_vol():
    results_headers = pd.DataFrame([["H2OT-eq_wtpc","CO2T-eq_ppmw","ST_ppmw","X_ppmw"]])
    return results_headers
# outputting model options used in the calculation
def results_table_model_options(models): 
    results_headers = pd.DataFrame([["setup opt","COH_species opt","H2S_m opt","species X opt","Hspeciation opt",
                 "fO2 opt","NNObuffer opt","FMQbuffer opt",
                 "carbon dioxide opt","water opt","hydrogen opt","sulfide opt","sulfate opt","hydrogen sulfide opt","methane opt","carbon monoxide opt","species X solubility opt","Cspeccomp opt","Hspeccomp opt",
                 "SCSS opt","SCAS opt","sulfur_saturation opt","sulfur_is_sat opt","graphite_saturation opt","ideal_gas opt",
                 "y_CO2 opt","y_SO2 opt","y_H2S opt","y_H2 opt","y_O2 opt","y_S2 opt","y_CO opt","y_CH4 opt","y_H2O opt","y_OCS opt","y_X opt",
                 "KHOg opt","KHOSg opt","KOSg opt","KOSg2 opt","KCOg opt","KCOHg opt","KOCSg opt","KCOs opt","carbonylsulfide opt",
                 "density opt","Date"]])
    results_values = pd.DataFrame([[models.loc["setup","option"],models.loc["COH_species","option"], models.loc["H2S_m","option"], models.loc["species X","option"],models.loc["Hspeciation","option"], 
                models.loc["fO2","option"], models.loc["NNObuffer","option"], models.loc["FMQbuffer","option"],
                 models.loc["carbon dioxide","option"], models.loc["water","option"], models.loc["hydrogen","option"], models.loc["sulfide","option"], models.loc["sulfate","option"], models.loc["hydrogen sulfide","option"], models.loc["methane","option"], models.loc["carbon monoxide","option"], models.loc["species X solubility","option"], models.loc["Cspeccomp","option"], models.loc["Hspeccomp","option"],
                 models.loc["SCSS","option"], models.loc["SCAS","option"], models.loc["sulfur_saturation","option"], models.loc["sulfur_is_sat","option"], models.loc["graphite_saturation","option"], models.loc["ideal_gas","option"],
                 models.loc["y_CO2","option"], models.loc["y_SO2","option"], models.loc["y_H2S","option"], models.loc["y_H2","option"], models.loc["y_O2","option"], models.loc["y_S2","option"], models.loc["y_CO","option"], models.loc["y_CH4","option"], models.loc["y_H2O","option"],models.loc["y_OCS","option"], models.loc["y_X","option"],
                 models.loc["KHOg","option"], models.loc["KHOSg","option"], models.loc["KOSg","option"], models.loc["KOSg2","option"], models.loc["KCOg","option"], models.loc["KCOHg","option"],models.loc["KOCSg","option"], models.loc["KCOs","option"],models.loc["carbonylsulfide","option"],
                 models.loc["density","option"],datetime.datetime.now()]])
    return results_headers, results_values
# outputting fugacities, partial pressures, gas mole fraction, fugacity coefficients, molecular masses, solubility constants, equilibrium constants, melt density
def results_table_f_p_xg_y_M_C_K_d(PT,melt_wf,models): 
    results_headers = pd.DataFrame([["fO2_DNNO","fO2_DFMQ",
                "fO2_bar","fH2_bar","fH2O_bar","fS2_bar","fSO2_bar","fH2S_bar","fCO2_bar","fCO_bar","fCH4_bar","fOCS_bar","fX_bar",
                "pO2_bar","pH2_bar","pH2O_bar","pS2_bar","pSO2_bar","pH2S_bar","pCO2_bar","pCO_bar","pCH4_bar","pOCS_bar","pX_bar",
                "xgO2_mf","xgH2_mf","xgH2O_mf","xgS2_mf","xgSO2_mf","xgH2S_mf","xgCO2_mf","xgCO_mf","xgCH4_mf","xgOCS_mf","xgX_mf","xgC_S_mf",
                "yO2","yH2","yH2O","yS2","ySO2","yH2S","yCO2","yCO","yCH4","yOCS","yX",
                "M_m_SO","M_m_ox",
                "C_H2O_mf_bar","C_H2_ppm_bar","C_CO2T_mf_bar","C_CO_ppm_bar","C_CH4_ppm_bar","C_S_ppm","C_SO4_ppm_bar","C_H2S_ppm_bar","C_X_ppm_bar",
                "KHOg","KHOSg","KCOg","KCOHg","KOCSg","KSOg","KSOg2","KHOm","KCOm","KCOs",
                "density_gcm3"]])
    results_values = pd.DataFrame([[mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,models),"NNO",models),mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,models),"FMQ",models),
                mdv.f_O2(PT,melt_wf,models),mg.f_H2(PT,melt_wf,models),mg.f_H2O(PT,melt_wf,models),mg.f_S2(PT,melt_wf,models),mg.f_SO2(PT,melt_wf,models),mg.f_H2S(PT,melt_wf,models),mg.f_CO2(PT,melt_wf,models),mg.f_CO(PT,melt_wf,models),mg.f_CH4(PT,melt_wf,models),mg.f_OCS(PT,melt_wf,models),mg.f_X(PT,melt_wf,models),
                mg.p_O2(PT,melt_wf,models),mg.p_H2(PT,melt_wf,models),mg.p_H2O(PT,melt_wf,models),mg.p_S2(PT,melt_wf,models),mg.p_SO2(PT,melt_wf,models),mg.p_H2S(PT,melt_wf,models),mg.p_CO2(PT,melt_wf,models),mg.p_CO(PT,melt_wf,models),mg.p_CH4(PT,melt_wf,models),mg.p_OCS(PT,melt_wf,models),mg.p_X(PT,melt_wf,models),
                mg.xg_O2(PT,melt_wf,models),mg.xg_H2(PT,melt_wf,models),mg.xg_H2O(PT,melt_wf,models),mg.xg_S2(PT,melt_wf,models),mg.xg_SO2(PT,melt_wf,models),mg.xg_H2S(PT,melt_wf,models),mg.xg_CO2(PT,melt_wf,models),mg.xg_CO(PT,melt_wf,models),mg.xg_CH4(PT,melt_wf,models),mg.xg_OCS(PT,melt_wf,models),mg.xg_X(PT,melt_wf,models),mg.gas_CS(PT,melt_wf,models),
                mdv.y_O2(PT,models),mdv.y_H2(PT,models),mdv.y_H2O(PT,models),mdv.y_S2(PT,models),mdv.y_SO2(PT,models),mdv.y_H2S(PT,models),mdv.y_CO2(PT,models),mdv.y_CO(PT,models),mdv.y_CH4(PT,models),mdv.y_OCS(PT,models),mdv.y_X(PT,models),
                mg.M_m_SO(melt_wf),mg.M_m_ox(melt_wf,models),
                mdv.C_H2O(PT,melt_wf,models),mdv.C_H2(PT,melt_wf,models),mdv.C_CO3(PT,melt_wf,models),mdv.C_CO(PT,melt_wf,models),mdv.C_CH4(PT,melt_wf,models),mdv.C_S(PT,melt_wf,models),mdv.C_SO4(PT,melt_wf,models),mdv.C_H2S(PT,melt_wf,models),mdv.C_X(PT,melt_wf,models),
                mdv.KHOg(PT,models),mdv.KHOSg(PT,models),mdv.KCOg(PT,models),mdv.KCOHg(PT,models),mdv.KOCSg(PT,models),mdv.KOSg(PT,models),mdv.KOSg2(PT,models),mdv.KHOm(PT,melt_wf,models),mdv.KCOm(PT,melt_wf,models),mdv.KCOs(PT,models),
                mdv.melt_density(PT,melt_wf,models)]])
    return results_headers, results_values
# headers for open system degassing all gas
def results_table_open_all_gas():
    results_headers = pd.DataFrame([["xgO2_all_mf","xgH2_all_mf","xgH2O_all_mf","xgS2_all_mf","xgSO2_all_mf","xgH2S_all_mf","xgCO2_all_mf","xgCO_all_mf","xgCH4_all_mf","xgOCS_all_mf","xgX_all_mf","xgC_S_all_mf"]])
    return results_headers
# saturation conditions
def results_table_sat(sulf_sat_result,PT,melt_wf,models):
    results_headers = pd.DataFrame([["SCSS_ppm","sulfide saturated","SCAS_ppm","anhydrite saturated","ST melt if sat","graphite saturated"]])
    results_values = pd.DataFrame([[sulf_sat_result["SCSS"],sulf_sat_result["sulfide_sat"],sulf_sat_result["SCAS"],sulf_sat_result["sulfate_sat"],sulf_sat_result["ST"],c.graphite_saturation(PT,melt_wf,models)]])
    return results_headers, results_values
# isotopes
def results_table_isotope_R(R,R_all_species_S,R_m_g_S,R_all_species_C,R_m_g_C,R_all_species_H,R_m_g_H):
    headers = pd.DataFrame([["R_ST","R_S_m","R_S_g","R_S_S2-","R_S_S2","R_S_OCS","R_S_H2S","R_S_SO2","R_S_S6+","R_S_H2Smol","a_S_g_m",
                            "R_CT","R_C_m","R_C_g","R_C_CO2","R_C_CO","R_C_CH4",'R_C_OCS',"R_C_COmol","R_C_CH4mol","R_C_CO2mol","R_C_CO32-","a_C_g_m",
                            "R_HT","R_H_m","R_H_g","R_H_H2O","R_H_H2","R_H_CH4","R_H_H2S","R_H_H2mol","R_H_CH4mol","R_H_H2Smol","R_H_H2Omol","R_H_OH-","a_H_g_m"]])
    values = pd.DataFrame([[R['S'],R_m_g_S["R_m"],R_m_g_S["R_g"],R_all_species_S["A"],R_all_species_S["B"],R_all_species_S["C"],R_all_species_S["D"],R_all_species_S["E"],R_all_species_S["F"],R_all_species_S["G"],R_m_g_S["R_g"]/R_m_g_S["R_m"],
                           R['C'],R_m_g_C["R_m"],R_m_g_C["R_g"],R_all_species_C["A"],R_all_species_C["B"],R_all_species_C["C"],R_all_species_C["D"],R_all_species_C["E"],R_all_species_C["F"],R_all_species_C["G"],R_all_species_C["H"],R_m_g_C["R_g"]/R_m_g_C["R_m"],
                           R['H'],R_m_g_H["R_m"],R_m_g_H["R_g"],R_all_species_H["A"],R_all_species_H["B"],R_all_species_H["C"],R_all_species_H["D"],R_all_species_H["E"],R_all_species_H["F"],R_all_species_H["G"],R_all_species_H["H"],R_all_species_H["I"],R_m_g_H["R_g"]/R_m_g_H["R_m"]]])
    return headers, values
#def results_table_isotope_a_D():
#    return headers, values
def results_table_isotope_d(R,R_all_species_S,R_m_g_S,R_all_species_C,R_m_g_C,R_all_species_H,R_m_g_H):
    headers = pd.DataFrame([["d_ST","d_S_m","d_S_g","d_S_S2-","d_S_S2","d_S_OCS","d_S_H2S","d_S_SO2","d_S_S6+","d_S_H2Smol","D_S_g_m",
                            "d_CT","d_C_m","d_C_g","d_C_CO2","d_C_CO","d_C_CH4",'d_C_OCS',"d_C_COmol","d_C_CH4mol","d_C_CO2mol","d_C_CO32-","D_C_g_m",
                            "d_HT","d_H_m","d_H_g","d_H_H2O","d_H_H2","d_H_CH4","d_H_H2S","d_H_H2mol","d_H_CH4mol","d_H_H2Smol","d_H_H2Omol","d_H_OH-","D_H_g_m"]])
    values = pd.DataFrame([[iso.ratio2delta("VCDT",34,'S',R['S']),iso.ratio2delta("VCDT",34,'S',R_m_g_S["R_m"]),iso.ratio2delta("VCDT",34,'S',R_m_g_S["R_g"]),iso.ratio2delta("VCDT",34,'S',R_all_species_S["A"]),iso.ratio2delta("VCDT",34,'S',R_all_species_S["B"]),iso.ratio2delta("VCDT",34,'S',R_all_species_S["C"]),iso.ratio2delta("VCDT",34,'S',R_all_species_S["D"]),iso.ratio2delta("VCDT",34,'S',R_all_species_S["E"]),iso.ratio2delta("VCDT",34,'S',R_all_species_S["F"]),iso.ratio2delta("VCDT",34,'S',R_all_species_S["G"]),iso.alpha2Delta((R_m_g_S["R_g"]/R_m_g_S["R_m"])),
                            iso.ratio2delta("VPDB",13,'C',R['C']),iso.ratio2delta("VPDB",13,'C',R_m_g_C["R_m"]),iso.ratio2delta("VPDB",13,'C',R_m_g_C["R_g"]),iso.ratio2delta("VPDB",13,'C',R_all_species_C["A"]),iso.ratio2delta("VPDB",13,'C',R_all_species_C["B"]),iso.ratio2delta("VPDB",13,'C',R_all_species_C["C"]),iso.ratio2delta("VPDB",13,'C',R_all_species_C["D"]),iso.ratio2delta("VPDB",13,'C',R_all_species_C["E"]),iso.ratio2delta("VPDB",13,'C',R_all_species_C["F"]),iso.ratio2delta("VPDB",13,'C',R_all_species_C["G"]),iso.ratio2delta("VPDB",13,'C',R_all_species_C["H"]),iso.alpha2Delta((R_m_g_C["R_g"]/R_m_g_C["R_m"])),
                            iso.ratio2delta("VSMOW",2,'H',R['H']),iso.ratio2delta("VSMOW",2,'H',R_m_g_H["R_m"]),iso.ratio2delta("VSMOW",2,'H',R_m_g_H["R_g"]),iso.ratio2delta("VSMOW",2,'H',R_all_species_H["A"]),iso.ratio2delta("VSMOW",2,'H',R_all_species_H["B"]),iso.ratio2delta("VSMOW",2,'H',R_all_species_H["C"]),iso.ratio2delta("VSMOW",2,'H',R_all_species_H["D"]),iso.ratio2delta("VSMOW",2,'H',R_all_species_H["E"]),iso.ratio2delta("VSMOW",2,'H',R_all_species_H["F"]),iso.ratio2delta("VSMOW",2,'H',R_all_species_H["G"]),iso.ratio2delta("VSMOW",2,'H',R_all_species_H["H"]),iso.ratio2delta("VSMOW",2,'H',R_all_species_H["I"]),iso.alpha2Delta((R_m_g_H["R_g"]/R_m_g_H["R_m"]))]])
    return headers, values


###############################
### options from setup file ###
###############################
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
    if models.loc["setup","option"] == "False":
        return models
    elif models.loc["setup","option"] == "True":
        # species
        if models.loc["COH_species","option"] == "setup":
            models.loc["COH_species","option"] = setup.loc[run,"COH_species"]
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

##################################################
### calculate the pressure of vapor saturation ###
##################################################
def calc_Pvsat(setup,models=mdv.default_models,first_row=0,last_row=None,p_tol=1.e-1,nr_step=1.,nr_tol=1.e-9):
    
    """ 
    Calculates the pressure of vapor saturation for multiple melt compositions given volatile-free melt composition, volatile content, temperature, and an fO2 estimate.


    Parameters
    ----------
    setup: pandas.DataFrame
        Dataframe with melt compositions to be used, requires following headers: 
        Sample, T_C, 
        DNNO or DFMQ or logfO2 or (Fe2O3 and FeO) or Fe3FeT or S6ST
        SiO2, TiO2, Al2O3, (Fe2O3T or FeOT unless Fe2O3 and FeO given), MnO, MgO, CaO, Na2O, K2O, P2O5, 
        H2O and/or CO2ppm and/or STppm and/or Xppm
        Note: concentrations (unless otherwise stated) are in wt%
    
    Optional:
    models: pandas.DataFrame
        Dataframe of options for different models.
    first_row: float
        Integer of the first row in the setup file to run (note the first row under the headers is row 0). Default = 0  
    last_row: float
        Integer of the last row in the setup file to run (note the first row under the headers is row 0). Default = length of setup
    p_tol: float
        Required tolerance for convergence of Pvsat in bars. Default = 1.e-1
    nr_step: float
        Step size for Newton-Raphson solver for melt speciation (this can be made smaller if there are problems with convergence.). Default = 1
    nr_tol: float
        Tolerance for the Newton-Raphson solver for melt speciation in weight fraction (this can be made larger if there are problems with convergence). Default = 1.e-9

    Returns
    -------
    results: pandas.DataFrame

    Outputs
    -------
    results_saturation_pressures: csv file (if output csv = yes in models)

    """
    if last_row == None:
        last_row = len(setup)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        PT={"T":setup.loc[run,"T_C"]}
        melt_wf=mg.melt_comp(run,setup)
        melt_wf['CO2'] = setup.loc[run,"CO2ppm"]/1000000.
        melt_wf["H2OT"] = setup.loc[run,"H2O"]/100.
        if "sulf_XFe" in setup:
            melt_wf["sulf_XFe"] = setup.loc[run,"sulf_XFe"]
        if "sulf_XCu" in setup:
            melt_wf["sulf_XCu"] = setup.loc[run,"sulf_XCu"]
        if "sulf_XNi" in setup:
            melt_wf["sulf_XNi"] = setup.loc[run,"sulf_XNi"]

        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)
 
        # calculate Pvsat assuming only H2O CO2 in vapour and melt
        #if setup.loc[run,"Fe3FeT"] > 0.:
        #    melt_wf['Fe3FeT'] = setup.loc[run,"Fe3FeT"]
        #else:
        #    melt_wf['Fe3FeT'] = 0.
        P_sat_H2O_CO2_only, P_sat_H2O_CO2_result = c.P_sat_H2O_CO2(PT,melt_wf,models,p_tol,nr_step,nr_tol)

        if models.loc["calc_sat","option"] == "fO2_fX":
            P_sat_fO2_fS2_result = c.P_sat_fO2_fS2(PT,melt_wf,models,p_tol)
            PT["P"] = P_sat_fO2_fS2_result["P_tot"]
        else:
            wm_ST = setup.loc[run,"STppm"]/1000000.
        melt_wf['ST'] = wm_ST
        melt_wf['CT'] = (melt_wf['CO2']/mdv.species.loc['CO2','M'])*mdv.species.loc['C','M']
        melt_wf['HT'] = (melt_wf['H2OT']/mdv.species.loc['H2O','M'])*(2.*mdv.species.loc['H','M'])
        wm_X = setup.loc[run,"Xppm"]/1000000.
        melt_wf['XT'] = wm_X
        if models.loc["bulk_composition","option"] == "melt-only":
            bulk_wf = {"H":(2.*mdv.species.loc["H","M"]*melt_wf["H2OT"])/mdv.species.loc["H2O","M"],"C":(mdv.species.loc["C","M"]*melt_wf["CO2"])/mdv.species.loc["CO2","M"],"S":wm_ST, "X":wm_X}
        else:
            raise TypeError('This is not currently possible')
        if models.loc["sulfur_is_sat","option"] == "yes":
            if melt_wf["XT"] > 0.:
                raise TypeError('This is not currently possible')
            P_sat_, conc, frac  = c.fO2_P_VSA(PT,melt_wf,models,nr_step,nr_tol,p_tol)
        elif models.loc["sulfur_saturation","option"] == "False":
            P_sat_, conc, frac = c.P_sat(PT,melt_wf,models,p_tol,nr_step,nr_tol)
        elif models.loc["sulfur_saturation","option"] == "True":
            if melt_wf["XT"] > 0.:
                raise TypeError('This is not currently possible')
            P_sat_, conc, frac = c.P_VSA(PT,melt_wf,models,nr_step,nr_tol,p_tol)
        PT["P"] = P_sat_
        melt_wf["H2OT"] = conc["wm_H2O"]
        melt_wf["CO2"] = conc["wm_CO2"]
        melt_wf["S2-"] = conc["wm_S2m"]
        melt_wf["Fe3FeT"] = conc["Fe3FeT"]
        #if models.loc["sulfur_is_sat","option"] == "yes":
        #    melt_wf["Fe3FeT"] = frac["Fe3FeT"]
        #else:
        #    melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,models)
        
        sulf_sat_result = c.sulfur_saturation(PT,melt_wf,models)
        # gas_mf = {"O2":mg.xg_O2(PT,melt_wf,models),"CO":mg.xg_CO(PT,melt_wf,models),"CO2":mg.xg_CO2(PT,melt_wf,models),"H2":mg.xg_H2(PT,melt_wf,models),"H2O":mg.xg_H2O(PT,melt_wf,models),"CH4":mg.xg_CH4(PT,melt_wf,models),"S2":mg.xg_S2(PT,melt_wf,models),"SO2":mg.xg_SO2(PT,melt_wf,models),"H2S":mg.xg_H2S(PT,melt_wf,models),"OCS":mg.xg_OCS(PT,melt_wf,models),"X":mg.xg_X(PT,melt_wf,models),"Xg_t":mg.Xg_tot(PT,melt_wf,models),"wt_g":0.}     
        melt_comp = mg.melt_normalise_wf(melt_wf,"yes","no")  
        
        # create results
        results_headers_table_sample_name, results_values_table_sample_name = results_table_sample_name(setup,run)
        results_headers_table_melt_comp_etc, results_values_table_melt_comp_etc = results_table_melt_comp_etc(PT,melt_comp,conc,frac,melt_wf)
        results_headers_table_model_options, results_values_table_model_options = results_table_model_options(models)    
        results_headers_table_f_p_xg_y_M_C_K_d, results_values_table_f_p_xg_y_M_C_K_d = results_table_f_p_xg_y_M_C_K_d(PT,melt_wf,models)
        results_headers_table_sat, results_values_table_sat = results_table_sat(sulf_sat_result,PT,melt_wf,models)
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
        
        if models.loc["print status","option"] == "True":
            print(n, setup.loc[run,"Sample"],PT["P"])
    
    results.columns = results.iloc[0]
    results = results[1:]  
    if models.loc["output csv","option"] == "True":
        results.to_csv('results_saturation_pressures.csv', index=False, header=True)
    
    return results

###################################
### cacluate re/degassing paths ###
###################################
def calc_gassing(setup,models=mdv.default_models,run=0,nr_step=1.,nr_tol=1.e-9,dp_step="auto",psat_tol=0.1,dwtg=1.e-6,i_nr_step=1.e-1,i_nr_tol=1.-9,nr_step_eq=1.):
     
    """ 
    Calculates the pressure of vapor saturation for multiple melt compositions given volatile-free melt composition, volatile content, temperature, and an fO2 estimate.


    Parameters
    ----------
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

    Optional:
    run: float
        Integer of the row in the setup file to run (note the first row under the headers is row 0). Default = 0
    nr_step: float
        Step size for Newton-Raphson solver for melt speciation (typically 1 is fine, but this can be made smaller if there are problems with convergence).
    nr_tol: float
        Tolerance for the Newton-Raphson solver for melt speciation in weight fraction (can be increased if there are problems with convergence). Default = 1.e-6
    dp_step: float
        Pressure step size for gassing calculation in bars. Default = 10    
    psat_tol: float
        Required tolerance for convergence of Pvsat in bars. Default = 0.1
    dwtg: float
        Amount of gas to add at each step if regassing in an open-system in wt fraction total system. Default = 1.e-6
    i_nr_step: float
        Step-size for newton-raphson convergence for isotopes (can be increased if there are problems with convergence). Default = 1e.-1
    i_nr_tol: float
        Tolerance for newton-raphson convergence for isotopes (can be increased if there are problems with convergence). Default = 1.e-9

    Returns
    -------
    results: pandas.DataFrame

    Outputs
    -------
    If output csv = yes in models
    results_gassing_chemistry: csv file

    """

    if models.loc["print status","option"] == "True":
        print(setup.loc[run,"Sample"])

    # check if any options need to be read from the setup file rather than the models file
    models = options_from_setup(run,models,setup)

    if models.loc["fO2","option"] != "Kress91A":
        raise TypeError("Change 'fO2' option in models to 'Kress91A' (other fO2 options are not currently supported)")

    # set T and volatile composition of the melt
    PT={"T":setup.loc[run,"T_C"]}
    melt_wf = mg.melt_comp(run,setup)
    melt_wf['CO2'] = setup.loc[run,"CO2ppm"]/1000000.
    melt_wf["H2OT"] = setup.loc[run,"H2O"]/100.
    melt_wf["ST"] = setup.loc[run,"STppm"]/1000000.
    melt_wf["H2"] = 0.
    melt_wf["XT"] = setup.loc[run,"Xppm"]/1000000.
    melt_wf["CT"] = (melt_wf["CO2"]/mdv.species.loc["CO2","M"])*mdv.species.loc["C","M"]
    melt_wf["HT"] = (2.*melt_wf["H2OT"]/mdv.species.loc["H2O","M"])*mdv.species.loc["H","M"]
    melt_wf["ST"] = melt_wf["ST"]
    if "S6ST" in setup:
        melt_wf["S6ST"] = setup.loc[run,"S6ST"]
    if "sulf_XFe" in setup:
        melt_wf["sulf_XFe"] = setup.loc[run,"sulf_XFe"]
    if "sulf_XCu" in setup:
        melt_wf["sulf_XCu"] = setup.loc[run,"sulf_XCu"]
    if "sulf_XNi" in setup:
        melt_wf["sulf_XNi"] = setup.loc[run,"sulf_XNi"]

    # Calculate saturation pressure for composition given in setup file
    if models.loc["COH_species","option"] == "H2O-CO2 only":  
        P_sat_, P_sat_H2O_CO2_result = c.P_sat_H2O_CO2(PT,melt_wf,models,psat_tol,nr_step,nr_tol)
        conc = {"wm_H2O":P_sat_H2O_CO2_result["wm_H2O"], "wm_CO2":P_sat_H2O_CO2_result["wm_CO2"], "wm_H2":0., "wm_CO":0., "wm_CH4":0., "wm_H2S":0., "wm_S2m":0., "wm_S6p":0., "ST": 0.}
        frac = c.melt_species_ratios(conc)
    else:
        P_sat_, conc, frac = c.P_sat(PT,melt_wf,models,psat_tol,nr_step,nr_tol)
    PT["P"] = P_sat_
    if models.loc["print status","option"] == "True":
        print("T=",PT["T"],"P=",PT["P"],datetime.datetime.now())

    # update melt composition at saturation pressure, check for sulfur saturation, and calculate some things
    melt_wf["H2OT"] = conc["wm_H2O"]
    melt_wf["CO2"] = conc["wm_CO2"]
    melt_wf["CO"] = conc["wm_CO"]
    melt_wf["CH4"] = conc["wm_CH4"]
    melt_wf["H2"] = conc["wm_H2"]
    melt_wf["S2-"] = conc["wm_S2m"]
    melt_wf["S6+"] = conc["wm_S6p"]
    melt_wf["H2S"] = conc["wm_H2S"]
    melt_wf["Fe3FeT"] = conc['Fe3FeT']
    melt_wf["S6ST"] = mg.S6ST(PT,melt_wf,models)
    sulf_sat_result = c.sulfur_saturation(PT,melt_wf,models)    
    wm_CO2eq, wm_H2Oeq = mg.melt_H2O_CO2_eq(melt_wf)
    melt_comp = mg.melt_normalise_wf(melt_wf,"yes","no")
    
    # Set bulk composition
    bulk_comp = c.bulk_composition(run,PT,melt_wf,setup,models)
    bulk_wf = {"C":bulk_comp["wt_C"],"H":bulk_comp["wt_H"],"O":bulk_comp["wt_O"],"S":bulk_comp["wt_S"],"Fe":bulk_comp["wt_Fe"],"Wt":bulk_comp["Wt"],"X":bulk_comp["wt_X"]}

    # set system and initial guesses
    system = eq.set_system(melt_wf,models)
    guesses = eq.initial_guesses(run,PT,melt_wf,setup,models,system)

    # create results
    results_headers_table_sample_name, results_values_table_sample_name = results_table_sample_name(setup,run)
    results_headers_table_melt_comp_etc, results_values_table_melt_comp_etc = results_table_melt_comp_etc(PT,melt_comp,conc,frac,melt_wf)
    results_headers_table_model_options, results_values_table_model_options = results_table_model_options(models)    
    results_headers_table_f_p_xg_y_M_C_K_d, results_values_table_f_p_xg_y_M_C_K_d = results_table_f_p_xg_y_M_C_K_d(PT,melt_wf,models)
    results_headers_table_sat, results_values_table_sat = results_table_sat(sulf_sat_result,PT,melt_wf,models)
    results_headers_table_melt_vol = results_table_melt_vol() # "H2OT-eq_wtpc","CO2T-eq_ppmw","ST_ppmw","X_ppmw"
    results_values_table_melt_vol = pd.DataFrame([[wm_H2Oeq*100.,wm_CO2eq*1000000.,conc["wm_ST"]*1000000.,melt_wf["XT"]*1000000.]])
    results_headers_table_wtg_etc = pd.DataFrame([["wt_g_wtpc","wt_g_O_wtf","wt_g_C_wtf","wt_g_H_wtf","wt_g_S_wtf","wt_g_X_wtf","wt_O_wtpc","wt_C_wtpc","wt_H_wtpc","wt_S_wtpc","wt_X_wtpc","Solving species",'mass balance C','mass balance O','mass balance H','mass balance S']])
    results_values_table_wtg_etc = pd.DataFrame([[bulk_comp["wt_g"]*100.,"","","","","",bulk_wf["O"]*100.,bulk_wf["C"]*100.,bulk_wf["H"]*100.,bulk_wf["S"]*100.,bulk_wf["X"]*100.,"","","","",""]])
    if models.loc["gassing_style","option"] == "open" and models.loc["gassing_direction","option"] == "degas":
        results_headers_table_open_all_gas = results_table_open_all_gas()
        results_values_table_open_all_gas = pd.DataFrame([[mg.xg_O2(PT,melt_wf,models),mg.xg_H2(PT,melt_wf,models),mg.xg_H2O(PT,melt_wf,models),mg.xg_S2(PT,melt_wf,models),mg.xg_SO2(PT,melt_wf,models),mg.xg_H2S(PT,melt_wf,models),mg.xg_CO2(PT,melt_wf,models),mg.xg_CO(PT,melt_wf,models),mg.xg_CH4(PT,melt_wf,models),mg.xg_OCS(PT,melt_wf,models),mg.xg_X(PT,melt_wf,models),mg.gas_CS(PT,melt_wf,models)]])
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
        results1 = pd.DataFrame([[PT["P"],PT["T"],mg.xg_O2(PT,melt_wf,models),mg.xg_CO(PT,melt_wf,models),mg.xg_CO2(PT,melt_wf,models),mg.xg_H2(PT,melt_wf,models),mg.xg_H2O(PT,melt_wf,models),mg.xg_CH4(PT,melt_wf,models),mg.xg_S2(PT,melt_wf,models),mg.xg_SO2(PT,melt_wf,models),mg.xg_H2S(PT,melt_wf,models),mg.xg_OCS(PT,melt_wf,models),wt_g_,
melt_wf["CO2"],melt_wf["H2OT"],0,(mg.wm_S(PT,melt_wf,models)/100),(mg.wm_SO3(PT,melt_wf,models)/100),melt_wf["ST"],melt_wf["Fe3FeT"],mg.S6ST(PT,melt_wf,models),
mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,models),"FMQ"),mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,models),"NNO"),SCSS_,sulfide_sat,SCAS_,sulfate_sat,
R_S_S2_,R_S_SO4_,"","","","",R_i["S"],"",ratio2delta("VCDT",R_S_S2_),ratio2delta("VCDT",R_S_SO4_),"","","","",ratio2delta("VCDT",R_i["S"]),"",
a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_,""]])
        results_isotopes = pd.concat([results_isotopes, results1], ignore_index=True) 
        if models.loc["output csv","option"] == "True":
            results_isotopes.to_csv('results_gassing_isotopes.csv', index=False, header=False)
    
    if dp_step == "auto":
        dp_step_choice = "auto"
        if models.loc["gassing_style","option"] == "open":
            dp_step = 1.
        else:
            if PT["P"] > 5000.:
                dp_step = 500.
            elif PT["P"] > 200.:
                dp_step = 100.
            elif PT["P"] > 50.:
                dp_step = 10.
            else:
                dp_step = 1.
    else:
        dp_step_choice = "user"

    if models.loc["P_variation","option"] == "polybaric":
        # pressure ranges and options
        starting_P = models.loc["starting_P","option"]
        if starting_P == "set":
            initial = int(setup.loc[run,"P_bar"])
        else:
            if models.loc["gassing_direction","option"] == "degas":
                answer = math.floor(PT["P"]/dp_step)
                initial = round(answer*dp_step)
            elif models.loc["gassing_direction","option"] == "regas":
                answer = math.ceil(PT["P"]/dp_step)
                initial = round(answer*dp_step)
        if models.loc["gassing_direction","option"] == "degas":
            #step = int(-1*dp_step) # pressure step in bars
            if "final_P" in setup:
                final = int(setup.loc[run,"final_P"])
            else:
                final = 1.
        elif models.loc["gassing_direction","option"] == "regas":
            #step = int(dp_step)
            final = int(setup.loc[run,"final_P"])
    elif models.loc["T_variation","option"] == "polythermal": # temperature ranges and options
        PT["P"] = setup.loc[run,"P_bar"]
        final = int(setup.loc[run,"final_T"])
        if setup.loc[run,"final_T"] > setup.loc[run,"T_C"]:
            initial = int(round(PT["T"])) 
            #step = int(dp_step) # temperature step in 'C
        elif setup.loc[run,"final_T"] < setup.loc[run,"T_C"]:
            initial = int(round(PT["T"]))
            #step = int(-1.*dp_step) # temperature step in 'C
    
    # add some gas to the system if doing open-system regassing
    if models.loc["gassing_direction","option"] == "regas" and models.loc["gassing_style","option"] == "open":
        gas_mf = {"O2":mg.xg_O2(PT,melt_wf,models),"CO":mg.xg_CO(PT,melt_wf,models),"CO2":mg.xg_CO2(PT,melt_wf,models),"H2":mg.xg_H2(PT,melt_wf,models),"H2O":mg.xg_H2O(PT,melt_wf,models),"CH4":mg.xg_CH4(PT,melt_wf,models),"S2":mg.xg_S2(PT,melt_wf,models),"SO2":mg.xg_SO2(PT,melt_wf,models),"H2S":mg.xg_H2S(PT,melt_wf,models),"OCS":mg.xg_OCS(PT,melt_wf,models),"X":mg.xg_X(PT,melt_wf,models),"Xg_t":mg.Xg_tot(PT,melt_wf,models),"wt_g":0.}
        new_comp = c.new_bulk_regas_open(PT,melt_wf,bulk_wf,gas_mf,dwtg,models)
        bulk_wf = {"C":new_comp['wt_C'],"H":new_comp['wt_H'],"O":new_comp['wt_O'],"S":new_comp['wt_S'],"Fe":new_comp['wt_Fe'],"X":new_comp['wt_X'],"Wt":new_comp['Wt']}
    
    # run over different pressures #
    number_of_step = 0.
    
    PT['P'] = initial
    while PT["P"] > 1.:
    #for i in range(initial,final,step): # P is pressure in bars or T is temperature in 'C
        number_of_step = number_of_step + 1.
        eq_Fe = models.loc["eq_Fe","option"]
        guesses_original = guesses # store original guesses in case the calculation needs to be restarted

        if dp_step_choice == "auto":
            if models.loc["gassing_style","option"] == "open":
                dp_step = 1.
            else:
                if PT["P"] > 5000.:
                    dp_step = 500.
                elif PT["P"] > 200.:
                    dp_step = 100.
                elif PT["P"] > 50.:
                    dp_step = 10.
                else:
                    dp_step = 1.
        
        if number_of_step == 1. and dp_step == 1.:
            dp_step = 0.
        
        if models.loc["gassing_direction","option"] == "regas":
            dp_step = -1.*dp_step
        
        if models.loc["P_variation","option"] == "polybaric": 
            #P = i - dp_step
            P = PT["P"] - dp_step
            if P < dp_step or P < 1.:
                P = 1.
            PT["P"] = P
        elif models.loc["T_variation","option"] == "polythermal":
            T = i - dp_step
            PT["T"] = T

        if models.loc["gassing_style","option"] == "open": # check melt is still vapor-saturated
            PT_ = {'P':PT['P'],'T':PT['T']}
            if models.loc["COH_species","option"] == "H2O-CO2 only":  
                P_sat_, P_sat_H2O_CO2_result = c.P_sat_H2O_CO2(PT_,melt_wf,models,psat_tol,nr_step,nr_tol)
                conc = {"wm_H2O":P_sat_H2O_CO2_result["wm_H2O"], "wm_CO2":P_sat_H2O_CO2_result["wm_CO2"], "wm_H2":0., "wm_CO":0., "wm_CH4":0., "wm_H2S":0., "wm_S2m":0., "wm_S6p":0., "ST": 0.}
                frac = c.melt_species_ratios(conc)
            else:
                P_sat_, conc, frac = c.P_sat(PT_,melt_wf,models,psat_tol,nr_step,nr_tol)
            if models.loc["gassing_direction","option"] == "degas":
                checkingP = PT['P']
                while P_sat_ < checkingP:
                    checkingP = checkingP - dp_step
                    PT_['P'] = checkingP
                    if models.loc["COH_species","option"] == "H2O-CO2 only":  
                        P_sat_, P_sat_H2O_CO2_result = c.P_sat_H2O_CO2(PT_,melt_wf,models,psat_tol,nr_step,nr_tol)
                        conc = {"wm_H2O":P_sat_H2O_CO2_result["wm_H2O"], "wm_CO2":P_sat_H2O_CO2_result["wm_CO2"], "wm_H2":0., "wm_CO":0., "wm_CH4":0., "wm_H2S":0., "wm_S2m":0., "wm_S6p":0., "ST": 0.}
                        frac = c.melt_species_ratios(conc)
                    else:
                        P_sat_, conc, frac = c.P_sat(PT_,melt_wf,models,psat_tol,nr_step,nr_tol)
                PT['P'] = checkingP

        if P_sat_ > PT["P"] or models.loc["gassing_direction","option"] == "regas":  
            # work out equilibrium partitioning between melt and gas phase
            xg, conc, melt_and_gas, guesses, new_models, solve_species, mass_balance = eq.mg_equilibrium(PT,melt_wf,bulk_wf,models,nr_step_eq,nr_tol,guesses)
            models = new_models
            #if xg["xg_O2"] == 1.0:
            #    print('tried resetting guesses')
            #    guesses = eq.initial_guesses(run,PT,melt_wf,setup,models,system)
            #    xg, conc, melt_and_gas, guesses, new_models, solve_species, mass_balance = eq.mg_equilibrium(PT,melt_wf,bulk_wf,models,nr_step_eq,nr_tol,guesses)
            #    models = new_models
            if models.loc["gassing_style","option"] == "closed":
                if xg["xg_O2"] == 1.0:
                    guesses = guesses_original
                    oldP = P + dp_step
                    if dp_step < 1. or dp_step == 1.:
                        if PT['P'] <= 10.:
                            PT["P"] = 1.
                            xg, conc, melt_and_gas, guesses, new_models, solve_species, mass_balance = eq.mg_equilibrium(PT,melt_wf,bulk_wf,models,nr_step_eq,nr_tol,guesses)
                            models = new_models
                            if xg["xg_O2"] == 1.0:
                                results.columns = results.iloc[0]
                                results = results[1:]  
                                results.reset_index(drop=True,inplace=True)
                                if models.loc["output csv","option"] == "True":
                                    results.to_csv('results_gassing_chemistry.csv', index=False, header=True) 
                                print("solver failed, calculation aborted at P = ", PT["P"], datetime.datetime.now())
                                return results                          
                        print('solver failed, increasing step size by factor 10')
                        dp_step = dp_step*10.
                    else:
                        print('solver failed, decreasing step size by factor 10')
                        dp_step = dp_step/10.
                    newP = oldP - dp_step
                    if newP < 1.:
                        newP = 1.
                    PT["P"] = newP
                    xg, conc, melt_and_gas, guesses, new_models, solve_species, mass_balance = eq.mg_equilibrium(PT,melt_wf,bulk_wf,models,nr_step_eq,nr_tol,guesses)
                    models = new_models
            if xg["xg_O2"] == 1.0:
                results.columns = results.iloc[0]
                results = results[1:]  
                results.reset_index(drop=True,inplace=True)
                if models.loc["output csv","option"] == "True":
                    results.to_csv('results_gassing_chemistry.csv', index=False, header=True)
                print("solver failed, calculation aborted at P = ", PT["P"], datetime.datetime.now())
                return results
            # gas composition
            gas_mf = {"O2":xg["xg_O2"],"CO":xg["xg_CO"],"S2":xg["xg_S2"],"CO2":xg["xg_CO2"],"H2O":xg["xg_H2O"],"H2":xg["xg_H2"],"CH4":xg["xg_CH4"],"SO2":xg["xg_SO2"],"H2S":xg["xg_H2S"],"OCS":xg["xg_OCS"],"X":xg["xg_X"],"Xg_t":xg["Xg_t"],"wt_g":melt_and_gas["wt_g"]}
        #else: # NEEDS SORTING ###
            #conc = eq.melt_speciation(PT,melt_wf,models,nr_step,nr_tol)
            #frac = c.melt_species_ratios(conc)
            #wm_ST_ = wm_S_ + wm_S6p_
            #S62 = S6T/S2m_ST
            #Fe3T = melt_wf["Fe3FeT"]
            #Fe32 = mg.overtotal2ratio(Fe3T)
            #xg = {"xg_O2":0., "xg_H2":0., "xg_S2":0., "xg_H2O":0., "xg_CO":0., "xg_CO2":0., "xg_SO2":0., "xg_CH4":0., "xg_H2S":0., "xg_OCS":0., "xg_X":0., "Xg_t":0.}
            #if number_of_step == 1:
            #    melt_and_gas = {}
            #melt_and_gas["wt_g_O"],melt_and_gas["wt_g_C"],melt_and_gas["wt_g_H"],melt_and_gas["wt_g_S"],melt_and_gas["wt_g_X"],melt_and_gas["wt_g"] = 0.,0.,0.,0.,0.
            #guesses = eq.initial_guesses(run,PT,melt_wf,setup,models,system)
            #solve_species = "na"
            #gas_mf = {"O2":xg["xg_O2"],"CO":xg["xg_CO"],"S2":xg["xg_S2"],"CO2":xg["xg_CO2"],"H2O":xg["xg_H2O"],"H2":xg["xg_H2"],"CH4":xg["xg_CH4"],"SO2":xg["xg_SO2"],"H2S":xg["xg_H2S"],"OCS":xg["xg_OCS"],"X":xg["xg_X"],"Xg_t":xg["Xg_t"],"wt_g":melt_and_gas["wt_g"]}
        
        if P_sat_ > PT["P"] or models.loc["gassing_direction","option"] == "regas": 
            if models.loc["gassing_style","option"] == "open" and models.loc["gassing_direction","option"] == "degas": 
                if number_of_step == 1:
                    gas_mf_all = gas_mf
                else:
                    gas_mf_all = c.gas_comp_all_open(gas_mf,gas_mf_all,models)
            if models.loc["COH_species","option"] == "H2O-CO2 only":
                Fe3T = melt_wf["Fe3FeT"]
                Fe32 = mg.overtotal2ratio(Fe3T)

        # set melt composition for forward calculation
        melt_wf["CO2"] = conc["wm_CO2"]
        melt_wf["H2OT"] = conc["wm_H2O"]
        melt_wf["H2"] = conc["wm_H2"]
        melt_wf["CO"] = conc["wm_CO"]
        melt_wf["CH4"] = conc["wm_CH4"]
        melt_wf["H2S"] = conc["wm_H2S"]
        melt_wf["S6+"] = (conc["wm_SO3"]/mdv.species.loc["SO3","M"])*mdv.species.loc["S","M"]
        melt_wf["S2-"] = conc["wm_S2m"]
        melt_wf["ST"] = conc["wm_ST"]
        melt_wf["XT"] = conc["wm_X"]
        melt_wf["Fe3FeT"] = conc["Fe3T"]
        if P_sat_ < PT["P"]:  
            bulk_comp = c.bulk_composition(run,PT,melt_wf,setup,models)
    
        # check for sulfur saturation and display warning in outputs
        sulf_sat_result = c.sulfur_saturation(PT,melt_wf,models)
        if sulf_sat_result["sulfide_sat"] == "yes":
            warning = "WARNING: sulfide-saturated"
        elif sulf_sat_result["sulfate_sat"] == "yes":
            warning = "WARNING: sulfate-saturated"
        else:
            warning = ""
        
        # calculate fO2
        if eq_Fe == "yes":
            fO2_ = mdv.f_O2(PT,melt_wf,models)
        elif eq_Fe == "no":
            fO2_ = (xg_O2_*mdv.y_O2(PT,models)*PT["P"])
        
        wm_CO2eq, wm_H2Oeq = mg.melt_H2O_CO2_eq(melt_wf)
        melt_comp = mg.melt_normalise_wf(melt_wf,"yes","no")
        frac = c.melt_species_ratios(conc)

        # store results
        results_headers_table_sample_name, results_values_table_sample_name = results_table_sample_name(setup,run)
        results_headers_table_melt_comp_etc, results_values_table_melt_comp_etc = results_table_melt_comp_etc(PT,melt_comp,conc,frac,melt_wf)
        results_headers_table_model_options, results_values_table_model_options = results_table_model_options(models)    
        results_headers_table_f_p_xg_y_M_C_K_d, results_values_table_f_p_xg_y_M_C_K_d = results_table_f_p_xg_y_M_C_K_d(PT,melt_wf,models)
        results_headers_table_sat, results_values_table_sat = results_table_sat(sulf_sat_result,PT,melt_wf,models)
        results_values_table_melt_vol = pd.DataFrame([[wm_H2Oeq*100.,wm_CO2eq*1000000.,conc["wm_ST"]*1000000.,melt_wf["XT"]*1000000.]])
        results_values_table_wtg_etc = pd.DataFrame([[melt_and_gas["wt_g"]*100.,melt_and_gas["wt_g_O"],melt_and_gas["wt_g_C"],melt_and_gas["wt_g_H"],melt_and_gas["wt_g_S"],melt_and_gas["wt_g_X"],melt_and_gas["wt_O"]*100.,melt_and_gas["wt_C"]*100.,melt_and_gas["wt_H"]*100.,melt_and_gas["wt_S"]*100.,melt_and_gas["wt_X"]*100.,solve_species,mass_balance['C'],mass_balance['O'],mass_balance['H'],mass_balance['S']]])
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
            A, B = iso.i2s6("S",PT,R_i,melt_wf,gas_mf,i_nr_step,i_nr_tol,guessx)
            RS_Sm, RS_H2S, RS_SO4, RS_S2, RS_SO2, RS_OCS = A
            RS_m, RS_g = B
            a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_ = iso.i2s6_S_alphas(PT)
            xg_SO3_ = 0.
            results2 = pd.DataFrame([[PT["P"],PT["T"],xg_O2_,xg_CO_,xg_CO2_,xg_H2_,xg_H2O_,xg_CH4_,xg_S2_,xg_SO2_,xg_SO3_,xg_H2S_,xg_OCS_,wt_g,
                        wm_CO2_,wm_H2O_,wm_H2_,wm_S_,wm_SO3_,wm_ST_,Fe3T,S6T,
                        mg.fO22Dbuffer(PT,fO2_,"FMQ"),mg.fO22Dbuffer(PT,fO2_,"NNO"),SCSS_,sulfide_sat,SCAS_,sulfate_sat,
                        RS_Sm, RS_SO4, RS_H2S, RS_SO2, RS_S2, RS_OCS, RS_m, RS_g, ratio2delta("VCDT",RS_Sm),ratio2delta("VCDT",RS_SO4),ratio2delta("VCDT",RS_H2S),ratio2delta("VCDT",RS_SO2),ratio2delta("VCDT",RS_S2),ratio2delta("VCDT",RS_OCS),ratio2delta("VCDT",RS_m),ratio2delta("VCDT",RS_g),a_H2S_S_,a_SO4_S_,a_S2_S_,a_SO2_S_,a_OCS_S_,RS_g/RS_m]])
            results_isotopes = pd.concat([results_isotopes, results2], ignore_index=True)
            if models.loc["output csv","option"] == "True":
                results_isotopes.to_csv('results_gassing_isotopes.csv', index=False, header=False)
        
        if models.loc["print status","option"] == "True":
            if(number_of_step % 100==0):
                print(PT["T"],PT["P"],mg.fO22Dbuffer(PT,fO2_,"FMQ",models),warning,datetime.datetime.now())

        # recalculate bulk composition if needed
        if models.loc["gassing_style","option"] == "open":
            results_me = mg.melt_elements(PT,melt_wf,bulk_wf,gas_mf,models)
            if models.loc["gassing_direction","option"] == "degas":
                Wt_ = bulk_wf['Wt']
                if results_me["wm_C"] < 1.e-6: # 1 ppm C
                    results_me["wm_C"] = 0.
                if results_me["wm_H"] < 1.e-6: # 1 ppm H
                    results_me["wm_H"] = 0.
                if results_me["wm_S"] < 1.e-6: # 1 ppm S
                    results_me["wm_S"] = 0.
                if results_me["wm_X"] < 1.e-6: # 1 ppm X
                    results_me["wm_X"] = 0.
                bulk_wf = {"C":results_me["wm_C"],"H":results_me["wm_H"],"O":results_me["wm_O"],"S":results_me["wm_S"],"X":results_me["wm_X"],"Fe":results_me["wm_Fe"],"Wt":(Wt_*(1. - melt_and_gas["wt_g"]))}
                melt_wf["CT"] = results_me["wm_C"]
                melt_wf["HT"] = results_me["wm_H"]
                melt_wf["ST"] = results_me["wm_S"]
                melt_wf["XT"] = results_me["wm_X"]
                system = eq.set_system(melt_wf,models)
            elif models.loc["gassing_direction","option"] == "regas":
                results_nbro = c.new_bulk_regas_open(PT,melt_wf,bulk_wf,gas_mf,dwtg,models)
                bulk_wf = {"C":results_nbro["wt_C"],"H":results_nbro["wt_H"],"O":results_nbro["wt_O"],"S":results_nbro["wt_S"],"X":results_nbro["wt_X"],"Fe":results_nbro["wt_Fe"],"Wt":results_nbro["Wt"]}
                #melt_wf["CT"] = results_nbro["wm_C"]
                #melt_wf["HT"] = results_nbro["wm_H"]
                #melt_wf["ST"] = results_nbro["wm_S"]
                #melt_wf["XT"] = results_nbro["wm_X"]
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
            
        if models.loc["gassing_direction","option"] == "regas":
            if (PT["P"] + dp_step) >= final:
                break

    results.columns = results.iloc[0]
    results = results[1:]
    results.reset_index(drop=True,inplace=True)  
    if models.loc["output csv","option"] == "True":
        results.to_csv('results_gassing_chemistry.csv', index=False, header=True)
    
    if models.loc["print status","option"] == "True":
        print("done", datetime.datetime.now())

    return results

#########################
### calculate isobars ###
#########################
def calc_isobar(setup,run=0,models=mdv.default_models,initial_P=1000.,final_P=10000.,step_P=1000.):
    
    if models.loc["COH_species","option"] == "H2O-CO2 only":
        PT={"T":setup.loc[run,"T_C"]}
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        # set up results table
        results = pd.DataFrame([["P_bar","H2O_wtpc","CO2_ppm"]])
    
        initial_P = int(initial_P)
        final_P = int(final_P+1)
        step_P = int(step_P)
        melt_wf=mg.melt_comp(run,setup)

        for n in range(initial_P,final_P,step_P):
            PT["P"] = n # pressure in bars
            results1 = c.calc_isobar_CO2H2O(PT,melt_wf,models)
            results = pd.concat([results, results1], ignore_index=True)
            if models.loc["print status","option"] == "True":
                print(setup.loc[run,"Sample"],n)
        results.columns = results.iloc[0]
        results = results[1:]

    else:
        raise TypeError("COH_species option must be H2O-CO2 only")
    if models.loc["output csv","option"] == "True":
        results.to_csv('results_isobars.csv', index=False, header=False)

    return results

#################################
### calculate pure solubility ###
#################################   
def calc_pure_solubility(setup,run=0,models=mdv.default_models,initial_P=5000.):
    if models.loc["print status","option"] == "True":
        print(setup.loc[run,"Sample"],initial_P)
    PT={"T":setup.loc[run,"T_C"]}

    # check if any options need to be read from the setup file rather than the models file
    models = options_from_setup(run,models,setup)

    # set up results table
    results = pd.DataFrame([["P_bar","H2O_wtpc","CO2_ppmw"]])
    
    initial_P = int(initial_P)
        
    for n in range(initial_P,1,-1):
        PT["P"] = n # pressure in bars
        melt_wf=mg.melt_comp(run,setup)
        results1 = c.calc_pure_solubility(PT,melt_wf,models)
        results = pd.concat([results, results1], ignore_index=True)
    
    results.columns = results.iloc[0]
    results = results[1:]
    if models.loc["output csv","option"] == "True":    
        results.to_csv('results_pure_solubility.csv', index=False, header=False)
    if models.loc["print status","option"] == "True":
        print("done")

    return results

######################################
### calculate solubility constants ###
######################################
# print capacities for multiple melt compositions in input file
def calc_sol_consts(setup,first_row=0,last_row=None,models=mdv.default_models):
    
    # set up results table
    results_headers_models = pd.DataFrame([["species X opt","Hspeciation opt",
                 "fO2 opt","NNObuffer opt","FMQbuffer opt",
                 "carbon dioxide opt","water opt","hydrogen opt","sulfide opt","sulfate opt","hydrogen sulfide opt","methane opt","carbon monoxide opt","species X solubility opt","Cspeccomp opt","Hspeccomp opt","Date"]])
    results_headers_values = pd.DataFrame([["Sample","Pressure (bar)","T ('C)","SiO2","TiO2","Al2O3","FeOT","MnO","MgO","CaO","Na2O","K2O","P2O5",
                "H2O","CO2 (ppm)","ST (ppm)","Fe3+/FeT","fO2 DFMQ","ln[C_CO2T]","ln[C_H2OT]","ln[C_S2-]","ln[C_S6+]","ln[C_H2S]","ln[C_H2]","ln[C_CO]","ln[C_CH4]","ln[C_X]","M_m_SO"]])
    results_headers = pd.concat([results_headers_values,results_headers_models],axis=1)
    
    if last_row == None:
        last_row = len(setup)
    
    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        PT={"T":setup.loc[run,"T_C"]}
        melt_wf=mg.melt_comp(run,setup)
        CO2 = setup.loc[run,"CO2ppm"]/1000000.
        CT = (CO2/mdv.species.loc["CO2","M"])*mdv.species.loc["C","M"]
        H2O = setup.loc[run,"H2O"]/100.
        HT = (H2O/mdv.species.loc["H2O","M"])*mdv.species.loc["H2","M"]
        XT = setup.loc[run,"Xppm"]/1000000.
        ST = setup.loc[run,"STppm"]/1000000.
        melt_wf['CO2']=CO2
        melt_wf["H2OT"]=H2O
        melt_wf["ST"]=ST
        melt_wf["X"]=XT
        melt_wf["XT"]=XT
        melt_wf["HT"]=HT
        melt_wf["CT"]=CT
        PT["P"] = setup.loc[run,"P_bar"]
        melt_wf["Fe3FeT"] = mg.Fe3FeT_i(PT,melt_wf,models)
        C_CO32 = mdv.C_CO3(PT,melt_wf,models)
        C_H2OT = mdv.C_H2O(PT,melt_wf,models)
        C_S2 = mdv.C_S(PT,melt_wf,models)
        C_S6 = mdv.C_SO4(PT,melt_wf,models)
        C_H2S = mdv.C_H2S(PT,melt_wf,models)
        C_H2 = mdv.C_H2(PT,melt_wf,models)
        C_CO = mdv.C_CO(PT,melt_wf,models)
        C_CH4 = mdv.C_CH4(PT,melt_wf,models)
        C_X = mdv.C_X(PT,melt_wf,models)
        fO2_ = mg.fO22Dbuffer(PT,mdv.f_O2(PT,melt_wf,models),"FMQ",models)
        M_m = mg.M_m_SO(melt_wf)
        melt_comp = mg.melt_normalise_wf(melt_wf,"yes","no")
                
        ### store results ###
        results_values_models = pd.DataFrame([[models.loc["species X","option"],models.loc["Hspeciation","option"], 
                models.loc["fO2","option"], models.loc["NNObuffer","option"], models.loc["FMQbuffer","option"],
                 models.loc["carbon dioxide","option"], models.loc["water","option"], models.loc["hydrogen","option"], models.loc["sulfide","option"], models.loc["sulfate","option"], models.loc["hydrogen sulfide","option"], models.loc["methane","option"], models.loc["carbon monoxide","option"], models.loc["species X solubility","option"], models.loc["Cspeccomp","option"], models.loc["Hspeccomp","option"],datetime.datetime.now()]])
        results_values_values = pd.DataFrame([[setup.loc[run,"Sample"],PT["P"],PT["T"],melt_comp["SiO2"]*100., melt_comp["TiO2"]*100., melt_comp["Al2O3"]*100., melt_comp["FeOT"]*100., melt_comp["MnO"]*100., melt_comp["MgO"]*100., melt_comp["CaO"]*100., melt_comp["Na2O"]*100., melt_comp["K2O"]*100., melt_comp["P2O5"]*100.,setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"STppm"],melt_wf["Fe3FeT"],fO2_,math.log(C_CO32),math.log(C_H2OT),math.log(C_S2),math.log(C_S6),math.log(C_H2S),math.log(C_H2),math.log(C_CO),math.log(C_CH4),math.log(C_X),M_m]])
        results1 = pd.concat([results_values_values,results_values_models],axis=1)
    
        if n == first_row:
            results = pd.concat([results_headers, results1])
        else:                         
            results = pd.concat([results, results1])
        
        if models.loc["print status","option"] == "True":
            print(n, setup.loc[run,"Sample"],math.log(C_CO32),math.log(C_H2OT),math.log(C_S2),math.log(C_S6),math.log(C_H2S),math.log(C_H2),math.log(C_CO),math.log(C_CH4),M_m)
    
    results.columns = results.iloc[0]
    results = results[1:]
                  
    if models.loc["output csv","option"] == "True":
        results.to_csv('capacities.csv', index=False, header=False)

    return results

#######################################      
### calculate fugacity coefficients ###
#######################################
def calc_fugacity_coefficients(setup,first_row=0,last_row=None,models=mdv.default_models):

    # set up results table
    results_headers_models = pd.DataFrame([["y_CO2 opt","y_SO2 opt","y_H2S opt","y_H2 opt","y_O2 opt","y_S2 opt","y_CO opt","y_CH4 opt","y_H2O opt","y_OCS opt","y_X opt","Date"]])
    results_headers_values = pd.DataFrame([["Sample","P_bar","T_C","yO2","yH2","yH2O","yS2","ySO2","yH2S","yCO2","yCO","yCH4","yOCS","yX"]])
    results_headers = pd.concat([results_headers_values,results_headers_models],axis=1)
    
    if last_row == None:
        last_row = len(setup)

    for n in range(first_row,last_row,1): # n is number of rows of data in conditions file
        run = n
        
        # check if any options need to be read from the setup file rather than the models file
        models = options_from_setup(run,models,setup)

        PT={"T":setup.loc[run,"T_C"]}
        PT["P"] = setup.loc[run,"P_bar"]
        
        ### store results ###
        results_values_models = pd.DataFrame([[models.loc["y_CO2","option"], models.loc["y_SO2","option"], models.loc["y_H2S","option"], models.loc["y_H2","option"], models.loc["y_O2","option"], models.loc["y_S2","option"], models.loc["y_CO","option"], models.loc["y_CH4","option"], models.loc["y_H2O","option"],models.loc["y_OCS","option"], models.loc["y_X","option"],datetime.datetime.now()]])
        results_values_values = pd.DataFrame([[setup.loc[run,"Sample"],PT["P"],PT["T"],mdv.y_O2(PT,models),mdv.y_H2(PT,models),mdv.y_H2O(PT,models),mdv.y_S2(PT,models),mdv.y_SO2(PT,models),mdv.y_H2S(PT,models),mdv.y_CO2(PT,models),mdv.y_CO(PT,models),mdv.y_CH4(PT,models),mdv.y_OCS(PT,models),mdv.y_X(PT,models)]])
        results1 = pd.concat([results_values_values,results_values_models],axis=1)
    
        if n == first_row:
            results = pd.concat([results_headers, results1])
        else:                         
            results = pd.concat([results, results1])
        
        if models.loc["print status","option"] == "True":
            print(n, setup.loc[run,"Sample"],PT["P"])
    
    results.columns = results.iloc[0]
    results = results[1:]  
    
    if models.loc["output csv","option"] == "True":
        results.to_csv('results_fugacity_coefficients.csv', index=False, header=True)

    return results

###############################               
### Use melt S oxybarometer ###
###############################        
def calc_melt_S_oxybarometer(setup,first_row=0,last_row=None,models=mdv.default_models,p_tol=0.1,nr_step=1.,nr_tol=1.e-9):
    """ 
    Calculates the range in oxygen fugacity based on the melt sulfur content for multiple melt compositions given volatile-free melt composition, volatile content, temperature, and either pressure or assumes Pvsat.


    Parameters
    ----------
    setup: pandas.DataFrame
        Dataframe with melt compositions to be used, requires following headers: 
        Sample, T_C, 
        SiO2, TiO2, Al2O3, (Fe2O3T or FeOT unless Fe2O3 and FeO given), MnO, MgO, CaO, Na2O, K2O, P2O5, 
        H2O and/or CO2ppm and/or STppm and/or Xppm
        Note: concentrations (unless otherwise stated) are in wt%
        Optional
        P_bar is pressure is given (otherwise calculation is at Pvsat)
        Fe3FeT is P_bar is specified
    
    Optional:
    models: pandas.DataFrame
        Dataframe of options for different models.
    first_row: float
        Integer of the first row in the setup file to run (note the first row under the headers is row 0). Default = 0  
    last_row: float
        Integer of the last row in the setup file to run (note the first row under the headers is row 0). Default = length of setup
    p_tol: float
        Required tolerance for convergence of Pvsat in bars. Default = 1.e-1
    nr_step: float
        Step size for Newton-Raphson solver for melt speciation (this can be made smaller if there are problems with convergence.). Default = 1
    nr_tol: float
        Tolerance for the Newton-Raphson solver for melt speciation in weight fraction (this can be made larger if there are problems with convergence). Default = 1.e-9

    Returns
    -------
    results: pandas.DataFrame

    Outputs
    -------
    fO2_range_from_S: csv file (if output csv = yes in models)

    """

    if last_row == None:
        last_row = len(setup)

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
        melt_wf["HT"] = ((setup.loc[run, "H2O"]/100.)/mdv.species.loc["H2O","M"])*mdv.species.loc["H2","M"]
        melt_wf["ST"] = setup.loc[run, "STppm"]/1000000.
        melt_wf["XT"] = setup.loc[run, "Xppm"]/1000000.
        melt_wf["CT"] = ((setup.loc[run, "CO2ppm"]/1000000.)/mdv.species.loc["CO2","M"])*mdv.species.loc["C","M"]
        if "sulf_XFe" in setup:
            melt_wf["sulf_XFe"] = setup.loc[run,"sulf_XFe"]
        if "sulf_XCu" in setup:
            melt_wf["sulf_XCu"] = setup.loc[run,"sulf_XCu"]
        if "sulf_XNi" in setup:
            melt_wf["sulf_XNi"] = setup.loc[run,"sulf_XNi"]

        if 'P_bar' in setup:    
            if setup.loc[run,"P_bar"] > 0.:
                PT["P"] = setup.loc[run,"P_bar"]
                melt_wf["Fe3FeT"] = setup.loc[run, "Fe3FeT"]
                sulfsat_results = c.fO2_range_from_S(PT,melt_wf,models)
                sulfsat_results["P_sat_sulf"] = setup.loc[run,"P_bar"]
                sulfsat_results["P_sat_anh"] = setup.loc[run,"P_bar"]
            else: 
                sulfsat_results = c.P_sat_sulf_anh(PT,melt_wf,models,p_tol,nr_step,nr_tol)
        else:
            sulfsat_results = c.P_sat_sulf_anh(PT,melt_wf,models,p_tol,nr_step,nr_tol)

        # create results
        results_headers_table_sample_name, results_values_table_sample_name = results_table_sample_name(setup,run)
        results_headers_table_model_options, results_values_table_model_options = results_table_model_options(models)  
        results_headers_T, results_values_T = pd.DataFrame([["T ('C)"]]), pd.DataFrame([[PT["T"]]])
        results_headers_table_melt_vol = results_table_melt_vol() # "H2OT-eq_wtpc","CO2T-eq_ppmw","ST_ppmw","X_ppmw"
        results_values_table_melt_vol = pd.DataFrame([[setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"STppm"],setup.loc[run,"Xppm"]]])
        results_headers_table_sulfsat = pd.DataFrame([["P (bar) sulf","S2- SCSS","sulfide saturated?","DFMQ-sulfide","fO2-sulfide","Fe3FeT-sulfide","S6ST-sulfide","P (bar) anh","S6+ SCAS","sulfate saturated?","DFMQ-sulfate","fO2-sulfate","Fe3FeT-sulfate","S6ST-sulfate"]])
        results_values_table_sulfsat = pd.DataFrame([[sulfsat_results["P_sat_sulf"],sulfsat_results["SCSS"],sulfsat_results["sulf_sat"],sulfsat_results["DFMQ_sulf"],sulfsat_results["fO2_sulf"],sulfsat_results["Fe3T_sulf"],sulfsat_results["S6T_sulf"],sulfsat_results["P_sat_anh"],sulfsat_results["SCAS"],sulfsat_results["anh_sat"],sulfsat_results["DFMQ_anh"],sulfsat_results["fO2_anh"],sulfsat_results["Fe3T_anh"],sulfsat_results["S6T_anh"]]])
        results_headers = pd.concat([results_headers_table_sample_name,results_headers_T,results_headers_table_melt_vol,results_headers_table_sulfsat,results_headers_table_model_options],axis=1)
        results1 = pd.concat([results_values_table_sample_name,results_values_T,results_values_table_melt_vol,results_values_table_sulfsat,results_values_table_model_options],axis=1)
    
        if n == first_row:
            results = pd.concat([results_headers, results1])
        else:                         
            results = pd.concat([results, results1])
        
        if models.loc["print status","option"] == "True":
            print(n, setup.loc[run,"Sample"],sulfsat_results["sulf_sat"],sulfsat_results["anh_sat"])
    
    results.columns = results.iloc[0]
    results = results[1:]  
        
    if models.loc["output csv","option"] == "True":
        results.to_csv('fO2_range_from_S.csv', index=False, header=True)

    return results

########################################
### measured parameters within error ### 
########################################
def calc_comp_error(setup,run,iterations,models=mdv.default_models):
    
    # set up results table
    results = pd.DataFrame([["Sample","T_C",
                  "SiO2","TiO2","Al2O3","FeOT","MnO","MgO","CaO","Na2O","K2O","P2O5",
                "H2O","CO2ppm","Xppm","STppm","Fe3FeT"]])
    results1 = pd.DataFrame([[setup.loc[run,"Sample"],setup.loc[run,"T_C"],setup.loc[run,"SiO2"],setup.loc[run,"TiO2"],setup.loc[run,"Al2O3"],setup.loc[run,"FeOT"],setup.loc[run,"MnO"],setup.loc[run,"MgO"],setup.loc[run,"CaO"],setup.loc[run,"Na2O"],setup.loc[run,"K2O"],setup.loc[run,"P2O5"],
                setup.loc[run,"H2O"],setup.loc[run,"CO2ppm"],setup.loc[run,"Xppm"],setup.loc[run,"STppm"], setup.loc[run,"Fe3FeT"]]])
                             
    results = pd.concat([results, results1], ignore_index=True)
    results1 = pd.DataFrame([["sds","",setup.loc[run,"SiO2_sd"],setup.loc[run,"TiO2_sd"],setup.loc[run,"Al2O3_sd"],setup.loc[run,"FeOT_sd"],setup.loc[run,"MnO_sd"],setup.loc[run,"MgO_sd"],setup.loc[run,"CaO_sd"],setup.loc[run,"Na2O_sd"],setup.loc[run,"K2O_sd"],setup.loc[run,"P2O5_sd"],
                setup.loc[run,"H2O_sd"],setup.loc[run,"CO2ppm_sd"],setup.loc[run,"Xppm_sd"],setup.loc[run,"STppm_sd"], setup.loc[run,"Fe3FeT_sd"]]])
                             
    results = pd.concat([results, results1], ignore_index=True)
    results1 = pd.DataFrame([["sd types","",setup.loc[run,"SiO2_sd_type"],setup.loc[run,"TiO2_sd_type"],setup.loc[run,"Al2O3_sd_type"],setup.loc[run,"FeOT_sd_type"],setup.loc[run,"MnO_sd_type"],setup.loc[run,"MgO_sd_type"],setup.loc[run,"CaO_sd_type"],setup.loc[run,"Na2O_sd_type"],setup.loc[run,"K2O_sd_type"],setup.loc[run,"P2O5_sd_type"],
                setup.loc[run,"H2O_sd_type"],setup.loc[run,"CO2ppm_sd_type"],setup.loc[run,"Xppm_sd_type"],setup.loc[run,"STppm_sd_type"], setup.loc[run,"Fe3FeT_sd_type"]]])
                             
    results = pd.concat([results, results1], ignore_index=True)
    for n in range(0,iterations,1): # n is number of rows of data in conditions file
        results1 = c.compositions_within_error(run,setup)
        results1 = pd.DataFrame([[run,setup.loc[run,"T_C"],results1["SiO2"],results1["TiO2"],results1["Al2O3"],results1["FeOT"],results1["MnO"],results1["MgO"],results1["CaO"],results1["Na2O"],results1["K2O"],results1["P2O5"],results1["H2O"],results1["CO2ppm"],results1["Xppm"],results1["STppm"],results1["Fe3FeT"]]])
        results = pd.concat([results, results1], ignore_index=True)
    
    results.columns = results.iloc[0]
    results = results[1:]
    if models.loc["output csv","option"] == "True":
        results.to_csv('random_compositions.csv', index=False, header=False)
    if models.loc["print status","option"] == "True":
        print(n, setup.loc[run,"Sample"],SiO2)
    
    return results