# model_dependent_variables.py

import pandas as pd
import numpy as np
import gmpy2 as gp
import math
import densityx as dx
import PySulfSat as ss

#import VolFe.melt_gas as mg

#from VolFe.melt_gas import *
import VolFe.melt_gas as mg

################
### Contents ###
################
# Solubility constants
# Solid/liquid saturation
# Equilibrium constants
# Fugacity coefficients
# Speciation 
# Isotope fractionation factors
# Constants

def make_models_df(models): 
    # Create the pandas DataFrame
    models = pd.DataFrame(models, columns=['type', 'option'])
    models = models.set_index('type')
    return models

# define default models
default_models = [['COH_species','yes_H2_CO_CH4_melt'],['H2S_m','True'],['species X','Ar'],['Hspeciation','none'],
              ['fO2','Kress91A'],['NNObuffer','Frost91'],['FMQbuffer','Frost91'],
              ['melt composition','Basalt'],['carbon dioxide','MORB_Dixon95'],['water','Basalt_Hughes24'],['hydrogen','Basalt_Hughes24'],['sulfide','ONeill21dil'],['sulfate','ONeill22dil'],['hydrogen sulfide','Basalt_Hughes24'],['methane','Basalt_Ardia13'],['carbon monoxide','Basalt_Hughes24'],['species X solubility','Ar_Basalt_HughesIP'],['Cspeccomp','Basalt'],['Hspeccomp','MORB_HughesIP'],
              ['SCSS','ONeill21hyd'],['SCAS','Zajacz19'],['sulfur_saturation','False'],['sulfur_is_sat','no'],['graphite_saturation','False'],
              ['ideal_gas','False'],['y_CO2','Shi92'],['y_SO2','Shi92_Hughes23'],['y_H2S','Shi92_Hughes24'],['y_H2','Shaw64'],['y_O2','Shi92'],['y_S2','Shi92'],['y_CO','Shi92'],['y_CH4','Shi92'],['y_H2O','Holland91'],['y_OCS','Shi92'],['y_X','ideal'],
              ['KHOg','Ohmoto97'],['KHOSg','Ohmoto97'],['KOSg','Ohmoto97'],['KOSg2','ONeill22'], ['KCOg','Ohmoto97'],['KCOHg','Ohmoto97'],['KOCSg','Moussallam19'],['KCOs','Holloway92'],['carbonylsulfide','COS'],
              ['bulk_composition','melt-only'],['starting_P','Pvsat'],['gassing_style','closed'],['gassing_direction','degas'],['P_variation','polybaric'],['eq_Fe','yes'],['solve_species','OCS'],
              ['beta_factors','Richet77'],['alpha_H_CH4v_CH4m','no fractionation'],['alpha_H_H2v_H2m','no fractionation'],['alpha_H_H2Ov_OHmm','Rust04'],['alpha_H_H2Ov_H2Om','Rust04'],['alpha_H_H2Sv_H2Sm','no fractionation'],
              ['alpha_C_CH4v_CH4m','no fractionation'],['alpha_C_COv_COm','no fractionation'],['alpha_C_CO2v_CO2T','LeePP'],['alpha_C_CO2v_CO2m','Blank93'],['alpha_C_CO2v_CO32mm','LeePP'],
              ['alpha_S_H2Sv_H2Sm','no fractionation'],['alpha_SO2_SO4','Fiege15'],['alpha_H2S_S','Fiege15'],
              ['density','DensityX'],['isotopes','no'],['T_variation','isothermal'],['crystallisation','no'],['mass_volume','mass'],['calc_sat','fO2_melt'],['bulk_O','exc_S'],['error',0.1],
              ['print status','False'],['output csv','True'],['setup','False'],['high precision','False']]
# Create the pandas DataFrame
default_models = pd.DataFrame(default_models, columns=['type', 'option'])
default_models = default_models.set_index('type')

# define default models for rhyolite
default_models_rhyolite = [['carbon dioxide','Rhyolite_Blank93'],['water','Rhyolite_HughesIP'],['hydrogen','Andesite_Hughes24'],
                  ['hydrogen sulfide','BasalticAndesite_Hughes24'],['species X solubility','Ar_Rhyolite_HughesIP'],
                  ['Cspeccomp','Rhyolite'],['Hspeccomp','Rhyolite_Zhang97']]

def check_default_options(models):
    
    def return_options(default,name,models):
        if name in models.index:
            variable = models.loc[name,'option']
            if variable == "default":
                    variable = default
        else:
            variable = default
        return variable

    # species
    COH_species = return_options(default_models.loc['COH_species','option'],'COH_species',models)
    H2S_m = return_options(default_models.loc['H2S_m','option'],'H2S_m',models)
    species_X = return_options(default_models.loc['species X','option'],'species X',models)
    Hspeciation = return_options(default_models.loc['Hspeciation','option'],'Hspeciation',models)
    # oxygen fugacity
    fO2 = return_options(default_models.loc['fO2','option'],'fO2',models)
    NNObuffer = return_options(default_models.loc['NNObuffer','option'],'NNObuffer',models)
    FMQbuffer = return_options(default_models.loc['FMQbuffer','option'],'FMQbuffer',models)
    # solubility constants
    melt_comp = return_options(default_models.loc['melt composition','option'],'melt composition',models)
    S2m = return_options(default_models.loc['sulfide','option'],'sulfide',models)
    S6p = return_options(default_models.loc['sulfate','option'],'sulfate',models)
    CH4 = return_options(default_models.loc['methane','option'],'methane',models)
    CO = return_options(default_models.loc['carbon monoxide','option'],'carbon monoxide',models)
    if melt_comp == 'Basalt':
        default_models_MC = default_models
    elif melt_comp == 'Rhyolite':
        default_models_MC = default_models_rhyolite
    CO2 = return_options(default_models_MC.loc['carbon dioxide','option'],'carbon dioxide',models)
    H2O = return_options(default_models_MC.loc['water','option'],'water',models)
    H2 = return_options(default_models_MC.loc['hydrogen','option'],'hydrogen',models)
    H2S = return_options(default_models_MC.loc['hydrogen sulfide','option'],'hydrogen sulfide',models)
    X = return_options(default_models_MC.loc['species X solubility','option'],'species X solubility',models)
    Cspec = return_options(default_models_MC.loc['Cspeccomp','option'],'Cspeccomp',models)
    Hspec = return_options(default_models_MC.loc['Hspeccomp','option'],'Hspeccomp',models)
    # saturation conditions
    SCSS = return_options(default_models.loc['SCSS','option'],'SCSS',models)
    SCAS = return_options(default_models.loc['SCAS','option'],'SCAS',models)
    sulfur_saturation = return_options(default_models.loc['sulfur_saturation','option'],'sulfur_saturation',models)
    sulfur_is_sat = return_options(default_models.loc['sulfur_is_sat','option'],'sulfur_is_sat',models)
    graphite_saturation = return_options(default_models.loc['graphite_saturation','option'],'graphite_saturation',models)
    # fugacity coefficients
    ideal_gas = return_options(default_models.loc['ideal_gas','option'],'ideal_gas',models)
    yCO2 = return_options(default_models.loc['y_CO2','option'],'y_CO2',models)
    ySO2 = return_options(default_models.loc['y_SO2','option'],'y_SO2',models)
    yH2S = return_options(default_models.loc['y_H2S','option'],'y_H2S',models)
    yH2 = return_options(default_models.loc['y_H2','option'],'y_H2',models)
    yO2 = return_options(default_models.loc['y_O2','option'],'y_O2',models)
    yS2 = return_options(default_models.loc['y_S2','option'],'y_S2',models)
    yCO = return_options(default_models.loc['y_CO','option'],'y_CO',models)
    yCH4 = return_options(default_models.loc['y_CH4','option'],'y_CH4',models)
    yH2O = return_options(default_models.loc['y_H2O','option'],'y_H2O',models)
    yOCS = return_options(default_models.loc['y_OCS','option'],'y_OCS',models)
    yX = return_options(default_models.loc['y_X','option'],'y_X',models)
    # equilibrium constants
    KHOg = return_options(default_models.loc['KHOg','option'],'KHOg',models)
    KHOSg = return_options(default_models.loc['KHOSg','option'],'KHOSg',models)
    KOSg = return_options(default_models.loc['KOSg','option'],'KOSg',models)
    KOSg2 = return_options(default_models.loc['KOSg2','option'],'KOSg2',models)
    KCOg = return_options(default_models.loc['KCOg','option'],'KCOg',models)
    KCOHg = return_options(default_models.loc['KCOHg','option'],'KCOHg',models)
    KOCSg = return_options(default_models.loc['KOCSg','option'],'KOCSg',models)
    KCOs = return_options(default_models.loc['KCOs','option'],'KCOs',models)
    OCS = return_options(default_models.loc['carbonylsulfide','option'],'carbonlysulfide',models)
    # degassing calculation
    bulk_composition = return_options(default_models.loc['bulk_composition','option'],'bulk_composition',models)
    starting_P = return_options(default_models.loc['starting_P','option'],'starting_P',models)
    gassing_style = return_options(default_models.loc['gassing_style','option'],'gassing_style',models)
    gassing_direction = return_options(default_models.loc['gassing_direction','option'],'gassing_direction',models)
    P_variation = return_options(default_models.loc['P_variation','option'],'P_variation',models)
    eq_Fe = return_options(default_models.loc['eq_Fe','option'],'eq_Fe',models)
    solve_species = return_options(default_models.loc['solve_species','option'],'solve_species',models)
    # other
    density = return_options(default_models.loc['density','option'],'density',models)
    isotopes = return_options(default_models.loc['isotopes','option'],'isotopes',models)
    T_variation = return_options(default_models.loc['T_variation','option'],'T_variation',models)
    crystallisation = return_options(default_models.loc['crystallisation','option'],'cystallisation',models)
    mass_volume = return_options(default_models.loc['mass_volume','option'],'mass_volume',models)
    calc_sat = return_options(default_models.loc['calc_sat','option'],'calc_sat',models)
    bulk_O = return_options(default_models.loc['bulk_O','option'],'bulk_O',models)
    error = return_options(default_models.loc['error','option'],'error',models)
    print_status = return_options(default_models.loc['print status','option'],'print status',models)
    output_csv = return_options(default_models.loc['output csv','option'],'output csv',models)
    setup = return_options(default_models.loc['setup','option'],'setup',models)
    precision = return_options(default_models.loc['high precision','option'],'high precision',models)

    models = [['COH_species',COH_species],['H2S_m',H2S_m],['species X',species_X],['Hspeciation',Hspeciation],
              ['fO2',fO2],['NNObuffer',NNObuffer],['FMQbuffer',FMQbuffer],
              ['melt composition',melt_comp],['carbon dioxide',CO2],['water',H2O],['hydrogen',H2],['sulfide',S2m],['sulfate',S6p],['hydrogen sulfide',H2S],['methane',CH4],['carbon monoxide',CO],['species X solubility',X],['Cspeccomp',Cspec],['Hspeccomp',Hspec],
              ['SCSS',SCSS],['SCAS',SCAS],['sulfur_saturation',sulfur_saturation],['sulfur_is_sat',sulfur_is_sat],['graphite_saturation',graphite_saturation],
              ['ideal_gas',ideal_gas],['y_CO2',yCO2],['y_SO2',ySO2],['y_H2S',yH2S],['y_H2',yH2],['y_O2',yO2],['y_S2',yS2],['y_CO',yCO],['y_CH4',yCH4],['y_H2O',yH2O],['y_OCS',yOCS],['y_X',yX],
              ['KHOg',KHOg],['KHOSg',KHOSg],['KOSg',KOSg],['KOSg2',KOSg2], ['KCOg',KCOg],['KCOHg',KCOHg],['KOCSg',KOCSg],['KCOs',KCOs],['carbonylsulfide',OCS],
              ['bulk_composition',bulk_composition],['starting_P',starting_P],['gassing_style',gassing_style],['gassing_direction',gassing_direction],['P_variation',P_variation],['eq_Fe',eq_Fe],['solve_species',solve_species],
              ['density',density],['isotopes',isotopes],['T_variation',T_variation],['crystallisation',crystallisation],['mass_volume',mass_volume],['calc_sat',calc_sat],['bulk_O',bulk_O],['error',error],
              ['print status',print_status],['output csv',output_csv],['setup',setup],['high precision',precision]]
    
    # Create the pandas DataFrame
    models = pd.DataFrame(models, columns=['type', 'option'])
    models = models.set_index('type')
    return models

def make_df_and_add_model_defaults(models):
    """

    Converts user-provided model configurations (e.g. ['carbon dioxide','MORB_Dixon95'],['hydrogen sulfide','basaltic andesite']
    into a structured pandas DataFrame, combined with default options for anything not specified
    

    Parameters
    ----------
    models : list of [str, str]
        A list of lists, where each inner list contains two elements: the model type (str)
        and the user-specified option (str) for that model type.

    Returns
    -------
    pandas.DataFrame
        A DataFrame where the first column is 'type', set as the index, and the second column
        is 'option', containing the user-specified option or the default option if none is provided.
        
    
    Model Parameters and Options
    ---------------------------------
    The following parameters can be overridden in models. 
    Any parameter can be set to 'setup', in which case the parameter is specified in the setup dataframe instead.
    
    
    ### Specifying species ###
    
    COH_species: Specifying what COH species are present in the melt and vapor.
        default: 'yes_H2_CO_CH4_melt' Include H2mol (if H present), COmol (if C present) and/or CH4mol (if H and C present) as dissolved melt species.
        Other options: 
        'no_H2_CO_CH4_melt' H2, CO and/or CH4 are insoluble in the melt but they are still present in the vapor (H2 in the vapor if H present, CO in the vapor if C present, CH4 in the vapor if both H and C present).
        'H2O-CO2 only' The only species present in the vapor are H2O and CO2 and in the melt are H2OT and CO2T (i.e., no CO, H2, and/or CH4 in the melt or vapor).
        
    H2S_m: Is H2S a dissolved melt species.
        default 'True' Include H2Smol as a dissolved melt species. 
        Other options:  
        'False' H2Smol is insoluble in the melt.
        
    species X: Chemical identity of species X, which defines its atomic mass.
        default 'Ar' Species X is argon (i.e., atomic mass of ~40).
        Other options:
        'Ne' Species X is Ne (i.e., atomic mass of ~20).
        Other noble gases not currently supported, but we can add them if you get in touch!
    
    ### Hspeciation: default 'none' Oxidised H in the melt only occurs as H2O species (i.e., no OH-).
        Other options:
        WORK IN PROGRESS

    
    ### Oxygen fugacity ###

    fO2: Model for parameterisation of relationship between fO2 and Fe3+/FeT
        default: 'Kress91A' Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
        Other options:
        'Kress91' Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
        'ONeill18' Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162 doi:10.1016/j.epsl.2018.10.0020012-821X
        'Borisov18' Eq. (4) from Borisov et al. (2018) CMP 173:98 doi:10.1007/s00410-018-1524-8  

    NNObuffer: Model for the parameterisation for the fO2 value of the NNO buffer.
        default: 'Frost91' Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance" doi:10.1515/9781501508684-004
        Only one option available currently, included for future development.

    FMQbuffer: Model for the parameterisation for the fO2 value of the FMQ buffer.
        default: 'Frost91' FM[beta]Q in Table 1 of Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance" doi:10.1515/9781501508684-004
        Other options:
        'ONeill87' O'Neill (1897) AmMin 72(1-2):67-75

    
    ### Models for solubility and speciation constants ###

    carbon dioxide: Model for the parameterisation of the CO2T solubility constant.
        default: 'MORB_Dixon95' Bullet (5) of summary from Dixon et al. (1995) JPet 36(6):1607-1631 doi:10.1093/oxfordjournals.petrology.a037267
        Other options:
        'Basalt_Dixon97' Eq. (7) from Dixon et al. (1997) AmMin 82(3-4)368-378 doi:10.2138/am-1997-3-415
        'NorthArchBasalt_Dixon97' Eq. (8) from Dixon et al. (1997) AmMin 82(3-4)368-378 doi:10.2138/am-1997-3-415
        'Basalt_Lesne11' Eq. (25,26) from Lesne et al. (2011) CMP 162:153-168 doi:10.1007/s00410-010-0585-0
        'VesuviusAlkaliBasalt_Lesne11' VES-9 in Table 4 from Lesne et al. (2011) CMP 162:153-168 doi:10.1007/s00410-010-0585-0
        'EtnaAlkaliBasalt_Lesne11' ETN-1 in Table 4 from Lesne et al. (2011) CMP 162:153-168 doi:10.1007/s00410-010-0585-0
        'StromboliAlkaliBasalt_Lense11' PST-9 in Table 4 from Lesne et al. (2011) CMP 162:153-168 doi:10.1007/s00410-010-0585-0
        'SunsetCraterAlkaliBasalt_Allison19' Sunset Crater in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
        'SVFVBasalticAndesite_Allison19' SVFV in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
        'ErebusPhonotephrite_Allison19' Erebus in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
        'VesuviusPhonotephrite_Allison19' Vesuvius in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
        'EtnaTrachybasalt_Allison19' Etna in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
        'StromboliAlkaliBasalt_Allison19' Stromboli in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
        'Basanite_Holloway94' Basanite in Table 5 from Holloway and Blank (1994) RiMG 30:187-230 doi:10.1515/9781501509674-012
        'Leucitite_Thibault94' Leucitite from Thibault & Holloway (1994) CMP 116:216-224 doi:10.1007/BF00310701
        'TholeiiteBasalt_Allison22' N72 basalt in Table 2 from Allison et al. (2022) CMP 177:40 doi:10.1007/s00410-022-01903-y
        'Rhyolite_Blank93' Fig.2 caption from Blank et al. (1993) EPSL 119:27-36 doi:10.1016/0012-821X(93)90004-S

    ### water: Model for the parameterisation for the H2O solubility constant.
        default: 'Basalt_Hughes24' Hughes et al. (2024) AmMin 109(3):422-438 based on data compiliation from Allison et al. (2022) CMP 177(3):40 doi:10.1007/s00410-022-01903-y
        Other options:
        'Rhyolite_HughesIP' Fig.SX from Hughes et al. (in prep) based on data in Fig. 3 of Blank et al. (1993)
    
    hydrogen: Model for the parameterisation of the H2 solubility constant.
        default: 'Basalt_Hughes24' Basalt H2 in Table S4 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
        Other options:
        'Andesite_Hughes24' Andesite H2 in Table S4 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
    
    sulfide: Model for the parameterisation for the *S2- solubility constant (all calibrated over wide range of silicate melt compositions).
        default: 'ONeill21dil' Eq. (10.34) inc. H2O dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        Other options:
        'ONeill21' Eq. (10.34) ex. H2O dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        'ONeill21hyd' (hydrous) Eq. (10.34, 10.49) from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        'Boulliung23eq6' Eq. (6) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9
        'Boulliung23eq7' Eq. (7) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9 
    
    sulfate: Model for the parameterisation of the S6+ solubility constant (all calibrated over wide range of silicate melt compositions).
        default: 'ONeill22dil' Eq. (12a) inc. H2O dilution from O'Neill & Mavrogenes (2022) GCA 334:368-382 10.1016/j.gca.2022.06.020
        Other options:
        'ONeill22' Eq. (12a) without H2O dilution from O'Neill & Mavrogenes (2022) GCA 334:368-382 doi:10.1016/j.gca.2022.06.020
        'Boulliung22nP' (no P-dependence) Eq. (5) from Boulliung & Wood (2023) GCA 343:420 doi:10.1016/j.gca.2022.11.025
        'Boulliung22wP' (inc. P-dependece) Eq. (5) from Boulliung & Wood (2023) GCA 343:420 doi:10.1016/j.gca.2022.11.025 and Eq. (8) for P from Boulliung & Wood (2022) GCA 336:150-164 doi:10.1016/j.gca.2022.08.032
        'Boulliung23eq9' Eq. (9) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9
        'Boulliung23eq11' Eq. (11) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9
    
    hydrogen sulfide: Model for the parameterisation for the H2S solubility constant.
        default 'Basalt_Hughes24' Fig.S6 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
        Other options:
        'BasalticAndesite_Hughes24' Fig.S6 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
    
    methane: Model for the parameterisation of the CH4 solubility constant.
        default: 'Basalt_Ardia13' Eq. (7a) from Ardia et al. (2013) GCA 114:52-71 doi:10.1016/j.gca.2013.03.028
        Only one option available currently, included for future development.

    carbon monoxide: Model for the parameterisation of the CO solubility constant.
        default: 'Basalt_Hughes24' CO in Table S4 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
        Only one option available currently, included for future development.

    species X solubility: Model for the parameterisation of the X solubility constant. 
        default: 'Ar_Basalt_HughesIP' Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        Other options:
        Ar_Rhyolite_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        Ne_Basalt_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        Ne_Rhyolite_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
    
    Cspeccomp: Model for the parameterisation of the speciation constant for CO2mol and CO32- in the melt.
        default: 'Basalt' Assume all oxidised carbon in the melt is present as carbonate ions.
        Other options:
        Andesite_Botcharnikov06: Eq. (8) from Botcharnikov et al. (2006) Chem.Geol. 229(1-3)125-143 doi:10.1016/j.chemgeo.2006.01.016
        Dacite_Botcharnikov06: Botcharnikov et al. (2006) Chem.Geol. 229(1-3)125-143 doi:10.1016/j.chemgeo.2006.01.016
        Rhyolite: Assume all oxidised carbon in the melt is present as molecular CO2.
    
    Hspeccomp: Model for the parameterisation of the speciation constant for H2Omol and OH- in the melt, either assuming ideal or regular solution models.
        default: 'MORB_HughesIP' [ideal solution only] Eq. SX in Hughes et al. (in prep) 
        Other options:
        MORB_Dixon95: [regular solution only] Table 5 of Dixon et al. (1995) JPet 36(6):1607-1631 doi:10.1093/oxfordjournals.petrology.a037267
        AlkaliBasalt_Lesne10: [regular solution only] Eq. (24-27) Lesne et al. (2010) CMP 162:133-151 doi:10.1007/s00410-010-0588-x
        StromboliAlkaliBasalt_Lesne10: Eq. (15) [ideal solution] or PST-9 in Table 5 [regular solution] from Lesne et al. (2010) CMP 162:133-151 doi:10.1007/s00410-010-0588-x
        VesuviusAlkaliBasalt_Lesne10: Eq. (16) [ideal solution] or VES-9 in Table 5 [regular solution] from Lesne et al. (2010) CMP 162:133-151 doi:10.1007/s00410-010-0588-x
        EtnaAlkaliBasalt_Lesne10: Eq. (17) [ideal solution] or ETN-1 in Table 5 [regular solution] from Lesne et al. (2010) CMP 162:133-151 doi:10.1007/s00410-010-0588-x
        Andesite_Botcharnikov06: [ideal solution only] Eq (7) from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143 doi:10.1016/j.chemgeo.2006.01.016
        Albite_Silver89: Fig. 8 [ideal solution only] or in the text [regular solution] from Silver & Stolper (1989) J.Pet 30(3)667-709 doi:10.1093/petrology/30.3.667
        Rhyolite_Zhang97: [ideal solution only] Eq. (9) from Zhang et al. (1997) GCA 61(15):3089-3100 doi:10.1016/S0016-7037(97)00151-8

            
    ### Saturation conditions ###

    SCSS: Model for parameterisation of the sulfide content at sulfide saturation (S2-CSS).
        default: 'ONeill21hyd' Eq. (10.34, 10.43, 10.45, 10.46, 10.49) from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        Other options:
        'ONeill21' Eq. (10.34, 10.43, 10.45, 10.46) excluding water dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        'ONeill21dil' Eq. (10.34, 10.43, 10.45, 10.46) including water dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        'Liu07' Eq. (9) in Liu et al. (2007) GCA 71:1783-1799 doi:10.1016/j.gca.2007.01.004
        'Fortin15' Eq. (7) Fortin et al. (2015) GCA 160:100-116 doi:10.1016/j.gca.2015.03.022
        'Liu21' Eq. (2) Liu et al. (2021) Chem.Geol. 559:119913 doi:10.1016.j.chemgeo.2020.119913
        'Fortin15_pss' Fortin et al. (2015) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        'Liu21_pss' Liu et al. (2021) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        'ONeill21_pss' O'Neill (2021) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        'ONeill22_pss' O'Neill & Mavrogenes (2022) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        'Smythe17_pss' Smythe et al. (2017) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    
    SCAS: Model for parameterisation of the sulfate content at anhydrite saturation (S6+CAS).
        default: 'Chowdhury19' Eq. (8) using Table 5 in Chowdhury & Dasgupta (2019) Chem.Geol. 522:162-174 doi:10.1016/j.chemgeo.2019.05.020
        Other options:
        'Zajacz19' Eq. (8-14) Zajacz & Tsay (2019) GCA 261:288-304 doi:10.1016/j.gca.2019.07.007
        'Liu23' Eq. (4) Liu et al. (2023) GCA 349:135-145 doi:10.1016/j.gca.2023.04.007
        'Zajacz19_pss' Zajacz and Tsay (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        'Chowdhury19_pss' Chowdhury & Dasgupta (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    
    sulfur_saturation: Is sulfur allowed to form sulfide or anhydrite if sulfur content of the melt reaches saturation levels for these phases?
        default: 'False' melt ± vapor are the only phases present - results are metastable with respect to sulfide and anhydrite if they could saturate.
        Other options:
        'True' If saturation conditions for sulfide or anhydrite are met, melt sulfur content reflects this.
    
    graphite_saturation: Is graphite allowed to form if the carbon content of the melt reaches saturation levels for graphite?
        default: 'False' melt ± vapor are the only phases present - results are metastable with respect to graphite if it could saturate.
        Other options:
        'True' If saturation conditions for graphite are met, melt carbon content reflects this.
           
    ### Fugacity coefficients ###
          
    ideal_gas: Treat all vapor species as ideal gases (i.e., all fugacity coefficients = 1 at all P).
        default: 'False' At least some of the vapor species are not treated as ideal gases. 
        Other options:
        'True' All fugacity coefficients = 1 at all P.
    
    y_CO2: Model for the parameterisation of the CO2 fugacity coefficient.
        default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        Other options:
        'Holland91' Holland & Powell (1991) CMP 109:265-273 10.1007/BF00306484
        'ideal' Treat CO2 as ideal gas species, fugacity coefficient = 1 at all P.

    y_SO2: Model for the parameterisation of the SO2 fugacity coefficient.
        default: 'Shi92_Hughes23' Fig.S1 Hughes et al. (2023) JGSL 180(3) doi:10.1144/jgs2021-12
        Other options:
        'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        'ideal' Treat SO2 as ideal gas species, fugacity coefficient = 1 at all P.

    y_H2S: Model for the parameterisation of the H2S fugacity coefficient.
        default: 'Shi92_Hughes24' Fig.S1 Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
        Other options:
        'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        'ideal' Treat H2S as ideal gas species, fugacity coefficient = 1 at all P.

    y_H2: Model for the parameterisation of the H2 fugacity coefficient.
        default: 'Shaw64' Eq. (4) from Shaw & Wones (1964) AmJSci 262:918-929
        Other options:
        'ideal' Treat H2 as ideal gas species, fugacity coefficient = 1 at all P.

    y_O2: Model for the parameterisation of the O2 fugacity coefficient.
        default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        Other options:
        'ideal' Treat O2 as ideal gas species, fugacity coefficient = 1 at all P.

    y_S2: Model for the parameterisation of the O2 fugacity coefficient.
        default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        Other options:
        'ideal' Treat S2 as ideal gas species, fugacity coefficient = 1 at all P.

    y_CO: Model for the parameterisation of the CO fugacity coefficient.
        default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        Other options:
        'ideal' Treat CO as ideal gas species, fugacity coefficient = 1 at all P.        

    y_CH4: Model for the parameterisation of the CH4 fugacity coefficient.
        default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        Other options:
        'ideal' Treat CH4 as ideal gas species, fugacity coefficient = 1 at all P.    

    y_H2O: Model for the parameterisation of the H2O fugacity coefficient.
        default: 'Holland91' Holland & Powell (1991) CMP 109:265-273 10.1007/BF00306484
        Other options:
        'ideal' Treat H2O as ideal gas species, fugacity coefficient = 1 at all P.    
    
    y_OCS: Model for the parameterisation of the OCS fugacity coefficient.
        default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        Other options:
        'ideal' Treat OCS as ideal gas species, fugacity coefficient = 1 at all P.            

    y_X: Model for the parameterisation of the X fugacity coefficient.
        default: 'ideal' Treat X as ideal gas species, fugacity coefficient = 1 at all P.
        Only one option available currently, included for future development.  
          
        
    ### Equilibrium constants ###
    
    KHOg: Model for the parameterisation of the equilibiurm constant for H2 + 0.5O2 = H2O.
        default: 'Ohmoto97' Reaction (d) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Only one option available currently, included for future development.

    KHOSg: Model for the parameterisation of the equilibiurm constant for 0.5S2 + H2O = H2S + 0.5O2.
        default: 'Ohmoto97' Reaction (h) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Other options:
        'no H2S' Stops H2S forming in the vapor (K = 0).
    
    KOSg: Model for the parameterisation of the equilibiurm constant for 0.5S2 + O2 = SO2.
        default: 'Ohmoto97' Reaction (f) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Other options:
        'no SO2' Stops SO2 forming in the vapor (K = 0). As a by-product, OCS will also stop forming.

    KOSg2: Model for the parameterisation of the equilibiurm constant for 0.5S2 + 1.5O2 = SO3.
        default: 'ONeill2' Eq. (6b) from O'Neill & Mavrogenes (2022) GCA 334:368-382 doi:10.1016/j.gca.2022.06.020
        Only one option available currently, included for future development.

    KOCg: Model for the parameterisation of the equilibiurm constant for CO + 0.5O2 = CO2.
        default: 'Ohmoto97' Reaction (c) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Only one option available currently, included for future development. 

    KCOHg: Model for the parameterisation of the equilibiurm constant for CH4 + 2O2 = CO2 + 2H2O.
        default: 'Ohmoto97' Reaction (e) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Other options:
        'no CH4' Stops CH4 forming in the vapor (K = large number).

    KOCSg: Model for the parameterisation of the equilibiurm constant for OCS.
        default: 'Moussallam19' Eq. (8) for 2CO2 + OCS ⇄ 3CO + SO2 in Moussallam et al. (2019) EPSL 520:260-267 doi:10.1016/j.epsl.2019.05.036 for 
        Other options:
        'no OCS' Stops OCS forming in the vapor (K = large number).  

    KCOs: Model for the parameterisation of the equilibiurm constant for Cgrahite + O2 = CO2.
        default: 'Holloway92' Eq. (3) KI in Holloway et al. (1992) EuropeanJ.Mineralogy 4(1):105-114.
        Only one option available currently, included for future development.

    carbonylsulfide: Reaction equilibrium KOCSg is for. 
        default: 'COS' 2CO2 + OCS ⇄ 3CO + SO2
        Only one option available currently, included for future development.

        
    ### Degassing calculation ###

    bulk_composition: Specifying what the inputted melt composition (i.e., dissolved volatiles and fO2-estimate) corresponds to for the degassing calculation
        default: 'melt-only' The inputted melt composition (i.e., dissolved volatiles) represents the bulk system - there is no vapor present. 
        The fO2-estimate is calculated at Pvsat for this melt composition.
        Other options:
        'melt+vapor_wtg' The inputted melt composition (i.e., dissolved volatiles) is in equilibrium with a vapor phase. 
        The amount of vapor as weight fraction gas (wtg) is specified in the inputs. The bulk system composition will be calculated by calculating Pvsat and the vapor composition given the input composition.
        'melt+vapor_initialCO2' The inputted melt composition (i.e., dissolved volatiles) is in equilibrium with a vapor phase. 
        The initial CO2 content of the melt (i.e., before degassing) is specified in the inputs. 
        The bulk system composition will be calculated by calculating Pvsat and the vapor composition given the input 
        composition. The amount of vapor present is calculated using the given initial CO2.
    
    starting_P: Determing the starting pressure for a degassing calculation.
        default: 'Pvsat' Calculation starts at Pvsat for the inputted melt composition (i.e., dissolved volatiles).
        Other options:
        'set' Calculation starts at the pressure specified in the inputs (using P_bar, pressure in bars).
    
    gassing_style: Does the bulk composition of the system (including oxygen) remain constant during the re/degassing 
    calculation.
        default: 'closed' The bulk composition of the system (inc. oxygen) is constant during re/degassing calculation - vapor and melt remain in chemical equilibrium throughout.
        Other options:
        'open' At each pressure-step, the vapor in equilibrium with the melt is removed (or added for regassing), such that the bulk composition of the system changes. This does not refer to being buffered in terms of fO2.
    
    gassing_direction: Is pressure increasing or decreasing from the starting perssure.
        default: 'degas' Pressure progressively decreases from starting pressure for isothermal, polybaric calculations (i.e., degassing).
        Other options:
        'regas' Pressure progressively increases from starting pressure for isothermal, polybaric calculations (i.e., regassing).
    
    P_variation: Is pressure varying during the calculation?
        default: 'polybaric' Pressure progressively changes during the calculation.
        Only one option available currently, included for future development.
    
    T_variation: Is temperature varying during the calculation?
        default: 'isothermal' Temperature is constant during the calculation.
        Only one option available currently, included for future development.
    
    # WHY DID I THINK THIS WAS USEFUL #
    eq_Fe: Does iron in the melt equilibrate with fO2.
        default: 'yes' Iron equilibrates with fO2
        Only one option available currently, included for future development.
    
    solve_species: What species are used to solve the equilibrium equations? 
    This should not need to be changed unless the solver is struggling.
        default: 'OCS' Guess mole fractions of O2, CO, and S2 in the vapor to solve the equilibrium equations.
        'OHS' Guess mole fractions of O2, H2, and S2 in the vapor to solve the equilibrium equations.
        'OCH' Guess mole fractions of O2, CO, and H2 in the vapor to solve the equilibrium equations.

        
    ### Other ###

    density: Model for parameterisation of melt density
        default: 'DensityX' DensityX from Iacovino & Till (2019) Volcanica 2(1):1-10 doi:10.30909/vol.02.01.0110
        Only one option available currently, included for future development.

    setup: Specifies whether model options are specified in the models or setup dataframe. 
        default: 'False' All model options are specified in the models dataframe.
        Other options:
        'True' Some of the model options are specified in the setup dataframe.
    
    print status: Specifies whether some sort of status information during the calculation is outputted to let you know progress.
        default: 'False' No information about calculation progress is printed.
        Other options:
        'True': Some information about calculation progress is printed.

    output csv: Specicies whether a csv of the outputted dataframe is saved at the end of the calculation.
        default: 'True' csv is outputted
        'False' csv is not outputted

    high precision: Is high preicision used for calculations?
        TRUE OR FALSE WHAT PRECISION IS IT USING
        default: 'False' normal precision used for calculations
        'True' high precision used

            
    ### In development ###

    For now, just leave them as their default option and everything should work fine!
          
    isotopes
        default: 'no'
        
    crystallisation
        default: 'no'
        
    mass_volume
        default: 'mass'
        
    calc_sat
        default: 'fO2_melt'
    
    bulk_O
        default: 'exc_S'
    
    error
        default: 0.1

    sulfur_is_sat
        default: 'no'
    """
    
    
    df_models = make_models_df(models)
    added_defaults = check_default_options(df_models)
    return added_defaults

#################################################################################################################################
##################################################### SOLUBILITY CONSTANTS ######################################################
#################################################################################################################################

###################################
### Solubility constant for H2O ###
###################################
def C_H2O(PT,melt_wf,models=default_models):
    
    """ 
    NOT COMPLETE

    Solubility constant for disolving H2O in the melt.

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions:
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
        
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "Hspeciation" and "water" and column label of "option"

    Model options
    -------------
    If "Hspeciation" = "none"
        default: 'Basalt_Hughes24' Fig.S2 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
        Other options:
        'Rhyolite_HughesIP' Fig.SX from Hughes et al. (in prep) based on data in Fig. 3 of Blank et al. (1993)
    

    Returns
    -------
    Solubility constant for H2O as <class 'mpfr'>

    """
    model_speciation = models.loc["Hspeciation","option"]
    model_solubility = models.loc["water","option"]

    if model_speciation in ["none",'none+ideal','none+regular']: ### C_H2O = (xmH2O)^2/fH2O ### (mole fraction)
        if model_solubility == "Rhyolite_HughesIP": # Fig.SX from Hughes et al. (in prep) based on data in Fig. 3 of Blank et al. (1993)
            C = 5.3851E-06 
        elif model_solubility == "Basalt_Hughes24": # Fig.S2 from Hughes et al. (2024) based on data compilation from Allison et al. (2022) for basalts with H2O < 6 wt%
            C = 4.6114e-6

        ### WORK IN PROGRESS ###
        elif model_solubility == "Lesne11mod": # modified general model from Lesne et al. (2011) 162:133-151
            C = 5.62316e-6
        elif model_solubility == "ETN-1" or model_solubility == "PST-9": # Fitted to ETN-1 and PST-9 from Lesne et al. (2011) 162:133-151
            C = 4.77591e-6
        elif model_solubility == "VES-9": # Fitted to VES-9 from Lesne et al. (2011) 162:133-151
            C = 5.46061e-6
        elif model_solubility ==  "evo": # Fitted to match EVo
            C = 2.782e-6
        elif model_solubility == "test": #test
            R_ = 83.144621 # cm3 bar K−1 mol−1
            DV = 12 # cm3/mol
            P0 = 1.0 # bar
            A = 4.6114e-6
            B = -((DV/(R_*T_K))*(P-P0))
            if models.loc["high precision","option"] == "True":
                C = A*gp.exp(B)
            else:
                C = A*math.exp(B)
        elif model_solubility == "test2": # for Ptot paper
            if models.loc["high precision","option"] == "True":
                C = gp.exp(-12.29)
            else:
                C = math.exp(-12.29)
        elif model_solubility == "carbon":
            C = 1.5e-9
        elif model_solubility == "test3":
            C = 6.22885E-09 # like 1000 ppm H2O at 730 bar
    
    elif model_speciation == "linear": ### C_H2O = xmH2O/fH2O ### (mole fraction)
        C = 0.00007925494 # like AllisonDataComp... I think.
            
    else: ### C_H2O = xmH2Omol/fH2O ### (mole fraction)
        P = PT['P']
        P0 = 1.0 # bar
        R = 83.15 # cm3 etc.
        T0 = 1473.15 # K
        if model_solubility == "Dixon95": # Dixon et al. (1995) - no compositional dependence
            DV = 12.
            A = 3.28e-5
            if models.loc["high precision","option"] == "True":
                C = A*gp.exp((-DV*(P-P0))/(R*T0))
            else:
                C = A*math.exp((-DV*(P-P0))/(R*T0))
        elif model_solubility == "alkali basalt": # Lesne et al. (2011) 162:133-151 eqn 31 with added RT term otherwise it will not work
            A = 5.71e-5 # XmH2Om0
            DV = 26.9 # VH2Om0 in cm3/mol
            if models.loc["high precision","option"] == "True":
                C = A*gp.exp((-DV*(P-P0))/(R*T0)) 
            else:
                C = A*math.exp((-DV*(P-P0))/(R*T0)) 
        
        ### Work in progress ###
        elif model_solubility == "ETN-1": # Fitted to ETN-1 and VES-9 Xm_H2Omol calculated at 1200 'C data from Lesne et al. (2011) 162:133-151
            C = 3.3989655e-6 
        elif model_solubility == "VES-9": # Fitted to ETN-1 and VES-9 Xm_H2Omol calculated at 1200 'C data from Lesne et al. (2011) 162:133-151
            C = 3.3989655e-6
        elif model_solubility == "PST-9": # Fitted to PST-9 Xm_H2Omol calculated at 1200 'C data from Lesne et al. (2011) 162:133-151
            C = 1.7022269e-6
    
    return C

          
##############################################
### Solubility constant for carbon dioxide ###
##############################################
def C_CO3(PT,melt_wf,models=default_models): ### C_CO2,T = xmCO2,T/fCO2 ### (mole fraction) ***except Shishkina14 - wmCO2 ppm***
    
    """ 
    Solubility constant for disolving CO2 in the melt: C_CO2,T = xmCO2,T/fCO2 [all oxidised carbon - i.e., CO2mol and CO32- - as CO2,T]


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "carbon dioxide" and column label of "option"

    Returns
    -------
    Solubility constant for CO2 as <class 'mpfr'>

    Model options
    -------------
    default: 'MORB_Dixon95' Bullet (5) of summary from Dixon et al. (1995) JPet 36(6):1607-1631 doi:10.1093/oxfordjournals.petrology.a037267
    Other options:
    'Basalt_Dixon97' Eq. (7) from Dixon et al. (1997) AmMin 82(3-4)368-378 doi:10.2138/am-1997-3-415
    'NorthArchBasalt_Dixon97' Eq. (8) from Dixon et al. (1997) AmMin 82(3-4)368-378 doi:10.2138/am-1997-3-415
    'Basalt_Lesne11' Eq. (25,26) from Lesne et al. (2011) CMP 162:153-168 doi:10.1007/s00410-010-0585-0
    'VesuviusAlkaliBasalt_Lesne11' VES-9 in Table 4 from Lesne et al. (2011) CMP 162:153-168 doi:10.1007/s00410-010-0585-0
    'EtnaAlkaliBasalt_Lesne11' ETN-1 in Table 4 from Lesne et al. (2011) CMP 162:153-168 doi:10.1007/s00410-010-0585-0
    'StromboliAlkaliBasalt_Lense11' PST-9 in Table 4 from Lesne et al. (2011) CMP 162:153-168 doi:10.1007/s00410-010-0585-0
    'SunsetCraterAlkaliBasalt_Allison19' Sunset Crater in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
    'SVFVBasalticAndesite_Allison19' SVFV in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
    'ErebusPhonotephrite_Allison19' Erebus in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
    'VesuviusPhonotephrite_Allison19' Vesuvius in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
    'EtnaTrachybasalt_Allison19' Etna in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
    'StromboliAlkaliBasalt_Allison19' Stromboli in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4 
    'Basanite_Holloway94' Basanite in Table 5 from Holloway and Blank (1994) RiMG 30:187-230 doi:10.1515/9781501509674-012
    'Leucitite_Thibault94' Leucitite from Thibault & Holloway (1994) CMP 116:216-224 doi:10.1007/BF00310701
    'TholeiiteBasalt_Allison22' N72 basalt in Table 2 from Allison et al. (2022) CMP 177:40 doi:10.1007/s00410-022-01903-y
    'Rhyolite_Blank93' Fig.2 caption from Blank et al. (1993) EPSL 119:27-36 doi:10.1016/0012-821X(93)90004-S

    """

    model = models.loc["carbon dioxide","option"]

    P = PT['P']
    T_K = PT['T']+273.15     
    # Calculate cation proportions with no volatiles but correct Fe speciation if available (a la Dixon 1997)
    melt_comp = mg.melt_cation_proportion(melt_wf,"no","no")

    R = 83.15
    T0 = 1473.15 # K
    PI = -6.5*(melt_comp["Si"]+melt_comp["Al"]) + 20.17*(melt_comp["Ca"] + 0.8*melt_comp["K"] + 0.7*melt_comp["Na"] + 0.4*melt_comp["Mg"] + 0.4*melt_comp["FeT"]) # Eq. (7) from Dixon (1997) Am. Min. 82:368-378
    PI_ = (melt_comp["Ca"] + 0.8*melt_comp["K"] + 0.7*melt_comp["Na"] + 0.4*melt_comp["Mg"] + 0.4*melt_comp["FeT"])/(melt_comp["Si"] + melt_comp["Al"]) # Shishkina et al. (2014) Chem. Geol. 388:112-129
 
    if model == "MORB_Dixon95": # Bullet (5) of summary from Dixon et al. (1995), which includes values from Pan et al. (1991)
        DV = 23. # cm3/mol
        P0 = 1.0 # bar
        A = 3.8e-7
        B = (-DV*(P-P0))/(R*T0)
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "Basalt_Dixon97": # Compositional dependence of Eq. (7) from Dixon (1997) Am. Min. 82:368-378 as shown in Eq. (1,5) from Witham et al. (2012)
        DV = 23 # cm3/mol
        P0 = 1.0 # bar
        A = (7.94e-7)*(PI+0.762)
        B = (-DV*(P-P0))/(R*T0)
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "NorthArchBasalt_Dixon97": # Eq. (8) from Dixon (1997) Am. Min. 82:368-378
        melt_comp_ox = mg.melt_normalise_wf(melt_wf,"no","no")
        DV = 23 # cm3/mol
        P0 = 1.0 # bar
        A = (8.70e-6)-((1.7e-7)*(melt_comp_ox["SiO2"]*100.))
        B = (-DV*(P-P0))/(R*T0)
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "Basalt_Lesne11": # Eq. (25, 26) from Lesne et al. (2011) based on Dixon (1997)
        DV = 25 # cm3/mol ±3
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(0.893*PI - 15.247) # Eq. (25)
        else:
            A = math.exp(0.893*PI - 15.247) # Eq. (25)
        B = (-DV*(P-P0))/(R*T0)
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "VesuviusAlkaliBasalt_Lesne11": # VES-9 in Table 4 of Lesne et al. (2011) CMP 162:153-168
        DV = 31.0 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.10) # ±0.03
        else:
            A = math.exp(-14.10) # ±0.03
        B = -((DV/(R*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "EtnaAlkaliBasalt_Lesne11": # ETN-1 in Table 4 of Lesne et al. (2011) CMP 162:153-168
        DV = 23.0 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.55) # ±0.00
        else:
            A = math.exp(-14.55) # ±0.00
        B = -((DV/(R*T_K))*(P-P0)) 
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "StromboliAlkaliBasalt_Lesne11": # PST-9 in Table 4 of Lesne et al. (2011) CMP 162:153-168
        DV = 6.0 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.74) # ±0.01
        else:
            A = math.exp(-14.74) # ±0.01
        B = -((DV/(R*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "Shishkina14": # modified from Shishkina et al. (2014) Chem. Geol. 388:112-129
        A = 1.164 # modified by converting P^A to APyCO2 but only including data up to and including 400 MPa
        B = 6.71*PI_-1.345
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "SunsetCraterAlkaliBasalt_Allison19": # Sunset Crater in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 16.40 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.67)
        else:
            A = math.exp(-14.67)
        B = -((DV/(R_*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "SFVFBasalticAndesite_Allison19": # SFVF in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 15.02 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.87)
        else:
            A = math.exp(-14.87)
        B = -((DV/(R_*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "ErebusPhonotephrite_Allison19": # Erebus in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = -14.65 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.65)
        else:
            A = math.exp(-14.65)
        B = -((DV/(R_*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "VesuviusPhonotephrite_Allison19": # Vesuvius in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 24.42 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.04)
        else:
            A = math.exp(-14.04)
        B = -((DV/(R_*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "EtnaTrachybasalt_Allison19": # Etna in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 21.59 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.28)
        else:
            A = math.exp(-14.28)
        B = -((DV/(R_*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "StromboliAlkaliBasalt_Allison19": # Stromboli in Table 4 from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 14.93 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.68)
        else:
            A = math.exp(-14.68)
        B = -((DV/(R_*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "Basanite_Holloway94": # Basanite in Table 5 from Holloway and Blank (1994) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 21.72 # cm3/mol ± 1.27
        DH = -13.1 # kJmol ± 13.9
        T0 = 1200. + 273.15 # K 
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.32)
        else:
            A = math.exp(-14.32)
        B = -((DV/(R_*T_K))*(P-P0)) + (DH/R)*((1.0/T0) - (1.0/T_K))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "Leucitite_Thibault94": # Leucitite from Thibault & Holloway (1994) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 21.53 # cm3/mol ± 0.42
        DH = -28.15 #kJ/mol ±4.24
        T0 = 1200. + 273.15 # K
        P0 = 1000.0 # bar
        A = gp.exp(-13.36)
        B = -((DV/(R_*T_K))*(P-P0)) + (DH/R)*((1.0/T0) - (1.0/T_K))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "TholeiiteBasalt_Allison22": # N72 basalt in Table 2 from Allison et al. (2022) CMP 177:40, based on experiments from Shishkina et al. (2010) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 19.05 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.86)
        B = -((DV/(R_*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "Rhyolite_Blank93": # Fig. 2 caption from Blank et al. (1993) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 28. # cm3/mol ± 2
        DH = -27.2 # kJ/mole ±2.1 
        P0 = 1.0 # bar
        T0 = 850. + 273.15 # K
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.45) # ±0.02
        else:
            A = math.exp(-14.45) # ±0.02
        B = -((DV/(R_*T_K))*(P-P0)) + (DH/R)*((1.0/T0) - (1.0/T_K))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    
    ### WORK IN PROGRESS BELOW HERE ###
    elif model == "Phonotephrite_Allison22": # AH3 Phonotephrite in Table 2 from Allison et al. (2022) CMP 177:40, based on experiments from Vetere et al. (2014) 
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 30.45 # cm3/mol
        P0 = 1000.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-13.26)
        else:
            A = math.exp(-13.26)
        B = -((DV/(R_*T_K))*(P-P0))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "Allison22mod": # modified from Allison et al. (2022) CMP 177:40
        P0 = 1000. # bars
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = -3350.650 + 3375.552*(Si+Na) + 2625.385*Ti + 3105.426*Al + 3628.018*Fe2 + 3323.320*(Mg+Ca) + 3795.115*K + 47.004*(Na/(Na+K)) # cm/mol
        lnK0 = -128.365 + 114.098*Si + 92.263*(Ti+Al) + 122.644*(Fe2+Ca+Na) + 111.549*Mg + 138.855*K + 2.239*(Na/(Na+K))
        if models.loc["high precision","option"] == "True":
            A = gp.exp(lnK0)
        else:
            A = math.exp(lnK0)
        B = ((-1.*DV)*(P-P0))/(R_*T_K)
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "scaledCsulfate": # O'Neill & Mavrogenes (2022) GCA 334:368-382 eq[12a]
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        lnC = -8.02 + ((21100. + 44000.*Na + 18700.*Mg + 4300.*Al + 44200.*K + 35600.*Ca + 12600.*Mn + 16500.*FeT)/T_K) #CS6+ = [S6+, ppm]/fSO3 
        Csulfate = gp.exp(lnC)*KOSg2(PT,models) # ppm S
        lnCsulfate = math.log(Csulfate)
        lnC=-0.46*lnCsulfate
        if models.loc["high precision","option"] == "True":
            A = gp.exp(lnC)
        else:
            A = math.exp(lnC)
        DV = 23.
        P0 = 1000. # bars
        R_ = 83.144621 # cm3 bar K−1 mol−1
        B = ((-1.*DV)*(P-P0))/(R_*T_K)
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "Behrens04fit": # Fit to Behrens et al. (2004) - tried for workshop
        DV = 41.8 # cm3/mol
        P0 = 1.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.2)
        else:
            A = math.exp(-14.2)
        B = (-DV*(P-P0))/(R*(1250.+273.15))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "dacite": # Fit to Behrens et al. (2004) using Ptot
        DV = 36.5 # cm3/mol
        P0 = 1.0 # bar
        if models.loc["high precision","option"] == "True":
            A = gp.exp(-14.3)
        else:
            A = math.exp(-14.3)
        B = (-DV*(P-P0))/(R*(1250.+273.15))
        if models.loc["high precision","option"] == "True":
            C = A*gp.exp(B)
        else:
            C = A*math.exp(B)
    elif model == "test": # for Ptot paper!!!
        P0 = 1.0 # bar
        # "CO2mol"
        #DV1 = 24.2340838328 # cm3/mol
        #A1 = gp.exp(-14.92978383)
        #B1 = (-DV1*(P-P0))/(R*T_K)
        # "CO32-"
        #DV2 = 8.912862511 # cm3/mol
        #A2 =gp.exp(])
        #B2 = (-DV2*(P-P0))/(R*T_K)
        #C = A1*gp.exp(B1) + A2*gp.exp(B2)
        # "average"
        DV2 = 16.57 # cm3/mol
        if models.loc["high precision","option"] == "True":
            A2 =gp.exp(-15.275)
        else:
            A2 =math.exp(-15.275)
        B2 = (-DV2*(P-P0))/(R*T_K)
        if models.loc["high precision","option"] == "True":
            C = A2*gp.exp(B2)
        else:
            C = A2*math.exp(B2)
    elif model == "water":
        C = 27.e-6
    return C


########################################
### solubility constant for sulfide ###
########################################
def C_S(PT,melt_wf,models=default_models): ### C_S = wmS2-*(fO2/fS2)^0.5 ### (weight ppm)
    
    """ 
    Solubility constant for disolving S in the melt as S2- in ppmw: C_S = wmS2-*(fO2/fS2)^0.5


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "sulfide" and column label of "option"
    
    Model options
    -------------
    default: 'ONeill21dil' Eq. (10.34) inc. H2O dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    Other options:
    'ONeill21' Eq. (10.34) ex. H2O dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    'ONeill21hyd' (hydrous) Eq. (10.34, 10.49) from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    'Boulliung23eq6' Eq. (6) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9
    'Boulliung23eq7' Eq. (7) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9
    
    Returns
    -------
    Solubility constant for S2- as <class 'mpfr'>

    """
    model = models.loc["sulfide","option"]
    
    T = PT['T'] + 273.15 # T in K
    
    def ONeill21(T,melt_comp): # Eq. (10.34) in O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        lnC = (8.77 - (23590.0/T) + (1673.0/T)*(6.7*(melt_comp["Na"]+melt_comp["K"]) + 4.9*melt_comp["Mg"] + 8.1*melt_comp["Ca"] + 8.9*(melt_comp["FeT"]+melt_comp["Mn"]) + 5.0*melt_comp["Ti"] + 1.8*melt_comp["Al"] - 22.2*melt_comp["Ti"]*(melt_comp["FeT"]+melt_comp["Mn"]) + 7.2*melt_comp["FeT"]*melt_comp["Si"]) - 2.06*math.erf(-7.2*(melt_comp["FeT"]+melt_comp["Mn"])))
        return lnC
    
    if model == "ONeill21": # Eq. (10.34) in O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","no")
        lnC = ONeill21(T,melt_comp)
        C = math.exp(lnC) 
     
    if model == "ONeill21dil": # Eq. (10.34) O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf,"water","no")        
        lnC = ONeill21(T,melt_comp)
        C = math.exp(lnC) 

    if model == "ONeill21hyd": # Eq. (10.34, 10.49) in O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","no")
        lnC_nondil = ONeill21(T,melt_comp)
        melt_comp = mg.melt_cation_proportion(melt_wf,"water","no")
        lnCdil = ONeill21(T,melt_comp)
        lnCH = (melt_comp["H"]*(6.4 + 12.4*melt_comp["H"] - 20.3*melt_comp["Si"] + 73.0*(melt_comp["Na"]+melt_comp["K"]))) # Eq. (10.49)
        lnC = lnCdil+lnCH
        C = math.exp(lnC) 

    if model == "Boulliung23eq6": # Eq. (6) from Boulliung, J., Wood, B.J. Sulfur oxidation state and solubility in silicate melts. Contrib Mineral Petrol 178, 56 (2023). https://doi.org/10.1007/s00410-023-02033-9
        # Mole fractions in the melt on cationic lattice with no volatiles and Fe speciated
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","yes")
        logC = 0.338 + (24328.*melt_comp["Fe2"] + 5411.*melt_comp["Ca"] + 15872.*melt_comp["Mn"] - 9697.)/T
        C = 10.**(logC) 
    
    if model == "Boulliung23eq7": # Eq. (7) from Boulliung, J., Wood, B.J. Sulfur oxidation state and solubility in silicate melts. Contrib Mineral Petrol 178, 56 (2023). https://doi.org/10.1007/s00410-023-02033-9
        # Mole fractions in the melt on cationic lattice with no volatiles and Fe speciated
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","yes")
        logC = 0.225 + (25237.*melt_comp["Fe2"] + 5214.*melt_comp["Ca"] + 12705.*melt_comp["Mn"] + 19829.*melt_comp["K"] - 1109.*melt_comp["Si"] - 8879.)/T
        C = 10.**(logC)

    #elif model == "FR54-S1":
    #    lnC = math.log(((1.3e-4)*10000.))

    return C



########################################
### solubility constant for sulfate ###
########################################
def C_SO4(PT,melt_wf,models=default_models): ### C_SO4 = wmS6+*(fS2*fO2^3)^-0.5 ### (weight ppm)
    
    """ 
    Solubility constant for disolving S6+ in the melt: C_SO4 = wmS6+(fS2*fO2^3)^-0.5 (in ppmw and bar)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "sulfate" and column label of "option"

    Returns
    -------
    Solubility constant for S6+ as <class 'mpfr'>

    Model options
    -------------
    default: 'ONeill22dil' Eq. (12a) inc. H2O dilution from O'Neill & Mavrogenes (2022) GCA 334:368-382 10.1016/j.gca.2022.06.020
    Other options:
    'ONeill22' Eq. (12a) without H2O dilution from O'Neill & Mavrogenes (2022) GCA 334:368-382 doi:10.1016/j.gca.2022.06.020
    'Boulliung22nP' Eq. (5) from Boulliung & Wood (2023) GCA 343:420 doi:10.1016/j.gca.2022.11.025
    'Boulliung22wP' Eq. (5) from Boulliung & Wood (2023) GCA 343:420 doi:10.1016/j.gca.2022.11.025 and Eq. (8) for P from Boulliung & Wood (2022) GCA 336:150-164 doi:10.1016/j.gca.2022.08.032
    'Boulliung23eq9' Eq. (9) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9
    'Boulliung23eq11' Eq. (11) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9

    """

    model = models.loc["sulfate","option"]

    T = PT['T'] + 273.15 # T in Kelvin
    P = PT['P'] # P in bars
    #slope = 115619.707 # slope for T-dependence for melt inclusion fits
    
    if model in ["Boulliung22nP","Boulliung22wP"]: # Eq. (5) in Boullioung & Wood (2022) GCA 336:150-164 - corrected!
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","no")
        logCS6 = -12.948 + ((15602.*melt_comp["Ca"] + 28649.*melt_comp["Na"] - 9596.*melt_comp["Mg"] + 4194.*melt_comp["Al"] +16016.*melt_comp["Mn"] + 29244.)/T) # wt% S
        if model == "Boulliung22wP": 
            logCS6 = logCS6 - ((0.1*((10.*P)-0.1))*1.5237)/T # wt% S, Eq. (8) Boullioung & Wood (2022) GCA 336:150-164
        Csulfate = (10.**logCS6)*10000. # ppm S
    elif model in ["ONeill22","ONeill22dil"]: 
        if model == "ONeill22": # Eq. (12a) in O'Neill & Mavrogenes (2022) GCA 334:368-382
            melt_comp = mg.melt_cation_proportion(melt_wf,"no","yes") # Mole fractions in the melt on cationic lattice (Fe as Fe2 and Fe3) no volatiles   
        elif model == "ONeill22dil": # Eq. (12a) in O'Neill & Mavrogenes (2022) GCA 334:368-382
            melt_comp = mg.melt_cation_proportion(melt_wf,"water","yes") # Mole fractions in the melt on cationic lattice (Fe as Fe2 and Fe3) includes water
        lnC = -8.02 + ((21100. + 44000.*melt_comp["Na"] + 18700.*melt_comp["Mg"] + 4300.*melt_comp["Al"] + 44200.*melt_comp["K"] + 35600.*melt_comp["Ca"] + 12600.*melt_comp["Mn"] + 16500.*melt_comp["Fe2"])/T) #CS6+ = [S6+, ppm]/fSO3
        if models.loc["high precision","option"] == "True":
            Csulfate = gp.exp(lnC)*KOSg2(PT,models) # ppm S
        else:
            Csulfate = math.exp(lnC)*KOSg2(PT,models) # ppm S
    elif model == "Boulliung23eq9": # Eq. (9) from Boulliung, J., Wood, B.J. Sulfur oxidation state and solubility in silicate melts. Contrib Mineral Petrol 178, 56 (2023). https://doi.org/10.1007/s00410-023-02033-9
        # Mole fractions in the melt on cationic lattice with no volatiles and Fe speciated
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","yes")
        logC = -12.948 + (28649.*melt_comp["Na"] + 15602.*melt_comp["Ca"] + 9496.*melt_comp["Mg"] + 16016.*melt_comp["Mn"] + 4194.*melt_comp["Al"] + 29244.)/T
        Csulfate = 10.**(logC)
    elif model == "Boulliung23eq11": # Eq. (11) from Boulliung, J., Wood, B.J. Sulfur oxidation state and solubility in silicate melts. Contrib Mineral Petrol 178, 56 (2023). https://doi.org/10.1007/s00410-023-02033-9
        # Mole fractions in the melt on cationic lattice with no volatiles and Fe speciated
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","yes")
        logC = -213.65 + (25696.*melt_comp["Na"] + 15076.*melt_comp["Ca"] + 9543.*melt_comp["Mg"] + 16158.*melt_comp["Mn"] + 4316.*melt_comp["Al"] + 68254.)/T + 55.03*math.log10(T)
        Csulfate = 10.**(logC)
    
    ### OLD ###
    #elif model == "Nash19": # Nash et al. (2019) EPSL 507:187-198
    #    S = 1. # S6+/S2- ratio of S6+/S2- of 0.5
    #    Csulfide = C_S(PT,melt_wf,models)
    #    A = PT_KCterm(PT,melt_wf,models) # P, T, compositional term from Kress & Carmicheal (1991)
    #    B = (8743600/T**2) - (27703/T) + 20.273 # temperature dependence from Nash et al. (2019)
    #    a = 0.196 # alnfO2 from Kress & Carmicheal (1991)
    #    F = 10**(((math.log10(S))-B)/8.)
    #    fO2 = math.exp(((math.log(0.5*F))-A)/a)
    #    Csulfate = (S*Csulfide)/(fO2**2)
    #elif model == "S6ST":
    #    Csulfide = C_S(PT,melt_wf,models)
    #    fO2 = f_O2(PT,melt_wf,models)
    #    S6ST_ = melt_wf["S6ST"]
    #    S = overtotal2ratio(S6ST_)
    #    Csulfate = (S*Csulfide)/(fO2**2)
    #elif model == "Hawaii":
    #    Csulfate = gp.exp(30.4) # Using Brounce et al. (2017) dataset at 1200 'C
    #    Csulfate = math.exp(slope*(1./T) -48.)
    #elif model == "Etna":
    #    Csulfate = math.exp(slope*(1./T) -50.15)
    #elif model == "Fuego":
    #    Csulfate = math.exp(slope*(1./T) -48.5)
    #elif model == "Erta Ale":
    #    Csulfate = math.exp(slope*(1./T) -45.5)
    #elif model == "FR54-S1":
    #    Csulfate = ((67.e6)*10000.)
    #elif model == "JdF": # 1100 'C ONLY
    #    Csulfate = 10.**17.
           
    return Csulfate



###################################
### solubility constant for H2S ###
###################################
def C_H2S(PT,melt_wf,models=default_models): # C_H2S = wmH2S/fH2S (ppm H2S, fH2S bar)
    
    """ 
    Solubility constant for disolving H2S in the melt: C_H2S = wmH2S/fH2S (ppmw/bar)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "hydrogen sulfide" and column label of "option"

    Returns
    -------
    Solubility constant for H2S as <class 'mpfr'>

    Model options
    -------------
    default: 'Basalt_Hughes24' Fig.S6 from Hughes et al. (2024) based on experimental data Moune et al. (2009) and calculations in Lesne et al. (2011)
    Other options:
    'BasalticAndesite_Hughes24' Fig.S6 from Hughes et al. (2024) based on experimental data Moune et al. (2009) and calculations in Lesne et al. (2011)

    """

    model = models.loc["hydrogen sulfide","option"]
    if model == "Basalt_Hughes24":
        K = 10.23 # fitted to basalt data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015 ChemGeol 418:104–116
    elif model == "BasalticAndesite_Hughes24":
        K = 6.82 # fitted to basaltic andesite data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015 ChemGeol 418:104–116 
    return K



########################################
### solubility constant for hydrogen ###
########################################
def C_H2(PT,melt_wf,models=default_models): # C_H2 = wmH2/fH2 (wtppm)
    
    """ 
    Solubility constant for disolving H2 in the melt: C_H2 = wmH2/fH2 (ppmw/bar)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "hydrogen" and column label of "option"

    Returns
    -------
    Solubility constant for H2 as <class 'mpfr'>

    Model options
    -------------
    default: 'Basalt_Hughes24' Basalt in Table S4 from Hughes et al. (2024) based on experimetnal data from Hirschmann et al. (2012)
    Other options:
    'Andesite_Hughes24' Andesite in Table S4 from Hughes et al. (2024) based on experimental data from Hirschmann et al. (2012)

    """

    model = models.loc["hydrogen","option"] 
    
    # Hirchmann et al. (2012) EPSL 345-348:38-48
    R = 83.144598 # bar cm3 /mol /K
    P = PT['P'] # pressure in bars
    T = PT['T'] + 273.15 # T in Kelvin SHOULD BE T0
    P0 = 100.*0.01 # kPa to bars
    if model == "Basalt_Hughes24": # Basalt in Table S4 from Hughes et al. (2024) based on experimetnal data from Hirschmann et al. (2012)
        #lnK0 = -11.4 # T0 = 1400 'C, P0 = 100 kPa for mole fraction H2
        lnK0 = -0.9624 # for ppm H2 (fitted in excel)
        DV = 10.6 # cm3/mol
    elif model == "Basalt_Hughes24": # Andesite in Table S4 from Hughes et al. (2024) based on experimental data from Hirschmann et al. (2012)
        #lnK0 = -10.6 # T0 = 1400 'C, P0 = 100 kPa for mole fraction H2
        lnK0 = -0.1296 # for ppm H2 (fitted in excel)
        DV = 11.3 # cm3/mol
    lnK = lnK0 - (DV*(P-P0))/(R*T) # = ln(XH2/fH2) in ppm/bar
    if models.loc["high precision","option"] == "True":
        C = gp.exp(lnK) 
    else:
        C = math.exp(lnK) 
    return C



######################################
### solubility constant for methane ##
######################################
def C_CH4(PT,melt_wf,models=default_models): # C_CH4 = wmCH4/fCH4 (ppm)
    
    """ 
    Solubility constant for disolving CH4 in the melt: C_CH4 = wmCH4/fCH4 (ppmw/bar)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "methane" and column label of "option"

    Returns
    -------
    Solubility constant for CH4 as <class 'mpfr'>

    Model options
    -------------
    default: 'Basalt_Ardia13' Eq. (7a) from Ardia et al. (2013)
    Only one option available currently, included for future development.

    """

    model = models.loc["methane","option"]

    if model == "Basalt_Ardia13": # Eq. (7a) from Ardia et al. (2013) GCA 114:52-71
        R = 83.144598 # bar cm3 /mol /K 
        P = PT['P'] # pressure in bars
        T = PT['T'] + 273.15 # T in Kelvin SHOULD BE T0
        P0 = 100.*0.01 # kPa to bars
        lnK0 = 4.93 # ppm CH4 
        #lnK0 = -7.63 # mole fraction CH4
        DV = 26.85 # cm3/mol
        lnK = lnK0 - (DV*(P-P0))/(R*T) 
        if models.loc["high precision","option"] == "True":
            K_ = gp.exp(lnK) # for fCH4 in GPa
        else:
            K_ = math.exp(lnK) # for fCH4 in GPa
        K = 0.0001*K_ # for fCH4 in bars 
    return K



#################################
### solubility constant for CO ##
#################################
def C_CO(PT,melt_wf,models=default_models): # C_CO = wmCO/fCO (ppm)
    
    """ 
    Solubility constant for disolving CO in the melt: C_CO = wmCO/fCO (ppmw/bar)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "carbon monoxide" and column label of "option"

    Returns
    -------
    Solubility constant for CH4 as <class 'mpfr'>

    Model options
    -------------
    default: 'Basalt_Hughes24' CO in Table S4 from Hughes et al. (2024) based on data from Armstrong et al. (2015), Stanley et al., (2014), and Wetzel et al., (2013)
    Only one option available currently, included for future development.

    """

    model = models.loc["carbon monoxide","option"]

    if model == "Basalt_Hughes24": # from fitting Armstrong et al. (2015) GCA 171:283-302; Stanley+2014, and Wetzel+13 thermodynamically
        R = 83.144598 # bar cm3 /mol /K 
        P = PT['P'] # pressure in bars
        T = PT['T'] +273.15 # T in Kelvin
        P0 = 1. # in bars
        lnK0 = -2.11 # ppm CO
        DV = 15.20 # cm3/mol
        lnK = lnK0 - (DV*(P-P0))/(R*T) 
        if models.loc["high precision","option"] == "True":
            K = gp.exp(lnK) # CO(ppm)/fCO(bars)
        else:
            K = math.exp(lnK) # CO(ppm)/fCO(bars)
    return K


#################################
### solubility constant for X ###
#################################
def C_X(PT,melt_wf,models=default_models): # C_X = wmX/fX (ppm)
    
    """ 
    Solubility constant for disolving X in the melt: C_X = wmX/fX (ppmw/bar)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "species X solubility" and column label of "option"

    Returns
    -------
    Solubility constant for X as <class 'mpfr'>

    Model options
    -------------
    default: 'Ar_Basalt_HughesIP' Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
    Other options:
    Ar_Rhyolite_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
    Ne_Basalt_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
    Ne_Rhyolite_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157

    """

    model = models.loc["species X solubility","option"]
        
    if model == "Ar_Basalt_HughesIP": # Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        K = 0.0799 # fitted assuming Ar is an ideal gas... i.e. yAr = 1.
    elif model == "Ar_Rhyolite_HughesIP": # Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        K = 0.4400 # fitted assuming Ar is an ideal gas... i.e. yAr = 1.
    elif model == "Ne_Basalt_HughesIP": # Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        K = 0.1504 # fitted assuming Ne is an ideal gas... i.e. yNe = 1.
    elif model == "Ne_Rhyolite_HughesIP": # Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        K = 0.8464 # fitted assuming Ne is an ideal gas... i.e. yNe = 1.
    
    ### WORK IN PROGRESS ###
    elif model == "test": 
        #K = 40. # similar to H2O
        #K = 6. # similar to S @ DFMQ+1.25
        #K = 21. # similar to S @ DFMQ+3
        #K = 155 # similar to S @ DFMQ0
        #K = 918005 # similar to S @DFMQ-3
        #K = 10.23 # similar to H2S
        #K = 0.51 # similar to CO32-
        #K = 1.37 # degassed at a similar depth to H2OT at 3wt%
        #K = 100.
        K = 35.

    return K

########################
### solubility of Cl ###    
########################
def C_Cl(PT,melt_wf):
    melt_comp = vf.melt_cation_proportion(melt_wf,"no","no")
    P = PT['P']/10000. # bar to GPa
    T = PT['T'] + 273.15 # 'C to 'K
    logC = 1.601 + (4470*melt_comp['Ca'] - 3430*melt_comp['Si'] + 2592*melt_comp['FeT'] - 4092*melt_comp['K'] - 894*P)/T
    C = float(10.**logC)
    return C

#################################################################################################################################
################################################ solid/liquid volatile saturation ###############################################
#################################################################################################################################

################################################
### sulfate content at anhydrite saturation ###
################################################
def SCAS(PT,melt_wf,models=default_models): 
    """ 
    Sulfate content at anhydrite saturation (S6+CAS)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "SCAS" and column label of "option"

    Returns
    -------
    S6+CAS in ppm as <class 'mpfr'>

    Model options
    -------------
    default: 'Zajacz19' Eq. (8-14) Zajacz & Tsay (2019) GCA 261:288-304 doi:10.1016/j.gca.2019.07.007
    Other options:
    'Chowdhury19' Eq. (8) using Table 5 in Chowdhury & Dasgupta (2019) Chem.Geol. 522:162-174 doi:10.1016/j.chemgeo.2019.05.020
    'Liu23' Eq. (4) Liu et al. (2023) GCA 349:135-145 doi:10.1016/j.gca.2023.04.007
    'Chowdhury19_pss' Chowdhury & Dasgupta (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    'Zajacz19_pss' Zajacz and Tsay (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    """
    model = models.loc["SCAS","option"]
    
    T = PT['T'] +273.15
    
    comp = mg.melt_pysulfsat(melt_wf)
 
    if model == "Chowdhury19": # Eq. (8) using Table 5 in Chowdhury & Dasgupta (2019) Chem.Geol. 522:162-174 doi:10.1016/j.chemgeo.2019.05.020
        # sulfate content (ppm) at anhydrite saturation [T in K]
        # mole fraction melt composition including water but all Fe as FeOT
        melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no") # different molecular weights used to original paper
        tot = 100.*melt_comp["mol_tot"]
        a = -13.23
        b = -0.50
        dSi = 3.02
        dCa = 36.70
        dMg = 2.84
        dFe = 10.14
        dAl = 44.28
        dNa = 26.27
        dK = -25.77
        e = 0.09
        f = 0.54
        wm_H2OT = 100.*melt_wf['H2OT']
        dX = dSi*melt_comp["SiO2"] + dCa*melt_comp["CaO"] + dMg*melt_comp["MgO"] + dFe*melt_comp["FeOT"] + dAl*melt_comp["Al2O3"] + dNa*melt_comp["Na2O"] + dK*melt_comp["K2O"]
        if models.loc["high precision","option"] == "True":
            lnxm_SO4 = a + b*((10.0**4.0)/T) + dX + e*wm_H2OT - f*gp.log(melt_comp["CaO"])                                                                                  
            xm_SO4 = gp.exp(lnxm_SO4)
        else:
            lnxm_SO4 = a + b*((10.0**4.0)/T) + dX + e*wm_H2OT - f*math.log(melt_comp["CaO"])                                                                                  
            xm_SO4 = math.exp(lnxm_SO4)
        Xm_SO4 = xm_SO4*(xm_SO4 + tot)
        S6CAS = Xm_SO4*species.loc["S","M"]*10000.0
                           
    elif model == "Zajacz19": # Eq. (8-14) Zajacz & Tsay (2019) GCA 261:288-304 doi:10.1016/j.gca.2019.07.007
        # mole fraction melt composition including water but all Fe as FeOT
        melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no") # different molecular weights used to original paper
        tot = 100.*melt_comp["mol_tot"]
        if melt_comp["Na2O"] + melt_comp["K2O"] + melt_comp["CaO"] >= melt_comp["Al2O3"]:
            P_Rhyo = 3.11*(melt_comp["Na2O"]+melt_comp["K2O"]+melt_comp["CaO"]-melt_comp["Al2O3"])
        else:
            P_Rhyo = 1.54*(melt_comp["Al2O3"]-(melt_comp["Na2O"]+melt_comp["K2O"]+melt_comp["CaO"]))
        NBOT = (2.*melt_comp["Na2O"]+2.*melt_comp["K2O"]+2.*(melt_comp["CaO"]+melt_comp["MgO"]+melt_comp["FeOT"])-melt_comp["Al2O3"]*2.)/(melt_comp["SiO2"]+2.*melt_comp["Al2O3"]) # according to spreadsheet not paper
        P_H2O = melt_comp["H2O"]*(2.09 - 1.65*NBOT) + 0.42*NBOT + 0.23
        P_C = ((P_Rhyo + 251.*melt_comp["CaO"]**2. + 57.*melt_comp["MgO"]**2. + 154.*melt_comp["FeOT"]**2.)/(2.*melt_comp["Al2O3"] + melt_comp["SiO2"]))/(1. + 4.8*NBOT)
        if models.loc["high precision","option"] == "True":
            P_T = gp.exp(-7890./T)
            Ksm_SPAnh = gp.exp(1.226*gp.log(P_C*P_T*P_H2O) + 0.079)         
        else:
            P_T = math.exp(-7890./T)
            Ksm_SPAnh = math.exp(1.226*math.log(P_C*P_T*P_H2O) + 0.079)                    
        Xsm_S = Ksm_SPAnh/melt_comp["CaO"]  
        S6CAS = Xsm_S*tot*32.07*10000.
        
    elif model == "Liu23": # Eq. (4) Liu et al. (2023) GCA 349:135-145 doi:10.1016/j.gca.2023.04.007
        melt_comp = mg.melt_mole_fraction(melt_wf,models,"no","no")
        NBOT = (2.*melt_comp["Na2O"]+2.*melt_comp["K2O"]+2.*(melt_comp["CaO"]+melt_comp["MgO"]+melt_comp["FeOT"])-melt_comp["Al2O3"]*2.)/(melt_comp["SiO2"]+2.*melt_comp["Al2O3"]) 
        melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no")
        lnSCAS = 12.185 - (8541./T) + (1.409*NBOT) + 9.984*melt_comp["CaO"] + melt_wf["H2OT"]*100.
        if models.loc["high precision","option"] == "True":
            S6CAS = gp.exp(lnSCAS)
        else:
            S6CAS = math.exp(lnSCAS)        
    
    elif model == "Chowdhury19_pss": # Chowdhury & Dasgupta (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_CD2019_SCAS(df=comp, T_K=T, H2O_Liq=float(100.*melt_wf['H2OT']), Fe3Fet_Liq=None, P_kbar=None)
        S6CAS = float(output["SCAS6_ppm"].iloc[0])
    elif model == "Zajacz19_pss": # Zajacz and Tsay (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_ZT2022_SCAS(df=comp, T_K=T, H2O_Liq=float(100.*melt_wf['H2OT']), Fe3Fet_Liq=None, P_kbar=None)
        S6CAS = float(output["SCAS6_ppm"].iloc[0])
    elif model == "Masotta15_pss": # Masotta and Kepler (2015) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_MK2015_SCAS(df=comp, T_K=T, H2O_Liq=float(100.*melt_wf['H2OT']), Fe3Fet_Liq=None, P_kbar=None)
        S6CAS = float(output["SCAS6_ppm"].iloc[0])
    
    return S6CAS

###############################################
### sulfide content at sulfide saturation ###
###############################################
def SCSS(PT,melt_wf,models=default_models): # sulfide content (ppm) at sulfide saturation from O'Neill (2020) [P in bar,T in K]
    """ 
    Sulfide content at sulfide saturation (S2-CSS)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "SCSS" and column label of "option"

    Returns
    -------
    S2-CCS in ppm as <class 'mpfr'>

    Model options
    -------------
    default: 'ONeill21hyd' Eq. (10.34, 10.43, 10.45, 10.46, 10.49) from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    Other options:
    'ONeill21' Eq. (10.34, 10.43, 10.45, 10.46) excluding water dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    'ONeill21dil' Eq. (10.34, 10.43, 10.45, 10.46) including water dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    'Liu07' Eq. (9) in Liu et al. (2007) GCA 71:1783-1799 doi:10.1016/j.gca.2007.01.004
    'Fortin15' Eq. (7) Fortin et al. (2015) GCA 160:100-116 doi:10.1016/j.gca.2015.03.022
    'Liu21' Eq. (2) Liu et al. (2021) Chem.Geol. 559:119913 doi:10.1016.j.chemgeo.2020.119913
    'Fortin15_pss' Fortin et al. (2015) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    'Liu21_pss' Liu et al. (2021) assuming pure FeS using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    'ONeill22_pss' O'Neill & Mavrogenes (2022) assuming pure FeS using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    'ONeill21_pss' O'Neill (2021) assuming pure FeS using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    'Smythe17_pss' Smythe et al. (2017) assuming pure FeS using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127

    """
    
    model = models.loc["SCSS","option"]

    P_bar = PT['P']
    T = PT['T'] + 273.15
    Fe3FeT = melt_wf["Fe3FeT"]
    comp = mg.melt_pysulfsat(melt_wf)

    if "sulf_XFe" in melt_wf:
        sulf_XFe = melt_wf["sulf_XFe"]
    else:
        sulf_XFe = 1.
    if "sulf_XCu" in melt_wf:
        sulf_XCu = melt_wf["sulf_XCu"]       
    else:
        sulf_XCu = 0.
    if "sulf_XNi" in melt_wf:
        sulf_XNi = melt_wf["sulf_XNi"]       
    else:
        sulf_XNi = None

    if model in ["ONeill21","ONeill21dil","ONeill21hyd"]: # O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        R = 8.31441
        P = (1.0e-4)*P_bar # pressure in GPa
        if models.loc["high precision","option"] == "True":
            D = (137778.0 - 91.666*T + 8.474*T*gp.log(T)) # J/mol Eq. (10.45)
        else:
            D = (137778.0 - 91.666*T + 8.474*T*math.log(T)) # J/mol Eq. (10.45)
        sulfide_comp = sulf_XFe
        if model == "ONeill21": # Eq. (10.34, 10.43, 10.45, 10.46) 
            # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) no volatiles
            melt_comp = mg.melt_cation_proportion(melt_wf,"no","yes")
        elif model == "ONeill21dil": # Eq. (10.34, 10.43, 10.45, 10.46)
            # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) and water
            melt_comp = mg.melt_cation_proportion(melt_wf,"water","yes")
        elif model == "ONeill21hyd": # Eq. (10.34, 10.43, 10.45, 10.46, 10.49)
            # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) and water
            melt_comp = mg.melt_cation_proportion(melt_wf,"water","yes")
        if models.loc["high precision","option"] == "True":
            lnaFeS = gp.log((1.0 - melt_comp["Fe2"])*sulfide_comp)
        else:
            lnaFeS = math.log((1.0 - melt_comp["Fe2"])*sulfide_comp)
        # lnyFe2 from Eq. (10.46)
        lnyFe2 = (((1.0-melt_comp["Fe2"])**2.0)*(28870.0 - 14710.0*melt_comp["Mg"] + 1960.0*melt_comp["Ca"] + 43300.0*melt_comp["Na"] + 95380.0*melt_comp["K"] - 76880.0*melt_comp["Ti"]) + (1.0-melt_comp["Fe2"])*(-62190.0*melt_comp["Si"] + 31520.0*melt_comp["Si"]**2.0))/(R*T)
        # lnS from Eq. (10.43)    
        if models.loc["high precision","option"] == "True":
            lnS = D/(R*T) + gp.log(C_S(PT,melt_wf,models)) - gp.log(melt_comp["Fe2"]) - lnyFe2 + lnaFeS + (-291.0*P + 351.0*gp.erf(P))/T
            SCSS = gp.exp(lnS)  
        else:
            lnS = D/(R*T) + math.log(C_S(PT,melt_wf,models)) - math.log(melt_comp["Fe2"]) - lnyFe2 + lnaFeS + (-291.0*P + 351.0*math.erf(P))/T
            SCSS = math.exp(lnS)  
    
    elif model == "Liu07": # Eq. (9) in Liu et al. (2007) GCA 71:1783-1799 doi:10.1016/j.gca.2007.01.004
        # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) and water
        melt_comp = mg.melt_cation_proportion(melt_wf,"water","yes")
        MFM = (melt_comp["Na"]+melt_comp["K"]+2.*(melt_comp["Ca"]+melt_comp["Mg"]+melt_comp["Fe2"]))/(melt_comp["Si"]*(melt_comp["Al"]+melt_comp["Fe3"]))
        if models.loc["high precision","option"] == "True":
            lnS = 11.35251 - (4454.6/T) - 0.03190*(PT["P"]/T) + 0.71006*gp.log(MFM) - 1.98063*(MFM*melt_comp["H"]) + 0.21867*gp.log(melt_comp["H"]) + 0.36192*gp.log(melt_comp["Fe2"])
            SCSS = gp.exp(lnS)
        else:
            lnS = 11.35251 - (4454.6/T) - 0.03190*(PT["P"]/T) + 0.71006*math.log(MFM) - 1.98063*(MFM*melt_comp["H"]) + 0.21867*math.log(melt_comp["H"]) + 0.36192*math.log(melt_comp["Fe2"])
            SCSS = math.exp(lnS)
        
    elif model == "Fortin15": # Eq. (7) Fortin et al. (2015) GCA 160:100-116 doi:10.1016/j.gca.2015.03.022
        # Mole fractions in the melt on cationic lattice (all Fe as FeOT) and water - molecular masses used are different to spreadsheet attached to paper
        melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no")
        lnS = 34.7837 - (5772.322/T) - 346.5377*((0.0001*PT["P"])/T) - 20.3934*melt_comp["H2O"] - 25.4986*melt_comp["SiO2"] - 18.3435*melt_comp["TiO2"] - 27.3807*melt_comp["Al2O3"] - 17.2752*melt_comp["FeOT"] - 22.3975*melt_comp["MgO"] - 20.3778*melt_comp["CaO"] - 18.9539*melt_comp["Na2O"] - 32.1944*melt_comp["K2O"]
        if models.loc["high precision","option"] == "True":
            SCSS = gp.exp(lnS)
        else:
            SCSS = math.exp(lnS)
    
    elif model == "Liu21": # Eq. (2) Liu et al. (2021) Chem.Geol. 559:119913 doi:10.1016.j.chemgeo.2020.119913
        XFeS = sulf_XFe
        H2O = melt_wf["H2OT"]*100.
        if models.loc["high precision","option"] == "True":
            SCSS = (XFeS*gp.exp(13.88 - (9744./T) - (328.*(0.0001*PT["P"])/T))) + 104.*H2O
        else:
            SCSS = (XFeS*math.exp(13.88 - (9744./T) - (328.*(0.0001*PT["P"])/T))) + 104.*H2O
    
    elif model == "pss_Li22": # Li and Zhang (2022) assuming pure FeS using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_LiZhang2022_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), H2O_Liq=float(100.*melt_wf['H2OT']), Fe_FeNiCu_Sulf=sulf_XFe, Cu_FeNiCu_Sulf=sulf_XCu, Ni_FeNiCu_Sulf=sulf_XNi, Fe3Fet_Liq=Fe3FeT)
        SCSS = float(output["SCSS2_ppm"].iloc[0])
    elif model == "pss_Blanchard21": # Blanchard et al. (2021) assuming pure FeS using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_B2021_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), H2O_Liq=float(100.*melt_wf['H2OT']), Fe_FeNiCu_Sulf=sulf_XFe, Cu_FeNiCu_Sulf=sulf_XCu, Ni_FeNiCu_Sulf=sulf_XNi,  Fe3Fet_Liq=Fe3FeT)
        SCSS = float(output["SCSS2_ppm"].iloc[0])
    elif model == "Fortin15_pss": # Fortin et al. (2015) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_F2015_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), H2O_Liq=float(100.*melt_wf['H2OT']))
        SCSS = float(output["SCSS2_ppm"].iloc[0])
    elif model == "Liu21_pss": # Liu et al. (2021) FeS using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_Liu2021_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), H2O_Liq=float(100.*melt_wf['H2OT']), Fe_FeNiCu_Sulf=sulf_XFe, Cu_FeNiCu_Sulf=sulf_XCu, Ni_FeNiCu_Sulf=sulf_XNi, Fe3Fet_Liq=Fe3FeT)
        SCSS = float(output["SCSS2_ppm"].iloc[0])
    elif model == "ONeill22_pss": # O'Neill & Mavrogenes (2022) FeS using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_OM2022_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), Fe3Fet_Liq=Fe3FeT, Fe_FeNiCu_Sulf=sulf_XFe, Cu_FeNiCu_Sulf=sulf_XCu, Ni_FeNiCu_Sulf=sulf_XNi)
        SCSS = float(output["SCSS2_ppm"].iloc[0])
    elif model == "ONeill21_pss": # O'Neill (2021) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_O2021_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), Fe3Fet_Liq=Fe3FeT, Fe_FeNiCu_Sulf=sulf_XFe, Cu_FeNiCu_Sulf=sulf_XCu, Ni_FeNiCu_Sulf=sulf_XNi)
        SCSS = float(output["SCSS2_ppm"].iloc[0])
    elif model == "Smythe17_pss": # Smythe et al. (2017) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
        output = ss.calculate_S2017_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), Fe3Fet_Liq=Fe3FeT, Fe_FeNiCu_Sulf=sulf_XFe, Cu_FeNiCu_Sulf=sulf_XCu, Ni_FeNiCu_Sulf=sulf_XNi, H2O_Liq=float(100.*melt_wf['H2OT']))
        SCSS = float(output["SCSS2_ppm_ideal_Smythe2017"].iloc[0])
    return SCSS

#################################################################################################################################
#################################### EQUILIBRIUM CONSTANTS FOR HOMOGENEOUS VAPOR EQUILIBRIA #####################################
#################################################################################################################################

# H2 + 0.5O2 = H2O
# K = fH2O/(fH2*(fO2)^0.5)
def KHOg(PT,models=default_models):
    
    """ 
    Equilibrium constant for H2 + 0.5O2 = H2O, K = fH2O/(fH2*(fO2)^0.5)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KHOg" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'Ohmoto97' Reaction (d) in Table 1 of Ohmoto & Kerrick (1997)
    Only one option available currently, included for future development.

    """
    
    model = models.loc["KHOg","option"]

    T_K = PT['T']+273.15
    if model == "Ohmoto97": # Reaction (d) in Table 1 of Ohmoto & Kerrick (1997)
        if models.loc["high precision","option"] == "True":
            K = 10.**((12510.0/T_K)-0.979*(gp.log10(T_K))+0.483)
        else:
            K = 10.**((12510.0/T_K)-0.979*(math.log10(T_K))+0.483)
    return K

# H2O + 0.5S2 = H2S + 0.5O2
# K = (fH2S*(fO2)^0.5)/((fS2^0.5)*fH2O)
def KHOSg(PT,models=default_models):
    """ 
    Equilibrium constant for H2O + 0.5S2 = H2S + 0.5O2, K = (fH2S*(fO2)^0.5)/((fS2^0.5)*fH2O)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KHOSg" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'Ohmoto97' Reaction (h) in Table 1 of Ohmoto & Kerrick (1997)
    Other options:
    'noH2S' Stops H2S forming in the vapor, K = 0.

    """
    
    model = models.loc["KHOSg","option"]

    T_K = PT['T']+273.15
    if model == "Ohmoto97": # Reaction (h) in Table 1 of Ohmoto & Kerrick (1997)
        if models.loc["high precision","option"] == "True":
            K = 10.**((-8117.0/T_K)+0.188*gp.log10(T_K)-0.352)
        else:
            K = 10.**((-8117.0/T_K)+0.188*math.log10(T_K)-0.352)
    elif model == "noH2S": # H2S doesn't form in the gas...
        K = 0.
    return K

# 0.5S2 + O2 = SO2
# K = fSO2/((fS2^0.5)*fO2)
def KOSg(PT,models=default_models):
    
    """ 
    Equilibrium constant for 0.5S2 + O2 = SO2, K = fSO2/((fS2^0.5)*fO2)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KOSg" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'Ohmoto97' Reaction (f) in Table 1 of Ohmoto & Kerrick (1997)
    Other options:
    'noSO2' Stops SO2 forming in the vapor, K = 0.

    """
    
    model = models.loc["KOSg","option"]

    T_K = PT['T']+273.15
    if model == "Ohmoto97": # Reaction (f) in Table 1 of Ohmoto & Kerrick (1997)
        K = 10.**((18929.0/T_K)-3.783)
    if model == 'noSO2':
        K = 0.
    return K

# 0.5S2 + 1.5O2 = SO3
# K = fSO3/((fS2^0.5)*(fO2^1.5)
def KOSg2(PT,models=default_models):
    
    """ 
    Equilibrium constant for 0.5S2 + 1.5O2 = SO3, K = fSO3/((fS2^0.5)*(fO2^1.5)


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KOsg2" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'ONeill22' Eq (6b) in O’Neill and Mavrogenes (2022)
    Only one option available currently, included for future development.

    """

    model = models.loc["KOSg2","option"]

    T_K = PT['T']+273.15
    if model == "ONeill22": # Eq (6b) in O’Neill and Mavrogenes (2022)
        if models.loc["high precision","option"] == "True":
            lnK = (55921./T_K) - 25.07 + 0.6465*gp.log(T_K)
            K = gp.exp(lnK) 
        else:
            lnK = (55921./T_K) - 25.07 + 0.6465*math.log(T_K)
            K = math.exp(lnK)
    return K

# CO + 0.5O = CO2
# K = fCO2/(fCO*(fO2^0.5))
def KCOg(PT,models=default_models):
    
    """ 
    Equilibrium constant for CO + 0.5O = CO2, K = fCO2/(fCO*(fO2^0.5))


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KCOg" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'Ohmoto97' Reaction (c) in Table 1 of Ohmoto & Kerrick (1997)
    Only one option available currently, included for future development.

    """

    model = models.loc["KCOg","option"]

    T_K = PT['T']+273.15
    if model == "Ohmoto97": # Reaction (c) in Table 1 of Ohmoto & Kerrick (1997)
        K = 10.**((14751.0/T_K)-4.535)
    return K

# CH4 + 2O2 = CO2 + 2H2O
# K = (fCO2*(fH2O^2))/(fCH4*(fO2^2))
def KCOHg(PT,models=default_models): 
    
    """ 
    Equilibrium constant for CH4 + 2O2 = CO2 + 2H2O, K = (fCO2*(fH2O^2))/(fCH4*(fO2^2))


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KCOHg" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'Ohmoto97' Reaction (e) in Table 1 of Ohmoto & Kerrick (1997)
    Other options:
    'noCH4' Almost stops CH4 forming in the vapor, K = very large.

    """

    model = models.loc["KCOHg","option"]

    T_K = PT['T']+273.15
    if model == "Ohmoto97": # Reaction (e) in Table 1 of Ohmoto & Kerrick (1997)
        if models.loc["high precision","option"] == "True":
            K = 10.**((41997.0/T_K)+0.719*gp.log10(T_K)-2.404)
        else:
            K = 10.**((41997.0/T_K)+0.719*math.log10(T_K)-2.404)
    if model == "noCH4":
        K = 1.e40 # really big
    return K

def KOCSg(PT,models=default_models): # OCS - depends on system
    
    """ 
    Equilibrium constant for OCS (either K = (fCO2*fH2S)/(fOCS*fH2O) or (fCO^3*fSO2)/(fCO2^2*fOCS))


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KOCSg" and "carbonylsulfide" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'Moussallam19' Eq. (8) in Moussallam et al. (2019) for KOCSg AND 'COS' for carbonlysulfide
    Other options:
    Other options:
    'noOCS' Almost stops OCS forming in the vapor, K = very large.

    """

    reaction = models.loc["carbonylsulfide","option"]
    model = models.loc["KOCSg","option"]

    T = PT['T']+273.15
    if reaction == "COS":
    # 2CO2 + OCS = 3CO + SO2 - 
    # K = (fCO^3*fSO2)/(fCO2^2*fOCS)    
        if model == "Moussallam19": # Eq. (8) in Moussallam et al. (2019)
            K = 10.**(9.24403 - (15386.45/T)) # P and f in bars, T in K 
        if model == "noOCS":
            K = 1.e30 # really big
        return K
    
    ### WORK IN PROGRESS ###
    if reaction == "COHS":
    # OCS + H2O = CO2 + H2S
    # K = (fCO2*fH2S)/(fOCS*fH2O)
        if models == "EVo": 
            if models.loc["high precision","option"] == "True":
                K = gp.exp(0.482 + (16.166e-2/T) + 0.081e-3*T - (5.715e-3/T**2) - 2.224e-1*gp.log(T))
            else:
                K = math.exp(0.482 + (16.166e-2/T) + 0.081e-3*T - (5.715e-3/T**2) - 2.224e-1*math.log(T))
            return K
        
# Cgraphite + O2 = CO2
def KCOs(PT,models=default_models): 
    
    """ 
    Equilibrium constant for Cgraphite + O2 = CO2


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KCOs" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'Holloway92' Eq (3) KI in Holloway et al. (1992) Eur J. Mineral. 4:105-114
    Only one option available currently, included for future development.

    """

    model = models.loc["KCOs","option"]

    T_K = PT['T']+273.15
    P = PT['P']
    if model == "Holloway92": # Eq (3) KI in Holloway et al. (1992) Eur J. Mineral. 4:105-114
        a = 40.07639
        b = -2.5392e-2
        c = 5.27096e-6
        d = 0.0267
        log10K = a + (b*T_K) + (c*T_K**2) + ((P-1)/T_K)
        K = 10.**log10K     
    return K
    
################################################################################################
##################################### SPECIATION CONSTANTS #####################################
################################################################################################
       
# H2Omol + O = 2OH
# K = xOH*2/(xH2Omol*xO)
def KHOm(PT,melt_wf,models=default_models):
    
    """ 
    Speciation constant for H2Omol + O = 2OH- assuming ideal mixing


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "Hspeccomp" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'MORB_HughesIP' Eq. SX in Hughes et al. (in prep) based on data from Dixon et al. (1995)
    Other options:
    StromboliAlkaliBasalt_Lesne10": Eq. (15) Lesne et al. (2010) CMP 162:133-151
    VesuviusAlkaliBasalt_Lesne10: Eq. (16) Lesne et al. (2010) CMP 162:133-151
    EtnaAlkaliBasalt_Lesne10: Eq. (17) Lesne et al. (2010) CMP 162:133-151
    Andesite_Botcharnikov06: Eq (7) from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143
    Albite_Silver89: Fig. 8 from Silver & Stolper (1989) J.Pet 30(3)667-709
    Rhyolite_Zhang97: Eq. (9) from Zhang et al. (1997) GCA 61(15):3089-3100

    """

    Hspeccomp = models.loc["Hspeccomp","option"]
    
    T_K = PT['T']+273.15
    
    if Hspeccomp == "Rhyolite_Zhang97": # Eq. (9) from Zhang et al. (1997) GCA 61(15):3089-3100
        a = -3110.
        b = 1.876
        if models.loc["high precision","option"] == "True":
            K = gp.exp((a/T_K) + b)
        else:
            K = math.exp((a/T_K) + b)
    elif Hspeccomp == "StromboliAlkaliBasalt_Lesne10": # Eq. (15) Lesne et al. (2010) CMP 162:133-151
        a = -8710.
        b = 8.5244
        if models.loc["high precision","option"] == "True":
            K = gp.exp((a/T_K) + b)
        else:
            K = math.exp((a/T_K) + b)
    elif Hspeccomp == "VesuviusAlkaliBasalt_Lesne10": # Eq. (16) Lesne et al. (2010) CMP 162:133-151
        a = -8033.
        b = 7.4222
        if models.loc["high precision","option"] == "True":
            K = gp.exp((a/T_K) + b)
        else:
            K = math.exp((a/T_K) + b)
    elif Hspeccomp == "EtnaAlkaliBasalt_Lesne10": # Eq. (17) Lesne et al. (2010) CMP 162:133-151
        a = -8300.
        b = 7.4859
        if models.loc["high precision","option"] == "True":
            K = gp.exp((a/T_K) + b)
        else:
            K = math.exp((a/T_K) + b)
    elif Hspeccomp == "Andesite_Botcharnikov06": # Eq (7) from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143
        a = -3650.
        b = 2.99
        if models.loc["high precision","option"] == "True":
            K = gp.exp((a/T_K) + b)
        else:
            K = math.exp((a/T_K) + b)
    elif Hspeccomp == "MORB_HughesIP": # fit to Dixon et al. (1995) data digitised from Lesne et al. (2010) CMP 162:133-151 in Hughes et al. (in prep)
        a = -2204.99
        b = 1.2600
        if models.loc["high precision","option"] == "True":
            K = gp.exp((a/T_K) + b)
        else:
            K = math.exp((a/T_K) + b)
    elif Hspeccomp == "Albite_Silver89": # Fig. 8 from Silver & Stolper (1989) J.Pet 30(3)667-709
        K = 0.17
    
    ### Work in progress ###
    elif Hspeccomp == "AlkaliBasalt": # average of eqn-15-17 from Lesne et al. (2010) CMP 162:133-151
        a = -8348.0 # VES-9 = -8033.0, ETN-1 = -8300.0, and PST-9 = -8710.0
        b = 7.8108 # VES-9 = 7.4222, ETN-1 = 7.4859, and PEST-9 = 8.5244
        if models.loc["high precision","option"] == "True":
            K = gp.exp((a/T_K) + b)
        else:
            K = math.exp((a/T_K) + b)    
    
    return K

def KregH2O(PT,melt_wf,models=default_models):
    
    """ 
    Speciation constant for H2Omol + O = 2OH- assuming regular mixing


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "Hspeccomp" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'MORB_HughesIP' WHICH IS NOT USABLE WITH THIS FUNCTION
    Other options:
    MORB_Dixon95: Table 5 of Dixon et al. (1995)
    AlkaliBasalt_Lesne10": Eq. (24-27) Lesne et al. (2010) CMP 162:133-151
    StromboliAlkaliBasalt_Lesne10": PST-9 in Table 5 from Lesne et al. (2010) 162:133-151
    VesuviusAlkaliBasalt_Lesne10: VES-9 in Table 5 from Lesne et al. (2010) 162:133-151
    EtnaAlkaliBasalt_Lesne10: ETN-1 in Table 5 from Lesne et al. (2010) 162:133-151
    Albite_Silver89: Silver & Stolper (1989) J.Pet 30(3)667-709

    """

    Hspeccomp = models.loc["Hspeccomp","option"]

    if Hspeccomp in ["MORB_Dixon95","Albite_Silver89"]: # Table 5 of Dixon et al. (1995); Silver & Stolper (1989) J.Pet 30(3)667-709
        A = 0.403
        B = 15.333
        C = 10.894
    elif Hspeccomp == "AlkaliBasalt_Lesne10": # Eq (24-27) from Lesne et al. (2010) CMP 162:133-151
        lnK = ((-2704.4/T_K) + 0.641)
        A = lnK + 49016.0/(R*T_K)
        B = -2153326.51/(R*T_K)
        C = 1.965495217/(R*T_K)
    elif Hspeccomp == "VesuviusAlkaliBasalt_Lesne10": # VES-9 in Table 5 from Lesne et al. (2010) 162:133-151
        A = 3.139
        B = -29.555
        C = 20.535
    elif Hspeccomp == "EtnaAlkaliBasalt_Lesne10": # ETN-1 in Table 5 from Lesne et al. (2010) 162:133-151
        A = 4.128
        B = -45.905
        C = 21.311
    elif Hspeccomp == "StromboliAlkaliBasalt_Lesne10": # PST-9 in Table 5 from Lesne et al. (2010) 162:133-151
        A = 2.6
        B = -22.476
        C = 22.295
    
    ### Work in progress ###
    elif Hspeccomp == "alkali basalt XT": # No T-dependence, hence its the speciation frozen in the glass. Eqn 7-10 from Lesne et al. (2010) CMP 162:133-151 (eqn 7 is wrong)
        # wt% normalised including H2O, all Fe as FeOT
        melt_comp = melt_normalise_wf(melt_wf,"volatiles","Fe speciation")
        Na = melt_comp["Na2O"]*100.
        K = melt_comp["K2O"]*100.
        A = 0.5761*(Na+K) - 0.2884 # eqn-8
        B = -8.9589*(Na+K) + 24.65 # eqn-9
        C = 1.7013*(Na+K) + 9.6481 # eqn-1
       
    
# CO2 + O = CO3
def KCOm(PT,melt_wf,models=default_models): # K = 
    """ 
    Speciation constant for CO2 + O = CO3


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "Cspeccomp" and column label of "option"

    Returns
    -------
    Equilibrium constant as <class 'mpfr'>

    Model options
    -------------
    default: 'Basalt' Assume all oxidised carbon in the melt is present as carbonate ions.
    Other options:
    Andesite_Botcharnikov06: Eq. (8) from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143
    Dacite_Botcharnikov06: Eq. in the text from Botcharnikov et al. (2006), based on data from Behrens et al. (2004)
    Rhyolite: Assume all oxidised carbon in the melt is present as molecular CO2.

    """
    
    Cspeccomp = models.loc["Cspeccomp","option"]

    T_K = PT['T']+273.15
    
    if Cspeccomp == "Andesite_Botcharnikov06": # Eq. (8) from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143
        a = 8665.0
        b = -5.11
        if models.loc["high precision","option"] == "True":
            value = gp.exp((a/T_K) + b)  
        else:
            value = math.exp((a/T_K) + b)  
    elif Cspeccomp == "Dacite_Botcharnikov06": # Eq. in the text from Botcharnikov et al. (2006), based on data from Behrens et al. (2004)
        a = 9787.0
        b = -7.69
        if models.loc["high precision","option"] == "True":
            value = gp.exp((a/T_K) + b)  
        else:
            value = math.exp((a/T_K) + b)
    elif Cspeccomp == "Basalt": # all oxidised carbon is CO32-
        value = "infinite"
    elif Cspeccomp == "Rhyolite": # all oxidised carbon is CO2,mol
        value = 0.
    return value


################################################################################################################################# 
##################################################### FUGACITY COEFFICIENTS #####################################################
#################################################################################################################################

# all fugacity coefficients are assumed to equal 1 below 1 bar.

##########################################
### CORK using Holland & Powell (1991) ###
##########################################

def CORK(PT,p0,a,b,c,d,models):
    P = PT['P']
    T_K = PT['T']+273.15
    def MRK(P_kb,VMRK,R,T_K,a,b): # MRK volume equation rearranged to equal 0
        return P_kb*pow(VMRK,3.0) - R*T_K*pow(VMRK,2.0) - (b*R*T_K + pow(b,2.0)*P_kb - a*pow(T_K,-0.5))*VMRK - (a*b)*pow(T_K,-0.5)

    def dMRK(P_kb,VMRK,R,T_K,a,b): # derivative of above
        return 3.0*P_kb*pow(VMRK,2.0) - 2.0*R*T_K*VMRK - (b*R*T_K + pow(b,2.0)*P_kb - a*pow(T_K,-0.5))

    def dVMRK(MRK,P_kb,VMRK,R,T_K,a,b):
        return abs(0-MRK(P_kb,VMRK,R,T_K,a,b))

    def NR_VMRK(MRK, dMRK, VMRK0, e1, P_kb,R,T_K,a,b):
        delta1 = dVMRK(MRK,P_kb,VMRK0,R,T_K,a,b)
        while delta1 > e1:
            VMRK0 = VMRK0 - MRK(P_kb,VMRK0,R,T_K,a,b)/dMRK(P_kb,VMRK0,R,T_K,a,b)
            delta1 = dVMRK(MRK,P_kb,VMRK0,R,T_K,a,b)
        return VMRK0
    
    R = 8.314e-3 # in kJ/mol/K
    P_kb = P/1000.0
    
    Vi = ((R*T_K)/P_kb) + b
        
    VMRK = NR_VMRK(MRK, dMRK, Vi, 1E20, P_kb,R,T_K,a,b)
        
    if P_kb > p0:
        V = VMRK + c*pow((P_kb-p0),0.5) + d*(P_kb-p0)
        ln_y_virial = (1/(R*T_K))*((2./3.)*c*pow((P_kb-p0),1.5) + (d/2.0)*pow((P_kb-p0),2.0))
    else:
        V = VMRK
        ln_y_virial = 0.0
        
    z = (P_kb*V)/(R*T_K)
    A = a/(b*R*pow(T_K,1.5))
    B = (b*P_kb)/(R*T_K)
        
    if models.loc["high precision","option"] == "True":
        ln_y = z - 1.0 - gp.log(z-B) - A*gp.log(1.0 + (B/z)) + ln_y_virial
        value = gp.exp(ln_y)
    else:
        ln_y = z - 1.0 - math.log(z-B) - A*math.log(1.0 + (B/z)) + ln_y_virial
        value = math.exp(ln_y)
    return value
    
###########################
### Shi & Saxena (1992) ###
###########################

def lny_SS(PT,Pcr,Tcr,models):
    P = PT['P']
    T_K = PT['T']+273.15
    Tr = T_K/Tcr
    A, B, C, D, P0, integral0 = Q_SS(PT,Tr,Pcr,models)
    Pr = P/Pcr
    P0r = P0/Pcr
    if models.loc["high precision","option"] == "True":
        integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
    else:
        integral = A*math.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
    integral_total = integral + integral0
    return integral_total

def Q_SS(PT,Tr,Pcr,models):
    P = PT['P']
    def Q1000(Pcr):
        Pr_ = 1000.0/Pcr
        P0r_ = 1.0/Pcr
        A0 = 1.0
        B0 = 0.9827e-1*pow(Tr,-1.0) + -0.2709*pow(Tr,-3.0)
        C0 = -0.1030e-2*pow(Tr,-1.5) + 0.1427e-1*pow(Tr,-4.0)
        D0 = 0.0
        if models.loc["high precision","option"] == "True":
            value = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
        else:
            value = A0*math.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
        return value
    
    def Q5000(Pcr):
        Pr_ = 5000.0/Pcr
        P0r_ = 1000.0/Pcr
        A0 = 1.0 + -5.917e-1*pow(Tr,-2.0)
        B0 = 9.122e-2*pow(Tr,-1.0)
        D0 = 0.0
        if models.loc["high precision","option"] == "True":
            C0 = -1.416e-4*pow(Tr,-2.0) + -2.835e-6*gp.log(Tr)
            value = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
        else:
            C0 = -1.416e-4*pow(Tr,-2.0) + -2.835e-6*math.log(Tr)
            value = A0*math.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
        return value
    if P > 5000.0:
        if models.loc["high precision","option"] == "True":
            A = 2.0614 + -2.235*pow(Tr,-2.0) + -3.941e-1*gp.log(Tr)
        else:
            A = 2.0614 + -2.235*pow(Tr,-2.0) + -3.941e-1*math.log(Tr)
        B = 5.513e-2*pow(Tr,-1.0) + 3.934e-2*pow(Tr,-2.0)
        C = -1.894e-6*pow(Tr,-1.0) + -1.109e-5*pow(Tr,-2.0) + -2.189e-5*pow(Tr,-3.0)
        D = 5.053e-11*pow(Tr,-1.0) + -6.303e-21*pow(Tr,3.0)
        P0 = 5000.0
        integral0 = Q1000(Pcr) + Q5000(Pcr)
        return A, B, C, D, P0, integral0
    elif P == 5000.0:
        A = 0
        B = 0
        C = 0
        D = 0
        P0 = 5000.0
        integral0 = Q1000(Pcr) + Q5000(Pcr)
        return A, B, C, D, P0, integral0
    elif P > 1000.0 and P < 5000.0:
        A = 1.0 + -5.917e-1*pow(Tr,-2.0)
        B = 9.122e-2*pow(Tr,-1.0)
        if models.loc["high precision","option"] == "True":
            C = -1.416e-4*pow(Tr,-2.0) + -2.835e-6*gp.log(Tr)
        else:
            C = -1.416e-4*pow(Tr,-2.0) + -2.835e-6*math.log(Tr)
        D = 0.0
        P0 = 1000.0
        integral0 = Q1000(Pcr)
        return A, B, C, D, P0, integral0
    elif P == 1000.0:
        A = 0
        B = 0
        C = 0
        D = 0.0
        P0 = 1000.0
        integral0 = Q1000(Pcr)
        return A, B, C, D, P0, integral0
    else:
        A = 1.0
        B = 0.9827e-1*pow(Tr,-1.0) + -0.2709*pow(Tr,-3.0)
        C = -0.1030e-2*pow(Tr,-1.5) + 0.1427e-1*pow(Tr,-4.0)
        D = 0.0
        P0 = 1.0
        integral0 = 0.0
        return A, B, C, D, P0, integral0
    
def y_SS(gas_species,PT,models=default_models):
    P = PT['P']
    T_K = PT['T']+273.15
    
    ideal_gas = models.loc["ideal_gas","option"]

    if ideal_gas == "True":
        return 1.0
    elif P < 1.: # ideal gas below 1 bar
        return 1.
    else:    
        Tcr = species.loc[gas_species,"Tcr"]
        Pcr = species.loc[gas_species,"Pcr"]
        if models.loc["high precision","option"] == "True":
            value = gp.exp(lny_SS(PT,Pcr,Tcr,models))/P    
        else:
            value = math.exp(lny_SS(PT,Pcr,Tcr,models))/P    
        return value

##############################
### for each vapor species ###    
##############################

def y_H2(PT,models=default_models):
    """ 
    Fugacity coefficient for H2

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_H2" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_H2
    ----------------------
    default: 'Shaw64' Eq. (4) from Shaw & Wones (1964) AmJSci 262:918-929
    Other options:
    ideal: Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    P = PT['P']
    T_K = PT['T']+273.15

    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_H2","option"]

    if ideal_gas == "True" or model == "ideal":
        return 1.0
    elif P < 1.: # ideal gas below 1 bar
        return 1.
    elif model == "Shaw64": # Eq. (4) from Shaw & Wones (1964) AmJSci 262:918-929
        P_atm = 0.986923*P
        if models.loc["high precision","option"] == "True":
            SW1 = gp.exp(-3.8402*pow(T_K,0.125)+0.5410)
            SW2 = gp.exp(-0.1263*pow(T_K,0.5)-15.980)
            SW3 = 300*gp.exp((-0.011901*T_K)-5.941) # NB used a value of -0.011901 instead of -0.11901 as reported to match data in Table 2
            ln_y = SW1*P_atm - SW2*pow(P_atm,2.0) + SW3*gp.exp((-P_atm/300.0)-1.0)
            value = gp.exp(ln_y)
        else:
            SW1 = math.exp(-3.8402*pow(T_K,0.125)+0.5410)
            SW2 = math.exp(-0.1263*pow(T_K,0.5)-15.980)
            SW3 = 300*math.exp((-0.011901*T_K)-5.941) # NB used a value of -0.011901 instead of -0.11901 as reported to match data in Table 2
            ln_y = SW1*P_atm - SW2*pow(P_atm,2.0) + SW3*math.exp((-P_atm/300.0)-1.0)
            value = math.exp(ln_y)
        return value
    ### work in progress ###
    elif model == "Shi92": # Shi & Saxena (1992) NOT WORKING
        Tcr = 33.25 # critical temperature in K 
        Pcr = 12.9696 # critical temperature in bar
        Tr = T_K/Tcr
        # Q for 1-1000 bar
        Q1_A_LP, Q2_A_LP, Q3_A_LP, Q4_A_LP, Q5_A_LP = 1., 0., 0., 0., 0.
        Q1_B_LP, Q2_B_LP, Q3_B_LP, Q4_B_LP, Q5_B_LP = 0., 0.9827e-1, 0., -0.2709, 0.
        Q1_C_LP, Q2_C_LP, Q3_C_LP, Q4_C_LP, Q5_C_LP = 0., 0., -0.1030e-2, 0., 0.1427e-1
        # Q for 1000-10000 bar
        Q1_A_HP, Q2_A_HP, Q3_A_HP, Q4_A_HP, Q5_A_HP, Q6_A_HP, Q7_A_HP, Q8_A_HP = 2.2615, 0., -6.8712e1, 0., -1.0573e4, 0., 0., -1.6936e-1
        Q1_B_HP, Q2_B_HP, Q3_B_HP, Q4_B_HP, Q5_B_HP, Q6_B_HP, Q7_B_HP, Q8_B_HP = -2.6707e-4, 0., 2.0173e-1, 0., 4.5759, 0., 0., 3.1452e-5
        Q1_C_HP, Q2_C_HP, Q3_C_HP, Q4_C_HP, Q5_C_HP, Q6_C_HP, Q7_C_HP, Q8_C_HP = -2.3376e-9, 0., 3.4091e-7, 0., -1.4188e-3, 0., 0., 3.0117e-10
        Q1_D_HP, Q2_D_HP, Q3_D_HP, Q4_D_HP, Q5_D_HP, Q6_D_HP, Q7_D_HP, Q8_D_HP = -3.2606e-15, 0., 2.4402e-12, 0., -2.4027e-9, 0., 0., 0.
        if P < 1000.:
            A = Q1_A_LP + Q2_A_LP*Tr**(-1.) + Q3_A_LP*Tr**(-1.5) + Q4_A_LP*Tr**(-3.) + Q5_A_LP*Tr**(-4.)
            B = Q1_B_LP + Q2_B_LP*Tr**(-1.) + Q3_B_LP*Tr**(-1.5) + Q4_B_LP*Tr**(-3.) + Q5_B_LP*Tr**(-4.)
            C = Q1_C_LP + Q2_C_LP*Tr**(-1.) + Q3_C_LP*Tr**(-1.5) + Q4_C_LP*Tr**(-3.) + Q5_C_LP*Tr**(-4.)
            D = 0.0
            P0 = 1.0
            integral0 = 0.
        elif P == 1000.:
            A = 0.0
            B = 0.0
            C = 0.0
            D = 0.0
            P0 = 1000.0
            Pr_ = 1000.0/Pcr
            P0r_ = 1.0/Pcr
            A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.)
            B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.)
            C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.)
            D0 = 0.0
            if models.loc["high precision","option"] == "True":
                integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))   
            else:
                integral0 = A0*math.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))   
        elif P > 1000.:
            if models.loc["high precision","option"] == "True":
                A = Q1_A_HP + Q2_A_HP*Tr + Q3_A_HP*Tr**(-1.) + Q4_A_HP*Tr**2. + Q5_A_HP*Tr**(-2.) + Q6_A_HP*Tr**3. + Q7_A_HP*Tr**(-3.0) + Q8_A_HP*gp.log(Tr)
                B = Q1_B_HP + Q2_B_HP*Tr + Q3_B_HP*Tr**(-1.) + Q4_B_HP*Tr**2. + Q5_B_HP*Tr**(-2.) + Q6_B_HP*Tr**3. + Q7_B_HP*Tr**(-3.0) + Q8_B_HP*gp.log(Tr)
                C = Q1_C_HP + Q2_C_HP*Tr + Q3_C_HP*Tr**(-1.) + Q4_C_HP*Tr**2. + Q5_C_HP*Tr**(-2.) + Q6_C_HP*Tr**3. + Q7_C_HP*Tr**(-3.0) + Q8_C_HP*gp.log(Tr)
                D = Q1_D_HP + Q2_D_HP*Tr + Q3_D_HP*Tr**(-1.) + Q4_D_HP*Tr**2. + Q5_D_HP*Tr**(-2.) + Q6_D_HP*Tr**3. + Q7_D_HP*Tr**(-3.0) + Q8_D_HP*gp.log(Tr)
            else:
                A = Q1_A_HP + Q2_A_HP*Tr + Q3_A_HP*Tr**(-1.) + Q4_A_HP*Tr**2. + Q5_A_HP*Tr**(-2.) + Q6_A_HP*Tr**3. + Q7_A_HP*Tr**(-3.0) + Q8_A_HP*math.log(Tr)
                B = Q1_B_HP + Q2_B_HP*Tr + Q3_B_HP*Tr**(-1.) + Q4_B_HP*Tr**2. + Q5_B_HP*Tr**(-2.) + Q6_B_HP*Tr**3. + Q7_B_HP*Tr**(-3.0) + Q8_B_HP*math.log(Tr)
                C = Q1_C_HP + Q2_C_HP*Tr + Q3_C_HP*Tr**(-1.) + Q4_C_HP*Tr**2. + Q5_C_HP*Tr**(-2.) + Q6_C_HP*Tr**3. + Q7_C_HP*Tr**(-3.0) + Q8_C_HP*math.log(Tr)
                D = Q1_D_HP + Q2_D_HP*Tr + Q3_D_HP*Tr**(-1.) + Q4_D_HP*Tr**2. + Q5_D_HP*Tr**(-2.) + Q6_D_HP*Tr**3. + Q7_D_HP*Tr**(-3.0) + Q8_D_HP*math.log(Tr)
            P0 = 1000.0
            Pr_ = 1000.0/Pcr
            P0r_ = 1.0/Pcr
            A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.)
            B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.)
            C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.)
            D0 = 0.0
            if models.loc["high precision","option"] == "True":
                integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
            else:
                integral0 = A0*math.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
        P0r = P0/Pcr
        Pr = P/Pcr
        if models.loc["high precision","option"] == "True":
            integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
            value = gp.exp(integral + integral0)/P
        else:
            integral = A*math.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
            value = math.exp(integral + integral0)/P
        return value

def y_H2O(PT,models=default_models):
    """ 
    Fugacity coefficient for H2O

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_H2O" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_H2O
    -----------------------
    default: 'Holland91' Holland & Powell (1991) CMP 109:265-273 10.1007/BF00306484
    Other options:
    ideal: Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_H2O","option"]

    P = PT['P']
    T_K = PT['T']+273.15

    if ideal_gas == "True" or model == "ideal":
        return 1.
    elif P < 1.: # ideal gas below 1 bar
        return 1.
    if model == "Holland91": # Holland & Powell (1991) CMP 109:265-273 10.1007/BF00306484
    # (T > 673 K only) - using Holland & Powell (1991) CORK
        p0 = 2.00 # in kb
        a = 1113.4 + -0.22291*(T_K - 673.0) + -3.8022e-4*pow((T_K-673.0),2.0) + 1.7791e-7*pow((T_K-673.0),3.0)
        b = 1.465
        c = -3.025650e-2 + -5.343144e-6*T_K
        d = -3.2297554e-3 + 2.2215221e-6*T_K
        y = CORK(PT,p0,a,b,c,d,models)
        return y

def y_CO2(PT,models=default_models):
    """ 
    Fugacity coefficient for CO2

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_CO2" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_CO2
    -----------------------
    default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    Other options:
    ideal: Treat as ideal gas, y = 1 at all P.
    Holland91: Holland & Powell (1991) CMP 109:265-273 10.1007/BF00306484
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_CO2","option"]

    P = PT['P']
    T_K = PT['T']+273.15

    if ideal_gas == "True" or model == "ideal":
        return 1.0
    elif P < 1.: # ideal gas below 1 bar
        return 1.
    else:
        if model == "Holland91": # Holland & Powell (1991) CMP 109:265-273 10.1007/BF00306484
            p0 = 5.00 # in kb
            a = 741.2 + -0.10891*(T_K) + -3.4203e-4*pow(T_K,2.0)
            b = 3.057
            c = -2.26924e-1 + -7.73793e-5*T_K
            d = 1.33790e-2 + -1.1740e-5*T_K
            y = CORK(PT,p0,a,b,c,d,models)
        elif model == "Shi92": # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
            gas_species = "CO2"
            y = y_SS(gas_species,PT,models)
        return y
    
def y_O2(PT,models=default_models):
    """ 
    Fugacity coefficient for O2

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_O2" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_O2
    ----------------------
    default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    Other options:
    ideal: Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_O2","option"]

    if model == "Shi92": # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "O2"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y
    
def y_S2(PT,models=default_models):
    """ 
    Fugacity coefficient for S2

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_S2" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_S2
    ----------------------
    default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    Other options:
    ideal: Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_S2","option"]

    if model == "Shi92": # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "S2"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y

def y_CO(PT,models=default_models):
    """ 
    Fugacity coefficient for CO

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_CO" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_CO
    ----------------------
    default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    Other options:
    ideal: Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_CO","option"]

    if model == "Shi92": # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "CO"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y
    
def y_CH4(PT,models=default_models):
    """ 
    Fugacity coefficient for CH4

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_CH4" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_CH4
    -----------------------
    default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    Other options:
    ideal: Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_CH4","option"]

    if model == "Shi92": # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "CH4"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y
    
def y_OCS(PT,models=default_models):
    """ 
    Fugacity coefficient for OCS

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_OCS" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_OCS
    -----------------------
    default: 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    Other options:
    ideal: Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_OCS","option"]

    if model == "Shi92": # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "OCS"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y

def y_X(PT,models=default_models): # species X fugacity coefficient
    """ 
    Fugacity coefficient for X

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_X" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_X
    -----------------------
    default: "ideal" Treat as ideal gas, y = 1 at all P.
    Only one option available currently, included for future development.        
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_X","option"]

    if model == "ideal":  # ideal gas
        y = 1.
    return y

#################################################################################
### SO2 and H2S from Shi & Saxena (1992) with option to modify below 500 bars ###
#################################################################################

def y_SO2(PT,models=default_models):
    """ 
    Fugacity coefficient for SO2

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_SO2" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_SO2
    -----------------------
    default: 'Shi92_Hughes23' Fig.S1 modified from Shi & Saxena (1992) from Hughes et al. (2023) JGSL 180(3) doi:10.1144/jgs2021-12
    Other options:
    Shi92: Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    ideal: Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_SO2","option"]

    P = PT['P']
    T_K = PT['T']+273.15

    gas_species = "SO2"
    if ideal_gas == "True" or model == "ideal":
        return 1.
    elif P < 1.: # ideal gas below 1 bar
        return 1.
    else: # 1-10000 bar
        Tcr = species.loc[gas_species,"Tcr"] # critical temperature in K
        Pcr = species.loc[gas_species,"Pcr"] # critical temperature in bar
        P0 = 1.0
        P0r = P0/Pcr
        Tr = T_K/Tcr
        Q1_A, Q2_A, Q3_A, Q4_A, Q5_A, Q6_A, Q7_A, Q8_A  = 0.92854, 0.43269e-1, -0.24671, 0., 0.24999, 0., -0.53182, -0.16461e-1
        Q1_B, Q2_B, Q3_B, Q4_B, Q5_B, Q6_B, Q7_B, Q8_B  = 0.84866e-3, -0.18379e-2, 0.66787e-1, 0., -0.29427e-1, 0., 0.29003e-1, 0.54808e-2
        Q1_C, Q2_C, Q3_C, Q4_C, Q5_C, Q6_C, Q7_C, Q8_C  = -0.35456e-3, 0.23316e-4, 0.94159e-3, 0., -0.81653e-3, 0., 0.23154e-3, 0.55542e-4
        if models.loc["high precision","option"] == "True":
            A = Q1_A + Q2_A*Tr + Q3_A*Tr**(-1.) + Q4_A*Tr**2. + Q5_A*Tr**(-2.) + Q6_A*Tr**3. + Q7_A*Tr**(-3.0) + Q8_A*gp.log(Tr)
            B = Q1_B + Q2_B*Tr + Q3_B*Tr**(-1.) + Q4_B*Tr**2. + Q5_B*Tr**(-2.) + Q6_B*Tr**3. + Q7_B*Tr**(-3.0) + Q8_B*gp.log(Tr)
            C = Q1_C + Q2_C*Tr + Q3_C*Tr**(-1.) + Q4_C*Tr**2. + Q5_C*Tr**(-2.) + Q6_C*Tr**3. + Q7_C*Tr**(-3.0) + Q8_C*gp.log(Tr)
        else:
            A = Q1_A + Q2_A*Tr + Q3_A*Tr**(-1.) + Q4_A*Tr**2. + Q5_A*Tr**(-2.) + Q6_A*Tr**3. + Q7_A*Tr**(-3.0) + Q8_A*math.log(Tr)
            B = Q1_B + Q2_B*Tr + Q3_B*Tr**(-1.) + Q4_B*Tr**2. + Q5_B*Tr**(-2.) + Q6_B*Tr**3. + Q7_B*Tr**(-3.0) + Q8_B*math.log(Tr)
            C = Q1_C + Q2_C*Tr + Q3_C*Tr**(-1.) + Q4_C*Tr**2. + Q5_C*Tr**(-2.) + Q6_C*Tr**3. + Q7_C*Tr**(-3.0) + Q8_C*math.log(Tr)
        D = 0.0
        if P >= 500.: # above 500 bar Shi & Saxena (1992) AmMin 77(9-10):1038-1049
            Pr = P/Pcr
            if models.loc["high precision","option"] == "True":
                integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
                y = (gp.exp(integral))/P
            else:
                integral = A*math.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
                y = (math.exp(integral))/P
        elif models.loc["y_SO2","option"] == "Shi92": # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
            Pr = P/Pcr
            if models.loc["high precision","option"] == "True":
                integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
                y = (gp.exp(integral))/P
            else:
                integral = A*math.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
                y = (math.exp(integral))/P
        elif models.loc["y_SO2","option"] == "Shi92_Hughes23": # Fig.S1 Hughes et al. (2023) JGSL 180(3) doi:10.1144/jgs2021-12
            Pr = 500./Pcr # calculate y at 500 bar
            if models.loc["high precision","option"] == "True":
                integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
                y_500 = (gp.exp(integral))/500.
            else:
                integral = A*math.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
                y_500 = (math.exp(integral))/500.
            y = ((y_500 - 1.)*(P/500.)) + 1. # linear extrapolation to P of interest
        return y       
            
def y_H2S(PT,models=default_models):
    """ 
    Fugacity coefficient for H2S

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_H2S" and "ideal_gas" and column label of "option"

    Returns
    -------
    Fugacity coefficient as <class 'mpfr'>

    Model options for y_H2S
    -----------------------
    default: 'Shi92_Hughes23' Fig.S1 modified from Shi & Saxena (1992) Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
    Other options:
    Shi92: Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    ideal: Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_H2S","option"]

    P = PT['P']
    T_K = PT['T']+273.15

    gas_species = "H2S"
    if ideal_gas == "True" or model == "ideal":
        return 1.0
    elif ideal_gas == "False":
        Tcr = species.loc[gas_species,"Tcr"] # critical temperature in K 
        Pcr = species.loc[gas_species,"Pcr"] # critical temperature in bar
        Tr = T_K/Tcr
        # Q for 1-500 bar
        Q1_A_LP, Q2_A_LP, Q3_A_LP, Q4_A_LP, Q5_A_LP, Q6_A_LP, Q7_A_LP, Q8_A_LP = 0.14721e1, 0.11177e1, 0.39657e1, 0., -0.10028e2, 0., 0.45484e1, -0.382e1
        Q1_B_LP, Q2_B_LP, Q3_B_LP, Q4_B_LP, Q5_B_LP, Q6_B_LP, Q7_B_LP, Q8_B_LP = 0.16066, 0.10887, 0.29014, 0., -0.99593, 0., -0.18627, -0.45515
        Q1_C_LP, Q2_C_LP, Q3_C_LP, Q4_C_LP, Q5_C_LP, Q6_C_LP, Q7_C_LP, Q8_C_LP = -0.28933, -0.70522e-1, 0.39828, 0., -0.50533e-1, 0., 0.1176, 0.33972
        # Q for 500-10000 bar
        Q1_A_HP, Q2_A_HP, Q3_A_HP, Q4_A_HP, Q5_A_HP, Q6_A_HP, Q7_A_HP, Q8_A_HP = 0.59941, -0.1557e-2, 0.4525e-1, 0., 0.36687, 0., -0.79248, 0.26058
        Q1_B_HP, Q2_B_HP, Q3_B_HP, Q4_B_HP, Q5_B_HP, Q6_B_HP, Q7_B_HP, Q8_B_HP = 0.22545e-1, 0.17473e-2, 0.48253e-1, 0., -0.1989e-1, 0., 0.32794e-1, -0.10985e-1
        Q1_C_HP, Q2_C_HP, Q3_C_HP, Q4_C_HP, Q5_C_HP, Q6_C_HP, Q7_C_HP, Q8_C_HP = 0.57375e-3, -0.20944e-5, -0.11894e-2, 0., 0.14661e-2, 0., -0.75605e-3, -0.27985e-3
        if P < 1.:
            return 1. # ideal gas below 1 bar
        elif P < 500.:
            if models.loc["y_H2S","option"] == "Shi92": # as is Shi and Saxena (1992) 
                if models.loc["high precision","option"] == "True":
                    A = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*gp.log(Tr)
                    B = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*gp.log(Tr)
                    C = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*gp.log(Tr)
                else:
                    A = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*math.log(Tr)
                    B = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*math.log(Tr)
                    C = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*math.log(Tr)
                D = 0.0
                P0 = 1.0
                integral0 = 0.
            elif models.loc["y_H2S","option"] == "Shi92_Hughes24": # Fig.S1 modified from Shi & Saxena (1992) Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739 
                P0 = 500.0 # calculate y at 500 bars
                Pr_ = 500.0/Pcr
                P0r_ = 1.0/Pcr
                D0 = 0.0
                if models.loc["high precision","option"] == "True":
                    A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*gp.log(Tr)
                    B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*gp.log(Tr)
                    C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*gp.log(Tr)
                    integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))            
                    y_500 = gp.exp(integral0)/500.
                else:
                    A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*math.log(Tr)
                    B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*math.log(Tr)
                    C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*math.log(Tr)
                    integral0 = A0*math.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))            
                    y_500 = math.exp(integral0)/500.
                y = ((y_500 - 1.)*(P/500.)) + 1. # linear extrapolation to P of interest
                return y
        elif P == 500.:
            A = 0.0
            B = 0.0
            C = 0.0
            D = 0.0
            P0 = 500.0
            Pr_ = 500.0/Pcr
            P0r_ = 1.0/Pcr
            D0 = 0.0
            if models.loc["high precision","option"] == "True":
                A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*gp.log(Tr)
                B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*gp.log(Tr)
                C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*gp.log(Tr)
                integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))        
            else:
                A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*math.log(Tr)
                B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*math.log(Tr)
                C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*math.log(Tr)
                integral0 = A0*math.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))     
        elif P > 500.:
            
            D = 0.0
            P0 = 500.0
            Pr_ = 500.0/Pcr
            P0r_ = 1.0/Pcr
            D0 = 0.0
            if models.loc["high precision","option"] == "True":
                A = Q1_A_HP + Q2_A_HP*Tr + Q3_A_HP*Tr**(-1.) + Q4_A_HP*Tr**2. + Q5_A_HP*Tr**(-2.) + Q6_A_HP*Tr**3. + Q7_A_HP*Tr**(-3.0) + Q8_A_HP*gp.log(Tr)
                B = Q1_B_HP + Q2_B_HP*Tr + Q3_B_HP*Tr**(-1.) + Q4_B_HP*Tr**2. + Q5_B_HP*Tr**(-2.) + Q6_B_HP*Tr**3. + Q7_B_HP*Tr**(-3.0) + Q8_B_HP*gp.log(Tr)
                C = Q1_C_HP + Q2_C_HP*Tr + Q3_C_HP*Tr**(-1.) + Q4_C_HP*Tr**2. + Q5_C_HP*Tr**(-2.) + Q6_C_HP*Tr**3. + Q7_C_HP*Tr**(-3.0) + Q8_C_HP*gp.log(Tr)
                A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*gp.log(Tr)
                B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*gp.log(Tr)
                C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*gp.log(Tr)
                integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
            else:
                A = Q1_A_HP + Q2_A_HP*Tr + Q3_A_HP*Tr**(-1.) + Q4_A_HP*Tr**2. + Q5_A_HP*Tr**(-2.) + Q6_A_HP*Tr**3. + Q7_A_HP*Tr**(-3.0) + Q8_A_HP*math.log(Tr)
                B = Q1_B_HP + Q2_B_HP*Tr + Q3_B_HP*Tr**(-1.) + Q4_B_HP*Tr**2. + Q5_B_HP*Tr**(-2.) + Q6_B_HP*Tr**3. + Q7_B_HP*Tr**(-3.0) + Q8_B_HP*math.log(Tr)
                C = Q1_C_HP + Q2_C_HP*Tr + Q3_C_HP*Tr**(-1.) + Q4_C_HP*Tr**2. + Q5_C_HP*Tr**(-2.) + Q6_C_HP*Tr**3. + Q7_C_HP*Tr**(-3.0) + Q8_C_HP*math.log(Tr)
                A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*math.log(Tr)
                B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*math.log(Tr)
                C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*math.log(Tr)
                integral0 = A0*math.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
        P0r = P0/Pcr
        Pr = P/Pcr
        if models.loc["high precision","option"] == "True":
            integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
            value = gp.exp(integral + integral0)/P
        else:
            integral = A*math.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
            value = math.exp(integral + integral0)/P
        return value
    
def y_SO3(PT,models=default_models):
    return 1.

#######################
### oxygen fugacity ###
#######################

# buffers
def NNO(PT,models=default_models):
    """ 
    Value of NNO buffer

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "NNObuffer" and column label of "option"

    Model options
    -------------
    default: 'Frost91' Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance" doi:10.1515/9781501508684-004
    Only one option available currently, included for future development.
        
    Returns
    -------
    log10fO2 value of NNO buffer as <class 'mpfr'>

    """

    model = models.loc["NNObuffer","option"]

    P = PT['P']
    T_K = PT["T"]+273.15
    if model == "Frost91": # Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance" doi:10.1515/9781501508684-004
        buffer = (-24930/T_K + 9.36 + 0.046*(P-1.0)/T_K)
    return buffer

def FMQ(PT,models=default_models):
    """ 
    Value of FMQ buffer

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "FMQbuffer" and column label of "option"

    Model options
    -------------
    default: 'Frost91' Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance" doi:10.1515/9781501508684-004
    Other options:
    'ONeill87' O'Neill (1897) AmMin 72(1-2):67-75
        
    Returns
    -------
    log10fO2 value of FMQ buffer as <class 'mpfr'>

    """

    model = models.loc["FMQbuffer","option"]

    P = PT['P']
    T_K = PT["T"]+273.15
    if model == "Frost91": # Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance" doi:10.1515/9781501508684-004
        buffer = (-25096.3/T_K + 8.735 + 0.11*(P-1.0)/T_K)
    elif model == "ONeill87": # O'Neill (1897) AmMin 72(1-2):67-75
        buffer = (8.58 - (25050/T_K))
    return buffer


# terms for different equations

def FefO2_KC91_Eq7_terms(PT,melt_wf,models=default_models): # terms for Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    # ln(XFe2O3/XFeO) = alnfO2 + (b/T) + c + sum(dX) + e[1 - (T0/T) = ln(T/T0)] + f(P/T) + g(((T-T0)P)/T) + h(P2/T)
    # ln(XFe2O3/XFeO) = alnfO2 + B
    # terms
    a = 0.196
    # sum(dX)
    # mole frations in the melt based on oxide components (all Fe as FeO) with no volatiles
    melt_comp = mg.melt_mole_fraction(melt_wf,models,"no","no")
    DAl = -2.243
    DFe = -1.828
    DCa = 3.201
    DNa = 5.854
    DK = 6.215
    D4X = DAl*melt_comp["Al2O3"] + DFe*melt_comp["FeOT"] + DCa*melt_comp["CaO"] + DNa*melt_comp["Na2O"] + DK*melt_comp["K2O"]
    # PT term
    P = PT['P']
    T_K = PT['T']+273.15
    b = 1.1492e4 # K
    c = -6.675
    e = -3.36
    f = -7.01e-7 # K/Pa
    g = -1.54e-10 # /Pa
    h = 3.85e-17 # K/Pa2
    T0 = 1673.0 # K
    P_Pa = P*1.0e5 # converts bars to pascals
    B = (b/T_K) + c + D4X + e*(1.0 - (T0/T_K) - math.log(T_K/T0)) + f*(P_Pa/T_K) + g*(((T_K-T0)*P_Pa)/T_K) + h*((P_Pa**2.0)/T_K) 
    return a, B

def FefO2_KC91_EqA_terms(PT,melt_wf,models=default_models): # terms for Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    # XFeO1.5/XFeO = (KD1*fO2**0.25 + 2y*KD2*KD1**2y*fO2**0.5y)/(1 + (1-2y)KD2*KD1**2y*fO2**0.5y)
    KD2 = 0.4
    y = 0.3
    # compositional term
    # mole frations in the melt based on oxide components (all Fe as FeO) with no volatiles
    melt_comp = mg.melt_mole_fraction(melt_wf,models,"no","no")
    DWAl = 39.86e3             #J
    DWCa = -62.52e3            #J
    DWNa = -102.0e3            #J
    DWK = -119.0e3             #J
    D4X = DWAl*melt_comp["Al2O3"]+DWCa*melt_comp["CaO"]+DWNa*melt_comp["Na2O"]+DWK*melt_comp["K2O"]
    # KD1
    T_K = PT['T']+273.15
    P = PT['P']
    DH = -106.2e3               #J
    DS = -55.10                 #J/K
    DCp = 31.86                 #J/K
    DV = 7.42e-6                #m3
    DVdot = 1.63e-9             #m3/K
    DVdash = -8.16e-16          #m3/Pa
    T0 = 1673.0                 # K
    P0 = 1.0e5                  # Pa 
    R = 8.3144598               # J/K/mol
    P_Pa = P*1.0e5
    if models.loc["high precision","option"] == "True":
        KD1 = math.exp((-DH/(R*T_K)) + (DS/R) - (DCp/R)*(1.0 - (T0/T_K) - gp.log(T_K/T0)) - (1.0/(R*T_K))*D4X - ((DV*(P_Pa-P0))/(R*T_K)) - ((DVdot*(T_K-T0)*(P_Pa-P0))/(R*T_K)) - (DVdash/(2.0*R*T_K))*pow((P_Pa-P0),2.0))
    else:
        KD1 = math.exp((-DH/(R*T_K)) + (DS/R) - (DCp/R)*(1.0 - (T0/T_K) - math.log(T_K/T0)) - (1.0/(R*T_K))*D4X - ((DV*(P_Pa-P0))/(R*T_K)) - ((DVdot*(T_K-T0)*(P_Pa-P0))/(R*T_K)) - (DVdash/(2.0*R*T_K))*pow((P_Pa-P0),2.0))
    return KD1, KD2, y

def FefO2_ONeill18_terms(PT,melt_wf,models=default_models): # Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162 doi:10.1016/j.epsl.2018.10.0020012-821X
    # 1n(Fe3Fe2) = a*DFMQ + B
    a = 0.125
    # mole fractions on a single cation basis in the melt based on oxide components (all Fe as FeO) with no volatiles
    melt_comp = mg.melt_cation_proportion(melt_wf,"volatiles","Fe speciation")
    B = - 1.36 + 2.4*melt_comp["Ca"] + 2.0*melt_comp["Na"] + 3.7*melt_comp["K"] - 2.4*melt_comp["P"]
    # FMQ
    T_K = PT['T']+273.15
    FMQ = 8.58 - (25050/T_K) # O'Neill (1987)
    return a, B, FMQ

def FefO2_Borisov18_terms(PT,melt_wf,models=default_models): # Eq. (4) from Borisov et al. (2018) CMP 173:98 doi:10.1007/s00410-018-1524-8
    T_K = PT['T']+273.15
    # Borisov et al. (2018) CMP 173
    a = 0.207
    # melt mole fraction with no volatiles and all Fe as FeOT
    melt_comp = mg.melt_mole_fraction(melt_wf,models,"no","no")  
    B = (4633.3/T_K -0.445*melt_comp["SiO2"] - 0.900*melt_comp["TiO2"] + 1.532*melt_comp["MgO"] + 0.314*melt_comp["CaO"] + 2.030*melt_comp["Na2O"] + 3.355*melt_comp["K2O"] - 4.851*melt_comp["P2O5"] - 3.081*melt_comp["SiO2"]*melt_comp["Al2O3"] -  4.370*melt_comp["SiO2"]*melt_comp["MgO"] - 1.852)
    return a, B

def fO22Fe3FeT(fO2,PT,melt_wf,models=default_models): # converting fO2 to Fe3/FeT
    """ 
    Fe3+/FeT in the melt from fO2

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "fO2" and column label of "option"
    
    Returns
    -------
    Fe3+/FeT in the melt <class 'mpfr'>

    Model options
    -------------
    default: 'Kress91A' Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    Other options:
    'Kress91' Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    'ONeill18' Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162 doi:10.1016/j.epsl.2018.10.0020012-821X
    'Borisov18' Eq. (4) from Borisov et al. (2018) CMP 173:98 doi:10.1007/s00410-018-1524-8   

    """
    model = models.loc["fO2","option"]

    T_K = PT['T']+273.15
    
    if model == "Kress91": # Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
        a, PTterm = FefO2_KC91_Eq7_terms(PT,melt_wf,models)
        if models.loc["high precision","option"] == "True":
            lnXFe2O3XFeO = a*gp.log(fO2) + PTterm
            XFe2O3XFeO = gp.exp(lnXFe2O3XFeO)
        else:
            lnXFe2O3XFeO = a*math.log(fO2) + PTterm
            XFe2O3XFeO = math.exp(lnXFe2O3XFeO)
        return (2.0*XFe2O3XFeO)/((2.0*XFe2O3XFeO)+1.0)
    
    elif model == "Kress91A": # Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
        kd1, KD2, y = FefO2_KC91_EqA_terms(PT,melt_wf,models)
        XFeO15XFeO = ((kd1*fO2**0.25)+(2.0*y*KD2*(kd1**(2.0*y))*(fO2**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(kd1**(2.0*y))*(fO2**(0.5*y)))
        return XFeO15XFeO/(XFeO15XFeO+1.0)  
    
    elif model == "ONeill18": # Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162 doi:10.1016/j.epsl.2018.10.0020012-821X
        a,B,FMQ = FefO2_ONeill18_terms(PT,melt_wf,models)
        if models.loc["high precision","option"] == "True":
            DQFM = gp.log10(fO2) - FMQ
        else:
            DQFM = math.log10(fO2) - FMQ
        lnFe3Fe2 = a*DQFM + B
        if models.loc["high precision","option"] == "True":
            Fe3Fe2 =  gp.exp(lnFe3Fe2)
        else:
            Fe3Fe2 =  math.exp(lnFe3Fe2)
        return Fe3Fe2/(Fe3Fe2 + 1.0)
    
    elif model == "Borisov18": # Eq. (4) from Borisov et al. (2018) CMP 173:98 doi:10.1007/s00410-018-1524-8
        a,B = FefO2_Borisov18_terms(PT,melt_wf,models)
        if models.loc["high precision","option"] == "True":
            Fe3Fe2 = 10.**(a*gp.log10(fO2) + B)
        else:
            Fe3Fe2 = 10.**(a*math.log10(fO2) + B)
        return Fe3Fe2/(Fe3Fe2 + 1.0)

def f_O2(PT,melt_wf,models=default_models):
    
    """ 
    fO2 from Fe3+/FeT in the melt

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "fO2" and column label of "option"
    
    Returns
    -------
    fO2 in bars <class 'mpfr'>

    Model options
    -------------
    default: 'Kress91A' Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    Other options:
    'Kress91' Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    'ONeill18' Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162 doi:10.1016/j.epsl.2018.10.0020012-821X
    'Borisov18' Eq. (4) from Borisov et al. (2018) CMP 173:98 doi:10.1007/s00410-018-1524-8   

    """

    model = models.loc["fO2","option"]
    
    def KC91(PT,melt_wf,models):
        a, PTterm = FefO2_KC91_Eq7_terms(PT,melt_wf,models)
        F = 0.5*mg.Fe3Fe2(melt_wf) # XFe2O3/XFeO
        alnfO2 = math.log(F) - PTterm
        fO2 = math.exp(alnfO2/a)
        return fO2
    
    #if model == "yes":
    #    return 10.0**(setup.loc[run,"logfO2"]) 
    
    if model == "Kress91": # Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
        fO2 = KC91(PT,melt_wf,models)
        return fO2
    
    elif model == "Kress91A": # Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
        F = mg.Fe3Fe2(melt_wf) # XFeO1.5/XFeO
        kd1, KD2, y = FefO2_KC91_EqA_terms(PT,melt_wf,models)
            
        def f(y,F,KD2,kd1,x): # KC91A rearranged to equal 0
            f = ((2.0*y - F + 2.0*y*F)*KD2*kd1**(2.0*y)*x**(0.5*y) + kd1*x**0.25 - F)
            return f

        def df(y,F,KD2,kd1,x): # derivative of above
            df = (0.5*y)*(2.0*y - F +2.0*y*F)*KD2*kd1**(2.0*y)*x**((0.5*y)-1.0) + 0.25*kd1*x**-0.75
            return df

        def dx(x):
            diff = abs(0-f(y,F,KD2,kd1,x))
            return diff
 
        def nr(x0, e1):
            delta1 = dx(x0)
            while delta1 > e1:
                x0 = x0 - f(y,F,KD2,kd1,x0)/df(y,F,KD2,kd1,x0)
                delta1 = dx(x0)
            return x0
            
        x0 = KC91(PT,melt_wf,models)
    
        fO2 = nr(x0, 1.e-15)
        return fO2
        
    elif model == "ONeill18": # Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162 doi:10.1016/j.epsl.2018.10.0020012-821X
        F = mg.Fe3Fe2(melt_wf) # Fe3+/Fe2+
        a,B,FMQ = FefO2_ONeill18_terms(PT,melt_wf,models)
        DQFM = (math.log(F) - B)/a
        logfO2 = DQFM + FMQ
        return 10.0**logfO2
    
    #elif model == "S6ST": # remove?!?!
    #    S6T = melt_wf['S6ST']
    #    S62 = mg.overtotal2ratio(S6T)
    #    fO2 = mg.S6S2_2_fO2(S62,melt_wf,run,PT,setup,models)
    #    return fO2
    
    elif model == "Borisov18": # Eq. (4) from Borisov et al. (2018) CMP 173:98 doi:10.1007/s00410-018-1524-8
        F = mg.Fe3Fe2(melt_wf)
        a,B = FefO2_Borisov18_terms(PT,melt_wf,models)
        if models.loc["high precision","option"] == "True":
            fO2 = 10.**((gp.log10(F) - B)/a)
        else:
            fO2 = 10.**((math.log10(F) - B)/a)
        return fO2
    
def S_Nash19_terms(PT): # Nash et al. 2019
    T_K = PT['T']+273.15
    A = 8.
    B = ((8.7436e6)/pow(T_K,2.0)) - (27703.0/T_K) + 20.273
    return A, B

def melt_density(PT,melt_wf,models=default_models): # g/cm3
    """ 
    Melt density

    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"
    
    melt_wf: pandas.DataFrame
        Dataframe of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.
    
    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "density" and column label of "option"

    Model options
    -------------
    default: 'DensityX' DensityX from Iacovino & Till (2019) Volcanica 2(1):1-10 doi:10.30909/vol.02.01.0110 
        
    Returns
    -------
    melt density in g/cm3 <class 'mpfr'>

    """
    if models.loc["density","option"] == "DensityX": # DensityX from Iacovino & Till (2019) Volcanica 2(1):1-10 doi:10.30909/vol.02.01.0110
        melt_comp = mg.melt_normalise_wf(melt_wf,"water","yes")
        P = PT["P"]
        T = PT["T"]
        melt_dx = pd.DataFrame([["sample", melt_comp["SiO2"], melt_comp["TiO2"], melt_comp["Al2O3"], melt_comp["FeO"], melt_comp["Fe2O3"], melt_comp["MgO"], melt_comp["CaO"], melt_comp["Na2O"], melt_comp["K2O"], melt_comp["H2O"],P,T]])
        melt_dx.columns = ["Sample_ID","SiO2","TiO2","Al2O3","FeO","Fe2O3","MgO","CaO","Na2O","K2O","H2O","P","T"]
        output = dx.Density(melt_dx)
        density = output.loc[0,"Density_g_per_cm3"]
    return density

#################################################################################################################################
########################################################### CONSTANTS ###########################################################
#################################################################################################################################

species = [['H',1.008,1.,0.,1.,'','','',''],
        ['C',12.011,'',0.,'','','',''],
        ['O',15.999,-2.,0.,'','','',''],				
        ['Na',22.99,'',0.,'','','','',''],
        ['Mg',24.305,'',0.,'','','','',''],
        ['Al',26.982,'',0.,'','','','',''],
        ['Si',28.085,'',0.,'','','','',''],
        ['P',30.974,'',0.,'','','','',''],
        ['S',32.06,'',0.,'','','','',''],
        ['K',39.098,'',0.,'','','','',''],
        ['Ca',40.078,'',0.,'','','','',''],
        ['Ti',47.867,'',0.,'','','','',''],
        ['Mn',54.938,'',0.,'','','','',''],
        ['Fe',55.845,'',0.,'','','','',''],
        ['SiO2',60.083,0.,2.,4.,1.,2.,'',''],
        ['TiO2',79.865,0.,2.,3.,1.,2.,'',''],
        ['Al2O3',101.961,0.,3.,3.,2.,3.,'',''],
        ['Fe2O3',159.687,0.,3.,3.,2.,3.,'',''],
        ['FeO1.5',79.8435,0.,1.5,3.,1.,1.5,'',''],
        ['FeO',71.844,0.,1.,2.,1.,1.,'',''],
        ['MnO',70.937,0.,1.,2.,1.,1.,'',''],
        ['MgO',40.304,0.,1.,2.,1.,1.,'',''],
        ['CaO',56.077,0.,1.,2.,1.,1.,'',''],
        ['Na2O',61.979,0.,1.,1.,2.,1.,'',''],
        ['K2O',94.195,0.,1.,1.,2.,1.,'',''],
        ['P2O5',141.943,0.,1.,5.,2.,1.,'',''],
        ['OH',17.007,-1.,1.,1.,1.,1.,'',''],
        ['H2O',18.015,0.,1.,1.,1.,1.,647.25,221.1925],
        ['H2S',34.076,0.,0.,1.,1.,1.,373.55,90.0779],
        ['CO',28.01,0.,1.,2.,1.,1.,133.15,34.9571],
        ['CO2',44.009,0.,2.,4.,1.,2.,304.15,73.8659],
        ['CO3',60.008,-2.,3.,4.,1.,3.,'',''],
        ['S2',64.12,0.,0.,'','','',208.15,72.954],
        ['SO2',64.058,0.,2.,4.,1.,2.,430.95,78.7295],
        ['SO3',80.057,0.,3.,6.,1.,3.,'',''],
        ['SO4',96.056,-2.,4.,6,'',4.,'',''],
        ['OCS',60.07,0.,1.,'','','',377.55,65.8612],
        ['O2',31.998,0.,2.,'','','',154.75,50.7638],
        ['H2',2.016,0.,0.,'','','',33.25,12.9696],
        ['CH4',16.043,0.,0.,'','','',191.05,46.4069],
        ['Ar',39.948,'','','','','','',''],
        ['Ne',20.1797,'','','','','','','']]    
# Create the pandas DataFrame
species = pd.DataFrame(species,columns=['species','M','overall_charge','no_O','cat_charge','no_cat','no_an','Tcr','Pcr'])
species = species.set_index('species')