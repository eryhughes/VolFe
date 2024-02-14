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

def check_default_options(models):
    
    def return_options(default,name,models):
        if models == "default_options":
            variable = default
        elif name in models.index:
            variable = models.loc[name,'option']
            if variable == "default":
                variable = default
        else:
            variable = default
        return variable

    # species
    insolubles = return_options('yes','insolubles',models)
    H2S_m = return_options('yes','H2S_m',models)
    species_X = return_options('Ar','species X',models)
    Hspeciation = return_options('none','Hspeciation',models)
    # oxygen fugacity
    fO2 = return_options('Kress91A','fO2',models)
    NNObuffer = return_options('Frost91','NNObuffer',models)
    FMQbuffer = return_options('Frost91','FMQbuffer',models)
    # solubility constants
    CO2 = return_options('Dixon95','carbon dioxide',models)
    H2O = return_options('AllisonDataComp','water',models)
    H2 = return_options('basalt','hydrogen',models)
    S2m = return_options('ONeill21dil','sulfide',models)
    S6p = return_options('ONeill22dil','sulfate',models)
    H2S = return_options('basalt','hydrogen sulfide',models)
    CH4 = return_options('Ardia13','methane',models)
    CO = return_options('basalt','carbon monoxide',models)
    X = return_options('Iacono-Marziano10_Ar_basalt','species X solubility',models)
    Cspec = return_options('basalt','Cspeccomp',models)
    Hspec = return_options('MORB','Hspeccomp',models)
    # saturation conditions
    SCSS = return_options('ONeill21hyd','SCSS',models)
    SCAS = return_options('Zajacz19','SCAS',models)
    sulfur_saturation = return_options('no','sulfur_saturation',models)
    sulfur_is_sat = return_options('no','sulfur_is_sat',models)
    graphite_saturation = return_options('no','graphite_saturation',models)
    # fugacity coefficients
    ideal_gas = return_options('no','ideal_gas',models)
    yCO2 = return_options('SS92','y_CO2',models)
    ySO2 = return_options('SS92_modified','y_SO2',models)
    yH2S = return_options('SS92_modified','y_H2S',models)
    yH2 = return_options('SW64','y_H2',models)
    yO2 = return_options('SS92','y_O2',models)
    yS2 = return_options('SS92','y_S2',models)
    yCO = return_options('SS92','y_CO',models)
    yCH4 = return_options('SS92','y_CH4',models)
    yH2O = return_options('HP91','y_H2O',models)
    yOCS = return_options('SS92','y_OCS',models)
    yX = return_options('ideal','y_X',models)
    # equilibrium constants
    KHOg = return_options('KO97','KOHg',models)
    KHOSg = return_options('KO97','KHOSg',models)
    KOSg = return_options('KO97','KOSg',models)
    KOSg2 = return_options('OM22','KOSg2',models)
    KCOg = return_options('KO97','KCOg',models)
    KCOHg = return_options('KO97','KCOHg',models)
    KOCSg = return_options('Moussallam19','KOCSg',models)
    KCOs = return_options('Holloway92','KCOs',models)
    OCS = return_options('COS','carbonlysulfide',models)
    # degassing calculation
    bulk_composition = return_options('yes','bulk_composition',models)
    starting_P = return_options('bulk','starting_P',models)
    gassing_style = return_options('closed','gassing_style',models)
    gassing_direction = return_options('degas','gassing_direction',models)
    P_variation = return_options('polybaric','P_variation',models)
    eq_Fe = return_options('yes','eq_Fe',models)
    solve_species = return_options('OCS','solve_species',models)
    # other
    density = return_options('DensityX','density',models)
    isotopes = return_options('no','isotopes',models)
    T_variation = return_options('isothermal','T_variation',models)
    crystallisation = return_options('no','cystallisation',models)
    mass_volume = return_options('mass','mass_volume',models)
    calc_sat = return_options('fO2_melt','calc_sat',models)
    bulk_O = return_options('exc S','bulk_O',models)
    error = return_options(0.1,'error',models)
    print_status = return_options('yes','print status',models)
    output_csv = return_options('yes','output csv',models)
    setup = return_options('no','setup',models)

    models = [['insolubles',insolubles],['H2S_m',H2S_m],['species X',species_X],['Hspeciation',Hspeciation],
              ['fO2',fO2],['NNObuffer',NNObuffer],['FMQbuffer',FMQbuffer],
              ['carbon dioxide',CO2],['water',H2O],['hydrogen',H2],['sulfide',S2m],['sulfate',S6p],['hydrogen sulfide',H2S],['methane',CH4],['carbon monoxide',CO],['species X solubility',X],['Cspeccomp',Cspec],['Hspeccomp',Hspec],
              ['SCSS',SCSS],['SCAS',SCAS],['sulfur_saturation',sulfur_saturation],['sulfur_is_sat',sulfur_is_sat],['graphite_saturation',graphite_saturation],
              ['ideal_gas',ideal_gas],['y_CO2',yCO2],['y_SO2',ySO2],['y_H2S',yH2S],['y_H2',yH2],['y_O2',yO2],['y_S2',yS2],['y_CO',yCO],['y_CH4',yCH4],['y_H2O',yH2O],['y_OCS',yOCS],['y_X',yX],
              ['KHOg',KHOg],['KHOSg',KHOSg],['KOSg',KOSg],['KOSg2',KOSg2], ['KCOg',KCOg],['KCOHg',KCOHg],['KOCSg',KOCSg],['KCOs',KCOs],['carbonylsulfide',OCS],
              ['bulk_composition',bulk_composition],['starting_P',starting_P],['gassing_style',gassing_style],['gassing_direction',gassing_direction],['P_variation',P_variation],['eq_Fe',eq_Fe],['solve_species',solve_species],
              ['density',density],['isotopes',isotopes],['T_variation',T_variation],['crystallisation',crystallisation],['mass_volume',mass_volume],['calc_sat',calc_sat],['bulk_O',bulk_O],['error',error],
              ['print status',print_status],['output csv',output_csv],['setup',setup]]
    
    # Create the pandas DataFrame
    models = pd.DataFrame(models, columns=['type', 'option'])
    models = models.set_index('type')
    return models

#################################################################################################################################
##################################################### SOLUBILITY CONSTANTS ######################################################
#################################################################################################################################

###################################
### Solubility constant for H2O ###
###################################
def C_H2O(PT,melt_wf,models):
    
    """ 
    Solubility constant for disolving H2O in the melt.


    Parameters
    ----------
    PT: pandas.DataFrame
        Dataframe of pressure-temperature conditions:
        pressure (bars) as "P"
        temperature ('C) as "T"
        
    melt_wf: pandas.DataFrame
        Dataframe of melt composition, not used but kept in case dependence on melt composition is required.
        
    species: pandas.DataFrame
        Dataframe of species.csv file.
    
    models: pandas.DataFrame
        Typically dataframe of models.csv file
        Minimum requirement is dataframe with row labels of "Hspeciation" and "water" and column label of "option"
        See models.csv for available options

    Returns
    -------
    Solubility constant for H2O as <class 'mpfr'>

    """
    model_speciation = models.loc["Hspeciation","option"]
    model_solubility = models.loc["water","option"]

    if model_speciation == "none": ### C_H2O = (xmH2O)^2/fH2O ### (mole fraction)
        if model_solubility == "ETN-1" or model_solubility == "PST-9": # Fitted to ETN-1 and PST-9 from Lesne et al. (2011) 162:133-151
            C = 4.77591e-6
        elif model_solubility == "VES-9": # Fitted to VES-9 from Lesne et al. (2011) 162:133-151
            C = 5.46061e-6
        elif model_solubility == "rhyolite": # Fitted to Blank et al. (1993) and Silver et al. (1990) datasets
            C = 5.13488743E-06 
        elif model_solubility ==  "evo": # Fitted to match EVo
            C = 2.782e-6
        elif model_solubility == "Lesne11mod": # modified general model from Lesne et al. (2011) 162:133-151
            C = 5.62316e-6
        elif model_solubility == "AllisonDataComp": # fitted to experimental data compilation from Allison et al. (2022) for H2O < 6 wt%
            C = 4.6114e-6
        elif model_solubility == "test": #test
            R_ = 83.144621 # cm3 bar K−1 mol−1
            DV = 12 # cm3/mol
            P0 = 1.0 # bar
            A = 4.6114e-6
            B = -((DV/(R_*T_K))*(P-P0))
            C = A*gp.exp(B)
        elif model_solubility == "test2": # for Ptot paper
            C = gp.exp(-12.29)
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
            C = A*gp.exp((-DV*(P-P0))/(R*T0))
        elif model_solubility == "alkali basalt": # Lesne et al. (2011) 162:133-151 eqn 31 with added RT term otherwise it will not work
            A = 5.71e-5 # XmH2Om0
            DV = 26.9 # VH2Om0 in cm3/mol
            C = A*gp.exp((-DV*(P-P0))/(R*T0)) 
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
def C_CO3(PT,melt_wf,models): ### C_CO2,T = xmCO2,T/fCO2 ### (mole fraction) ***except Shishkina14 - wmCO2 ppm***
    
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
        Not normally required unless "carbon dioxide" option requires melt composition.
        
    species: pandas.DataFrame
        Dataframe of species.csv file.
    
    models: pandas.DataFrame
        Typically dataframe of models.csv file
        Minimum requirement is dataframe with row label of "carbon dioxide" and column label of "option"
        See models.csv for available options

    Returns
    -------
    Solubility constant for CO2 as <class 'mpfr'>

    """
    model = models.loc["carbon dioxide","option"]

    P = PT['P']
    T_K = PT['T']+273.15     
    # Calculate cation proportions with no volatiles but correct Fe speciation if available (a la Dixon 1997)
    melt_comp = mg.melt_cation_proportion(melt_wf,"no","no")

    R = 83.15
    T0 = 1473.15 # K
    PI = -6.5*(melt_comp["Si"]+melt_comp["Al"]) + 20.17*(melt_comp["Ca"] + 0.8*melt_comp["K"] + 0.7*melt_comp["Na"] + 0.4*melt_comp["Mg"] + 0.4*melt_comp["FeT"]) # Dixon (1997) Am. Min. 82:368-378
    PI_ = (melt_comp["Ca"] + 0.8*melt_comp["K"] + 0.7*melt_comp["Na"] + 0.4*melt_comp["Mg"] + 0.4*melt_comp["FeT"])/(melt_comp["Si"] + melt_comp["Al"]) # Shishkina et al. (2014) Chem. Geol. 388:112-129
    DH = -13.1 # kJ/mol # Lesne et al. (2011) CMP 162:153-168 from basanite of Holloway & Blank (1994)
 
    if model == "Dixon95": # Dixon et al. (1995)
        DV = 23. # cm3/mol
        P0 = 1.0 # bar
        A = 3.8e-7
        B = (-DV*(P-P0))/(R*T0)
        C = A*gp.exp(B)
    elif model == "Dixon97": # Compositional dependence from Dixon (1997) Am. Min. 82:368-378 as shown by Witham et al. (2012) [assumes PI-SiO2 relationship in caption of figre 2 is 10.19 instead of 10.9 - if 10.9 is assumed you get negative C_CO3]
        DV = 23 # cm3/mol
        P0 = 1.0 # bar
        A = (7.94e-7)*(PI+0.762)
        B = (-DV*(P-P0))/(R*T0)
        C = A*gp.exp(B)
    elif model == "Lesne11": # Lesne et al. (2011)
        DV = 23 # cm3/mol
        P0 = 1.0 # bar
        A = 7.94e-7*((((871*PI)+93.0)/1000.0)+0.762)
        B = (-DV*(P-P0))/(R*T0)
        C = A*gp.exp(B)
    elif model == "VES-9": # Lesne et al. (2011) CMP 162:153-168
        DV = 31.0 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.10)
        B = -((DV/(R*T_K))*(P-P0)) + (DH/R)*((1.0/T0) - (1.0/T_K))
        C = A*gp.exp(B)
    elif model == "ETN-1": # Lesne et al. (2011) CMP 162:153-168
        DV = 23.0 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.55)
        B = -((DV/(R*T_K))*(P-P0)) + (DH/R)*((1.0/T0) - (1.0/T_K))    
        C = A*gp.exp(B)
    elif model == "PST-9": # Lesne et al. (2011) CMP 162:153-168
        DV = 6.0 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.74)
        B = -((DV/(R*T_K))*(P-P0)) + (DH/R)*((1.0/T0) - (1.0/T_K))
        C = A*gp.exp(B)
    elif model == "Shishkina14": # modified from Shishkina et al. (2014) Chem. Geol. 388:112-129
        A = 1.164 # modified by converting P^A to APyCO2 but only including data up to and including 400 MPa
        B = 6.71*PI_-1.345
        C = A*gp.exp(B)
    elif model == "Sunset Crater": # Sunset Crater from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 16.40 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.67)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "SFVF": # SFVF from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 15.02 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.87)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "Erebus": # Erebus from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = -14.65 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.65)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "Vesuvius": # Vesuvius from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 24.42 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.04)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "Etna": # Etna from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 21.59 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.28)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "Stromboli": # Stromboli from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 14.93 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.68)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "Basanite": # Basanite composition from Holloway and Blank (1994), data from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 21.72 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.32)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "Leucite": # Leucite composition from Thibault and Holloway (1994), data from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 21.53 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-13.36)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "AH3 phonotephrite": # AH3 Phonotephrite composition from Vetere et al. (2014), data from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 30.45 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-13.26)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "N72 Basalt": # N72 basalt composition from Shishkina et al. (2010), data from Allison et al. (2022) CMP 177:40
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = 19.05 # cm3/mol
        P0 = 1000.0 # bar
        A = gp.exp(-14.86)
        B = -((DV/(R_*T_K))*(P-P0))
        C = A*gp.exp(B)
    elif model == "Allison22mod": # modified from Allison et al. (2022) CMP 177:40
        P0 = 1000. # bars
        R_ = 83.144621 # cm3 bar K−1 mol−1
        DV = -3350.650 + 3375.552*(Si+Na) + 2625.385*Ti + 3105.426*Al + 3628.018*Fe2 + 3323.320*(Mg+Ca) + 3795.115*K + 47.004*(Na/(Na+K)) # cm/mol
        lnK0 = -128.365 + 114.098*Si + 92.263*(Ti+Al) + 122.644*(Fe2+Ca+Na) + 111.549*Mg + 138.855*K + 2.239*(Na/(Na+K))
        A = gp.exp(lnK0)
        B = ((-1.*DV)*(P-P0))/(R_*T_K)
        C = A*gp.exp(B)
    elif model == "scaledCsulfate": # O'Neill & Mavrogenes (2022) GCA 334:368-382 eq[12a]
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        lnC = -8.02 + ((21100. + 44000.*Na + 18700.*Mg + 4300.*Al + 44200.*K + 35600.*Ca + 12600.*Mn + 16500.*FeT)/T_K) #CS6+ = [S6+, ppm]/fSO3 
        Csulfate = gp.exp(lnC)*KOSg2(PT,models) # ppm S
        lnCsulfate = math.log(Csulfate)
        lnC=-0.46*lnCsulfate
        A = gp.exp(lnC)
        DV = 23.
        P0 = 1000. # bars
        R_ = 83.144621 # cm3 bar K−1 mol−1
        B = ((-1.*DV)*(P-P0))/(R_*T_K)
        C = A*gp.exp(B)
    elif model == "Blank93": # Blank et al. (1993) - rhyolite - tried for workshop
        DV = 28 # cm3/mol
        P0 = 1.0 # bar
        A = gp.exp(-14.45)
        B = (-DV*(P-P0))/(R*(850.+273.15))
        C = A*gp.exp(B)
    elif model == "Behrens04fit": # Fit to Behrens et al. (2004) - tried for workshop
        DV = 41.8 # cm3/mol
        P0 = 1.0 # bar
        A = gp.exp(-14.2)
        B = (-DV*(P-P0))/(R*(1250.+273.15))
        C = A*gp.exp(B)
    elif model == "dacite": # Fit to Behrens et al. (2004) using Ptot
        DV = 36.5 # cm3/mol
        P0 = 1.0 # bar
        A = gp.exp(-14.3)
        B = (-DV*(P-P0))/(R*(1250.+273.15))
        C = A*gp.exp(B)
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
        A2 =gp.exp(-15.275)
        B2 = (-DV2*(P-P0))/(R*T_K)
        C = A2*gp.exp(B2)
    elif model == "water":
        C = 27.e-6
    return C


########################################
### solubility constant for sulfide ###
########################################
def C_S(PT,melt_wf,models): ### C_S = wmS2-*(fO2/fS2)^0.5 ### (weight ppm)
    
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
        
    species: pandas.DataFrame
        Dataframe of species.csv file.
    
    models: pandas.DataFrame
        Typically dataframe of models.csv file
        Minimum requirement is dataframe with row label of "sulfide" and column label of "option"
        See models.csv for available options

    Returns
    -------
    Solubility constant for CO2 as <class 'mpfr'>

    """
    model = models.loc["sulfide","option"]
    
    T = PT['T'] + 273.15 # T in K
    
    def ONeill21(T,melt_comp):
        lnC = (8.77 - (23590.0/T) + (1673.0/T)*(6.7*(melt_comp["Na"]+melt_comp["K"]) + 4.9*melt_comp["Mg"] + 8.1*melt_comp["Ca"] + 8.9*(melt_comp["FeT"]+melt_comp["Mn"]) + 5.0*melt_comp["Ti"] + 1.8*melt_comp["Al"] - 22.2*melt_comp["Ti"]*(melt_comp["FeT"]+melt_comp["Mn"]) + 7.2*melt_comp["FeT"]*melt_comp["Si"]) - 2.06*math.erf(-7.2*(melt_comp["FeT"]+melt_comp["Mn"])))
        return lnC
    
    if model == "ONeill21":
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","no")
        lnC = ONeill21(T,melt_comp)
     
    if model == "ONeill21dil":
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf,"water","no")        
        lnC = ONeill21(T,melt_comp)

    if model == "ONeill21hyd":
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","no")
        lnC_nondil = ONeill21(T,melt_comp)
        melt_comp = mg.melt_cation_proportion(melt_wf,"water","no")
        lnCdil = ONeill21(T,melt_comp)
        lnCH = (melt_comp["H"]*(6.4 + 12.4*melt_comp["H"] - 20.3*melt_comp["Si"] + 73.0*(melt_comp["Na"]+melt_comp["K"])))
        lnC = lnCdil+lnCH
        
    elif model == "FR54-S1":
        lnC = math.log(((1.3e-4)*10000.))

    C = math.exp(lnC) 
    return C



########################################
### solubility constant for sulfate ###
########################################
def C_SO4(PT,melt_wf,models): ### C_SO4 = wmS6+*(fS2*fO2^3)^-0.5 ### (weight ppm)
    model = models.loc["sulfate","option"]

    T = PT['T'] + 273.15 # T in Kelvin
    P = PT['P'] # P in bars
    slope = 115619.707 # slope for T-dependence for melt inclusion fits
    
    if model == "Nash19": # Nash et al. (2019) EPSL 507:187-198
        S = 1. # S6+/S2- ratio of S6+/S2- of 0.5
        Csulfide = C_S(PT,melt_wf,models)
        A = PT_KCterm(PT,melt_wf,models) # P, T, compositional term from Kress & Carmicheal (1991)
        B = (8743600/T**2) - (27703/T) + 20.273 # temperature dependence from Nash et al. (2019)
        a = 0.196 # alnfO2 from Kress & Carmicheal (1991)
        F = 10**(((math.log10(S))-B)/8.)
        fO2 = math.exp(((math.log(0.5*F))-A)/a)
        Csulfate = (S*Csulfide)/(fO2**2)
    elif model == "S6ST":
        Csulfide = C_S(PT,melt_wf,models)
        fO2 = f_O2(PT,melt_wf,models)
        S6ST_ = melt_wf["S6ST"]
        S = overtotal2ratio(S6ST_)
        Csulfate = (S*Csulfide)/(fO2**2)
    elif model == "Hawaii":
        #Csulfate = gp.exp(30.4) # Using Brounce et al. (2017) dataset at 1200 'C
        Csulfate = math.exp(slope*(1./T) -48.)
    elif model == "Etna":
        Csulfate = math.exp(slope*(1./T) -50.15)
    elif model == "Fuego":
        Csulfate = math.exp(slope*(1./T) -48.5)
    elif model == "Erta Ale":
        Csulfate = math.exp(slope*(1./T) -45.5)
    elif model == "FR54-S1":
        Csulfate = ((67.e6)*10000.)
    elif model == "JdF": # 1100 'C ONLY
        Csulfate = 10.**17.
    elif model == "Boulliung22nP": # Boullioung & Wood (2022) GCA 336:150-164 [eq5] - corrected!
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","no")
        logCS6 = -12.948 + ((15602.*melt_comp["Ca"] + 28649.*melt_comp["Na"] - 9596.*melt_comp["Mg"] + 4194.*melt_comp["Al"] +16016.*melt_comp["Mn"] + 29244.)/T) # wt% S
        Csulfate = (10.**logCS6)*10000. # ppm S
    elif model == "Boulliung22wP": # Boullioung & Wood (2022) GCA 336:150-164 [eq5] - corrected!
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf,"no","no")
        logCS6 = -12.659 + ((3692.*melt_comp["Ca"] - 7592.*melt_comp["Si"] - 13736.*melt_comp["Ti"] + 3762.*melt_comp["Al"] + 34483)/T) - (0.1*P*1.5237)/T # wt% S
        Csulfate = (10.**logCS6)*10000. # ppm S
    elif model in ["ONeill22","ONeill22dil"]: 
        if model == "ONeill22": # O'Neill & Mavrogenes (2022) GCA 334:368-382 eq[12a]
            melt_comp = mg.melt_cation_proportion(melt_wf,"no","yes") # Mole fractions in the melt on cationic lattice (Fe as Fe2 and Fe3) no volatiles   
        elif model == "ONeill22dil": # O'Neill & Mavrogenes (2022) GCA 334:368-382 eq[12a]
            melt_comp = mg.melt_cation_proportion(melt_wf,"water","yes") # Mole fractions in the melt on cationic lattice (Fe as Fe2 and Fe3) includes water
        lnC = -8.02 + ((21100. + 44000.*melt_comp["Na"] + 18700.*melt_comp["Mg"] + 4300.*melt_comp["Al"] + 44200.*melt_comp["K"] + 35600.*melt_comp["Ca"] + 12600.*melt_comp["Mn"] + 16500.*melt_comp["Fe2"])/T) #CS6+ = [S6+, ppm]/fSO3
        Csulfate = gp.exp(lnC)*KOSg2(PT,models) # ppm S
        
    return Csulfate



###################################
### solubility constant for H2S ###
###################################
def C_H2S(PT,melt_wf,models): # C_H2S = wmH2S/fH2S (ppm H2S, fH2S bar)
    model = models.loc["hydrogen sulfide","option"]
    if model == "basalt":
        K = 10.23 # fitted to basalt data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015 ChemGeol 418:104–116
    elif model == "basaltic andesite":
        K = 6.82 # fitted to basaltic andesite data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015 ChemGeol 418:104–116 
    return K



########################################
### solubility constant for hydrogen ###
########################################
def C_H2(PT,melt_wf,models): # C_H2 = wmH2/fH2 (wtppm)
    model = models.loc["hydrogen","option"] 
    
    # Hirchmann et al. (2012) EPSL 345-348:38-48
    R = 83.144598 # bar cm3 /mol /K
    P = PT['P'] # pressure in bars
    T = PT['T'] + 273.15 # T in Kelvin SHOULD BE T0
    P0 = 100.*0.01 # kPa to bars
    if model == "basalt":
        #lnK0 = -11.4 # T0 = 1400 'C, P0 = 100 kPa for mole fraction H2
        lnK0 = -0.9624 # for ppm H2 (fitted in excel)
        DV = 10.6 # cm3/mol
    elif model == "andesite":
        #lnK0 = -10.6 # T0 = 1400 'C, P0 = 100 kPa for mole fraction H2
        lnK0 = -0.1296 # for ppm H2 (fitted in excel)
        DV = 11.3 # cm3/mol
    lnK = lnK0 - (DV*(P-P0))/(R*T) # = ln(XH2/fH2) in ppm/bar
    C = gp.exp(lnK) 
    return C



######################################
### solubility constant for methane ##
######################################
def C_CH4(PT,melt_wf,models): # C_CH4 = wmCH4/fCH4 (ppm)
    model = models.loc["methane","option"]

    if model == "Ardia13": # Ardia et al. (2013) GCA 114:52-71
        R = 83.144598 # bar cm3 /mol /K 
        P = PT['P'] # pressure in bars
        T = PT['T'] + 273.15 # T in Kelvin SHOULD BE T0
        P0 = 100.*0.01 # kPa to bars
        lnK0 = 4.93 # ppm CH4 
        #lnK0 = -7.63 # mole fraction CH4
        DV = 26.85 # cm3/mol
        lnK = lnK0 - (DV*(P-P0))/(R*T) 
        K_ = gp.exp(lnK) # for fCH4 in GPa
        K = 0.0001*K_ # for fCH4 in bars 
    return K



#################################
### solubility constant for CO ##
#################################
def C_CO(PT,melt_wf,models): # C_CO = wmCO/fCO (ppm)
    model = models.loc["carbon monoxide","option"]

    if model == "basalt": # from fitting Armstrong et al. (2015) GCA 171:283-302; Stanley+2014, and Wetzel+13 thermodynamically
        R = 83.144598 # bar cm3 /mol /K 
        P = PT['P'] # pressure in bars
        T = PT['T'] +273.15 # T in Kelvin
        P0 = 1. # in bars
        lnK0 = -2.11 # ppm CO
        DV = 15.20 # cm3/mol
        lnK = lnK0 - (DV*(P-P0))/(R*T) 
        K = gp.exp(lnK) # CO(ppm)/fCO(bars)
    return K


#################################
### solubility constant for X ###
#################################
def C_X(PT,melt_wf,models): # C_X = wmX/fX (ppm)
    species = models.loc["species X","option"]
    model = models.loc["species X solubility","option"]
        
    if species == "Ar":
        if model == "Iacono-Marziano10_Ar_basalt": # Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
            K = 0.0799 # fitted assuming Ar is an ideal gas... i.e. yAr = 1.
        elif model == "Iacono-Marziano10_Ar_rhyolite": # Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
            K = 0.4400 # fitted assuming Ar is an ideal gas... i.e. yAr = 1.
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
    if species == "Ne":
        if model == "Iacono-Marziano10_Ne_basalt": # Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
            K = 0.1504 # fitted assuming Ne is an ideal gas... i.e. yNe = 1.
        elif model == "Iacono-Marziano10_Ne_rhyolite": # Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
            K = 0.8464 # fitted assuming Ne is an ideal gas... i.e. yNe = 1.
    return K


#################################################################################################################################
################################################ solid/liquid volatile saturation ###############################################
#################################################################################################################################

################################################
### sulfate content at anhydrite saturation ###
################################################
def SCAS(PT,melt_wf,models): 
    model = models.loc["SCAS","option"]
    
    T = PT['T'] +273.15
    
    comp = mg.melt_pysulfsat(melt_wf)
 
    if model == "Chowdhury18": # sulfate content (ppm) at anhydrite saturation from Chowdhury & Dasgupta (2018) [T in K]
        # mole fraction melt composition including water but all Fe as FeOT
        melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no")
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
        lnxm_SO4 = a + b*((10.0**4.0)/T) + dX + e*wm_H2OT - f*gp.log(melt_comp["CaO"])                                                                                  
        xm_SO4 = gp.exp(lnxm_SO4) 
        Xm_SO4 = xm_SO4*(xm_SO4 + tot)
        S6CAS = Xm_SO4*species.loc["S","M"]*10000.0
                           
    elif model == "Zajacz19":
        # mole fraction melt composition including water but all Fe as FeOT
        melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no")
        tot = 100.*melt_comp["mol_tot"]
        if melt_comp["Na2O"] + melt_comp["K2O"] + melt_comp["CaO"] >= melt_comp["Al2O3"]:
            P_Rhyo = 3.11*(melt_comp["Na2O"]+melt_comp["K2O"]+melt_comp["CaO"]-melt_comp["Al2O3"])
        else:
            P_Rhyo = 1.54*(melt_comp["Al2O3"]-(melt_comp["Na2O"]+melt_comp["K2O"]+melt_comp["CaO"]))
        NBOT = (2.*melt_comp["Na2O"]+2.*melt_comp["K2O"]+2.*(melt_comp["CaO"]+melt_comp["MgO"]+melt_comp["FeOT"])-melt_comp["Al2O3"]*2.)/(melt_comp["SiO2"]+2.*melt_comp["Al2O3"]) # according to spreadsheet not paper
        P_T = gp.exp(-7890./T)
        P_H2O = melt_comp["H2O"]*(2.09 - 1.65*NBOT) + 0.42*NBOT + 0.23
        P_C = ((P_Rhyo + 251.*melt_comp["CaO"]**2. + 57.*melt_comp["MgO"]**2. + 154.*melt_comp["FeOT"]**2.)/(2.*melt_comp["Al2O3"] + melt_comp["SiO2"]))/(1. + 4.8*NBOT)
        Ksm_SPAnh = gp.exp(1.226*gp.log(P_C*P_T*P_H2O) + 0.079)                         
        Xsm_S = Ksm_SPAnh/melt_comp["CaO"]  
        S6CAS = Xsm_S*tot*32.07*10000.
        
    elif model == "Liu23": # Liu et al. (2023) GCA 349:135-145 eq. 4
        melt_comp = mg.melt_mole_fraction(melt_wf,models,"no","no")
        NBOT = (2.*melt_comp["Na2O"]+2.*melt_comp["K2O"]+2.*(melt_comp["CaO"]+melt_comp["MgO"]+melt_comp["FeOT"])-melt_comp["Al2O3"]*2.)/(melt_comp["SiO2"]+2.*melt_comp["Al2O3"]) 
        melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no")
        lnSCAS = 12.185 - (8541./T) + (1.409*NBOT) + 9.984*melt_comp["CaO"] + melt_wf["H2OT"]*100.
        S6CAS = gp.exp(lnSCAS)
    
    elif model == "pss_CD2019": # Chowdhury & Dasgupta (2019) using PySulfSat by Wieser and Gleeson (2023)
        output = ss.calculate_CD2019_SCAS(df=comp, T_K=T, H2O_Liq=float(100.*melt_wf['H2OT']), Fe3Fet_Liq=None, P_kbar=None)
        S6CAS = float(output["SCAS6_ppm"])
    elif model == "pss_ZT2022": # Zajacz and Tsay (2022) using PySulfSat by Wieser and Gleeson (2023)
        output = ss.calculate_ZT2022_SCAS(df=comp, T_K=T, H2O_Liq=float(100.*melt_wf['H2OT']), Fe3Fet_Liq=None, P_kbar=None)
        S6CAS = float(output["SCAS6_ppm"])
    #elif model == "pss_MK2015": # NOT WORKING: Masotta and Kepler (2015) using PySulfSat by Wieser and Gleeson (2023)
    #    output = ss.calculate_MK2015_SCAS(df=comp, T_K=T, H2O_Liq=float(100.*melt_wf['H2OT']), Fe3Fet_Liq=None, P_kbar=None)
    #    S6CAS = float(output["SCAS6_ppm"])
    
    return S6CAS

###############################################
### sulfide content at sulfide saturation ###
###############################################
def SCSS(PT,melt_wf,models): # sulfide content (ppm) at sulfide saturation from O'Neill (2020) [P in bar,T in K]
    model = models.loc["SCSS","option"]

    P_bar = PT['P']
    T = PT['T'] + 273.15
    Fe3FeT = melt_wf["Fe3FeT"]
    comp = mg.melt_pysulfsat(melt_wf)
    
    if model in ["ONeill21","ONeill21dil","ONeill21hyd"]:
        R = 8.31441
        P = (1.0e-4)*P_bar # pressure in GPa
        D = (137778.0 - 91.66*T + 8.474*T*gp.log(T)) # J/mol
        sulfide_comp = 1.0 # assumes the sulfide is pure FeS (no Ni, Cu, etc.)
        if model == "ONeill21":
            # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) no volatiles
            melt_comp = mg.melt_cation_proportion(melt_wf,"no","yes")
        elif model == "ONeill21dil":
            # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) and water
            melt_comp = mg.melt_cation_proportion(melt_wf,"water","yes")
        elif model == "ONeill21hyd":
            # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) and water
            melt_comp = mg.melt_cation_proportion(melt_wf,"water","yes")
        lnaFeS = gp.log((1.0 - melt_comp["Fe2"])*sulfide_comp)
        lnyFe2 = (((1.0-melt_comp["Fe2"])**2.0)*(28870.0 - 14710.0*melt_comp["Mg"] + 1960.0*melt_comp["Ca"] + 43300.0*melt_comp["Na"] + 95380.0*melt_comp["K"] - 76880.0*melt_comp["Ti"]) + (1.0-melt_comp["Fe2"])*(-62190.0*melt_comp["Si"] + 31520.0*melt_comp["Si"]**2.0))/(R*T)
        lnS = D/(R*T) + gp.log(C_S(PT,melt_wf,models)) - gp.log(melt_comp["Fe2"]) - lnyFe2 + lnaFeS + (-291.0*P + 351.0*gp.erf(P))/T
        SCSS = gp.exp(lnS)  
    
    elif model == "Liu07":
        # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) and water
        melt_comp = mg.melt_cation_proportion(melt_wf,"water","yes")
        MFM = (melt_comp["Na"]+melt_comp["K"]+2.*(melt_comp["Ca"]+melt_comp["Mg"]+melt_comp["Fe2"]))/(melt_comp["Si"]*(melt_comp["Al"]+melt_comp["Fe3"]))
        lnS = 11.35251 - (4454.6/T) - 0.03190*(PT["P"]/T) + 0.71006*gp.log(MFM) - 1.98063*(MFM*melt_comp["H"]) + 0.21867*gp.log(melt_comp["H"]) + 0.36192*gp.log(melt_comp["Fe2"])
        SCSS = gp.exp(lnS)
        
    elif model == "Fortin15":
        # Mole fractions in the melt on cationic lattice (all Fe as FeOT) and water
        melt_comp = mg.melt_cation_proportion(melt_wf,"water","no")
        lnS = 34.784 - (5772.3/T) - 346.54*((0.0001*PT["P"])/T) - 20.393*melt_comp["H"] - 25.499*melt_comp["Si"] - 18.344*melt_comp["Ti"] - 27.381*melt_comp["Al"] - 17.275*melt_comp["Fe"] - 22.398*melt_comp["Mg"] - 20.378*melt_comp["Ca"] - 18.954*melt_comp["Na"] - 32.195*melt_comp["K"]
        SCSS = gp.exp(lnS)
    
    elif model == "Liu21":
        XFeS = 1.
        H2O = melt_wf["H2OT"]*100.
        SCSS = (XFeS*gp.exp(13.88 - (9744./T) - (328.*(0.0001*PT["P"])/T))) + 104.*H2O
    
    elif model == "pss_LiZhang2022": # NOT WORKING Li and Zhang (2022) using PySulfSat by Wieser and Gleeson (2023) - assuming pure FeS
        output = ss.calculate_LiZhang2022_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), H2O_Liq=float(100.*melt_wf['H2OT']), Fe_FeNiCu_Sulf=1., Cu_FeNiCu_Sulf=0., Ni_FeNiCu_Sulf=0., Fe3Fet_Liq=Fe3FeT, logfo2=None,Ni_Liq=None, Cu_Liq=None, Ni_Sulf_init=5, Cu_Sulf_init=5, Fe_Sulf=None, Cu_Sulf=None, Ni_Sulf=None, T_K_transition=True,highT=False, lowT=False)
        SCSS = float(output["SCSS_Tot"])
    elif model == "pss_F2015": # Fortin et al. (2015) using PySulfSat by Wieser and Gleeson (2023)
        output = ss.calculate_F2015_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), H2O_Liq=float(100.*melt_wf['H2OT']), Fe3Fet_Liq=None,Fe_FeNiCu_Sulf=None)
        SCSS = float(output["SCSS2_ppm"])
    elif model == "pss_Liu2021": # Liu et al. (2021) using PySulfSat by Wieser and Gleeson (2023) - assuming pure FeS
        output = ss.calculate_Liu2021_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), Fe_FeNiCu_Sulf=1., H2O_Liq=float(100.*melt_wf['H2OT']), Cu_FeNiCu_Sulf=0., Ni_FeNiCu_Sulf=0.,
Ni_Liq=None, Cu_Liq=None, Fe_Sulf=None, Cu_Sulf=None, Ni_Sulf=None, Ni_Sulf_init=5, Cu_Sulf_init=5, Fe3Fet_Liq=Fe3FeT)
        SCSS = float(output["SCSS2_ppm"])
    elif model == "pss_OM2022": # O'Neill & Mavrogenes (2022) using PySulfSat by Wieser and Gleeson (2023) - assuming pure FeS
        output = ss.calculate_OM2022_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), Fe3Fet_Liq=Fe3FeT, Fe_FeNiCu_Sulf=1., Cu_FeNiCu_Sulf=0., Ni_FeNiCu_Sulf=0., Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None, Ni_Liq=None, Cu_Liq=None, Ni_Sulf_init=5, Cu_Sulf_init=5)
        SCSS = float(output["SCSS2_ppm"])
    elif model == "pss_O2021": # O'Neill (2021) using PySulfSat by Wieser and Gleeson (2023) - assuming pure FeS
        output = ss.calculate_O2021_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), Fe3Fet_Liq=Fe3FeT, Fe_FeNiCu_Sulf=1., Cu_FeNiCu_Sulf=0., Ni_FeNiCu_Sulf=0.,Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None,Ni_Liq=None, Cu_Liq=None, Ni_Sulf_init=5, Cu_Sulf_init=5)
        SCSS = float(output["SCSS2_ppm"])
    elif model == "pss_S2017": # Smythe et al. (2017) using PySulfSat by Wieser and Gleeson (2023) - assuming pure FeS    
        output = ss.calculate_S2017_SCSS(df=comp, T_K=T, P_kbar=(P_bar/1000.), Fe3Fet_Liq=Fe3FeT, Fe_FeNiCu_Sulf=0.999999, Cu_FeNiCu_Sulf=0.000001, Ni_FeNiCu_Sulf=0.,Fe_Sulf=None, Ni_Sulf=None, Cu_Sulf=None, Ni_Liq=None, Cu_Liq=None, H2O_Liq=float(100.*melt_wf['H2OT']),Ni_Sulf_init=5, Cu_Sulf_init=5)
        SCSS = float(output["SCSS2_ppm_ideal_Smythe2017"])
    return SCSS

#################################################################################################################################
#################################### EQUILIBRIUM CONSTANTS FOR HOMOGENEOUS VAPOR EQUILIBRIA #####################################
#################################################################################################################################

# H2 + 0.5O2 = H2O
# K = fH2O/(fH2*(fO2)^0.5)
def KHOg(PT,models):
    model = models.loc["KHOg","option"]

    T_K = PT['T']+273.15
    if model == "KO97":
        K = 10.**((12510.0/T_K)-0.979*(gp.log10(T_K))+0.483)
    return K

# H2O + 0.5S2 = H2S + 0.5O2
# K = (fH2S*(fO2)^0.5)/((fS2^0.5)*fH2O)
def KHOSg(PT,models):
    model = models.loc["KHOSg","option"]

    T_K = PT['T']+273.15
    if model == "KO97": # Kerrick & Ohmoto (1997)
        K = 10.**((-8117.0/T_K)+0.188*gp.log10(T_K)-0.352)
    elif model == "noH2S": # H2S doesn't form in the gas...
        K = 0.
    return K

# 0.5S2 + O2 = SO2
# K = fSO2/((fS2^0.5)*fO2)
def KOSg(PT,models):
    model = models.loc["KOSg","option"]

    T_K = PT['T']+273.15
    if model == "KO97": # Kerrick & Ohmoto (1997)
        K = 10.**((18929.0/T_K)-3.783)
    return K

# 0.5S2 + 1.5O2 = SO3
# K = fSO3/((fS2^0.5)*(fO2^1.5)
def KOSg2(PT,models):
    model = models.loc["KOSg2","option"]

    T_K = PT['T']+273.15
    if model == "OM22": # O'Neill+Mavrogenes2022 from JANF 
        lnK = (55921./T_K) - 25.07 + 0.6465*gp.log(T_K)
        K = gp.exp(lnK) 
    return K

# CO + 0.5O = CO2
# K = fCO2/(fCO*(fO2^0.5))
def KCOg(PT,models):
    model = models.loc["KCOg","option"]

    T_K = PT['T']+273.15
    if model == "KO97": # Kerrick & Ohmoto (1997)
        K = 10.**((14751.0/T_K)-4.535)
    return K

# CH4 + 2O2 = CO2 + 2H2O
# K = (fCO2*(fH2O^2))/(fCH4*(fO2^2))
def KCOHg(PT,models): 
    model = models.loc["KCOHg","option"]

    T_K = PT['T']+273.15
    if model == "KO97": # Kerrick & Ohmoto (1997)
        K = 10.**((41997.0/T_K)+0.719*gp.log10(T_K)-2.404)
    return K

def KOCSg(PT,models): # OCS - depends on system
    reaction = models.loc["carbonylsulfide","option"]
    model = models.loc["KOCSg","option"]

    T = PT['T']+273.15
    if reaction == "COHS":
    # OCS + H2O = CO2 + H2S
    # K = (fCO2*fH2S)/(fOCS*fH2O)
        if models == "EVo": 
            K = gp.exp(0.482 + (16.166e-2/T) + 0.081e-3*T - (5.715e-3/T**2) - 2.224e-1*gp.log(T))
            return K
    if reaction == "COS":
    # 2CO2 + OCS = 3CO + SO2 - 
    # K = (fCO^3*fSO2)/(fCO2^2*fOCS)    
        if model == "Moussallam19": # Moussallam et al. (2019) EPSL 520:260-267
            K = 10.**(9.24403 - (15386.45/T)) # P and f in bars, T in K 
        return K

# Cgraphite + O2 = CO2
def KCOs(PT,models): 
    model = models.loc["KCOs","option"]

    T_K = PT['T']+273.15
    P = PT['P']
    if model == "Holloway92": # Holloway et al. (1992) Eur J. Mineral. 4:105-114 equation (3) KI
        a = 40.07639
        b = -2.5392e-2
        c = 5.27096e-6
        d = 0.0267
        log10K = a + (b*T_K) + (c*T_K**2) + ((P-1)/T_K)
        K = 10.**log10K     
    return K
    
#################################################################################################################################
##################################### EQUILIBRIUM CONSTANTS FOR HOMOGENEOUS MELT EQUILIBRIA #####################################
#################################################################################################################################
       
# H2Omol + O = 2OH
# K = xOH*2/(xH2Omol*xO)
def KHOm(PT,melt_wf,models):
    Hspeccomp = models.loc["Hspeccomp","option"]
    
    T_K = PT['T']+273.15
    
    if Hspeccomp == "rhyolite": # Zhang (1999) Reviews in Geophysics 37(4):493-516
        a = -3120.0
        b = 1.89
        K = gp.exp((a/T_K) + b)
    elif Hspeccomp == "alkali basalt": # average of eqn-15-17 from Lesne et al. (2010) CMP 162:133-151
        a = -8348.0 # VES-9 = -8033.0, ETN-1 = -8300.0, and PST-9 = -8710.0
        b = 7.8108 # VES-9 = 7.4222, ETN-1 = 7.4859, and PEST-9 = 8.5244
        K = gp.exp((a/T_K) + b)
    elif Hspeccomp == "PST-9": # Lesne et al. (2010) CMP 162:133-151 eqn 15
        a = -8710.0
        b = 8.5244
        K = gp.exp((a/T_K) + b)
    elif Hspeccomp == "VES-9": # Lesne et al. (2010) CMP 162:133-151 eqn 16
        a = -8033.0
        b = 7.4222
        K = gp.exp((a/T_K) + b)
    elif Hspeccomp == "ETN-1": # Lesne et al. (2010) CMP 162:133-151 eqn 17
        a = -8300.0
        b = 7.4859
        K = gp.exp((a/T_K) + b)
    elif Hspeccomp == "andesite": # eqn-7 from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143
        a = -3650.0
        b = 2.99
        K = gp.exp((a/T_K) + b)
    elif Hspeccomp == "MORB": # fit to Dixon et al. (1995) data digitised from Lesne et al. (2010) CMP 162:133-151
        a = -2204.99
        b = 1.2600
        K = gp.exp((a/T_K) + b)
    elif Hspeccomp == "albite": # Silver & Stolper (1989) J.Pet 30(3)667-709
        K = 0.17
    return K

def KregH2O(PT,melt_wf,models):
    Hspeccomp = models.loc["Hspeccomp","option"]

    if Hspeccomp == "MORB": # Dixon et al. (1995)
        A = 0.403
        B = 15.333
        C = 10.894
    elif Hspeccomp == "alkali basalt XT": # No T-dependence, hence its the speciation frozen in the glass. Eqn 7-10 from Lesne et al. (2010) CMP 162:133-151 (eqn 7 is wrong)
        # wt% normalised including H2O, all Fe as FeOT
        melt_comp = melt_normalise_wf(melt_wf,"volatiles","Fe speciation")
        Na = melt_comp["Na2O"]*100.
        K = melt_comp["K2O"]*100.
        A = 0.5761*(Na+K) - 0.2884 # eqn-8
        B = -8.9589*(Na+K) + 24.65 # eqn-9
        C = 1.7013*(Na+K) + 9.6481 # eqn-1
    elif Hspeccomp == "alkali basalt": # Includes T-dependence, hence speciation in the melt. Eqn 24-27 from Lesne et al. (2010) CMP 162:133-151
        lnK = ((-2704.4/T_K) + 0.641)
        A = lnK + 49016.0/(R*T_K)
        B = -2153326.51/(R*T_K)
        C = 1.965495217/(R*T_K)
    elif Hspeccomp == "VES-9": # Lesne et al. (2011) 162:133-151 Table 5
        A = 3.139
        B = -29.555
        C = 20.535
    elif Hspeccomp == "ETN-1": # Lesne et al. (2011) 162:133-151 Table 5
        A = 4.128
        B = -45.905
        C = 21.311
    elif Hspeccomp == "PST-9": # Lesne et al. (2011) 162:133-151 Table 5
        A = 2.6
        B = -22.476
        C = 22.295
    elif Hspeccomp == "albite": # Silver & Stolper (1989) J.Pet 30(3)667-709
        A = 0.403
        B = 15.333
        C = 10.894
    
# CO2 + O = CO3
def KCOm(PT,melt_wf,models): # K = 
    Cspeccomp = models.loc["Cspeccomp","option"]

    T_K = PT['T']+273.15
    
    if Cspeccomp == "andesite": # eqn-8 from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143
        a = 8665.0
        b = -5.11
        value = gp.exp((a/T_K) + b)  
    elif Cspeccomp == "dacite": # from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143
        a = 9787.0
        b = -7.69
        value = gp.exp((a/T_K) + b)
    elif Cspeccomp == "basalt": # all oxidised carbon is CO32-
        value = "infinite"
    elif Cspeccomp == "rhyolite": # all oxidised carbon is CO2,mol
        value = 0.
    return value


################################################################################################################################# 
##################################################### FUGACITY COEFFICIENTS #####################################################
#################################################################################################################################

# all fugacity coefficients are assumed to equal 1 below 1 bar.

##########################################
### CORK using Holland & Powell (1991) ###
##########################################

def CORK(PT,p0,a,b,c,d):
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
        
    ln_y = z - 1.0 - gp.log(z-B) - A*gp.log(1.0 + (B/z)) + ln_y_virial
    return gp.exp(ln_y)
    
###########################
### Shi & Saxena (1992) ###
###########################

def lny_SS(PT,Pcr,Tcr):
    P = PT['P']
    T_K = PT['T']+273.15
    Tr = T_K/Tcr
    A, B, C, D, P0, integral0 = Q_SS(PT,Tr,Pcr)
    Pr = P/Pcr
    P0r = P0/Pcr
    integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
    integral_total = integral + integral0
    return integral_total

def Q_SS(PT,Tr,Pcr):
    P = PT['P']
    def Q1000(Pcr):
        Pr_ = 1000.0/Pcr
        P0r_ = 1.0/Pcr
        A0 = 1.0
        B0 = 0.9827e-1*pow(Tr,-1.0) + -0.2709*pow(Tr,-3.0)
        C0 = -0.1030e-2*pow(Tr,-1.5) + 0.1427e-1*pow(Tr,-4.0)
        D0 = 0.0
        return A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
    def Q5000(Pcr):
        Pr_ = 5000.0/Pcr
        P0r_ = 1000.0/Pcr
        A0 = 1.0 + -5.917e-1*pow(Tr,-2.0)
        B0 = 9.122e-2*pow(Tr,-1.0)
        C0 = -1.416e-4*pow(Tr,-2.0) + -2.835e-6*gp.log(Tr)
        D0 = 0.0
        return A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
    if P > 5000.0:
        A = 2.0614 + -2.235*pow(Tr,-2.0) + -3.941e-1*gp.log(Tr)
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
        C = -1.416e-4*pow(Tr,-2.0) + -2.835e-6*gp.log(Tr)
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
    
def y_SS(gas_species,PT,models):
    P = PT['P']
    T_K = PT['T']+273.15
    
    ideal_gas = models.loc["ideal_gas","option"]

    if ideal_gas == "yes":
        return 1.0
    elif P < 1.: # ideal gas below 1 bar
        return 1.
    else:    
        Tcr = species.loc[gas_species,"Tcr"]
        Pcr = species.loc[gas_species,"Pcr"]
        return gp.exp(lny_SS(PT,Pcr,Tcr))/P    

##############################
### for each vapor species ###    
##############################

def y_H2(PT,models):
    P = PT['P']
    T_K = PT['T']+273.15

    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_H2","option"]

    if ideal_gas == "yes" or model == "ideal":
        return 1.0
    elif P < 1.: # ideal gas below 1 bar
        return 1.
    elif model == "SW64": # Shaw & Wones (1964)
        SW1 = gp.exp(-3.8402*pow(T_K,0.125)+0.5410)
        SW2 = gp.exp(-0.1263*pow(T_K,0.5)-15.980)
        SW3 = 300*gp.exp((-0.011901*T_K)-5.941) # NB used a value of -0.011901 instead of -0.11901 as reported to match data in Table 2
        P_atm = 0.986923*P
        ln_y = SW1*P_atm - SW2*pow(P_atm,2.0) + SW3*gp.exp((-P_atm/300.0)-1.0)
        return gp.exp(ln_y)
    elif model == "SS92": # Shi & Saxena (1992) NOT WORKING
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
            integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))            
        elif P > 1000.:
            A = Q1_A_HP + Q2_A_HP*Tr + Q3_A_HP*Tr**(-1.) + Q4_A_HP*Tr**2. + Q5_A_HP*Tr**(-2.) + Q6_A_HP*Tr**3. + Q7_A_HP*Tr**(-3.0) + Q8_A_HP*gp.log(Tr)
            B = Q1_B_HP + Q2_B_HP*Tr + Q3_B_HP*Tr**(-1.) + Q4_B_HP*Tr**2. + Q5_B_HP*Tr**(-2.) + Q6_B_HP*Tr**3. + Q7_B_HP*Tr**(-3.0) + Q8_B_HP*gp.log(Tr)
            C = Q1_C_HP + Q2_C_HP*Tr + Q3_C_HP*Tr**(-1.) + Q4_C_HP*Tr**2. + Q5_C_HP*Tr**(-2.) + Q6_C_HP*Tr**3. + Q7_C_HP*Tr**(-3.0) + Q8_C_HP*gp.log(Tr)
            D = Q1_D_HP + Q2_D_HP*Tr + Q3_D_HP*Tr**(-1.) + Q4_D_HP*Tr**2. + Q5_D_HP*Tr**(-2.) + Q6_D_HP*Tr**3. + Q7_D_HP*Tr**(-3.0) + Q8_D_HP*gp.log(Tr)
            P0 = 1000.0
            Pr_ = 1000.0/Pcr
            P0r_ = 1.0/Pcr
            A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.)
            B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.)
            C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.)
            D0 = 0.0
            integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
        P0r = P0/Pcr
        Pr = P/Pcr
        integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
        return gp.exp(integral + integral0)/P

def y_H2O(PT,models):
    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_H2O","option"]

    P = PT['P']
    T_K = PT['T']+273.15

    if ideal_gas == "yes" or model == "ideal":
        return 1.
    elif P < 1.: # ideal gas below 1 bar
        return 1.
    if model == "HP91": 
    # (T > 673 K only) - using Holland & Powell (1991) CORK
        p0 = 2.00 # in kb
        a = 1113.4 + -0.22291*(T_K - 673.0) + -3.8022e-4*pow((T_K-673.0),2.0) + 1.7791e-7*pow((T_K-673.0),3.0)
        b = 1.465
        c = -3.025650e-2 + -5.343144e-6*T_K
        d = -3.2297554e-3 + 2.2215221e-6*T_K
        y = CORK(PT,p0,a,b,c,d)
        return y

def y_CO2(PT,models):
    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_CO2","option"]

    P = PT['P']
    T_K = PT['T']+273.15

    if ideal_gas == "yes" or model == "ideal":
        return 1.0
    elif P < 1.: # ideal gas below 1 bar
        return 1.
    else:
        if model == "HP91": # use Holland & Powell (1991)
            p0 = 5.00 # in kb
            a = 741.2 + -0.10891*(T_K) + -3.4203e-4*pow(T_K,2.0)
            b = 3.057
            c = -2.26924e-1 + -7.73793e-5*T_K
            d = 1.33790e-2 + -1.1740e-5*T_K
            y = CORK(PT,p0,a,b,c,d)
        elif model == "SS92": # use Shi & Saxena (1992)
            gas_species = "CO2"
            y = y_SS(gas_species,PT,models)
        return y
    
def y_O2(PT,models):
    model = models.loc["y_O2","option"]

    if model == "SS92":
        gas_species = "O2"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y
    
def y_S2(PT,models):
    model = models.loc["y_S2","option"]

    if model == "SS92":
        gas_species = "S2"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y

def y_CO(PT,models):
    model = models.loc["y_CO","option"]

    if model == "SS92":
        gas_species = "CO"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y
    
def y_CH4(PT,models):
    model = models.loc["y_CH4","option"]

    if model == "SS92":
        gas_species = "CH4"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y
    
def y_OCS(PT,models):
    model = models.loc["y_OCS","option"]

    if model == "SS92":
        gas_species = "OCS"
        y = y_SS(gas_species,PT,models)
    elif model == "ideal":
        y = 1.
    return y

def y_X(PT,models): # species X fugacity coefficient
    model = models.loc["y_X","option"]

    if model == "ideal":  # ideal gas
        y = 1.
    return y

#################################################################################
### SO2 and H2S from Shi & Saxena (1992) with option to modify below 500 bars ###
#################################################################################

def y_SO2(PT,models):
    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_SO2","option"]

    P = PT['P']
    T_K = PT['T']+273.15

    gas_species = "SO2"
    if ideal_gas == "yes" or model == "ideal":
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
        A = Q1_A + Q2_A*Tr + Q3_A*Tr**(-1.) + Q4_A*Tr**2. + Q5_A*Tr**(-2.) + Q6_A*Tr**3. + Q7_A*Tr**(-3.0) + Q8_A*gp.log(Tr)
        B = Q1_B + Q2_B*Tr + Q3_B*Tr**(-1.) + Q4_B*Tr**2. + Q5_B*Tr**(-2.) + Q6_B*Tr**3. + Q7_B*Tr**(-3.0) + Q8_B*gp.log(Tr)
        C = Q1_C + Q2_C*Tr + Q3_C*Tr**(-1.) + Q4_C*Tr**2. + Q5_C*Tr**(-2.) + Q6_C*Tr**3. + Q7_C*Tr**(-3.0) + Q8_C*gp.log(Tr)
        D = 0.0
        if P >= 500.: # above 500 bar using Shi and Saxena (1992) as is
            Pr = P/Pcr
            integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
            y = (gp.exp(integral))/P
        elif models.loc["y_SO2","option"] == "SS92": # as is Shi and Saxena (1992)
            Pr = P/Pcr
            integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
            y = (gp.exp(integral))/P
        elif models.loc["y_SO2","option"] == "SS92_modified": # below 500 bar linear fit between the value at 500 bar and y = 1 at 1 bar to avoid weird behaviour...
            Pr = 500./Pcr # calculate y at 500 bar
            integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
            y_500 = (gp.exp(integral))/500.
            y = ((y_500 - 1.)*(P/500.)) + 1. # linear extrapolation to P of interest
        return y       
            
def y_H2S(PT,models):
    ideal_gas = models.loc["ideal_gas","option"]
    model = models.loc["y_H2S","option"]

    P = PT['P']
    T_K = PT['T']+273.15

    gas_species = "H2S"
    if ideal_gas == "yes" or model == "ideal":
        return 1.0
    elif ideal_gas == "no":
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
            if models.loc["y_H2S","option"] == "SS92": # as is Shi and Saxena (1992) 
                A = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*gp.log(Tr)
                B = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*gp.log(Tr)
                C = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*gp.log(Tr)
                D = 0.0
                P0 = 1.0
                integral0 = 0.
            elif models.loc["y_SO2","option"] == "SS92_modified": # below 500 bar linear fit between the value at 500 bar and y = 1 at 1 bar to avoid weird behaviour... 
                P0 = 500.0 # calculate y at 500 bars
                Pr_ = 500.0/Pcr
                P0r_ = 1.0/Pcr
                A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*gp.log(Tr)
                B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*gp.log(Tr)
                C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*gp.log(Tr)
                D0 = 0.0
                integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))            
                y_500 = gp.exp(integral0)/500.
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
            A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*gp.log(Tr)
            B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*gp.log(Tr)
            C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*gp.log(Tr)
            D0 = 0.0
            integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))            
        elif P > 500.:
            A = Q1_A_HP + Q2_A_HP*Tr + Q3_A_HP*Tr**(-1.) + Q4_A_HP*Tr**2. + Q5_A_HP*Tr**(-2.) + Q6_A_HP*Tr**3. + Q7_A_HP*Tr**(-3.0) + Q8_A_HP*gp.log(Tr)
            B = Q1_B_HP + Q2_B_HP*Tr + Q3_B_HP*Tr**(-1.) + Q4_B_HP*Tr**2. + Q5_B_HP*Tr**(-2.) + Q6_B_HP*Tr**3. + Q7_B_HP*Tr**(-3.0) + Q8_B_HP*gp.log(Tr)
            C = Q1_C_HP + Q2_C_HP*Tr + Q3_C_HP*Tr**(-1.) + Q4_C_HP*Tr**2. + Q5_C_HP*Tr**(-2.) + Q6_C_HP*Tr**3. + Q7_C_HP*Tr**(-3.0) + Q8_C_HP*gp.log(Tr)
            D = 0.0
            P0 = 500.0
            Pr_ = 500.0/Pcr
            P0r_ = 1.0/Pcr
            A0 = Q1_A_LP + Q2_A_LP*Tr + Q3_A_LP*Tr**(-1.) + Q4_A_LP*Tr**2. + Q5_A_LP*Tr**(-2.) + Q6_A_LP*Tr**3. + Q7_A_LP*Tr**(-3.0) + Q8_A_LP*gp.log(Tr)
            B0 = Q1_B_LP + Q2_B_LP*Tr + Q3_B_LP*Tr**(-1.) + Q4_B_LP*Tr**2. + Q5_B_LP*Tr**(-2.) + Q6_B_LP*Tr**3. + Q7_B_LP*Tr**(-3.0) + Q8_B_LP*gp.log(Tr)
            C0 = Q1_C_LP + Q2_C_LP*Tr + Q3_C_LP*Tr**(-1.) + Q4_C_LP*Tr**2. + Q5_C_LP*Tr**(-2.) + Q6_C_LP*Tr**3. + Q7_C_LP*Tr**(-3.0) + Q8_C_LP*gp.log(Tr)
            D0 = 0.0
            integral0 = A0*gp.log(Pr_/P0r_) + B0*(Pr_ - P0r_) + (C0/2.0)*(pow(Pr_,2.0) - pow(P0r_,2.0)) + (D0/3.0)*(pow(Pr_,3.0) - pow(P0r_,3.0))
        P0r = P0/Pcr
        Pr = P/Pcr
        integral = A*gp.log(Pr/P0r) + B*(Pr - P0r) + (C/2.0)*(pow(Pr,2.0) - pow(P0r,2.0)) + (D/3.0)*(pow(Pr,3.0) - pow(P0r,3.0))
        return gp.exp(integral + integral0)/P
    
def y_SO3(PT,models):
    return 1.

#######################
### oxygen fugacity ###
#######################

# buffers
def NNO(PT,models):
    model = models.loc["NNObuffer","option"]

    P = PT['P']
    T_K = PT["T"]+273.15
    if model == "Frost91":
        buffer = (-24930/T_K + 9.36 + 0.046*(P-1.0)/T_K) # Frost (1991)
    return buffer

def FMQ(PT,models):
    model = models.loc["FMQbuffer","option"]

    P = PT['P']
    T_K = PT["T"]+273.15
    if model == "Frost91":
        buffer = (-25096.3/T_K + 8.735 + 0.11*(P-1.0)/T_K) # Frost (1991)
    elif model == "ONeill87":
        buffer = (8.58 - (25050/T_K)) # O'Neill (1987)
    return buffer


# terms for different equations

def FefO2_KC91_Eq7_terms(PT,melt_wf,models): # terms for Kress & Carmichael (1991) Equation 7
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

def FefO2_KC91_EqA_terms(PT,melt_wf,models): # terms for Kress & Carmichael (1991) Appendix Equations A1-6
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
    KD1 = math.exp((-DH/(R*T_K)) + (DS/R) - (DCp/R)*(1.0 - (T0/T_K) - gp.log(T_K/T0)) - (1.0/(R*T_K))*D4X - ((DV*(P_Pa-P0))/(R*T_K)) - ((DVdot*(T_K-T0)*(P_Pa-P0))/(R*T_K)) - (DVdash/(2.0*R*T_K))*pow((P_Pa-P0),2.0))
    return KD1, KD2, y

def FefO2_ONeill18_terms(PT,melt_wf,models):
    # O'Neill et al. (2018) EPSL 504:152-162
    # 1n(Fe3Fe2) = a*DFMQ + B
    a = 0.125
    # mole fractions on a single cation basis in the melt based on oxide components (all Fe as FeO) with no volatiles
    melt_comp = mg.melt_cation_proportion(melt_wf,"volatiles","Fe speciation")
    B = - 1.36 + 2.4*melt_comp["Ca"] + 2.0*melt_comp["Na"] + 3.7*melt_comp["K"] - 2.4*melt_comp["P"]
    # FMQ
    T_K = PT['T']+273.15
    FMQ = 8.58 - (25050/T_K) # O'Neill (1987)
    return a, B, FMQ

def FefO2_Borisov18_terms(PT,melt_wf,models):
    T_K = PT['T']+273.15
    # Borisov et al. (2018) CMP 173
    a = 0.207
    # melt mole fraction with no volatiles and all Fe as FeOT
    melt_comp = mg.melt_mole_fraction(melt_wf,models,"no","no")  
    B = (4633.3/T_K -0.445*melt_comp["SiO2"] - 0.900*melt_comp["TiO2"] + 1.532*melt_comp["MgO"] + 0.314*melt_comp["CaO"] + 2.030*melt_comp["Na2O"] + 3.355*melt_comp["K2O"] - 4.851*melt_comp["P2O5"] - 3.081*melt_comp["SiO2"]*melt_comp["Al2O3"] -  4.370*melt_comp["SiO2"]*melt_comp["MgO"] - 1.852)
    return a, B

def fO22Fe3FeT(fO2,PT,melt_wf,models): # converting fO2 to Fe3/FeT
    model = models.loc["fO2","option"]

    T_K = PT['T']+273.15
    
    if model == "Kress91":
        a, PTterm = FefO2_KC91_Eq7_terms(PT,melt_wf,models)
        lnXFe2O3XFeO = a*gp.log(fO2) + PTterm
        XFe2O3XFeO = gp.exp(lnXFe2O3XFeO)
        return (2.0*XFe2O3XFeO)/((2.0*XFe2O3XFeO)+1.0)
    
    elif model == "Kress91A": 
        kd1, KD2, y = FefO2_KC91_EqA_terms(PT,melt_wf,models)
        XFeO15XFeO = ((kd1*fO2**0.25)+(2.0*y*KD2*(kd1**(2.0*y))*(fO2**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(kd1**(2.0*y))*(fO2**(0.5*y)))
        return XFeO15XFeO/(XFeO15XFeO+1.0)  
    
    elif model == "ONeill18": # O'Neill et al. (2018) EPSL 504:152-162
        a,B,FMQ = FefO2_ONeill18_terms(PT,melt_wf,models)
        DQFM = gp.log10(fO2) - FMQ
        lnFe3Fe2 = a*DQFM + B
        Fe3Fe2 =  gp.exp(lnFe3Fe2)
        return Fe3Fe2/(Fe3Fe2 + 1.0)
    
    elif model == "Borisov18": # Borisov et al. (2018) CMP 173:
        a,B = FefO2_Borisov18_terms(PT,melt_wf,models)
        Fe3Fe2 = 10.**(a*gp.log10(fO2) + B)
        return Fe3Fe2/(Fe3Fe2 + 1.0)

def f_O2(PT,melt_wf,models):
    model = models.loc["fO2","option"]
    
    def KC91(PT,melt_wf,models):
        a, PTterm = FefO2_KC91_Eq7_terms(PT,melt_wf,models)
        F = 0.5*mg.Fe3Fe2(melt_wf) # XFe2O3/XFeO
        alnfO2 = math.log(F) - PTterm
        fO2 = math.exp(alnfO2/a)
        return fO2
    
    #if model == "yes":
    #    return 10.0**(setup.loc[run,"logfO2"]) 
    
    if model == "Kress91":
        fO2 = KC91(PT,melt_wf,models)
        return fO2
    
    elif model == "Kress91A": 
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
        
    elif model == "ONeill18": # O'Neill et al. (2018) EPSL 504:152-162
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
    
    elif model == "Borisov18": # Borisov et al. (2018) CMP 173
        F = mg.Fe3Fe2(melt_wf)
        a,B = FefO2_Borisov18_terms(PT,melt_wf,models)
        fO2 = 10.**((gp.log10(F) - B)/a)
        return fO2
    
def S_Nash19_terms(PT): # Nash et al. 2019
    T_K = PT['T']+273.15
    A = 8.
    B = ((8.7436e6)/pow(T_K,2.0)) - (27703.0/T_K) + 20.273
    return A, B

# density of the melt in g/cm3 using DensityX (Iacovino & Till 2019 Volcanica 2(1))
def melt_density(PT,melt_wf,models):
    if models.loc["density","option"] == "DensityX":
        melt_comp = mg.melt_normalise_wf(melt_wf,"water","yes")
        P = PT["P"]
        T = PT["T"]
        melt_dx = pd.DataFrame([["sample", melt_comp["SiO2"], melt_comp["TiO2"], melt_comp["Al2O3"], melt_comp["FeO"], melt_comp["Fe2O3"], melt_comp["MgO"], melt_comp["CaO"], melt_comp["Na2O"], melt_comp["K2O"], melt_comp["H2O"],P,T]])
        melt_dx.columns = ["Sample_ID","SiO2","TiO2","Al2O3","FeO","Fe2O3","MgO","CaO","Na2O","K2O","H2O","P","T"]
        output = dx.Density(melt_dx)
        density = output.loc[0,"Density_g_per_cm3"]
    return density

#################################################################################################################################
################################################# ISOTOPE FRACTIONATION FACTORS #################################################
#################################################################################################################################

# beta factors from Richet et al. (1977) fitted to quadratic equation for 600 < T'C < 1300

def beta_gas(PT,element):
    t = 1./(PT["T"]+273.15)
    if element == "S":
        if species == "SO2":
            a, b, c = 4872.56428, 0.76400, 0.99975
        elif species == "S2":
            a, b, c = 1708.22425, -0.76202, 1.00031
        elif species == "OCS":
            a, b, c = 980.75175, 1.74954, 0.99930 
        elif species == "H2S":
            a, b, c = 935.84901, 1.29355, 0.99969
    if models.loc["beta_factors","option"] == "Richet77":
        value = a*t**2 + b*t + c
    return value

def alpha_gas(element,A,B,PT):
    beta_A = beta_gas(PT,element,A)
    beta_B = beta_gas(PT,element,B) 
    result = beta_A/beta_B
    return result

def alpha_H2S_S(PT,models): # Fiege et al. (2015) Chemical Geology equation 8 - H2S fluid and S2- melt
    model = models.loc["alpha_H2S_S","option"]

    if model == "Fiege15":
        T_K = PT["T"] + 273.15
        lna103 = (10.84*((1000./T_K)**2)) - 2.5
        a = gp.exp(lna103/1000.)
    return a

def alpha_SO2_SO4(PT,models): # Fiege et al. (2015) Chemical Geology equation 9 - SO2 fluid and SO4 melt
    model = models.loc["alpha_SO2_SO4","option"]

    if model == "Fiege15":
        T_K = PT["T"] + 273.15
        lna103 = (-0.42*((1000./T_K)**3)) - (2.133*((1000./T_K)**3)) - (0.105*(1000./T_K)) - 0.41
        a = gp.exp(lna103/1000.)
    return a

def alpha_A_B(element,A,B,PT,models):
    if A == "SO2" and B == "H2S":
        a = alpha_gas(element,A,B,PT)
    return a

#################
### constants ###
#################

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