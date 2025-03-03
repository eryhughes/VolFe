# model_dependent_variables.py

import pandas as pd
import numpy as np
import gmpy2 as gp
import math
import densityx as dx
import PySulfSat as ss

# import VolFe.melt_gas as mg

# from VolFe.melt_gas import *
import VolFe.melt_gas as mg

# CONTENTS
# Models
# Solubility constants
# Solid/liquid saturation
# Equilibrium constants
# Fugacity coefficients
# Speciation
# Isotope fractionation factors
# Constants

########################################################################################
# MODELS ###############################################################################
########################################################################################


def make_models_df(models):
    """
    Creates the models options DataFrame

    Parameters
    ----------
    models: list of [str, str]
        A list of lists, where each inner list contains two elements: the model type
        (str)
        and the user-specified option (str) for that model type.

    Returns
    -------
    pandas.DataFrame
        A DataFrame where the first column is 'type', set as the index, and the second
        column
        is 'option', containing the user-specified option. Default options are not
        added if no option is provided.
    """
    # Create the pandas DataFrame
    models = pd.DataFrame(models, columns=["type", "option"])
    models = models.set_index("type")
    return models


# define default models
default_models = [
    ["COH_species", "yes_H2_CO_CH4_melt"],
    ["H2S_m", "True"],
    ["species X", "Ar"],
    ["Hspeciation", "none"],
    ["fO2", "Kress91A"],
    ["NNObuffer", "Frost91"],
    ["FMQbuffer", "Frost91"],
    ["melt composition", "Basalt"],
    ["carbon dioxide", "MORB_Dixon95"],
    ["water", "Basalt_Hughes24"],
    ["hydrogen", "Basalt_Hughes24"],
    ["sulfide", "ONeill21dil"],
    ["sulfate", "ONeill22dil"],
    ["hydrogen sulfide", "Basalt_Hughes24"],
    ["methane", "Basalt_Ardia13"],
    ["carbon monoxide", "Basalt_Hughes24"],
    ["species X solubility", "Ar_Basalt_HughesIP"],
    ["Cspeccomp", "Basalt"],
    ["Hspeccomp", "MORB_HughesIP"],
    ["SCSS", "ONeill21hyd"],
    ["SCAS", "Zajacz19_pss"],
    ["sulfur_saturation", "False"],
    ["sulfur_is_sat", "no"],
    ["graphite_saturation", "False"],
    ["ideal_gas", "False"],
    ["y_CO2", "Shi92"],
    ["y_SO2", "Shi92_Hughes23"],
    ["y_H2S", "Shi92_Hughes24"],
    ["y_H2", "Shaw64"],
    ["y_O2", "Shi92"],
    ["y_S2", "Shi92"],
    ["y_CO", "Shi92"],
    ["y_CH4", "Shi92"],
    ["y_H2O", "Holland91"],
    ["y_OCS", "Shi92"],
    ["y_X", "ideal"],
    ["KHOg", "Ohmoto97"],
    ["KHOSg", "Ohmoto97"],
    ["KOSg", "Ohmoto97"],
    ["KOSg2", "ONeill22"],
    ["KCOg", "Ohmoto97"],
    ["KCOHg", "Ohmoto97"],
    ["KOCSg", "Moussallam19"],
    ["KCOs", "Holloway92"],
    ["carbonylsulfide", "COS"],
    ["bulk_composition", "melt-only"],
    ["starting_P", "Pvsat"],
    ["gassing_style", "closed"],
    ["gassing_direction", "degas"],
    ["P_variation", "polybaric"],
    ["eq_Fe", "yes"],
    ["solve_species", "OCS"],
    ["beta_factors", "Richet77"],
    ["alpha_H_CH4v_CH4m", "no fractionation"],
    ["alpha_H_H2v_H2m", "no fractionation"],
    ["alpha_H_H2Ov_OHmm", "Rust04"],
    ["alpha_H_H2Ov_H2Om", "Rust04"],
    ["alpha_H_H2Sv_H2Sm", "no fractionation"],
    ["alpha_C_CH4v_CH4m", "no fractionation"],
    ["alpha_C_COv_COm", "no fractionation"],
    ["alpha_C_CO2v_CO2T", "LeePP"],
    ["alpha_C_CO2v_CO2m", "Blank93"],
    ["alpha_C_CO2v_CO32mm", "LeePP"],
    ["alpha_S_H2Sv_H2Sm", "no fractionation"],
    ["alpha_SO2_SO4", "Fiege15"],
    ["alpha_H2S_S", "Fiege15"],
    ["density", "DensityX"],
    ["isotopes", "no"],
    ["T_variation", "isothermal"],
    ["crystallisation", "no"],
    ["mass_volume", "mass"],
    ["calc_sat", "fO2_melt"],
    ["bulk_O", "exc_S"],
    ["error", 0.1],
    ["print status", "False"],
    ["output csv", "True"],
    ["setup", "False"],
    ["high precision", "False"],
]
# Create the pandas DataFrame
default_models = pd.DataFrame(default_models, columns=["type", "option"])
default_models = default_models.set_index("type")

# define default models for rhyolite
default_models_rhyolite = [
    ["carbon dioxide", "Rhyolite_Blank93"],
    ["water", "Rhyolite_HughesIP"],
    ["hydrogen", "Andesite_Hughes24"],
    ["hydrogen sulfide", "BasalticAndesite_Hughes24"],
    ["species X solubility", "Ar_Rhyolite_HughesIP"],
    ["Cspeccomp", "Rhyolite"],
    ["Hspeccomp", "Rhyolite_Zhang97"],
]


def check_default_options(models):
    """
    Adds default options to models DataFrame if no user-option is given.

    Parameters
    ----------
    models: pd.DataFrame
        Datarame with user specified model options

    Returns
    -------
    pandas.DataFrame
        A DataFrame where the first column is 'type', set as the index, and
        the second column
        is 'option', containing the user-specified option and default options
        where none are speficied.
    """

    def return_options(default, name, models):
        if name in models.index:
            variable = models.loc[name, "option"]
            if variable == "default":
                variable = default
        else:
            variable = default
        return variable

    # species
    COH_species = return_options(
        default_models.loc["COH_species", "option"], "COH_species", models
    )
    H2S_m = return_options(default_models.loc["H2S_m", "option"], "H2S_m", models)
    species_X = return_options(
        default_models.loc["species X", "option"], "species X", models
    )
    Hspeciation = return_options(
        default_models.loc["Hspeciation", "option"], "Hspeciation", models
    )
    # oxygen fugacity
    fO2 = return_options(default_models.loc["fO2", "option"], "fO2", models)
    NNObuffer = return_options(
        default_models.loc["NNObuffer", "option"], "NNObuffer", models
    )
    FMQbuffer = return_options(
        default_models.loc["FMQbuffer", "option"], "FMQbuffer", models
    )
    # solubility constants
    melt_comp = return_options(
        default_models.loc["melt composition", "option"], "melt composition", models
    )
    S2m = return_options(default_models.loc["sulfide", "option"], "sulfide", models)
    S6p = return_options(default_models.loc["sulfate", "option"], "sulfate", models)
    CH4 = return_options(default_models.loc["methane", "option"], "methane", models)
    CO = return_options(
        default_models.loc["carbon monoxide", "option"], "carbon monoxide", models
    )
    if melt_comp == "Basalt":
        default_models_MC = default_models
    elif melt_comp == "Rhyolite":
        default_models_MC = default_models_rhyolite
    CO2 = return_options(
        default_models_MC.loc["carbon dioxide", "option"], "carbon dioxide", models
    )
    H2O = return_options(default_models_MC.loc["water", "option"], "water", models)
    H2 = return_options(default_models_MC.loc["hydrogen", "option"], "hydrogen", models)
    H2S = return_options(
        default_models_MC.loc["hydrogen sulfide", "option"], "hydrogen sulfide", models
    )
    X = return_options(
        default_models_MC.loc["species X solubility", "option"],
        "species X solubility",
        models,
    )
    Cspec = return_options(
        default_models_MC.loc["Cspeccomp", "option"], "Cspeccomp", models
    )
    Hspec = return_options(
        default_models_MC.loc["Hspeccomp", "option"], "Hspeccomp", models
    )
    # saturation conditions
    SCSS = return_options(default_models.loc["SCSS", "option"], "SCSS", models)
    SCAS = return_options(default_models.loc["SCAS", "option"], "SCAS", models)
    sulfur_saturation = return_options(
        default_models.loc["sulfur_saturation", "option"], "sulfur_saturation", models
    )
    sulfur_is_sat = return_options(
        default_models.loc["sulfur_is_sat", "option"], "sulfur_is_sat", models
    )
    graphite_saturation = return_options(
        default_models.loc["graphite_saturation", "option"],
        "graphite_saturation",
        models,
    )
    # fugacity coefficients
    ideal_gas = return_options(
        default_models.loc["ideal_gas", "option"], "ideal_gas", models
    )
    yCO2 = return_options(default_models.loc["y_CO2", "option"], "y_CO2", models)
    ySO2 = return_options(default_models.loc["y_SO2", "option"], "y_SO2", models)
    yH2S = return_options(default_models.loc["y_H2S", "option"], "y_H2S", models)
    yH2 = return_options(default_models.loc["y_H2", "option"], "y_H2", models)
    yO2 = return_options(default_models.loc["y_O2", "option"], "y_O2", models)
    yS2 = return_options(default_models.loc["y_S2", "option"], "y_S2", models)
    yCO = return_options(default_models.loc["y_CO", "option"], "y_CO", models)
    yCH4 = return_options(default_models.loc["y_CH4", "option"], "y_CH4", models)
    yH2O = return_options(default_models.loc["y_H2O", "option"], "y_H2O", models)
    yOCS = return_options(default_models.loc["y_OCS", "option"], "y_OCS", models)
    yX = return_options(default_models.loc["y_X", "option"], "y_X", models)
    # equilibrium constants
    KHOg = return_options(default_models.loc["KHOg", "option"], "KHOg", models)
    KHOSg = return_options(default_models.loc["KHOSg", "option"], "KHOSg", models)
    KOSg = return_options(default_models.loc["KOSg", "option"], "KOSg", models)
    KOSg2 = return_options(default_models.loc["KOSg2", "option"], "KOSg2", models)
    KCOg = return_options(default_models.loc["KCOg", "option"], "KCOg", models)
    KCOHg = return_options(default_models.loc["KCOHg", "option"], "KCOHg", models)
    KOCSg = return_options(default_models.loc["KOCSg", "option"], "KOCSg", models)
    KCOs = return_options(default_models.loc["KCOs", "option"], "KCOs", models)
    OCS = return_options(
        default_models.loc["carbonylsulfide", "option"], "carbonlysulfide", models
    )
    # degassing calculation
    bulk_composition = return_options(
        default_models.loc["bulk_composition", "option"], "bulk_composition", models
    )
    starting_P = return_options(
        default_models.loc["starting_P", "option"], "starting_P", models
    )
    gassing_style = return_options(
        default_models.loc["gassing_style", "option"], "gassing_style", models
    )
    gassing_direction = return_options(
        default_models.loc["gassing_direction", "option"], "gassing_direction", models
    )
    P_variation = return_options(
        default_models.loc["P_variation", "option"], "P_variation", models
    )
    eq_Fe = return_options(default_models.loc["eq_Fe", "option"], "eq_Fe", models)
    solve_species = return_options(
        default_models.loc["solve_species", "option"], "solve_species", models
    )
    # other
    density = return_options(default_models.loc["density", "option"], "density", models)
    isotopes = return_options(
        default_models.loc["isotopes", "option"], "isotopes", models
    )
    T_variation = return_options(
        default_models.loc["T_variation", "option"], "T_variation", models
    )
    crystallisation = return_options(
        default_models.loc["crystallisation", "option"], "cystallisation", models
    )
    mass_volume = return_options(
        default_models.loc["mass_volume", "option"], "mass_volume", models
    )
    calc_sat = return_options(
        default_models.loc["calc_sat", "option"], "calc_sat", models
    )
    bulk_O = return_options(default_models.loc["bulk_O", "option"], "bulk_O", models)
    error = return_options(default_models.loc["error", "option"], "error", models)
    print_status = return_options(
        default_models.loc["print status", "option"], "print status", models
    )
    output_csv = return_options(
        default_models.loc["output csv", "option"], "output csv", models
    )
    setup = return_options(default_models.loc["setup", "option"], "setup", models)
    precision = return_options(
        default_models.loc["high precision", "option"], "high precision", models
    )

    models = [
        ["COH_species", COH_species],
        ["H2S_m", H2S_m],
        ["species X", species_X],
        ["Hspeciation", Hspeciation],
        ["fO2", fO2],
        ["NNObuffer", NNObuffer],
        ["FMQbuffer", FMQbuffer],
        ["melt composition", melt_comp],
        ["carbon dioxide", CO2],
        ["water", H2O],
        ["hydrogen", H2],
        ["sulfide", S2m],
        ["sulfate", S6p],
        ["hydrogen sulfide", H2S],
        ["methane", CH4],
        ["carbon monoxide", CO],
        ["species X solubility", X],
        ["Cspeccomp", Cspec],
        ["Hspeccomp", Hspec],
        ["SCSS", SCSS],
        ["SCAS", SCAS],
        ["sulfur_saturation", sulfur_saturation],
        ["sulfur_is_sat", sulfur_is_sat],
        ["graphite_saturation", graphite_saturation],
        ["ideal_gas", ideal_gas],
        ["y_CO2", yCO2],
        ["y_SO2", ySO2],
        ["y_H2S", yH2S],
        ["y_H2", yH2],
        ["y_O2", yO2],
        ["y_S2", yS2],
        ["y_CO", yCO],
        ["y_CH4", yCH4],
        ["y_H2O", yH2O],
        ["y_OCS", yOCS],
        ["y_X", yX],
        ["KHOg", KHOg],
        ["KHOSg", KHOSg],
        ["KOSg", KOSg],
        ["KOSg2", KOSg2],
        ["KCOg", KCOg],
        ["KCOHg", KCOHg],
        ["KOCSg", KOCSg],
        ["KCOs", KCOs],
        ["carbonylsulfide", OCS],
        ["bulk_composition", bulk_composition],
        ["starting_P", starting_P],
        ["gassing_style", gassing_style],
        ["gassing_direction", gassing_direction],
        ["P_variation", P_variation],
        ["eq_Fe", eq_Fe],
        ["solve_species", solve_species],
        ["density", density],
        ["isotopes", isotopes],
        ["T_variation", T_variation],
        ["crystallisation", crystallisation],
        ["mass_volume", mass_volume],
        ["calc_sat", calc_sat],
        ["bulk_O", bulk_O],
        ["error", error],
        ["print status", print_status],
        ["output csv", output_csv],
        ["setup", setup],
        ["high precision", precision],
    ]

    # Create the pandas DataFrame
    models = pd.DataFrame(models, columns=["type", "option"])
    models = models.set_index("type")
    return models


def make_df_and_add_model_defaults(models):
    """

    Converts user-provided model configurations (e.g. ['carbon dioxide','MORB_Dixon95'],
    ['hydrogen sulfide','basaltic andesite'] into a structured pandas DataFrame,
    combined with default options for anything not specified


    Parameters
    ----------
    models : list of [str, str]
        A list of lists, where each inner list contains two elements: the model type
        (str)
        and the user-specified option (str) for that model type.

    Returns
    -------
    pandas.DataFrame
        A DataFrame where the first column is 'type', set as the index, and the second
        column is 'option', containing the user-specified option or the default option
        if none is provided.


    Model Parameters and Options
    ----------------------------
    The following parameters can be overridden in models.
    Any parameter can be set to 'setup', in which case the parameter is specified in
    the setup dataframe instead.


    ### Specifying species ###

    COH_species: Specifying what COH species are present in the melt and vapor.
        - 'yes_H2_CO_CH4_melt' [default] Include H2mol (if H present), COmol (if C
        present) and/or CH4mol (if H and C present) as dissolved melt species.
        - 'no_H2_CO_CH4_melt' H2, CO and/or CH4 are insoluble in the melt but they are
        still present in the vapor (H2 in the vapor if H present, CO in the vapor if C
        present, CH4 in the vapor if both H and C present).
        - 'H2O-CO2 only' The only species present in the vapor are H2O and CO2 and in
        the melt are H2OT and CO2T (i.e., no CO, H2, and/or CH4 in the melt or vapor).

    H2S_m: Is H2S a dissolved melt species.
        - 'True' [default] Include H2Smol as a dissolved melt species.
        - 'False' H2Smol is insoluble in the melt.

    species X: Chemical identity of species X, which defines its atomic mass.
        - 'Ar' [default] Species X is argon (i.e., atomic mass of ~40).
        - 'Ne' Species X is Ne (i.e., atomic mass of ~20).
        Other noble gases not currently supported, but we can add them if you get in
        touch!

    Hspeciation:
        - 'none' [default] Oxidised H in the melt only occurs as H2O species (i.e., no
        OH-).
        Only one option available currently, included for future development.


    ### Oxygen fugacity ###

    fO2: Model for parameterisation of relationship between fO2 and Fe3+/FeT
        See function fO2 for options.

    NNObuffer: Model for the parameterisation for the fO2 value of the NNO buffer.
        See function NNO for options.

    FMQbuffer: Model for the parameterisation for the fO2 value of the FMQ buffer.
        See function FMQ for options.


    ### Models for solubility and speciation constants ###

    carbon dioxide: Model for the parameterisation of the CO2T solubility constant.
        See function C_CO3 for options.

    water: Model for the parameterisation for the H2O solubility constant.
        See function C_H2O for options.

    hydrogen: Model for the parameterisation of the H2 solubility constant.
        See function C_H2 for options.

    sulfide: Model for the parameterisation for the *S2- solubility constant (all
    calibrated over wide range of silicate melt compositions).
        See function C_S for options.

    sulfate: Model for the parameterisation of the S6+ solubility constant (all
    calibrated over wide range of silicate melt compositions).
        See function C_SO4 for options.

    hydrogen sulfide: Model for the parameterisation for the H2S solubility constant.
        See function C_H2S for options.

    methane: Model for the parameterisation of the CH4 solubility constant.
        See function C_CH4 for options.

    carbon monoxide: Model for the parameterisation of the CO solubility constant.
        See function C_CO for options.

    species X solubility: Model for the parameterisation of the X solubility constant.
        See function C_X for options.

    Cspeccomp: Model for the parameterisation of the speciation constant for CO2mol and
    CO32- in the melt.
        See function K_COm for options.

    Hspeccomp: Model for the parameterisation of the speciation constant for H2Omol and
    OH- in the melt, either assuming ideal or regular solution models.
        See function K_HOm for options.


    ### Saturation conditions ###

    SCSS: Model for parameterisation of the sulfide content at sulfide saturation
    (S2-CSS).
        See function SCSS for options.

    SCAS: Model for parameterisation of the sulfate content at anhydrite saturation
    (S6+CAS).
        See function SCAS for options.

    sulfur_saturation: Is sulfur allowed to form sulfide or anhydrite if sulfur content
    of the melt reaches saturation levels for these phases?
        - 'False' [default] melt ± vapor are the only phases present - results are
        metastable with respect to sulfide and anhydrite if they could saturate.
        - 'True' If saturation conditions for sulfide or anhydrite are met, melt sulfur
        content reflects this.

    graphite_saturation: Is graphite allowed to form if the carbon content of the melt
    reaches saturation levels for graphite?
        - 'False' [default] melt ± vapor are the only phases present - results are
        metastable with respect to graphite if it could saturate.
        - 'True' If saturation conditions for graphite are met, melt carbon content
        reflects this.

    ### Fugacity coefficients ###

    ideal_gas: Treat all vapor species as ideal gases (i.e., all fugacity coefficients =
      1 at all P).
        - 'False' [default] At least some of the vapor species are not treated as ideal
        gases.
        - 'True' All fugacity coefficients = 1 at all P.

    y_CO2: Model for the parameterisation of the CO2 fugacity coefficient.
        See function y_CO2 for options.

    y_SO2: Model for the parameterisation of the SO2 fugacity coefficient.
        See function y_SO2 for options.

    y_H2S: Model for the parameterisation of the H2S fugacity coefficient.
        See function y_H2S for options.

    y_H2: Model for the parameterisation of the H2 fugacity coefficient.
        See function y_H2 for options.

    y_O2: Model for the parameterisation of the O2 fugacity coefficient.
        See function y_O2 for options.

    y_S2: Model for the parameterisation of the O2 fugacity coefficient.
        See function y_S2 for options.

    y_CO: Model for the parameterisation of the CO fugacity coefficient.
        See function y_CO for options.

    y_CH4: Model for the parameterisation of the CH4 fugacity coefficient.
        See function y_CH4 for options.

    y_H2O: Model for the parameterisation of the H2O fugacity coefficient.
        See function y_H2O for options.

    y_OCS: Model for the parameterisation of the OCS fugacity coefficient.
        See function y_OCS for options.

    y_X: Model for the parameterisation of the X fugacity coefficient.
        See function y_X for options.


    ### Equilibrium constants ###

    KHOg: Model for the parameterisation of the equilibiurm constant for
    H2 + 0.5O2 = H2O.
        See function KHOg for options.

    KHOSg: Model for the parameterisation of the equilibiurm constant for
    0.5S2 + H2O = H2S + 0.5O2.
        See function KHOSg for options.

    KOSg: Model for the parameterisation of the equilibiurm constant for
    0.5S2 + O2 = SO2.
        See function KOSg for options.

    KOSg2: Model for the parameterisation of the equilibiurm constant for
    0.5S2 + 1.5O2 = SO3.
        See function KOSg2 for options.

    KOCg: Model for the parameterisation of the equilibiurm constant for
    CO + 0.5O2 = CO2.
        See function KOCg for options.

    KCOHg: Model for the parameterisation of the equilibiurm constant for
    CH4 + 2O2 = CO2 + 2H2O.
        See function KCOHg for options.

    KOCSg: Model for the parameterisation of the equilibiurm constant for OCS.
        See function KOCSg for options.

    KCOs: Model for the parameterisation of the equilibiurm constant for
    Cgrahite + O2 = CO2.
        See function KCOs for options.

    carbonylsulfide: Reaction equilibrium KOCSg is for.
        - 'COS' [default] 2CO2 + OCS ⇄ 3CO + SO2
        Only one option available currently, included for future development.


    ### Degassing calculation ###

    bulk_composition: Specifying what the inputted melt composition (i.e., dissolved
    volatiles and fO2-estimate) corresponds to for the degassing calculation
        -'melt-only' [default] The inputted melt composition (i.e., dissolved volatiles)
        represents the bulk system - there is no vapor present. The fO2-estimate is
        calculated at Pvsat for this melt composition.
        - 'melt+vapor_wtg' The inputted melt composition (i.e., dissolved volatiles) is
        in equilibrium with a vapor phase. The amount of vapor as weight fraction gas
        (wtg) is specified in the inputs. The bulk system composition will be calculated
        by calculating Pvsat and the vapor composition given the input composition.
        - 'melt+vapor_initialCO2' The inputted melt composition (i.e., dissolved
        volatiles) is in equilibrium with a vapor phase. The initial CO2 content of the
        melt (i.e., before degassing) is specified in the inputs. The bulk system
        composition will be calculated by calculating Pvsat and the vapor composition
        given the input composition. The amount of vapor present is calculated using the
        given initial CO2.

    starting_P: Determing the starting pressure for a degassing calculation.
        - 'Pvsat' [default] Calculation starts at Pvsat for the inputted melt
        composition (i.e., dissolved volatiles).
        - 'set' Calculation starts at the pressure specified in the inputs (using P_bar,
        pressure in bars).

    gassing_style: Does the bulk composition of the system (including oxygen) remain
    constant during the re/degassing
    calculation.
        - 'closed' [default] The bulk composition of the system (inc. oxygen) is
        constant during re/degassing calculation - vapor and melt remain in chemical
        equilibrium throughout.
        - 'open' At each pressure-step, the vapor in equilibrium with the melt is
        removed (or added for regassing), such that the bulk composition of the system
        changes. This does not refer to being buffered in terms of fO2.

    gassing_direction: Is pressure increasing or decreasing from the starting perssure.
        - 'degas' [default] Pressure progressively decreases from starting pressure for
        isothermal, polybaric calculations (i.e., degassing).
        - 'regas' Pressure progressively increases from starting pressure for
        isothermal, polybaric calculations (i.e., regassing).

    P_variation: Is pressure varying during the calculation?
        - 'polybaric' [default] Pressure progressively changes during the calculation.
        Only one option available currently, included for future development.

    T_variation: Is temperature varying during the calculation?
        - 'isothermal' [default] Temperature is constant during the calculation.
        Only one option available currently, included for future development.

    # WHY DID I THINK THIS WAS USEFUL #
    eq_Fe: Does iron in the melt equilibrate with fO2.
        - 'yes' [default] Iron equilibrates with fO2
        Only one option available currently, included for future development.

    solve_species: What species are used to solve the equilibrium equations? This should
    not need to be changed unless the solver is struggling.
        - 'OCS' [default] Guess mole fractions of O2, CO, and S2 in the vapor to solve
        the equilibrium equations.
        - 'OHS' Guess mole fractions of O2, H2, and S2 in the vapor to solve the
        equilibrium equations.
        - 'OCH' Guess mole fractions of O2, CO, and H2 in the vapor to solve the
        equilibrium equations.


    ### Other ###

    density: Model for parameterisation of melt density
        See function melt_density for options.

    setup: Specifies whether model options are specified in the models or setup
    dataframe.
        - 'False' [default] All model options are specified in the models dataframe.
        - 'True' Some of the model options are specified in the setup dataframe.

    print status: Specifies whether some sort of status information during the
    calculation is outputted to let you know progress.
        - 'False' [default] No information about calculation progress is printed.
        - 'True' Some information about calculation progress is printed.

    output csv: Specicies whether a csv of the outputted dataframe is saved at the end
    of the calculation.
        - 'True' [default] csv is outputted
        - 'False' csv is not outputted

    high precision: Is high preicision used for calculations?
        TRUE OR FALSE WHAT PRECISION IS IT USING
        - 'False' [default] normal precision used for calculations
        - 'True' high precision used


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


########################################################################################
# SOLUBILITY CONSTANTS #################################################################
########################################################################################


###################################
# Solubility constant for H2O #####
###################################
def C_H2O(PT, melt_wf, models=default_models):
    """
    Solubility constant for disolving H2O in the melt.

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions:
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with indexes of "Hspeciation" and "water" and
        column label of "option"

    Returns
    -------
    float
        Solubility constant for H2O as <class 'mpfr'>


    Model options for Hspeciation
    -------------
    - "none" [default]
    - Only one option available currently, included for future development.


    Model options for water
    ------------------------
    - 'Basalt_Hughes24' [default] Fig.S2 from Hughes et al. (2024) AmMin 109(3):422-438
    doi:10.2138/am-2023-8739
    - 'Rhyolite_HughesIP' Fig.SX from Hughes et al. (in prep) based on data in Fig. 3 of
    Blank et al. (1993)

    """
    model_speciation = models.loc["Hspeciation", "option"]
    model_solubility = models.loc["water", "option"]

    # C_H2O = (xmH2O)^2/fH2O (mole fraction)
    if model_speciation in [
        "none",
        "none+ideal",
        "none+regular",
    ]:
        # Fig.SX from Hughes et al. (in prep) based on data in Fig. 3 of Blank et al.
        # (1993)
        if model_solubility == "Rhyolite_HughesIP":
            C = 5.3851e-06

        # Fig.S2 from Hughes et al. (2024) based on data compilation from Allison et
        # al. (2022) for basalts with H2O < 6 wt%
        elif model_solubility == "Basalt_Hughes24":
            C = 4.6114e-6

        # WORK IN PROGRESS

        # modified general model from Lesne et al. (2011) 162:133-151
        elif model_solubility == "Lesne11mod":
            C = 5.62316e-6

        # Fitted to ETN-1 and PST-9 from Lesne et al. (2011) 162:133-151
        elif model_solubility == "ETN-1" or model_solubility == "PST-9":
            C = 4.77591e-6

        # Fitted to VES-9 from Lesne et al. (2011) 162:133-151
        elif model_solubility == "VES-9":
            C = 5.46061e-6

        # Fitted to match EVo
        elif model_solubility == "evo":
            C = 2.782e-6

        # test
        elif model_solubility == "test":
            R_ = 83.144621  # cm3 bar K−1 mol−1
            DV = 12  # cm3/mol
            P0 = 1.0  # bar
            A = 4.6114e-6
            T_K = PT["T"] + 273.15
            B = -((DV / (R_ * T_K)) * (PT["P"] - P0))
            if models.loc["high precision", "option"] == "True":
                C = A * gp.exp(B)
            else:
                C = A * math.exp(B)

        # for Ptot paper
        elif model_solubility == "test2":
            if models.loc["high precision", "option"] == "True":
                C = gp.exp(-12.29)
            else:
                C = math.exp(-12.29)

        elif model_solubility == "carbon":
            C = 1.5e-9

        # like 1000 ppm H2O at 730 bar
        elif model_solubility == "test3":
            C = 6.22885e-09

    # C_H2O = xmH2O/fH2O (mole fraction)
    # like AllisonDataComp... I think.
    elif model_speciation == "linear":
        C = 0.00007925494

    # C_H2O = xmH2Omol/fH2O (mole fraction)
    else:
        P = PT["P"]
        P0 = 1.0  # bar
        R = 83.15  # cm3 etc.
        T0 = 1473.15  # K

        # Dixon et al. (1995) - no compositional dependence
        if model_solubility == "Dixon95":
            DV = 12.0
            A = 3.28e-5
            if models.loc["high precision", "option"] == "True":
                C = A * gp.exp((-DV * (P - P0)) / (R * T0))
            else:
                C = A * math.exp((-DV * (P - P0)) / (R * T0))

        # Lesne et al. (2011) 162:133-151 eqn 31 with added RT term otherwise it will
        # not work
        elif model_solubility == "alkali basalt":
            A = 5.71e-5  # XmH2Om0
            DV = 26.9  # VH2Om0 in cm3/mol
            if models.loc["high precision", "option"] == "True":
                C = A * gp.exp((-DV * (P - P0)) / (R * T0))
            else:
                C = A * math.exp((-DV * (P - P0)) / (R * T0))

        # Eq. (9) from Dixon (1997) Am. Min. 82:368-378
        elif model_solubility == "NorthArchBasalt_Dixon97":
            A = (-3.4e-5) - ((1.29e-6) * (mg.melt_comp_ox["SiO2"] * 100.0))  # XmH2Om0
            DV = 12.0  # VH2Om0 in cm3/mol
            if models.loc["high precision", "option"] == "True":
                C = A * gp.exp((-DV * (P - P0)) / (R * T0))
            else:
                C = A * math.exp((-DV * (P - P0)) / (R * T0))

        # WORK IN PROGRESS

        # Fitted to ETN-1 and VES-9 Xm_H2Omol calculated at 1200 'C data from Lesne et
        # al. (2011) 162:133-151
        elif model_solubility == "ETN-1":
            C = 3.3989655e-6

        # Fitted to ETN-1 and VES-9 Xm_H2Omol calculated at 1200 'C data from Lesne et
        # al. (2011) 162:133-151
        elif model_solubility == "VES-9":
            C = 3.3989655e-6

        # Fitted to PST-9 Xm_H2Omol calculated at 1200 'C data from Lesne et al. (2011)
        # 162:133-151
        elif model_solubility == "PST-9":
            C = 1.7022269e-6

    return C


##############################################
# Solubility constant for carbon dioxide #####
##############################################
# C_CO2,T = xmCO2,T/fCO2 (mole fraction) ***except Shishkina14 - wmCO2 ppm***
def C_CO3(PT, melt_wf, models=default_models):
    """
    Solubility constant for disolving CO2 as CO2,T (all oxidised carbon, i.e., CO2mol
    and CO32-, as CO2,T) in the melt: C_CO2,T = xmCO2,T/fCO2 (mole fraction/bar)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "carbon dioxide" and column label
        of "option"

    Returns
    -------
    float
        Solubility constant for CO2 in mole fraction/bar


    Model options for 'carbon dioxide'
    -------------
    - 'MORB_Dixon95' [default] Bullet (5) of summary from Dixon et al. (1995) JPet
    36(6):1607-1631 https://doi.org/10.1093/oxfordjournals.petrology.a037267
    - 'Basalt_Dixon97' Eq. (7) from Dixon et al. (1997) AmMin 82(3-4)368-378
    https://doi.org/10.2138/am-1997-3-415
    - 'NorthArchBasalt_Dixon97' Eq. (8) from Dixon et al. (1997) AmMin 82(3-4)368-378
    https://doi.org/10.2138/am-1997-3-415
    - 'Basalt_Lesne11' Eq. (25,26) from Lesne et al. (2011) CMP 162:153-168
    https://doi.org/10.1007/s00410-010-0585-0
    - 'VesuviusAlkaliBasalt_Lesne11' VES-9 in Table 4 from Lesne et al. (2011) CMP
    162:153-168 https://doi.org/10.1007/s00410-010-0585-0
    - 'EtnaAlkaliBasalt_Lesne11' ETN-1 in Table 4 from Lesne et al. (2011) CMP
    162:153-168 https://doi.org/10.1007/s00410-010-0585-0
    - 'StromboliAlkaliBasalt_Lense11' PST-9 in Table 4 from Lesne et al. (2011)
    CMP 162:153-168 https://doi.org/10.1007/s00410-010-0585-0
    - 'Basanite_Holloway94' Basanite in Table 5 from Holloway and Blank (1994) RiMG
    30:187-230 https://doi.org/10.1515/9781501509674-012
    - 'Leucitite_Thibault94' Leucitite from Thibault & Holloway (1994) CMP 116:216-224
    https://doi.org/10.1007/BF00310701
    - 'Rhyolite_Blank93' Fig.2 caption from Blank et al. (1993) EPSL 119:27-36
    https://doi.org/10.1016/0012-821X(93)90004-S

    """

    model = models.loc["carbon dioxide", "option"]

    P = PT["P"]
    T_K = PT["T"] + 273.15

    # Calculate cation proportions with no volatiles but correct Fe speciation if
    # available (a la Dixon 1997)
    melt_comp = mg.melt_cation_proportion(melt_wf, "no", "no")

    R = 83.15
    T0 = 1473.15  # K

    # Eq. (7) from Dixon (1997) Am. Min. 82:368-378
    PI = -6.5 * (melt_comp["Si"] + melt_comp["Al"]) + 20.17 * (
        melt_comp["Ca"]
        + 0.8 * melt_comp["K"]
        + 0.7 * melt_comp["Na"]
        + 0.4 * melt_comp["Mg"]
        + 0.4 * melt_comp["FeT"]
    )

    # Shishkina et al. (2014) Chem. Geol. 388:112-129
    # PI_ = (
    #    melt_comp["Ca"]
    #    + 0.8 * melt_comp["K"]
    #    + 0.7 * melt_comp["Na"]
    #    + 0.4 * melt_comp["Mg"]
    #    + 0.4 * melt_comp["FeT"]
    # ) / (melt_comp["Si"] + melt_comp["Al"])

    # Bullet (5) of summary from Dixon et al. (1995), which includes values from Pan et
    # al. (1991)
    if model == "MORB_Dixon95":
        DV = 23.0  # cm3/mol
        P0 = 1.0  # bar
        A = 3.8e-7
        B = (-DV * (P - P0)) / (R * T0)
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)

    # Compositional dependence of Eq. (7) from Dixon (1997) Am. Min. 82:368-378 as shown
    # in Eq. (1,5) from Witham et al. (2012)
    elif model == "Basalt_Dixon97":
        DV = 23  # cm3/mol
        P0 = 1.0  # bar
        A = (7.94e-7) * (PI + 0.762)
        B = (-DV * (P - P0)) / (R * T0)
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)

    # Eq. (8) from Dixon (1997) Am. Min. 82:368-378
    elif model == "NorthArchBasalt_Dixon97":
        melt_comp_ox = mg.melt_normalise_wf(melt_wf, "no", "no")
        DV = 23  # cm3/mol
        P0 = 1.0  # bar
        A = (8.70e-6) - ((1.7e-7) * (melt_comp_ox["SiO2"] * 100.0))
        B = (-DV * (P - P0)) / (R * T0)
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)

    # Eq. (25, 26) from Lesne et al. (2011) based on Dixon (1997)
    elif model == "Basalt_Lesne11":
        DV = 25  # cm3/mol ±3
        P0 = 1000.0  # bar
        if models.loc["high precision", "option"] == "True":
            A = gp.exp(0.893 * PI - 15.247)  # Eq. (25)
        else:
            A = math.exp(0.893 * PI - 15.247)  # Eq. (25)
        B = (-DV * (P - P0)) / (R * T0)
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)

    # VES-9 in Table 4 of Lesne et al. (2011) CMP 162:153-168
    elif model == "VesuviusAlkaliBasalt_Lesne11":
        DV = 31.0  # cm3/mol
        P0 = 1000.0  # bar
        if models.loc["high precision", "option"] == "True":
            A = gp.exp(-14.10)  # ±0.03
        else:
            A = math.exp(-14.10)  # ±0.03
        B = -((DV / (R * T_K)) * (P - P0))
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)

    # ETN-1 in Table 4 of Lesne et al. (2011) CMP 162:153-168
    elif model == "EtnaAlkaliBasalt_Lesne11":
        DV = 23.0  # cm3/mol
        P0 = 1000.0  # bar
        if models.loc["high precision", "option"] == "True":
            A = gp.exp(-14.55)  # ±0.00
        else:
            A = math.exp(-14.55)  # ±0.00
        B = -((DV / (R * T_K)) * (P - P0))
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)

    # PST-9 in Table 4 of Lesne et al. (2011) CMP 162:153-168
    elif model == "StromboliAlkaliBasalt_Lense11":
        DV = 6.0  # cm3/mol
        P0 = 1000.0  # bar
        if models.loc["high precision", "option"] == "True":
            A = gp.exp(-14.74)  # ±0.01
        else:
            A = math.exp(-14.74)  # ±0.01
        B = -((DV / (R * T_K)) * (P - P0))
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)
    # elif model == "Shishkina14": # modified from Shishkina et al. (2014) Chem. Geol.
    # 88:112-129
    #    A = 1.164 # modified by converting P^A to APyCO2 but only including data up to
    # and including 400 MPa
    #    B = 6.71*PI_-1.345
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    # elif model == "SunsetCraterAlkaliBasalt_Allison19": # Sunset Crater in Table 4
    # from Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4)
    #    R_ = 83.144621 # cm3 bar K−1 mol−1
    #    DV = 16.40 # cm3/mol
    #    P0 = 1000.0 # bar
    #    if models.loc["high precision","option"] == "True":
    #        A = gp.exp(-14.67)
    #    else:
    #        A = math.exp(-14.67)
    #    B = -((DV/(R_*T_K))*(P-P0))
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    # elif model == "SFVFBasalticAndesite_Allison19": # SFVF in Table 4 from Allison et
    # al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4)
    #    R_ = 83.144621 # cm3 bar K−1 mol−1
    #    DV = 15.02 # cm3/mol
    #    P0 = 1000.0 # bar
    #    if models.loc["high precision","option"] == "True":
    #        A = gp.exp(-14.87)
    #    else:
    #        A = math.exp(-14.87)
    #    B = -((DV/(R_*T_K))*(P-P0))
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    # elif model == "ErebusPhonotephrite_Allison19": # Erebus in Table 4 from Allison
    # et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4)
    #    R_ = 83.144621 # cm3 bar K−1 mol−1
    #    DV = -14.65 # cm3/mol
    #    P0 = 1000.0 # bar
    #    if models.loc["high precision","option"] == "True":
    #        A = gp.exp(-14.65)
    #    else:
    #        A = math.exp(-14.65)
    #    B = -((DV/(R_*T_K))*(P-P0))
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    # elif model == "VesuviusPhonotephrite_Allison19": # Vesuvius in Table 4 from
    # Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4)
    #    R_ = 83.144621 # cm3 bar K−1 mol−1
    #    DV = 24.42 # cm3/mol
    #    P0 = 1000.0 # bar
    #    if models.loc["high precision","option"] == "True":
    #        A = gp.exp(-14.04)
    #    else:
    #        A = math.exp(-14.04)
    #    B = -((DV/(R_*T_K))*(P-P0))
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    # elif model == "EtnaTrachybasalt_Allison19": # Etna in Table 4 from Allison et al.
    # (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4)
    #    R_ = 83.144621 # cm3 bar K−1 mol−1
    #    DV = 21.59 # cm3/mol
    #    P0 = 1000.0 # bar
    #    if models.loc["high precision","option"] == "True":
    #        A = gp.exp(-14.28)
    #    else:
    #        A = math.exp(-14.28)
    #    B = -((DV/(R_*T_K))*(P-P0))
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    # elif model == "StromboliAlkaliBasalt_Allison19": # Stromboli in Table 4 from
    # Allison et al. (2019) CMP 174:58 doi:10.1007/s00410-019-1592-4)
    #   R_ = 83.144621 # cm3 bar K−1 mol−1
    #    DV = 14.93 # cm3/mol
    #    P0 = 1000.0 # bar
    #    if models.loc["high precision","option"] == "True":
    #        A = gp.exp(-14.68)
    #    else:
    #        A = math.exp(-14.68)
    #    B = -((DV/(R_*T_K))*(P-P0))
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    elif (
        model == "Basanite_Holloway94"
    ):  # Basanite in Table 5 from Holloway and Blank (1994)
        R_ = 83.144621  # cm3 bar K−1 mol−1
        DV = 21.72  # cm3/mol ± 1.27
        DH = -13.1  # kJmol ± 13.9
        T0 = 1200.0 + 273.15  # K
        P0 = 1000.0  # bar
        if models.loc["high precision", "option"] == "True":
            A = gp.exp(-14.32)
        else:
            A = math.exp(-14.32)
        B = -((DV / (R_ * T_K)) * (P - P0)) + (DH / R) * ((1.0 / T0) - (1.0 / T_K))
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)
    elif model == "Leucitite_Thibault94":  # Leucitite from Thibault & Holloway (1994)
        R_ = 83.144621  # cm3 bar K−1 mol−1
        DV = 21.53  # cm3/mol ± 0.42
        DH = -28.15  # kJ/mol ±4.24
        T0 = 1200.0 + 273.15  # K
        P0 = 1000.0  # bar
        A = gp.exp(-13.36)
        B = -((DV / (R_ * T_K)) * (P - P0)) + (DH / R) * ((1.0 / T0) - (1.0 / T_K))
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)
    # elif model == "TholeiiteBasalt_Allison22": # N72 basalt in Table 2 from Allison et
    # al. (2022) CMP 177:40, based on experiments from Shishkina et al. (2010)
    #    R_ = 83.144621 # cm3 bar K−1 mol−1
    #    DV = 19.05 # cm3/mol
    #    P0 = 1000.0 # bar
    #    A = gp.exp(-14.86)
    #    B = -((DV/(R_*T_K))*(P-P0))
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    elif model == "Rhyolite_Blank93":  # Fig. 2 caption from Blank et al. (1993)
        R_ = 83.144621  # cm3 bar K−1 mol−1
        DV = 28.0  # cm3/mol ± 2
        DH = -27.2  # kJ/mole ±2.1
        P0 = 1.0  # bar
        T0 = 850.0 + 273.15  # K
        if models.loc["high precision", "option"] == "True":
            A = gp.exp(-14.45)  # ±0.02
        else:
            A = math.exp(-14.45)  # ±0.02
        B = -((DV / (R_ * T_K)) * (P - P0)) + (DH / R) * ((1.0 / T0) - (1.0 / T_K))
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)

    # WORK IN PROGRESS BELOW HERE #
    # elif model == "Phonotephrite_Allison22": # AH3 Phonotephrite in Table 2 from
    # Allison et al. (2022) CMP 177:40, based on experiments from Vetere et al. (2014)
    #    R_ = 83.144621 # cm3 bar K−1 mol−1
    #    DV = 30.45 # cm3/mol
    #    P0 = 1000.0 # bar
    #    if models.loc["high precision","option"] == "True":
    #        A = gp.exp(-13.26)
    #    else:
    #        A = math.exp(-13.26)
    #    B = -((DV/(R_*T_K))*(P-P0))
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    # elif model == "Allison22mod": # modified from Allison et al. (2022) CMP 177:40
    #    P0 = 1000. # bars
    #    R_ = 83.144621 # cm3 bar K−1 mol−1
    #    DV = -3350.650 + 3375.552*(Si+Na) + 2625.385*Ti + 3105.426*Al + 3628.018*Fe2 +
    # 3323.320*(Mg+Ca) + 3795.115*K + 47.004*(Na/(Na+K)) # cm/mol
    #    lnK0 = -128.365 + 114.098*Si + 92.263*(Ti+Al) + 122.644*(Fe2+Ca+Na) +
    # 111.549*Mg + 138.855*K + 2.239*(Na/(Na+K))
    #    if models.loc["high precision","option"] == "True":
    #        A = gp.exp(lnK0)
    #    else:
    #        A = math.exp(lnK0)
    #    B = ((-1.*DV)*(P-P0))/(R_*T_K)
    #    if models.loc["high precision","option"] == "True":
    #        C = A*gp.exp(B)
    #    else:
    #        C = A*math.exp(B)
    elif model == "Behrens04fit":  # Fit to Behrens et al. (2004) - tried for workshop
        DV = 41.8  # cm3/mol
        P0 = 1.0  # bar
        if models.loc["high precision", "option"] == "True":
            A = gp.exp(-14.2)
        else:
            A = math.exp(-14.2)
        B = (-DV * (P - P0)) / (R * (1250.0 + 273.15))
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)
    elif model == "dacite":  # Fit to Behrens et al. (2004) using Ptot
        DV = 36.5  # cm3/mol
        P0 = 1.0  # bar
        if models.loc["high precision", "option"] == "True":
            A = gp.exp(-14.3)
        else:
            A = math.exp(-14.3)
        B = (-DV * (P - P0)) / (R * (1250.0 + 273.15))
        if models.loc["high precision", "option"] == "True":
            C = A * gp.exp(B)
        else:
            C = A * math.exp(B)
    elif model == "test":  # for Ptot paper!!!
        P0 = 1.0  # bar
        # "CO2mol"
        # DV1 = 24.2340838328 # cm3/mol
        # A1 = gp.exp(-14.92978383)
        # B1 = (-DV1*(P-P0))/(R*T_K)
        # "CO32-"
        # DV2 = 8.912862511 # cm3/mol
        # A2 =gp.exp(])
        # B2 = (-DV2*(P-P0))/(R*T_K)
        # C = A1*gp.exp(B1) + A2*gp.exp(B2)
        # "average"
        DV2 = 16.57  # cm3/mol
        if models.loc["high precision", "option"] == "True":
            A2 = gp.exp(-15.275)
        else:
            A2 = math.exp(-15.275)
        B2 = (-DV2 * (P - P0)) / (R * T_K)
        if models.loc["high precision", "option"] == "True":
            C = A2 * gp.exp(B2)
        else:
            C = A2 * math.exp(B2)
    elif model == "water":
        C = 27.0e-6
    return C


########################################
# solubility constant for sulfide ######
########################################
# C_S = wmS2-*(fO2/fS2)^0.5 (weight ppm)
def C_S(PT, melt_wf, models=default_models):
    """
    Solubility constant for disolving S in the melt as S2- in ppmw:
    C_S = wmS2-*(fO2/fS2)^0.5


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "sulfide" and column label of
        "option"

    Returns
    -------
    float
        Solubility constant for S2- as <class 'mpfr'>


    Model options for sulfide
    -------------
    - 'ONeill21dil' [default] Eq. (10.34) inc. H2O dilution from O'Neill (2021) in
    "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    - 'ONeill21' Eq. (10.34) ex. H2O dilution from O'Neill (2021) in "Magma Redox
    Geochemistry" doi:10.1002/9781119473206.ch10
    - 'ONeill21hyd' (hydrous) Eq. (10.34, 10.49) from O'Neill (2021) in "Magma Redox
    Geochemistry" doi:10.1002/9781119473206.ch10
    - 'Boulliung23_eq6' Eq. (6) from Boulliung & Wood (2023) CMP 178:56
    doi:10.1007/s00410-023-02033-9
    - 'Boulliung23_eq7' Eq. (7) from Boulliung & Wood (2023) CMP 178:56
    doi:10.1007/s00410-023-02033-9
    """
    model = models.loc["sulfide", "option"]

    T = PT["T"] + 273.15  # T in K

    # Eq. (10.34) in O'Neill (2021) in "Magma Redox Geochemistry"
    # doi:10.1002/9781119473206.ch10
    def ONeill21(T, melt_comp):
        lnC = (
            8.77
            - (23590.0 / T)
            + (1673.0 / T)
            * (
                6.7 * (melt_comp["Na"] + melt_comp["K"])
                + 4.9 * melt_comp["Mg"]
                + 8.1 * melt_comp["Ca"]
                + 8.9 * (melt_comp["FeT"] + melt_comp["Mn"])
                + 5.0 * melt_comp["Ti"]
                + 1.8 * melt_comp["Al"]
                - 22.2 * melt_comp["Ti"] * (melt_comp["FeT"] + melt_comp["Mn"])
                + 7.2 * melt_comp["FeT"] * melt_comp["Si"]
            )
            - 2.06 * math.erf(-7.2 * (melt_comp["FeT"] + melt_comp["Mn"]))
        )
        return lnC

    # Eq. (10.34) in O'Neill (2021) in "Magma Redox Geochemistry"
    # doi:10.1002/9781119473206.ch10
    if model == "ONeill21":
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(
            melt_wf, "no", "no", molmass="M_ONeill21", majors="majors_ONeill21"
        )
        lnC = ONeill21(T, melt_comp)
        C = math.exp(lnC)

    # Eq. (10.34) O'Neill (2021) in "Magma Redox Geochemistry"
    # doi:10.1002/9781119473206.ch10
    if model == "ONeill21dil":
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(
            melt_wf, "water", "no", molmass="M_ONeill21", majors="majors_ONeill21"
        )
        lnC = ONeill21(T, melt_comp)
        C = math.exp(lnC)

    # Eq. (10.34, 10.49) in O'Neill (2021) in "Magma Redox Geochemistry"
    # doi:10.1002/9781119473206.ch10
    if model == "ONeill21hyd":
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(
            melt_wf, "no", "no", molmass="M_ONeill21", majors="majors_ONeill21"
        )
        melt_comp = mg.melt_cation_proportion(
            melt_wf, "water", "no", molmass="M_ONeill21", majors="majors_ONeill21"
        )
        lnCdil = ONeill21(T, melt_comp)
        lnCH = melt_comp["H"] * (
            6.4
            + 12.4 * melt_comp["H"]
            - 20.3 * melt_comp["Si"]
            + 73.0 * (melt_comp["Na"] + melt_comp["K"])
        )  # Eq. (10.49)
        lnC = lnCdil + lnCH
        C = math.exp(lnC)

    # Eq. (6) from Boulliung, J., Wood, B.J. Sulfur oxidation state and solubility in
    # silicate melts. Contrib Mineral Petrol 178, 56 (2023).
    # https://doi.org/10.1007/s00410-023-02033-9
    if model == "Boulliung23_eq6":
        # Mole fractions in the melt on cationic lattice with no volatiles and Fe
        # speciated
        melt_comp = mg.melt_single_O(
            melt_wf, "no", "yes", molmass="M_Boulliung23", majors="majors_Boulliung23"
        )
        logC = (
            0.338
            + (
                24328.0 * melt_comp["FeO"]
                + 5411.0 * melt_comp["CaO"]
                + 15872.0 * melt_comp["MnO"]
                - 9697.0
            )
            / T
        )
        C = 10.0 ** (logC)

    # Eq. (7) from Boulliung, J., Wood, B.J. Sulfur oxidation state and solubility in
    # silicate melts. Contrib Mineral Petrol 178, 56 (2023).
    # https://doi.org/10.1007/s00410-023-02033-9
    if model == "Boulliung23_eq7":
        # Mole fractions in the melt on cationic lattice with no volatiles and Fe
        # speciated
        melt_comp = mg.melt_single_O(
            melt_wf, "no", "yes", molmass="M_Boulliung23", majors="majors_Boulliung23"
        )
        # 8879.108 used rather than 8879 to match spreadsheet
        logC = (
            0.225
            + (
                25237.0 * melt_comp["FeO"]
                + 5214.0 * melt_comp["CaO"]
                + 12705.0 * melt_comp["MnO"]
                + 19829.0 * melt_comp["K2O"]
                - 1109.0 * melt_comp["SiO2"]
                - 8879.108
            )
            / T
        )
        C = 10.0 ** (logC)

    # elif model == "FR54-S1":
    #    lnC = math.log(((1.3e-4)*10000.))

    return C


########################################
# solubility constant for sulfate ######
########################################
# C_SO4 = wmS6+*(fS2*fO2^3)^-0.5 (weight ppm)
def C_SO4(PT, melt_wf, models=default_models):
    """
    Solubility constant for disolving S6+ in the melt: C_SO4 = wmS6+(fS2*fO2^3)^-0.5
    (in ppmw and bar)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "sulfate" and column label of
        "option"

    Returns
    -------
    float
        Solubility constant for S6+ as <class 'mpfr'>


    Model options for sulfate
    -------------
    - 'ONeill22dil' [default] Eq. (12a) inc. H2O dilution from O'Neill & Mavrogenes
    (2022) GCA 334:368-382 10.1016/j.gca.2022.06.020
    - 'ONeill22' Eq. (12a) without H2O dilution from O'Neill & Mavrogenes (2022) GCA
    334:368-382 doi:10.1016/j.gca.2022.06.020
    - 'Boulliung22nP' Eq. (5) from Boulliung & Wood (2023) GCA 343:420
    doi:10.1016/j.gca.2022.11.025
    - 'Boulliung22wP' Eq. (5) from Boulliung & Wood (2023) GCA 343:420
    doi:10.1016/j.gca.2022.11.025 and Eq. (8) for P from Boulliung & Wood (2022) GCA
    336:150-164 doi:10.1016/j.gca.2022.08.032
    - 'Boulliung23_eq9' Eq. (9) from Boulliung & Wood (2023) CMP 178:56
    doi:10.1007/s00410-023-02033-9
    - 'Boulliung23_eq11' Eq. (11) from Boulliung & Wood (2023) CMP 178:56
    doi:10.1007/s00410-023-02033-9

    """

    model = models.loc["sulfate", "option"]

    T = PT["T"] + 273.15  # T in Kelvin
    P = PT["P"]  # P in bars
    # slope = 115619.707 # slope for T-dependence for melt inclusion fits

    # Eq. (5) in Boullioung & Wood (2022) GCA 336:150-164 - corrected!
    if model in [
        "Boulliung22nP",
        "Boulliung22wP",
    ]:
        # Mole fractions in the melt on cationic lattice (all Fe as FeO) no volatiles
        melt_comp = mg.melt_cation_proportion(melt_wf, "no", "no")
        logCS6 = -12.948 + (
            (
                15602.0 * melt_comp["Ca"]
                + 28649.0 * melt_comp["Na"]
                - 9596.0 * melt_comp["Mg"]
                + 4194.0 * melt_comp["Al"]
                + 16016.0 * melt_comp["Mn"]
                + 29244.0
            )
            / T
        )  # wt% S

        # wt% S, Eq. (8) Boullioung & Wood (2022) GCA 336:150-164
        if model == "Boulliung22wP":
            logCS6 = logCS6 - ((0.1 * ((10.0 * P) - 0.1)) * 1.5237) / T
        Csulfate = (10.0**logCS6) * 10000.0  # ppm S
    elif model in ["ONeill22", "ONeill22dil"]:
        # Eq. (12a) in O'Neill & Mavrogenes (2022) GCA 334:368-382
        if model == "ONeill22":
            # Mole fractions in the melt on cationic lattice (Fe as Fe2 and Fe3)
            # volatiles
            melt_comp = mg.melt_cation_proportion(
                melt_wf, "no", "no", molmass="M_ONeill21", majors="majors_ONeill21"
            )

        # Eq. (12a) in O'Neill & Mavrogenes (2022) GCA 334:368-382
        elif model == "ONeill22dil":
            # Mole fractions in the melt on cationic lattice (Fe as Fe2 and Fe3)
            # includes water
            melt_comp = mg.melt_cation_proportion(
                melt_wf, "water", "no", molmass="M_ONeill21", majors="majors_ONeill21"
            )
        Fe2 = melt_comp["FeT"] * (1.0 - melt_wf["Fe3FeT"])
        lnC = -8.02 + (
            (
                21100.0
                + 44000.0 * melt_comp["Na"]
                + 18700.0 * melt_comp["Mg"]
                + 4300.0 * melt_comp["Al"]
                + 44200.0 * melt_comp["K"]
                + 35600.0 * melt_comp["Ca"]
                + 12600.0 * melt_comp["Mn"]
                + 16500.0 * Fe2
            )
            / T
        )  # CS6+ = [S6+, ppm]/fSO3
        if models.loc["high precision", "option"] == "True":
            Csulfate = gp.exp(lnC) * KOSg2(PT, models)  # ppm S
        else:
            Csulfate = math.exp(lnC) * KOSg2(PT, models)  # ppm S

    # Eq. (9) from Boulliung, J., Wood, B.J. Sulfur oxidation state and solubility in
    # silicate melts. Contrib Mineral Petrol 178, 56 (2023).
    # https://doi.org/10.1007/s00410-023-02033-9
    elif model == "Boulliung23_eq9":
        # Mole fractions in the melt on cationic lattice with no volatiles and Fe
        # speciated
        # Used 29244.299 instead of 292544 to match spreadsheet
        melt_comp = mg.melt_single_O(
            melt_wf, "no", "yes", molmass="M_Boulliung23", majors="majors_Boulliung23"
        )
        logC = (
            -12.948
            + (
                28649.0 * melt_comp["Na2O"]
                + 15602.0 * melt_comp["CaO"]
                + 9496.0 * melt_comp["MgO"]
                + 16016.0 * melt_comp["MnO"]
                + 4194.0 * melt_comp["Al2O3"]
                + 29244.229
            )
            / T
        )
        Csulfate = 10.0 ** (logC)

    # Eq. (11) from Boulliung, J., Wood, B.J. Sulfur oxidation state and solubility in
    # silicate melts. Contrib Mineral Petrol 178, 56 (2023).
    # https://doi.org/10.1007/s00410-023-02033-9
    elif model == "Boulliung23_eq11":
        # Mole fractions in the melt on cationic lattice with no volatiles and Fe
        # speciated
        melt_comp = mg.melt_single_O(
            melt_wf, "no", "yes", molmass="M_Boulliung23", majors="majors_Boulliung23"
        )
        # Used -213.645 instead of -213.65, 55.029 instead of 55.03 to match spreadsheet
        logC = (
            -213.645
            + (
                (
                    25696.0 * melt_comp["Na2O"]
                    + 15076.0 * melt_comp["CaO"]
                    + 9543.0 * melt_comp["MgO"]
                    + 16158.0 * melt_comp["MnO"]
                    + 4316.0 * melt_comp["Al2O3"]
                    + 68254.0
                )
                / T
            )
            + 55.029 * math.log10(T)
        )
        Csulfate = 10.0 ** (logC)

    # OLD #
    # elif model == "Nash19": # Nash et al. (2019) EPSL 507:187-198
    #    S = 1. # S6+/S2- ratio of S6+/S2- of 0.5
    #    Csulfide = C_S(PT,melt_wf,models)
    # P, T, compositional term from Kress &Carmicheal (1991)
    #    A = PT_KCterm(PT,melt_wf,models)
    # temperature dependence from Nash et al. (2019)
    #    B = (8743600/T**2) - (27703/T) + 20.273
    #    a = 0.196 # alnfO2 from Kress & Carmicheal (1991)
    #    F = 10**(((math.log10(S))-B)/8.)
    #    fO2 = math.exp(((math.log(0.5*F))-A)/a)
    #    Csulfate = (S*Csulfide)/(fO2**2)
    # elif model == "S6ST":
    #    Csulfide = C_S(PT,melt_wf,models)
    #    fO2 = f_O2(PT,melt_wf,models)
    #    S6ST_ = melt_wf["S6ST"]
    #    S = overtotal2ratio(S6ST_)
    #    Csulfate = (S*Csulfide)/(fO2**2)
    # elif model == "Hawaii":
    #    Csulfate = gp.exp(30.4) # Using Brounce et al. (2017) dataset at 1200 'C
    #    Csulfate = math.exp(slope*(1./T) -48.)
    # elif model == "Etna":
    #    Csulfate = math.exp(slope*(1./T) -50.15)
    # elif model == "Fuego":
    #    Csulfate = math.exp(slope*(1./T) -48.5)
    # elif model == "Erta Ale":
    #    Csulfate = math.exp(slope*(1./T) -45.5)
    # elif model == "FR54-S1":
    #    Csulfate = ((67.e6)*10000.)
    # elif model == "JdF": # 1100 'C ONLY
    #    Csulfate = 10.**17.

    return Csulfate


###################################
# solubility constant for H2S #####
###################################
# C_H2S = wmH2S/fH2S (ppm H2S, fH2S bar)
def C_H2S(PT, melt_wf, models=default_models):
    """
    Solubility constant for disolving H2S in the melt: C_H2S = wmH2S/fH2S (ppmw/bar)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "hydrogen sulfide" and column
        label of "option"

    Returns
    -------
    float
        Solubility constant for H2S as <class 'mpfr'>


    Model options for hydrogen sulfide
    -------------
    - 'Basalt_Hughes24' [default] Fig.S6 from Hughes et al. (2024) based on experimental
    data Moune et al. (2009) and calculations in Lesne et al. (2011)
    - 'BasalticAndesite_Hughes24' Fig.S6 from Hughes et al. (2024) based on experimental
    data Moune et al. (2009) and calculations in Lesne et al. (2011)

    """

    model = models.loc["hydrogen sulfide", "option"]
    # fitted to basalt data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015 ChemGeol
    # 418:104–116
    if model == "Basalt_Hughes24":
        K = 10.23
    # fitted to basaltic andesite data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015
    # ChemGeol 418:104–116
    elif model == "BasalticAndesite_Hughes24":
        K = 6.82
    return K


########################################
# solubility constant for hydrogen #####
########################################
# C_H2 = wmH2/fH2 (wtppm)
def C_H2(PT, melt_wf, models=default_models):
    """
    Solubility constant for disolving H2 in the melt: C_H2 = wmH2/fH2 (ppmw/bar)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "hydrogen" and column label of
        "option"

    Returns
    -------
    float
        Solubility constant for H2 as <class 'mpfr'>


    Model options for hydrogen
    -------------
    - 'Basalt_Hughes24' [default] Basalt in Table S4 from Hughes et al. (2024) based on
    experimetnal data from Hirschmann et al. (2012)
    - 'Andesite_Hughes24' Andesite in Table S4 from Hughes et al. (2024) based on
    experimental data from Hirschmann et al. (2012)

    """

    model = models.loc["hydrogen", "option"]

    # Hirchmann et al. (2012) EPSL 345-348:38-48
    R = 83.144598  # bar cm3 /mol /K
    P = PT["P"]  # pressure in bars
    T = PT["T"] + 273.15  # T in Kelvin SHOULD BE T0
    P0 = 100.0 * 0.01  # kPa to bars

    # Basalt in Table S4 from Hughes et al. (2024) based on experimetnal data from
    # Hirschmann et al. (2012)
    if model == "Basalt_Hughes24":
        # lnK0 = -11.4 # T0 = 1400 'C, P0 = 100 kPa for mole fraction H2
        lnK0 = -0.9624  # for ppm H2 (fitted in excel)
        DV = 10.6  # cm3/mol

    # Andesite in Table S4 from Hughes et al. (2024) based on experimental data from
    # Hirschmann et al. (2012)
    elif model == "Andesite_Hughes24":
        # lnK0 = -10.6 # T0 = 1400 'C, P0 = 100 kPa for mole fraction H2
        lnK0 = -0.1296  # for ppm H2 (fitted in excel)
        DV = 11.3  # cm3/mol
    lnK = lnK0 - (DV * (P - P0)) / (R * T)  # = ln(XH2/fH2) in ppm/bar
    if models.loc["high precision", "option"] == "True":
        C = gp.exp(lnK)
    else:
        C = math.exp(lnK)
    return C


######################################
# solubility constant for methane ####
######################################
# C_CH4 = wmCH4/fCH4 (ppm)
def C_CH4(PT, melt_wf, models=default_models):
    """
    Solubility constant for disolving CH4 in the melt: C_CH4 = wmCH4/fCH4 (ppmw/bar)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "methane" and column label of
        "option"

    Returns
    -------
    float
        Solubility constant for CH4 in ppmw/bar


    Model options for "methane"
    -------------
    - 'Basalt_Ardia13' [default] Eq. (7a) from Ardia et al. (2013) GCA 114:52-71
    https://doi.org/10.1016/j.gca.2013.03.028
    - Only one option available currently, included for future development.

    """

    model = models.loc["methane", "option"]

    # Eq. (7a) from Ardia et al. (2013) GCA 114:52-71
    # https://doi.org/10.1016/j.gca.2013.03.028
    if model == "Basalt_Ardia13":
        R = 83.144598  # bar cm3 /mol /K
        P = PT["P"]  # pressure in bars
        T = PT["T"] + 273.15  # T in Kelvin
        P0 = 100.0 * 0.01  # kPa to bars
        lnK0 = 4.93  # ppm CH4
        # lnK0 = -7.63 # mole fraction CH4
        DV = 26.85  # cm3/mol
        lnK = lnK0 - (DV * (P - P0)) / (R * T)
        if models.loc["high precision", "option"] == "True":
            K_ = gp.exp(lnK)  # for fCH4 in GPa
        else:
            K_ = math.exp(lnK)  # for fCH4 in GPa
        K = 0.0001 * K_  # for fCH4 in bars
    return K


#################################
# solubility constant for CO ####
#################################
# C_CO = wmCO/fCO (ppm)
def C_CO(PT, melt_wf, models=default_models):
    """
    Solubility constant for disolving CO in the melt: C_CO = wmCO/fCO (ppmw/bar)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "carbon monoxide" and column
        label of "option"

    Returns
    -------
    float
        Solubility constant for CH4 in ppmw/bar


    Model options for 'carbon monoxide'
    -------------
    - 'Basalt_Hughes24' [default] CO in Table S4 from Hughes et al. (2024) AmMin
    109:422-438 https://doi.org/10.2138/am-2023-8739. Based on data from Armstrong et
    al. (2015), Stanley et al., (2014), and Wetzel et al., (2013)
    - Only one option available currently, included for future development.

    """

    model = models.loc["carbon monoxide", "option"]

    # from fitting Armstrong et al. (2015) GCA 171:283-302; Stanley+2014, andWetzel+13
    # thermodynamically
    if model == "Basalt_Hughes24":
        R = 83.144598  # bar cm3 /mol /K
        P = PT["P"]  # pressure in bars
        T = PT["T"] + 273.15  # T in Kelvin
        P0 = 1.0  # in bars
        lnK0 = -2.11  # ppm CO
        DV = 15.20  # cm3/mol
        lnK = lnK0 - (DV * (P - P0)) / (R * T)
        if models.loc["high precision", "option"] == "True":
            K = gp.exp(lnK)  # CO(ppm)/fCO(bars)
        else:
            K = math.exp(lnK)  # CO(ppm)/fCO(bars)
    return K


#################################
# solubility constant for X #####
#################################
# C_X = wmX/fX (ppm)
def C_X(PT, melt_wf, models=default_models):
    """
    Solubility constant for disolving X in the melt: C_X = wmX/fX (ppmw/bar)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "species X solubility" and column
        label of "option"

    Returns
    -------
    float
        Solubility constant for X as <class 'mpfr'>


    Model options for species X solubility
    -------------
    - 'Ar_Basalt_HughesIP' [default] Hughes et al. (in prep) based on data from
    Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
    - Ar_Rhyolite_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano
    et al. (2010) Chemical Geology 279(3–4):145-157
    - Ne_Basalt_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et
    al. (2010) Chemical Geology 279(3–4):145-157
    - Ne_Rhyolite_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano
    et al. (2010) Chemical Geology 279(3–4):145-157
    - [user specified number]: User can type a number that will be used instead (i.e., a
    constant value)

    """

    model = models.loc["species X solubility", "option"]

    # Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical
    # Geology 279(3–4):145-157
    if model == "Ar_Basalt_HughesIP":
        K = 0.0799  # fitted assuming Ar is an ideal gas... i.e. yAr = 1.

    # Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical
    # Geology 279(3–4):145-157
    elif model == "Ar_Rhyolite_HughesIP":
        K = 0.4400  # fitted assuming Ar is an ideal gas... i.e. yAr = 1.

    # Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical
    # Geology 279(3–4):145-157
    elif model == "Ne_Basalt_HughesIP":
        K = 0.1504  # fitted assuming Ne is an ideal gas... i.e. yNe = 1.

    # Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical
    # Geology 279(3–4):145-157
    elif model == "Ne_Rhyolite_HughesIP":
        K = 0.8464  # fitted assuming Ne is an ideal gas... i.e. yNe = 1.

    # WORK IN PROGRESS
    elif model == "test":
        # K = 40. # similar to H2O
        # K = 6. # similar to S @ DFMQ+1.25
        # K = 21. # similar to S @ DFMQ+3
        # K = 155 # similar to S @ DFMQ0
        # K = 918005 # similar to S @DFMQ-3
        # K = 10.23 # similar to H2S
        # K = 0.51 # similar to CO32-
        # K = 1.37 # degassed at a similar depth to H2OT at 3wt%
        # K = 100.
        K = 35.0
    else:
        K = float(model)

    return K


########################
# solubility of Cl #####
########################
def C_Cl(PT, melt_wf):
    # WORK IN PROGRESS

    melt_comp = mg.melt_cation_proportion(melt_wf, "no", "no")
    P = PT["P"] / 10000.0  # bar to GPa
    T = PT["T"] + 273.15  # 'C to 'K
    logC = (
        1.601
        + (
            4470 * melt_comp["Ca"]
            - 3430 * melt_comp["Si"]
            + 2592 * melt_comp["FeT"]
            - 4092 * melt_comp["K"]
            - 894 * P
        )
        / T
    )
    C = float(10.0**logC)
    return C


########################################################################################
# solid/liquid volatile saturation #####################################################
########################################################################################


################################################
# sulfate content at anhydrite saturation ######
################################################
def SCAS(PT, melt_wf, models=default_models):
    """
    Sulfate content at anhydrite saturation (S6+CAS)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "SCAS" and column label of
        "option"

    Returns
    -------
    float
        S6+CAS in ppm as <class 'mpfr'>


    Model options for SCAS
    -------------
    - 'Liu23' Eq. (4) Liu et al. (2023) GCA 349:135-145 doi:10.1016/j.gca.2023.04.007
    - 'Chowdhury19_pss' Chowdhury & Dasgupta (2019) using PySulfSat by Wieser and
    Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    - 'Zajacz19_pss' Zajacz and Tsay (2019) using PySulfSat by Wieser and Gleeson (2023)
    Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    - 'Masotta15_pss' NEEDS FIXING Masotta and Kepler (2015) using PySulfSat by Wieser
    and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    """
    model = models.loc["SCAS", "option"]

    T = PT["T"] + 273.15

    comp = mg.melt_pysulfsat(melt_wf)

    # Eq. (8) using Table 5 in Chowdhury & Dasgupta (2019) Chem.Geol. 522:162-174
    # doi:10.1016/j.chemgeo.2019.05.020
    # if model == "Chowdhury19":
    #    # sulfate content (ppm) at anhydrite saturation [T in K]
    #    # mole fraction melt composition including water but all Fe as FeOT
    #    # different molecular weights used to original paper
    # melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no")
    #    tot = 100.*melt_comp["mol_tot"]
    #    a = -13.23
    #    b = -0.50
    #    dSi = 3.02
    #    dCa = 36.70
    #    dMg = 2.84
    #    dFe = 10.14
    #    dAl = 44.28
    #    dNa = 26.27
    #    dK = -25.77
    #    e = 0.09
    #    f = 0.54
    #    wm_H2OT = 100.*melt_wf['H2OT']
    #    dX = dSi*melt_comp["SiO2"] + dCa*melt_comp["CaO"] + dMg*melt_comp["MgO"] +
    #           dFe*melt_comp["FeOT"] + dAl*melt_comp["Al2O3"] + dNa*melt_comp["Na2O"] +
    #           dK*melt_comp["K2O"]
    #    if models.loc["high precision","option"] == "True":
    #        lnxm_SO4 = a + b*((10.0**4.0)/T) + dX + e*wm_H2OT -
    #        f*gp.log(melt_comp["CaO"])
    #        xm_SO4 = gp.exp(lnxm_SO4)
    #    else:
    #        lnxm_SO4 = a + b*((10.0**4.0)/T) + dX + e*wm_H2OT -
    #        f*math.log(melt_comp["CaO"])
    #        xm_SO4 = math.exp(lnxm_SO4)
    #    Xm_SO4 = xm_SO4*(xm_SO4 + tot)
    #    S6CAS = Xm_SO4*species.loc["S","M"]*10000.0

    # Eq. (8-14) Zajacz & Tsay (2019) GCA 261:288-304 doi:10.1016/j.gca.2019.07.007
    # elif model == "Zajacz19":
    #    # mole fraction melt composition including water but all Fe as FeOT
    #    # different molecular weights used to original paper
    #    melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no")
    #    tot = 100.*melt_comp["mol_tot"]
    #    if melt_comp["Na2O"] + melt_comp["K2O"] + melt_comp["CaO"] >=
    #   melt_comp["Al2O3"]:
    #        P_Rhyo = 3.11*(melt_comp["Na2O"]+melt_comp["K2O"]+melt_comp["CaO"]-
    #           melt_comp["Al2O3"])
    #    else:
    #        P_Rhyo = 1.54*(melt_comp["Al2O3"]-(melt_comp["Na2O"]+melt_comp["K2O"]+
    #       melt_comp["CaO"]))
    #    NBOT = (2.*melt_comp["Na2O"]+2.*melt_comp["K2O"]+2.*(melt_comp["CaO"]+
    #   melt_comp["MgO"]+melt_comp["FeOT"])-melt_comp["Al2O3"]*2.)/(melt_comp["SiO2"]+
    #   2.*melt_comp["Al2O3"]) # according to spreadsheet not paper
    #    P_H2O = melt_comp["H2O"]*(2.09 - 1.65*NBOT) + 0.42*NBOT + 0.23
    #    P_C = ((P_Rhyo + 251.*melt_comp["CaO"]**2. + 57.*melt_comp["MgO"]**2. +
    #   154.*melt_comp["FeOT"]**2.)/(2.*melt_comp["Al2O3"] + melt_comp["SiO2"]))/(1. +
    #   4.8*NBOT)
    #    if models.loc["high precision","option"] == "True":
    #        P_T = gp.exp(-7890./T)
    #        Ksm_SPAnh = gp.exp(1.226*gp.log(P_C*P_T*P_H2O) + 0.079)
    #    else:
    #        P_T = math.exp(-7890./T)
    #        Ksm_SPAnh = math.exp(1.226*math.log(P_C*P_T*P_H2O) + 0.079)
    #    Xsm_S = Ksm_SPAnh/melt_comp["CaO"]
    #    S6CAS = Xsm_S*tot*32.07*10000.

    # Eq. (4) Liu et al. (2023) GCA 349:135-145 doi:10.1016/j.gca.2023.04.007
    if model == "Liu23":
        melt_comp = mg.melt_mole_fraction(melt_wf, models, "no", "no")
        NBOT = (
            2.0 * melt_comp["Na2O"]
            + 2.0 * melt_comp["K2O"]
            + 2.0 * (melt_comp["CaO"] + melt_comp["MgO"] + melt_comp["FeOT"])
            - melt_comp["Al2O3"] * 2.0
        ) / (melt_comp["SiO2"] + 2.0 * melt_comp["Al2O3"])
        melt_comp = mg.melt_mole_fraction(melt_wf, models, "water", "no")
        lnSCAS = (
            12.185
            - (8541.0 / T)
            + (1.409 * NBOT)
            + 9.984 * melt_comp["CaO"]
            + melt_wf["H2OT"] * 100.0
        )
        if models.loc["high precision", "option"] == "True":
            S6CAS = gp.exp(lnSCAS)
        else:
            S6CAS = math.exp(lnSCAS)

    # Chowdhury & Dasgupta (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica
    # 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "Chowdhury19_pss":
        output = ss.calculate_CD2019_SCAS(
            df=comp,
            T_K=T,
            H2O_Liq=float(100.0 * melt_wf["H2OT"]),
            Fe3Fet_Liq=None,
            P_kbar=None,
        )
        S6CAS = float(output["SCAS6_ppm"].iloc[0])

    # Zajacz and Tsay (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica
    # 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "Zajacz19_pss":
        output = ss.calculate_ZT2019_SCAS(
            df=comp,
            T_K=T,
            H2O_Liq=float(100.0 * melt_wf["H2OT"]),
            Fe3Fet_Liq=None,
            P_kbar=None,
        )
        S6CAS = float(output["SCAS6_ppm"].iloc[0])

    # Masotta and Kepler (2015) using PySulfSat by Wieser and Gleeson (2023) Volcanica
    # 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "Masotta15_pss":
        output = ss.calculate_MK2015_SCAS(
            df=comp,
            T_K=T,
            H2O_Liq=(float(100.0 * melt_wf["H2OT"])),
            Fe3Fet_Liq=None,
            P_kbar=None,
        )
        S6CAS = float(output["SCAS6_ppm"].iloc[0])

    return S6CAS


###############################################
# sulfide content at sulfide saturation #######
###############################################
# sulfide content (ppm) at sulfide saturation
def SCSS(PT, melt_wf, models=default_models):
    """
    Sulfide content at sulfide saturation (S2-CSS)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.). Assumes sulfide is pure FeS
        unless specified in melt_wf

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "SCSS" and column label of
        "option"

    Returns
    -------
    float
        S2-CCS in ppm as <class 'mpfr'>


    Model options for SCSS
    -------------
    - 'ONeill21hyd' [default] Eq. (10.34, 10.43, 10.45, 10.46, 10.49) from O'Neill
    (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    - 'ONeill21' Eq. (10.34, 10.43, 10.45, 10.46) excluding water dilution from O'Neill
    (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    - 'ONeill21dil' Eq. (10.34, 10.43, 10.45, 10.46) including water dilution from
    O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    - 'Liu07' Eq. (9) in Liu et al. (2007) GCA 71:1783-1799
    doi:10.1016/j.gca.2007.01.004
    - 'Fortin15_pss' Fortin et al. (2015) using PySulfSat by Wieser & Gleeson (2023)
    Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    - 'Liu21_pss' Liu et al. (2021) using PySulfSat by Wieser & Gleeson (2023) Volcanica
    6(1):107-127 doi:10.30909/vol.06.01.107127
    - 'ONeill22_pss' O'Neill & Mavrogenes (2022) using PySulfSat by Wieser & Gleeson
    (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    - 'ONeill21_pss' O'Neill (2021) using PySulfSat by Wieser & Gleeson (2023) Volcanica
    6(1):107-127 doi:10.30909/vol.06.01.107127
    - 'Smythe17_pss' Smythe et al. (2017) using PySulfSat by Wieser & Gleeson (2023)
    Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    - 'Li22_pss' Li and Zhang (2022) using PySulfSat by Wieser and Gleeson (2023)
    Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    - 'Blanchard21_eq11_pss' Eq. (11) from Blanchard et al. (2021) using PySulfSat by
    Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    - 'Blanchard21_eq12_pss' Eq. (12) from Blanchard et al. (2021) using PySulfSat by
    Wieser and Gleeson (2023) Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127

    """

    model = models.loc["SCSS", "option"]

    P_bar = PT["P"]
    T = PT["T"] + 273.15
    Fe3FeT = melt_wf["Fe3FeT"]
    comp = mg.melt_pysulfsat(melt_wf)

    if "sulf_XFe" in melt_wf:
        sulf_XFe = melt_wf["sulf_XFe"]
    else:
        sulf_XFe = 1.0
    if "sulf_XCu" in melt_wf:
        sulf_XCu = melt_wf["sulf_XCu"]
    else:
        sulf_XCu = 0.0
    if "sulf_XNi" in melt_wf:
        sulf_XNi = melt_wf["sulf_XNi"]
    else:
        sulf_XNi = None

    # O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
    if model in [
        "ONeill21",
        "ONeill21dil",
        "ONeill21hyd",
    ]:
        R = 8.31441
        P = (1.0e-4) * P_bar  # pressure in GPa
        if models.loc["high precision", "option"] == "True":
            D = 137778.0 - 91.666 * T + 8.474 * T * gp.log(T)  # J/mol Eq. (10.45)
        else:
            D = 137778.0 - 91.666 * T + 8.474 * T * math.log(T)  # J/mol Eq. (10.45)
        sulfide_comp = sulf_XFe
        if model == "ONeill21":  # Eq. (10.34, 10.43, 10.45, 10.46)
            # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) no volatiles
            melt_comp = mg.melt_cation_proportion(
                melt_wf, "no", "no", molmass="M_ONeill21", majors="majors_ONeill21"
            )
        elif model == "ONeill21dil":  # Eq. (10.34, 10.43, 10.45, 10.46)
            # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) and water
            melt_comp = mg.melt_cation_proportion(
                melt_wf, "water", "no", molmass="M_ONeill21", majors="majors_ONeill21"
            )
        elif model == "ONeill21hyd":  # Eq. (10.34, 10.43, 10.45, 10.46, 10.49)
            # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) and water
            melt_comp = mg.melt_cation_proportion(
                melt_wf, "water", "no", molmass="M_ONeill21", majors="majors_ONeill21"
            )
        Fe2 = melt_comp["FeT"] * (1.0 - Fe3FeT)
        if models.loc["high precision", "option"] == "True":
            lnaFeS = gp.log((1.0 - Fe2) * sulfide_comp)
        else:
            lnaFeS = math.log((1.0 - Fe2) * sulfide_comp)
        # lnyFe2 from Eq. (10.46)
        lnyFe2 = (
            ((1.0 - Fe2) ** 2.0)
            * (
                28870.0
                - 14710.0 * melt_comp["Mg"]
                + 1960.0 * melt_comp["Ca"]
                + 43300.0 * melt_comp["Na"]
                + 95380.0 * melt_comp["K"]
                - 76880.0 * melt_comp["Ti"]
            )
            + (1.0 - Fe2)
            * (-62190.0 * melt_comp["Si"] + 31520.0 * melt_comp["Si"] ** 2.0)
        ) / (R * T)
        # lnS from Eq. (10.43)
        if models.loc["high precision", "option"] == "True":
            lnS = (
                D / (R * T)
                + gp.log(C_S(PT, melt_wf, models))
                - gp.log(Fe2)
                - lnyFe2
                + lnaFeS
                + (-291.0 * P + 351.0 * gp.erf(P)) / T
            )
            SCSS = gp.exp(lnS)
        else:
            lnS = (
                D / (R * T)
                + math.log(C_S(PT, melt_wf, models))
                - math.log(Fe2)
                - lnyFe2
                + lnaFeS
                + (-291.0 * P + 351.0 * math.erf(P)) / T
            )
            SCSS = math.exp(lnS)

    # Eq. (9) in Liu et al. (2007) GCA 71:1783-1799 doi:10.1016/j.gca.2007.01.004
    elif model == "Liu07":
        # Mole fractions in the melt on cationic lattice (Fe2 and Fe3) and water
        melt_comp = mg.melt_cation_proportion(melt_wf, "water", "yes")
        MFM = (
            melt_comp["Na"]
            + melt_comp["K"]
            + 2.0 * (melt_comp["Ca"] + melt_comp["Mg"] + melt_comp["Fe2"])
        ) / (melt_comp["Si"] * (melt_comp["Al"] + melt_comp["Fe3"]))
        if models.loc["high precision", "option"] == "True":
            lnS = (
                11.35251
                - (4454.6 / T)
                - 0.03190 * (PT["P"] / T)
                + 0.71006 * gp.log(MFM)
                - 1.98063 * (MFM * melt_comp["H"])
                + 0.21867 * gp.log(melt_comp["H"])
                + 0.36192 * gp.log(melt_comp["Fe2"])
            )
            SCSS = gp.exp(lnS)
        else:
            lnS = (
                11.35251
                - (4454.6 / T)
                - 0.03190 * (PT["P"] / T)
                + 0.71006 * math.log(MFM)
                - 1.98063 * (MFM * melt_comp["H"])
                + 0.21867 * math.log(melt_comp["H"])
                + 0.36192 * math.log(melt_comp["Fe2"])
            )
            SCSS = math.exp(lnS)

    # Eq. (7) Fortin et al. (2015) GCA 160:100-116 doi:10.1016/j.gca.2015.03.022
    # elif model == "Fortin15":
    #    # Mole fractions in the melt on cationic lattice (all Fe as FeOT) and water -
    # molecular masses used are different to spreadsheet attached to paper
    #    melt_comp = mg.melt_mole_fraction(melt_wf,models,"water","no")
    #    lnS = 34.7837 - (5772.322/T) - 346.5377*((0.0001*PT["P"])/T) -
    # 20.3934*melt_comp["H2O"] - 25.4986*melt_comp["SiO2"] - 18.3435*melt_comp["TiO2"] -
    # 27.3807*melt_comp["Al2O3"] - 17.2752*melt_comp["FeOT"] - 22.3975*melt_comp["MgO"]
    # - 20.3778*melt_comp["CaO"] - 18.9539*melt_comp["Na2O"] - 32.1944*melt_comp["K2O"]
    #    if models.loc["high precision","option"] == "True":
    #        SCSS = gp.exp(lnS)
    #    else:
    #        SCSS = math.exp(lnS)

    # Eq. (2) Liu et al. (2021) Chem.Geol. 559:119913 doi:10.1016.j.chemgeo.2020.119913
    # elif model == "Liu21":
    #    XFeS = sulf_XFe
    #    H2O = melt_wf["H2OT"]*100.
    #    if models.loc["high precision","option"] == "True":
    #        SCSS = (XFeS*gp.exp(13.88 - (9744./T) - (328.*(0.0001*PT["P"])/T))) +
    #   104.*H2O
    #    else:
    #        SCSS = (XFeS*math.exp(13.88 - (9744./T) - (328.*(0.0001*PT["P"])/T))) +
    #   104.*H2O

    # Li and Zhang (2022) using PySulfSat by Wieser and Gleeson (2023) Volcanica
    # 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "Li22_pss":
        output = ss.calculate_LiZhang2022_SCSS(
            df=comp,
            T_K=T,
            P_kbar=(P_bar / 1000.0),
            H2O_Liq=float(100.0 * melt_wf["H2OT"]),
            Fe_FeNiCu_Sulf=sulf_XFe,
            Cu_FeNiCu_Sulf=sulf_XCu,
            Ni_FeNiCu_Sulf=sulf_XNi,
            Fe3Fet_Liq=Fe3FeT,
        )
        SCSS = float(output["SCSS_Tot"].iloc[0])

    # Eq. (11) from Blanchard et al. (2021) using PySulfSat by Wieser and Gleeson (2023)
    # Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "Blanchard21_eq11_pss":
        output = ss.calculate_B2021_SCSS(
            df=comp,
            T_K=T,
            P_kbar=(P_bar / 1000.0),
            H2O_Liq=float(100.0 * melt_wf["H2OT"]),
            Fe_FeNiCu_Sulf=sulf_XFe,
            Cu_FeNiCu_Sulf=sulf_XCu,
            Ni_FeNiCu_Sulf=sulf_XNi,
            Fe3Fet_Liq=Fe3FeT,
        )
        SCSS = float(output["SCSS2_ppm_eq11"].iloc[0])

    # Eq. (12) from Blanchard et al. (2021) using PySulfSat by Wieser and Gleeson (2023)
    # Volcanica 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "Blanchard21_eq12_pss":
        output = ss.calculate_B2021_SCSS(
            df=comp,
            T_K=T,
            P_kbar=(P_bar / 1000.0),
            H2O_Liq=float(100.0 * melt_wf["H2OT"]),
            Fe_FeNiCu_Sulf=sulf_XFe,
            Cu_FeNiCu_Sulf=sulf_XCu,
            Ni_FeNiCu_Sulf=sulf_XNi,
            Fe3Fet_Liq=Fe3FeT,
        )
        SCSS = float(output["SCSS2_ppm_eq12"].iloc[0])

    # Fortin et al. (2015) using PySulfSat by Wieser & Gleeson (2023) Volcanica
    # 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "Fortin15_pss":
        output = ss.calculate_F2015_SCSS(
            df=comp,
            T_K=T,
            P_kbar=(P_bar / 1000.0),
            H2O_Liq=float(100.0 * melt_wf["H2OT"]),
        )
        SCSS = float(output["SCSS2_ppm"].iloc[0])

    # Liu et al. (2021) using PySulfSat by Wieser & Gleeson (2023) Volcanica
    # 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "Liu21_pss":
        output = ss.calculate_Liu2021_SCSS(
            df=comp,
            T_K=T,
            P_kbar=(P_bar / 1000.0),
            H2O_Liq=float(100.0 * melt_wf["H2OT"]),
            Fe_FeNiCu_Sulf=sulf_XFe,
            Cu_FeNiCu_Sulf=sulf_XCu,
            Ni_FeNiCu_Sulf=sulf_XNi,
            Fe3Fet_Liq=Fe3FeT,
        )
        SCSS = float(output["SCSS2_ppm"].iloc[0])

    # O'Neill & Mavrogenes (2022) using PySulfSat by Wieser & Gleeson (2023) Volcanica
    # 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "ONeill22_pss":
        output = ss.calculate_OM2022_SCSS(
            df=comp,
            T_K=T,
            P_kbar=(P_bar / 1000.0),
            Fe3Fet_Liq=Fe3FeT,
            Fe_FeNiCu_Sulf=sulf_XFe,
            Cu_FeNiCu_Sulf=sulf_XCu,
            Ni_FeNiCu_Sulf=sulf_XNi,
        )
        SCSS = float(output["SCSS2_ppm"].iloc[0])

    # O'Neill (2021) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127
    # doi:10.30909/vol.06.01.107127
    elif model == "ONeill21_pss":
        output = ss.calculate_O2021_SCSS(
            df=comp,
            T_K=T,
            P_kbar=(P_bar / 1000.0),
            Fe3Fet_Liq=Fe3FeT,
            Fe_FeNiCu_Sulf=sulf_XFe,
            Cu_FeNiCu_Sulf=sulf_XCu,
            Ni_FeNiCu_Sulf=sulf_XNi,
        )
        SCSS = float(output["SCSS2_ppm"].iloc[0])

    # Smythe et al. (2017) using PySulfSat by Wieser & Gleeson (2023) Volcanica
    # 6(1):107-127 doi:10.30909/vol.06.01.107127
    elif model == "Smythe17_pss":
        output = ss.calculate_S2017_SCSS(
            df=comp,
            T_K=T,
            P_kbar=(P_bar / 1000.0),
            Fe3Fet_Liq=Fe3FeT,
            Fe_FeNiCu_Sulf=sulf_XFe,
            Cu_FeNiCu_Sulf=sulf_XCu,
            Ni_FeNiCu_Sulf=sulf_XNi,
            H2O_Liq=float(100.0 * melt_wf["H2OT"]),
        )
        SCSS = float(output["SCSS2_ppm_ideal_Smythe2017"].iloc[0])
    return SCSS


########################################################################################
# EQUILIBRIUM CONSTANTS FOR HOMOGENEOUS VAPOR EQUILIBRIA ###############################
########################################################################################


# H2 + 0.5O2 = H2O
# K = fH2O/(fH2*(fO2)^0.5)
def KHOg(PT, models=default_models):
    """
    Equilibrium constant for H2 + 0.5O2 = H2O, K = fH2O/(fH2*(fO2)^0.5)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KHOg" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for KHOg
    -------------
    - 'Ohmoto97' [default] Reaction (d) in Table 1 of Ohmoto & Kerrick (1997)
    - Only one option available currently, included for future development.

    """

    model = models.loc["KHOg", "option"]

    T_K = PT["T"] + 273.15
    if model == "Ohmoto97":  # Reaction (d) in Table 1 of Ohmoto & Kerrick (1997)
        if models.loc["high precision", "option"] == "True":
            K = 10.0 ** ((12510.0 / T_K) - 0.979 * (gp.log10(T_K)) + 0.483)
        else:
            K = 10.0 ** ((12510.0 / T_K) - 0.979 * (math.log10(T_K)) + 0.483)
    return K


# H2O + 0.5S2 = H2S + 0.5O2
# K = (fH2S*(fO2)^0.5)/((fS2^0.5)*fH2O)
def KHOSg(PT, models=default_models):
    """
    Equilibrium constant for H2O + 0.5S2 = H2S + 0.5O2,
    K = (fH2S*(fO2)^0.5)/((fS2^0.5)*fH2O)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KHOSg" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for KHOSg
    -------------
    - 'Ohmoto97' [default] Reaction (h) in Table 1 of Ohmoto & Kerrick (1997)
    - 'noH2S' Stops H2S forming in the vapor, K = 0.

    """

    model = models.loc["KHOSg", "option"]

    T_K = PT["T"] + 273.15
    if model == "Ohmoto97":  # Reaction (h) in Table 1 of Ohmoto & Kerrick (1997)
        if models.loc["high precision", "option"] == "True":
            K = 10.0 ** ((-8117.0 / T_K) + 0.188 * gp.log10(T_K) - 0.352)
        else:
            K = 10.0 ** ((-8117.0 / T_K) + 0.188 * math.log10(T_K) - 0.352)
    elif model == "noH2S":  # H2S doesn't form in the gas...
        K = 0.0
    return K


# 0.5S2 + O2 = SO2
# K = fSO2/((fS2^0.5)*fO2)
def KOSg(PT, models=default_models):
    """
    Equilibrium constant for 0.5S2 + O2 = SO2, K = fSO2/((fS2^0.5)*fO2)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KOSg" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for KOSg
    -------------
    - 'Ohmoto97' [default] Reaction (f) in Table 1 of Ohmoto & Kerrick (1997)
    - 'noSO2' Stops SO2 forming in the vapor, K = 0.

    """

    model = models.loc["KOSg", "option"]

    T_K = PT["T"] + 273.15
    if model == "Ohmoto97":  # Reaction (f) in Table 1 of Ohmoto & Kerrick (1997)
        K = 10.0 ** ((18929.0 / T_K) - 3.783)
    if model == "noSO2":
        K = 0.0
    return K


# 0.5S2 + 1.5O2 = SO3
# K = fSO3/((fS2^0.5)*(fO2^1.5)
def KOSg2(PT, models=default_models):
    """
    Equilibrium constant for 0.5S2 + 1.5O2 = SO3, K = fSO3/((fS2^0.5)*(fO2^1.5)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KOsg2" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for KOSg2
    -------------
    - 'ONeill22' [default] Eq (6b) in O’Neill and Mavrogenes (2022)
    - Only one option available currently, included for future development.

    """

    model = models.loc["KOSg2", "option"]

    T_K = PT["T"] + 273.15
    if model == "ONeill22":  # Eq (6b) in O’Neill and Mavrogenes (2022)
        if models.loc["high precision", "option"] == "True":
            lnK = (55921.0 / T_K) - 25.07 + 0.6465 * gp.log(T_K)
            K = gp.exp(lnK)
        else:
            lnK = (55921.0 / T_K) - 25.07 + 0.6465 * math.log(T_K)
            K = math.exp(lnK)
    return K


# CO + 0.5O = CO2
# K = fCO2/(fCO*(fO2^0.5))
def KCOg(PT, models=default_models):
    """
    Equilibrium constant for CO + 0.5O = CO2, K = fCO2/(fCO*(fO2^0.5))


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KCOg" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for KCOg
    -------------
    - 'Ohmoto97' [default] Reaction (c) in Table 1 of Ohmoto & Kerrick (1997)
    - Only one option available currently, included for future development.

    """

    model = models.loc["KCOg", "option"]

    T_K = PT["T"] + 273.15
    if model == "Ohmoto97":  # Reaction (c) in Table 1 of Ohmoto & Kerrick (1997)
        K = 10.0 ** ((14751.0 / T_K) - 4.535)
    return K


# CH4 + 2O2 = CO2 + 2H2O
# K = (fCO2*(fH2O^2))/(fCH4*(fO2^2))
def KCOHg(PT, models=default_models):
    """
    Equilibrium constant for CH4 + 2O2 = CO2 + 2H2O, K = (fCO2*(fH2O^2))/(fCH4*(fO2^2))


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KCOHg" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for KCOHg
    -------------
    - 'Ohmoto97' [default] Reaction (e) in Table 1 of Ohmoto & Kerrick (1997)
    - 'noCH4' Almost stops CH4 forming in the vapor, K = very large.

    """

    model = models.loc["KCOHg", "option"]

    T_K = PT["T"] + 273.15
    if model == "Ohmoto97":  # Reaction (e) in Table 1 of Ohmoto & Kerrick (1997)
        if models.loc["high precision", "option"] == "True":
            K = 10.0 ** ((41997.0 / T_K) + 0.719 * gp.log10(T_K) - 2.404)
        else:
            K = 10.0 ** ((41997.0 / T_K) + 0.719 * math.log10(T_K) - 2.404)
    if model == "noCH4":
        K = 1.0e40  # really big
    return K


def KOCSg(PT, models=default_models):  # OCS - depends on system
    """
    Equilibrium constant for OCS (either K = (fCO2*fH2S)/(fOCS*fH2O) or
    (fCO^3*fSO2)/(fCO2^2*fOCS))


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KOCSg" and "carbonylsulfide" and
        column label of "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for KOCSg
    -------------
    - 'Moussallam19' [default] Eq. (8) in Moussallam et al. (2019) for KOCSg AND 'COS'
    for carbonlysulfide
    - 'noOCS' Almost stops OCS forming in the vapor, K = very large.

    """

    reaction = models.loc["carbonylsulfide", "option"]
    model = models.loc["KOCSg", "option"]

    T = PT["T"] + 273.15
    if reaction == "COS":
        # 2CO2 + OCS = 3CO + SO2 -
        # K = (fCO^3*fSO2)/(fCO2^2*fOCS)
        if model == "Moussallam19":  # Eq. (8) in Moussallam et al. (2019)
            K = 10.0 ** (9.24403 - (15386.45 / T))  # P and f in bars, T in K
        if model == "noOCS":
            K = 1.0e30  # really big
        return K

    # WORK IN PROGRESS
    if reaction == "COHS":
        # OCS + H2O = CO2 + H2S
        # K = (fCO2*fH2S)/(fOCS*fH2O)
        if models == "EVo":
            if models.loc["high precision", "option"] == "True":
                K = gp.exp(
                    0.482
                    + (16.166e-2 / T)
                    + 0.081e-3 * T
                    - (5.715e-3 / T**2)
                    - 2.224e-1 * gp.log(T)
                )
            else:
                K = math.exp(
                    0.482
                    + (16.166e-2 / T)
                    + 0.081e-3 * T
                    - (5.715e-3 / T**2)
                    - 2.224e-1 * math.log(T)
                )
            return K


# Cgraphite + O2 = CO2
def KCOs(PT, models=default_models):
    """
    Equilibrium constant for Cgraphite + O2 = CO2


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "KCOs" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for KCOs
    -------------
    - 'Holloway92' [default] Eq (3) KI in Holloway et al. (1992) Eur J. Mineral.
    4:105-114
    Only one option available currently, included for future development.

    """

    model = models.loc["KCOs", "option"]

    T_K = PT["T"] + 273.15
    P = PT["P"]
    if (
        model == "Holloway92"
    ):  # Eq (3) KI in Holloway et al. (1992) Eur J. Mineral. 4:105-114
        a = 40.07639
        b = -2.5392e-2
        c = 5.27096e-6
        d = 0.0267
        log10K = a + (b * T_K) + (c * T_K**2) + d * ((P - 1.0) / T_K)
        K = 10.0**log10K
    return K


########################################################################################
# SPECIATION CONSTANTS #################################################################
########################################################################################


# H2Omol + O = 2OH
# K = xOH*2/(xH2Omol*xO)
def KHOm(PT, melt_wf, models=default_models):
    """
    Speciation constant for H2Omol + O = 2OH- assuming ideal mixing


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "Hspeccomp" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for Hspeccomp
    -------------
    - 'MORB_HughesIP' [default] Eq. SX in Hughes et al. (in prep) based on data from
    Dixon et al. (1995)
    - StromboliAlkaliBasalt_Lesne10": Eq. (15) Lesne et al. (2010) CMP 162:133-151
    - VesuviusAlkaliBasalt_Lesne10: Eq. (16) Lesne et al. (2010) CMP 162:133-151
    - EtnaAlkaliBasalt_Lesne10: Eq. (17) Lesne et al. (2010) CMP 162:133-151
    - Andesite_Botcharnikov06: Eq (7) from Botcharnikov et al. (2006) Chem. Geol.
    229(1-3)125-143
    - Albite_Silver89: Fig. 8 from Silver & Stolper (1989) J.Pet 30(3)667-709
    - Rhyolite_Zhang97: Eq. (9) from Zhang et al. (1997) GCA 61(15):3089-3100

    """

    Hspeccomp = models.loc["Hspeccomp", "option"]

    T_K = PT["T"] + 273.15

    if (
        Hspeccomp == "Rhyolite_Zhang97"
    ):  # Eq. (9) from Zhang et al. (1997) GCA 61(15):3089-3100
        a = -3110.0
        b = 1.876
        if models.loc["high precision", "option"] == "True":
            K = gp.exp((a / T_K) + b)
        else:
            K = math.exp((a / T_K) + b)
    elif (
        Hspeccomp == "StromboliAlkaliBasalt_Lesne10"
    ):  # Eq. (15) Lesne et al. (2010) CMP 162:133-151
        a = -8710.0
        b = 8.5244
        if models.loc["high precision", "option"] == "True":
            K = gp.exp((a / T_K) + b)
        else:
            K = math.exp((a / T_K) + b)
    elif (
        Hspeccomp == "VesuviusAlkaliBasalt_Lesne10"
    ):  # Eq. (16) Lesne et al. (2010) CMP 162:133-151
        a = -8033.0
        b = 7.4222
        if models.loc["high precision", "option"] == "True":
            K = gp.exp((a / T_K) + b)
        else:
            K = math.exp((a / T_K) + b)
    elif (
        Hspeccomp == "EtnaAlkaliBasalt_Lesne10"
    ):  # Eq. (17) Lesne et al. (2010) CMP 162:133-151
        a = -8300.0
        b = 7.4859
        if models.loc["high precision", "option"] == "True":
            K = gp.exp((a / T_K) + b)
        else:
            K = math.exp((a / T_K) + b)
    elif (
        Hspeccomp == "Andesite_Botcharnikov06"
    ):  # Eq (7) from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143
        a = -3650.0
        b = 2.99
        if models.loc["high precision", "option"] == "True":
            K = gp.exp((a / T_K) + b)
        else:
            K = math.exp((a / T_K) + b)

    # fit to Dixon et al. (1995) data digitised from Lesne et al. (2010) CMP 162:133-151
    # in Hughes et al. (in prep)
    elif Hspeccomp == "MORB_HughesIP":
        a = -2204.99
        b = 1.2600
        if models.loc["high precision", "option"] == "True":
            K = gp.exp((a / T_K) + b)
        else:
            K = math.exp((a / T_K) + b)

    # Fig. 8 from Silver & Stolper (1989) J.Pet 30(3)667-709
    elif Hspeccomp == "Albite_Silver89":
        K = 0.17
    else:
        K = 0.17

    # Work in progress
    if (
        Hspeccomp == "AlkaliBasalt"
    ):  # average of eqn-15-17 from Lesne et al. (2010) CMP 162:133-151
        a = -8348.0  # VES-9 = -8033.0, ETN-1 = -8300.0, and PST-9 = -8710.0
        b = 7.8108  # VES-9 = 7.4222, ETN-1 = 7.4859, and PEST-9 = 8.5244
        if models.loc["high precision", "option"] == "True":
            K = gp.exp((a / T_K) + b)
        else:
            K = math.exp((a / T_K) + b)

    return K


def KregH2O(PT, melt_wf, models=default_models):
    """
    Speciation constant for H2Omol + O = 2OH- assuming regular mixing


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "Hspeccomp" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for Hspeccomp
    -------------
    - 'MORB_HughesIP' [default] WHICH IS NOT USABLE WITH THIS FUNCTION
    - MORB_Dixon95: Table 5 of Dixon et al. (1995)
    - AlkaliBasalt_Lesne10": Eq. (24-27) Lesne et al. (2010) CMP 162:133-151
    - StromboliAlkaliBasalt_Lesne10": PST-9 in Table 5 from Lesne et al. (2010)
    162:133-151
    - VesuviusAlkaliBasalt_Lesne10: VES-9 in Table 5 from Lesne et al. (2010)
    162:133-151
    - EtnaAlkaliBasalt_Lesne10: ETN-1 in Table 5 from Lesne et al. (2010) 162:133-151
    - Albite_Silver89: Silver & Stolper (1989) J.Pet 30(3)667-709

    """

    Hspeccomp = models.loc["Hspeccomp", "option"]

    if Hspeccomp in [
        "MORB_Dixon95",
        "Albite_Silver89",
    ]:  # Table 5 of Dixon et al. (1995); Silver & Stolper (1989) J.Pet 30(3)667-709
        A = 0.403
        B = 15.333
        C = 10.894
    elif (
        Hspeccomp == "AlkaliBasalt_Lesne10"
    ):  # Eq (24-27) from Lesne et al. (2010) CMP 162:133-151
        T_K = PT["T"] + 273.15
        R = 83.15  # cm3bar/molK
        lnK = (-2704.4 / T_K) + 0.641
        A = lnK + 49016.0 / (R * T_K)
        B = -2153326.51 / (R * T_K)
        C = 1.965495217 / (R * T_K)
    elif (
        Hspeccomp == "VesuviusAlkaliBasalt_Lesne10"
    ):  # VES-9 in Table 5 from Lesne et al. (2010) 162:133-151
        A = 3.139
        B = -29.555
        C = 20.535
    elif (
        Hspeccomp == "EtnaAlkaliBasalt_Lesne10"
    ):  # ETN-1 in Table 5 from Lesne et al. (2010) 162:133-151
        A = 4.128
        B = -45.905
        C = 21.311
    elif (
        Hspeccomp == "StromboliAlkaliBasalt_Lesne10"
    ):  # PST-9 in Table 5 from Lesne et al. (2010) 162:133-151
        A = 2.6
        B = -22.476
        C = 22.295

    # Work in progress
    # No T-dependence, hence its the speciation frozen in the glass. Eqn 7-10 from Lesne
    # et al. (2010) CMP 162:133-151 (eqn 7 is wrong)
    elif Hspeccomp == "alkali basalt XT":
        # wt% normalised including H2O, all Fe as FeOT
        melt_comp = mg.melt_normalise_wf(melt_wf, "volatiles", "Fe speciation")
        Na = melt_comp["Na2O"] * 100.0
        K = melt_comp["K2O"] * 100.0
        A = 0.5761 * (Na + K) - 0.2884  # eqn-8
        B = -8.9589 * (Na + K) + 24.65  # eqn-9
        C = 1.7013 * (Na + K) + 9.6481  # eqn-1
    return A, B, C


# CO2 + O = CO3
def KCOm(PT, melt_wf, models=default_models):  # K =
    """
    Speciation constant for CO2 + O = CO3


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)
        Not normally used unless model option requires melt composition.

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "Cspeccomp" and column label of
        "option"

    Returns
    -------
    float
        Equilibrium constant as <class 'mpfr'>


    Model options for Cspeccomp
    -------------
    - 'Basalt' [default] Assume all oxidised carbon in the melt is present as carbonate
    ions.
    - Andesite_Botcharnikov06: Eq. (8) from Botcharnikov et al. (2006) Chem. Geol.
    229(1-3)125-143
    - Dacite_Botcharnikov06: Eq. in the text from Botcharnikov et al. (2006), based on
    data from Behrens et al. (2004)
    - Rhyolite: Assume all oxidised carbon in the melt is present as molecular CO2.

    """

    Cspeccomp = models.loc["Cspeccomp", "option"]

    T_K = PT["T"] + 273.15

    # Eq. (8) from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143
    if Cspeccomp == "Andesite_Botcharnikov06":
        a = 8665.0
        b = -5.11
        if models.loc["high precision", "option"] == "True":
            value = gp.exp((a / T_K) + b)
        else:
            value = math.exp((a / T_K) + b)

    # Eq. in the text from Botcharnikov et al. (2006), based on data from Behrens et al.
    # (2004)
    elif Cspeccomp == "Dacite_Botcharnikov06":
        a = 9787.0
        b = -7.69
        if models.loc["high precision", "option"] == "True":
            value = gp.exp((a / T_K) + b)
        else:
            value = math.exp((a / T_K) + b)
    elif Cspeccomp == "Basalt":  # all oxidised carbon is CO32-
        value = "infinite"
    elif Cspeccomp == "Rhyolite":  # all oxidised carbon is CO2,mol
        value = 0.0
    return value


########################################################################################
# FUGACITY COEFFICIENTS ################################################################
########################################################################################

# all fugacity coefficients are assumed to equal 1 below 1 bar.

##########################################
# CORK using Holland & Powell (1991) #####
##########################################


def y_CORK(species, PT, models):
    """
    Fugacity coefficient using eq. (4,A1-3) from Holland & Powell (1991) CMP 109:265-273
    10.1007/BF00306484


    Parameters
    ----------
    species: str
        species of interest (e.g., 'H2O', 'CO2')

    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Dataframe detailing model options

    Returns
    -------
    float
        Fugacity coefficient
    """
    P = PT["P"]
    T_K = PT["T"] + 273.15
    R = 8.3144598e-3  # in kJ/mol/K
    P_kb = P / 1000.0

    # Appendix: Calculation of CORK volumes
    V = vol_CORK(species, PT, models)

    # Appendix: Calculation of CORK fugacities
    a, b, c, d, p0 = parameters_Holland91(species, PT, models)
    if P_kb > p0:
        # Eq. (A.3)
        ln_y_virial = (1 / (R * T_K)) * (
            (2.0 / 3.0) * c * pow((P_kb - p0), 1.5) + (d / 2.0) * pow((P_kb - p0), 2.0)
        )
    else:
        ln_y_virial = 0.0

    z = (P_kb * V) / (R * T_K)
    A = a / (b * R * pow(T_K, 1.5))
    B = (b * P_kb) / (R * T_K)

    if z < B:
        value = 1.0
    elif models.loc["high precision", "option"] == "True":
        # Eq. (A.2)
        ln_y = z - 1.0 - gp.log(z - B) - A * gp.log(1.0 + (B / z)) + ln_y_virial
        value = gp.exp(ln_y)
    else:
        # Eq. (A.2)
        ln_y = z - 1.0 - math.log(z - B) - A * math.log(1.0 + (B / z)) + ln_y_virial
        value = math.exp(ln_y)

    return value


def y_sCORK(species, PT, models):
    """Fugacity coefficient using eq. (8) from Holland & Powell (1991) CMP 109:265-273
    10.1007/BF00306484

    Args:
        species (str): species of interest (e.g., 'H2O', 'CO2')
        PT (dict): Dictionary of pressure-temperature conditions with pressure (bars) as
        "P" and temperature ('C) as "T"

    Returns:
        float: fugacity coefficient
    """
    R = 8.3144598e-3
    P = PT["P"]  # P in bar
    T_K = PT["T"] + 273.15  # T in K
    P_kb = P / 1000.0  # P in kb

    a, b, c, d, p0 = parameters_Holland91(species, PT, models)  # noqa

    # Eq. (8) rearranged for lnf
    lnf = (
        (R * T_K * math.log(1000 * P_kb))
        + (b * P_kb)
        + (a / (b * T_K**0.5))
        * ((math.log(R * T_K + b * P_kb)) - (math.log(R * T_K + 2 * b * P_kb)))
        + ((2.0 / 3.0) * c * P_kb * P_kb**0.5)
        + ((d / 2.0) * P_kb**2.0)
    ) / (R * T_K)
    y = (math.exp(lnf)) / P
    return y


def parameters_Holland91(species, PT, models):
    """Parameters for (simplified) CORK equations in Holland & Powell (1991) CMP =
    109:265-273 10.1007/BF00306484

    Args:
        species (str): species of interest (e.g., 'H2O', 'CO2')
        PT (dict): Dictionary of pressure-temperature conditions with pressure (bars) as
        "P" and temperature ('C) as "T"
        models (pandas.DataFrame): Dataframe detailing model options

    Returns:
        tuple(float,float,float,float,float): a, b, c, d, p0
    """

    # Parameters for gas species of interest using corresponding states using eq. (9)
    # and Table 2
    def corresponding_states_Holland91(PT, Tc, Pc):
        T = PT["T"] + 273.1  # T in K
        # Table 2
        a0 = 5.45963e-5
        a1 = -8.63920e-6
        b0 = 9.18301e-4
        c0 = -3.30558e-5
        c1 = 2.30524e-6
        d0 = 6.93054e-7
        d1 = -8.38293e-8
        # Eq. (9)
        a = (
            a0 * (Tc ** (5 / 2) / Pc) + a1 * (Tc * (3 / 2) / Pc) * T
        )  # Kj2 kbar^-1 K^(1/2) mol^(-2)
        b = b0 * (Tc / Pc)  # kJ kbar^-1 mol^-1
        c = c0 * (Tc / Pc ** (3 / 2)) + (c1 / Pc ** (3 / 2)) * T
        d = d0 * (Tc / Pc**2) + (d1 / Pc**2) * T
        return a, b, c, d

    def constants_CO2_Holland91(PT, models):
        model = models.loc["y_CO2", "option"]
        T_K = PT["T"] + 273.15

        if model in ["Holland91_eq4,A1-3_tab1", "Holland91_eq8_tab1"]:
            # Table 1
            a = 741.2 + -0.10891 * (T_K) + -3.4203e-4 * pow(T_K, 2.0)
            b = 3.057
            c = -2.26924e-1 + 7.73793e-5 * T_K  # Eq. (4)
            d = 1.33790e-2 + -1.01740e-5 * T_K  # Eq. (4)
            p0 = 5.0  # kbar
        elif model == "Holland91_eq8,9_tab2":
            Tc = 304.2  # Critical temperature in K
            Pc = 0.0738  # Critical pressure in kbar
            a, b, c, d = corresponding_states_Holland91(PT, Tc, Pc)  # Eq. (9), Table 2
            p0 = ""

        return a, b, c, d, p0

    def constants_H2O_Holland91(PT):
        T_K = PT["T"] + 273.15

        # Table 1
        p0 = 2.00  # in kb
        # Eq. (6) T >673 K
        a = (
            1113.4
            + -0.22291 * (T_K - 673.0)
            + -3.8022e-4 * pow((T_K - 673.0), 2.0)
            + 1.7791e-7 * pow((T_K - 673.0), 3.0)
        )
        b = 1.465
        c = -3.025650e-2 + -5.343144e-6 * T_K
        d = -3.2297554e-3 + 2.2215221e-6 * T_K

        return a, b, c, d, p0

    if species == "H2O":
        a, b, c, d, p0 = constants_H2O_Holland91(PT)
    elif species == "CO2":
        a, b, c, d, p0 = constants_CO2_Holland91(PT, models)
    elif species in ["CH4", "CO", "H2"]:
        if species == "CH4":
            Tc = 190.6  # Critical temperature in K
            Pc = 0.0460  # Critical pressure in kbar
        elif species == "H2":
            Tc = 41.26  # Critical temperature in K
            Pc = 0.0211  # Critical pressure in kbar
        elif species == "CO":
            Tc = 132.9  # Critical temperature in K
            Pc = 0.0350  # Critical pressure in kbar
        a, b, c, d = corresponding_states_Holland91(PT, Tc, Pc)
        p0 = ""

    return a, b, c, d, p0


# Flowers (1979) modified from code from MIMiC (Rasmussen et al., 2021:
# https://github.com/DJRgeoscience/MIMiC), originally from VolatileCalc (Newman &
# Lowenstern, 2001)
def MRK(PT, X_1):  # Redlich-Kwong routine to estimate endmember H2O and CO2 fugacities
    # CO2, X_1 = 0.
    # H2O, X_1 = 1.
    P = PT["P"]
    TK = PT["T"] + 273.15

    def FNA(TK):
        return (
            166800000
            - 193080 * (TK - 273.15)
            + 186.4 * (TK - 273.15) ** 2
            - 0.071288 * ((TK - 273.15) ** 3)
        ) * 1.01325

    def FNB(TK):
        return 1.01325 * (73030000 - 71400 * (TK - 273.15) + 21.57 * (TK - 273.15) ** 2)

    def FNC(TK):
        R = 83.14321
        return 1.01325 * (
            np.exp(-11.071 + 5953 / TK - 2746000 / TK**2 + 464600000 / TK**3)
            * 0.5
            * R
            * R
            * TK**2.5
            / 1.02668
            + 40123800
        )

    def FNF(V, TK, A, B, P):
        R = 83.14321
        return R * TK / (V - B) - A / ((V * V + B * V) * TK**0.5) - P

    R = 83.14321
    B_1 = 14.6
    B_2 = 29.7

    B = X_1 * B_1 + (1 - X_1) * B_2
    A = X_1**2 * FNA(TK) + 2 * X_1 * (1 - X_1) * FNC(TK) + (1 - X_1) ** 2 * FNB(TK)
    Temp2 = B + 5
    Q = 1
    Temp1 = 0
    while abs(Temp2 - Temp1) >= 0.00001:
        Temp1 = Temp2
        F_1 = (FNF(Temp1 + 0.01, TK, A, B, P) - FNF(Temp1, TK, A, B, P)) / 0.01
        Temp2 = Temp1 - Q * FNF(Temp1, TK, A, B, P) / F_1
        F_2 = (FNF(Temp2 + 0.01, TK, A, B, P) - FNF(Temp2, TK, A, B, P)) / 0.01
        if F_2 * F_1 <= 0:
            Q = Q / 2.0
        if abs(Temp2 - Temp1) > 0.00001:
            F_1 = F_2
    V = Temp2
    G = (
        np.log(V / (V - B))
        + B_1 / (V - B)
        - 2
        * (X_1 * FNA(TK) + (1 - X_1) * FNC(TK))
        * np.log((V + B) / V)
        / (R * TK**1.5 * B)
    )
    G = (
        G
        + (np.log((V + B) / V) - B / (V + B)) * A * B_1 / (R * TK**1.5 * B**2)
        - np.log(P * V / (R * TK))
    )
    G = np.exp(G)
    return G


###########################
# Shi & Saxena (1992) #####
###########################


def lny_SS(PT, Pcr, Tcr, models):
    """
    Natural log of the fugacity coefficient using Shi & Saxena (1992)


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    Pcr: float
        Critical pressure

    Tcr: float
        Critical temperature

    models: pandas.DataFrame
        Dataframe detailing model options

    Returns
    -------
    float
        Natural log of the fugacity coefficient
    """
    P = PT["P"]
    T_K = PT["T"] + 273.15
    Tr = T_K / Tcr
    A, B, C, D, P0, integral0 = Q_SS(PT, Tr, Pcr, models)
    Pr = P / Pcr
    P0r = P0 / Pcr
    if models.loc["high precision", "option"] == "True":
        integral = (
            A * gp.log(Pr / P0r)
            + B * (Pr - P0r)
            + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
            + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
        )
    else:
        integral = (
            A * math.log(Pr / P0r)
            + B * (Pr - P0r)
            + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
            + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
        )
    integral_total = integral + integral0
    return integral_total


def Q_SS(PT, Tr, Pcr, models):
    """
    Modelling constants for Shi & Saxena (1992) - Table 1.


    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    Pcr: float
        Critical pressure

    Tcr: float
        Critical temperature

    models: pandas.DataFrame
        Dataframe detailing model options

    Returns
    -------
    str
        Calculated values for A, B, C, D, P0, integral0
    """
    P = PT["P"]

    def Q1000(Pcr):
        Pr_ = 1000.0 / Pcr
        P0r_ = 1.0 / Pcr
        A0 = 1.0
        B0 = 0.9827e-1 * pow(Tr, -1.0) + -0.2709 * pow(Tr, -3.0)
        C0 = -0.1030e-2 * pow(Tr, -1.5) + 0.1427e-1 * pow(Tr, -4.0)
        D0 = 0.0
        if models.loc["high precision", "option"] == "True":
            value = (
                A0 * gp.log(Pr_ / P0r_)
                + B0 * (Pr_ - P0r_)
                + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
            )
        else:
            value = (
                A0 * math.log(Pr_ / P0r_)
                + B0 * (Pr_ - P0r_)
                + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
            )
        return value

    def Q5000(Pcr):
        Pr_ = 5000.0 / Pcr
        P0r_ = 1000.0 / Pcr
        A0 = 1.0 + -5.917e-1 * pow(Tr, -2.0)
        B0 = 9.122e-2 * pow(Tr, -1.0)
        D0 = 0.0
        if models.loc["high precision", "option"] == "True":
            C0 = -1.416e-4 * pow(Tr, -2.0) + -2.835e-6 * gp.log(Tr)
            value = (
                A0 * gp.log(Pr_ / P0r_)
                + B0 * (Pr_ - P0r_)
                + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
            )
        else:
            C0 = -1.416e-4 * pow(Tr, -2.0) + -2.835e-6 * math.log(Tr)
            value = (
                A0 * math.log(Pr_ / P0r_)
                + B0 * (Pr_ - P0r_)
                + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
            )
        return value

    if P > 5000.0:
        if models.loc["high precision", "option"] == "True":
            A = 2.0614 + -2.235 * pow(Tr, -2.0) + -3.941e-1 * gp.log(Tr)
        else:
            A = 2.0614 + -2.235 * pow(Tr, -2.0) + -3.941e-1 * math.log(Tr)
        B = 5.513e-2 * pow(Tr, -1.0) + 3.934e-2 * pow(Tr, -2.0)
        C = (
            -1.894e-6 * pow(Tr, -1.0)
            + -1.109e-5 * pow(Tr, -2.0)
            + -2.189e-5 * pow(Tr, -3.0)
        )
        D = 5.053e-11 * pow(Tr, -1.0) + -6.303e-21 * pow(Tr, 3.0)
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
        A = 1.0 + -5.917e-1 * pow(Tr, -2.0)
        B = 9.122e-2 * pow(Tr, -1.0)
        if models.loc["high precision", "option"] == "True":
            C = -1.416e-4 * pow(Tr, -2.0) + -2.835e-6 * gp.log(Tr)
        else:
            C = -1.416e-4 * pow(Tr, -2.0) + -2.835e-6 * math.log(Tr)
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
        B = 0.9827e-1 * pow(Tr, -1.0) + -0.2709 * pow(Tr, -3.0)
        C = -0.1030e-2 * pow(Tr, -1.5) + 0.1427e-1 * pow(Tr, -4.0)
        D = 0.0
        P0 = 1.0
        integral0 = 0.0
        return A, B, C, D, P0, integral0


def y_SS(gas_species, PT, models=default_models):
    """
    Fugacity coefficient using Shi & Saxena (1992).


    Parameters
    ----------
    gas_species: string
        Name of gas species

    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Dataframe detailing model options

    Returns
    -------
    float
        Fugacity coefficient
    """
    P = PT["P"]

    ideal_gas = models.loc["ideal_gas", "option"]

    if ideal_gas == "True":
        return 1.0
    elif P < 1.0:  # ideal gas below 1 bar
        return 1.0
    else:
        Tcr = species.loc[gas_species, "Tcr"]
        Pcr = species.loc[gas_species, "Pcr"]
        if models.loc["high precision", "option"] == "True":
            value = gp.exp(lny_SS(PT, Pcr, Tcr, models)) / P
        else:
            value = math.exp(lny_SS(PT, Pcr, Tcr, models)) / P
        return value


##############################
# for each vapor species #####
##############################


def y_H2(PT, models=default_models):
    """
    Fugacity coefficient for H2

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_H2" and "ideal_gas" and column
        label of "option"

    Returns
    -------
    float
        Fugacity coefficient for H2


    Model options for y_H2
    ----------------------
    - 'Shaw64' [default] Eq. (4) from Shaw & Wones (1964) AmJSci 262:918-929
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    P = PT["P"]
    T_K = PT["T"] + 273.15

    ideal_gas = models.loc["ideal_gas", "option"]
    model = models.loc["y_H2", "option"]

    if ideal_gas == "True" or model == "ideal":
        return 1.0
    elif P < 1.0:  # ideal gas below 1 bar
        return 1.0
    elif model == "Shaw64":  # Eq. (4) from Shaw & Wones (1964) AmJSci 262:918-929
        P_atm = 0.986923 * P
        if models.loc["high precision", "option"] == "True":
            SW1 = gp.exp(-3.8402 * pow(T_K, 0.125) + 0.5410)
            SW2 = gp.exp(-0.1263 * pow(T_K, 0.5) - 15.980)
            # NB used a value of -0.011901 instead of -0.11901 as reported to match data
            # in Table 2
            SW3 = 300 * gp.exp((-0.011901 * T_K) - 5.941)
            ln_y = (
                SW1 * P_atm
                - SW2 * pow(P_atm, 2.0)
                + SW3 * gp.exp((-P_atm / 300.0) - 1.0)
            )
            value = gp.exp(ln_y)
        else:
            SW1 = math.exp(-3.8402 * pow(T_K, 0.125) + 0.5410)
            SW2 = math.exp(-0.1263 * pow(T_K, 0.5) - 15.980)
            # NB used a value of -0.011901 instead of -0.11901 as reported to match data
            # in Table 2
            SW3 = 300 * math.exp((-0.011901 * T_K) - 5.941)
            ln_y = (
                SW1 * P_atm
                - SW2 * pow(P_atm, 2.0)
                + SW3 * math.exp((-P_atm / 300.0) - 1.0)
            )
            value = math.exp(ln_y)
        return value

    # WORK IN PROGRESS
    elif model == "Shi92":  # Shi & Saxena (1992) NOT WORKING
        Tcr = 33.25  # critical temperature in K
        Pcr = 12.9696  # critical temperature in bar
        Tr = T_K / Tcr
        # Q for 1-1000 bar
        Q1_A_LP, Q2_A_LP, Q3_A_LP, Q4_A_LP, Q5_A_LP = 1.0, 0.0, 0.0, 0.0, 0.0
        Q1_B_LP, Q2_B_LP, Q3_B_LP, Q4_B_LP, Q5_B_LP = 0.0, 0.9827e-1, 0.0, -0.2709, 0.0
        Q1_C_LP, Q2_C_LP, Q3_C_LP, Q4_C_LP, Q5_C_LP = (
            0.0,
            0.0,
            -0.1030e-2,
            0.0,
            0.1427e-1,
        )
        # Q for 1000-10000 bar
        Q1_A_HP, Q2_A_HP, Q3_A_HP, Q4_A_HP, Q5_A_HP, Q6_A_HP, Q7_A_HP, Q8_A_HP = (
            2.2615,
            0.0,
            -6.8712e1,
            0.0,
            -1.0573e4,
            0.0,
            0.0,
            -1.6936e-1,
        )
        Q1_B_HP, Q2_B_HP, Q3_B_HP, Q4_B_HP, Q5_B_HP, Q6_B_HP, Q7_B_HP, Q8_B_HP = (
            -2.6707e-4,
            0.0,
            2.0173e-1,
            0.0,
            4.5759,
            0.0,
            0.0,
            3.1452e-5,
        )
        Q1_C_HP, Q2_C_HP, Q3_C_HP, Q4_C_HP, Q5_C_HP, Q6_C_HP, Q7_C_HP, Q8_C_HP = (
            -2.3376e-9,
            0.0,
            3.4091e-7,
            0.0,
            -1.4188e-3,
            0.0,
            0.0,
            3.0117e-10,
        )
        Q1_D_HP, Q2_D_HP, Q3_D_HP, Q4_D_HP, Q5_D_HP, Q6_D_HP, Q7_D_HP, Q8_D_HP = (
            -3.2606e-15,
            0.0,
            2.4402e-12,
            0.0,
            -2.4027e-9,
            0.0,
            0.0,
            0.0,
        )
        if P < 1000.0:
            A = (
                Q1_A_LP
                + Q2_A_LP * Tr ** (-1.0)
                + Q3_A_LP * Tr ** (-1.5)
                + Q4_A_LP * Tr ** (-3.0)
                + Q5_A_LP * Tr ** (-4.0)
            )
            B = (
                Q1_B_LP
                + Q2_B_LP * Tr ** (-1.0)
                + Q3_B_LP * Tr ** (-1.5)
                + Q4_B_LP * Tr ** (-3.0)
                + Q5_B_LP * Tr ** (-4.0)
            )
            C = (
                Q1_C_LP
                + Q2_C_LP * Tr ** (-1.0)
                + Q3_C_LP * Tr ** (-1.5)
                + Q4_C_LP * Tr ** (-3.0)
                + Q5_C_LP * Tr ** (-4.0)
            )
            D = 0.0
            P0 = 1.0
            integral0 = 0.0
        elif P == 1000.0:
            A = 0.0
            B = 0.0
            C = 0.0
            D = 0.0
            P0 = 1000.0
            Pr_ = 1000.0 / Pcr
            P0r_ = 1.0 / Pcr
            A0 = (
                Q1_A_LP
                + Q2_A_LP * Tr
                + Q3_A_LP * Tr ** (-1.0)
                + Q4_A_LP * Tr**2.0
                + Q5_A_LP * Tr ** (-2.0)
            )
            B0 = (
                Q1_B_LP
                + Q2_B_LP * Tr
                + Q3_B_LP * Tr ** (-1.0)
                + Q4_B_LP * Tr**2.0
                + Q5_B_LP * Tr ** (-2.0)
            )
            C0 = (
                Q1_C_LP
                + Q2_C_LP * Tr
                + Q3_C_LP * Tr ** (-1.0)
                + Q4_C_LP * Tr**2.0
                + Q5_C_LP * Tr ** (-2.0)
            )
            D0 = 0.0
            if models.loc["high precision", "option"] == "True":
                integral0 = (
                    A0 * gp.log(Pr_ / P0r_)
                    + B0 * (Pr_ - P0r_)
                    + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                    + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                )
            else:
                integral0 = (
                    A0 * math.log(Pr_ / P0r_)
                    + B0 * (Pr_ - P0r_)
                    + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                    + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                )
        elif P > 1000.0:
            if models.loc["high precision", "option"] == "True":
                A = (
                    Q1_A_HP
                    + Q2_A_HP * Tr
                    + Q3_A_HP * Tr ** (-1.0)
                    + Q4_A_HP * Tr**2.0
                    + Q5_A_HP * Tr ** (-2.0)
                    + Q6_A_HP * Tr**3.0
                    + Q7_A_HP * Tr ** (-3.0)
                    + Q8_A_HP * gp.log(Tr)
                )
                B = (
                    Q1_B_HP
                    + Q2_B_HP * Tr
                    + Q3_B_HP * Tr ** (-1.0)
                    + Q4_B_HP * Tr**2.0
                    + Q5_B_HP * Tr ** (-2.0)
                    + Q6_B_HP * Tr**3.0
                    + Q7_B_HP * Tr ** (-3.0)
                    + Q8_B_HP * gp.log(Tr)
                )
                C = (
                    Q1_C_HP
                    + Q2_C_HP * Tr
                    + Q3_C_HP * Tr ** (-1.0)
                    + Q4_C_HP * Tr**2.0
                    + Q5_C_HP * Tr ** (-2.0)
                    + Q6_C_HP * Tr**3.0
                    + Q7_C_HP * Tr ** (-3.0)
                    + Q8_C_HP * gp.log(Tr)
                )
                D = (
                    Q1_D_HP
                    + Q2_D_HP * Tr
                    + Q3_D_HP * Tr ** (-1.0)
                    + Q4_D_HP * Tr**2.0
                    + Q5_D_HP * Tr ** (-2.0)
                    + Q6_D_HP * Tr**3.0
                    + Q7_D_HP * Tr ** (-3.0)
                    + Q8_D_HP * gp.log(Tr)
                )
            else:
                A = (
                    Q1_A_HP
                    + Q2_A_HP * Tr
                    + Q3_A_HP * Tr ** (-1.0)
                    + Q4_A_HP * Tr**2.0
                    + Q5_A_HP * Tr ** (-2.0)
                    + Q6_A_HP * Tr**3.0
                    + Q7_A_HP * Tr ** (-3.0)
                    + Q8_A_HP * math.log(Tr)
                )
                B = (
                    Q1_B_HP
                    + Q2_B_HP * Tr
                    + Q3_B_HP * Tr ** (-1.0)
                    + Q4_B_HP * Tr**2.0
                    + Q5_B_HP * Tr ** (-2.0)
                    + Q6_B_HP * Tr**3.0
                    + Q7_B_HP * Tr ** (-3.0)
                    + Q8_B_HP * math.log(Tr)
                )
                C = (
                    Q1_C_HP
                    + Q2_C_HP * Tr
                    + Q3_C_HP * Tr ** (-1.0)
                    + Q4_C_HP * Tr**2.0
                    + Q5_C_HP * Tr ** (-2.0)
                    + Q6_C_HP * Tr**3.0
                    + Q7_C_HP * Tr ** (-3.0)
                    + Q8_C_HP * math.log(Tr)
                )
                D = (
                    Q1_D_HP
                    + Q2_D_HP * Tr
                    + Q3_D_HP * Tr ** (-1.0)
                    + Q4_D_HP * Tr**2.0
                    + Q5_D_HP * Tr ** (-2.0)
                    + Q6_D_HP * Tr**3.0
                    + Q7_D_HP * Tr ** (-3.0)
                    + Q8_D_HP * math.log(Tr)
                )
            P0 = 1000.0
            Pr_ = 1000.0 / Pcr
            P0r_ = 1.0 / Pcr
            A0 = (
                Q1_A_LP
                + Q2_A_LP * Tr
                + Q3_A_LP * Tr ** (-1.0)
                + Q4_A_LP * Tr**2.0
                + Q5_A_LP * Tr ** (-2.0)
            )
            B0 = (
                Q1_B_LP
                + Q2_B_LP * Tr
                + Q3_B_LP * Tr ** (-1.0)
                + Q4_B_LP * Tr**2.0
                + Q5_B_LP * Tr ** (-2.0)
            )
            C0 = (
                Q1_C_LP
                + Q2_C_LP * Tr
                + Q3_C_LP * Tr ** (-1.0)
                + Q4_C_LP * Tr**2.0
                + Q5_C_LP * Tr ** (-2.0)
            )
            D0 = 0.0
            if models.loc["high precision", "option"] == "True":
                integral0 = (
                    A0 * gp.log(Pr_ / P0r_)
                    + B0 * (Pr_ - P0r_)
                    + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                    + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                )
            else:
                integral0 = (
                    A0 * math.log(Pr_ / P0r_)
                    + B0 * (Pr_ - P0r_)
                    + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                    + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                )
        P0r = P0 / Pcr
        Pr = P / Pcr
        if models.loc["high precision", "option"] == "True":
            integral = (
                A * gp.log(Pr / P0r)
                + B * (Pr - P0r)
                + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
            )
            value = gp.exp(integral + integral0) / P
        else:
            integral = (
                A * math.log(Pr / P0r)
                + B * (Pr - P0r)
                + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
            )
            value = math.exp(integral + integral0) / P
        return value


def y_H2O(PT, models=default_models):
    """
    Fugacity coefficient for H2O

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_H2O" and "ideal_gas" and
        column label of "option"

    Returns
    -------
    float
        Fugacity coefficient for H2O


    Model options for y_H2O
    -----------------------
    - 'Holland91' [default] Eq. (4,6,A1-3) and Table 1 (T > 673 K only) from Holland &
    Powell (1991) CMP 109:265-273 10.1007/BF00306484
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    ideal_gas = models.loc["ideal_gas", "option"]
    model = models.loc["y_H2O", "option"]

    P = PT["P"]

    if ideal_gas == "True" or model == "ideal":
        return 1.0
    elif P < 1.0:  # ideal gas below 1 bar
        return 1.0
    else:
        # Eq. (4,6,A1-3) and Table 1 (T > 673 K only) from Holland & Powell (1991) CMP
        # 109:265-273 10.1007/BF00306484
        if model == "Holland91":
            # Eq. (4,A1-3)
            y = y_CORK("H2O", PT, models)
        elif model == "Flowers79":
            y = MRK(PT, 1.0)
        return y


def y_CO2(PT, models=default_models):
    """
    Fugacity coefficient for CO2

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_CO2" and "ideal_gas" and
        column label of "option"

    Returns
    -------
    float
        Fugacity coefficient for CO2


    Model options for y_CO2
    -----------------------
    - 'Holland91_eq8_tab1' [default] Eq. (8) and Table 1 from Holland & Powell (1991)
    CMP 109:265-273 10.1007/BF00306484
    - 'Holland91_eq4,A1-3_tab1' Eq. (4,A1-3) and Table 1 from Holland & Powell (1991)
    CMP 109:265-273 10.1007/BF00306484
    - 'Holland91_eq8,9_tab2' Eq. (8,9) and Table 2 from Holland & Powell (1991) CMP
    109:265-273 10.1007/BF00306484
    - 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    ideal_gas = models.loc["ideal_gas", "option"]
    model = models.loc["y_CO2", "option"]

    P = PT["P"]

    if ideal_gas == "True" or model == "ideal":
        return 1.0
    elif P < 1.0:  # ideal gas below 1 bar
        return 1.0
    else:
        # Eq. (4,A1-3) and Table 1 from Holland & Powell (1991) CMP 109:265-273
        if model == "Holland91_eq4,A1-3_tab1":
            y = y_CORK("CO2", PT, models)  # Eq. (4,A1-3)
        # Eq. (8,9) and Table 2 from Holland & Powell (1991) CMP 109:265-273
        elif model == "Holland91_eq8,9_tab2":
            y = y_sCORK("CO2", PT, models)  # Eq. (8)
        # Eq. (8) and Table 1 from Holland & Powell (1991) CMP 109:265-273
        elif model == "Holland91_eq8_tab1":
            y = y_sCORK("CO2", PT, models)  # Eq. (8)
        # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        elif model == "Shi92":
            gas_species = "CO2"
            y = y_SS(gas_species, PT, models)
        if model == "Flowers79":
            y = MRK(PT, 0.0)
        return y


def y_O2(PT, models=default_models):
    """
    Fugacity coefficient for O2

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_O2" and "ideal_gas" and column
        label of "option"

    Returns
    -------
    float
        Fugacity coefficient for O2


    Model options for y_O2
    ----------------------
    - 'Shi92' [default] Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_O2", "option"]

    if model == "Shi92":  # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "O2"
        y = y_SS(gas_species, PT, models)
    elif model == "ideal":
        y = 1.0
    return y


def y_S2(PT, models=default_models):
    """
    Fugacity coefficient for S2

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_S2" and "ideal_gas" and column
        label of "option"

    Returns
    -------
    float
        Fugacity coefficient for S2


    Model options for y_S2
    ----------------------
    - 'Shi92' [default] Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_S2", "option"]

    if model == "Shi92":  # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "S2"
        y = y_SS(gas_species, PT, models)
    elif model == "ideal":
        y = 1.0
    return y


def y_CO(PT, models=default_models):
    """
    Fugacity coefficient for CO

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_CO" and "ideal_gas" and column
        label of "option"

    Returns
    -------
    float
        Fugacity coefficient for CO


    Model options for y_CO
    ----------------------
    - 'Shi92' [default] Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_CO", "option"]

    if model == "Shi92":  # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "CO"
        y = y_SS(gas_species, PT, models)
    elif model == "ideal":
        y = 1.0
    return y


def y_CH4(PT, models=default_models):
    """
    Fugacity coefficient for CH4

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_CH4" and "ideal_gas" and
        column label of "option"

    Returns
    -------
    float
        Fugacity coefficient for CH4


    Model options for y_CH4
    -----------------------
    - 'Shi92' [default] Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_CH4", "option"]

    if model == "Shi92":  # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "CH4"
        y = y_SS(gas_species, PT, models)
    elif model == "ideal":
        y = 1.0
    return y


def y_OCS(PT, models=default_models):
    """
    Fugacity coefficient for OCS

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_OCS" and "ideal_gas" and
        column label of "option"

    Returns
    -------
    float
        Fugacity coefficient for OCS


    Model options for y_OCS
    -----------------------
    - 'Shi92' [default] Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_OCS", "option"]

    if model == "Shi92":  # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
        gas_species = "OCS"
        y = y_SS(gas_species, PT, models)
    elif model == "ideal":
        y = 1.0
    return y


def y_X(PT, models=default_models):  # species X fugacity coefficient
    """
    Fugacity coefficient for X

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_X" and "ideal_gas" and column
        label of "option"

    Returns
    -------
    float
        Fugacity coefficient for X


    Model options for y_X
    -----------------------
    - "ideal" Treat as ideal gas, y = 1 at all P.
    - Only one option available currently, included for future development.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    model = models.loc["y_X", "option"]

    if model == "ideal":  # ideal gas
        y = 1.0
    return y


#################################################################################
# SO2 and H2S from Shi & Saxena (1992) with option to modify below 500 bars #####
#################################################################################


def y_SO2(PT, models=default_models):
    """
    Fugacity coefficient for SO2

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_SO2" and "ideal_gas" and
        column label of "option"

    Returns
    -------
    float
        Fugacity coefficient for SO2


    Model options for y_SO2
    -----------------------
    - 'Shi92_Hughes23' [default] Fig.S1 modified from Shi & Saxena (1992) from Hughes et
    al. (2023) JGSL 180(3) doi:10.1144/jgs2021-12
    - 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    ideal_gas = models.loc["ideal_gas", "option"]
    model = models.loc["y_SO2", "option"]

    P = PT["P"]
    T_K = PT["T"] + 273.15

    gas_species = "SO2"
    if ideal_gas == "True" or model == "ideal":
        return 1.0
    elif P < 1.0:  # ideal gas below 1 bar
        return 1.0
    else:  # 1-10000 bar
        Tcr = species.loc[gas_species, "Tcr"]  # critical temperature in K
        Pcr = species.loc[gas_species, "Pcr"]  # critical temperature in bar
        P0 = 1.0
        P0r = P0 / Pcr
        Tr = T_K / Tcr
        Q1_A, Q2_A, Q3_A, Q4_A, Q5_A, Q6_A, Q7_A, Q8_A = (
            0.92854,
            0.43269e-1,
            -0.24671,
            0.0,
            0.24999,
            0.0,
            -0.53182,
            -0.16461e-1,
        )
        Q1_B, Q2_B, Q3_B, Q4_B, Q5_B, Q6_B, Q7_B, Q8_B = (
            0.84866e-3,
            -0.18379e-2,
            0.66787e-1,
            0.0,
            -0.29427e-1,
            0.0,
            0.29003e-1,
            0.54808e-2,
        )
        Q1_C, Q2_C, Q3_C, Q4_C, Q5_C, Q6_C, Q7_C, Q8_C = (
            -0.35456e-3,
            0.23316e-4,
            0.94159e-3,
            0.0,
            -0.81653e-3,
            0.0,
            0.23154e-3,
            0.55542e-4,
        )
        if models.loc["high precision", "option"] == "True":
            A = (
                Q1_A
                + Q2_A * Tr
                + Q3_A * Tr ** (-1.0)
                + Q4_A * Tr**2.0
                + Q5_A * Tr ** (-2.0)
                + Q6_A * Tr**3.0
                + Q7_A * Tr ** (-3.0)
                + Q8_A * gp.log(Tr)
            )
            B = (
                Q1_B
                + Q2_B * Tr
                + Q3_B * Tr ** (-1.0)
                + Q4_B * Tr**2.0
                + Q5_B * Tr ** (-2.0)
                + Q6_B * Tr**3.0
                + Q7_B * Tr ** (-3.0)
                + Q8_B * gp.log(Tr)
            )
            C = (
                Q1_C
                + Q2_C * Tr
                + Q3_C * Tr ** (-1.0)
                + Q4_C * Tr**2.0
                + Q5_C * Tr ** (-2.0)
                + Q6_C * Tr**3.0
                + Q7_C * Tr ** (-3.0)
                + Q8_C * gp.log(Tr)
            )
        else:
            A = (
                Q1_A
                + Q2_A * Tr
                + Q3_A * Tr ** (-1.0)
                + Q4_A * Tr**2.0
                + Q5_A * Tr ** (-2.0)
                + Q6_A * Tr**3.0
                + Q7_A * Tr ** (-3.0)
                + Q8_A * math.log(Tr)
            )
            B = (
                Q1_B
                + Q2_B * Tr
                + Q3_B * Tr ** (-1.0)
                + Q4_B * Tr**2.0
                + Q5_B * Tr ** (-2.0)
                + Q6_B * Tr**3.0
                + Q7_B * Tr ** (-3.0)
                + Q8_B * math.log(Tr)
            )
            C = (
                Q1_C
                + Q2_C * Tr
                + Q3_C * Tr ** (-1.0)
                + Q4_C * Tr**2.0
                + Q5_C * Tr ** (-2.0)
                + Q6_C * Tr**3.0
                + Q7_C * Tr ** (-3.0)
                + Q8_C * math.log(Tr)
            )
        D = 0.0
        if P >= 500.0:  # above 500 bar Shi & Saxena (1992) AmMin 77(9-10):1038-1049
            Pr = P / Pcr
            if models.loc["high precision", "option"] == "True":
                integral = (
                    A * gp.log(Pr / P0r)
                    + B * (Pr - P0r)
                    + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                    + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
                )
                y = (gp.exp(integral)) / P
            else:
                integral = (
                    A * math.log(Pr / P0r)
                    + B * (Pr - P0r)
                    + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                    + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
                )
                y = (math.exp(integral)) / P
        elif (
            models.loc["y_SO2", "option"] == "Shi92"
        ):  # Shi & Saxena (1992) AmMin 77(9-10):1038-1049
            Pr = P / Pcr
            if models.loc["high precision", "option"] == "True":
                integral = (
                    A * gp.log(Pr / P0r)
                    + B * (Pr - P0r)
                    + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                    + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
                )
                y = (gp.exp(integral)) / P
            else:
                integral = (
                    A * math.log(Pr / P0r)
                    + B * (Pr - P0r)
                    + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                    + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
                )
                y = (math.exp(integral)) / P
        elif (
            models.loc["y_SO2", "option"] == "Shi92_Hughes23"
        ):  # Fig.S1 Hughes et al. (2023) JGSL 180(3) doi:10.1144/jgs2021-12
            Pr = 500.0 / Pcr  # calculate y at 500 bar
            if models.loc["high precision", "option"] == "True":
                integral = (
                    A * gp.log(Pr / P0r)
                    + B * (Pr - P0r)
                    + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                    + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
                )
                y_500 = (gp.exp(integral)) / 500.0
            else:
                integral = (
                    A * math.log(Pr / P0r)
                    + B * (Pr - P0r)
                    + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                    + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
                )
                y_500 = (math.exp(integral)) / 500.0
            y = (
                (y_500 - 1.0) * (P / 500.0)
            ) + 1.0  # linear extrapolation to P of interest
        return y


def y_H2S(PT, models=default_models):
    """
    Fugacity coefficient for H2S

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "y_H2S" and "ideal_gas" and
        column label of "option"

    Returns
    -------
    float
        Fugacity coefficient for H2S


    Model options for y_H2S
    -----------------------
    - 'Shi92_Hughes23' [default] Fig.S1 modified from Shi & Saxena (1992) Hughes et al.
    (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
    - 'Shi92' Shi & Saxena (1992) AmMin 77(9-10):1038-1049
    - 'ideal' Treat as ideal gas, y = 1 at all P.
    (Note: "ideal_gas" = "True" overides chosen option)
    """
    ideal_gas = models.loc["ideal_gas", "option"]
    model = models.loc["y_H2S", "option"]

    P = PT["P"]
    T_K = PT["T"] + 273.15

    gas_species = "H2S"
    if ideal_gas == "True" or model == "ideal":
        return 1.0
    elif ideal_gas == "False":
        Tcr = species.loc[gas_species, "Tcr"]  # critical temperature in K
        Pcr = species.loc[gas_species, "Pcr"]  # critical temperature in bar
        Tr = T_K / Tcr
        # Q for 1-500 bar
        Q1_A_LP, Q2_A_LP, Q3_A_LP, Q4_A_LP, Q5_A_LP, Q6_A_LP, Q7_A_LP, Q8_A_LP = (
            0.14721e1,
            0.11177e1,
            0.39657e1,
            0.0,
            -0.10028e2,
            0.0,
            0.45484e1,
            -0.382e1,
        )
        Q1_B_LP, Q2_B_LP, Q3_B_LP, Q4_B_LP, Q5_B_LP, Q6_B_LP, Q7_B_LP, Q8_B_LP = (
            0.16066,
            0.10887,
            0.29014,
            0.0,
            -0.99593,
            0.0,
            -0.18627,
            -0.45515,
        )
        Q1_C_LP, Q2_C_LP, Q3_C_LP, Q4_C_LP, Q5_C_LP, Q6_C_LP, Q7_C_LP, Q8_C_LP = (
            -0.28933,
            -0.70522e-1,
            0.39828,
            0.0,
            -0.50533e-1,
            0.0,
            0.1176,
            0.33972,
        )
        # Q for 500-10000 bar
        Q1_A_HP, Q2_A_HP, Q3_A_HP, Q4_A_HP, Q5_A_HP, Q6_A_HP, Q7_A_HP, Q8_A_HP = (
            0.59941,
            -0.1557e-2,
            0.4525e-1,
            0.0,
            0.36687,
            0.0,
            -0.79248,
            0.26058,
        )
        Q1_B_HP, Q2_B_HP, Q3_B_HP, Q4_B_HP, Q5_B_HP, Q6_B_HP, Q7_B_HP, Q8_B_HP = (
            0.22545e-1,
            0.17473e-2,
            0.48253e-1,
            0.0,
            -0.1989e-1,
            0.0,
            0.32794e-1,
            -0.10985e-1,
        )
        Q1_C_HP, Q2_C_HP, Q3_C_HP, Q4_C_HP, Q5_C_HP, Q6_C_HP, Q7_C_HP, Q8_C_HP = (
            0.57375e-3,
            -0.20944e-5,
            -0.11894e-2,
            0.0,
            0.14661e-2,
            0.0,
            -0.75605e-3,
            -0.27985e-3,
        )
        if P < 1.0:
            return 1.0  # ideal gas below 1 bar
        elif P < 500.0:
            if models.loc["y_H2S", "option"] == "Shi92":  # as is Shi and Saxena (1992)
                if models.loc["high precision", "option"] == "True":
                    A = (
                        Q1_A_LP
                        + Q2_A_LP * Tr
                        + Q3_A_LP * Tr ** (-1.0)
                        + Q4_A_LP * Tr**2.0
                        + Q5_A_LP * Tr ** (-2.0)
                        + Q6_A_LP * Tr**3.0
                        + Q7_A_LP * Tr ** (-3.0)
                        + Q8_A_LP * gp.log(Tr)
                    )
                    B = (
                        Q1_B_LP
                        + Q2_B_LP * Tr
                        + Q3_B_LP * Tr ** (-1.0)
                        + Q4_B_LP * Tr**2.0
                        + Q5_B_LP * Tr ** (-2.0)
                        + Q6_B_LP * Tr**3.0
                        + Q7_B_LP * Tr ** (-3.0)
                        + Q8_B_LP * gp.log(Tr)
                    )
                    C = (
                        Q1_C_LP
                        + Q2_C_LP * Tr
                        + Q3_C_LP * Tr ** (-1.0)
                        + Q4_C_LP * Tr**2.0
                        + Q5_C_LP * Tr ** (-2.0)
                        + Q6_C_LP * Tr**3.0
                        + Q7_C_LP * Tr ** (-3.0)
                        + Q8_C_LP * gp.log(Tr)
                    )
                else:
                    A = (
                        Q1_A_LP
                        + Q2_A_LP * Tr
                        + Q3_A_LP * Tr ** (-1.0)
                        + Q4_A_LP * Tr**2.0
                        + Q5_A_LP * Tr ** (-2.0)
                        + Q6_A_LP * Tr**3.0
                        + Q7_A_LP * Tr ** (-3.0)
                        + Q8_A_LP * math.log(Tr)
                    )
                    B = (
                        Q1_B_LP
                        + Q2_B_LP * Tr
                        + Q3_B_LP * Tr ** (-1.0)
                        + Q4_B_LP * Tr**2.0
                        + Q5_B_LP * Tr ** (-2.0)
                        + Q6_B_LP * Tr**3.0
                        + Q7_B_LP * Tr ** (-3.0)
                        + Q8_B_LP * math.log(Tr)
                    )
                    C = (
                        Q1_C_LP
                        + Q2_C_LP * Tr
                        + Q3_C_LP * Tr ** (-1.0)
                        + Q4_C_LP * Tr**2.0
                        + Q5_C_LP * Tr ** (-2.0)
                        + Q6_C_LP * Tr**3.0
                        + Q7_C_LP * Tr ** (-3.0)
                        + Q8_C_LP * math.log(Tr)
                    )
                D = 0.0
                P0 = 1.0
                integral0 = 0.0

            # Fig.S1 modified from Shi & Saxena (1992) Hughes et al. (2024)
            # AmMin 109(3):422-438 doi:10.2138/am-2023-8739
            elif models.loc["y_H2S", "option"] == "Shi92_Hughes24":
                P0 = 500.0  # calculate y at 500 bars
                Pr_ = 500.0 / Pcr
                P0r_ = 1.0 / Pcr
                D0 = 0.0
                if models.loc["high precision", "option"] == "True":
                    A0 = (
                        Q1_A_LP
                        + Q2_A_LP * Tr
                        + Q3_A_LP * Tr ** (-1.0)
                        + Q4_A_LP * Tr**2.0
                        + Q5_A_LP * Tr ** (-2.0)
                        + Q6_A_LP * Tr**3.0
                        + Q7_A_LP * Tr ** (-3.0)
                        + Q8_A_LP * gp.log(Tr)
                    )
                    B0 = (
                        Q1_B_LP
                        + Q2_B_LP * Tr
                        + Q3_B_LP * Tr ** (-1.0)
                        + Q4_B_LP * Tr**2.0
                        + Q5_B_LP * Tr ** (-2.0)
                        + Q6_B_LP * Tr**3.0
                        + Q7_B_LP * Tr ** (-3.0)
                        + Q8_B_LP * gp.log(Tr)
                    )
                    C0 = (
                        Q1_C_LP
                        + Q2_C_LP * Tr
                        + Q3_C_LP * Tr ** (-1.0)
                        + Q4_C_LP * Tr**2.0
                        + Q5_C_LP * Tr ** (-2.0)
                        + Q6_C_LP * Tr**3.0
                        + Q7_C_LP * Tr ** (-3.0)
                        + Q8_C_LP * gp.log(Tr)
                    )
                    integral0 = (
                        A0 * gp.log(Pr_ / P0r_)
                        + B0 * (Pr_ - P0r_)
                        + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                        + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                    )
                    y_500 = gp.exp(integral0) / 500.0
                else:
                    A0 = (
                        Q1_A_LP
                        + Q2_A_LP * Tr
                        + Q3_A_LP * Tr ** (-1.0)
                        + Q4_A_LP * Tr**2.0
                        + Q5_A_LP * Tr ** (-2.0)
                        + Q6_A_LP * Tr**3.0
                        + Q7_A_LP * Tr ** (-3.0)
                        + Q8_A_LP * math.log(Tr)
                    )
                    B0 = (
                        Q1_B_LP
                        + Q2_B_LP * Tr
                        + Q3_B_LP * Tr ** (-1.0)
                        + Q4_B_LP * Tr**2.0
                        + Q5_B_LP * Tr ** (-2.0)
                        + Q6_B_LP * Tr**3.0
                        + Q7_B_LP * Tr ** (-3.0)
                        + Q8_B_LP * math.log(Tr)
                    )
                    C0 = (
                        Q1_C_LP
                        + Q2_C_LP * Tr
                        + Q3_C_LP * Tr ** (-1.0)
                        + Q4_C_LP * Tr**2.0
                        + Q5_C_LP * Tr ** (-2.0)
                        + Q6_C_LP * Tr**3.0
                        + Q7_C_LP * Tr ** (-3.0)
                        + Q8_C_LP * math.log(Tr)
                    )
                    integral0 = (
                        A0 * math.log(Pr_ / P0r_)
                        + B0 * (Pr_ - P0r_)
                        + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                        + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                    )
                    y_500 = math.exp(integral0) / 500.0
                y = (
                    (y_500 - 1.0) * (P / 500.0)
                ) + 1.0  # linear extrapolation to P of interest
                return y
        elif P == 500.0:
            A = 0.0
            B = 0.0
            C = 0.0
            D = 0.0
            P0 = 500.0
            Pr_ = 500.0 / Pcr
            P0r_ = 1.0 / Pcr
            D0 = 0.0
            if models.loc["high precision", "option"] == "True":
                A0 = (
                    Q1_A_LP
                    + Q2_A_LP * Tr
                    + Q3_A_LP * Tr ** (-1.0)
                    + Q4_A_LP * Tr**2.0
                    + Q5_A_LP * Tr ** (-2.0)
                    + Q6_A_LP * Tr**3.0
                    + Q7_A_LP * Tr ** (-3.0)
                    + Q8_A_LP * gp.log(Tr)
                )
                B0 = (
                    Q1_B_LP
                    + Q2_B_LP * Tr
                    + Q3_B_LP * Tr ** (-1.0)
                    + Q4_B_LP * Tr**2.0
                    + Q5_B_LP * Tr ** (-2.0)
                    + Q6_B_LP * Tr**3.0
                    + Q7_B_LP * Tr ** (-3.0)
                    + Q8_B_LP * gp.log(Tr)
                )
                C0 = (
                    Q1_C_LP
                    + Q2_C_LP * Tr
                    + Q3_C_LP * Tr ** (-1.0)
                    + Q4_C_LP * Tr**2.0
                    + Q5_C_LP * Tr ** (-2.0)
                    + Q6_C_LP * Tr**3.0
                    + Q7_C_LP * Tr ** (-3.0)
                    + Q8_C_LP * gp.log(Tr)
                )
                integral0 = (
                    A0 * gp.log(Pr_ / P0r_)
                    + B0 * (Pr_ - P0r_)
                    + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                    + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                )
            else:
                A0 = (
                    Q1_A_LP
                    + Q2_A_LP * Tr
                    + Q3_A_LP * Tr ** (-1.0)
                    + Q4_A_LP * Tr**2.0
                    + Q5_A_LP * Tr ** (-2.0)
                    + Q6_A_LP * Tr**3.0
                    + Q7_A_LP * Tr ** (-3.0)
                    + Q8_A_LP * math.log(Tr)
                )
                B0 = (
                    Q1_B_LP
                    + Q2_B_LP * Tr
                    + Q3_B_LP * Tr ** (-1.0)
                    + Q4_B_LP * Tr**2.0
                    + Q5_B_LP * Tr ** (-2.0)
                    + Q6_B_LP * Tr**3.0
                    + Q7_B_LP * Tr ** (-3.0)
                    + Q8_B_LP * math.log(Tr)
                )
                C0 = (
                    Q1_C_LP
                    + Q2_C_LP * Tr
                    + Q3_C_LP * Tr ** (-1.0)
                    + Q4_C_LP * Tr**2.0
                    + Q5_C_LP * Tr ** (-2.0)
                    + Q6_C_LP * Tr**3.0
                    + Q7_C_LP * Tr ** (-3.0)
                    + Q8_C_LP * math.log(Tr)
                )
                integral0 = (
                    A0 * math.log(Pr_ / P0r_)
                    + B0 * (Pr_ - P0r_)
                    + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                    + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                )
        elif P > 500.0:
            D = 0.0
            P0 = 500.0
            Pr_ = 500.0 / Pcr
            P0r_ = 1.0 / Pcr
            D0 = 0.0
            if models.loc["high precision", "option"] == "True":
                A = (
                    Q1_A_HP
                    + Q2_A_HP * Tr
                    + Q3_A_HP * Tr ** (-1.0)
                    + Q4_A_HP * Tr**2.0
                    + Q5_A_HP * Tr ** (-2.0)
                    + Q6_A_HP * Tr**3.0
                    + Q7_A_HP * Tr ** (-3.0)
                    + Q8_A_HP * gp.log(Tr)
                )
                B = (
                    Q1_B_HP
                    + Q2_B_HP * Tr
                    + Q3_B_HP * Tr ** (-1.0)
                    + Q4_B_HP * Tr**2.0
                    + Q5_B_HP * Tr ** (-2.0)
                    + Q6_B_HP * Tr**3.0
                    + Q7_B_HP * Tr ** (-3.0)
                    + Q8_B_HP * gp.log(Tr)
                )
                C = (
                    Q1_C_HP
                    + Q2_C_HP * Tr
                    + Q3_C_HP * Tr ** (-1.0)
                    + Q4_C_HP * Tr**2.0
                    + Q5_C_HP * Tr ** (-2.0)
                    + Q6_C_HP * Tr**3.0
                    + Q7_C_HP * Tr ** (-3.0)
                    + Q8_C_HP * gp.log(Tr)
                )
                A0 = (
                    Q1_A_LP
                    + Q2_A_LP * Tr
                    + Q3_A_LP * Tr ** (-1.0)
                    + Q4_A_LP * Tr**2.0
                    + Q5_A_LP * Tr ** (-2.0)
                    + Q6_A_LP * Tr**3.0
                    + Q7_A_LP * Tr ** (-3.0)
                    + Q8_A_LP * gp.log(Tr)
                )
                B0 = (
                    Q1_B_LP
                    + Q2_B_LP * Tr
                    + Q3_B_LP * Tr ** (-1.0)
                    + Q4_B_LP * Tr**2.0
                    + Q5_B_LP * Tr ** (-2.0)
                    + Q6_B_LP * Tr**3.0
                    + Q7_B_LP * Tr ** (-3.0)
                    + Q8_B_LP * gp.log(Tr)
                )
                C0 = (
                    Q1_C_LP
                    + Q2_C_LP * Tr
                    + Q3_C_LP * Tr ** (-1.0)
                    + Q4_C_LP * Tr**2.0
                    + Q5_C_LP * Tr ** (-2.0)
                    + Q6_C_LP * Tr**3.0
                    + Q7_C_LP * Tr ** (-3.0)
                    + Q8_C_LP * gp.log(Tr)
                )
                integral0 = (
                    A0 * gp.log(Pr_ / P0r_)
                    + B0 * (Pr_ - P0r_)
                    + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                    + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                )
            else:
                A = (
                    Q1_A_HP
                    + Q2_A_HP * Tr
                    + Q3_A_HP * Tr ** (-1.0)
                    + Q4_A_HP * Tr**2.0
                    + Q5_A_HP * Tr ** (-2.0)
                    + Q6_A_HP * Tr**3.0
                    + Q7_A_HP * Tr ** (-3.0)
                    + Q8_A_HP * math.log(Tr)
                )
                B = (
                    Q1_B_HP
                    + Q2_B_HP * Tr
                    + Q3_B_HP * Tr ** (-1.0)
                    + Q4_B_HP * Tr**2.0
                    + Q5_B_HP * Tr ** (-2.0)
                    + Q6_B_HP * Tr**3.0
                    + Q7_B_HP * Tr ** (-3.0)
                    + Q8_B_HP * math.log(Tr)
                )
                C = (
                    Q1_C_HP
                    + Q2_C_HP * Tr
                    + Q3_C_HP * Tr ** (-1.0)
                    + Q4_C_HP * Tr**2.0
                    + Q5_C_HP * Tr ** (-2.0)
                    + Q6_C_HP * Tr**3.0
                    + Q7_C_HP * Tr ** (-3.0)
                    + Q8_C_HP * math.log(Tr)
                )
                A0 = (
                    Q1_A_LP
                    + Q2_A_LP * Tr
                    + Q3_A_LP * Tr ** (-1.0)
                    + Q4_A_LP * Tr**2.0
                    + Q5_A_LP * Tr ** (-2.0)
                    + Q6_A_LP * Tr**3.0
                    + Q7_A_LP * Tr ** (-3.0)
                    + Q8_A_LP * math.log(Tr)
                )
                B0 = (
                    Q1_B_LP
                    + Q2_B_LP * Tr
                    + Q3_B_LP * Tr ** (-1.0)
                    + Q4_B_LP * Tr**2.0
                    + Q5_B_LP * Tr ** (-2.0)
                    + Q6_B_LP * Tr**3.0
                    + Q7_B_LP * Tr ** (-3.0)
                    + Q8_B_LP * math.log(Tr)
                )
                C0 = (
                    Q1_C_LP
                    + Q2_C_LP * Tr
                    + Q3_C_LP * Tr ** (-1.0)
                    + Q4_C_LP * Tr**2.0
                    + Q5_C_LP * Tr ** (-2.0)
                    + Q6_C_LP * Tr**3.0
                    + Q7_C_LP * Tr ** (-3.0)
                    + Q8_C_LP * math.log(Tr)
                )
                integral0 = (
                    A0 * math.log(Pr_ / P0r_)
                    + B0 * (Pr_ - P0r_)
                    + (C0 / 2.0) * (pow(Pr_, 2.0) - pow(P0r_, 2.0))
                    + (D0 / 3.0) * (pow(Pr_, 3.0) - pow(P0r_, 3.0))
                )
        P0r = P0 / Pcr
        Pr = P / Pcr
        if models.loc["high precision", "option"] == "True":
            integral = (
                A * gp.log(Pr / P0r)
                + B * (Pr - P0r)
                + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
            )
            value = gp.exp(integral + integral0) / P
        else:
            integral = (
                A * math.log(Pr / P0r)
                + B * (Pr - P0r)
                + (C / 2.0) * (pow(Pr, 2.0) - pow(P0r, 2.0))
                + (D / 3.0) * (pow(Pr, 3.0) - pow(P0r, 3.0))
            )
            value = math.exp(integral + integral0) / P
        return value


def y_SO3(PT, models=default_models):
    return 1.0


#######################
# oxygen fugacity #####
#######################


# buffers
def NNO(PT, models=default_models):
    """
    Value of NNO buffer

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "NNObuffer" and column label of
        "option"

    Returns
    -------
    float
        log10fO2 value of NNO buffer as <class 'mpfr'>


    Model options for NNObuffer
    -------------
    - 'Frost91' [default] Frost (1991) in "Oxide Minerals: Petrologic and Magnetic
    Significance" doi:10.1515/9781501508684-004
    - Only one option available currently, included for future development.

    """

    model = models.loc["NNObuffer", "option"]

    P = PT["P"]
    T_K = PT["T"] + 273.15

    # Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance"
    # doi:10.1515/9781501508684-004
    if model == "Frost91":
        buffer = -24930 / T_K + 9.36 + 0.046 * (P - 1.0) / T_K
    return buffer


def FMQ(PT, models=default_models):
    """
    Value of FMQ buffer

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "FMQbuffer" and column label of
        "option"

    Returns
    -------
    float
        log10fO2 value of FMQ buffer as <class 'mpfr'>


    Model options for FMQbuffer
    -------------
    - 'Frost91' Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance"
    doi:10.1515/9781501508684-004
    - 'ONeill87' O'Neill (1897) AmMin 72(1-2):67-75

    """

    model = models.loc["FMQbuffer", "option"]

    P = PT["P"]
    T_K = PT["T"] + 273.15

    # Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance"
    # doi:10.1515/9781501508684-004
    if model == "Frost91":
        buffer = -25096.3 / T_K + 8.735 + 0.11 * (P - 1.0) / T_K

    # O'Neill (1897) AmMin 72(1-2):67-75
    elif model == "ONeill87":
        buffer = 8.58 - (25050 / T_K)
    return buffer


# terms for different equations


# terms for Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
def FefO2_KC91_Eq7_terms(PT, melt_wf, models=default_models):
    """
    Terms for Kress & Carmichael (1991) equation (7) to calculate fO2-Fe3+/FeT

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Dataframe of model options

    Returns
    -------
    str
        Terms a, B
    """

    # ln(XFe2O3/XFeO) = alnfO2 + (b/T) + c + sum(dX) + e[1 - (T0/T) = ln(T/T0)] + f(P/T)
    # + g(((T-T0)P)/T) + h(P2/T)
    # ln(XFe2O3/XFeO) = alnfO2 + B
    # terms
    a = 0.196
    # sum(dX)
    # mole frations in the melt based on oxide components (all Fe as FeO) with no
    # volatiles
    melt_comp = mg.melt_mole_fraction(melt_wf, models, "no", "no")
    DAl = -2.243
    DFe = -1.828
    DCa = 3.201
    DNa = 5.854
    DK = 6.215
    D4X = (
        DAl * melt_comp["Al2O3"]
        + DFe * melt_comp["FeOT"]
        + DCa * melt_comp["CaO"]
        + DNa * melt_comp["Na2O"]
        + DK * melt_comp["K2O"]
    )
    # PT term
    P = PT["P"]
    T_K = PT["T"] + 273.15
    b = 1.1492e4  # K
    c = -6.675
    e = -3.36
    f = -7.01e-7  # K/Pa
    g = -1.54e-10  # /Pa
    h = 3.85e-17  # K/Pa2
    T0 = 1673.0  # K
    P_Pa = P * 1.0e5  # converts bars to pascals
    B = (
        (b / T_K)
        + c
        + D4X
        + e * (1.0 - (T0 / T_K) - math.log(T_K / T0))
        + f * (P_Pa / T_K)
        + g * (((T_K - T0) * P_Pa) / T_K)
        + h * ((P_Pa**2.0) / T_K)
    )
    return a, B


# terms for Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92
# doi:10.1007/BF00307328
def FefO2_KC91_EqA_terms(PT, melt_wf, models=default_models):
    """
    Terms for Kress & Carmichael (1991) equation (A-5,6) to calculate fO2-Fe3+/FeT

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Dataframe of model options

    Returns
    -------
    str
        Terms KD1, KD2, y
    """
    # XFeO1.5/XFeO = (KD1*fO2**0.25 + 2y*KD2*KD1**2y*fO2**0.5y)/(1 +
    # (1-2y)KD2*KD1**2y*fO2**0.5y)
    KD2 = 0.4
    y = 0.3
    # compositional term
    # mole frations in the melt based on oxide components (all Fe as FeO) with no
    # volatiles
    melt_comp = mg.melt_mole_fraction(melt_wf, models, "no", "no")
    DWAl = 39.86e3  # J
    DWCa = -62.52e3  # J
    DWNa = -102.0e3  # J
    DWK = -119.0e3  # J
    D4X = (
        DWAl * melt_comp["Al2O3"]
        + DWCa * melt_comp["CaO"]
        + DWNa * melt_comp["Na2O"]
        + DWK * melt_comp["K2O"]
    )
    # KD1
    T_K = PT["T"] + 273.15
    P = PT["P"]
    DH = -106.2e3  # J
    DS = -55.10  # J/K
    DCp = 31.86  # J/K
    DV = 7.42e-6  # m3
    DVdot = 1.63e-9  # m3/K
    DVdash = -8.16e-16  # m3/Pa
    T0 = 1673.0  # K
    P0 = 1.0e5  # Pa
    R = 8.3144598  # J/K/mol
    P_Pa = P * 1.0e5
    if models.loc["high precision", "option"] == "True":
        KD1 = math.exp(
            (-DH / (R * T_K))
            + (DS / R)
            - (DCp / R) * (1.0 - (T0 / T_K) - gp.log(T_K / T0))
            - (1.0 / (R * T_K)) * D4X
            - ((DV * (P_Pa - P0)) / (R * T_K))
            - ((DVdot * (T_K - T0) * (P_Pa - P0)) / (R * T_K))
            - (DVdash / (2.0 * R * T_K)) * pow((P_Pa - P0), 2.0)
        )
    else:
        KD1 = math.exp(
            (-DH / (R * T_K))
            + (DS / R)
            - (DCp / R) * (1.0 - (T0 / T_K) - math.log(T_K / T0))
            - (1.0 / (R * T_K)) * D4X
            - ((DV * (P_Pa - P0)) / (R * T_K))
            - ((DVdot * (T_K - T0) * (P_Pa - P0)) / (R * T_K))
            - (DVdash / (2.0 * R * T_K)) * pow((P_Pa - P0), 2.0)
        )
    return KD1, KD2, y


# Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162
# doi:10.1016/j.epsl.2018.10.0020012-821X
def FefO2_ONeill18_terms(PT, melt_wf, models=default_models):
    """
    Terms for O'Neill et al. (2018) equation (9a) to calculate fO2-Fe3+/FeT

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Dataframe of model options

    Returns
    -------
    str
        Terms a, B, FMQ
    """
    # 1n(Fe3Fe2) = a*DFMQ + B
    a = 0.25
    # mole fractions on a single cation basis in the melt based on oxide components (all
    # Fe as FeO) with no volatiles
    melt_comp = mg.melt_cation_proportion(melt_wf, "no", "no")
    B = (
        -1.36
        + 2.4 * melt_comp["Ca"]
        + 2.0 * melt_comp["Na"]
        + 3.7 * melt_comp["K"]
        - 2.4 * melt_comp["P"]
    )
    # FMQ
    T_K = PT["T"] + 273.15
    FMQ = 8.58 - (25050 / T_K)  # O'Neill (1987)
    return a, B, FMQ


def FefO2_Borisov18_terms(
    PT, melt_wf, models=default_models
):  # Eq. (4) from Borisov et al. (2018) CMP 173:98 doi:10.1007/s00410-018-1524-8
    T_K = PT["T"] + 273.15
    # Borisov et al. (2018) CMP 173
    a = 0.207
    # melt mole fraction with no volatiles and all Fe as FeOT
    melt_comp = mg.melt_mole_fraction(
        melt_wf, models, "no", "no", molmass="ONeill22", majors="ONeill22"
    )
    B = (
        4633.3 / T_K
        - 0.445 * melt_comp["SiO2"]
        - 0.900 * melt_comp["TiO2"]
        + 1.532 * melt_comp["MgO"]
        + 0.314 * melt_comp["CaO"]
        + 2.030 * melt_comp["Na2O"]
        + 3.355 * melt_comp["K2O"]
        - 4.851 * melt_comp["P2O5"]
        - 3.081 * melt_comp["SiO2"] * melt_comp["Al2O3"]
        - 4.370 * melt_comp["SiO2"] * melt_comp["MgO"]
        - 1.852
    )
    return a, B


def fO22Fe3FeT(fO2, PT, melt_wf, models=default_models):  # converting fO2 to Fe3/FeT
    """
    Fe3+/FeT in the melt from fO2

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "fO2" and column label of
        "option"

    Returns
    -------
    float
        Fe3+/FeT in the melt <class 'mpfr'>


    Model options for fO2
    -------------
    - 'Kress91A' [default] Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92
    doi:10.1007/BF00307328
    - 'Kress91' Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92
    doi:10.1007/BF00307328
    - 'ONeill18' Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162
    doi:10.1016/j.epsl.2018.10.0020012-821X
    - 'Borisov18' Eq. (4) from Borisov et al. (2018) CMP 173:98
    doi:10.1007/s00410-018-1524-8

    """
    model = models.loc["fO2", "option"]

    # Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    if model == "Kress91":
        a, PTterm = FefO2_KC91_Eq7_terms(PT, melt_wf, models)
        if models.loc["high precision", "option"] == "True":
            lnXFe2O3XFeO = a * gp.log(fO2) + PTterm
            XFe2O3XFeO = gp.exp(lnXFe2O3XFeO)
        else:
            lnXFe2O3XFeO = a * math.log(fO2) + PTterm
            XFe2O3XFeO = math.exp(lnXFe2O3XFeO)
        return (2.0 * XFe2O3XFeO) / ((2.0 * XFe2O3XFeO) + 1.0)

    # Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    elif model == "Kress91A":
        kd1, KD2, y = FefO2_KC91_EqA_terms(PT, melt_wf, models)
        XFeO15XFeO = (
            (kd1 * fO2**0.25)
            + (2.0 * y * KD2 * (kd1 ** (2.0 * y)) * (fO2 ** (0.5 * y)))
        ) / (1.0 + (1.0 - 2.0 * y) * KD2 * (kd1 ** (2.0 * y)) * (fO2 ** (0.5 * y)))
        return XFeO15XFeO / (XFeO15XFeO + 1.0)

    # Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162
    # doi:10.1016/j.epsl.2018.10.0020012-821X
    elif model == "ONeill18":
        a, B, FMQ = FefO2_ONeill18_terms(PT, melt_wf, models)
        if models.loc["high precision", "option"] == "True":
            DQFM = gp.log10(fO2) - FMQ
        else:
            DQFM = math.log10(fO2) - FMQ
        lnFe3Fe2 = a * DQFM + B
        Fe3Fe2 = 10.0 ** (lnFe3Fe2)
        return Fe3Fe2 / (Fe3Fe2 + 1.0)

    # Eq. (4) from Borisov et al. (2018) CMP 173:98 doi:10.1007/s00410-018-1524-8
    elif model == "Borisov18":
        a, B = FefO2_Borisov18_terms(PT, melt_wf, models)
        if models.loc["high precision", "option"] == "True":
            Fe3Fe2 = 10.0 ** (a * gp.log10(fO2) + B)
        else:
            Fe3Fe2 = 10.0 ** (a * math.log10(fO2) + B)
        return Fe3Fe2 / (Fe3Fe2 + 1.0)


def f_O2(PT, melt_wf, models=default_models):
    """
    fO2 from Fe3+/FeT in the melt

    Parameters
    ----------
    PT: dict
        Dictionary of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "fO2" and column label of
        "option"

    Returns
    -------
    float
        fO2 in bars <class 'mpfr'>


    Model options for fO2
    -------------
    - 'Kress91A' [default] Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92
    doi:10.1007/BF00307328
    - 'Kress91' Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92
    doi:10.1007/BF00307328
    - 'ONeill18' Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162
    doi:10.1016/j.epsl.2018.10.0020012-821X
    - 'Borisov18' Eq. (4) from Borisov et al. (2018) CMP 173:98
    doi:10.1007/s00410-018-1524-8

    """

    model = models.loc["fO2", "option"]

    def KC91(PT, melt_wf, models):
        a, PTterm = FefO2_KC91_Eq7_terms(PT, melt_wf, models)
        F = 0.5 * mg.Fe3Fe2(melt_wf)  # XFe2O3/XFeO
        alnfO2 = math.log(F) - PTterm
        fO2 = math.exp(alnfO2 / a)
        return fO2

    # if model == "yes":
    #    return 10.0**(setup.loc[run,"logfO2"])

    # Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    if model == "Kress91":
        fO2 = KC91(PT, melt_wf, models)
        return fO2

    # Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92 doi:10.1007/BF00307328
    elif model == "Kress91A":
        F = mg.Fe3Fe2(melt_wf)  # XFeO1.5/XFeO
        kd1, KD2, y = FefO2_KC91_EqA_terms(PT, melt_wf, models)

        def f(y, F, KD2, kd1, x):  # KC91A rearranged to equal 0
            f = (
                (2.0 * y - F + 2.0 * y * F) * KD2 * kd1 ** (2.0 * y) * x ** (0.5 * y)
                + kd1 * x**0.25
                - F
            )
            return f

        def df(y, F, KD2, kd1, x):  # derivative of above
            df = (0.5 * y) * (2.0 * y - F + 2.0 * y * F) * KD2 * kd1 ** (
                2.0 * y
            ) * x ** ((0.5 * y) - 1.0) + 0.25 * kd1 * x**-0.75
            return df

        def dx(x):
            diff = abs(0 - f(y, F, KD2, kd1, x))
            return diff

        def nr(x0, e1):
            delta1 = dx(x0)
            while delta1 > e1:
                x0 = x0 - f(y, F, KD2, kd1, x0) / df(y, F, KD2, kd1, x0)
                delta1 = dx(x0)
            return x0

        x0 = KC91(PT, melt_wf, models)

        fO2 = nr(x0, 1.0e-15)
        return fO2

    # Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162
    # doi:10.1016/j.epsl.2018.10.0020012-821X
    elif model == "ONeill18":
        F = mg.Fe3Fe2(melt_wf)  # Fe3+/Fe2+
        a, B, FMQ = FefO2_ONeill18_terms(PT, melt_wf, models)
        DQFM = (math.log10(F) - B) / a
        logfO2 = DQFM + FMQ
        return 10.0**logfO2

    # elif model == "S6ST": # remove?!?!
    #    S6T = melt_wf['S6ST']
    #    S62 = mg.overtotal2ratio(S6T)
    #    fO2 = mg.S6S2_2_fO2(S62,melt_wf,run,PT,setup,models)
    #    return fO2

    # Eq. (4) from Borisov et al. (2018) CMP 173:98 doi:10.1007/s00410-018-1524-8
    elif model == "Borisov18":
        F = mg.Fe3Fe2(melt_wf)
        a, B = FefO2_Borisov18_terms(PT, melt_wf, models)
        if models.loc["high precision", "option"] == "True":
            fO2 = 10.0 ** ((gp.log10(F) - B) / a)
        else:
            fO2 = 10.0 ** ((math.log10(F) - B) / a)
        return fO2


def S_Nash19_terms(PT):  # Nash et al. 2019
    # WORK IN PROGRESS
    T_K = PT["T"] + 273.15
    A = 8.0
    B = ((8.7436e6) / pow(T_K, 2.0)) - (27703.0 / T_K) + 20.273
    return A, B


def melt_density(PT, melt_wf, models=default_models):  # g/cm3
    """
    Melt density

    Parameters
    ----------
    PT: dict
        Dataframe of pressure-temperature conditions
        pressure (bars) as "P"
        temperature ('C) as "T"

    melt_wf: dict
        Dictionary of melt composition (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Minimum requirement is dataframe with index of "density" and column label of
        "option"

    Returns
    -------
    float
        melt density in g/cm3 <class 'mpfr'>


    Model options for density
    -------------
    - 'DensityX' [default] DensityX from Iacovino & Till (2019) Volcanica 2(1):1-10
    doi:10.30909/vol.02.01.0110
    - Only one option available currently, included for future development.

    """
    # DensityX from Iacovino & Till (2019) Volcanica 2(1):1-10
    # doi:10.30909/vol.02.01.0110
    if models.loc["density", "option"] == "DensityX":
        melt_comp = mg.melt_normalise_wf(melt_wf, "water", "yes")
        P = PT["P"]
        T = PT["T"]
        melt_dx = pd.DataFrame(
            [
                [
                    "sample",
                    melt_comp["SiO2"],
                    melt_comp["TiO2"],
                    melt_comp["Al2O3"],
                    melt_comp["FeO"],
                    melt_comp["Fe2O3"],
                    melt_comp["MgO"],
                    melt_comp["CaO"],
                    melt_comp["Na2O"],
                    melt_comp["K2O"],
                    melt_comp["H2O"],
                    P,
                    T,
                ]
            ]
        )
        melt_dx.columns = [
            "Sample_ID",
            "SiO2",
            "TiO2",
            "Al2O3",
            "FeO",
            "Fe2O3",
            "MgO",
            "CaO",
            "Na2O",
            "K2O",
            "H2O",
            "P",
            "T",
        ]
        output = dx.Density(melt_dx)
        density = output.loc[0, "density_g_per_cm"]
    return density


########################################################################################
# VOLUME ###############################################################################
########################################################################################


def vol_CORK(species, PT, models):
    """Appendix: Calculation of CORK volumes using eq. (4, A1) from Holland & Powell
    (1991) CMP 109:265-273 10.1007/BF00306484

    Args:
        species (str): gas species of interest (e.g., 'H2O', 'CO2')
        PT (dict): Dictionary of pressure-temperature conditions: pressure (bars) as
        "P", temperature ('C) as "T"
        models (pandas.DataFrame): Dataframe of the model options

    Returns:
        float: volume in kJ/kbar
    """

    P = PT["P"]
    T_K = PT["T"] + 273.15

    # Eq. (A.1)
    def MRK(P_kb, VMRK, R, T_K, a, b):  # MRK volume equation rearranged to equal 0
        return (
            P_kb * pow(VMRK, 3.0)
            - R * T_K * pow(VMRK, 2.0)
            - (b * R * T_K + pow(b, 2.0) * P_kb - a * pow(T_K, -0.5)) * VMRK
            - (a * b) * pow(T_K, -0.5)
        )

    def dMRK(P_kb, VMRK, R, T_K, a, b):  # derivative of above
        return (
            3.0 * P_kb * pow(VMRK, 2.0)
            - 2.0 * R * T_K * VMRK
            - (b * R * T_K + pow(b, 2.0) * P_kb - a * pow(T_K, -0.5))
        )

    def dVMRK(MRK, P_kb, VMRK, R, T_K, a, b):
        return abs(0 - MRK(P_kb, VMRK, R, T_K, a, b))

    def NR_VMRK(MRK, dMRK, VMRK0, e1, P_kb, R, T_K, a, b):
        delta1 = dVMRK(MRK, P_kb, VMRK0, R, T_K, a, b)
        while delta1 > e1:
            VMRK0 = VMRK0 - MRK(P_kb, VMRK0, R, T_K, a, b) / dMRK(
                P_kb, VMRK0, R, T_K, a, b
            )
            delta1 = dVMRK(MRK, P_kb, VMRK0, R, T_K, a, b)
        return VMRK0

    R = 8.3144598e-3  # in kJ/mol/K
    P_kb = P / 1000.0

    a, b, c, d, p0 = parameters_Holland91(species, PT, models)

    Vi = ((R * T_K) / P_kb) + b

    VMRK = NR_VMRK(MRK, dMRK, Vi, 1e-5, P_kb, R, T_K, a, b)

    if P_kb > p0:
        # Eq. (4)
        V = VMRK + c * pow((P_kb - p0), 0.5) + d * (P_kb - p0)
    else:
        V = VMRK

    return V


def vol_sCORK(species, PT, models):
    """Volume using eq. (7a) from Holland & Powell (1991) CMP 109:265-273
    10.1007/BF00306484

    Args:
        species (str): gas species of interest (e.g., 'H2O', 'CO2')
        PT (dict): Dictionary of pressure-temperature conditions: pressure (bars) as
        "P", temperature ('C) as "T"
        models (pandas.DataFrame): Dataframe of the model options


    Returns:
        float: volume in kJ/kbar
    """
    R = 8.3144598e-3
    P_bar = PT["P"]  # P in bar
    T = PT["T"] + 273.15  # T in K
    P = P_bar / 1000.0  # P in kb

    a, b, c, d, p0 = parameters_Holland91(species, PT, models)  # noqa

    # Eq. (7a)
    V = (
        ((R * T) / P)
        + b
        - ((a * R * T**0.5) / (((R * T) + (b * P)) * ((R * T) + (2.0 * b * P))))
        + (c * P**0.5)
        + (d * P)
    )

    return V


def gas_molar_volume(species, PT, models):
    """Calculates molar volume of a given gas species

    Args:
        species (str): gas species of interest (e.g., 'H2O', 'CO2')
        PT (dict): Dictionary of pressure-temperature conditions: pressure (bars) as
        "P", temperature ('C) as "T"
        models (pandas.DataFrame): Dataframe of the model options

    Returns:
        float: volume in kJ/kbar
    """
    if species == "H2O":
        model = models.loc["y_H2O", "option"]
        if model == "Holland91":
            V = vol_CORK(species, PT, models)
    if species == "CO2":
        model = models.loc["y_CO2", "option"]
        if model == "Holland91_eq4,A1-3_tab1":
            V = vol_CORK(species, PT, models)
        elif model == "Holland91_eq8,9_tab2":
            V = vol_sCORK(species, PT, models)
        elif model == "Holland91_eq8_tab1":
            V = vol_sCORK(species, PT, models)
    if species in "CH4":
        model = models.loc["y_CH4", "option"]
        if model == "Holland91_eq8,9_tab2":
            V = vol_sCORK(species, PT, models)
    if species == "H2":
        model = models.loc["y_H2", "option"]
        if model == "Holland91_eq8,9_tab2":
            V = vol_sCORK(species, PT, models)
    if species == "CO":
        model = models.loc["y_CO", "option"]
        if model == "Holland91_eq8,9_tab2":
            V = vol_sCORK(species, PT, models)

    return V


########################################################################################
# CONSTANTS ############################################################################
########################################################################################

species = [
    ["H", "", 1.008, 1.0, 0.0, 1.0, "", "", "", "", "", "", "", ""],
    ["C", "", 12.011, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["O", "", 15.999, -2.0, 0.0, "", "", "", "", "", 16.0, "", "", ""],
    ["Na", "", 22.99, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["Mg", "", 24.305, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["Al", "", 26.982, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["Si", "", 28.085, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["P", "", 30.974, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["S", "", 32.06, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["K", "", 39.098, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["Ca", "", 40.078, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["Ti", "", 47.867, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["Mn", "", 54.938, "", 0.0, "", "", "", "", "", "", "", "", ""],
    ["Fe", "", 55.845, "", 0.0, "", "", "", "", "", 55.85, "", 55.85, ""],
    ["SiO2", "Y", 60.083, 0.0, 2.0, 4.0, 1.0, 2.0, "", "", 60.0855, "Y", 60.08, "Y"],
    ["TiO2", "Y", 79.865, 0.0, 2.0, 3.0, 1.0, 2.0, "", "", 79.867, "Y", 79.9, "Y"],
    [
        "Al2O3",
        "Y",
        101.961,
        0.0,
        3.0,
        3.0,
        2.0,
        3.0,
        "",
        "",
        101.9633078,
        "Y",
        101.96,
        "Y",
    ],
    [
        "Fe2O3",
        "",
        159.687,
        0.0,
        3.0,
        3.0,
        2.0,
        3.0,
        "",
        "",
        143.77185,
        "Y",
        159.687,
        "",
    ],
    ["FeO1.5", "", 79.8435, 0.0, 1.5, 3.0, 1.0, 1.5, "", "", "", "", "", ""],
    ["FeO", "Y", 71.844, 0.0, 1.0, 2.0, 1.0, 1.0, "", "", 71.85, "Y", 71.85, "Y"],
    ["MnO", "Y", 70.937, 0.0, 1.0, 2.0, 1.0, 1.0, "", "", 70.94, "Y", 70.94, "Y"],
    ["MgO", "Y", 40.304, 0.0, 1.0, 2.0, 1.0, 1.0, "", "", 40.32, "Y", 40.32, "Y"],
    ["CaO", "Y", 56.077, 0.0, 1.0, 2.0, 1.0, 1.0, "", "", 56.06, "Y", 56.08, "Y"],
    ["Na2O", "Y", 61.979, 0.0, 1.0, 1.0, 2.0, 1.0, "", "", 61.88, "Y", 61.98, "Y"],
    ["K2O", "Y", 94.195, 0.0, 1.0, 1.0, 2.0, 1.0, "", "", 94.2, "Y", 94.2, "Y"],
    ["P2O5", "Y", 141.943, 0.0, 1.0, 5.0, 2.0, 1.0, "", "", 141.943, "N", 141.943, "N"],
    ["OH", "", 17.007, -1.0, 1.0, 1.0, 1.0, 1.0, "", "", "", "", "", ""],
    [
        "H2O",
        "",
        18.015,
        0.0,
        1.0,
        1.0,
        1.0,
        1.0,
        647.25,
        221.1925,
        18.015,
        "",
        18.015,
        "",
    ],
    ["H2S", "", 34.076, 0.0, 0.0, 1.0, 1.0, 1.0, 373.55, 90.0779, "", "", "", ""],
    ["CO", "", 28.01, 0.0, 1.0, 2.0, 1.0, 1.0, 133.15, 34.9571, "", "", "", ""],
    [
        "CO2",
        "",
        44.009,
        0.0,
        2.0,
        4.0,
        1.0,
        2.0,
        304.15,
        73.8659,
        44.009,
        "",
        44.009,
        "",
    ],
    ["CO3", "", 60.008, -2.0, 3.0, 4.0, 1.0, 3.0, "", "", "", "", "", ""],
    ["S2", "", 64.12, 0.0, 0.0, "", "", "", 208.15, 72.954, "", "", "", ""],
    ["SO2", "", 64.058, 0.0, 2.0, 4.0, 1.0, 2.0, 430.95, 78.7295, "", "", "", ""],
    ["SO3", "", 80.057, 0.0, 3.0, 6.0, 1.0, 3.0, "", "", "", "", "", ""],
    ["SO4", "", 96.056, -2.0, 4.0, 6, "", 4.0, "", "", "", "", "", ""],
    ["OCS", "", 60.07, 0.0, 1.0, "", "", "", 377.55, 65.8612, "", "", "", ""],
    ["O2", "", 31.998, 0.0, 2.0, "", "", "", 154.75, 50.7638, "", "", "", ""],
    ["H2", "", 2.016, 0.0, 0.0, "", "", "", 33.25, 12.9696, 2.016, "", 2.016, ""],
    ["CH4", "", 16.043, 0.0, 0.0, "", "", "", 191.05, 46.4069, "", "", "", ""],
    ["Ar", "", 39.948, "", "", "", "", "", "", "", "", "", "", ""],
    ["Ne", "", 20.1797, "", "", "", "", "", "", "", "", "", "", ""],
]
# If a paper doesn't have a molcular mass for a certain species, the molecular mass in
# "M" is used.
# Create the pandas DataFrame
species = pd.DataFrame(
    species,
    columns=[
        "species",
        "majors",
        "M",
        "overall_charge",
        "no_O",
        "cat_charge",
        "no_cat",
        "no_an",
        "Tcr",
        "Pcr",
        "M_Boulliung23",
        "majors_Boulliung23",
        "M_ONeill21",
        "majors_ONeill21",
    ],
)
species = species.set_index("species")

########################################################################################
########################################################################################
# IN DEVELOPMENT BELOW HERE ############################################################
########################################################################################
########################################################################################

########################################################################################
# ISOTOPE FRACTIONATION FACTORS ########################################################
########################################################################################

#############
# vapor #####S
#############


def beta_gas(PT, element, species, models):  # beta factors
    # from Richet et al. (1977)fitted to quadratic equation for 600 < T'C < 1300
    if models.loc["beta_factors", "option"] == "Richet77":
        t = 1.0 / (PT["T"] + 273.15)
        if element == "S":
            if species == "SO2":
                a, b, c = 4872.56428, 0.76400, 0.99975
            elif species == "S2":
                a, b, c = 1708.22425, -0.76202, 1.00031
            elif species == "OCS":
                a, b, c = 980.75175, 1.74954, 0.99930
            elif species == "H2S":
                a, b, c = 935.84901, 1.29355, 0.99969
        if element == "H":
            if species == "H2O":
                a, b, c = 635425.0, -61.78, 1.0197
            elif species == "H2":
                a, b, c = 150565.0, 211.22, 0.9355
            elif species == "CH4":
                a, b, c = 571003.0, -129.28, 1.0386
            elif species == "H2S":
                a, b, c = 327635, -9.8154, 1.0004
        if element == "C":
            if species == "CO2":
                a, b, c = 17637.0, 13.968, 0.9955
            elif species == "CO":
                a, b, c = 10080.0, 7.4637, 0.9978
            elif species == "CH4":
                a, b, c = 8527.5, 14.052, 0.9955
            elif species == "OCS":
                a, b, c = 14572.0, 8.359, 0.9973
        value = a * t**2 + b * t + c
    return value


######################
# Sulfur species #####
######################


def alpha_S_H2Sv_S2mm(PT, comp, models):  # alpha for 32/34S between H2S(v) and S2-(m)
    model = models.loc["alpha_H2S_S", "option"]
    if model == "Fiege15":  # Fiege et al. (2015) Chemical Geology eq. 8
        T_K = PT["T"] + 273.15
        lna103 = (10.84 * ((1000.0 / T_K) ** 2)) - 2.5
        if models.loc["high precision", "option"] == "True":
            a = gp.exp(lna103 / 1000.0)
        else:
            a = math.exp(lna103 / 1000.0)
    elif model == "no fractionation":
        a = 1.0
    return a


def alpha_S_SO2v_S6pm(PT, comp, models):  # alpha for 32/34S between SO2(v) and S6+(m)
    model = models.loc["alpha_SO2_SO4", "option"]
    if model == "Fiege15":  # Fiege et al. (2015) Chemical Geology eq. 9
        T_K = PT["T"] + 273.15
        lna103 = (
            (-0.42 * ((1000.0 / T_K) ** 3))
            - (2.133 * ((1000.0 / T_K) ** 3))
            - (0.105 * (1000.0 / T_K))
            - 0.41
        )
        if models.loc["high precision", "option"] == "True":
            a = gp.exp(lna103 / 1000.0)
        else:
            a = math.exp(lna103 / 1000.0)
    elif model == "no fractionation":
        a = 1.0
    return a


def alpha_S_H2Sv_H2Sm(PT, comp, models):  # alpha for 32/34S between H2S(v) and H2S(m)
    model = models.loc["alpha_S_H2Sv_H2Sm", "option"]
    if model == "no fractionation":  #
        a = 1.0
    return a


######################
# Carbon species #####
######################


def alpha_C_CO2v_CO32mm(
    PT, comp, models
):  # alpha for 13/12C between CO2(v) and CO32-(m)
    model = models.loc["alpha_C_CO2v_CO32mm", "option"]
    if model == "LeePP":
        a = math.exp(2.9 / 1000.0)
    elif model == "no fractionation":
        a = 1.0
    return a


def alpha_C_CO2v_CO2m(
    PT, comp, models
):  # alpha for 13/12C between CO2(v) and CO2mol(m)
    model = models.loc["alpha_C_CO2v_CO2m", "option"]
    if model == "Blank93":
        a = 1.0
    elif model == "no fractionation":
        a = 1.0
    return a


def alpha_C_CO2v_CO2Tm(PT, comp, models):  # alpha for 13/12C between CO2(v) and CO2T(m)
    model = models.loc["alpha_C_CO2v_CO2T", "option"]
    if model == "LeePP":
        a = math.exp(2.9 / 1000.0)  # FIX THIS
    elif model == "no fractionation":
        a = 1.0
    return a


def alpha_C_COv_COm(PT, comp, models):  # alpha for 13/12C between CO(v) and COmol(m)
    model = models.loc["alpha_C_COv_COm", "option"]
    if model == "no fractionation":
        a = 1.0
    return a


def alpha_C_CH4v_CH4m(
    PT, comp, models
):  # alpha for 13/12C between CH4(v) and CH4mol(m)
    model = models.loc["alpha_C_CH4v_CH4m", "option"]
    if model == "no fractionation":
        a = 1.0
    return a


######################
# Hydrogen species ###
######################


def alpha_H_H2Ov_H2Om(PT, comp, models):
    model = models.loc["alpha_H_H2Ov_H2Om", "option"]
    if model == "no fractionation":
        a = 1.0
    elif model == "Rust04":
        a = 0.9896
    return a


def alpha_H_H2Ov_OHmm(PT, comp, models):
    model = models.loc["alpha_H_H2Ov_OHmm", "option"]
    if model == "no fractionation":
        a = 1.0
    elif model == "Rust04":
        a = 1.0415
    return a


def alpha_H_H2v_H2m(PT, comp, models):  # alpha for D/H between H2(v) and H2(m)
    model = models.loc["alpha_H_H2v_H2m", "option"]
    if model == "no fractionation":
        a = 1.0
    return a


def alpha_H_CH4v_CH4m(PT, comp, models):  # alpha for D/H between CH4(v) and CH4(m)
    model = models.loc["alpha_H_CH4v_CH4m", "option"]
    if model == "no fractionation":
        a = 1.0
    return a


def alpha_H_H2Sv_H2Sm(PT, comp, models):  # alpha for D/H between H2S(v) and H2S(m)
    model = models.loc["alpha_H_H2Sv_H2Sm", "option"]
    if model == "no fractionation":  #
        a = 1.0
    return a
