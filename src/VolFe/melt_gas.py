# melt_gas.py

import pandas as pd
import math
import warnings as w

# import model_dependent_variables as mdv
import VolFe.model_dependent_variables as mdv
import VolFe.equilibrium_equations as eq

########################################################################################
# Contents #############################################################################
########################################################################################

# Fugacity, mole fraction, partial pressure
# Speciation
# Converting gas and melt compositions
# Volume and density

########################################################################################
# FUGACITY, MOLE FRACTION, PARTIAL PRESSURE ############################################
########################################################################################


def f_H2O(PT, melt_wf, models):
    """
    Calculates fugacity of water in the vapor from concentration (mole fraction) of
    water in the melt.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (e.g., SiO2, TiO2, etc.), particularly H2O.

    models: pandas.DataFrame
        Model options: minimum requirement is a DataFrame with index of “Hspeciation”
        and column label of “option”.

    Returns
    -------
    float
        fH2O in bars


    See C_H2O for model options for 'Hspeciation'.

    """
    Hspeciation = models.loc["Hspeciation", "option"]

    # fH2O = xmH2OT^2/CH2O
    if Hspeciation in ["none", "none+ideal", "none+regular"]:
        value = ((xm_H2OT_so(melt_wf)) ** 2.0) / mdv.C_H2O(PT, melt_wf, models)
        return value

    # fH2O = xmH2OT/CH2O
    elif Hspeciation == "linear":
        value = xm_H2OT_so(melt_wf) / mdv.C_H2O(PT, melt_wf, models)
        return value

    # regular or ideal: fH2O = xmH2Omol/CH2O
    else:
        value = xm_H2Omol_so(PT, melt_wf, models) / mdv.C_H2O(PT, melt_wf, models)
        return value


def f_CO2(PT, melt_wf, models):
    """
    Calculates fugacity of CO2 in the vapor from concentration (mole fraction) of total
    CO2 (i.e., CO2,mol + CO32-) in the melt.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.), in particular CO2.

    models: pandas.DataFrame
        Model options: minimum requirement is dataframe with index of “carbon dioxide”
        and column label of “option”.

    Returns
    -------
    float
        fCO2 in bars


    See C_CO2 for model options for "carbon dioxide".

    """
    CO3model = models.loc["carbon dioxide", "option"]

    # wtppm CO2 modified from Shishkina et al. (2014) Chem. Geol. 388:112-129
    if CO3model == "Shishkina14":
        f = (melt_wf["CO2"] * 1000000.0) / mdv.C_CO3(PT, melt_wf, models)

    # fCO2 = xmCO2/C_CO3
    else:
        f = xm_CO2_so(melt_wf) / mdv.C_CO3(PT, melt_wf, models)

    return f


def f_S2(PT, melt_wf, models):
    """
    Calculates fugacity of S2 in the vapor from concentration (weight fraction) of *S2-
    in the melt.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.), in particular *S2-.

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fS2 in bars

    """
    K = mdv.C_S(PT, melt_wf, models) / 1000000.0
    return ((melt_wf["S2-"] / K) ** 2.0) * mdv.f_O2(PT, melt_wf, models)


def f_H2(PT, melt_wf, models):
    """
    Calculates fugacity of H2 in the vapor from fH2O and fO2 based on melt composition.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fH2 in bars

    """
    K = mdv.KHOg(PT, models)
    return f_H2O(PT, melt_wf, models) / (K * mdv.f_O2(PT, melt_wf, models) ** 0.5)


def f_CO(PT, melt_wf, models):
    """
    Calculates fugacity of CO in the vapor from fCO2 and fO2 based on melt composition.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fCO in bars

    """
    K = mdv.KCOg(PT, models)
    return f_CO2(PT, melt_wf, models) / (K * mdv.f_O2(PT, melt_wf, models) ** 0.5)


def f_H2S(PT, melt_wf, models):
    """
    Calculates fugacity of H2S in the vapor from fS2, fH2O, and fO2 based on melt
    composition.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fH2S in bars

    """
    K = mdv.KHOSg(PT, models)
    return (
        K * f_S2(PT, melt_wf, models) ** 0.5 * f_H2O(PT, melt_wf, models)
    ) / mdv.f_O2(PT, melt_wf, models) ** 0.5


def f_SO2(PT, melt_wf, models):
    """
    Calculates fugacity of SO2 in the vapor from fS2 and fO2 based on melt composition.

    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fSO2 in bars

    """
    K = mdv.KOSg(PT, models)
    return K * mdv.f_O2(PT, melt_wf, models) * f_S2(PT, melt_wf, models) ** 0.5


def f_SO3(PT, melt_wf, models):
    """
    Calculates fugacity of SO3 in the vapor from fS2 and fO2 based on melt composition.

    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fSO2 in bars

    """
    K = mdv.KOSg2(PT, models)
    return (
        K * (mdv.f_O2(PT, melt_wf, models)) ** 1.5 * (f_S2(PT, melt_wf, models)) ** 0.5
    )


def f_CH4(PT, melt_wf, models):
    """
    Calculates fugacity of CH4 in the vapor from fCO2, fH2O, and fO2 based on melt
    composition.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
       Melt composition in weight fraction (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Models options.

    Returns
    -------
    float
        fCH4 in bars

    """
    K = mdv.KCOHg(PT, models)
    return (f_CO2(PT, melt_wf, models) * f_H2O(PT, melt_wf, models) ** 2.0) / (
        K * mdv.f_O2(PT, melt_wf, models) ** 2.0
    )


def f_OCS(PT, melt_wf, models):
    """
    Calculates fugacity of OCS in the vapor from fCO2, fH2O, and fH2S or fSO2 based on
    melt composition.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fOCS in bars

    """
    OCSmodel = models.loc["carbonylsulfide", "option"]
    K = mdv.KOCSg(PT, models)
    if OCSmodel == "COHS":
        if f_H2O(PT, melt_wf, models) > 0.0:
            return (f_CO2(PT, melt_wf, models) * f_H2S(PT, melt_wf, models)) / (
                f_H2O(PT, melt_wf, models) * K
            )
        else:
            return 0.0
    else:
        if f_CO2(PT, melt_wf, models) > 0.0:
            return ((f_CO(PT, melt_wf, models) ** 3.0) * f_SO2(PT, melt_wf, models)) / (
                (f_CO2(PT, melt_wf, models) ** 2.0) * K
            )
        else:
            return 0.0


def f_X(PT, melt_wf, models):
    """
    Calculates fugacity of X in the vapor from concentration (weight fraction) of X in
    the melt.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fX in bars

    """
    K = mdv.C_X(PT, melt_wf, models) / 1000000.0
    return melt_wf["XT"] / K


#######################
# oxygen fugacity #####
#######################


def fO22Dbuffer(PT, fO2, buffer, models):
    """
    Converts absolute fO2 to value relative-to-buffer.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    fO2: float
        fO2 in bars

    buffer: str
        Buffer of interest, currently available options are Ni-NiO (NNO) or
        Fayalite-Magnetite-Quartz (FMQ).

    models: pandas.DataFrame
        Models options.

    Returns
    -------
    float
        fO2 value relative to buffer

    """
    # Ni-NiO
    if buffer == "NNO":
        return math.log10(fO2) - mdv.NNO(PT, models)

    # Fayalite-Magnetite-Quartz
    elif buffer == "FMQ":
        return math.log10(fO2) - mdv.FMQ(PT, models)


def Dbuffer2fO2(PT, D, buffer, model):
    """
    Converts value relative-to-buffer to absolute fO2.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    D: float
        fO2 value relative to buffer.

    buffer: str
        Buffer of interest, current options are Ni-NiO (NNO) or
        Fayalite-Magnetite-Quartz (FMQ).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fO2 in bars

    """
    # Ni-NiO
    if buffer == "NNO":
        return 10.0 ** (D + mdv.NNO(PT, model))

    # Fayalite-Magnetite-Quartz
    elif buffer == "FMQ":
        return 10.0 ** (D + mdv.FMQ(PT, model))


def S6S2_2_fO2(S62, melt_wf, PT, models):
    """
    Converts S6+/S2- in the melt to fO2.


    Parameters
    ----------
    S62: float
        S6+/S2- in melt.

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        fO2 in bars

    """

    def calc_fO2(S62, melt_wf, PT, models):
        CSO4 = mdv.C_SO4(PT, melt_wf, models) / 1000000.0
        CS = mdv.C_S(PT, melt_wf, models) / 1000000.0
        if models.loc["H2S_m", "option"] == "False":
            fO2 = ((S62 * CS) / CSO4) ** 0.5
        elif models.loc["H2S_m", "option"] == "True":
            KHS = mdv.KHOSg(PT, models)
            CH2S = (
                (mdv.C_H2S(PT, melt_wf, models) / 1000000.0)
                / mdv.species.loc["H2S", "M"]
            ) * mdv.species.loc["S", "M"]
            CH2OT = mdv.C_H2O(PT, melt_wf, models)
            xmH2O = xm_H2OT_so(melt_wf)
            W = CSO4 / ((KHS * CH2S * (xmH2O**2.0 / CH2OT)) + CS)
            fO2 = (S62 / W) ** 0.5
        return fO2, xmH2O

    melt_wf["Fe3FeT"] = 0.1
    fO2, xmH2O = calc_fO2(S62, melt_wf, PT, models)
    Fe3FeT = mdv.fO22Fe3FeT(fO2, PT, melt_wf, models)
    delta1 = abs(melt_wf["Fe3FeT"] - Fe3FeT)
    conc = eq.melt_speciation(PT, melt_wf, models, 1.0, 1.0e-9)
    melt_wf["H2OT"] = conc["wm_H2O"]
    melt_wf["CO2"] = conc["wm_CO2"]
    delta2 = abs(xmH2O - xm_H2OT_so(melt_wf))
    while delta1 > 0.0001 or delta2 > 1.0e-9:
        melt_wf["Fe3FeT"] = Fe3FeT
        fO2, xmH2O = calc_fO2(S62, melt_wf, PT, models)
        Fe3FeT = mdv.fO22Fe3FeT(fO2, PT, melt_wf, models)
        delta1 = abs(melt_wf["Fe3FeT"] - Fe3FeT)
        conc = eq.melt_speciation(PT, melt_wf, models, 1.0, 1.0e-9)
        melt_wf["H2OT"] = conc["wm_H2O"]
        melt_wf["CO2"] = conc["wm_CO2"]
        delta2 = abs(xmH2O - xm_H2OT_so(melt_wf))
    else:
        fO2, xmH2O = calc_fO2(S62, melt_wf, PT, models)

    # S6p = melt_wf["ST"]*(S62/(1.+S62))
    # S2mT = melt_wf["ST"] - S6p
    # S_H2S_S2m = (mdv.KHOSg(PT,models)*(mdv.C_H2S(PT,melt_wf,models)/1000000.)*
    # (xmH2O**2./mdv.C_H2O(PT,melt_wf,models)))/(mdv.C_S(PT,melt_wf,models)/1000000.)
    # S_H2S = S2mT*(S_H2S_S2m/(1.+S_H2S_S2m))
    # H2S = (S_H2S/mdv.species.loc['S','M'])*mdv.species.loc['H2S','M']
    # S2m = S2mT - S_H2S
    return fO2


########################
# partial pressure #####
########################


def p_H2(PT, melt_wf, models):
    """
    Converts fH2 based on melt composition to partial pressure of H2 (pH2).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pH2 in bars

    """
    return f_H2(PT, melt_wf, models) / mdv.y_H2(PT, models)


def p_H2O(PT, melt_wf, models):
    """
    Converts fH2O based on melt composition to partial pressure of H2O (pH2O).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pH2O in bars

    """
    return f_H2O(PT, melt_wf, models) / mdv.y_H2O(PT, models)


def p_O2(PT, melt_wf, models):
    """
    Converts fO2 based on melt composition to partial pressure of O2 (pO2).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pO2 in bars

    """
    return mdv.f_O2(PT, melt_wf, models) / mdv.y_O2(PT, models)


def p_SO2(PT, melt_wf, models):
    """
    Converts fSO2 based on melt composition to partial pressure of SO2 (pSO2).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pSO2 in bars

    """
    return f_SO2(PT, melt_wf, models) / mdv.y_SO2(PT, models)


def p_S2(PT, melt_wf, models):
    """
    Converts fS2 based on melt composition to partial pressure of S2 (pS2).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.)

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pS2 in bars

    """
    return f_S2(PT, melt_wf, models) / mdv.y_S2(PT, models)


def p_H2S(PT, melt_wf, models):
    """
    Converts fH2S based on melt composition to partial pressure of H2S (pH2S).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pH2S in bars

    """
    return f_H2S(PT, melt_wf, models) / mdv.y_H2S(PT, models)


def p_CO2(PT, melt_wf, models):
    """
    Converts fCO2 based on melt composition to partial pressure of CO2 (pCO2).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pCO2 in bars

    """
    return f_CO2(PT, melt_wf, models) / mdv.y_CO2(PT, models)


def p_CO(PT, melt_wf, models):
    """
    Converts fCO based on melt composition to partial pressure of CO (pCO).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pCO in bars

    """
    return f_CO(PT, melt_wf, models) / mdv.y_CO(PT, models)


def p_CH4(PT, melt_wf, models):
    """
    Converts fCH4 based on melt composition to partial pressure of CH4 (pCH4).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pCH4 in bars

    """
    return f_CH4(PT, melt_wf, models) / mdv.y_CH4(PT, models)


def p_OCS(PT, melt_wf, models):
    """
    Converts fOCS based on melt composition to partial pressure of OCS (pOCS).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pOCS in bars

    """
    return f_OCS(PT, melt_wf, models) / mdv.y_OCS(PT, models)


def p_X(PT, melt_wf, models):
    """
    Converts fX based on melt composition to partial pressure of X (pX).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        pX in bars

    """
    return f_X(PT, melt_wf, models) / mdv.y_X(PT, models)


def p_tot(PT, melt_wf, models):
    """
    Calculate total pressure by summing partial pressures of all vapor species based on
    melt composition.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        total P in bars

    """
    return (
        p_H2(PT, melt_wf, models)
        + p_H2O(PT, melt_wf, models)
        + p_O2(PT, melt_wf, models)
        + p_SO2(PT, melt_wf, models)
        + p_S2(PT, melt_wf, models)
        + p_H2S(PT, melt_wf, models)
        + p_CO2(PT, melt_wf, models)
        + p_CO(PT, melt_wf, models)
        + p_CH4(PT, melt_wf, models)
        + p_OCS(PT, melt_wf, models)
        + p_X(PT, melt_wf, models)
    )


############################
# vapor molar fraction #####
############################


def xg_H2(PT, melt_wf, models):
    """
    Converts pH2 based on melt composition to mole fraction of H2 in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of H2 in the vapor

    """
    return p_H2(PT, melt_wf, models) / PT["P"]


def xg_H2O(PT, melt_wf, models):
    """
    Converts pH2O based on melt composition to mole fraction of H2O in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of H2O in the vapor

    """
    return p_H2O(PT, melt_wf, models) / PT["P"]


def xg_O2(PT, melt_wf, models):
    """
    Converts pO2 based on melt composition to mole fraction of O2 in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of O2 in the vapor

    """
    return p_O2(PT, melt_wf, models) / PT["P"]


def xg_SO2(PT, melt_wf, models):
    """
    Converts pSO2 based on melt composition to mole fraction of SO2 in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of SO2 in the vapor

    """
    return p_SO2(PT, melt_wf, models) / PT["P"]


def xg_S2(PT, melt_wf, models):
    """
    Converts pS2 based on melt composition to mole fraction of S2 in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of S2 in the vapor

    """
    return p_S2(PT, melt_wf, models) / PT["P"]


def xg_H2S(PT, melt_wf, models):
    """
    Converts pH2S based on melt composition to mole fraction of H2S in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of H2S in the vapor

    """
    return p_H2S(PT, melt_wf, models) / PT["P"]


def xg_CO2(PT, melt_wf, models):
    """
    Converts pCO2 based on melt composition to mole fraction of CO2 in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of CO2 in the vapor

    """
    return p_CO2(PT, melt_wf, models) / PT["P"]


def xg_CO(PT, melt_wf, models):
    """
    Converts pCO based on melt composition to mole fraction of CO in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of CO in the vapor

    """
    return p_CO(PT, melt_wf, models) / PT["P"]


def xg_CH4(PT, melt_wf, models):
    """
    Converts pCH4 based on melt composition to mole fraction of CH4 in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of CH4 in the vapor

    """
    return p_CH4(PT, melt_wf, models) / PT["P"]


def xg_OCS(PT, melt_wf, models):
    """
    Converts pOCS based on melt composition to mole fraction of OCS in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of OCS in the vapor

    """
    return p_OCS(PT, melt_wf, models) / PT["P"]


def xg_X(PT, melt_wf, models):
    """
    Converts pX based on melt composition to mole fraction of X in the vapor.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of X in the vapor

    """
    return p_X(PT, melt_wf, models) / PT["P"]


def Xg_tot(PT, melt_wf, models):
    """
    Sum over over all vapor species of (mole fraction)*(molecular mass)


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        Sum over over all vapor species of (mole fraction)*(molecular mass)

    """
    species_X = models.loc["species X", "option"]
    Xg_t = (
        xg_CO2(PT, melt_wf, models) * mdv.species.loc["CO2", "M"]
        + xg_CO(PT, melt_wf, models) * mdv.species.loc["CO", "M"]
        + xg_O2(PT, melt_wf, models) * mdv.species.loc["O2", "M"]
        + xg_H2O(PT, melt_wf, models) * mdv.species.loc["H2O", "M"]
        + xg_H2(PT, melt_wf, models) * mdv.species.loc["H2", "M"]
        + xg_CH4(PT, melt_wf, models) * mdv.species.loc["CH4", "M"]
        + xg_SO2(PT, melt_wf, models) * mdv.species.loc["SO2", "M"]
        + xg_S2(PT, melt_wf, models) * mdv.species.loc["S2", "M"]
        + xg_H2S(PT, melt_wf, models) * mdv.species.loc["H2S", "M"]
        + xg_OCS(PT, melt_wf, models) * mdv.species.loc["OCS", "M"]
        + xg_X(PT, melt_wf, models) * mdv.species.loc[species_X, "M"]
    )
    return Xg_t


##########################
# melt mole fraction #####
##########################


def wm_vol(melt_wf):
    """
    Calculate wt% of H2OT-eq + CO2-eq in the melt.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        wt% of [H2OT-eq + CO2-eq] in the melt

    """
    wm_H2OT = 100.0 * melt_wf["H2OT"]
    wm_CO2 = 100.0 * melt_wf["CO2"]
    return wm_H2OT + wm_CO2  # + wm_S(wm_ST) + wm_SO3(wm_ST)


def wm_nvol(melt_wf):
    """
    Calculate wt% of everything except H2OT-eq + CO2-eq in the melt (excluding S and X).


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        wt% of everything except [H2OT-eq + CO2-eq] in the melt

    """
    return 100.0 - wm_vol(melt_wf)


def M_m_SO(melt_wf):
    """
    Calculate molecular mass of the melt on a singular oxygen basis with no volatiles
    and all Fe as FeOT.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        molecular mass of melt on singular oxygen basis

    """
    # single oxide, no volatiles, all Fe as FeOT
    melt_comp = melt_single_O(melt_wf, "no", "no")
    M_m = 1.0 / melt_comp["Xmtot"]
    return M_m


def M_m_ox(melt_wf, models):  # no volatiles, all Fe as FeOT
    """
    Calculates molecular mass of the melt on an oxide basis with no volatiles and all Fe
    as FeOT.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        molecular mass of melt on oxide basis

    """
    melt_comp = melt_mole_fraction(melt_wf, models, "no", "no")
    M_m = 1.0 / melt_comp["mol_tot"]
    return M_m


def Xm_H2OT(melt_wf):
    """
    Calculates number of moles of H2OT in the melt.


    Parameters
    ----------
    melt_wf: dict
        H2OT mass concentration in weight fraction

    Returns
    -------
    float
        number of moles of H2OT in melt

    """
    wm_H2OT = 100.0 * melt_wf["H2OT"]
    return wm_H2OT / mdv.species.loc["H2O", "M"]


def Xm_CO2(melt_wf):
    """
    Calculates number of moles of CO2 in the melt from mass of CO2.


    Parameters
    ----------
    melt_wf: dict
        CO2 mass concentration in weight fraction

    Returns
    -------
    float
        number of moles of CO2 in melt

    """
    wm_CO2 = 100.0 * melt_wf["CO2"]
    return wm_CO2 / mdv.species.loc["CO2", "M"]


# Mole fraction in the melt based on mixing between volatile-free melt on a singular
# oxygen basis and volatiles
def Xm_m_so(melt_wf):  # singular oxygen basis
    """
    Calculates number of moles of melt in the melt on a singular oxygen basis.


    Parameters
    ----------
    melt_wf: dictionary
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        number of moles of melt in melt

    """
    return wm_nvol(melt_wf) / M_m_SO(melt_wf)


def Xm_tot_so(melt_wf):
    """
    Calculates total moles in the melt on a singular oxygen basis including H2OT, CO2,
    and melt.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        total moles of melt

    """
    return Xm_H2OT(melt_wf) + Xm_CO2(melt_wf) + Xm_m_so(melt_wf)


def xm_H2OT_so(melt_wf):
    """
    Calculates mole fraction of H2OT in the melt on a singular oxygen basis including
    H2OT, CO2, and melt.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        mole fraction of H2OT in the melt

    """
    return Xm_H2OT(melt_wf) / Xm_tot_so(melt_wf)


def xm_CO2_so(melt_wf):
    """
    Calculates mole fraction of CO2 in the melt on a singular oxygen basis including
    H2OT, CO2, and melt.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        mole fraction of CO2 in the melt

    """
    return Xm_CO2(melt_wf) / Xm_tot_so(melt_wf)


def xm_melt_so(melt_wf):
    """
    Calculates mole fraction of melt in the melt on a singular oxygen basis including
    H2OT, CO2, and melt.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        mole fraction of melt in the melt

    """
    return Xm_m_so(melt_wf) / Xm_tot_so(melt_wf)


def Xm_t_so(melt_wf):
    """
    Calculates total "mass" in the melt on a singular oxygen basis including H2OT, CO2,
    and melt from sum of (mole fraction)*(molecular mass).


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        total "mass" in the melt

    """
    return (
        xm_H2OT_so(melt_wf) * mdv.species.loc["H2O", "M"]
        + xm_CO2_so(melt_wf) * mdv.species.loc["CO2", "M"]
        + xm_melt_so(melt_wf) * M_m_SO(melt_wf)
    )


########################################################################################
# SPECIATION ###########################################################################
########################################################################################


def ratio2overtotal(x):
    """
    Converts ratio of two species (e.g., Fe3+/Fe2+) to ratio of numerator species over
    total (e.g., Fe3+/FeT).


    Parameters
    ----------
    x: float
        Ratio of species 1 over species 2 (i.e., species1/species2).

    Returns
    -------
    float
        ratio of species 1 over total of both species (i.e.,
        species1/[species1+species2)

    """
    return x / (x + 1.0)


def overtotal2ratio(x):
    """
    Converts ratio of numerator species over total (e.g., Fe3+/FeT) to ratio of two
    species (e.g., Fe3+/Fe2+).


    Parameters
    ----------
    x: float
        ratio of species 1 over total of both species (i.e.,
        species1/[species1+species2])

    Returns
    -------
    float
        ratio of species 1 over species 2 (i.e., species1/species2)

    """
    return x / (1.0 - x)


###########################
# hydrogen speciation #####
###########################


def xm_OH_so(PT, melt_wf, models):
    """
    Calculate mole fraction of OH- in the melt on a singular oxygen basis.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of OH- in the melt


    See C_H2O for model options for 'Hspeciation'.

    """
    Hspeciation = models.loc["Hspeciation", "option"]

    Z = xm_H2OT_so(melt_wf)

    def f(
        A, B, C, Z, x
    ):  # regular solution model rearranged to equal 0 to solve for xm_OH
        return (
            A
            + B * x
            + C * (Z - 0.5 * x)
            + math.log((x**2.0) / ((Z - 0.5 * x) * (1.0 - Z - 0.5 * x)))
        )

    def df(B, C, Z, x):  # derivative of above
        return (
            B
            - 0.5 * C
            - (2.0 / x)
            + (0.5 / (Z - 0.5 * x) + (0.5 / (1.0 - Z - 0.5 * x)))
        )

    def dx(A, B, C, x):
        return abs(0 - f(A, B, C, Z, x))

    def nr(A, B, C, x0, e1):
        delta1 = dx(A, B, C, x0)
        while delta1 > e1:
            x0 = x0 - f(A, B, C, Z, x0) / df(B, C, Z, x0)
            delta1 = dx(A, B, C, x0)
        return x0

    if Z == 0.0:
        return 0.0

    # ideal mixture from Silver & Stolper (1985) J. Geol. 93(2) 161-177
    elif Hspeciation in [
        "ideal",
        "none+ideal",
    ]:
        K = mdv.KHOm(PT, melt_wf, models)
        return (0.5 - (0.25 - ((K - 4.0) / K) * (Z - Z**2.0)) ** 0.5) / (
            (K - 4.0) / (2.0 * K)
        )

    # all OH-
    elif Hspeciation == "none":
        return 0.0

    # Regular solution model from Silver & Stolper (1989) J.Pet. 30(3)667-709
    elif Hspeciation in [
        "regular",
        "none+regular",
    ]:
        tol = 0.000001  # tolerance
        K = mdv.KHOm(PT, melt_wf, models)
        x0 = (0.5 - (0.25 - ((K - 4.0) / K) * (Z - Z**2.0)) ** 0.5) / (
            (K - 4.0) / (2.0 * K)
        )
        A, B, C = mdv.KregH2O(PT, melt_wf, models)
        return nr(A, B, C, x0, tol)


def xm_H2Omol_so(PT, melt_wf, models):
    """
    Calculate mole fraction of H2Omol in the melt on a singular oxygen basis.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        mole fraction of H2Omol in the melt

    """
    Z = xm_H2OT_so(melt_wf)
    return Z - 0.5 * xm_OH_so(PT, melt_wf, models)


def wm_H2Omol_OH(PT, melt_wf, models):  # wt fraction
    """
    Converts mole fractions of OH- and H2Omol on a singular oxygen basis to weight
    fractions of OH- and H2Omol.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    tuple(float,float)
        weight fractions in the melt of H2Omol, OH-

    """
    H2Omol = xm_H2Omol_so(PT, melt_wf, models) * mdv.species.loc["H2O", "M"]
    OH_H2O = xm_OH_so(PT, melt_wf, models) * 0.5 * mdv.species.loc["H2O", "M"]
    melt = xm_melt_so(melt_wf) * M_m_SO(melt_wf)
    CO2T = xm_CO2_so(melt_wf) * mdv.species.loc["CO2", "M"]
    wm_H2Omol = H2Omol / (H2Omol + OH_H2O + melt + CO2T)
    wm_OH_H2O = OH_H2O / (H2Omol + OH_H2O + melt + CO2T)
    wm_OH = ((wm_OH_H2O / mdv.species.loc["H2O", "M"]) * 2.0) * mdv.species.loc[
        "OH", "M"
    ]
    return wm_H2Omol, wm_OH


#########################
# carbon speciation #####
#########################


def xm_CO32_CO2mol(PT, melt_wf, models):  # mole fraction
    """
    Calculate mole fraction of carbonate and CO2,mol in the melt on a singular oxygen
    basis.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    tuple(float,float)
        mole fraction in the melt of carbonate, CO2mol


    See KCOm for model options for 'Cspeccomp'.

    """
    xm_CO2T = xm_CO2_so(melt_wf)
    xm_H2OT = xm_H2OT_so(melt_wf)
    if models.loc["Cspeccomp", "option"] == "Basalt":
        xm_CO32 = xm_CO2T
        xm_CO2mol = 0.0
    elif models.loc["Cspeccomp", "option"] == "Rhyolite":
        xm_CO32 = 0.0
        xm_CO2mol = xm_CO2T
    else:
        K = mdv.KCOm(PT, melt_wf, models)
        xm_CO32 = (K * xm_CO2T * (1.0 - xm_CO2T - xm_H2OT)) / (
            1.0 + K * (1.0 - xm_CO2T - xm_H2OT)
        )
        xm_CO2mol = xm_CO2T - xm_CO32
    return xm_CO32, xm_CO2mol


def wm_CO32_CO2mol(PT, melt_wf, models):  # wt fraction
    """
    Converts mole fractions of carbonate and CO2mol on a singular oxygen basis to weight
    fractions of carbonate and CO2mol.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    tuple(float,float)
        weight fractions in the melt of carbonate, CO2mol

    """
    xm_CO32, xm_CO2mol = xm_CO32_CO2mol(PT, melt_wf, models)
    CO2mol = xm_CO2mol * mdv.species.loc["CO2", "M"]
    CO2carb = xm_CO32 * mdv.species.loc["CO2", "M"]
    melt = xm_melt_so(melt_wf) * M_m_SO(melt_wf)
    H2OT = xm_H2OT_so(melt_wf) * mdv.species.loc["H2O", "M"]
    wm_CO2mol = CO2mol / (CO2mol + CO2carb + H2OT + melt)
    wm_CO2carb = CO2carb / (CO2mol + CO2carb + H2OT + melt)
    return wm_CO2carb, wm_CO2mol


##########################
# sulfur speciation ######
##########################


def S6S2(PT, melt_wf, models):
    """
    Calculates S6+/S2- in the melt (where S2- includes *S2- and H2S).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        S6+/S2- in the melt

    """
    model = models.loc["sulfate", "option"]
    if model == "Nash19":
        A, B = mdv.S_Nash19_terms(PT)
        result = 10.0 ** (A * math.log10(Fe3Fe2(melt_wf)) + B)
    else:
        CSO4 = mdv.C_SO4(PT, melt_wf, models)
        CS = mdv.C_S(PT, melt_wf, models)
        fO2 = mdv.f_O2(PT, melt_wf, models)
        if models.loc["H2S_m", "option"] == "False":
            result = (CSO4 / CS) * fO2**2.0
        elif models.loc["H2S_m", "option"] == "True":
            KHS = mdv.KHOSg(PT, models)
            CH2S = mdv.C_H2S(PT, melt_wf, models)
            CH2OT = mdv.C_H2O(PT, melt_wf, models)
            xmH2O = xm_H2OT_so(melt_wf)
            result = (CSO4 / ((KHS * CH2S * (xmH2O**2.0 / CH2OT)) + CS)) * fO2**2.0
    return result


def S6ST(PT, melt_wf, models):
    """
    Calculates S6+/ST in the melt (where ST includes all sulfur species in the melt:
    *S2-, H2S, and S6+).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        S6+/ST in the melt

    """
    S6S2_ = S6S2(PT, melt_wf, models)
    return S6S2_ / (S6S2_ + 1.0)


def wm_S(PT, melt_wf, models):
    """
    Calculates wt% total S2- as S in the melt (both *S2- and H2S).


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        total S2- as S in wt%

    """
    wm_ST = 100.0 * melt_wf["ST"]
    S6ST_ = S6ST(PT, melt_wf, models)
    return wm_ST * (1.0 - S6ST_)


def wm_SO3(PT, melt_wf, models):
    """
    Calculates wt% total S6+ in the melt as SO3.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.

    Returns
    -------
    float
        S6+ as SO3 in wt%

    """
    wm_ST = 100.0 * melt_wf["ST"]
    S6ST_ = S6ST(PT, melt_wf, models)
    return ((wm_ST * S6ST_) / mdv.species.loc["S", "M"]) * mdv.species.loc["SO3", "M"]


#######################
# iron speciation #####
#######################


def Fe3Fe2(melt_wf):
    """
    Converts Fe3+/FeT to Fe3+/Fe2+ in the melt.


    Parameters
    ----------
    melt_wf: dict
        Fe3+/FeT as "Fe3FeT"

    Returns
    -------
    float
        Fe3+/Fe2+ in the melt

    """
    Fe3FeT = melt_wf["Fe3FeT"]
    return Fe3FeT / (1.0 - Fe3FeT)


def Wm_FeT(melt_wf, molmass="M"):
    """
    Calculate weight of total Fe in the melt from FeOT or Fe2O3T or FeO+Fe2O3.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        weight total Fe in the melt

    """
    if melt_wf["FeOT"] > 0.0:
        return (melt_wf["FeOT"] / mdv.species.loc["FeO", molmass]) * mdv.species.loc[
            "Fe", molmass
        ]
    elif melt_wf["Fe2O3T"] > 0.0:
        return (
            (melt_wf["Fe2O3T"] / mdv.species.loc["Fe2O3", molmass])
            * 2.0
            * mdv.species.loc["Fe", molmass]
        )
    else:
        return (
            (
                (melt_wf["FeO"] / mdv.species.loc["FeO", molmass])
                + (melt_wf["Fe2O3"] / mdv.species.loc["Fe2O3", molmass])
            )
            * 2.0
            * mdv.species.loc["Fe", molmass]
        )


def Wm_FeO(melt_wf, molmass="M"):
    """
    Calculate weight of FeO in the melt from total Fe and Fe3+/FeT.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        weight FeO in the melt

    """
    Fe3FeT = melt_wf["Fe3FeT"]
    return (
        (Wm_FeT(melt_wf, molmass=molmass) / mdv.species.loc["Fe", molmass])
        * (1.0 - Fe3FeT)
        * mdv.species.loc["FeO", molmass]
    )


def Wm_Fe2O3(melt_wf, molmass="M"):
    """
    Calculate weight of Fe2O3 in the melt from total Fe and Fe3+/FeT.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        weight Fe2O3 in the melt

    """
    Fe3FeT = melt_wf["Fe3FeT"]
    return (
        ((Wm_FeT(melt_wf, molmass=molmass) / mdv.species.loc["Fe", molmass]) * Fe3FeT)
        / 2.0
    ) * mdv.species.loc["Fe2O3", molmass]


def Wm_FeOT(melt_wf, molmass="M"):
    """
    Calculate weight of FeOT in the melt from total Fe.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        weight FeOT in the melt

    """
    return (
        Wm_FeT(melt_wf, molmass=molmass) / mdv.species.loc["Fe", molmass]
    ) * mdv.species.loc["FeO", molmass]


def wm_Fe_nv(melt_wf):  # no volatiles
    """
    Calculate weight fraction of total Fe in the melt from all oxides, FeO, Fe2O3 but no
    volatiles.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    Returns
    -------
    float
        weight fraction of total Fe in the melt

    """
    Wm_tot = (
        melt_wf["SiO2"]
        + melt_wf["TiO2"]
        + melt_wf["Al2O3"]
        + melt_wf["MnO"]
        + melt_wf["MgO"]
        + melt_wf["MnO"]
        + melt_wf["CaO"]
        + melt_wf["Na2O"]
        + melt_wf["K2O"]
        + melt_wf["P2O5"]
        + Wm_FeO(melt_wf)
        + Wm_Fe2O3(melt_wf)
    )
    FeT = mdv.species.loc["Fe", "M"] * (
        (2.0 * Wm_Fe2O3(melt_wf) / mdv.species.loc["Fe2O3", "M"])
        + (Wm_FeO(melt_wf) / mdv.species.loc["FeO", "M"])
    )
    return 100.0 * FeT / Wm_tot


def Fe3FeT_i(PT, melt_wf, models):
    """
    Calculate initial Fe3+/FeT in the melt based on either: Fe3+/FeT, logfO2, NNO, FMQ,
    S6+/ST, or FeO+Fe2O3.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.


    Returns
    -------
    float
        Fe3+/FeT in the melt

    """
    model = models.loc["fO2", "option"]

    if model == "buffered":
        fO2 = 10 ** (melt_wf["logfO2_i"])
        return mdv.fO22Fe3FeT(fO2, PT, melt_wf, models)
    else:
        if pd.isnull(melt_wf["Fe3FeT_i"]) is False:
            if (
                pd.isnull(melt_wf["logfO2_i"]) is False
                or pd.isnull(melt_wf["DNNO"]) is False
                or pd.isnull(melt_wf["DFMQ"]) is False
                or pd.isnull(melt_wf["S6ST_i"]) is False
                or pd.isnull(melt_wf["Fe2O3"]) is False
                or pd.isnull(melt_wf["FeO"]) is False
            ):
                w.warn(
                    "you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered Fe3+/FeT"  # noqa E501
                )
            return melt_wf["Fe3FeT_i"]
        elif pd.isnull(melt_wf["logfO2_i"]) is False:
            if (
                pd.isnull(melt_wf["DNNO"]) is False
                or pd.isnull(melt_wf["DFMQ"]) is False
                or pd.isnull(melt_wf["S6ST_i"]) is False
                or pd.isnull(melt_wf["Fe2O3"]) is False
                or pd.isnull(melt_wf["FeO"]) is False
            ):
                w.warn(
                    "you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered log(fO2)"  # noqa E501
                )
            fO2 = 10.0 ** (melt_wf["logfO2_i"])
            return mdv.fO22Fe3FeT(fO2, PT, melt_wf, models)
        elif pd.isnull(melt_wf["DNNO"]) is False:
            if (
                pd.isnull(melt_wf["DFMQ"]) is False
                or pd.isnull(melt_wf["S6ST_i"]) is False
                or pd.isnull(melt_wf["Fe2O3"]) is False
                or pd.isnull(melt_wf["FeO"]) is False
            ):
                w.warn(
                    "you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered DNNO"  # noqa E501
                )
            D = melt_wf["DNNO"]
            fO2 = Dbuffer2fO2(PT, D, "NNO", models)
            return mdv.fO22Fe3FeT(fO2, PT, melt_wf, models)
        elif pd.isnull(melt_wf["DFMQ"]) is False:
            if (
                pd.isnull(melt_wf["S6ST_i"]) is False
                or pd.isnull(melt_wf["Fe2O3"]) is False
                or pd.isnull(melt_wf["FeO"]) is False
            ):
                w.warn(
                    "you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered DFMQ"  # noqa E501
                )
            D = melt_wf["DFMQ"]
            fO2 = Dbuffer2fO2(PT, D, "FMQ", models)
            return mdv.fO22Fe3FeT(fO2, PT, melt_wf, models)
        elif pd.isnull(melt_wf["S6ST_i"]) is False:
            if (
                pd.isnull(melt_wf["Fe2O3"]) is False
                or pd.isnull(melt_wf["FeO"]) is False
            ):
                w.warn(
                    "you entered more than one way to infer iron speciation, note that this calcualtion is only considering the entered S6+/ST"  # noqa E501
                )
            S6T = melt_wf["S6ST_i"]
            S62 = overtotal2ratio(S6T)
            fO2 = S6S2_2_fO2(S62, melt_wf, PT, models)
            return mdv.fO22Fe3FeT(fO2, PT, melt_wf, models)
        else:
            return ((2.0 * melt_wf["Fe2O3"]) / mdv.species.loc["Fe2O3", "M"]) / (
                ((2.0 * melt_wf["Fe2O3"]) / mdv.species.loc["Fe2O3", "M"])
                + (melt_wf["FeO"] / mdv.species.loc["FeO", "M"])
            )


########################################################################################
# converting gas and melt compositions #################################################
########################################################################################


# C/S ratio of the vapor
def gas_CS(PT, melt_wf, models):
    """
    Calculate C/S mole fraction ratio in the vapor including all C and S species based
    on the melt composition.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.


    Returns
    -------
    float
        C/S mole fraction ratio in the vapor

    """
    S_values = (
        (2.0 * xg_S2(PT, melt_wf, models))
        + xg_SO2(PT, melt_wf, models)
        + xg_OCS(PT, melt_wf, models)
        + xg_H2S(PT, melt_wf, models)
    )
    if S_values > 0:
        xgCS = (
            xg_CO(PT, melt_wf, models)
            + xg_CO2(PT, melt_wf, models)
            + xg_OCS(PT, melt_wf, models)
            + xg_CH4(PT, melt_wf, models)
        ) / S_values
    else:
        xgCS = ""
    return xgCS


def gas_CS_alt(xg):
    """
    Calculate C/S mole fraction ratio in the vapor including all C and S species based
    on gas composition.


    Parameters
    ----------
    xg: dict
        Mole fraction of vapor species.


    Returns
    -------
    float
        C/S mole fraction ratio in the vapor

    """
    S_values = (2.0 * xg["S2"]) + xg["SO2"] + xg["OCS"] + xg["H2S"]
    if S_values > 0:
        xgCS = (xg["CO"] + xg["CO2"] + xg["OCS"] + xg["CH4"]) / S_values
    else:
        xgCS = ""
    return xgCS


# all carbon as CO2 and all hydrogen as H2O
def melt_H2O_CO2_eq(melt_wf):
    """
    Calculate weight fraction of all C as CO2 (CO2-eq) and all H as H2O (H2O-eq) in the
    melt.


    Parameters
    ----------
    PT: dict
        Pressure (bars) as "P" and temperature ('C) as "T".

    melt_wf: dict
        Melt composition in weight fraction (SiO2, TiO2, etc.).

    models: pandas.DataFrame
        Model options.


    Returns
    -------
    tuple(float,float)
        weight fraction in the melt of CO2-eq, H2O-eq

    """
    wmCO2eq = melt_wf["CO2"] + mdv.species.loc["CO2", "M"] * (
        (melt_wf["CO"] / mdv.species.loc["CO", "M"])
        + (melt_wf["CH4"] / mdv.species.loc["CH4", "M"])
    )
    wmH2Oeq = melt_wf["H2OT"] + mdv.species.loc["H2O", "M"] * (
        (melt_wf["H2"] / mdv.species.loc["H2", "M"])
        + (2.0 * melt_wf["CH4"] / mdv.species.loc["CH4", "M"])
        + (melt_wf["H2S"] / mdv.species.loc["H2S", "M"])
    )
    return wmCO2eq, wmH2Oeq


# calculate weight fraction of species in the gas
def gas_wf(gas_mf, models):
    """
    Convert gas composition from mole fraction to weight fraction.


    Parameters
    ----------
    gas_mf: dict
        Mole fractions of gas composition and total "mass" of gas.

    models: pandas.DataFrame
        Model options.


    Returns
    -------
    dict
        gas composition in weight fraction

    """
    wg_O2 = (mdv.species.loc["O2", "M"] * gas_mf["O2"]) / gas_mf["Xg_t"]
    wg_H2 = (mdv.species.loc["H2", "M"] * gas_mf["H2"]) / gas_mf["Xg_t"]
    wg_H2O = (mdv.species.loc["H2O", "M"] * gas_mf["H2O"]) / gas_mf["Xg_t"]
    wg_H2S = (mdv.species.loc["H2S", "M"] * gas_mf["H2S"]) / gas_mf["Xg_t"]
    wg_S2 = (mdv.species.loc["S2", "M"] * gas_mf["S2"]) / gas_mf["Xg_t"]
    wg_SO2 = (mdv.species.loc["SO2", "M"] * gas_mf["SO2"]) / gas_mf["Xg_t"]
    wg_CO2 = (mdv.species.loc["CO2", "M"] * gas_mf["CO2"]) / gas_mf["Xg_t"]
    wg_CO = (mdv.species.loc["CO", "M"] * gas_mf["CO"]) / gas_mf["Xg_t"]
    wg_CH4 = (mdv.species.loc["CH4", "M"] * gas_mf["CH4"]) / gas_mf["Xg_t"]
    wg_OCS = (mdv.species.loc["OCS", "M"] * gas_mf["OCS"]) / gas_mf["Xg_t"]
    species_X = models.loc["species X", "option"]
    wg_X = (mdv.species.loc[species_X, "M"] * gas_mf["X"]) / gas_mf["Xg_t"]
    result = {
        "wg_O2": wg_O2,
        "wg_H2": wg_H2,
        "wg_H2O": wg_H2O,
        "wg_H2S": wg_H2S,
        "wg_S2": wg_S2,
        "wg_SO2": wg_SO2,
        "wg_CO2": wg_CO2,
        "wg_CO": wg_CO,
        "wg_CH4": wg_CH4,
        "wg_OCS": wg_OCS,
        "wg_X": wg_X,
    }
    return result


# calculate weight fraction of species in the gas relative to total system
def gas_wft(gas_mf, models):
    """
    Calculate gas composition in weight fraction of total system from mole fraction of
    gas.


    Parameters
    ----------
    gas_mf: dict
        Mole fractions of gas composition and total "mass" of gas.

    models: pandas.DataFrame
        Model options.


    Returns
    -------
    dict
        gas composition in weight fraction of total system

    """
    gaswf = gas_wf(gas_mf, models)
    wgt_O2 = gaswf["wg_O2"] * gas_mf["wt_g"]
    wgt_H2 = gaswf["wg_H2"] * gas_mf["wt_g"]
    wgt_H2O = gaswf["wg_H2O"] * gas_mf["wt_g"]
    wgt_H2S = gaswf["wg_H2S"] * gas_mf["wt_g"]
    wgt_S2 = gaswf["wg_S2"] * gas_mf["wt_g"]
    wgt_SO2 = gaswf["wg_SO2"] * gas_mf["wt_g"]
    wgt_CO2 = gaswf["wg_CO2"] * gas_mf["wt_g"]
    wgt_CO = gaswf["wg_CO"] * gas_mf["wt_g"]
    wgt_CH4 = gaswf["wg_CH4"] * gas_mf["wt_g"]
    wgt_OCS = gaswf["wg_OCS"] * gas_mf["wt_g"]
    wgt_X = gaswf["wg_X"] * gas_mf["wt_g"]
    result = {
        "wgt_O2": wgt_O2,
        "wgt_H2": wgt_H2,
        "wgt_H2O": wgt_H2O,
        "wgt_H2S": wgt_H2S,
        "wgt_S2": wgt_S2,
        "wgt_SO2": wgt_SO2,
        "wgt_CO2": wgt_CO2,
        "wgt_CO": wgt_CO,
        "wgt_CH4": wgt_CH4,
        "wgt_OCS": wgt_OCS,
        "wgt_X": wgt_X,
    }
    return result


def gas_weight(gas_mf, bulk_wf):
    """
    Calculates mass of each gas species.


    Parameters
    ----------
    gas_mf: dict
       Mole fractions of gas composition and total "mass" of gas.

    bulk_wf: dict
        Bulk composition of system in weight fraction.


    Returns
    -------
    dict
        mass of each gas species

    """
    gaswft_ = gas_wft(gas_mf)
    Wg_O2 = gaswft_["wgt_O2"] * bulk_wf["Wt"]
    Wg_H2 = gaswft_["wgt_H2"] * bulk_wf["Wt"]
    Wg_H2O = gaswft_["wgt_H2O"] * bulk_wf["Wt"]
    Wg_H2S = gaswft_["wgt_H2S"] * bulk_wf["Wt"]
    Wg_S2 = gaswft_["wgt_S2"] * bulk_wf["Wt"]
    Wg_SO2 = gaswft_["wgt_SO2"] * bulk_wf["Wt"]
    Wg_CO2 = gaswft_["wgt_CO2"] * bulk_wf["Wt"]
    Wg_CO = gaswft_["wgt_CO"] * bulk_wf["Wt"]
    Wg_CH4 = gaswft_["wgt_CH4"] * bulk_wf["Wt"]
    Wg_OCS = gaswft_["wgt_OCS"] * bulk_wf["Wt"]
    Wg_X = gaswft_["wgt_X"] * bulk_wf["Wt"]
    Wg_t = gas_mf["wt_g"] * bulk_wf["Wt"]
    result = {
        "Wg_O2": Wg_O2,
        "Wg_H2": Wg_H2,
        "Wg_H2O": Wg_H2O,
        "Wg_H2S": Wg_H2S,
        "Wg_S2": Wg_S2,
        "Wg_SO2": Wg_SO2,
        "Wg_CO2": Wg_CO2,
        "Wg_CO": Wg_CO,
        "Wg_CH4": Wg_CH4,
        "Wg_OCS": Wg_OCS,
        "Wg_X": Wg_X,
        "Wg_t": Wg_t,
    }
    return result


def gas_moles(gas_mf, bulk_wf, models):
    """
    Calculates moles of each gas species.


    Parameters
    ----------
    gas_mf: dict
        Mole fractions of gas composition and total "mass" of gas.

    bulk_wf: dict
        Bulk composition of system in weight fraction.

    models: pandas.DataFrame
        Model options.


    Returns
    -------
    dict
        moles of each gas species

    """
    gasw = gas_weight(gas_mf, bulk_wf)
    Xg_O2 = gasw["Wg_O2"] / mdv.species.loc["O2", "M"]
    Xg_H2O = gasw["Wg_H2O"] / mdv.species.loc["H2O", "M"]
    Xg_H2 = gasw["Wg_H2"] / mdv.species.loc["H2", "M"]
    Xg_H2S = gasw["Wg_H2S"] / mdv.species.loc["H2S", "M"]
    Xg_S2 = gasw["Wg_S2"] / mdv.species.loc["S2", "M"]
    Xg_SO2 = gasw["Wg_SO2"] / mdv.species.loc["SO2", "M"]
    Xg_CO2 = gasw["Wg_CO2"] / mdv.species.loc["CO2", "M"]
    Xg_CO = gasw["Wg_CO"] / mdv.species.loc["CO", "M"]
    Xg_CH4 = gasw["Wg_CH4"] / mdv.species.loc["CH4", "M"]
    Xg_OCS = gasw["g_OCS"] / mdv.species.loc["OCS", "M"]
    species_X = models.loc["species X", "option"]
    Xg_X = gasw["Wg_X"] / mdv.species.loc[species_X, "M"]
    Xt_g = (
        Xg_O2
        + Xg_H2
        + Xg_H2O
        + Xg_H2S
        + Xg_S2
        + Xg_SO2
        + Xg_CO2
        + Xg_CO
        + Xg_CH4
        + Xg_OCS
        + Xg_X
    )
    result = {
        "Xg_O2": Xg_O2,
        "Xg_H2": Xg_H2,
        "Xg_H2O": Xg_H2O,
        "Xg_H2S": Xg_H2S,
        "Xg_S2": Xg_S2,
        "Xg_SO2": Xg_SO2,
        "Xg_CO2": Xg_CO2,
        "Xg_CO": Xg_CO,
        "Xg_CH4": Xg_CH4,
        "Xg_OCS": Xg_OCS,
        "Xg_X": Xg_X,
        "Xg_t": Xt_g,
    }
    return result


# calculate weight fraction of elements in the gas
def gas_elements(gas_mf, models):
    """
    Convert gas composition in mole fraction of each species to weight fraction of each
    element.


    Parameters
    ----------
    gas_mf: dict
        Mole fractions of gas composition and total "mass" of gas.

    models: pandas.DataFrame
        Model options.


    Returns
    -------
    dict
        weight fraction of each element in the gas

    """
    gaswf = gas_wf(gas_mf, models)
    wg_O = gaswf["wg_O2"] + mdv.species.loc["O", "M"] * (
        (gaswf["wg_H2O"] / mdv.species.loc["H2O", "M"])
        + (2.0 * gaswf["wg_SO2"] / mdv.species.loc["SO2", "M"])
        + (2.0 * gaswf["wg_CO2"] / mdv.species.loc["CO2", "M"])
        + (gaswf["wg_CO"] / mdv.species.loc["CO", "M"])
        + (gaswf["wg_OCS"] / mdv.species.loc["OCS", "M"])
    )
    wg_H = gaswf["wg_H2"] + mdv.species.loc["H", "M"] * (
        (2.0 * gaswf["wg_H2O"] / mdv.species.loc["H2O", "M"])
        + (2.0 * gaswf["wg_H2S"] / mdv.species.loc["H2S", "M"])
        + (4.0 * gaswf["wg_CH4"] / mdv.species.loc["CH4", "M"])
    )
    wg_S = gaswf["wg_S2"] + mdv.species.loc["S", "M"] * (
        (gaswf["wg_H2S"] / mdv.species.loc["H2S", "M"])
        + (gaswf["wg_SO2"] / mdv.species.loc["SO2", "M"])
        + (gaswf["wg_OCS"] / mdv.species.loc["OCS", "M"])
    )
    wg_C = mdv.species.loc["C", "M"] * (
        (gaswf["wg_CO2"] / mdv.species.loc["CO2", "M"])
        + (gaswf["wg_CO"] / mdv.species.loc["CO", "M"])
        + (gaswf["wg_OCS"] / mdv.species.loc["OCS", "M"])
        + (gaswf["wg_CH4"] / mdv.species.loc["CH4", "M"])
    )
    wg_X = gaswf["wg_X"]
    result = {"wg_O": wg_O, "wg_H": wg_H, "wg_S": wg_S, "wg_C": wg_C, "wg_X": wg_X}
    return result


# calculate weight fraction of elements in the melt
def melt_elements(melt_wf, bulk_wf, gas_comp):
    """
    Convert melt composition in weight fraction of each species to weight fraction of
    each element.


    Parameters
    ----------
    melt_wf: dict
        Melt composition in weight fraction.

    bulk_wf: dict
        Bulk composition of system in weight fraction.

    gas_comp: dict
        Mole fractions of gas composition and total "mass" of gas.


    Returns
    -------
    dict
        weight fraction of each element in the melt

    """
    wm_C = mdv.species.loc["C", "M"] * (
        melt_wf["CO2"] / mdv.species.loc["CO2", "M"]
        + melt_wf["CO"] / mdv.species.loc["CO", "M"]
        + melt_wf["CH4"] / mdv.species.loc["CH4", "M"]
    )
    wm_H = (
        2.0
        * mdv.species.loc["H", "M"]
        * (
            melt_wf["H2OT"] / mdv.species.loc["H2O", "M"]
            + melt_wf["H2"] / mdv.species.loc["H2", "M"]
            + melt_wf["H2S"] / mdv.species.loc["H2S", "M"]
            + (2.0 * melt_wf["CH4"]) / mdv.species.loc["CH4", "M"]
        )
    )
    wm_S = mdv.species.loc["S", "M"] * (
        (melt_wf["S2-"] + melt_wf["S6+"]) / mdv.species.loc["S", "M"]
        + (melt_wf["H2S"] / mdv.species.loc["H2S", "M"])
    )
    wm_Fe = bulk_wf["Fe"] / (1.0 - gas_comp["wt_g"])
    Fe32 = overtotal2ratio(melt_wf["Fe3FeT"])
    wm_O = mdv.species.loc["O", "M"] * (
        (melt_wf["H2OT"] / mdv.species.loc["H2O", "M"])
        + ((2.0 * melt_wf["CO2"]) / mdv.species.loc["CO2", "M"])
        + (3.0 * melt_wf["S6+"] / mdv.species.loc["S", "M"])
        + (melt_wf["CO"] / mdv.species.loc["CO", "M"])
        + (wm_Fe / mdv.species.loc["Fe", "M"]) * ((1.5 * Fe32 + 1.0) / (Fe32 + 1.0))
    )
    wm_X = melt_wf["XT"]
    result = {
        "wm_C": wm_C,
        "wm_H": wm_H,
        "wm_S": wm_S,
        "wm_Fe": wm_Fe,
        "wm_O": wm_O,
        "wm_X": wm_X,
    }
    return result


########################################################################################
# WORK IN PROGRESS: VOLUME AND DENSITY #################################################
########################################################################################


# molar volume of the gas in J/bar/mol
def Vm_gas(gas_mf, PT, models):
    # work in progress
    Vm = (
        gas_mf["O2"] * mdv.gas_molar_volume("O2", PT, models)
        + gas_mf["H2"] * mdv.gas_molar_volume("H2", PT, models)
        + gas_mf["H2O"] * mdv.gas_molar_volume("H2O", PT, models)
        + gas_mf["CO2"] * mdv.gas_molar_volume("CO2", PT, models)
        + gas_mf["CO"] * mdv.gas_molar_volume("CO", PT, models)
        + gas_mf["CH4"] * mdv.gas_molar_volume("CH4", PT, models)
        + gas_mf["S2"] * mdv.gas_molar_volume("S2", PT, models)
        + gas_mf["H2S"] * mdv.gas_molar_volume("H2S", PT, models)
        + gas_mf["SO2"] * mdv.gas_molar_volume("SO2", PT, models)
        + gas_mf["OCS"] * mdv.gas_molar_volume("OCS", PT, models)
    )
    return Vm


# volume of the gas in cm3
def gas_volume(PT, gas_mf, bulk_wf, models):
    # work in progress
    Xg_O2, Xg_H2, Xg_H2O, Xg_H2S, Xg_S2, Xg_SO2, Xg_CO2, Xg_CO, Xg_CH4, Xg_OCS, Xt_g = (
        gas_moles(gas_mf, bulk_wf)
    )
    volume = Xt_g * Vm_gas(gas_mf, PT, models)  # in J/bar
    volume_cm3 = volume * 10.0
    return volume_cm3


# density of the gas in g/cm3
def gas_density(PT, gas_mf, bulk_wf, models):
    # work in progress
    volume = gas_volume(PT, gas_mf, bulk_wf, models)  # cm3
    mass = gas_mf["wt_g"] * bulk_wf["Wt"]  # g
    density = mass / volume
    return density


# volume of the melt in cm3
def melt_volume(run, PT, melt_wf, bulk_wf, gas_mf, setup):
    # work in progress
    density = mdv.melt_density(run, PT, melt_wf, setup)
    mass = bulk_wf["Wt"] * (1 - gas_mf["wt_g"])
    volume = mass / density
    return volume


# volume of the system (melt + gas)
def system_density(run, PT, melt_wf, gas_mf, bulk_wf, setup, models):
    # work in progress
    wt_g_ = gas_mf["wt_g"]
    wt_m_ = 1.0 - wt_g_
    density_m = mdv.melt_density(run, PT, melt_wf, setup)
    if wt_g_ > 0.0:
        density_g = gas_density(PT, gas_mf, bulk_wf, models)
    else:
        density_g = 0.0
    density_sys = density_m * wt_m_ + density_g * wt_g_
    return density_sys


########################################################################################
# melt composition #####################################################################
########################################################################################


def melt_comp(run, setup):
    """
    Convert input DataFrame of melt composition into dictionary.


    Parameters
    ----------
    run: float
        Index of row of interest.

    setup: pandas.DataFrame
        Melt composition.


    Returns
    -------
    dict
        melt composition and "initial fO2" definer

    """

    melt_wf = {}

    oxides = setup.columns.tolist()

    for x in [
        "SiO2",
        "TiO2",
        "Al2O3",
        "MgO",
        "MnO",
        "CaO",
        "Na2O",
        "K2O",
        "P2O5",
    ]:
        if x in oxides:
            value = setup.loc[run, x]
        else:
            value = 0.0
        melt_wf[x] = value

    for x in [
        "CO2",
        "X",
        "ST",
    ]:
        if x + "ppm" in oxides:
            value = setup.loc[run, x + "ppm"] / 1000000.0
        else:
            value = 0.0
        melt_wf[x] = value

    if "H2O" in oxides:
        value = setup.loc[run, "H2O"] / 100.0
    else:
        value = 0.0
    melt_wf["H2OT"] = value

    for x in ["FeOT", "Fe2O3T", "FeO", "Fe2O3", "DNNO", "DFMQ"]:
        if x in oxides:
            value = setup.loc[run, x]
        else:
            value = float("NaN")
        melt_wf[x] = value

    for x in ["logfO2", "Fe3FeT", "S6ST"]:
        if x in oxides:
            value = setup.loc[run, x]
        else:
            value = float("NaN")
        melt_wf[x + "_i"] = value

    for x in ["sulf_XFe", "sulf_XCu", "sulf_XNi"]:
        if x in setup:
            melt_wf[x] = setup.loc[run, x]

    melt_wf["XT"] = melt_wf["X"]
    melt_wf["HT"] = (melt_wf["H2OT"] / mdv.species.loc["H2O", "M"]) * mdv.species.loc[
        "H2", "M"
    ]
    melt_wf["CT"] = (melt_wf["CO2"] / mdv.species.loc["CO2", "M"]) * mdv.species.loc[
        "C", "M"
    ]
    melt_wf["H2"] = 0.0

    return melt_wf


# normalise melt composition in weight fraction
def melt_normalise_wf(melt_wf, volatiles, Fe_speciation, molmass="M", majors="majors"):
    """
    Normalise melt composition in weight fraction.


    Parameters
    ----------
    melt_wf: dict
        Weight fraction of melt composition.

    volatiles: str
        Whether to include volatiles in the normalisation:
        "yes" assumes all H is H2O, all C is CO2, all S is S, all X is X;
        "water" assumes all H is H2O;
        "no" volatiles excluded.

    Fe_speciation: str
        How iron speciation is treated in the normalisation:
        "yes" Fe is split into FeO and Fe2O3;
        "no" Fe is treated as FeOT.


    Returns
    -------
    dict
        weight fraction of each oxide in the melt

    """
    if volatiles == "yes":  # assumes all H is H2O, all C is CO2, all S is S, all X is X
        H2O = (melt_wf["HT"] / mdv.species.loc["H2", molmass]) * mdv.species.loc[
            "H2O", molmass
        ]
        CO2 = (melt_wf["CT"] / mdv.species.loc["C", molmass]) * mdv.species.loc[
            "CO2", molmass
        ]
        S = melt_wf["ST"]
        X = melt_wf["XT"]
    elif volatiles == "water":  # assumes all H is H2O
        H2O = (melt_wf["HT"] / mdv.species.loc["H2", molmass]) * mdv.species.loc[
            "H2O", molmass
        ]
        CO2, S, X = 0.0, 0.0, 0.0
    elif volatiles == "no":
        H2O, CO2, S, X = 0.0, 0.0, 0.0, 0.0
    volatiles = H2O + CO2 + S + X
    if Fe_speciation == "no":
        Wm_FeOT_ = Wm_FeOT(melt_wf, molmass=molmass)
        Wm_FeO_ = 0.0
        Wm_Fe2O3_ = 0.0
    elif Fe_speciation == "yes":
        Wm_FeOT_ = 0.0
        Wm_FeO_ = Wm_FeO(melt_wf, molmass=molmass)
        Wm_Fe2O3_ = Wm_Fe2O3(melt_wf, molmass=molmass)
    if mdv.species.loc["SiO2", majors] == "Y":
        Wm_SiO2 = melt_wf["SiO2"]
    else:
        Wm_SiO2 = 0.0
    if mdv.species.loc["TiO2", majors] == "Y":
        Wm_TiO2 = melt_wf["TiO2"]
    else:
        Wm_TiO2 = 0.0
    if mdv.species.loc["Al2O3", majors] == "Y":
        Wm_Al2O3 = melt_wf["Al2O3"]
    else:
        Wm_Al2O3 = 0.0
    if mdv.species.loc["MgO", majors] == "Y":
        Wm_MgO = melt_wf["MgO"]
    else:
        Wm_MgO = 0.0
    if mdv.species.loc["MnO", majors] == "Y":
        Wm_MnO = melt_wf["MnO"]
    else:
        Wm_MnO = 0.0
    if mdv.species.loc["CaO", majors] == "Y":
        Wm_CaO = melt_wf["CaO"]
    else:
        Wm_CaO = 0.0
    if mdv.species.loc["Na2O", majors] == "Y":
        Wm_Na2O = melt_wf["Na2O"]
    else:
        Wm_Na2O = 0.0
    if mdv.species.loc["K2O", majors] == "Y":
        Wm_K2O = melt_wf["K2O"]
    else:
        Wm_K2O = 0.0
    if mdv.species.loc["P2O5", majors] == "Y":
        Wm_P2O5 = melt_wf["P2O5"]
    else:
        Wm_P2O5 = 0.0

    tot = (
        Wm_SiO2
        + Wm_TiO2
        + Wm_Al2O3
        + Wm_FeOT_
        + Wm_FeO_
        + Wm_Fe2O3_
        + Wm_MgO
        + Wm_MnO
        + Wm_CaO
        + Wm_Na2O
        + Wm_K2O
        + Wm_P2O5
    )
    result = {"SiO2": (Wm_SiO2 / tot) * (1.0 - volatiles)}
    result["TiO2"] = (Wm_TiO2 / tot) * (1.0 - volatiles)
    result["Al2O3"] = (Wm_Al2O3 / tot) * (1.0 - volatiles)
    result["FeOT"] = (Wm_FeOT_ / tot) * (1.0 - volatiles)
    result["FeO"] = (Wm_FeO_ / tot) * (1.0 - volatiles)
    result["Fe2O3"] = (Wm_Fe2O3_ / tot) * (1.0 - volatiles)
    result["MgO"] = (Wm_MgO / tot) * (1.0 - volatiles)
    result["MnO"] = (Wm_MnO / tot) * (1.0 - volatiles)
    result["CaO"] = (Wm_CaO / tot) * (1.0 - volatiles)
    result["Na2O"] = (Wm_Na2O / tot) * (1.0 - volatiles)
    result["K2O"] = (Wm_K2O / tot) * (1.0 - volatiles)
    result["P2O5"] = (Wm_P2O5 / tot) * (1.0 - volatiles)
    result["H2O"] = H2O
    result["CO2"] = CO2
    result["S"] = S
    result["X"] = X
    return result


# calculate cation proportions
def melt_cation_proportion(
    melt_wf, volatiles, Fe_speciation, molmass="M", majors="majors"
):
    """
    Calculate cation proportion of melt composition.


    Parameters
    ----------
    melt_wf: dict
        Weight fraction of melt composition.

    volatiles: str
        Whether to include volatiles in the normalisation:
        "yes" assumes all H is H2O, all C is CO2, all S is S, all X is X;
        "water" assumes all H is H2O;
        "no" volatiles excluded.

    Fe_speciation: str
        How iron speciation is treated in the normalisation:
        "yes" Fe is split into FeO and Fe2O3;
        "no" Fe is treated as FeOT.


    Returns
    -------
    dict
        cation mole fraction in the melt

    """
    melt_comp = melt_normalise_wf(
        melt_wf, volatiles, Fe_speciation, molmass=molmass, majors=majors
    )
    tot = (
        (
            (mdv.species.loc["SiO2", "no_cat"] * melt_comp["SiO2"])
            / mdv.species.loc["SiO2", molmass]
        )
        + (
            (mdv.species.loc["TiO2", "no_cat"] * melt_comp["TiO2"])
            / mdv.species.loc["TiO2", molmass]
        )
        + (
            (mdv.species.loc["Al2O3", "no_cat"] * melt_comp["Al2O3"])
            / mdv.species.loc["Al2O3", molmass]
        )
        + (
            (mdv.species.loc["FeO", "no_cat"] * melt_comp["FeOT"])
            / mdv.species.loc["FeO", molmass]
        )
        + (
            (mdv.species.loc["FeO", "no_cat"] * melt_comp["FeO"])
            / mdv.species.loc["FeO", molmass]
        )
        + (
            (mdv.species.loc["Fe2O3", "no_cat"] * melt_comp["Fe2O3"])
            / mdv.species.loc["Fe2O3", molmass]
        )
        + (
            (mdv.species.loc["MgO", "no_cat"] * melt_comp["MgO"])
            / mdv.species.loc["MgO", molmass]
        )
        + (
            (mdv.species.loc["MnO", "no_cat"] * melt_comp["MnO"])
            / mdv.species.loc["MnO", molmass]
        )
        + (
            (mdv.species.loc["CaO", "no_cat"] * melt_comp["CaO"])
            / mdv.species.loc["CaO", molmass]
        )
        + (
            (mdv.species.loc["Na2O", "no_cat"] * melt_comp["Na2O"])
            / mdv.species.loc["Na2O", molmass]
        )
        + (
            (mdv.species.loc["K2O", "no_cat"] * melt_comp["K2O"])
            / mdv.species.loc["K2O", molmass]
        )
        + (
            (mdv.species.loc["P2O5", "no_cat"] * melt_comp["P2O5"])
            / mdv.species.loc["P2O5", molmass]
        )
        + (
            (mdv.species.loc["H2O", "no_cat"] * melt_comp["H2O"])
            / mdv.species.loc["H2O", molmass]
        )
        + (
            (mdv.species.loc["CO2", "no_cat"] * melt_comp["CO2"])
            / mdv.species.loc["CO2", molmass]
        )
    )
    result = {
        "Si": (
            (mdv.species.loc["SiO2", "no_cat"] * melt_comp["SiO2"])
            / mdv.species.loc["SiO2", molmass]
        )
        / tot
    }
    result["Ti"] = (
        (mdv.species.loc["TiO2", "no_cat"] * melt_comp["TiO2"])
        / mdv.species.loc["TiO2", molmass]
    ) / tot
    result["Al"] = (
        (mdv.species.loc["Al2O3", "no_cat"] * melt_comp["Al2O3"])
        / mdv.species.loc["Al2O3", molmass]
    ) / tot
    result["FeT"] = (
        (mdv.species.loc["FeO", "no_cat"] * melt_comp["FeOT"])
        / mdv.species.loc["FeO", molmass]
    ) / tot
    result["Fe2"] = (
        (mdv.species.loc["FeO", "no_cat"] * melt_comp["FeO"])
        / mdv.species.loc["FeO", molmass]
    ) / tot
    result["Fe3"] = (
        (mdv.species.loc["Fe2O3", "no_cat"] * melt_comp["Fe2O3"])
        / mdv.species.loc["Fe2O3", molmass]
    ) / tot
    result["Mg"] = (
        (mdv.species.loc["MgO", "no_cat"] * melt_comp["MgO"])
        / mdv.species.loc["MgO", molmass]
    ) / tot
    result["Mn"] = (
        (mdv.species.loc["MgO", "no_cat"] * melt_comp["MnO"])
        / mdv.species.loc["MnO", molmass]
    ) / tot
    result["Ca"] = (
        (mdv.species.loc["CaO", "no_cat"] * melt_comp["CaO"])
        / mdv.species.loc["CaO", molmass]
    ) / tot
    result["Na"] = (
        (mdv.species.loc["Na2O", "no_cat"] * melt_comp["Na2O"])
        / mdv.species.loc["Na2O", molmass]
    ) / tot
    result["K"] = (
        (mdv.species.loc["K2O", "no_cat"] * melt_comp["K2O"])
        / mdv.species.loc["K2O", molmass]
    ) / tot
    result["P"] = (
        (mdv.species.loc["P2O5", "no_cat"] * melt_comp["P2O5"])
        / mdv.species.loc["P2O5", molmass]
    ) / tot
    result["H"] = (
        (mdv.species.loc["H2O", "no_cat"] * melt_comp["H2O"])
        / mdv.species.loc["H2O", molmass]
    ) / tot
    result["C"] = (
        (mdv.species.loc["CO2", "no_cat"] * melt_comp["CO2"])
        / mdv.species.loc["CO2", molmass]
    ) / tot
    return result


# calculate mole fractions in the melt
def melt_mole_fraction(
    melt_wf, models, volatiles, Fe_speciation, molmass="M", majors="majors"
):
    """
    Calculate oxide mole fraction of the melt.


    Parameters
    ----------
    melt_wf: dict
        Weight fraction of melt composition.

    models: pandas.DataFrame
        Model options.

    volatiles: str
        Whether to include volatiles in the normalisation:
        "yes" assumes all H is H2O, all C is CO2, all S is S, all X is X;
        "water" assumes all H is H2O;
        "no" volatiles excluded.

    Fe_speciation: str
        Wow iron speciation is treated in the normalisation:
        "yes" Fe is split into FeO and Fe2O3;
        "no" Fe is treated as FeOT.


    Returns
    -------
    dict
        oxide mole fraction of the melt

    """
    melt_comp = melt_normalise_wf(
        melt_wf, volatiles, Fe_speciation, molmass=molmass, majors=majors
    )
    species_X = models.loc["species X", "option"]
    mol_tot = (
        (melt_comp["SiO2"] / mdv.species.loc["SiO2", molmass])
        + (melt_comp["TiO2"] / mdv.species.loc["TiO2", molmass])
        + (melt_comp["Al2O3"] / mdv.species.loc["Al2O3", molmass])
        + (melt_comp["FeOT"] / mdv.species.loc["FeO", molmass])
        + (melt_comp["FeO"] / mdv.species.loc["FeO", molmass])
        + (melt_comp["Fe2O3"] / mdv.species.loc["Fe2O3", molmass])
        + (melt_comp["MnO"] / mdv.species.loc["MnO", molmass])
        + (melt_comp["MgO"] / mdv.species.loc["MgO", molmass])
        + (melt_comp["CaO"] / mdv.species.loc["CaO", molmass])
        + (melt_comp["Na2O"] / mdv.species.loc["Na2O", molmass])
        + (melt_comp["K2O"] / mdv.species.loc["K2O", molmass])
        + (melt_comp["P2O5"] / mdv.species.loc["P2O5", molmass])
        + (melt_comp["H2O"] / mdv.species.loc["H2O", molmass])
        + (melt_comp["CO2"] / mdv.species.loc["CO2", molmass])
        + (melt_comp["S"] / mdv.species.loc["S", molmass])
        + (melt_comp["X"] / mdv.species.loc[species_X, molmass])
    )
    result = {"SiO2": (melt_comp["SiO2"] / mdv.species.loc["SiO2", molmass]) / mol_tot}
    result["TiO2"] = (melt_comp["TiO2"] / mdv.species.loc["TiO2", molmass]) / mol_tot
    result["Al2O3"] = (melt_comp["Al2O3"] / mdv.species.loc["Al2O3", molmass]) / mol_tot
    result["FeOT"] = (melt_comp["FeOT"] / mdv.species.loc["FeO", molmass]) / mol_tot
    result["FeO"] = (melt_comp["FeO"] / mdv.species.loc["FeO", molmass]) / mol_tot
    result["Fe2O3"] = (melt_comp["Fe2O3"] / mdv.species.loc["Fe2O3", molmass]) / mol_tot
    result["MnO"] = (melt_comp["MnO"] / mdv.species.loc["MnO", molmass]) / mol_tot
    result["MgO"] = (melt_comp["MgO"] / mdv.species.loc["MgO", molmass]) / mol_tot
    result["CaO"] = (melt_comp["CaO"] / mdv.species.loc["CaO", molmass]) / mol_tot
    result["Na2O"] = (melt_comp["Na2O"] / mdv.species.loc["Na2O", molmass]) / mol_tot
    result["K2O"] = (melt_comp["K2O"] / mdv.species.loc["K2O", molmass]) / mol_tot
    result["P2O5"] = (melt_comp["P2O5"] / mdv.species.loc["P2O5", molmass]) / mol_tot
    result["H2O"] = (melt_comp["H2O"] / mdv.species.loc["H2O", molmass]) / mol_tot
    result["CO2"] = (melt_comp["CO2"] / mdv.species.loc["CO2", molmass]) / mol_tot
    result["S"] = (melt_comp["S"] / mdv.species.loc["S", molmass]) / mol_tot
    result["X"] = (melt_comp["X"] / mdv.species.loc[species_X, molmass]) / mol_tot
    result["mol_tot"] = mol_tot
    return result


def melt_single_O(melt_wf, volatiles, Fe_speciation, molmass="M", majors="majors"):
    """
    Calculate oxide mole fraction on a single oxygen basis in the melt.


    Parameters
    ----------
    melt_wf: dict
        Weight fraction of melt composition.

    volatiles: str
        Whether to include volatiles in the normalisation:
        "yes" assumes all H is H2O, all C is CO2, all S is S, all X is X;
        "water" assumes all H is H2O;
        "no" volatiles excluded.

    Fe_speciation: str
        How iron speciation is treated in the normalisation:
        "yes" Fe is split into FeO and Fe2O3;
        "no" Fe is treated as FeOT.


    Returns
    -------
    dict
        oxide mole fraction on a single oxygen basis in the melt

    """
    melt_comp = melt_normalise_wf(
        melt_wf, volatiles, Fe_speciation, molmass=molmass, majors=majors
    )
    Xmtot = (
        (
            melt_comp["SiO2"]
            / (mdv.species.loc["SiO2", molmass] / mdv.species.loc["SiO2", "no_O"])
        )
        + (
            melt_comp["TiO2"]
            / (mdv.species.loc["TiO2", molmass] / mdv.species.loc["TiO2", "no_O"])
        )
        + (
            melt_comp["Al2O3"]
            / (mdv.species.loc["Al2O3", molmass] / mdv.species.loc["Al2O3", "no_O"])
        )
        + (
            melt_comp["MnO"]
            / (mdv.species.loc["MnO", molmass] / mdv.species.loc["MnO", "no_O"])
        )
        + (
            melt_comp["MgO"]
            / (mdv.species.loc["MgO", molmass] / mdv.species.loc["MgO", "no_O"])
        )
        + (
            melt_comp["CaO"]
            / (mdv.species.loc["CaO", molmass] / mdv.species.loc["CaO", "no_O"])
        )
        + (
            melt_comp["Na2O"]
            / (mdv.species.loc["Na2O", molmass] / mdv.species.loc["Na2O", "no_O"])
        )
        + (
            melt_comp["K2O"]
            / (mdv.species.loc["K2O", molmass] / mdv.species.loc["K2O", "no_O"])
        )
        + (
            melt_comp["P2O5"]
            / (mdv.species.loc["P2O5", molmass] / mdv.species.loc["P2O5", "no_O"])
        )
        + (
            melt_comp["FeOT"]
            / (mdv.species.loc["FeO", molmass] / mdv.species.loc["FeO", "no_O"])
        )
        + (
            melt_comp["FeO"]
            / (mdv.species.loc["FeO", molmass] / mdv.species.loc["FeO", "no_O"])
        )
        + (
            melt_comp["Fe2O3"]
            / (mdv.species.loc["Fe2O3", molmass] / mdv.species.loc["Fe2O3", "no_O"])
        )
        + (
            melt_comp["H2O"]
            / (mdv.species.loc["H2O", molmass] / mdv.species.loc["H2O", "no_O"])
        )
        + (
            melt_comp["CO2"]
            / (mdv.species.loc["CO2", molmass] / mdv.species.loc["CO2", "no_O"])
        )
    )
    result = {
        "SiO2": (
            melt_comp["SiO2"]
            / (mdv.species.loc["SiO2", molmass] / mdv.species.loc["SiO2", "no_O"])
        )
        / Xmtot
    }
    result["TiO2"] = (
        melt_comp["TiO2"]
        / (mdv.species.loc["TiO2", molmass] / mdv.species.loc["TiO2", "no_O"])
    ) / Xmtot
    result["Al2O3"] = (
        melt_comp["Al2O3"]
        / (mdv.species.loc["Al2O3", molmass] / mdv.species.loc["Al2O3", "no_O"])
    ) / Xmtot
    result["FeOT"] = (
        melt_comp["FeOT"]
        / (mdv.species.loc["FeO", molmass] / mdv.species.loc["FeO", "no_O"])
    ) / Xmtot
    result["FeO"] = (
        melt_comp["FeO"]
        / (mdv.species.loc["FeO", molmass] / mdv.species.loc["FeO", "no_O"])
    ) / Xmtot
    result["Fe2O3"] = (
        melt_comp["Fe2O3"]
        / (mdv.species.loc["Fe2O3", molmass] / mdv.species.loc["Fe2O3", "no_O"])
    ) / Xmtot
    result["MnO"] = (
        melt_comp["MnO"]
        / (mdv.species.loc["MnO", molmass] / mdv.species.loc["MnO", "no_O"])
    ) / Xmtot
    result["MgO"] = (
        melt_comp["MgO"]
        / (mdv.species.loc["MgO", molmass] / mdv.species.loc["MgO", "no_O"])
    ) / Xmtot
    result["CaO"] = (
        melt_comp["CaO"]
        / (mdv.species.loc["CaO", molmass] / mdv.species.loc["CaO", "no_O"])
    ) / Xmtot
    result["P2O5"] = (
        melt_comp["P2O5"]
        / (mdv.species.loc["P2O5", molmass] / mdv.species.loc["P2O5", "no_O"])
    ) / Xmtot
    result["Na2O"] = (
        melt_comp["Na2O"]
        / (mdv.species.loc["Na2O", molmass] / mdv.species.loc["Na2O", "no_O"])
    ) / Xmtot
    result["K2O"] = (
        melt_comp["K2O"]
        / (mdv.species.loc["K2O", molmass] / mdv.species.loc["K2O", "no_O"])
    ) / Xmtot
    result["H2O"] = (melt_comp["H2O"] / mdv.species.loc["H2O", molmass]) / Xmtot
    result["CO2"] = (melt_comp["CO2"] / mdv.species.loc["CO2", molmass]) / Xmtot
    result["Xmtot"] = Xmtot
    return result


def melt_pysulfsat(melt_wf):  # output dataframe for pysulfsat
    """
    Convert headers to be read by PySulfSat and makes the dictionary a DataFrame.


    Parameters
    ----------
    melt_wf: dict
        Weight fraction of melt composition.


    Returns
    -------
    pandas.DataFrame
        melt composition in weight fraction

    """
    comp = {
        "SiO2_Liq": [melt_wf["SiO2"]],
        "TiO2_Liq": [melt_wf["TiO2"]],
        "Al2O3_Liq": [melt_wf["Al2O3"]],
        "FeOt_Liq": [Wm_FeOT(melt_wf)],
        "MnO_Liq": [melt_wf["MnO"]],
        "MgO_Liq": [melt_wf["MgO"]],
        "CaO_Liq": [melt_wf["CaO"]],
        "Na2O_Liq": [melt_wf["Na2O"]],
        "K2O_Liq": [melt_wf["K2O"]],
        "P2O5_Liq": [melt_wf["P2O5"]],
    }
    comp = pd.DataFrame(comp)
    return comp
