# isotopes.py

import pandas as pd
import math

import VolFe.melt_gas as mg
import VolFe.calculations as c
import VolFe.model_dependent_variables as mdv

# this file contains functions for
# delta notation
# calculating fractionation factors

####################################
# delta notation ###################
####################################


def delta_standard(standard, isotope, element):
    if element == "S":
        if standard == "VCDT":
            if isotope == 34:
                reference = 1 / 22.6436  # 34S/32S Ding et al. (2001)
    elif element == "C":
        if standard == "VPDB":
            if isotope == 13:
                reference = (
                    0.01123720  # 13C/12C International Atomic Energy Agency (1995)
                )
    elif element == "H":
        if standard == "VSMOW":
            if isotope == 2:
                reference = 155.76 / 1.0e6  # 2H/1H Hagemann et al. (1970)
    return reference


def ratio2delta(standard, isotope, element, ratio):
    reference = delta_standard(standard, isotope, element)
    d = ((ratio - reference) / reference) * 1000.0
    return d


def delta2ratio(standard, isotope, element, d):
    reference = delta_standard(standard, isotope, element)
    ratio = ((d / 1000.0) * reference) + reference
    return ratio


def alpha2Delta(a):
    D = 1000.0 * (1.0 - a)
    return D


#######################################################
# calculating fractionation factors ###################
#######################################################


def alpha_gas_using_beta(element, A, B, PT, models):  # using beta values
    beta_A = mdv.beta_gas(PT, element, A, models)  # gas species A
    beta_B = mdv.beta_gas(PT, element, B, models)  # gas species B
    result = beta_A / beta_B
    return result


#########################
# consistant alphas #####
#########################


def alphas_C(PT, comp, models):
    if float(comp["wt_g_wtpc"].iloc[0]) > 0.0:  # all alphas against CO2 in the vapor
        A = alpha_gas_using_beta("C", "CO", "CO2", PT, models)  # CO(v)-CO2(v)
        B = alpha_gas_using_beta("C", "CH4", "CO2", PT, models)  # CH4(v)-CO2(v)
        C = alpha_gas_using_beta("C", "OCS", "CO2", PT, models)  # OCS(v)-CO2(v)
        D = A / mdv.alpha_C_COv_COm(PT, comp, models)  # COmol(m)-CO2(v)
        E = B / mdv.alpha_C_CH4v_CH4m(PT, comp, models)  # CH4mol(m)-CO2(v)
        F = 1.0 / mdv.alpha_C_CO2v_CO2m(PT, comp, models)  # CO2mol(m)-CO2(v)
        G = 1.0 / mdv.alpha_C_CO2v_CO32mm(PT, comp, models)  # CO32-(m)-CO2(v)
        values = {
            "CO": A,
            "CH4": B,
            "OCS": C,
            "COmol": D,
            "CH4mol": E,
            "CO2mol": F,
            "CO32-": G,
            "CO2": 1.0,
        }
    else:
        a = alpha_gas_using_beta("C", "CO", "CO2", PT, models)  # CO(v)
        b = alpha_gas_using_beta("C", "CH4", "CO2", PT, models)  # CH4(v)
        c = alpha_gas_using_beta("C", "OCS", "CO2", PT, models)  # OCS(v)
        if (
            float(comp["CO2carb_ppmw"].iloc[0]) > 0.0
        ):  # all alphas against CO32- in the melt
            A = mdv.alpha_C_CO2v_CO32mm(PT, comp, models)  # CO32-(m)
            B = c / mdv.alpha_C_CO2v_CO2m(PT, comp, models)  # CO2mol(m)
            C_species1, C_species2 = "CO2mol", "CO32-"
        else:  # all alphas against CO2mol in the melt
            A = mdv.alpha_C_CO2v_CO2m(PT, comp, models)  # CO2mol(m)
            B = c / mdv.alpha_C_CO2v_CO32mm(PT, comp, models)  # CO32-(m)
            C_species1, C_species2 = "CO32-", "CO2mol"
        C = (a * A) / mdv.alpha_C_COv_COm(PT, comp, models)  # COmol(m)
        D = (b * A) / mdv.alpha_C_CH4v_CH4m(PT, comp, models)  # CH4mol(m)
        E = a / A
        F = b / A
        G = c / A
        values = {
            C_species1: A,
            "COmol": B,
            "CH4mol": C,
            "CO2": A,
            "CO": E,
            "CH4": F,
            "OCS": G,
            C_species2: 1.0,
        }
    return values


def alphas_H(PT, comp, models):
    if float(comp["wt_g_wtpc"].iloc[0]) > 0.0:  # all alphas against H2O in the vapor
        A = alpha_gas_using_beta("H", "H2", "H2O", PT, models)  # H2(v)
        B = alpha_gas_using_beta("H", "CH4", "H2O", PT, models)  # CH4(v)
        C = alpha_gas_using_beta("H", "H2S", "H2O", PT, models)  # H2S(v)
        D = A / mdv.alpha_H_H2v_H2m(PT, comp, models)  # H2mol(m)
        E = B / mdv.alpha_H_CH4v_CH4m(PT, comp, models)  # CH4mol(m)
        F = C / mdv.alpha_H_H2Sv_H2Sm(PT, comp, models)  # H2Smol(m)
        G = 1.0 / mdv.alpha_H_H2Ov_H2Om(PT, comp, models)  # H2OT(m)
        H = 1.0 / mdv.alpha_H_H2Ov_OHmm(PT, comp, models)  # H2OT(m)
        values = {
            "H2": A,
            "CH4": B,
            "H2S": C,
            "H2mol": D,
            "CH4mol": E,
            "H2Smol": F,
            "H2Omol": G,
            "OH-": H,
            "H2O": 1.0,
        }
    else:
        a = alpha_gas_using_beta("H", "H2", "H2O", PT, models)  # H2(v)
        b = alpha_gas_using_beta("H", "CH4", "H2O", PT, models)  # CH4(v)
        c = alpha_gas_using_beta("H", "H2S", "H2O", PT, models)  # H2S(v)
        if (
            float(comp["H2Omol_wtpc"].iloc[0]) > 0.0
        ):  # all alphas against H2Omol in the melt
            A = mdv.alpha_H_H2Ov_H2Om(PT, comp, models)  # H2Omol(m)
            B = c / mdv.alpha_H_H2Ov_OHmm(PT, comp, models)  # OH-(m)
            H_species1, H_species2 = "OH-", "H2Omol"
        else:  # all alphas against CO2mol in the melt
            A = mdv.alpha_C_CO2v_CO2m(PT, comp, models)  # OH-(m)
            B = c / mdv.alpha_H_H2Ov_OHmm(PT, comp, models)  # H2Omol(m)
            H_species1, H_species2 = "H2Omol", "OH-"
        C = (a * A) / mdv.alpha_H_H2v_H2m(PT, comp, models)  # H2mol(m)
        D = (b * A) / mdv.alpha_H_CH4v_CH4m(PT, comp, models)  # CH4mol(m)
        E = (c * A) / mdv.alpha_H_CH4v_CH4m(PT, comp, models)  # H2Smol(m)
        F = a / A
        G = b / A
        H = c / A
        values = {
            H_species1: A,
            "H2mol": B,
            "CH4mol": C,
            "H2Smol": D,
            "H2O": A,
            "H2": F,
            "CH4": G,
            "H2S": H,
            H_species2: 1.0,
        }

    return values


def alphas_S(PT, comp, models):  # all alphas against S2-(m)
    C = mdv.alpha_S_H2Sv_S2mm(PT, comp, models)  # H2S(v)
    A = C * alpha_gas_using_beta("S", "S2", "H2S", PT, models)  # S2(v)
    B = C * alpha_gas_using_beta("S", "OCS", "H2S", PT, models)  # OCS(v)
    D = C * alpha_gas_using_beta("S", "SO2", "H2S", PT, models)  # SO2(v)
    E = (C * D) / mdv.alpha_S_SO2v_S6pm(PT, comp, models)  # S6+(m)
    F = mdv.alpha_S_H2Sv_H2Sm(PT, comp, models) / A  # H2S(m)
    values = {
        "S2": A,
        "OCS": B,
        "H2S": C,
        "SO2": D,
        "SO42-": E,
        "H2Smol": F,
        "S2-": 1.0,
    }
    return values


##############
# Simple #####
##############


def simple_isotope_fractionation(D, db):
    """
    Calculates isotopic composition of melt and vapor during closed- and open-system
    degassing for a constant fractionation factor.


    Parameters
    ----------
    D: float
        Float of cap-delta Fractionation factor between vapor and melt in per mil.

    db: float
        Float initial little-delta isotope value of bulk system in per mil.

    Returns
    -------
    Dataframe.

    """
    for n in range(0, 1000, 1):
        F = 1.0 - (n / 1000.0)
        dm_closed = db - D * (1.0 - F)
        dv_closed = db + D * F
        dm_open = db + D * math.log(F)
        dv_open_inst = dm_open + D
        if n == 0.0:
            dv_open = dv_open_inst
        else:
            dv_open = dv_open_inst * (1.0 / n) + ((n - 1.0) / n) * dv_open
        results1 = pd.DataFrame(
            [[F, dm_closed, dv_closed, dm_open, dv_open_inst, dv_open]]
        )
        if n == 0.0:
            results_headers = pd.DataFrame(
                [["F", "dm_closed", "dv_closed", "dm_open", "dv_open_inst", "dv_open"]]
            )
            results = pd.concat([results_headers, results1])
        else:
            results = pd.concat([results, results1])
    results.columns = results.iloc[0]
    results = results[1:]
    return results


#############################
# newton raphson solver #####
#############################


def newton_raphson(x0, constants, e1, step, eqs, deriv, maxiter=100):
    def dx(x, eqs):
        f_ = eqs(x, constants)
        result = abs(0 - f_)
        return result

    def nr(x0, eqs, deriv):
        f_ = eqs(x0, constants)
        df_ = deriv
        x0 = x0 - step * (f_ / df_)
        return x0

    # create results table
    delta1 = dx(x0, eqs)
    results = pd.DataFrame([["guessx", "diff", "step"]])
    results1 = pd.DataFrame([[x0, delta1, step]])
    results = pd.concat([results, results1], ignore_index=True)

    i = 0.0
    for iter in range(maxiter):
        i = i + 1
        f_ = eqs(x0, constants)
        df_ = deriv(x0, constants)
        x0 = x0 - step * (f_ / df_)
        # while x0 < 0.:
        #    step = step/10.
        #    x0 = x0 - step*(f_/df_)
        delta1 = dx(x0, eqs)
        if abs(delta1) < e1:
            return x0
        results1 = pd.DataFrame([[x0, delta1, step]])
        results = pd.concat([results, results1], ignore_index=True)
        if i % 50 == 0:
            results.to_csv("results_nr_isotopes.csv", index=False, header=False)


# two isotopes, nine species


def allocate_species(element, comp, alphas, species_distribution):
    if element == "S":
        species = "S2-"
        T_a = species_distribution[species]
        species = "S2"
        a_b, T_b = alphas[species], species_distribution[species]
        species = "OCS"
        a_c, T_c = alphas[species], species_distribution[species]
        species = "H2S"
        a_d, T_d = alphas[species], species_distribution[species]
        species = "SO2"
        a_e, T_e = alphas[species], species_distribution[species]
        species = "SO42-"
        a_f, T_f = alphas[species], species_distribution[species]
        species = "H2Smol"
        a_g, T_g = alphas[species], species_distribution[species]
        a_h, T_h = 1.0, 0.0
        a_i, T_i = 1.0, 0.0
    if element == "C":
        if float(comp["wt_g_wtpc"].iloc[0]) > 0.0:
            species = "CO2"
            T_a = species_distribution[species]
            species = "CO2mol"
            a_g, T_g = alphas[species], species_distribution[species]
            species = "CO32-"
            a_h, T_h = alphas[species], species_distribution[species]
        else:
            if float(comp["CO2carb_ppmw"].iloc[0]) > 0.0:
                species = "CO32-"
                T_a = species_distribution[species]
                species = "CO2mol"
                a_g, T_g = alphas[species], species_distribution[species]
                species = "CO2"
                a_h, T_h = alphas[species], species_distribution[species]
            else:
                species = "CO2mol"
                T_a = species_distribution[species]
                species = "CO2"
                a_g, T_g = alphas[species], species_distribution[species]
                species = "CO32-"
                a_h, T_h = alphas[species], species_distribution[species]
        species = "CO"
        a_b, T_b = alphas[species], species_distribution[species]
        species = "CH4"
        a_c, T_c = alphas[species], species_distribution[species]
        species = "OCS"
        a_d, T_d = alphas[species], species_distribution[species]
        species = "COmol"
        a_e, T_e = alphas[species], species_distribution[species]
        species = "CH4mol"
        a_f, T_f = alphas[species], species_distribution[species]
        a_i, T_i = 1.0, 0.0
    if element == "H":
        if float(comp["wt_g_wtpc"].iloc[0]) > 0.0:
            species = "H2O"
            T_a = species_distribution[species]
            species = "H2Omol"
            a_h, T_h = alphas[species], species_distribution[species]
            species = "OH-"
            a_i, T_i = alphas[species], species_distribution[species]
        else:
            if float(comp["H2Omol_wtpc"].iloc[0]) > 0.0:
                species = "H2Omol"
                T_a = species_distribution[species]
                species = "H2O"
                a_h, T_h = alphas[species], species_distribution[species]
                species = "OH-"
                a_i, T_i = alphas[species], species_distribution[species]
            else:
                species = "OH-"
                T_a = species_distribution[species]
                species = "H2Omol"
                a_h, T_h = alphas[species], species_distribution[species]
                species = "H2O"
                a_i, T_i = alphas[species], species_distribution[species]
        species = "H2"
        a_b, T_b = alphas[species], species_distribution[species]
        species = "CH4"
        a_c, T_c = alphas[species], species_distribution[species]
        species = "H2S"
        a_d, T_d = alphas[species], species_distribution[species]
        species = "H2mol"
        a_e, T_e = alphas[species], species_distribution[species]
        species = "CH4mol"
        a_f, T_f = alphas[species], species_distribution[species]
        species = "H2Smol"
        a_g, T_g = alphas[species], species_distribution[species]
    alphas_out = {
        "B": a_b,
        "C": a_c,
        "D": a_d,
        "E": a_e,
        "F": a_f,
        "G": a_g,
        "H": a_h,
        "I": a_i,
    }
    species_distribution_out = {
        "A": T_a,
        "B": T_b,
        "C": T_c,
        "D": T_d,
        "E": T_e,
        "F": T_f,
        "G": T_g,
        "H": T_h,
        "I": T_i,
    }
    return alphas_out, species_distribution_out


def rename_output(element, input, comp):
    output = {}
    if element == "S":
        output["m_S2-"] = input["A"]
        output["g_S2"] = input["B"]
        output["g_OCS"] = input["C"]
        output["g_H2S"] = input["D"]
        output["g_SO2"] = input["E"]
        output["m_SO42-"] = input["F"]
        output["m_H2Smol"] = input["G"]
    if element == "C":
        if float(comp["wt_g_wtpc"].iloc[0]) > 0.0:
            output["g_CO2"] = input["A"]
            output["m_CO2mol"] = input["G"]
            output["m_CO32-"] = input["H"]
        else:
            if float(comp["CO2carb_ppmw"].iloc[0]) > 0.0:
                output["m_CO32-"] = input["A"]
                output["m_CO2mol"] = input["G"]
                output["g_CO2"] = input["H"]
            else:
                output["m_CO2mol"] = input["A"]
                output["g_CO2"] = input["G"]
                output["m_CO32-"] = input["H"]
        output["g_CO"] = input["B"]
        output["g_CH4"] = input["C"]
        output["g_OCS"] = input["D"]
        output["m_COmol"] = input["E"]
        output["m_CH4mol"] = input["F"]
    if element == "H":
        if float(comp["wt_g_wtpc"].iloc[0]) > 0.0:
            output["g_H2O"] = input["A"]
            output["m_H2Omol"] = input["H"]
            output["m_OH-"] = input["I"]
        else:
            if float(comp["H2Omol_wtpc"].iloc[0]) > 0.0:
                output["m_H2Omol"] = input["A"]
                output["g_H2O"] = input["H"]
                output["m_OH-"] = input["I"]
            else:
                output["m_OH-"] = input["A"]
                output["m_H2Omol"] = input["H"]
                output["g_H2O"] = input["I"]
        output["g_H2"] = input["B"]
        output["g_CH4"] = input["C"]
        output["g_H2S"] = input["D"]
        output["m_H2mol"] = input["E"]
        output["m_CH4mol"] = input["F"]
        output["m_H2Smol"] = input["G"]

    return output


def i2s9(element, PT, comp, R, models, nr_step, nr_tol):
    if element == "S":
        alphas = alphas_S(PT, comp, models)
        species_distribution = c.mf_S_species(comp)
        R_i = R["S"]
        guessx = iso_initial_guesses(element, R, comp)
    elif element == "C":
        alphas = alphas_C(PT, comp, models)
        species_distribution = c.mf_C_species(comp)
        R_i = R["C"]
        guessx = iso_initial_guesses(element, R, comp)
    elif element == "H":
        alphas = alphas_H(PT, comp, models)
        species_distribution = c.mf_H_species(comp)
        R_i = R["H"]
        guessx = iso_initial_guesses(element, R, comp)

    alphas_, species_distribution_ = allocate_species(
        element, comp, alphas, species_distribution
    )

    constants = alphas_, species_distribution_, R_i

    def isotope_distribution(l_a, constants):
        alphas, species_distribution, R_initial = constants
        R_a = (species_distribution["A"] - l_a) / l_a
        R_b = alphas["B"] * R_a
        R_c = alphas["C"] * R_a
        R_d = alphas["D"] * R_a
        R_e = alphas["E"] * R_a
        R_f = alphas["F"] * R_a
        R_g = alphas["G"] * R_a
        R_h = alphas["H"] * R_a
        R_i = alphas["I"] * R_a
        ratio = {
            "A": R_a,
            "B": R_b,
            "C": R_c,
            "D": R_d,
            "E": R_e,
            "F": R_f,
            "G": R_g,
            "H": R_h,
            "I": R_i,
        }
        return ratio

    def f(l_a, constants):
        alphas, species_distribution, R_i = constants
        R_a = (species_distribution["A"] / l_a) - 1.0
        if species_distribution["B"] > 0.0:
            l_b = species_distribution["B"] / (1.0 + alphas["B"] * R_a)
        else:
            l_b = 0.0
        if species_distribution["C"] > 0.0:
            l_c = species_distribution["C"] / (1.0 + alphas["C"] * R_a)
        else:
            l_c = 0.0
        if species_distribution["D"] > 0.0:
            l_d = species_distribution["D"] / (1.0 + alphas["D"] * R_a)
        else:
            l_d = 0.0
        if species_distribution["E"] > 0.0:
            l_e = species_distribution["E"] / (1.0 + alphas["E"] * R_a)
        else:
            l_e = 0.0
        if species_distribution["F"] > 0.0:
            l_f = species_distribution["F"] / (1.0 + alphas["F"] * R_a)
        else:
            l_f = 0.0
        if species_distribution["G"] > 0.0:
            l_g = species_distribution["G"] / (1.0 + alphas["G"] * R_a)
        else:
            l_g = 0.0
        if species_distribution["H"] > 0.0:
            l_h = species_distribution["H"] / (1.0 + alphas["H"] * R_a)
        else:
            l_h = 0.0
        if species_distribution["I"] > 0.0:
            l_i = species_distribution["I"] / (1.0 + alphas["I"] * R_a)
        else:
            l_i = 0.0
        total = l_a + l_b + l_c + l_d + l_e + l_f + l_g + l_h + l_i
        R_i_ = 1.0 / R_i
        l_t = R_i_ / (R_i_ + 1.0)
        return l_t - total

    def df(l_a, constants):
        alphas, species_distribution, R_i = constants
        a_b, a_c, a_d, a_e, a_f, a_g, a_h, a_i = (
            alphas["B"],
            alphas["C"],
            alphas["D"],
            alphas["E"],
            alphas["F"],
            alphas["G"],
            alphas["H"],
            alphas["I"],
        )
        T_a, T_b, T_c, T_d, T_e, T_f, T_g, T_h, T_i = (
            species_distribution["A"],
            species_distribution["B"],
            species_distribution["C"],
            species_distribution["D"],
            species_distribution["E"],
            species_distribution["F"],
            species_distribution["G"],
            species_distribution["H"],
            species_distribution["I"],
        )
        result = -1.0
        if T_b > 0.0:
            result = result - T_a * T_b * a_b / (
                l_a**2 * (a_b * (T_a / l_a - 1.0) + 1.0) ** 2
            )
        if T_c > 0.0:
            result = result - T_a * T_c * a_c / (
                l_a**2 * (a_c * (T_a / l_a - 1.0) + 1.0) ** 2
            )
        if T_d > 0.0:
            result = result - T_a * T_d * a_d / (
                l_a**2 * (a_d * (T_a / l_a - 1.0) + 1.0) ** 2
            )
        if T_e > 0.0:
            result = result - T_a * T_e * a_e / (
                l_a**2 * (a_e * (T_a / l_a - 1.0) + 1.0) ** 2
            )
        if T_f > 0.0:
            result = result - T_a * T_f * a_f / (
                l_a**2 * (a_f * (T_a / l_a - 1.0) + 1.0) ** 2
            )
        if T_g > 0.0:
            result = result - T_a * T_g * a_g / (
                l_a**2 * (a_g * (T_a / l_a - 1.0) + 1.0) ** 2
            )
        if T_h > 0.0:
            result = result - T_a * T_h * a_h / (
                l_a**2 * (a_h * (T_a / l_a - 1.0) + 1.0) ** 2
            )
        if T_i > 0.0:
            result = result - T_a * T_i * a_i / (
                l_a**2 * (a_i * (T_a / l_a - 1.0) + 1.0) ** 2
            )
        return result

    l_a = newton_raphson(guessx, constants, nr_tol, nr_step, f, df)
    result1 = isotope_distribution(l_a, constants)
    result2 = av_m_g(element, result1, constants)
    if float(comp["wt_g_wtpc"].iloc[0]) == 0.0:
        if element == "C":
            A = result1["A"]
            G = result1["G"]
            H = result1["H"]
            if float(comp["CO2carb_ppmw"].iloc[0]) > 0.0:
                result1["A"] = H
                result1["H"] = A
            else:
                result1["A"] = G
                result1["G"] = A
        if element == "H":
            A = result1["A"]
            H = result1["H"]
            if float(comp["H2Omol_wtpc"].iloc[0]) > 0.0:
                result1["A"] = H
                result1["H"] = A
            else:
                result1["A"] = G
                result1["G"] = A
    renamed_result1 = rename_output(element, result1, comp)
    return renamed_result1, result2


def av_m_g(element, ratio, constants):
    alphas, species_distribution, R_i = constants

    # heavy/total isotope ratio
    R_a_ = ratio["A"] / (ratio["A"] + 1.0)
    R_b_ = ratio["B"] / (ratio["B"] + 1.0)
    R_c_ = ratio["C"] / (ratio["C"] + 1.0)
    R_d_ = ratio["D"] / (ratio["D"] + 1.0)
    R_e_ = ratio["E"] / (ratio["E"] + 1.0)
    R_f_ = ratio["F"] / (ratio["F"] + 1.0)
    R_g_ = ratio["G"] / (ratio["G"] + 1.0)
    R_h_ = ratio["H"] / (ratio["H"] + 1.0)
    R_i_ = ratio["I"] / (ratio["I"] + 1.0)

    if element == "S":
        h_m = (
            R_a_ * species_distribution["A"]
            + R_f_ * species_distribution["F"]
            + R_g_ * species_distribution["G"]
        )  # 34S melt (S2- + SO42- + H2Smol)
        l_m = (
            (1.0 - R_a_) * species_distribution["A"]
            + (1.0 - R_f_) * species_distribution["F"]
            + (1.0 - R_g_) * species_distribution["G"]
        )  # 32S melt (S2- + SO42- + H2Smol)
        h_g = (
            R_b_ * species_distribution["B"]
            + R_c_ * species_distribution["C"]
            + R_d_ * species_distribution["D"]
            + R_e_ * species_distribution["E"]
        )  # 32S gas (H2S + S2 + SO2 + OCS)
        l_g = (
            (1.0 - R_b_) * species_distribution["B"]
            + (1.0 - R_c_) * species_distribution["C"]
            + (1.0 - R_d_) * species_distribution["D"]
            + (1.0 - R_e_) * species_distribution["E"]
        )  # 32S gas (H2S + S2 + SO2 + OCS)
    if element == "C":
        h_m = (
            R_e_ * species_distribution["E"]
            + R_f_ * species_distribution["F"]
            + R_g_ * species_distribution["G"]
            + R_h_ * species_distribution["H"]
        )  # 13C melt (COmol + CH4mol + CO2mol + CO32-mol)
        l_m = (
            (1.0 - R_e_) * species_distribution["E"]
            + (1.0 - R_f_) * species_distribution["F"]
            + (1.0 - R_g_) * species_distribution["G"]
            + (1.0 - R_h_) * species_distribution["H"]
        )  # 12C melt (COmol + CH4mol + CO2mol + CO32-mol)
        h_g = (
            R_b_ * species_distribution["B"]
            + R_c_ * species_distribution["C"]
            + R_d_ * species_distribution["D"]
            + R_a_ * species_distribution["A"]
        )  # 13C gas (CO + CH4 + OCS + CO2)
        l_g = (
            (1.0 - R_b_) * species_distribution["B"]
            + (1.0 - R_c_) * species_distribution["C"]
            + (1.0 - R_d_) * species_distribution["D"]
            + (1.0 - R_a_) * species_distribution["A"]
        )  # 12C gas (CO + CH4 + OCS + CO2)
    if element == "H":
        h_m = (
            R_e_ * species_distribution["E"]
            + R_f_ * species_distribution["F"]
            + R_g_ * species_distribution["G"]
            + R_h_ * species_distribution["H"]
            + R_i_ * species_distribution["I"]
        )  # D melt (H2mol + CH4mol + H2Smol + H2Omol + OH-)
        l_m = (
            (1.0 - R_e_) * species_distribution["E"]
            + (1.0 - R_f_) * species_distribution["F"]
            + (1.0 - R_g_) * species_distribution["G"]
            + (1.0 - R_h_) * species_distribution["H"]
            + (1.0 - R_i_) * species_distribution["I"]
        )  # H melt (H2mol + CH4mol + H2Smol + H2Omol + OH-)
        h_g = (
            R_b_ * species_distribution["B"]
            + R_c_ * species_distribution["C"]
            + R_d_ * species_distribution["D"]
            + R_a_ * species_distribution["A"]
        )  # D gas (H2 + CH4 + H2S + H2O)
        l_g = (
            (1.0 - R_b_) * species_distribution["B"]
            + (1.0 - R_c_) * species_distribution["C"]
            + (1.0 - R_d_) * species_distribution["D"]
            + (1.0 - R_a_) * species_distribution["A"]
        )  # H gas (H2 + CH4 + H2S + H2O)

    R_m = h_m / l_m
    R_g = h_g / l_g
    ratio_g_m = {"R_m": R_m, "R_g": R_g}
    return ratio_g_m


def iso_initial_guesses(element, R, comp):
    if element == "S":
        species_distribution = c.mf_S_species(comp)
        R_i = R["S"]
        A = "S2-"
    elif element == "C":
        species_distribution = c.mf_C_species(comp)
        R_i = R["C"]
        A = "CO2"
    elif element == "H":
        species_distribution = c.mf_H_species(comp)
        R_i = R["H"]
        A = "H2O"
    l_a = species_distribution[A] / (R_i + 1.0)
    return l_a


##################
# OLD ALPHAS #####
##################


def i2s6_S_alphas(PT):  # all alphas against S2- in the melt
    a_b = mg.alpha_H2S_S(PT)  # H2S-S
    a_c = (
        (1.0 / mg.alpha_SO2_SO4(PT))
        * mg.alpha_H2S_S(PT)
        * mg.alpha_gas("S", "SO2", "H2S", PT)
    )  # SO4-S
    a_d = mg.alpha_H2S_S(PT) * mg.alpha_gas("S", "S2", "H2S", PT)  # S2-S
    a_e = mg.alpha_H2S_S(PT) * mg.alpha_gas("S", "SO2", "H2S", PT)  # SO2-S
    a_f = mg.alpha_H2S_S(PT) * mg.alpha_gas("S", "OCS", "H2S", PT)  # OCS-S
    return a_b, a_c, a_d, a_e, a_f


def i2s7_S_alphas(PT):  # all alphas against S2- in the melt
    a_b = mg.alpha_H2S_S(PT)  # H2S-S
    a_c = (
        (1.0 / mg.alpha_SO2_SO4(PT))
        * mg.alpha_H2S_S(PT)
        * mg.alpha_gas("S", "SO2", "H2S", PT)
    )  # SO4-S
    a_d = mg.alpha_H2S_S(PT) * mg.alpha_gas("S", "S2", "H2S", PT)  # S2-S
    a_e = mg.alpha_H2S_S(PT) * mg.alpha_gas("S", "SO2", "H2S", PT)  # SO2-S
    a_f = mg.alpha_H2S_S(PT) * mg.alpha_gas("S", "OCS", "H2S", PT)  # OCS-S
    return a_b, a_c, a_d, a_e, a_f


def alpha_A_B(element, A, B, PT, models):
    if A == "SO2" and B == "H2S":
        a = mg.alpha_gas(element, A, B, PT)
    return a


##########################################
# OLD: different numbers of isotopes #####
##########################################


# two isotopes, two species
def i2s2(element, PT, R_i, melt_wf):
    if element == "S":
        knowns = i2s2_S_melt(PT, R_i, melt_wf)
    # a = fractionation factor A-B, x = mole fraction of S in B, L = mole fraction of
    # light isotope total
    a, x, L = knowns
    A = a - 1.0
    B = (1.0 - a) * (L + x) - 1.0
    C = a * x * L
    l_b = (-B - ((B**2) - (4.0 * A * C)) ** 0.5) / (2.0 * A)
    h_b = x - l_b
    R_b = h_b / l_b
    R_a = a * R_b
    return R_a, R_b


def i2s2_S_melt(PT, R_i, melt_wf):
    a = (
        (1.0 / mg.alpha_SO2_SO4(PT))
        * mg.alpha_H2S_S(PT)
        * mg.alpha_gas("S", "SO2", "H2S", PT)
    )  # SO4-S2-
    x = 1.0 - melt_wf["S6ST"]  # mole fraction of S as S2- in the melt
    R_i_ = 1.0 / R_i["S"]  # 32S/34S
    L = R_i_ / (R_i_ + 1.0)  # mole fraction of 32S
    return a, x, L


# two isotopes, six species
def i2s6(
    element, PT, R, melt_wf, gas_mf, nr_step, nr_tol, guessx
):  # species distribution is mole fraction of S in each species
    if element == "S":
        a_b, a_c, a_d, a_e, a_f = i2s6_S_alphas(PT)
        species_distribution = c.mf_S_species(melt_wf, gas_mf)
        T_a = species_distribution["S2-"]
        T_b = species_distribution["H2S"]
        T_c = species_distribution["SO42-"]
        T_d = species_distribution["S2"]
        T_e = species_distribution["SO2"]
        T_f = species_distribution["OCS"]
        R_i = R["S"]

    constants = a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i

    def isotope_distribution(l_a, constants):
        a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i = constants
        R_a = (T_a - l_a) / l_a  # 34S/32S
        R_b = a_b * R_a
        R_c = a_c * R_a
        R_d = a_d * R_a
        R_e = a_e * R_a
        R_f = a_f * R_a
        return R_a, R_b, R_c, R_d, R_e, R_f

    def av_m_g(element, l_a, constants):
        a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i = constants
        R_a, R_b, R_c, R_d, R_e, R_f = isotope_distribution(l_a, constants)

        # heavy/total isotope ratio
        R_a_ = R_a / (R_a + 1.0)
        R_b_ = R_b / (R_b + 1.0)
        R_c_ = R_c / (R_c + 1.0)
        R_d_ = R_d / (R_d + 1.0)
        R_e_ = R_e / (R_e + 1.0)
        R_f_ = R_f / (R_f + 1.0)

        if element == "S":
            h_m = R_a_ * T_a + R_c_ * T_c  # 34S melt (S2- + SO42-)
            l_m = (1.0 - R_a_) * T_a + (1.0 - R_c_) * T_c
            h_g = (
                R_b_ * T_b + R_d_ * T_d + R_e_ * T_e + R_f_ * T_f
            )  # 34S gas (H2S + S2 + SO2 + OCS)
            l_g = (
                (1.0 - R_b_) * T_b
                + (1.0 - R_d_) * T_d
                + (1.0 - R_e_) * T_e
                + (1.0 - R_f_) * T_f
            )  # 34S gas (H2S + S2 + SO2 + OCS)

        R_m = h_m / l_m
        R_g = h_g / l_g
        return R_m, R_g

    def f(l_a, constants):
        a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i = constants
        R_a = (T_a / l_a) - 1.0
        l_b = T_b / (1.0 + a_b * R_a)
        l_c = T_c / (1.0 + a_c * R_a)
        l_d = T_d / (1.0 + a_d * R_a)
        l_e = T_e / (1.0 + a_e * R_a)
        l_f = T_f / (1.0 + a_f * R_a)
        total = l_a + l_b + l_c + l_d + l_e + l_f
        R_i_ = 1.0 / R_i
        l_t = R_i_ / (R_i_ + 1.0)
        return l_t - total

    def df(l_a, constants):
        a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i = constants
        result = (
            -T_a * T_b * a_b / (l_a**2 * (a_b * (T_a / l_a - 1.0) + 1.0) ** 2)
            - T_a * T_c * a_c / (l_a**2 * (a_c * (T_a / l_a - 1.0) + 1.0) ** 2)
            - T_a * T_d * a_d / (l_a**2 * (a_d * (T_a / l_a - 1.0) + 1.0) ** 2)
            - T_a * T_f * a_f / (l_a**2 * (a_f * (T_a / l_a - 1.0) + 1.0) ** 2)
            - 1.0
        )
        return result

    l_a = newton_raphson(guessx, constants, nr_tol, nr_step, f, df)
    result1 = isotope_distribution(l_a, constants)
    result2 = av_m_g(element, l_a, constants)
    return result1, result2
