import numpy as np
from scipy.optimize import newton

# Copied on 20/02/2025 from:
# https://github.com/pipliggins/EVo/blob/main/evo/solvgas.py
# https://github.com/pipliggins/EVo/blob/main/evo/constants.py

m = {}
"""Atomic and molecular masses, g/mol"""
m["h"] = 1.00794
m["o"] = 15.9994
m["c"] = 12.0000
m["s"] = 32.066
m["n"] = 14.0067
m["si"] = 28.0855
m["ti"] = 47.867
m["al"] = 26.981539
m["mn"] = 54.938044
m["mg"] = 24.305
m["ca"] = 40.078
m["na"] = 22.990
m["k"] = 39.0983
m["p"] = 30.97376
m["li"] = 6.941
m["h2"] = 2.0164
m["o2"] = 31.9988
m["h2o"] = 18.01528
m["co2"] = 44.0095
m["co"] = 28.0101
m["ch4"] = 16.0425
m["so2"] = 64.0638
m["h2s"] = 34.08088
m["s2"] = 64.13
m["ocs"] = 60.0751
m["n2"] = 28.0134
m["fe"] = 55.845
m["sio2"] = 60.0843
m["tio2"] = 79.8658
m["al2o3"] = 101.961278
m["fe2o3"] = 159.6882
m["feo"] = 71.8444
m["feot"] = 71.8444
m["mno"] = 70.937445
m["mgo"] = 40.3044
m["cao"] = 56.0774
m["na2o"] = 61.978938
m["k2o"] = 94.196
m["p2o5"] = 141.944524
m["li2o"] = 29.8814

PTcrit = {}
PTcrit["O2"] = [154.75, 50.7638]
PTcrit["H2O"] = [647.25, 221.1925]
PTcrit["H2"] = [33.25, 12.9696]
PTcrit["CO2"] = [304.15, 73.8659]
PTcrit["CH4"] = [191.05, 46.4069]
PTcrit["CO"] = [133.15, 34.9571]
PTcrit["S2"] = [208.15, 72.954]
PTcrit["SO2"] = [430.95, 78.7295]
PTcrit["OCS"] = [377.55, 65.8612]
PTcrit["H2S"] = [373.55, 90.0779]
PTcrit["N2"] = [126.2, 33.9]  # Roskosz 2006


def find_Y(P, T, species_list):
    """
    Calculates the fugacity coefficients for all relevant species.

    Parameters
    ----------
    P : float
        Pressure (bar)
    T : float
        Temperature (K)
    species_list : [strings]
        A list of each molecule in the system given as capitalised strings.

    Returns
    -------
    tuple
        a tuple of all the fugacity coefficients, as floats.

    References
    ----------
    Shaw, H.R., & Wones, D.R. (1964). Fugacity coefficients for hydrogen
    gas between 0 degrees and 1000C, for pressures to 3000 atm.
    American Journal of Science.

    Holland, T., & Powell, R. (1991). A Compensated-Redlich-Kwong (CORK)
    equation for volumes and fugacities of CO2 and H2O in the range 1 bar
    to 50 kbar and 100-1600C. Contributions to Mineralogy and Petrology.

    Shi, P., & Saxena, S. K. (1992). Thermodynamic modeling of the C-H-O-S
    fluid system. American Mineralogist.
    """

    # Universal gas constant
    R = 8.3144598

    h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, ocs_y, n2_y = (
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    )

    for species in species_list:
        if species == "H2":
            # Shaw and Wones 1964
            logyH2 = (
                P * np.exp((-3.8402 * (T ** (1 / 8))) + 0.541)
                - P**2 * np.exp(-0.1263 * (T**0.5) - 15.98)
                + 300 * np.exp(-0.11901 * T - 5.941) * np.exp(-(P / 300) - 1)
            )

            h2_y = np.exp(logyH2)

        elif species == "CO2":
            # Holland and Powell 1991

            P0 = 5.00  # kbar
            R = R / 1000  # kJ/mol.K
            a = 741.2 + (-0.10891 * T) + (-3.4203e-4 * T**2)
            b = 3.057  # kJ/kbar/mol
            c = -2.26924e-1 + 7.73793e-5 * T
            d = 1.33790e-2 + (-1.01740e-5 * T)
            P_kb = P / 1000.0  # kbar

            lnf = (
                R * T * np.log(1000.0 * P_kb)
                + b * P_kb
                + (a / (b * (T) ** 0.5))
                * (np.log(R * T + b * P_kb) - np.log(R * T + 2.0 * b * P_kb))
                + (2 / 3) * c * P_kb**1.5
                + (d / 2.0) * P_kb**2.0
            ) / (R * T)

            co2_y = np.exp(lnf) / (P_kb * 1000.0)

        elif species == "N2":
            # Holland and Powell 1991 extension to other species

            Tc = PTcrit[species][0]
            Pc = PTcrit[species][1]

            R = R / 1000  # kJ/mol.K
            a = 5.45963e-5 * (Tc ** (5 / 2) / Pc) - 8.6392e-6 * (Tc ** (3 / 2) / Pc) * T
            b = 9.18301e-4 * (Tc / Pc)  # kJ/kbar/mol
            c = -3.30558e-5 * (Tc / Pc ** (3 / 2)) + T * (2.30524e-6 / Pc ** (3 / 2))
            d = 6.93054e-7 * (Tc / Pc**2) + (-8.38293e-8 / Pc**2) * T
            P_kb = P / 1000.0  # kbar

            lnf = (
                R * T * np.log(1000.0 * P_kb)
                + b * P_kb
                + (a / (b * np.sqrt(T)))
                * (np.log(R * T + b * P_kb) - np.log(R * T + 2.0 * b * P_kb))
                + (2 / 3) * c * P_kb**1.5
                + (d / 2.0) * P_kb**2.0
            ) / (R * T)

            n2_y = np.exp(lnf) / (P_kb * 1000.0)

        elif species == "H2O":
            # Holland and Powell 1991 - above 673 K ONLY

            P0 = 2.0
            P_kb = P / 1000.0  # kbar
            R = R / 1000  # KJ/ mol.K
            a = (
                1113.4
                + (-0.22291 * (T - 673))
                + (-3.8022e-4 * (T - 673) ** 2)
                + 1.7791e-7 * (T - 673) ** 3
            )
            b = 1.465  # kJ/kbar/mol
            c = -3.025650e-2 + (-5.343144e-6 * T)
            d = -3.2297554e-3 + (2.2215221e-6 * T)

            def find_V(T, P_kb, P0, R, a, b, c, d, guess):
                def volMRK(
                    T, P_kb, R, a, b, guess
                ):  # Temp, Pressure, H2O/CO2, a or agas, initial guess for NR
                    # PL: having trimmed down, there may just be one real root for this,
                    # negating need to newton algoritham, just use a solve.
                    def f(x):
                        return (
                            P_kb * (x**3)
                            - R * T * x**2
                            - (b * R * T + (b**2) * P_kb - a / T**0.5) * x
                            - (a * b / T**0.5)
                        )

                    # where Vmrk is designated by x

                    def df(x):  # Defines the differentiated eq of Vmrk
                        return (
                            3 * P_kb * x**2
                            - 2 * R * T * x
                            - (b * R * T + b**2 * P_kb - a / T**0.5)
                        )

                    # Carries out NR on eq f, with initial guess 'guess'
                    # and to an allowable error of tol
                    return newton(f, guess, fprime=df)

                assert volMRK(T, P_kb, R, a, b, guess) > 0, "MRK Volume is negative"

                if P_kb > P0:
                    V = (
                        volMRK(T, P_kb, R, a, b, guess)
                        + c * ((P_kb - P0) ** 0.5)
                        + d * (P_kb - P0)
                    )
                else:
                    V = volMRK(T, P_kb, R, a, b, guess)

                assert V > 0, "Total volume is negative"
                return V

            guess = (R * T / P_kb) + b
            Z = (P_kb * find_V(T, P_kb, P0, R, a, b, c, d, guess)) / (R * T)

            if P_kb > P0:
                lnYvirial = (1 / (R * T)) * (
                    ((2 / 3) * c * (P_kb - P0) ** 1.5) + ((d / 2) * (P_kb - P0) ** 2)
                )  # EQ. A.3
            else:
                lnYvirial = 0

            MRKa = a / (b * R * T**1.5)
            MRKb = (b * P_kb) / (R * T)

            h2o_y = np.exp(
                Z - 1 - np.log(Z - MRKb) - MRKa * np.log(1 + (MRKb / Z)) + lnYvirial
            )  # EQ. A.2

        elif species == "SO2":
            # Shi and Saxena 1992

            Pr = P / PTcrit[species][1]
            P0r = 1.0 / PTcrit[species][1]
            Tr = T / PTcrit[species][0]
            A = (
                0.92854
                + 0.43269e-1 * Tr
                + -0.24671 * Tr**-1
                + 0.24999 * Tr**-2
                + -0.53182 * Tr**-3
                + -0.16461e-1 * np.log(Tr)
            )
            B = (
                0.84866e-3
                + -0.18379e-2 * Tr
                + 0.66787e-1 * Tr**-1
                + -0.29427e-1 * Tr**-2
                + 0.29003e-1 * Tr**-3
                + 0.54808e-2 * np.log(Tr)
            )
            C = (
                -0.35456e-3
                + 0.23316e-4 * Tr
                + 0.94159e-3 * Tr**-1
                + -0.81653e-3 * Tr**-2
                + 0.23154e-3 * Tr**-3
                + 0.55542e-4 * np.log(Tr)
            )
            D = 0.0
            integral = (
                A * np.log(Pr / P0r)
                + B * (Pr - P0r)
                + (C / 2.0) * (Pr**2 - P0r**2)
                + (D / 3.0) * (Pr**3 - P0r**3)
            )

            so2_y = (np.exp(integral) + 0) / P

        else:  # For O2, CO, CH4, H2S, S2 and OCS
            # Shi and Saxena 1992

            def q_less(a, b, c, d, e, Tr):
                # when P < 1000, except H2S
                return a + b * Tr**-1 + c * Tr**-1.5 + d * Tr**-3 + e * Tr**-4

            def q_geq(a, b, c, d, e, f, g, h, Tr):
                # when P>= 1000, except H2S
                return (
                    a
                    + b * Tr
                    + c * Tr**-1
                    + d * Tr**2
                    + e * Tr**-2
                    + f * Tr**3
                    + g * Tr**-3
                    + h * np.log(Tr)
                )

            def Z(a, b, c, d, pr, p0r):
                """Returns the SS(1991) integral Z(p, t)/P, eq. 11"""
                return (
                    a * np.log(pr / p0r)
                    + b * (pr - p0r)
                    + (c / 2.0) * (pr**2 - p0r**2)
                    + (d / 3.0) * (pr**3 - p0r**3)
                )

            def Z1000(pcrit, Tr):
                # z at 1000 bar, except H2S
                Pr = 1000 / pcrit
                P0r = 1 / pcrit

                A = 1
                B = (0.09827 / Tr) + (-0.2709 / Tr**3)
                C = (-0.00103 / Tr**1.5) + (0.01427 / Tr**4)
                D = 0.0

                return Z(A, B, C, D, Pr, P0r)  # noqa: B023

            def Z5000(pcrit, Tr):
                # z at 5000 bar, except H2S
                Pr = 5000 / pcrit
                P0r = 1000 / pcrit

                A = 1 + (-0.5917 / Tr**2)
                B = 0.09122 / Tr
                C = (-0.0001416 / Tr**2) + (-0.000002835 * np.log(Tr))
                D = 0.0

                return Z(A, B, C, D, Pr, P0r)  # noqa: B023

            def H2S500(pcrit, Tr):
                # z for H2S at 500 bar
                Pr = 500.0 / pcrit
                P0r = 1.0 / pcrit

                A = q_geq(1.4721, 1.1177, 3.9657, 0, -10.028, 0, 4.5484, -3.8200, Tr)
                B = q_geq(
                    0.16066, 0.10887, 0.29014, 0, -0.99593, 0, -0.18627, -0.45515, Tr
                )
                C = q_geq(
                    -0.28933,
                    -7.0522e-02,
                    0.39828,
                    0,
                    -5.0533e-02,
                    0,
                    0.1176,
                    0.33972,
                    Tr,
                )
                D = 0.0

                return Z(A, B, C, D, Pr, P0r)  # noqa: B023

            def find_Z0(Tr, species):
                # Returns the A, B, C and D parameters for Z, and Z0

                if species != "H2S":
                    if P < 1000.0:
                        P0r = 1.0 / PTcrit[species][1]

                        A = q_less(1, 0, 0, 0, 0, Tr)
                        B = q_less(0, 0.09827, 0, -0.2709, 0, Tr)
                        C = q_less(0, 0, -0.00103, 0, 0.01427, Tr)
                        D = 0
                        Z0 = np.log(1.0)

                    elif P == 1000.0:
                        P0r = 1.0 / PTcrit[species][1]

                        A = 0
                        B = 0
                        C = 0
                        D = 0
                        Z0 = Z1000(PTcrit[species][1], Tr)

                    elif P > 1000.0 and P < 5000.0:
                        P0r = 1000.0 / PTcrit[species][1]

                        A = q_geq(1, 0, 0, 0, -0.5917, 0, 0, 0, Tr)
                        B = q_geq(0, 0, 0.09122, 0, 0, 0, 0, 0, Tr)
                        C = q_geq(0, 0, 0, 0, -0.0001416, 0, 0, -2.835e-06, Tr)
                        D = 0
                        Z0 = Z1000(PTcrit[species][1], Tr)

                    elif P == 5000.0:
                        P0r = 1000.0 / PTcrit[species][1]

                        A = 0
                        B = 0
                        C = 0
                        D = 0
                        Z0 = Z5000(PTcrit[species][1], Tr) + Z1000(
                            PTcrit[species][1], Tr
                        )

                    else:  # P > 5000 bar
                        P0r = 5000.0 / PTcrit[species][1]

                        A = q_geq(2.0614, 0, 0, 0, -2.235, 0, 0, -0.3941, Tr)
                        B = q_geq(0, 0, 0.05513, 0, 0.03934, 0, 0, 0, Tr)
                        C = q_geq(0, 0, -1.894e-06, 0, -1.109e-05, 0, -2.189e-05, 0, Tr)
                        D = q_geq(0, 0, 5.053e-11, 0, 0, -6.303e-21, 0, 0, Tr)
                        Z0 = Z5000(PTcrit[species][1], Tr) + Z1000(
                            PTcrit[species][1], Tr
                        )

                elif species == "H2S":
                    if P < 500.0:
                        P0r = 1.0 / PTcrit[species][1]

                        A = q_geq(
                            1.4721, 1.1177, 3.9657, 0, -10.028, 0, 4.5484, -3.8200, Tr
                        )
                        B = q_geq(
                            0.16066,
                            0.10887,
                            0.29014,
                            0,
                            -0.99593,
                            0,
                            -0.18627,
                            -0.45515,
                            Tr,
                        )
                        C = q_geq(
                            -0.28933,
                            -7.0522e-02,
                            0.39828,
                            0,
                            -5.0533e-02,
                            0,
                            0.1176,
                            0.33972,
                            Tr,
                        )
                        D = 0.0
                        Z0 = np.log(1.0)

                    elif P == 500.0:
                        P0r = 1.0 / PTcrit[species][1]

                        A = 0.0
                        B = 0.0
                        C = 0.0
                        D = 0.0
                        Z0 = H2S500(PTcrit["H2S"][1], Tr)

                    elif P > 500.0:
                        P0r = 500.0 / PTcrit["H2S"][1]

                        A = q_geq(
                            0.59941,
                            -1.557e-03,
                            0.04525,
                            0,
                            0.36687,
                            0,
                            -0.79248,
                            0.26058,
                            Tr,
                        )
                        B = q_geq(
                            2.2545e-02,
                            1.7473e-03,
                            0.048253,
                            0,
                            -0.01989,
                            0,
                            0.032794,
                            -0.010985,
                            Tr,
                        )
                        C = q_geq(
                            5.7375e-04,
                            -2.0944e-06,
                            -0.0011894,
                            0,
                            0.0014661,
                            0,
                            -0.00075605,
                            -0.00027985,
                            Tr,
                        )
                        D = 0.0
                        Z0 = H2S500(PTcrit["H2S"][1], Tr)

                return A, B, C, D, P0r, Z0

            def y(P, T, species):
                """
                Main equation for calculating fugacity coefficients for a
                single `species` using Shi & Saxena 1992.

                Parameters
                ----------
                P : float
                    Pressure (bar)
                T : float
                    Temperature (K)
                species : string
                    The species name

                Returns
                -------
                float
                    fugacity coefficient for `species`.
                """

                # only calibrated to 1 bar, they start increasing gamma again < 1 bar
                # which should happen - gases become more ideal at low pressure.
                if P < 1.0:
                    return 1

                Tr = T / PTcrit[species][0]
                Pr = P / PTcrit[species][1]

                A, B, C, D, P0r, Z0 = find_Z0(Tr, species)

                intZ = Z(A, B, C, D, Pr, P0r)  # noqa: B023

                return np.exp(intZ + Z0) / P

            if species == "O2":
                o2_y = y(P, T, "O2")
            elif species == "CO":
                co_y = y(P, T, "CO")
            elif species == "CH4":
                ch4_y = y(P, T, "CH4")
            elif species == "S2":
                s2_y = y(P, T, "S2")
            elif species == "H2S":
                h2s_y = y(P, T, "H2S")
            elif species == "OCS":
                ocs_y = y(P, T, "OCS")
            # elif species == 'N2':
            #     n2_y = y(P, T, 'N2')

    return h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y, ocs_y


# Copied on 20/02/2025 from https://github.com/sdecho/Sulfur_X/blob/main/fugacity.py
# Input changed to PT and added required P0.
def phiso2(PT):
    P0 = 1.0
    # This function calculates the fugacity coefficient of SO2 following Shi & Saxena 1992. To avoid discontinuity
    # in the result, only the calculation at high pressure is used. Instead, fugacity coefficient is assumed to be 1
    # at pressure lower than 20bar.
    if PT["P"] < 20:
        lnphiSO2 = 0
    else:
        Tcr = 430.95  # critical temperature in K
        Pcr = 78.7295  # critical pressure in bar
        Pr = PT["P"] / Pcr  # reduced pressure
        P0r = P0 / Pcr
        Tr = PT["T"] / Tcr  # reduced temperature

        # Paramatersfor EOS
        AQ1 = 0.92854e00
        AQ2 = 0.43269e-1
        AQ3 = -0.24671e00
        AQ4 = 0
        AQ5 = 0.24999e00
        AQ6 = 0
        AQ7 = -0.53182e00
        AQ8 = -0.16461e-01
        BQ1 = 0.84866e-03
        BQ2 = -0.18379e-02
        BQ3 = 0.66787e-01
        BQ4 = 0
        BQ5 = -0.29427e-01
        BQ6 = 0
        BQ7 = 0.29003e-01
        BQ8 = 0.54808e-02
        CQ1 = -0.35456e-03
        CQ2 = 0.23316e-04
        CQ3 = 0.94159e-03
        CQ4 = 0
        CQ5 = -0.81653e-03
        CQ6 = 0
        CQ7 = 0.23154e-03
        CQ8 = 0.55542e-04

        # ---------EOS from Shi & Saxena 1992 - --------------------------------------
        # Equation (3a) from Shi & Saxena 1992
        A = (
            AQ1
            + AQ2 * Tr
            + AQ3 * (Tr**-1)
            + AQ4 * (Tr**2)
            + AQ5 * (Tr**-2)
            + AQ6 * (Tr**3)
            + AQ7 * (Tr**-3)
            + AQ8 * np.log(Tr)
        )
        B = (
            BQ1
            + BQ2 * Tr
            + BQ3 * (Tr**-1)
            + BQ4 * (Tr**2)
            + BQ5 * (Tr**-2)
            + BQ6 * (Tr**3)
            + BQ7 * (Tr**-3)
            + BQ8 * np.log(Tr)
        )
        C = (
            CQ1
            + CQ2 * Tr
            + CQ3 * (Tr**-1)
            + CQ4 * (Tr**2)
            + CQ5 * (Tr**-2)
            + CQ6 * (Tr**3)
            + CQ7 * (Tr**-3)
            + CQ8 * np.log(Tr)
        )
        # Equation (10) from Shi & Saxena 1992
        Zcomp = A * np.log(Pr / P0r) + B * (Pr - P0r) + (C / 2) * (Pr**2 - P0r**2)
        # Equation (9) from Shi & Saxena 1992
        lnphiSO2 = Zcomp - np.log(PT["P"])

    return np.exp(lnphiSO2)
