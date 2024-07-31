# Tests the degassing function works as expected

import VolFe as vf
import pandas as pd
import os


def test_degassing_simple():
    "simple test of degassing function"
    my_analysis = {
        "Sample": "Sari15-04-33",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 47.89,  # wt%
        "TiO2": 0.75,  # wt%
        "Al2O3": 16.74,  # wt%
        "FeOT": 9.43,  # wt%
        "MnO": 0.18,  # wt%
        "MgO": 5.92,  # wt%
        "CaO": 11.58,  # wt%
        "Na2O": 2.14,  # wt%
        "K2O": 0.63,  # wt%
        "P2O5": 0.17,  # wt%
        "H2O": 4.17,  # wt%
        "CO2ppm": 1487.0,  # ppm
        "STppm": 1343.5,  # ppm
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.195,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])
    degas1 = vf.calc_gassing(my_analysis)

    # check that the output is correct here
    # This fails on my laptop after 10 bar

    os.remove("results_gassing_chemistry.csv")
    os.remove("results_jacnewton3.csv")
