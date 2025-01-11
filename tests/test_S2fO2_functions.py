# Tests the calc_pvsat function works as expected

import VolFe as vf
import pandas as pd
import pytest


def test_S2fO2_df_pvsat():
    "simple test of calc_melt_S_oxybarometer function at Pvsat"

    # Define the melt composition and T as a dictionary.
    my_analysis = {
        "Sample": "Sari15-04-34",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 46.94,  # wt%
        "TiO2": 0.65,  # wt%
        "Al2O3": 16.0,  # wt%
        "FeOT": 8.92,  # wt%
        "MnO": 0.15,  # wt%
        "MgO": 7.3,  # wt%
        "CaO": 13.94,  # wt%
        "Na2O": 1.61,  # wt%
        "K2O": 0.26,  # wt%
        "P2O5": 0.07,  # wt%
        "H2O": 3.83,  # wt%
        "CO2ppm": 1109.0,  # ppm
        "STppm": 1614.12,  # ppm
        "Xppm": 0.0,
    }  # ppm

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result1 = vf.calc_melt_S_oxybarometer(my_analysis)

    assert result1.loc[0, "P (bar) sulf"] == pytest.approx(3209.679275636066)
    assert result1.loc[0, "DFMQ-sulfide"] == pytest.approx(0.9059440119597797)
    assert result1.loc[0, "P (bar) anh"] == ""
    assert result1.loc[0, "DFMQ-sulfate"] == ""


def test_S2fO2_df_P():
    "simple test of calc_melt_S_oxybarometer function at given P"

    # Define the melt composition, P, and T as a dictionary.
    my_analysis = {
        "Sample": "Sari15-04-34",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 46.94,  # wt%
        "TiO2": 0.65,  # wt%
        "Al2O3": 16.0,  # wt%
        "FeOT": 8.92,  # wt%
        "MnO": 0.15,  # wt%
        "MgO": 7.3,  # wt%
        "CaO": 13.94,  # wt%
        "Na2O": 1.61,  # wt%
        "K2O": 0.26,  # wt%
        "P2O5": 0.07,  # wt%
        "H2O": 3.83,  # wt%
        "CO2ppm": 1109.0,  # ppm
        "STppm": 1614.12,  # ppm
        "Xppm": 0.0,  # ppm
        "P_bar": 1000.0,  # bar
        "Fe3FeT": 0.1,
    }

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result = vf.calc_melt_S_oxybarometer(my_analysis)

    assert result.loc[0, "P (bar) sulf"] == pytest.approx(1000.0)
    assert result.loc[0, "DFMQ-sulfide"] == pytest.approx(1.1588609007298638)
    assert result.loc[0, "P (bar) anh"] == pytest.approx(1000.0)
    assert result.loc[0, "DFMQ-sulfate"] == ""


def test_S2fO2_df_Xsulf():
    "simple test of calc_melt_S_oxybarometer function with sulfide composition"

    # Define the melt composition, sulfide composition, and T as a dictionary.
    my_analysis = {
        "Sample": "Sari15-04-34",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 46.94,  # wt%
        "TiO2": 0.65,  # wt%
        "Al2O3": 16.0,  # wt%
        "FeOT": 8.92,  # wt%
        "MnO": 0.15,  # wt%
        "MgO": 7.3,  # wt%
        "CaO": 13.94,  # wt%
        "Na2O": 1.61,  # wt%
        "K2O": 0.26,  # wt%
        "P2O5": 0.07,  # wt%
        "H2O": 3.83,  # wt%
        "CO2ppm": 1109.0,  # ppm
        "STppm": 1614.12,  # ppm
        "Xppm": 0.0,  # ppm
        "sulf_XFe": 0.8,
        "sulf_XCu": 0.15,
        "sulf_XNi": 0.05,
    }

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result = vf.calc_melt_S_oxybarometer(my_analysis)

    assert result.loc[0, "P (bar) sulf"] == pytest.approx(3223.517523236548)
    assert result.loc[0, "DFMQ-sulfide"] == pytest.approx(1.0232307813783903)
    assert result.loc[0, "P (bar) anh"] == ""
    assert result.loc[0, "DFMQ-sulfate"] == ""


def test_S2fO2_df_useropt():
    "simple test of calc_melt_S_oxybarometer function with user options"

    # Define the melt composition, sulfide composition, and T as a dictionary.
    my_analysis = {
        "Sample": "Sari15-04-34",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 46.94,  # wt%
        "TiO2": 0.65,  # wt%
        "Al2O3": 16.0,  # wt%
        "FeOT": 8.92,  # wt%
        "MnO": 0.15,  # wt%
        "MgO": 7.3,  # wt%
        "CaO": 13.94,  # wt%
        "Na2O": 1.61,  # wt%
        "K2O": 0.26,  # wt%
        "P2O5": 0.07,  # wt%
        "H2O": 3.83,  # wt%
        "CO2ppm": 1109.0,  # ppm
        "STppm": 1614.12,  # ppm
        "Xppm": 0.0,  # ppm
        "sulf_XFe": 0.8,
        "sulf_XCu": 0.15,
        "sulf_XNi": 0.05,
    }

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["SCSS", "Fortin15_pss"], ["SCAS", "Chowdhury19"], ["y_S2", "ideal"]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    # runs the calculation
    result = vf.calc_melt_S_oxybarometer(my_analysis, models=my_models)

    assert result.loc[0, "P (bar) sulf"] == pytest.approx(3515.0324378885393)
    assert result.loc[0, "DFMQ-sulfide"] == ""
    assert result.loc[0, "P (bar) anh"] == ""
    assert result.loc[0, "DFMQ-sulfate"] == ""
