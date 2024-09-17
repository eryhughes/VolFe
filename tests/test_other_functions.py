# Tests the calc_pvsat function works as expected

import VolFe as vf
import pandas as pd
import pytest


def test_other_df_y():
    "simple test of calc_fugacity_coefficients function"

    # Define conditions T as a dictionary.
    my_analysis = {
        "Sample": "test",
        "T_C": 1200.0,  # Temperature in 'C
        "P_bar": 1000.0,
    }  # Pressure in bar

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_fugacity_coefficients(my_analysis)

    assert result.loc[0, "yO2"] == pytest.approx(1.209858)
    assert result.loc[0, "yS2"] == pytest.approx(1.191327)
    assert result.loc[0, "yCO2"] == pytest.approx(1.265765)


def test_other_df_C():
    "simple test of calc_sol_consts function"

    # Define the melt composition, fO2 estimate, and T as a dictionary.
    my_analysis = {
        "Sample": "Sari15-04-33",
        "T_C": 1200.0,  # Temperature in 'C
        "P_bar": 1000.0,  # Pressure in bar
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

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_sol_consts(my_analysis)

    assert result.loc[0, "ln[C_CH4]"] == pytest.approx(-4.499333)
    assert result.loc[0, "ln[C_S6+]"] == pytest.approx(29.742992)
    assert result.loc[0, "ln[C_H2OT]"] == pytest.approx(-12.286979)


def test_other_df_puresol():
    "simple test of calc_pure_solubility"

    # Define the melt composition, fO2 estimate, and T as a dictionary.
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
        "Fe3FeT": 0.195,
        "initial_P": 5000.0,
    }  # bar

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_pure_solubility(my_analysis)

    assert result.loc[1, "H2O_wtpc"] == pytest.approx(9.222062871596409)
    assert result.loc[500, "CO2_ppmw"] == pytest.approx(2606.5792533661706)


def test_other_df_puresol_useropt():
    "simple test of calc_pure_solubility with user defined options"

    # Define the melt composition, fO2 estimate, and T as a dictionary.
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
        "Fe3FeT": 0.195,
        "initial_P": 5000.0,
    }  # bar

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["water", "Rhyolite_HughesIP"], ["carbon dioxide", "Rhyolite_Blank93"]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_pure_solubility(my_analysis, models=my_models)

    assert result.loc[1, "H2O_wtpc"] == pytest.approx(10.04115868113018)
    assert result.loc[500, "CO2_ppmw"] == pytest.approx(3025.990363824278)


def test_other_df_isobar():
    "simple test of calc_isobar"

    # Define the melt composition, fO2 estimate, and T as a dictionary.
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
        "Fe3FeT": 0.195,
        "initial_P": 5000.0,
    }  # bar

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # change just the "COH_species" option to "H2O-CO2 only"
    my_models = [["COH_species", "H2O-CO2 only"]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_isobar(
        my_analysis, models=my_models, initial_P=1000.0, final_P=4000.0, step_P=1000.0
    )

    assert result.loc[80, "P_bar"] == pytest.approx(4000)
    assert result.loc[50, "H2O_wtpc"] == pytest.approx(2.1880898834157447)
    assert result.loc[30, "CO2_ppm"] == pytest.approx(870.0216854013381)
