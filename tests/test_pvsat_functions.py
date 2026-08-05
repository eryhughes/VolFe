# Tests the calc_pvsat function works as expected

import VolFe as vf
import pandas as pd
import pytest


def test_pvsat_df_FeOT_Fe3FeT():
    "simple test of calc_pvsat function using FeOT and Fe3FeT using example 1a"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 56.98,  # wt%
        "TiO2": 1.66,  # wt%
        "Al2O3": 15.52,  # wt%
        "FeOT": 9.47,  # wt%
        "MnO": 0.24,  # wt%
        "MgO": 2.96,  # wt%
        "CaO": 6.49,  # wt%
        "Na2O": 4.06,  # wt%
        "K2O": 0.38,  # wt%
        "P2O5": 0.22,  # wt%
        "H2O": 1.88,  # wt%
        "CO2ppm": 13.0,  # ppm
        "STppm": 362.83,  # ppm
        "Fe3FeT": 0.155,
    }  # mole or weight fraction (they're the same)

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0, "P_bar"] == pytest.approx(330.781731652385)
    assert result1.loc[0, "fO2_DFMQ"] == pytest.approx(0.28960414241780263)


def test_pvsat_df_Fe2O3T_DFMQ():
    "simple test of calc_pvsat function using Fe2O3T and DFMQ using example 1a"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,
        "SiO2": 56.98,
        "TiO2": 1.66,
        "Al2O3": 15.52,
        "Fe2O3T": 9.47,  # Fe2O3T instead of FeOT
        "MnO": 0.24,
        "MgO": 2.96,
        "CaO": 6.49,
        "Na2O": 4.06,
        "K2O": 0.38,
        "P2O5": 0.22,
        "H2O": 1.88,
        "CO2ppm": 13.0,
        "STppm": 362.83,
        "DFMQ": 0.0,
    }  # DFMQ instead of Fe3+/FeT

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0, "P_bar"] == pytest.approx(324.46984008025794)
    assert result1.loc[0, "Fe3+/FeT"] == pytest.approx(0.13813893016613535)


def test_pvsat_df_FeO_Fe2O3():
    "simple test of calc_pvsat function using FeO and Fe2O3 using example 1a"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,
        "SiO2": 56.98,
        "TiO2": 1.66,
        "Al2O3": 15.52,
        "FeO": 8.03,  # FeO and Fe2O3 to contraint total Fe and fO2
        "Fe2O3": 2.57,  # ^
        "MnO": 0.24,
        "MgO": 2.96,
        "CaO": 6.49,
        "Na2O": 4.06,
        "K2O": 0.38,
        "P2O5": 0.22,
        "H2O": 1.88,
        "CO2ppm": 13.0,
        "STppm": 362.83,
    }

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0, "P_bar"] == pytest.approx(421.138400557486)
    assert result1.loc[0, "fO2_DFMQ"] == pytest.approx(1.287568712477623)


def test_pvsat_df_FeOT_DNNO():
    "simple test of calc_pvsat function using FeOT and DNNO using example 1a"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,
        "SiO2": 56.98,
        "TiO2": 1.66,
        "Al2O3": 15.52,
        "FeOT": 9.47,
        "MnO": 0.24,
        "MgO": 2.96,
        "CaO": 6.49,
        "Na2O": 4.06,
        "K2O": 0.38,
        "P2O5": 0.22,
        "H2O": 1.88,
        "CO2ppm": 13.0,
        "STppm": 362.83,
        "DNNO": 1.0,
    }

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0, "P_bar"] == pytest.approx(477.3369430755751)
    assert result1.loc[0, "Fe3+/FeT"] == pytest.approx(0.26725365489217384)


def test_pvsat_df_FeOT_S6ST():
    "simple test of calc_pvsat function using FeOT and S6+/ST using example 1a"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,
        "SiO2": 56.98,
        "TiO2": 1.66,
        "Al2O3": 15.52,
        "FeOT": 9.47,
        "MnO": 0.24,
        "MgO": 2.96,
        "CaO": 6.49,
        "Na2O": 4.06,
        "K2O": 0.38,
        "P2O5": 0.22,
        "H2O": 1.88,
        "CO2ppm": 13.0,
        "STppm": 362.83,
        "S6ST": 0.23,
    }

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0, "P_bar"] == pytest.approx(412.80473362888443)
    assert result1.loc[0, "fO2_DFMQ"] == pytest.approx(1.0229078685870343)


def test_pvsat_df_FeOT_DNNO_Fe3FeT():
    "simple test of calc_pvsat function using FeOT, DNNO, and Fe3+/FeT"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,
        "SiO2": 56.98,
        "TiO2": 1.66,
        "Al2O3": 15.52,
        "FeOT": 9.47,
        "MnO": 0.24,
        "MgO": 2.96,
        "CaO": 6.49,
        "Na2O": 4.06,
        "K2O": 0.38,
        "P2O5": 0.22,
        "H2O": 1.88,
        "CO2ppm": 13.0,
        "STppm": 362.83,
        "DNNO": 1.0,
        "Fe3FeT": 0.155,
    }

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation, check a warning is generated
    with pytest.warns(
        UserWarning, match="you entered more than one way to infer iron speciation"
    ):
        result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0, "P_bar"] == pytest.approx(330.781731652385)
    assert result1.loc[0, "Fe3+/FeT"] == pytest.approx(0.155)


def test_pvsat_df_FeOT_Fe3FeT_useroptions():
    """simple test of calc_pvsat function using FeOT and Fe3FeT with user defined
    options using example 1c"""

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 56.98,  # wt%
        "TiO2": 1.66,  # wt%
        "Al2O3": 15.52,  # wt%
        "FeOT": 9.47,  # wt%
        "MnO": 0.24,  # wt%
        "MgO": 2.96,  # wt%
        "CaO": 6.49,  # wt%
        "Na2O": 4.06,  # wt%
        "K2O": 0.38,  # wt%
        "P2O5": 0.22,  # wt%
        "H2O": 1.88,  # wt%
        "CO2ppm": 13.0,  # ppm
        "STppm": 362.83,  # ppm
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.155,
    }  # mole or weight fraction (they're the same)

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [
        ["carbon dioxide", "Basalt_Dixon97"],
        ["hydrogen sulfide", "BasalticAndesite_Hughes24"],
        ["y_S2", "ideal"],
    ]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result1 = vf.calc_Pvsat(my_analysis, models=my_models)

    assert result1.loc[0, "P_bar"] == pytest.approx(279.88640811865804)
    assert result1.loc[0, "fO2_DFMQ"] == pytest.approx(0.2883015061880503)


def test_pvsat_df_X_Ar_bas():
    "simple test of calc_pvsat function including X as Ar in basalt using example 1d"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 56.98,  # wt%
        "TiO2": 1.66,  # wt%
        "Al2O3": 15.52,  # wt%
        "FeOT": 9.47,  # wt%
        "MnO": 0.24,  # wt%
        "MgO": 2.96,  # wt%
        "CaO": 6.49,  # wt%
        "Na2O": 4.06,  # wt%
        "K2O": 0.38,  # wt%
        "P2O5": 0.22,  # wt%
        "H2O": 1.88,  # wt%
        "CO2ppm": 13.0,  # ppm
        "STppm": 362.83,  # ppm
        "Xppm": 20.0,  # ppm *** 20 ppm "X" added**
        "Fe3FeT": 0.155,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0, "P_bar"] == pytest.approx(580.1382319929767)
    assert result1.loc[0, "fO2_DFMQ"] == pytest.approx(0.29594295301335194)


def test_pvsat_df_X_Ar_rhy():
    "simple test of calc_pvsat function including X as Ar in rhyolite using example 1d"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 56.98,  # wt%
        "TiO2": 1.66,  # wt%
        "Al2O3": 15.52,  # wt%
        "FeOT": 9.47,  # wt%
        "MnO": 0.24,  # wt%
        "MgO": 2.96,  # wt%
        "CaO": 6.49,  # wt%
        "Na2O": 4.06,  # wt%
        "K2O": 0.38,  # wt%
        "P2O5": 0.22,  # wt%
        "H2O": 1.88,  # wt%
        "CO2ppm": 13.0,  # ppm
        "STppm": 362.83,  # ppm
        "Xppm": 20.0,  # ppm *** 20 ppm "X" added**
        "Fe3FeT": 0.155,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["species X solubility", "Ar_Rhyolite_Hughes25"]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis, models=my_models)

    assert result1.loc[0, "P_bar"] == pytest.approx(376.0555328591023)
    assert result1.loc[0, "fO2_DFMQ"] == pytest.approx(0.2907603795983995)


def test_pvsat_df_X_Ne_bas():
    "simple test of calc_pvsat function including X as Ne in basalt using example 1d"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 56.98,  # wt%
        "TiO2": 1.66,  # wt%
        "Al2O3": 15.52,  # wt%
        "FeOT": 9.47,  # wt%
        "MnO": 0.24,  # wt%
        "MgO": 2.96,  # wt%
        "CaO": 6.49,  # wt%
        "Na2O": 4.06,  # wt%
        "K2O": 0.38,  # wt%
        "P2O5": 0.22,  # wt%
        "H2O": 1.88,  # wt%
        "CO2ppm": 13.0,  # ppm
        "STppm": 362.83,  # ppm
        "Xppm": 20.0,  # ppm *** 20 ppm "X" added**
        "Fe3FeT": 0.155,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["species X", "Ne"], ["species X solubility", "Ne_Basalt_Hughes25"]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis, models=my_models)

    assert result1.loc[0, "P_bar"] == pytest.approx(463.21887716871447)
    assert result1.loc[0, "fO2_DFMQ"] == pytest.approx(0.29297974431457785)


def test_pvsat_df_X_Ne_rhy():
    "simple test of calc_pvsat function including X as Ne in rhyolite using example 1d"

    my_analysis = {
        "Sample": "TN273-01D-01-01",
        "T_C": 1200.0,  # Temperature in 'C
        "SiO2": 56.98,  # wt%
        "TiO2": 1.66,  # wt%
        "Al2O3": 15.52,  # wt%
        "FeOT": 9.47,  # wt%
        "MnO": 0.24,  # wt%
        "MgO": 2.96,  # wt%
        "CaO": 6.49,  # wt%
        "Na2O": 4.06,  # wt%
        "K2O": 0.38,  # wt%
        "P2O5": 0.22,  # wt%
        "H2O": 1.88,  # wt%
        "CO2ppm": 13.0,  # ppm
        "STppm": 362.83,  # ppm
        "Xppm": 20.0,  # ppm *** 20 ppm "X" added**
        "Fe3FeT": 0.155,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["species X", "Ne"], ["species X solubility", "Ne_Rhyolite_Hughes25"]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis, models=my_models)

    assert result1.loc[0, "P_bar"] == pytest.approx(354.317755081519)
    assert result1.loc[0, "fO2_DFMQ"] == pytest.approx(0.2902055195533597)
