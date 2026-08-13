# Tests the degassing function works as expected

import VolFe as vf
import pandas as pd
import pytest
import platform
import numpy as np


# tests to complete
# 2d closed regas
# 2d open regas

options = vf.default_models.copy()
options.loc["output csv", "option"] = False


def test_degas_df_default():
    "simple test of calc_gassing function using example 2a"

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
        "H2O": 1.0,  # wt%
        "CO2ppm": 100.0,  # ppm
        "STppm": 1000.0,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(337.8, rel=1e-3)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.3903631833963219)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(98.79891110755987)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00901325404278971, rel=1e-3)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        -0.483916101878636, rel=1e-3
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.001833243745489121, rel=2e-3
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        0.0007060269775851043, rel=1e-3
    )


def test_degas_df_sat_sulf():
    "simple test of calc_gassing function with sulfur saturation using example 2a"

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
        "H2O": 1.0,  # wt%
        "CO2ppm": 100.0,  # ppm
        "STppm": 1000.0,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # Choose the options I want to change for the calculation
    # - everything else will use the default options
    my_models = [["sulfur_saturation", True]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(337.8089065669246, rel=1e-3)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.3903631833963219, rel=1e-3)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(98.79891110755987, rel=1e-3)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00901325404278971, rel=1e-3)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        -0.483916101878636, rel=1e-3
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.00183, rel=2e-3)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        0.0007060269775851043, rel=1e-3
    )


def test_degas_df_closed_CO2i():
    """simple test of calc_gassing function with closed-system degassing and initial CO2
    using example 2b"""

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
        "H2O": 1.0,  # wt%
        "CO2ppm": 100.0,  # ppm
        "STppm": 1000.0,  # ppm
        "Fe3FeT": 0.177,
        "initial_CO2wtpc": 4.0,
    }  # initial CO2 content of the system in wt%

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["bulk_composition", "melt+vapor_initialCO2"], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(337.8089065669246, rel=1e-3)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.3903631833963219, rel=1e-3)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(98.79891110755987, rel=1e-3)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00901325404278971, rel=1e-3)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        -0.6901667944130256, rel=1e-3
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.1940201296406095, rel=1e-3
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        0.0048362983610370355, rel=1e-3
    )


def test_degas_df_closed_wtg():
    """simple test of calc_gassing function with closed-system degassing and wtg using
    example 2b"""

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
        "H2O": 1.0,  # wt%
        "CO2ppm": 100.0,  # ppm
        "STppm": 1000.0,  # ppm
        "Fe3FeT": 0.177,
        "wt_g": 3.0,
    }  # wt% vapor in equilibrium with the melt

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["bulk_composition", "melt+vapor_wtg"]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(337.8089065669246, rel=1e-3)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.3903631833963219, rel=1e-3)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(98.79891110755987, rel=1e-3)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0090132540427897, rel=1e-3)
    idx = result.index[result["P_bar"] == 2.0]
    assert result.loc[idx, "P_bar"] == 2.0
    assert result.loc[idx, "fO2_DFMQ"] == pytest.approx(
        -0.5457313637076062, rel=1e-3
    )
    assert result.loc[idx, "CO2T_ppmw"] == pytest.approx(
        0.32690414099985654, rel=1e-3
    )
    assert result.loc[idx, "xgS2_mf"] == pytest.approx(
        0.004065465875554149, rel=1e-3
    )


def test_degas_df_open():
    """simple test of calc_gassing function for open-system degassing using example 2c
    but with lower initial volatile content for speed"""

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
        "H2O": 0.5,  # wt%
        "CO2ppm": 10.0,  # ppm
        "STppm": 100,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["gassing_style", "open"], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(46.30360083810482, rel=1e-3)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.3828640585336194, rel=1e-3)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(9.87748349636939, rel=1e-3)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0005557511517431813, rel=1e-3)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        0.3293057401161059, rel=1e-3
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == 0.0
    assert result.loc[len(result) - 1, "xgS2_mf"] == 0.0


def test_regas_df_closed():
    "simple test of calc_gassing function for closed-system regassing using example 2d"

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
        "CO2ppm": 1487,  # ppm
        "STppm": 1343.5,  # ppm
        "Fe3FeT": 0.177,
        "initial_CO2wtpc": 3.0,  # initial CO2 content of the system in wt%
        "final_P": 4000.0,
    }  # bar

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [
        ["gassing_direction", "regas"],
        ["bulk_composition", "melt+vapor_initialCO2"],
        ["output csv", False],
    ]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3869.8535236453786, rel=1e-3)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.4734120855770545, rel=1e-3)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4157582963953, rel=1e-3)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0006611502479165442, rel=1e-3)
    assert result.loc[len(result) - 1, "P_bar"] == 4100.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        0.47225372640648366, rel=1e-3
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        1619.0205394465231, rel=1e-3
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        0.0006938672105624531, rel=1e-3
    )


def test_regas_df_open():
    """simple test of calc_gassing function for open-system regassing using example 2d
    but to 3900 bar to save time"""

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
        "CO2ppm": 1487,  # ppm
        "STppm": 1343.5,  # ppm
        "Fe3FeT": 0.177,
        "initial_CO2wtpc": 3.0,  # initial CO2 content of the system in wt%
        "final_P": 3900.0,
    }  # bar

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [
        ["gassing_direction", "regas"],
        ["gassing_style", "open"],
        ["output csv", False],
    ]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3869.8535236453786, rel=1e-3)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.4734120855770545, rel=1e-3)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4157582963953, rel=1e-3)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0006611502479165442, rel=1e-3)
    assert result.loc[len(result) - 1, "P_bar"] == 3901.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        0.47296926980278187, rel=1e-3
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        1494.1277819055208, rel=1e-3
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        0.0006551773560027822, rel=1e-3
    )


def test_degas_df_CHOAr_basalt():
    "simple test of calc_gassing function for CHOAr system in basalt using example 2e"

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
        "H2O": 1.0,  # wt%
        "CO2ppm": 100.0,  # ppm
        "Xppm": 1.0,  # ppm <<< treating this as Ar
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(314.0767584547407)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.38975633889039596)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(98.79749988350174)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.039848999388785264)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.4882152820352932)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.0020680645692377525
    )
    assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(5.060650006392194e-05)


def test_degas_df_HSO():
    "simple test of calc_gassing function for HSO system using example 2e"

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
        "H2O": 1.0,  # wt%
        "CO2ppm": 0.0,  # ppm
        "STppm": 1000.0,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(121.0752719229832)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.38479697596682794)
    assert result.loc[0, "CO2T_ppmw"] == 0.0
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.02503594688216864)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.4892544872641871)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == 0.0
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(0.000726801036520788)


def test_degas_df_CSO():
    "simple test of calc_gassing function for CSO system using example 2e"

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
        "H2O": 0.0,  # wt%
        "CO2ppm": 100.0,  # ppm
        "STppm": 1000.0,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(241.9094651596429, rel=1e-3)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.3879069754214717, rel=1e-3)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(98.781009072327416, rel=1e-3)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.008259717364763537, rel=1e-3)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        -0.2796046953029556, rel=1e-3
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.04928397455659001, rel=1e-3
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        0.07090392314280476, rel=1e-3
    )
