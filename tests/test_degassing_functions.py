# Tests the degassing function works as expected

import VolFe as vf
import pandas as pd
import pytest


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
        "H2O": 4.17,  # wt%
        "CO2ppm": 1487.0,  # ppm
        "STppm": 1343.5,  # ppm
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(3862.9286026262585)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47326338580309724)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4121779037337)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.000471000320606681)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.3686544338466735)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.00624264256020586
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        4.663409944151297e-05
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
        "H2O": 4.17,  # wt%
        "CO2ppm": 1487.0,  # ppm
        "STppm": 1343.5,  # ppm
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # Choose the options I want to change for the calculation
    # - everything else will use the default options
    my_models = [["sulfur_saturation", "True"]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3851.0721730969226)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.4730086622181213)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4061072002758)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00021511891574073062)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        -0.09215694434160326
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.006396595215163005
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        6.473394019681534e-06
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
        "H2O": 4.17,  # wt%
        "CO2ppm": 1487.0,  # ppm
        "STppm": 1343.5,  # ppm
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.177,
        "initial_CO2wtpc": 4.0,
    }  # initial CO2 content of the system in wt%

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["bulk_composition", "melt+vapor_initialCO2"], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3862.9286026262585)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47326338580309724)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4121779037337)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.000471000320606681)
    assert result.loc[len(result) - 1, "P_bar"] == 3.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        -0.30023881436684974
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.33585401816719396
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        0.00010126638630802683
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
        "H2O": 4.17,  # wt%
        "CO2ppm": 1487.0,  # ppm
        "STppm": 1343.5,  # ppm
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.177,
        "wt_g": 3.0,
    }  # wt% vapor in equilibrium with the melt

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["bulk_composition", "melt+vapor_wtg"]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3862.9286026262585)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47326338580309724)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4121779037337)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.000471000320606681)
    assert result.loc[len(result) - 1, "P_bar"] == 2.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.328178420366644)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.16651526559195437
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        7.805550256915342e-05
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
        "H2O": 1.0,  # wt%
        "CO2ppm": 50.0,  # ppm
        "STppm": 100,  # ppm
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["gassing_style", "open"], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(197.81717165883674)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.3867740973142686)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(49.395163028754354)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00011877068634021532)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.4171085)
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
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.177,
        "initial_CO2wtpc": 3.0,  # initial CO2 content of the system in wt%
        "final_P": 5000.0,
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

    assert result.loc[0, "P_bar"] == pytest.approx(3862.9286026262585)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47326338580309724)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4121779037337)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.000471000320606681)
    assert result.loc[len(result) - 1, "P_bar"] == 5100.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.4594909149353912)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(2334.2395489689416)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        0.00040847652169019897
    )


def test_regas_df_open():
    """simple test of calc_gassing function for open-system regassing using example 2d
    but to 4000 bar to save time"""

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
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.177,
        "initial_CO2wtpc": 3.0,  # initial CO2 content of the system in wt%
        "final_P": 4000.0,
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

    assert result.loc[0, "P_bar"] == pytest.approx(3862.9286026262585)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47326338580309724)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4121779037337)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.000471000320606681)
    assert result.loc[len(result) - 1, "P_bar"] == 4001.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.4696895464709012)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(1569.180690253006)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(0.000450238483173885)


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
        "H2O": 2.0,  # wt%
        "CO2ppm": 500.0,  # ppm
        "STppm": 0.0,  # ppm
        "Xppm": 10.0,  # ppm <<< treating this as Ar
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(1466.7665158846582)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.4184782384942274)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(494.3819055099977)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.0853281310975796)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.5808773404714334)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.004933568211088988
    )
    assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(
        0.00024081249606827505
    )


def test_degas_df_CHONe_basalt():
    "simple test of calc_gassing function for CHONe system in basalt using example 2e"

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
        "H2O": 2.0,  # wt%
        "CO2ppm": 500.0,  # ppm
        "STppm": 0.0,  # ppm
        "Xppm": 10.0,  # ppm <<< treating this as Ar
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [
        ["species X", "Ne"],
        ["species X solubility", "Ne_Basalt_HughesIP"],
        ["output csv", False],
    ]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    # run calculation
    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(1411.3395702781336)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.4171323423753375)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(494.3664109356707)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.04711081804999244)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.5807683081320834)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.004932307918374528
    )
    assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(0.0004765928303183229)


def test_degas_df_CHOAr_rhyolite():
    "simple test of calc_gassing function for CHOAr system in rhyolite using example 2e"

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
        "H2O": 2.0,  # wt%
        "CO2ppm": 500.0,  # ppm
        "STppm": 0.0,  # ppm
        "Xppm": 10.0,  # ppm <<< treating this as Ar
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [
        ["species X solubility", "Ar_Rhyolite_HughesIP"],
        ["output csv", False],
    ]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(1369.96395460828)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.4161253277862338)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(494.3548019730171)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.01658968664892445)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.5806861266434034)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.004931358248440175
    )
    assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(
        0.00024070388121406465
    )


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
        "H2O": 2.0,  # wt%
        "CO2ppm": 0.0,  # ppm
        "STppm": 1000.0,  # ppm
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(380.16027020784526)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.39144451230686617)
    assert result.loc[0, "CO2T_ppmw"] == 0.0
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0055645408964144945)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.3838932210138193)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == 0.0
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        0.00012487825695314437
    )


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
        "CO2ppm": 500.0,  # ppm
        "STppm": 1000.0,  # ppm
        "Xppm": 0.0,  # ppm
        "Fe3FeT": 0.177,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(1076.8106983140926)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.40893371042459314)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(494.15715728726246)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0018696283187867953)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        -0.32333317432324193
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.15816438531073498
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(0.04512950369492971)
