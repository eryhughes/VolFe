# Tests the degassing function works as expected

import VolFe as vf
import pandas as pd
import pytest


# tests to complete
# 2a sulf sat
# 2b closed-system wtg
# 2d closed regas
# 2d open regas

options = vf.default_models.copy()
options.loc["output csv", "option"] = False


def test_degas_df_default():
    "simple test of calc_gassing function"

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

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(3863.58)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.720353)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1475.998)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(6.72298e-4)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.0872088)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.006361992)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(1.36331e-05)


# def test_degas_df_sat_sulf():
#    "simple test of calc_gassing function with sulfur saturation"

#    my_analysis = {'Sample':'Sari15-04-33',
#           'T_C': 1200., # Temperature in 'C
#           'SiO2': 47.89, # wt%
#           'TiO2': 0.75, # wt%
#           'Al2O3': 16.74, # wt%
#           'FeOT': 9.43, # wt%
#           'MnO': 0.18, # wt%
#           'MgO': 5.92, # wt%
#           'CaO': 11.58, # wt%
#           'Na2O': 2.14, # wt%
#           'K2O': 0.63, # wt%
#           'P2O5': 0.17, # wt%
#           'H2O': 4.17, # wt%
#           'CO2ppm': 1487., # ppm
#           'STppm': 1343.5, # ppm
#           'Xppm': 0., # ppm
#           'Fe3FeT': 0.195}

#    my_analysis = pd.DataFrame(my_analysis, index=[0])

# Choose the options I want to change for the calculation
# - everything else will use the default options
#    my_models = [['sulfur_saturation','True']]

# turn to dataframe with correct column headers and indexes
#    my_models = vf.make_df_and_add_model_defaults(my_models)

#    result = vf.calc_gassing(my_analysis,models=my_models)

#    assert result.loc[0,"P_bar"] == pytest.approx(3848.9041360851006)
#    assert result.loc[0,"fO2_DFMQ"] == pytest.approx(0.7200378919508879)
#    assert result.loc[0,"CO2T_ppmw"] == pytest.approx(1475.9913934875094)
#    assert result.loc[0,"xgS2_mf"] == pytest.approx(0.00038774672125415206)
#    assert result.loc[len(result)-1,'P_bar'] == 1.0
#    assert result.loc[len(result)-1,"fO2_DFMQ"] == pytest.approx(0.11170286344042957)
#    assert result.loc[len(result)-1,"CO2T_ppmw"] == pytest.approx(0.006450832514387365)
#    assert result.loc[len(result)-1,"xgS2_mf"] == pytest.approx(3.2326503829701168e-06)


def test_degas_df_closed_CO2i():
    "simple test of calc_gassing function with closed-system degassing and initial CO2"

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
        "initial_CO2wtpc": 4.0,
    }  # initial CO2 content of the system in wt%

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["bulk_composition", "melt+vapor_initialCO2"], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3863.5831334191394)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.7203532684511362)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1475.9979291041839)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0006722975341518086)
    assert result.loc[len(result) - 1, "P_bar"] == 3.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.09899649250168707)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.33989935002108057
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        5.67343267164686e-05
    )


# def test_degas_df_closed_wtg():
#    "simple test of calc_gassing function with closed-system degassing and wtg"

#    my_analysis = {'Sample':'Sari15-04-33',
#           'T_C': 1200., # Temperature in 'C
#           'SiO2': 47.89, # wt%
#           'TiO2': 0.75, # wt%
#           'Al2O3': 16.74, # wt%
#           'FeOT': 9.43, # wt%
#           'MnO': 0.18, # wt%
#           'MgO': 5.92, # wt%
#           'CaO': 11.58, # wt%
#           'Na2O': 2.14, # wt%
#           'K2O': 0.63, # wt%
#           'P2O5': 0.17, # wt%
#           'H2O': 4.17, # wt%
#           'CO2ppm': 1487., # ppm
#           'STppm': 1343.5, # ppm
#           'Xppm': 0., # ppm
#           'Fe3FeT': 0.195,
#           'wt_g': 3.} # wt% vapor in equilibrium with the melt

#    my_analysis = pd.DataFrame(my_analysis, index=[0])

# choose the options I want - everything else will use the default options
#    my_models = [['bulk_composition','melt+vapor_wtg']]

# turn to dataframe with correct column headers and indexes
#    my_models = vf.make_df_and_add_model_defaults(my_models)

#    result = vf.calc_gassing(my_analysis,models=my_models)

#    assert result.loc[0, "P_bar"] == pytest.approx(3863.5831334191394)
#    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.7203532684511362)
#    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1475.9979291041839)
#    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0006722975341518086)
#    assert result.loc[len(result) - 1, "P_bar"] == 1.0
#    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.7588141360313543)
#    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx()
#    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx()


def test_degas_df_open():
    "simple test of calc_gassing function for open-system degassing"

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
        "Fe3FeT": 0.195,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["gassing_style", "open"], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(199.78358397995052)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.63390049438501)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(49.54358631280869)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00018482167916692897)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.6315427057280116)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == 0.0
    assert result.loc[len(result) - 1, "xgS2_mf"] == 0.0


def test_degas_df_HSO():
    "simple test of calc_gassing function for HSO system"

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
        "Fe3FeT": 0.195,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(406.37228244941497)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.6391885514882905)
    assert result.loc[0, "CO2T_ppmw"] == 0.0
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.008314410934512158)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        -0.09058551042920904
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == 0.0
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(
        3.5233776095813517e-05
    )


def test_degas_df_CSO():
    "simple test of calc_gassing function for CSO system"

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
        "Fe3FeT": 0.195,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(1098.0598069382127)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.6565341610092759)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(495.5956612185359)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0026586642922226352)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
        -0.03611792004490333
    )
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
        0.13291287430638418
    )
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(0.01615449833076651)


def test_degas_df_CHOAr_basalt():
    "simple test of calc_gassing function for CHOAr system in basalt"

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
        "Fe3FeT": 0.195,
    }

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(1460.1133014759528)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.6653926960786709)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(495.7596295369173)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.08571694089111581)
    # assert result.loc[len(result) - 1, "P_bar"] == 70.0
    # assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
    #     0.7115340646298742
    # )
    # assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
    #     0.5743842239770661
    # )
    # assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(
    #     0.00039491054832845393
    # )


def test_degas_df_CHONe_basalt():
    "simple test of calc_gassing function for CHONe system in basalt"

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
        "Fe3FeT": 0.195,
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

    assert result.loc[0, "P_bar"] == pytest.approx(1404.7646337241717)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.6640482770575096)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(495.7479166646544)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.047331318077006045)
    # assert result.loc[len(result) - 1, "P_bar"] == 120.
    # assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
    # 0.6941924132550827)
    # assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
    # 1.3001260013356675)
    # assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(
    # 0.001026819006300385)


def test_degas_df_CHOAr_rhyolite():
    "simple test of calc_gassing function for CHOAr system in rhyolite"

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
        "Fe3FeT": 0.195,
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

    assert result.loc[0, "P_bar"] == pytest.approx(1363.4483591372611)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.6630423937635648)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(495.7391411775598)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.01666896481627928)
    # assert result.loc[len(result) - 1, "P_bar"] == 120.
    # assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(
    #   0.7004837361155927)
    # assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(
    #   0.9697530002782397)
    # assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(
    #   0.0004647839509440657)
