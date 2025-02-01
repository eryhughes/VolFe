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
        "Fe3FeT": 0.177}

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(3847.3361533976154)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47292836384512427)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4015295658357)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00048030277135318553)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.3400085464233662)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.0062599242785629046)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(4.1326131294704035e-05)


def test_degas_df_sat_sulf():
    "simple test of calc_gassing function with sulfur saturation using example 2a"

    my_analysis = {'Sample':'Sari15-04-33',
           'T_C': 1200., # Temperature in 'C
           'SiO2': 47.89, # wt%
           'TiO2': 0.75, # wt%
           'Al2O3': 16.74, # wt%
           'FeOT': 9.43, # wt%
           'MnO': 0.18, # wt%
           'MgO': 5.92, # wt%
           'CaO': 11.58, # wt%
           'Na2O': 2.14, # wt%
           'K2O': 0.63, # wt%
           'P2O5': 0.17, # wt%
           'H2O': 4.17, # wt%
           'CO2ppm': 1487., # ppm
           'STppm': 1343.5, # ppm
           'Xppm': 0., # ppm
           'Fe3FeT': 0.177}

    my_analysis = pd.DataFrame(my_analysis, index=[0])

# Choose the options I want to change for the calculation
# - everything else will use the default options
    my_models = [['sulfur_saturation','True']]

# turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis,models=my_models)

    assert result.loc[0,"P_bar"] == pytest.approx(3834.8735665073505)
    assert result.loc[0,"fO2_DFMQ"] == pytest.approx(0.4726603883083058)
    assert result.loc[0,"CO2T_ppmw"] == pytest.approx(1472.3951153749524)
    assert result.loc[0,"xgS2_mf"] == pytest.approx(0.00021302858508563988)
    assert result.loc[len(result)-1,'P_bar'] == 1.0
    assert result.loc[len(result)-1,"fO2_DFMQ"] == pytest.approx(-0.05880714327033232)
    assert result.loc[len(result)-1,"CO2T_ppmw"] == pytest.approx(0.006398173687476844)
    assert result.loc[len(result)-1,"xgS2_mf"] == pytest.approx(5.39915910546583e-06)


def test_degas_df_closed_CO2i():
    "simple test of calc_gassing function with closed-system degassing and initial CO2 using example 2b"

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
        "initial_CO2wtpc": 4.0}  # initial CO2 content of the system in wt%

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["bulk_composition", "melt+vapor_initialCO2"], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3847.3361533976154)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47292836384512427)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4015295658357)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00048030277135318553)
    assert result.loc[len(result) - 1, "P_bar"] == 3.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.28872864539710363)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.33668806961068326)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(9.77474064866415e-05)


def test_degas_df_closed_wtg():
    "simple test of calc_gassing function with closed-system degassing and wtg using example 2b"

    my_analysis = {'Sample':'Sari15-04-33',
           'T_C': 1200., # Temperature in 'C
           'SiO2': 47.89, # wt%
           'TiO2': 0.75, # wt%
           'Al2O3': 16.74, # wt%
           'FeOT': 9.43, # wt%
           'MnO': 0.18, # wt%
           'MgO': 5.92, # wt%
           'CaO': 11.58, # wt%
           'Na2O': 2.14, # wt%
           'K2O': 0.63, # wt%
           'P2O5': 0.17, # wt%
           'H2O': 4.17, # wt%
           'CO2ppm': 1487., # ppm
           'STppm': 1343.5, # ppm
           'Xppm': 0., # ppm
           'Fe3FeT': 0.177,
           'wt_g': 3.} # wt% vapor in equilibrium with the melt

    my_analysis = pd.DataFrame(my_analysis, index=[0])

# choose the options I want - everything else will use the default options
    my_models = [['bulk_composition','melt+vapor_wtg']]

# turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis,models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3847.3361533976154)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47292836384512427)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4015295658357)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00048030277135318553)
    assert result.loc[len(result) - 1, "P_bar"] == 2.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.31229980033722704)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.16718214175986407)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(7.370194166547196e-05)


def test_degas_df_open():
    "simple test of calc_gassing function for open-system degassing using example 2c"

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
        "Fe3FeT": 0.177}

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [["gassing_style", "open"], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(197.81558695022989)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.38677405655738273)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(49.395162960661466)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00012094116333795436)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.4116067548696618)
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
        'initial_CO2wtpc': 3., # initial CO2 content of the system in wt%
        'final_P':5000.} # bar

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [['gassing_direction','regas'],['bulk_composition','melt+vapor_initialCO2'], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3847.3361533976154)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47292836384512427)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4015295658357)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00048030277135318553)
    assert result.loc[len(result) - 1, "P_bar"] == 5100.
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.4630266222472912)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(2345.1186943946186)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(0.0004195571211495622)

def test_regas_df_open():
    "simple test of calc_gassing function for open-system regassing using example 2d but to 4000 bar to save time"

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
        'initial_CO2wtpc': 3., # initial CO2 content of the system in wt%
        'final_P':4000.} # bar

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [['gassing_direction','regas'],['gassing_style','open'], ["output csv", False]]

    # turn to dataframe with correct column headers and indexes
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis, models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3847.3361533976154)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.47292836384512427)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1472.4015295658357)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.00048030277135318553)
    assert result.loc[len(result) - 1, "P_bar"] == 4001.
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.46976292373420847)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(1579.963744721686)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(0.0004575936690526102)

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

    assert result.loc[0, "P_bar"] == pytest.approx(1463.9612748927782)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.41841020608934265)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(494.3811228660554)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.08549163676895262)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.5594127323997524)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.004903203107918983)
    assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(0.0002394718896877798)


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

    assert result.loc[0, "P_bar"] == pytest.approx(1408.6576220265883)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.41706712820600345)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(494.36565954920235)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.04720051250386283)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.5593157668590552)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.004902089137362658)
    assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(0.00047395291067061073)


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

    assert result.loc[0, "P_bar"] == pytest.approx(1367.3755625959768)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.4160622647228962)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(494.3540745261418)
    assert result.loc[0, "xgX_mf"] == pytest.approx(0.01662109032000306)
    assert result.loc[len(result) - 1, "P_bar"] == 1.
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(0.5592426877679912)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.00490124988736414)
    assert result.loc[len(result) - 1, "xgX_mf"] == pytest.approx(0.0002393756490633622)

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
        "Fe3FeT": 0.177}

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis, models=options)

    assert result.loc[0, "P_bar"] == pytest.approx(380.10070821002193)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.39144299300803187)
    assert result.loc[0, "CO2T_ppmw"] == 0.0
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.005648279548092143)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.33361632492418103)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == 0.0
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(0.00010120701139514083)


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

    assert result.loc[0, "P_bar"] == pytest.approx(1077.2077950994615)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.4089435192488029)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(494.15727361904624)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.001913668221326871)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.2855805216179945)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.1530199730355921)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(0.04026569803310184)