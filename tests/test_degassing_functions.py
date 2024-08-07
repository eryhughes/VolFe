# Tests the degassing function works as expected

import VolFe as vf
import pandas as pd
import pytest


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

    result = vf.calc_gassing(my_analysis)

    assert result.loc[0, "P_bar"] == pytest.approx(3863.58)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.720353)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1475.998)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(6.72298e-4)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.0872088)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.006361992)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(1.36331e-05)


#def test_degas_df_sat_sulf():
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

    # Choose the options I want to change for the calculation - everything else will use the default options
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
           'Fe3FeT': 0.195,
           'initial_CO2wtpc': 4.} # initial CO2 content of the system in wt%

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # choose the options I want - everything else will use the default options
    my_models = [['bulk_composition','melt+vapor_initialCO2']]

    # turn to dataframe with correct column headers and indexes    
    my_models = vf.make_df_and_add_model_defaults(my_models)

    result = vf.calc_gassing(my_analysis,models=my_models)

    assert result.loc[0, "P_bar"] == pytest.approx(3863.5831334191394)
    assert result.loc[0, "fO2_DFMQ"] == pytest.approx(0.7203532684511362)
    assert result.loc[0, "CO2T_ppmw"] == pytest.approx(1475.9979291041839)
    assert result.loc[0, "xgS2_mf"] == pytest.approx(0.0006722975341518086)
    assert result.loc[len(result) - 1, "P_bar"] == 1.0
    assert result.loc[len(result) - 1, "fO2_DFMQ"] == pytest.approx(-0.1334627142973961)
    assert result.loc[len(result) - 1, "CO2T_ppmw"] == pytest.approx(0.11187099769675689)
    assert result.loc[len(result) - 1, "xgS2_mf"] == pytest.approx(2.3080542012091533e-05)