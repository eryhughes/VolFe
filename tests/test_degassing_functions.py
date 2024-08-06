# Tests the degassing function works as expected

import VolFe as vf
import pandas as pd
import pytest
import os


def test_degassing_simple_fails(capsys):
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

    # commented out for now as this takes too long on the CI

    # degas1 = vf.calc_gassing(my_analysis)

    # assert degas1.shape == (386, 177)

    # # This reads the terminal output and checks the failed solver message is present
    # terminal_output = capsys.readouterr()

    # assert "solver failed, calculation aborted" in terminal_output.out

    # os.remove("results_gassing_chemistry.csv")
    # os.remove("results_jacnewton3.csv")

def test_degas_df_default(capsys):
    "simple test of calc_gassing function"
    
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
           'Fe3FeT': 0.195}

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result = vf.calc_gassing(my_analysis)

    assert result.loc[0,"P_bar"] == pytest.approx(3863.5831334191394)
    assert result.loc[0,"fO2_DFMQ"] == pytest.approx(0.7203532684511362)
    assert result.loc[0,"CO2T_ppmw"] == pytest.approx(1475.9979291041839)
    assert result.loc[0,"xgS2_mf"] == pytest.approx(0.0006722975341518086)
    assert result.loc[len(result)-1,'P_bar'] == 1.0
    assert result.loc[len(result)-1,"fO2_DFMQ"] == pytest.approx(-0.08720880727594427)
    assert result.loc[len(result)-1,"CO2T_ppmw"] == pytest.approx(0.006361991861169213)
    assert result.loc[len(result)-1,"xgS2_mf"] == pytest.approx(1.3633089495821218e-05)   