# Tests the calc_pvsat function works as expected

import VolFe as vf
import pandas as pd
import pytest
import os


def test_pvsat_df_FeOT_Fe3FeT(capsys):
    "simple test of calc_pvsat function using FeOT and Fe3FeT"
    
    my_analysis = {'Sample':'TN273-01D-01-01',
           'T_C': 1200., # Temperature in 'C
           'SiO2': 56.98, # wt%
           'TiO2': 1.66, # wt%
           'Al2O3': 15.52, # wt%
           'FeOT': 9.47, # wt%
           'MnO': 0.24, # wt%
           'MgO': 2.96, # wt%
           'CaO': 6.49, # wt%
           'Na2O': 4.06, # wt%
           'K2O': 0.38, # wt%
           'P2O5': 0.22, # wt%
           'H2O': 1.88, # wt%
           'CO2ppm': 13., # ppm
           'STppm': 362.83, # ppm
           'Xppm': 0., # ppm
           'Fe3FeT': 0.171} # mole or weight fraction (they're the same)

    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0,"P_bar"] == pytest.approx(337.80934489140867)
    assert result1.loc[0,"fO2_DFMQ"] == pytest.approx(0.5346132659807692)

def test_pvsat_df_Fe2O3T_DFMQ(capsys):
    "simple test of calc_pvsat function using Fe2O3T and DFMQ"

    my_analysis = {'Sample':'TN273-01D-01-01',
           'T_C': 1200.,
           'SiO2': 56.98,
           'TiO2': 1.66,
           'Al2O3': 15.52,
           'Fe2O3T': 9.47, # Fe2O3T instead of FeOT
           'MnO': 0.24,
           'MgO': 2.96,
           'CaO': 6.49,
           'Na2O': 4.06,
           'K2O': 0.38,
           'P2O5': 0.22,
           'H2O': 1.88,
           'CO2ppm': 13.,
           'STppm': 362.83,
           'Xppm': 0.,
           'DFMQ': 0.} # DFMQ instead of Fe3+/FeT

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0,"P_bar"] == pytest.approx(317.2180305487919)
    assert result1.loc[0,"Fe3+/FeT"] == pytest.approx(0.13998995667087746)


def test_pvsat_df_FeO_Fe2O3(capsys):
    "simple test of calc_pvsat function using FeO and Fe2O3"

    my_analysis = {'Sample':'TN273-01D-01-01',
           'T_C': 1200.,
           'SiO2': 56.98,
           'TiO2': 1.66,
           'Al2O3': 15.52,
           'FeO': 8.03, # FeO and Fe2O3 to contraint total Fe and fO2
           'Fe2O3': 2.57, # ^
           'MnO': 0.24,
           'MgO': 2.96,
           'CaO': 6.49,
           'Na2O': 4.06,
           'K2O': 0.38,
           'P2O5': 0.22,
           'H2O': 1.88,
           'CO2ppm': 13.,
           'STppm': 362.83,
           'Xppm': 0.} 

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0,"P_bar"] == pytest.approx(436.4296485406602)
    assert result1.loc[0,"fO2_DFMQ"] == pytest.approx(1.227099405389815)

def test_pvsat_df_FeOT_DNNO(capsys):
    "simple test of calc_pvsat function using FeOT and DNNO"

    my_analysis = {'Sample':'TN273-01D-01-01',
           'T_C': 1200.,
           'SiO2': 56.98,
           'TiO2': 1.66,
           'Al2O3': 15.52,
           'FeOT': 9.47, ###
           'MnO': 0.24,
           'MgO': 2.96,
           'CaO': 6.49,
           'Na2O': 4.06,
           'K2O': 0.38,
           'P2O5': 0.22,
           'H2O': 1.88,
           'CO2ppm': 13.,
           'STppm': 362.83,
           'Xppm': 0.,
           'DNNO': 1.} ### 

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0,"P_bar"] == pytest.approx(451.9452953846255)
    assert result1.loc[0,"Fe3+/FeT"] == pytest.approx(0.26742007392545164)

def test_pvsat_df_FeOT_S6ST(capsys):
    "simple test of calc_pvsat function using FeOT and S6+/ST"

    my_analysis = {'Sample':'TN273-01D-01-01',
           'T_C': 1200.,
           'SiO2': 56.98,
           'TiO2': 1.66,
           'Al2O3': 15.52,
           'FeOT': 9.47, 
           'MnO': 0.24,
           'MgO': 2.96,
           'CaO': 6.49,
           'Na2O': 4.06,
           'K2O': 0.38,
           'P2O5': 0.22,
           'H2O': 1.88,
           'CO2ppm': 13.,
           'STppm': 362.83,
           'Xppm': 0.,
           'S6ST': 0.23} ###

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0,"P_bar"] == pytest.approx(395.6435433141458)
    assert result1.loc[0,"fO2_DFMQ"] == pytest.approx(1.0129274022003445)

def test_pvsat_df_FeOT_DNNO_Fe3FeT(capsys):
    "simple test of calc_pvsat function using FeOT, DNNO, and Fe3+/FeT"

    my_analysis = {'Sample':'TN273-01D-01-01',
           'T_C': 1200.,
           'SiO2': 56.98,
           'TiO2': 1.66,
           'Al2O3': 15.52,
           'FeOT': 9.47, ###
           'MnO': 0.24,
           'MgO': 2.96,
           'CaO': 6.49,
           'Na2O': 4.06,
           'K2O': 0.38,
           'P2O5': 0.22,
           'H2O': 1.88,
           'CO2ppm': 13.,
           'STppm': 362.83,
           'Xppm': 0.,
           'DNNO': 1.,
           'Fe3FeT':0.171} ###

    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result1 = vf.calc_Pvsat(my_analysis)

    assert result1.loc[0,"P_bar"] == pytest.approx(337.80934489140867)
    assert result1.loc[0,"fO2_DFMQ"] == pytest.approx(0.5346132659807692)


