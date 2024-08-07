# Tests the calc_pvsat function works as expected

import VolFe as vf
import pandas as pd
import pytest

def test_S2fO2_df_pvsat():
    "simple test of calc_melt_S_oxybarometer function at Pvsat"

    # Define the melt composition and T as a dictionary.
    my_analysis = {'Sample':'Sari15-04-34',
           'T_C': 1200., # Temperature in 'C
           'SiO2': 46.94, # wt%
           'TiO2': 0.65, # wt%
           'Al2O3': 16.0, # wt%
           'FeOT': 8.92, # wt%
           'MnO': 0.15, # wt%
           'MgO': 7.3, # wt%
           'CaO': 13.94, # wt%
           'Na2O': 1.61, # wt%
           'K2O': 0.26, # wt%
           'P2O5': 0.07, # wt%
           'H2O': 3.83, # wt%
           'CO2ppm': 1109., # ppm
           'STppm': 1614.12, # ppm
           'Xppm': 0.} # ppm
													
    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    result1 = vf.calc_melt_S_oxybarometer(my_analysis)

    assert result1.loc[0, "P (bar) sulf"] == pytest.approx(3209.621469)
    assert result1.loc[0, "DFMQ-sulfide"] == pytest.approx(0.905499)
    assert result1.loc[0, "P (bar) anh"] == pytest.approx('')
    assert result1.loc[0, "DFMQ-sulfate"] == pytest.approx('')

def test_S2fO2_df_P():
    "simple test of calc_melt_S_oxybarometer function at given P"

    # Define the melt composition, P, and T as a dictionary.
    my_analysis = {'Sample':'Sari15-04-34',
           'T_C': 1200., # Temperature in 'C
           'SiO2': 46.94, # wt%
           'TiO2': 0.65, # wt%
           'Al2O3': 16.0, # wt%
           'FeOT': 8.92, # wt%
           'MnO': 0.15, # wt%
           'MgO': 7.3, # wt%
           'CaO': 13.94, # wt%
           'Na2O': 1.61, # wt%
           'K2O': 0.26, # wt%
           'P2O5': 0.07, # wt%
           'H2O': 3.83, # wt%
           'CO2ppm': 1109., # ppm
           'STppm': 1614.12, # ppm
           'Xppm': 0., # ppm
           'P_bar': 1000., # bar
           'Fe3FeT': 0.1} 
													
    # Turn the dictionary into a pandas dataframe, setting the index to 0.
    my_analysis = pd.DataFrame(my_analysis, index=[0])

    # runs the calculation
    result = vf.calc_melt_S_oxybarometer(my_analysis)

    assert result.loc[0, "P (bar) sulf"] == pytest.approx(1000.0)
    assert result.loc[0, "DFMQ-sulfide"] == pytest.approx(1.158481)
    assert result.loc[0, "P (bar) anh"] == pytest.approx(1000.0)
    assert result.loc[0, "DFMQ-sulfate"] == pytest.approx('')