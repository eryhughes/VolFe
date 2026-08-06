# Tests that the correct warnings are raised when the user tries to do
# something suspicious.

import pytest
import VolFe as vf
import pandas as pd
import warnings

ANALYSIS_ERROR = {
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
    "Xppm": 0.0,
    "DNNO": 1.0,  # fO2 option
    "Fe3FeT": 0.171,
}  # fO2 option


def test_too_many_iron_speciation_options_warning():
    """
    Checks a warning is raised when the user tries to specify more than one way
    to infer iron speciation.
    """
    my_analysis = pd.DataFrame(ANALYSIS_ERROR, index=[0])
    with pytest.warns(
        UserWarning, match="you entered more than one way to infer iron speciation"
    ):
        vf.calc_Pvsat(my_analysis)


ANALYSIS_RUNTIME_WARN = {
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


# def test_warnings_still_raise():
#    my_analysis = pd.DataFrame(ANALYSIS_RUNTIME_WARN, index=[0])
#    models = [["bulk_composition", "melt+vapor_initialCO2"]]
#    my_models = vf.make_df_and_add_model_defaults(models)
#
#    with pytest.warns(RuntimeWarning):
#        vf.calc_gassing(my_analysis, models=my_models, suppress_warnings=False)


def test_warnings_suppressed():
    my_analysis = pd.DataFrame(ANALYSIS_RUNTIME_WARN, index=[0])
    models = [["bulk_composition", "melt+vapor_initialCO2"]]
    my_models = vf.make_df_and_add_model_defaults(models)

    with warnings.catch_warnings(record=True) as record:
        vf.calc_gassing(my_analysis, models=my_models)
    assert not record, "Warnings were not supressed"
