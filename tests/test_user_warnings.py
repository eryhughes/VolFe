# Tests that the correct warnings are raised when the user tries to do
# something suspicious.

import pytest
import VolFe as vf
import pandas as pd

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
