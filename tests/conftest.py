# In your conftest.py or test file
import os


def pytest_sessionfinish():
    # List of files to remove at the end of the session
    files_to_remove = [
        "results_jacnewton2.csv",
        "results_jacnewton3.csv",
        "results_newtraph.csv",
        "capacities.csv",
        "fO2_range_from_S.csv",
        "results_fugacity_coefficients.csv",
        "results_isobars.csv",
        "results_pure_solubility.csv",
        "results_saturation_pressures.csv",
    ]

    # Cleanup action
    for file in files_to_remove:
        if os.path.exists(file):
            os.remove(file)
