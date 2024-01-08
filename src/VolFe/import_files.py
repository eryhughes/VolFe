# import_data.py

import pandas as pd

# loads species.csv, models.csv, and user setup file
def import_files(filename):
    species = pd.read_csv("species.csv", index_col = [0]) # attributes of the different species used in the system
    models = pd.read_csv("models.csv", index_col = [0]) # model options
    setup = pd.read_csv(filename) # csv for initial conditions of the system - typically inputs.csv
    return species, models, setup