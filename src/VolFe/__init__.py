__author__ = 'Ery Hughes'

import pandas as pd
from datetime import date
import gmpy2 as gp
import numpy as np
import datetime
import math as math
from scipy import optimize
import densityx as dx
import PySulfSat as ss

# functions for calculating properties of the melt and gas
from VolFe.melt_gas import *
# functions for variables that can be model depedent (e.g., solubility constants, fugacity coefficients)
from VolFe.model_dependent_variables import *
# functions for calculating equilibrium speciation and concentration between melt ± gas at given P and T
from VolFe.equilibrium_equations import *
from VolFe.differential_equations import *
# functions to calculate equilibrium isotope fractionation given melt ± gas composition
from VolFe.batch_calculations import *
# functions to run a calculations
from VolFe.calculations import *
# functions to run calculations in batch mode
from VolFe.batch_calculations import *

# version
from ._version import __version__