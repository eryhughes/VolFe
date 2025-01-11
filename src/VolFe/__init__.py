__author__ = "Ery Hughes"
__version__ = "0.3"

# functions for calculating properties of the melt and gas
from .melt_gas import *

# functions for variables that can be model depedent (e.g., solubility constants,
# fugacity coefficients)
from .model_dependent_variables import *

# functions for calculating equilibrium speciation and concentration between melt ± gas
# at given P and T
from .equilibrium_equations import *
from .differential_equations import *

# functions to calculate equilibrium isotope fractionation given melt ± gas composition
from .batch_calculations import *

# functions to run a calculations
from .calculations import *

# functions to run calculations in batch mode
from .batch_calculations import *
