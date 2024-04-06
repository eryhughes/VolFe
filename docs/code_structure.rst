==========
Code structure
==========

This is a brief overview of how the VolFe code is structured and broadly where to find different types of functions in VolFe.
You don't need to understand this to run VolFe, but popped it here in case people are interested! 
For instance, you could use VolFe functions in your code if they're helpful (particuarly model dependent variables).


Melt gas
---
"melt_gas.py" hosts functions for calculating properties of the melt and gas a set of conditions (e.g., *P* and *T*) - so these are the forward calculations. Broadly, this includes things like: 

- fugacities (f_i), partial pressures (p_i), and mole fraction in vapor (xg_i) of vapor species given *P*, *T* and volatile concentration in the melt.

- coverting melt composition between mass and mole fractions, based on cations, oxides, single oxygen, etc.

- calculating melt speciation for oxidised carbon, oxidised hydrogen, sulfur, and iron given *P*, *T*, dissovled volatile concentration, and *f*O2.


Model dependent variables
---
"model_dependent_variables.py" hosts the functions for calculating various model dependent variables given a set of (e.g., *P* and *T*) - so these are the forward calculations. 
More information can be found on the :doc:`current model dependent variables <current_mdv>` page, but the types of parameters calculated and stored here are:

- element and species constants (e.g., atomic and molecular mass)

- functions to build the models option dataframe, including the defaul options

- functions to calculate solubility, speciation, and equilibrium constants; fugacity coefficients; SCSS and SCAS, and melt density

- functions to convert between fO2 and Fe3+/FeT and fO2 buffers


Equilibrium equations and Differential equations
---

“equilibrium_equations.py” contains the functions to calculate the equilibrium concentration and speciation of volatiles in melt ± vapor and the proportions of melt and vapor given *P*, *T*, and the bulk composition of the system.
These are iterative calculations, which are solved using Jacobian matrix/Newton-Raphson approach (solvers and initial guesses are also in this file).
The associated differential equations required for the solver are contained in “differential_equations.py”, where the symbolic differentiation was done using SymPy (Meurer et al., 2017).
It is a long file because there are separate functions for different combinations of species in the melt and vapor (e.g., only carbon species in the system; the full COHS system; only H2O-CO2 in the melt and vapor, etc.)

Calculations
---
“calculations.py” contains broadly two types of functions. 

Firstly, there are functions for calculating properties about the system, such as:

- bulk composition given melt and vapor compositions

- is sulfur or graphite saturated at a given set of conditions

- Monte Carlo simulation of compositional errors

Secondly, its calculation types that don't necessarily have independent variables of *P*, *T*, and bulk composition (as in "equilibrium_equations.py"), such as:

- *T*, melt composition (inc. fO2 estimate), vapor-saturated: calculate Pvsat and the melt/vapor composition and speciation.

- *T*, melt composition, vapor+sulfide+anhydrite-saturated, ±P: calculate fO2 ± P and the melt/vapor composition and speciation.

- *T*, melt composition, vapor+(sulfide or anhydrite)-saturated, ± P: calculate fO2 ± P and the melt/vapor composition and speciation.

- *T*, *P*, *f*O2, CO2-eq, H2O-eq, vapor±(sulfide and/or anhydrite)-saturation = calculate how much sulfur can dissolve in the melt.

- *T*, *P*, volatile-free composition, assuming only H2O and CO2: calculate concentrations of H2O and CO2 in the melt (i.e., isobars).


Batch calculations
---
"batch_calculations.py" joins calculations together to enable them to be run for multiple conditions. 
This might be to create a degassing calculations (i.e., run over multiple P) and calculate Pvsat for a whole spreadsheet of calculations. 
More details can be found in the Worked Examples.