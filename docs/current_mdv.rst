=========================
Options
=========================

Various parameters in VolFe can be calculated using different models - these include parameters like solubility functions, equilibrium constants, fugacity coefficients, etc. 
For solubility functions, oxygen fugacity to Fe3+/FeT relationships, and sulfur saturation conditions, this is how VolFe includes the effect of melt compositions on calculations.
    
There are some models already available in VolFe, which are detailed below. 
Others can be added as they become available (see Worked Example in :doc:`add your own <add_your_own>`) - let us know if you have a new model to be added! 

Additionally, there are various options for how the calculations are done in VolFe. 
This varies from which species are treated as insoluble in the melt to whether a csv of the results is created at the end of the calculation.

To see which functions calculate the different model dependent variables and show the different model options currently available in VolFe, 
click on the "Model dependent variables" heading in the "API reference" section on the left-hand ribbon. The click on the required function name, which are detailed below.


Equilibrium constants
-------
- KCOHg(): equilibrium constant for CH4 + 2O2 = CO2 + 2H2O

- KCOg(): equilibrium constant for CO + 0.5O2 = CO

- KCOs(): equilibrium constant for Cgraphite + O2 = CO2

- KHOSg(): equilibrium constant for 0.5S2 + H2O = H2S + 0.5O2

- KHOg(): equilibrium constant for H2 + 0.5O2 = H2O

- KOCSg(): equilibrium constant for 2CO2 + OCS = 3CO + SO2

- KOSg(): equilibrium constant for 0.5S2 + O2 = SO2

- KOSg2(): equilibrium constant for 0.5S2 + 1.5O2 = SO3


- :doc:`Equilibrium constants <mdv/equilibrium_constants>`

- :doc:`Fugacity coefficients <mdv/fugacity_coefficients>`

- :doc:`Solubility functions <mdv/solubility_constants>`

- :doc:`Oxygen fugacity and Fe3+/FeT <mdv/oxygen_fugacity>`

- :doc:`Sulfide/sulfate content at sulfide/anhydrite saturation <mdv/saturation_conditions>`

- :doc:`Calculation options <mdv/calculation_options>`