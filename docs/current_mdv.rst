=========================
Options
=========================

Various parameters in VolFe can be calculated using different models - these include parameters like solubility constants, equilibrium constants, fugacity coefficients, etc. 
For solubility constants, oxygen fugacity to Fe3+/FeT relationships, and sulfur saturation conditions, this is how VolFe includes the effect of melt compositions on calculations.
    
There are some models already available in VolFe, which are detailed below. 
Others can be added as they become available (see Worked Example in :doc:`add your own <add_your_own>`) - let us know if you have a new model to be added! 
Additionally, there are various options for how the calculations are done in VolFe. 
This varies from which species are treated as insoluble in the melt to whether a csv of the results is created at the end of the calcualtion.

- :doc:`Equilibrium constants <mdv/equilibrium_constants>`

- :doc:`Fugacity coefficients <mdv/fugacity_coefficients>`

- :doc:`Solubility constants <mdv/solubility_constants>`

- :doc:`Oxygen fugacity and Fe3+/FeT <mdv/oxygen_fugacity>`

- :doc:`Sulfide/sulfate content at sulfide/anhydrite saturation <mdv/saturation_conditions>`

- :doc:`Calculation options <mdv/calculation_options>`