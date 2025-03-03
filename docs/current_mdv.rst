=========================
Options
=========================

Various parameters in VolFe can be calculated using different models - these include parameters like solubility functions, equilibrium constants, fugacity coefficients, etc.
For solubility functions, oxygen fugacity to Fe3+/FeT relationships, and sulfur saturation conditions, this is how VolFe includes the effect of melt compositions on calculations.

There are some models already available in VolFe, which are detailed below.
Others can be added as they become available (see Worked Example in :doc:`add your own <add_your_own>`) - let us know if you have a new model to be added!

Additionally, there are various options for how the calculations are done in VolFe.
This varies from which species are treated as insoluble in the melt to whether a csv of the results is created at the end of the calculation.

To see the different options available, click on the "make_df_and_add_model_defaults()" under "Model dependent variables" in the API Refence section on the left-hand ribbon.

Here are some hints to the naming convention used for functions, where specific model options for these model dependent variables can be found:

- Equilibrium constants: functions starting with K (e.g., KCOHg() is the function for calculating the equilibrium constants for CH4 + 2O2 = CO2 + 2H2O).

- Fugacity coefficients: functions starting with y_ (e.g., y_CH4() is the function for calculating the fugacity coefficient for CH4).

- Solubility functions: functions starting with C_ (e.g., C_CH4() is the function for calculating the solubility function for CH4).

- Oxygen fugacity and Fe3+/FeT: functions are FMQ(), NNO(), fO22Fe3FeT().

- Sulfide/sulfate content at sulfide/anhydrite saturation: functions are SCSS() and SCAS().