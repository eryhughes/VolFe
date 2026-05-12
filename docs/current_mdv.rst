=========================
Options
=========================

Various parameters in VolFe can be calculated using different models - these include parameters like solubility functions, equilibrium constants, fugacity coefficients, isotope fractionation factors, etc.
For solubility functions, oxygen fugacity to Fe\ :sup:`3+`/Fe\ :sub:`T` relationships, and sulfur saturation conditions, this is how VolFe includes the effect of melt compositions on calculations.

There are lots of models already available in VolFe, which can be found in the API reference section.
Others can be added as they become available (see :doc:`Add your own <add_your_own>` in the Worked Examples section) - let us know if you have a new model to be added!

Additionally, there are various options for how the calculations are done in VolFe.
This varies from which species are treated as insoluble in the melt to whether a csv of the results is created at the end of the calculation.

To see the different options available (both calculation options and models for different parameters), click on ``make_df_and_add_model_defaults()`` under "Model dependent variables" in the API Reference section on the left-hand ribbon.
For calculation options, these are listed here.
For models for different parameters, it will point you to the function to check that contains the model options.

If you're using VolFe, you can use ``help()`` to see the model options for each parameter by searching for the relevant function - at the bottom there will be a section called 'model options for...', which includes references.

Here are some hints to the naming convention used for functions, where specific model options for these model dependent variables can be found:

- Equilibrium constants: functions starting with ``K`` (e.g., ``KCOHg()`` is the function for calculating the equilibrium constants for CH\ :sub:`4` + 2O\ :sub:`2` = CO\ :sub:`2` + 2H\ :sub:`2` O - ``help(volfe.KCOHg)``).

- Fugacity coefficients: functions starting with ``y_`` (e.g., ``y_CH4()`` is the function for calculating the fugacity coefficient for CH\ :sub:`4` - ``help(volfe.y_CH4)``).

- Solubility functions: functions starting with ``C_`` (e.g., ``C_CH4()`` is the function for calculating the solubility function for CH\ :sub:`4` - ``help(volfe.C_CH4)``).

- Oxygen fugacity and Fe\ :sup:`3+`/Fe\ :sub:`T`: functions are ``FMQ()``, ``NNO()``, ``fO22Fe3FeT()``, and ``f_O2()``.

- Sulfide/sulfate content at sulfide/anhydrite saturation: functions are ``SCSS()`` and ``SCAS()``.

- Isotope fractionation factors: functions are ``alpha_`` (e.g., ``alpha_S_SO2v_S6pm`` is the function for calculating the isotopic fractionation factor between SO\ :sub:`2` in the vapor and S\ :sup:`6+` in the melt - ``help(volfe.alpha_S_SO2v_S6pm)``)