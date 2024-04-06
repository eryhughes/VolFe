===================================================================================
Calculating the pressure of vapor saturation (Pvsat)
===================================================================================

The dissolved volatile content of a melt in equilibrium with a vapor can be used as a barometer because the sum of the partial pressures of all the vapor species must equal the total pressure (REF). 
This is the pressure of vapor saturation (Pvsat) barometer and is often applied to melt inclusions to calculate magma storage depths (REF) and sub-aqueous matrix glasses to calculate eruption depths (REF). 

VolFe calculates Pvsat, the melt speciation, and the vapor composition for a melt of given temperature and melt composition (including volatiles and fO2). 

This calculation was outlined in detail Hughes et al. (2024) (schematic in figure below, Hughes et al. in prep).

:image:`<figures/pvsatcalc.png>` :width: 800

In this example we'll show you how to run this calculation for: 

- :doc:`Example 1a <Examples/1a. pvsat 1MI_df>`: One analysis entered as a dataframe using default options. 

- :doc:`Example 1b <Examples/1b. pvsat csv>`: Analyses in a csv file using default options. 

- :doc:`Example 1c <Examples/1c. pvsat user_opt>`: Analyses in a csv file using user specified options.