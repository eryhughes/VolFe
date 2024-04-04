===============================================
1. Calculating the pressure of vapor saturation
===============================================

The volatile content of a melt in equilibrium with a vapor can be used as a barometer because the sum of the partial pressures of all the vapor species must equal the total pressure (REF). This is often applied to melt inclusions to calculate magma storage depths (REF) and sub-aqueous matrix glasses to calculate eruption depths (REF).

VolFe can calculate the pressure of vapor saturation (Pvsat), melt speciation, and vapor speciation for a melt of given temperature and melt composition (including volatiles and oxygen fugacity). This calculation is described in detail in Hughes et al. (2024).

In this example we'll show you how to run this calculation for: 

:doc:`Example 1a <Examples/1a. pvsat 1MI_df>`: One analysis entered as a dataframe using default options. 

:doc:`Example 1b <Examples/1b. pvsat csv>`: Analyses in a csv file using default options. 

:doc:`Example 1c <Examples/1c. pvsat user_opts>`: Analyses in a csv file using user specified options.

Hughes, E.C., Liggins, P., Saper, L. and M. Stolper, E., 2024. The efects of oxygen fugacity and sulfur on the pressure of vapor-saturation of magma. American Mineralogist, 109(3), pp.422-438. https://doi.org/10.2138/am-2022-8739 