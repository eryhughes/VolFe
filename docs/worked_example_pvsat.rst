=============================================
1. Calculating the pressure of vapor saturation
=============================================

The volatile content of a melt in equilibrium with a vapor can be used as a barometer because the sum of the partial pressures of all the vapor species must equal the total pressure (REF). This is often applied to melt inclusions to calculate magma storage depths (REF) and sub-aqueous matrix glasses to calculate eruption depths (REF).

VolFe can calculate the pressure of vapor saturation (Pvsat), melt speciation, and vapor speciation for a melt of given temperature and melt composition (including volatiles and oxygen fugacity). 

In this example we'll show you how to run this calculation for: 

a) a single analysis entered as a dataframe using default options: :doc:`Example 1a <Examples/1a. pvsat 1MI_df>`

b) a single analysis from a csv file using default options

c) all analyses in a csv file using default options

d) a set of consequative rows from a csv file using options

e) all analyses in a csv file using one set of user specified options for all analyses

f) all analyses in a csv file using different model options for each analysis