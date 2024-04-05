===============================================================
3. Calculating the oxygen fugacity from the melt sulfur content
===============================================================

Oxygen fugacity is a key thermodynamic parameter to estimate in magmatic systems because of its effects on the chemical and physical properties of the melt (REF). 
In certain circumstances, the sulfur content of the melt can be used to place bounds on the oxygen fugacity based on sulfide and anhydrite saturation (e.g., Beerman et al., 2011; Muth and Wallace, 2022; Hughes et al., 2023). 

This calculation was outlined in detail in Section “Using wmST as an oxybarometer” in Hughes et al. (2023) (schematic in figure below, Hughes et al. in prep).

:image:`<figures/SfO2calc.png>` :width: 800

In this example we'll show you how to run this calculation for: 

:doc:`Example 3a <Examples/3a. SfO2 1MI_df>`: One analysis entered as a dataframe using default options. 

:doc:`Example 3b <Examples/3b. SfO2 csv>`: Analyses in a csv file using default options. 

:doc:`Example 3c <Examples/3c. SfO2 user_opt>`: Analyses in a csv file using user specified options.

Beermann, O., Botcharnikov, R.E., Holtz, F., Diedrich, O. and Nowak, M., 2011. Temperature dependence of sulfide and sulfate solubility in olivine-saturated basaltic magmas. Geochimica et Cosmochimica Acta, 75(23), pp.7612-7631. https://doi.org/10.1016/j.gca.2011.09.024 

Hughes, E.C., Saper, L.M., Liggins, P., O'Neill, H.S.C. and Stolper, E.M., 2023. The sulfur solubility minimum and maximum in silicate melt. Journal of the Geological Society, 180(3), pp.jgs2021-125. https://doi.org/10.1144/jgs2021-125 

Muth, M.J. and Wallace, P.J., 2022. Sulfur recycling in subduction zones and the oxygen fugacity of mafic arc magmas. Earth and Planetary Science Letters, 599, p.117836. https://doi.org/10.1016/j.epsl.2022.117836 
