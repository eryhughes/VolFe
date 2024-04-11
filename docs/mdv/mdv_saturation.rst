===================================================================================
Sulfide/sulfate content at sulfide/anhydrite saturation
===================================================================================

Below are the different models for the sulfide content an sulfide saturation and sulate content at anhydrite currently in VolFe, which can also be viewed in the :doc:`saturation conditions notebook <SatCond>`.  
The sulfide is assumed to be pure FeS unless specified in the input dataframe using "sulf_XFe", "sulf_XCu", and "sulf_XNi" (mole fractions of FeS, CuS, and NiS in the sulfide, respectively)
The model type is in **bold** (i.e., what option will effect this variable) and the options of models to use are in *italics* (default options are indicated).


**SCSS**

Model for parameterisation of the sulfide content at sulfide saturation (S2-CSS).
        
- *ONeill21hyd*: Eq. (10.34, 10.43, 10.45, 10.46, 10.49) from O'Neill (2021) [default]

- *ONeill21*: Eq. (10.34, 10.43, 10.45, 10.46) excluding water dilution from O'Neill (2021)

- *ONeill21dil*: Eq. (10.34, 10.43, 10.45, 10.46) including water dilution from O'Neill (2021)

- *Liu07*: Eq. (9) in Liu et al. (2007)

- *Fortin15*: Eq. (7) Fortin et al. (2015)

- *Liu21*: Eq. (2) Liu et al. (2021)

- *Fortin15_pss*: Fortin et al. (2015) using PySulfSat by Wieser & Gleeson (2023)

- *Liu21_pss*: Liu et al. (2021) using PySulfSat by Wieser & Gleeson (2023)

- *ONeill22_pss*: O'Neill & Mavrogenes (2022) using PySulfSat by Wieser & Gleeson (2023)

- *ONeill21_pss*: O'Neill (2021) using PySulfSat by Wieser & Gleeson (2023)

- *Smythe17_pss*: Smythe et al. (2017) using PySulfSat by Wieser & Gleeson (2023)


**SCAS**

Model for parameterisation of the sulfate content at anhydrite saturation (S6+CAS).

- *Zajacz19*: Eq. (8-14) Zajacz & Tsay (2019)

- *Chowdhury19*: Eq. (8) using Table 5 in Chowdhury & Dasgupta (2019)

- *Liu23*: Eq. (4) Liu et al. (2023)

- *Chowdhury19_pss*: Chowdhury & Dasgupta (2019) using PySulfSat by Wieser and Gleeson (2023)

- *Zajacz19_pss*: Zajacz and Tsay (2019) using PySulfSat by Wieser and Gleeson (2023)