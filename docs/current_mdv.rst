=========================
Options
=========================

Various parameters in VolFe can be calculated using different models - these include parameters like solubility functions, equilibrium constants, fugacity coefficients, isotope fractionation factors, etc.
For solubility functions, oxygen fugacity to Fe\ :sup:`3+`/Fe\ :sub:`T` relationships, and sulfur saturation conditions, this is how VolFe includes the effect of melt compositions on calculations.

There are lots of models already available in VolFe, which can be found in the API reference section.
Others can be added as they become available (see :doc:`Add your own <add_your_own>` in the Worked Examples section) - let us know if you have a new model to be added!

Additionally, there are various options for how the calculations are done in VolFe.
This varies from which species are treated as insoluble in the melt to whether a csv of the results is created at the end of the calculation.

How to set options
------------------

Options are passed to calculation functions as a pandas DataFrame created from a list of ``[parameter_name, option_value]`` pairs using :func:`~VolFe.model_dependent_variables.make_df_and_add_model_defaults`:

.. code-block:: python

   import VolFe as vf

   my_models = [
       ['carbon dioxide', 'Basalt_Dixon97'],
       ['gassing_style', 'open'],
   ]
   models = vf.make_df_and_add_model_defaults(my_models)

Any parameters not specified will use the default option. The ``models`` DataFrame is then passed to calculation functions (e.g., ``vf.calc_Pvsat(setup, models=models)``).

Naming conventions
------------------

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

- Isotope fractionation factors: functions are ``alpha_`` (e.g., ``alpha_S_SO2v_S6pm`` is the function for calculating the isotopic fractionation factor between SO\ :sub:`2` in the vapor and S\ :sup:`6+` in the melt).


All available options
---------------------

Below is a comprehensive list of all model parameters, organized by category.
For each parameter, the default option is shown along with all available alternatives.


Species
^^^^^^^

**COH_species** *(default: yes_H2_CO_CH4_melt)*
   Specifying what COH species are present in the melt and vapor.

   - ``'yes_H2_CO_CH4_melt'`` — Include H2mol (if H present), COmol (if C present), and/or CH4mol (if H and C present) as dissolved melt species.
   - ``'no_H2_CO_CH4_melt'`` — H2, CO, and/or CH4 are insoluble in the melt but they are still present in the vapor (H2 in the vapor if H present, CO in the vapor if C present, CH4 in the vapor if both H and C present).
   - ``'H2O-CO2 only'`` — The only species present in the vapor are H2O and CO2 and in the melt are H2OT and CO2T (i.e., no CO, H2, and/or CH4 in the melt or vapor).

**H2S_m** *(default: True)*
   Is H2S a dissolved melt species.

   - ``'True'`` — Include H2Smol as a dissolved melt species.
   - ``'False'`` — H2Smol is insoluble in the melt.

**species X** *(default: Ar)*
   Chemical identity of species X, which defines its atomic mass. Other noble gases not currently supported, but we can add them if you get in touch!

   - ``'Ar'`` — Species X is argon (i.e., atomic mass of ~40).
   - ``'Ne'`` — Species X is Ne (i.e., atomic mass of ~20).

**Hspeciation** *(default: none)*

   - ``'none'`` — Oxidised H in the melt only occurs as H2OT species (i.e., no OH-).

Oxygen fugacity
^^^^^^^^^^^^^^^

**fO2** *(default: Kress91A)*
   Model for parameterisation of relationship between fO2 and Fe3+/FeT See :func:`~VolFe.model_dependent_variables.f_O2`.

   - ``'Kress91A'`` — Eq. (A-5, A-6) in Kress and Carmichael (1991) CMP 108:82-92 https//doi.org/10.1007/BF00307328
   - ``'Kress91'`` — Eq. (7) in Kress and Carmichael (1991) CMP 108:82-92 https//doi.org/10.1007/BF00307328
   - ``'ONeill18'`` — Eq. (9a) in O'Neill et al. (2018) EPSL 504:152-162 https//doi.org/10.1016/j.epsl.2018.10.0020012-821X
   - ``'Borisov18'`` — Eq. (4) from Borisov et al. (2018) CMP 173:98 https//doi.org/10.1007/s00410-018-1524-8

**NNObuffer** *(default: Frost91)*
   Model for the parameterisation for the fO2 value of the NNO buffer. See :func:`~VolFe.model_dependent_variables.NNO`.

   - ``'Frost91'`` — Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance" https//doi.org/10.1515/9781501508684-004

**FMQbuffer** *(default: Frost91)*
   Model for the parameterisation for the fO2 value of the FMQ buffer. See :func:`~VolFe.model_dependent_variables.FMQ`.

   - ``'Frost91'`` — Frost (1991) in "Oxide Minerals: Petrologic and Magnetic Significance" https//doi.org/10.1515/9781501508684-004
   - ``'ONeill87'`` — O'Neill (1897) AmMin 72(1-2):67-75

Solubility constants
^^^^^^^^^^^^^^^^^^^^

**carbon dioxide** *(default: MORB_Dixon95)*
   Model for the parameterisation of the CO2T solubility constant. See :func:`~VolFe.model_dependent_variables.C_CO3`.

   - ``'MORB_Dixon95'`` — Bullet (5) of summary from Dixon et al. (1995) JPet 36(6):1607-1631 https://doi.org/10.1093/oxfordjournals.petrology.a037267
   - ``'Basalt_Dixon97'`` — Eq. (7) from Dixon (1997) AmMin 82(3-4)368-378 https://doi.org/10.2138/am-1997-3-415
   - ``'NorthArchBasalt_Dixon97'`` — Eq. (8) from Dixon (1997) AmMin 82(3-4)368-378 https://doi.org/10.2138/am-1997-3-415
   - ``'Basalt_Lesne11'`` — Eq. (25,26) from Lesne et al. (2011) CMP 162:153-168 https://doi.org/10.1007/s00410-010-0585-0
   - ``'VesuviusAlkaliBasalt_Lesne11'`` — VES-9 in Table 4 from Lesne et al. (2011) CMP 162:153-168 https://doi.org/10.1007/s00410-010-0585-0
   - ``'EtnaAlkaliBasalt_Lesne11'`` — ETN-1 in Table 4 from Lesne et al. (2011) CMP 162:153-168 https://doi.org/10.1007/s00410-010-0585-0
   - ``'StromboliAlkaliBasalt_Lense11'`` — PST-9 in Table 4 from Lesne et al. (2011) CMP 162:153-168 https://doi.org/10.1007/s00410-010-0585-0
   - ``'Basanite_Holloway94'`` — Basanite in Table 5 from Holloway and Blank (1994) RiMG 30:187-230 https://doi.org/10.1515/9781501509674-012
   - ``'Leucitite_Thibault94'`` — Leucitite from Thibault & Holloway (1994) CMP 116:216-224 https://doi.org/10.1007/BF00310701
   - ``'Rhyolite_Blank93'`` — Fig.2 caption from Blank et al. (1993) EPSL 119:27-36 https://doi.org/10.1016/0012-821X(93)90004-S

**water** *(default: Basalt_Hughes24)*
   Model for the parameterisation for the H2O solubility constant. See :func:`~VolFe.model_dependent_variables.C_H2O`.

   - ``'Basalt_Hughes24'`` — Fig.S2 from Hughes et al. (2024) AmMin 109(3):422-438 https://doi.org/10.2138/am-2023-8739
   - ``'Rhyolite_HughesIP'`` — Fig.SX from Hughes et al. (in prep) based on data in Fig. 3 of Blank et al. (1993)

**hydrogen** *(default: Basalt_Hughes24)*
   Model for the parameterisation of the H2 solubility constant. See :func:`~VolFe.model_dependent_variables.C_H2`.

   - ``'Basalt_Hughes24'`` — Basalt in Table S4 from Hughes et al. (2024) https://doi.org/10.2138/am-2023-8739, based on experimental data from Hirschmann et al. (2012).
   - ``'Andesite_Hughes24'`` — Andesite in Table S4 from Hughes et al. (2024) https://doi.org/10.2138/am-2023-8739, based on experimental data from Hirschmann et al. (2012).

**sulfide** *(default: ONeill21dil)*
   Model for the parameterisation for the \*S2- solubility constant. See :func:`~VolFe.model_dependent_variables.C_S`.

   - ``'ONeill21dil'`` — Eq. (10.34) inc. H2O dilution from O'Neill (2021) in "Magma Redox Geochemistry" https://doi.org/10.1002/9781119473206.ch10
   - ``'ONeill21'`` — Eq. (10.34) ex. H2O dilution from O'Neill (2021) in "Magma Redox Geochemistry" https://doi.org/10.1002/9781119473206.ch10
   - ``'ONeill21hyd'`` — (hydrous) Eq. (10.34, 10.49) from O'Neill (2021) in "Magma Redox Geochemistry" https://doi.org/10.1002/9781119473206.ch10
   - ``'Boulliung23_eq6'`` — Eq. (6) from Boulliung & Wood (2023) CMP 178:56 https://doi.org10.1007/s00410-023-02033-9
   - ``'Boulliung23_eq7'`` — Eq. (7) from Boulliung & Wood (2023) CMP 178:56 https://doi.org10.1007/s00410-023-02033-9

**sulfate** *(default: ONeill22dil)*
   Model for the parameterisation of the S6+ solubility constant. See :func:`~VolFe.model_dependent_variables.C_SO4`.

   - ``'ONeill22dil'`` — Eq. (12a) inc. H2O dilution from O'Neill & Mavrogenes (2022) GCA 334:368-382 https://doi.org/10.1016/j.gca.2022.06.020
   - ``'ONeill22'`` — Eq. (12a) without H2O dilution from O'Neill & Mavrogenes (2022) GCA 334:368-382 https://doi.org/10.1016/j.gca.2022.06.020
   - ``'Boulliung22nP'`` — Eq. (5) from Boulliung & Wood (2023) GCA 343:420 https://doi.org/10.1016/j.gca.2022.11.025
   - ``'Boulliung22wP'`` — Eq. (5) from Boulliung & Wood (2023) GCA 343:420 https://doi.org/10.1016/j.gca.2022.11.025 and Eq. (8) for P from Boulliung & Wood (2022) GCA 336:150-164 https://doi.org/10.1016/j.gca.2022.08.032
   - ``'Boulliung23_eq9'`` — Eq. (9) from Boulliung & Wood (2023) CMP 178:56 https://doi.org/10.1007/s00410-023-02033-9
   - ``'Boulliung23_eq11'`` — Eq. (11) from Boulliung & Wood (2023) CMP 178:56 https://doi.org/10.1007/s00410-023-02033-9

**hydrogen sulfide** *(default: Basalt_Hughes24)*
   Model for the parameterisation for the H2S solubility constant. See :func:`~VolFe.model_dependent_variables.C_H2S`.

   - ``'Basalt_Hughes24'`` — Fig.S6 from Hughes et al. (2024) https://doi.org/10.2138/am-2023-8739, based on experimental data Moune et al. (2009) and calculations in Lesne et al. (2011).
   - ``'BasalticAndesite_Hughes24'`` — Fig.S6 from Hughes et al. (2024) https://doi.org/10.2138/am-2023-8739, based on experimental data Moune et al. (2009) and calculations in Lesne et al. (2011).

**methane** *(default: Basalt_Ardia13)*
   Model for the parameterisation of the CH4 solubility constant. See :func:`~VolFe.model_dependent_variables.C_CH4`.

   - ``'Basalt_Ardia13'`` — Eq. (7a) from Ardia et al. (2013) GCA 114:52-71 https://doi.org/10.1016/j.gca.2013.03.028

**carbon monoxide** *(default: Basalt_Hughes24)*
   Model for the parameterisation of the CO solubility constant. See :func:`~VolFe.model_dependent_variables.C_CO`.

   - ``'Basalt_Hughes24'`` — CO in Table S4 from Hughes et al. (2024) https://doi.org/10.2138/am-2023-8739, based on data from Armstrong et al. (2015), Stanley et al. (2014), and Wetzel et al. (2013).

**species X solubility** *(default: Ar_Basalt_HughesIP)*
   Model for the parameterisation of the X solubility constant. See :func:`~VolFe.model_dependent_variables.C_X`.

   - ``'Ar_Basalt_HughesIP'`` — Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
   - ``'Ar_Rhyolite_HughesIP'`` — Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
   - ``'Ne_Basalt_HughesIP'`` — Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
   - ``'Ne_Rhyolite_HughesIP'`` — Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157

**Cspeccomp** *(default: Basalt)*
   Model for the parameterisation of the speciation constant for CO2mol and CO32- in the melt. See :func:`~VolFe.model_dependent_variables.KCOm`.

   - ``'Basalt'`` — Assume all oxidised carbon in the melt is present as carbonate ions.
   - ``'Andesite_Botcharnikov06'`` — Eq. (8) from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143 https://doi.org/10.1016/j.chemgeo.2006.01.016
   - ``'Dacite_Botcharnikov06'`` — Eq. in the text from Botcharnikov et al. (2006) https://doi.org/10.1016/j.chemgeo.2006.01.016, based on data from Behrens et al. (2004)
   - ``'Rhyolite'`` — Assume all oxidised carbon in the melt is present as molecular CO2.

**Hspeccomp** *(default: MORB_HughesIP)*
   Model for the parameterisation of the speciation constant for H2Omol and OH- in the melt. See :func:`~VolFe.model_dependent_variables.KHOm`.

   - ``'MORB_HughesIP'`` — Eq. SX in Hughes et al. (in prep) based on data from Dixon et al. (1995)
   - ``'StromboliAlkaliBasalt_Lesne10'`` — Eq. (15) Lesne et al. (2010) CMP 162:133-151 https://doi.org/10.1007/s00410-010-0588-x
   - ``'VesuviusAlkaliBasalt_Lesne10'`` — Eq. (16) Lesne et al. (2010) CMP 162:133-151 https://doi.org/10.1007/s00410-010-0588-x
   - ``'EtnaAlkaliBasalt_Lesne10'`` — Eq. (17) Lesne et al. (2010) CMP 162:133-151 https://doi.org/10.1007/s00410-010-0588-x
   - ``'Andesite_Botcharnikov06'`` — Eq (7) from Botcharnikov et al. (2006) Chem. Geol. 229(1-3)125-143 https://doi.org/10.1016/j.chemgeo.2006.01.016
   - ``'Albite_Silver89'`` — Fig. 8 from Silver & Stolper (1989) J.Pet 30(3)667-709 https://doi.org/10.1093/petrology/30.3.667
   - ``'Rhyolite_Zhang97'`` — Eq. (9) from Zhang et al. (1997) GCA 61(15):3089-3100 https://doi.org/10.1016/S0016-7037(97)00151-8

Saturation conditions
^^^^^^^^^^^^^^^^^^^^^

**SCSS** *(default: ONeill21hyd)*
   Model for parameterisation of the sulfide content at sulfide saturation (S2-CSS). See :func:`~VolFe.model_dependent_variables.SCSS`.

   - ``'ONeill21hyd'`` — Eq. (10.34, 10.43, 10.45, 10.46, 10.49) from O'Neill (2021) in "Magma Redox Geochemistry" https://doi.org/10.1002/9781119473206.ch10
   - ``'ONeill21'`` — Eq. (10.34, 10.43, 10.45, 10.46) excluding water dilution from O'Neill (2021) in "Magma Redox Geochemistry" https://doi.org/10.1002/9781119473206.ch10
   - ``'ONeill21dil'`` — Eq. (10.34, 10.43, 10.45, 10.46) including water dilution from O'Neill (2021) in "Magma Redox Geochemistry" https://doi.org/10.1002/9781119473206.ch10
   - ``'Liu07'`` — Eq. (9) in Liu et al. (2007) GCA 71:1783-1799 https://doi.org/10.1016/j.gca.2007.01.004
   - ``'Fortin15_pss'`` — Eq. (7) in Fortin et al. (2015) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127
   - ``'Liu21_pss'`` — Eq. (2) in Liu et al. (2021) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127
   - ``'ONeill22_pss'`` — O'Neill & Mavrogenes (2022) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127
   - ``'ONeill21_pss'`` — O'Neill (2021) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127
   - ``'Smythe17_pss'`` — Smythe et al. (2017) using PySulfSat by Wieser & Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127
   - ``'Li22_pss'`` — Eq. (19) from Li and Zhang (2022) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127
   - ``'Blanchard21_eq11_pss'`` — Eq. (11) from Blanchard et al. (2021) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127
   - ``'Blanchard21_eq12_pss'`` — Eq. (12) from Blanchard et al. (2021) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127

**SCAS** *(default: Zajacz19_pss)*
   Model for parameterisation of the sulfate content at anhydrite saturation (S6+CAS). See :func:`~VolFe.model_dependent_variables.SCAS`.

   - ``'Liu23'`` — Eq. (4) Liu et al. (2023) GCA 349:135-145 https://doi.org/10.1016/j.gca.2023.04.007
   - ``'Chowdhury19_pss'`` — Chowdhury & Dasgupta (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127
   - ``'Zajacz19_pss'`` — Zajacz and Tsay (2019) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127
   - ``'Masotta15_pss'`` — Masotta and Kepler (2015) using PySulfSat by Wieser and Gleeson (2023) Volcanica 6(1):107-127 https://doi.org/10.30909/vol.06.01.107127

**sulfur_saturation** *(default: False)*
   Is sulfur allowed to form sulfide or anhydrite if sulfur content of the melt reaches saturation levels for these phases?

   - ``'False'`` — melt ± vapor are the only phases present - results are metastable with respect to sulfide and anhydrite if they could saturate.
   - ``'True'`` — If saturation conditions for sulfide or anhydrite are met, melt sulfur content reflects this.

**graphite_saturation** *(default: False)*
   Is graphite allowed to form if the carbon content of the melt reaches saturation levels for graphite?

   - ``'False'`` — melt ± vapor are the only phases present - results are metastable with respect to graphite if it could saturate.
   - ``'True'`` — If saturation conditions for graphite are met, melt carbon content reflects this.

Fugacity coefficients
^^^^^^^^^^^^^^^^^^^^^

**ideal_gas** *(default: False)*
   Treat all vapor species as ideal gases (i.e., all fugacity coefficients = 1 at all P).

   - ``'False'`` — At least some of the vapor species are not treated as ideal gases.
   - ``'True'`` — All fugacity coefficients = 1 at all P.

**y_CO2** *(default: Shi92)*
   Model for the parameterisation of the CO2 fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_CO2`.

   - ``'Shi92'`` — Shi & Saxena (1992) AmMin 77(9-10):1038-1049
   - ``'Holland91_eq8_tab1'`` — Eq. (8) and Table 1 from Holland & Powell (1991) CMP 109:265-273 https://doi.org/10.1007/BF00306484
   - ``'Holland91_eq4,A1-3_tab1'`` — Eq. (4,A1-3) and Table 1 from Holland & Powell (1991) CMP 109:265-273 https://doi.org/10.1007/BF00306484
   - ``'Holland91_eq8,9_tab2'`` — Eq. (8,9) and Table 2 from Holland & Powell (1991) CMP 109:265-273 https://doi.org/10.1007/BF00306484
   - ``'Flowers79'`` — Flowers (1979) modified from code from MIMiC (Rasmussen et al., 2021: https://github.com/DJRgeoscience/MIMiC), originally from VolatileCalc (Newman & # Lowenstern, 2001)
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_SO2** *(default: Shi92_Hughes23)*
   Model for the parameterisation of the SO2 fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_SO2`.

   - ``'Shi92_Hughes23'`` — Fig.S1 modified from Shi & Saxena (1992) from Hughes et al. (2023) JGSL 180(3) https://doi.org/10.1144/jgs2021-12
   - ``'Shi92'`` — Shi & Saxena (1992) AmMin 77(9-10):1038-1049
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_H2S** *(default: Shi92_Hughes24)*
   Model for the parameterisation of the H2S fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_H2S`.

   - ``'Shi92_Hughes23'`` — Fig.S1 modified from Shi & Saxena (1992) Hughes et al. (2024) AmMin 109(3):422-438 https://doi.org/10.2138/am-2023-8739
   - ``'Shi92'`` — Shi & Saxena (1992) AmMin 77(9-10):1038-1049
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_H2** *(default: Shaw64)*
   Model for the parameterisation of the H2 fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_H2`.

   - ``'Shaw64'`` — Eq. (4) from Shaw & Wones (1964) AmJSci 262:918-929
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_O2** *(default: Shi92)*
   Model for the parameterisation of the O2 fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_O2`.

   - ``'Shi92'`` — Shi & Saxena (1992) AmMin 77(9-10):1038-1049
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_S2** *(default: Shi92)*
   Model for the parameterisation of the S2 fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_S2`.

   - ``'Shi92'`` — Shi & Saxena (1992) AmMin 77(9-10):1038-1049
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_CO** *(default: Shi92)*
   Model for the parameterisation of the CO fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_CO`.

   - ``'Shi92'`` — Shi & Saxena (1992) AmMin 77(9-10):1038-1049
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_CH4** *(default: Shi92)*
   Model for the parameterisation of the CH4 fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_CH4`.

   - ``'Shi92'`` — Shi & Saxena (1992) AmMin 77(9-10):1038-1049
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_H2O** *(default: Holland91)*
   Model for the parameterisation of the H2O fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_H2O`.

   - ``'Holland91'`` — Eq. (4,6,A1-3) and Table 1 (T > 673 K only) from Holland & Powell (1991) CMP 109:265-273 https://doi.org/10.1007/BF00306484
   - ``'Flowers79'`` — Flowers (1979) modified from code from MIMiC (Rasmussen et al., 2021: https://github.com/DJRgeoscience/MIMiC), originally from VolatileCalc (Newman & Lowenstern, 2001)
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_OCS** *(default: Shi92)*
   Model for the parameterisation of the OCS fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_OCS`.

   - ``'Shi92'`` — Shi & Saxena (1992) AmMin 77(9-10):1038-1049
   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

**y_X** *(default: ideal)*
   Model for the parameterisation of the X fugacity coefficient. See :func:`~VolFe.model_dependent_variables.y_X`.

   - ``'ideal'`` — Treat as ideal gas, y = 1 at all P.

Equilibrium constants
^^^^^^^^^^^^^^^^^^^^^

**KHOg** *(default: Ohmoto97)*
   Model for the parameterisation of the equilibiurm constant for H2 + 0.5O2 ⇄ H2O. See :func:`~VolFe.model_dependent_variables.KHOg`.

   - ``'Ohmoto97'`` — Reaction (d) in Table 1 of Ohmoto & Kerrick (1997).

**KHOSg** *(default: Ohmoto97)*
   Model for the parameterisation of the equilibiurm constant for 0.5S2 + H2O ⇄ H2S + 0.5O2. See :func:`~VolFe.model_dependent_variables.KHOSg`.

   - ``'Ohmoto97'`` — Reaction (h) in Table 1 of Ohmoto & Kerrick (1997).
   - ``'noH2S'`` — Stops H2S forming in the vapor, K = 0.

**KOSg** *(default: Ohmoto97)*
   Model for the parameterisation of the equilibiurm constant for 0.5S2 + O2 ⇄ SO2. See :func:`~VolFe.model_dependent_variables.KOSg`.

   - ``'Ohmoto97'`` — Reaction (f) in Table 1 of Ohmoto & Kerrick (1997).
   - ``'noSO2'`` — Stops SO2 forming in the vapor, K = 0.

**KOSg2** *(default: ONeill22)*
   Model for the parameterisation of the equilibiurm constant for 0.5S2 + 1.5O2 ⇄ SO3. See :func:`~VolFe.model_dependent_variables.KOSg2`.

   - ``'ONeill22'`` — Eq (6b) in O’Neill and Mavrogenes (2022) https://doi.org/10.1016/j.gca.2022.06.020

**KOCg** *(default: Ohmoto97)*
   Model for the parameterisation of the equilibiurm constant for CO + 0.5O2 ⇄ CO2. See :func:`~VolFe.model_dependent_variables.KCOg`.

   - ``'Ohmoto97'`` — Reaction (c) in Table 1 of Ohmoto & Kerrick (1997).

**KCOHg** *(default: Ohmoto97)*
   Model for the parameterisation of the equilibiurm constant for CH4 + 2O2 ⇄ CO2 + 2H2O. See :func:`~VolFe.model_dependent_variables.KCOHg`.

   - ``'Ohmoto97'`` — Reaction (e) in Table 1 of Ohmoto & Kerrick (1997).
   - ``'noCH4'`` — Almost stops CH4 forming in the vapor, K = very large.

**KOCSg** *(default: Moussallam19)*
   Model for the parameterisation of the equilibiurm constant for OCS. See :func:`~VolFe.model_dependent_variables.KOCSg`.

   - ``'Moussallam19'`` — Eq. (8) in Moussallam et al. (2019) https://doi.org/10.1016/j.epsl.2019.05.036 for KOCSg and 'COS' for carbonlysulfide
   - ``'noOCS'`` — Almost stops OCS forming in the vapor, K = very large.

**KCOs** *(default: Holloway92)*
   Model for the parameterisation of the equilibiurm constant for Cgrahite + O2 ⇄ CO2. See :func:`~VolFe.model_dependent_variables.KCOs`.

   - ``'Holloway92'`` — Eq (3) KI in Holloway et al. (1992) Eur J. Mineral. 4:105-114

**carbonylsulfide** *(default: COS)*
   Reaction equilibrium KOCSg is for.

   - ``'COS'`` — 2CO2 + OCS ⇄ 3CO + SO2

Degassing calculation
^^^^^^^^^^^^^^^^^^^^^

**bulk_composition** *(default: melt-only)*
   Specifying what the inputted melt composition (i.e., dissolved volatiles and fO2-estimate) corresponds to for the re/degassing calculation.

   - ``'melt-only'`` — The inputted melt composition (i.e., dissolved volatiles) represents the bulk system - there is no vapor present. The fO2-estimate is calculated at Pvsat for this melt composition.
   - ``'melt+vapor_wtg'`` — The inputted melt composition (i.e., dissolved volatiles) is in equilibrium with a vapor phase. The amount of vapor as weight fraction gas (wtg) is specified in the inputs. The bulk system composition will be calculated by calculating Pvsat and the vapor composition given the input composition.
   - ``'melt+vapor_initialCO2'`` — The inputted melt composition (i.e., dissolved volatiles) is in equilibrium with a vapor phase. The initial CO2 content of the melt (i.e., before degassing) is specified in the inputs. The bulk system composition will be calculated by calculating Pvsat and the vapor composition given the input composition. The amount of vapor present is calculated using the given initial CO2.

**starting_P** *(default: Pvsat)*
   Determining the starting pressure for a re/degassing calculation.

   - ``'Pvsat'`` — Calculation starts at Pvsat for the inputted melt composition (i.e., dissolved volatiles).

**gassing_style** *(default: closed)*
   Does the bulk composition of the system (including oxygen) remain constant during the re/degassing calculation.

   - ``'closed'`` — The bulk composition of the system (inc. oxygen) is constant during the re/degassing calculation - vapor and melt remain in chemical equilibrium throughout.
   - ``'open'`` — At each pressure-step, the vapor in equilibrium with the melt is removed (or added for regassing), such that the bulk composition of the system changes. This does not refer to being externally buffered in terms of fO2.

**gassing_direction** *(default: degas)*
   Is pressure increasing or decreasing from the starting pressure.

   - ``'degas'`` — Pressure progressively decreases for isothermal, polybaric calculations (i.e., degassing).
   - ``'regas'`` — Pressure progressively increases for isothermal, polybaric calculations (i.e., regassing).

**P_variation** *(default: polybaric)*
   Is pressure varying during the calculation?

   - ``'polybaric'`` — Pressure progressively changes during the calculation.

**T_variation** *(default: isothermal)*
   Is temperature varying during the calculation?

   - ``'isothermal'`` — Temperature is constant during the calculation.

**eq_Fe** *(default: yes)*
   Does iron in the melt equilibrate with fO2.

   - ``'yes'`` — Iron equilibrates with fO2

Other
^^^^^

**density** *(default: DensityX)*
   Model for parameterisation of melt density. See :func:`~VolFe.model_dependent_variables.melt_density`.

   - ``'DensityX'`` — DensityX from Iacovino & Till (2019) Volcanica 2(1):1-10 https//doi.org/10.30909/vol.02.01.0110

**print status** *(default: False)*
   Specifies whether some sort of status information during the calculation is outputted to let you know progress.

   - ``'False'`` — No information about calculation progress is printed.
   - ``'True'`` — Some information about calculation progress is printed.

**output csv** *(default: True)*
   Specicies whether a csv of the outputted dataframe is automatically saved at the end of the calculation.

   - ``'True'`` — csv is outputted.
   - ``'False'`` — csv is not outputted.

In development
^^^^^^^^^^^^^^

The following parameters are under active development. For now, leave them at their default values.

**isotopes** *(default: no)*

**crystallisation** *(default: no)*

**mass_volume** *(default: mass)*

**calc_sat** *(default: fO2_melt)*

**bulk_O** *(default: exc_S)*

**error** *(default: 0.1)*

**sulfur_is_sat** *(default: no)*

**solve_species** *(default: auto)*

**setup** *(default: False)*

**high precision** *(default: False)*
- Isotope fractionation factors: functions are ``alpha_`` (e.g., ``alpha_S_SO2v_S6pm`` is the function for calculating the isotopic fractionation factor between SO\ :sub:`2` in the vapor and S\ :sup:`6+` in the melt - ``help(volfe.alpha_S_SO2v_S6pm)``)
