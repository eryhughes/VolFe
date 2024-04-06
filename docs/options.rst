=========================
Options within VolFe
=========================

There are various options for how the calculations are done in VolFe.
This varies from which species are treated as insoluble in the melt to whether a csv of the results is created at the end of the calcualtion.

Below are descriptions of the various options (option type is in **bold**, the different options themselves in *italics*, and the default option is labelled [default]):



Volatile species
------

Broadly, options related to the types of species VolFe is considering.
This is mostly done automatically by the volatile elements included in the setup file, but there are a few options you can decide.


**insolubles**

Specifying if H2, CO, and CH4 are present in the melt and/or vapor.
        
- *yes*: Include H2mol, COmol and/or CH4mol as dissolved melt species (which species depends on whether H and/or C are volatile elements). [default]

- *no*: H2, CO and/or CH4 are insoluble in the melt.

- *H2O-CO2 only*: The only species present in the vapor are H2O and CO2 and in the melt are H2OT and CO2T (i.e., no CO, H2, and/or CH4 in the melt or vapor).


**H2S_m**

Specify if H2S is a dissolved melt species.

- *yes*: Include H2Smol as a dissolved melt species. [default]
- *no*: H2Smol is insoluble in the melt.


**species X**

Chemical identity of species X, which defines its atomic mass. 
Note that its solubility is defined by **species X solubility** - see :doc:`solubility and speciation constants <mdv/mdv_solubility_speciation_constants>` for more details.

- *Ar*: Species X is argon (i.e., atomic mass of ~40). [default]

- *Ne*: Species X is Ne (i.e., atomic mass of ~20).


**Hspeciation**

[In progress]
              


Saturation of sulfide, anhydrite, and graphite
-----

**sulfur_saturation**

Is sulfur allowed to form sulfide or anhydrite if sulfur content of the melt reaches saturation levels for these phases.

- *no*: melt ± vapor are the only phases present - results are metastable with respect to sulfide and anhydrite if they could saturate. [default]

- *yes*: If saturation conditions for sulfide or anhydrite are met, melt sulfur content reflects this.


**graphite_saturation**

Is graphite allowed to form if the carbon content of the melt reaches saturation levels for graphite.

- *no*: melt ± vapor are the only phases present - results are metastable with respect to graphite if it could saturate. [default]

- *yes*: If saturation conditions for graphite are met, melt carbon content reflects this.



Degassing calculations
-----

Options specific to how the degassing calculation is run.


**bulk_composition**

Specifying what the inputted melt composition (i.e., dissolved volatiles and fO2-estimate) correspond to for the degassing calculation

- *yes*: The inputted melt composition (i.e., dissolved volatiles) represents the bulk system - there is no vapor present. The fO2-estimate is calculated at Pvsat for this melt composition. [default]

- *wtg*: The inputted melt composition (i.e., dissolved volatiles) is in equilibrium with a vapor phase. The amount of vapor is specified in the inputs. The bulk system composition will be calculated by calculating Pvsat and the vapor composition given the input composition.

- *CO2*: The inputted melt composition (i.e., dissolved volatiles) is in equilibrium with a vapor phase. The initial CO2 content of the melt (i.e., before degassing) is specified in the inputs. The bulk system composition will be calculated by calculating Pvsat and the vapor composition given the input composition.


**starting_P**

Determing the starting pressure for a degassing calculation.

- *bulk*: Calculation starts at Pvsat for the inputted melt composition (i.e., dissolved volatiles), which has no vapor present. [default]

- *set*: Calculation starts at the pressure specified in the inputs.

- *measured*: Calculation starts at Pvsat for the inputted melt composition (i.e., dissolved volatiles), which has vapor present.


**gassing_style**

Does the bulk composition of the system (including oxygen) remain constant during the re/degassing calculation.

- *closed*: The bulk composition of the system (inc. oxygen) is constant during re/degassing calculation - vapor and melt remain in chemical equilibrium throughout. [default]

- *open*: At each pressure-step, the vapor in equilibrium with the melt is removed (or added for regassing), such that the bulk composition of the system changes. This does not refer to being buffered in terms of fO2.


**gassing_direction**

Is pressure increasing or decreasing from the starting perssure.

- *degas*: Pressure progressively decreases from starting pressure for isothermal, polybaric calculations (i.e., degassing). [default]

- *regas*: Pressure progressively increases from starting pressure for isothermal, polybaric calculations (i.e., regassing). 
    

**P_variation**

Is pressure varying during the calculation?

- *polybaric*: Pressure progressively changes during the calculation. [default]

- Only one option available currently, included for future development.
    

**T_variation**

Is temperature varying during the calculation?

- *isothermal*: Temperature is constant during the calculation. [default]

- Only one option available currently, included for future development.
     

**solve_species**

What species are used to solve the equilibrium equations? This should not need to be changed unless the solver is struggling.

- *OCS*: Guess mole fractions of O2, CO, and S2 in the vapor to solve the equilibrium equations. [default]

- *OHS*: Guess mole fractions of O2, H2, and S2 in the vapor to solve the equilibrium equations.

- *OCH*: Guess mole fractions of O2, CO, and H2 in the vapor to solve the equilibrium equations.             



Other
----

**setup**

Specifies whether model options are specified in the models or setup dataframe. 
- *no*: All model options are specified in the models dataframe. [default]

- *yes*: Some of the model options are specified in the setup dataframe.


**print status**

Specifies whether some sort of status information during the calculation is outputted to let you know progress.

- *no*: No information about calculation progress is printed. [default]

- *yes*: Some information about calculation progress is printed.


**output csv**

Specicies whether a csv of the outputted dataframe is saved at the end of the calculation. 

- *yes*: csv is outputted [default]

- *no*: csv is not outputted    



In development
----

The following options are in development.
For now, just leave them as their default option and everything should work fine!

- **isotopes**: default = *no*

- **crystallisation**: default = *no*

- **mass_volume**: default = *mass*

- **calc_sat**: default = *fO2_melt*

- **bulk_O**: default = *exc_S*

- **error**: default = *0.1*

- **eq_Fe**: default = *yes*

- **sulfur_is_sat**: default = *no*