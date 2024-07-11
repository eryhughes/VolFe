===================================================================================
Equilibrium constants
===================================================================================

Below are the different models for the equilibrium constants currently in VolFe, which can also be viewed in the :doc:`equilibrium constants notebook <EqConst>`. 
The calibration range for these parametrisations is not listed in the original references but is suitable for magmatic systems.
Currently there is only one option for each equilibrium constant (except H2S), but more could be added in the future. 

**KHOg:** Model for the parameterisation of the equilibiurm constant for H2 + 0.5O2 = H2O.
    
- default: 'Ohmoto97' Reaction (d) in Table 1 from Ohmoto & Kerrick (1977).
    
- Only one option available currently, included for future development.

**KHOSg**: Model for the parameterisation of the equilibiurm constant for 0.5S2 + H2O = H2S + 0.5O2.
        
- default: 'Ohmoto97' Reaction (h) in Table 1 from Ohmoto & Kerrick (1977).
    
- 'no H2S' Stops H2S forming in the vapor (K = 0).
    
**KOSg**: Model for the parameterisation of the equilibiurm constant for 0.5S2 + O2 = SO2.
    
- default: 'Ohmoto97' Reaction (f) in Table 1 from Ohmoto & Kerrick (1977).
    
- 'no SO2' Stops SO2 forming in the vapor (K = 0). As a by-product, OCS will also stop forming.

**KOSg2**: Model for the parameterisation of the equilibiurm constant for 0.5S2 + 1.5O2 = SO3.
    
- default: 'ONeill2' Eq. (6b) from O'Neill & Mavrogenes (2022).
    
- Only one option available currently, included for future development.

**KOCg**: Model for the parameterisation of the equilibiurm constant for CO + 0.5O2 = CO2.
        
- default: 'Ohmoto97' Reaction (c) in Table 1 from Ohmoto & Kerrick (1977).
    
- Only one option available currently, included for future development. 

**KCOHg**: Model for the parameterisation of the equilibiurm constant for CH4 + 2O2 = CO2 + 2H2O.
    
- default: 'Ohmoto97' Reaction (e) in Table 1 from Ohmoto & Kerrick (1977).
    
- 'no CH4' Stops CH4 forming in the vapor (K = large number).

**KOCSg**: Model for the parameterisation of the equilibiurm constant for OCS.
    
- default: 'Moussallam19' Eq. (8) for 2CO2 + OCS ⇄ 3CO + SO2 in Moussallam et al. (2019).
    
- 'no OCS' Stops OCS forming in the vapor (K = large number).  


**KCOs**: Model for the parameterisation of the equilibiurm constant for Cgrahite + O2 = CO2.
    
- default: 'Holloway92' Eq. (3) KI in Holloway et al. (1992).
    
- Only one option available currently, included for future development.

**carbonylsulfide**: Reaction equilibrium KOCSg is for. 
    
- default: 'COS' 2CO2 + OCS ⇄ 3CO + SO2
    
- Only one option available currently, included for future development.