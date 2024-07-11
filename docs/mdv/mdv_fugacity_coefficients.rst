===================================================================================
Fugacity coefficients
===================================================================================

Below are the different models for the fugacity coefficients currently in VolFe, which can also be viewed in the :doc:`fugacity coefficients notebook <FugCoeff>`. 
VolFe treats the vapor as an ideal mixture of non-ideal gases; hence, these parameterisations depend no *P* and *T* but *not* the vapor composition.

**ideal_gas**: Treat all vapor species as ideal gases (i.e., all fugacity coefficients = 1 at all P).
    
- default: 'False' At least some of the vapor species are not treated as ideal gases. 
        
- 'True' All fugacity coefficients = 1 at all P.
    
**y_CO2**: Model for the parameterisation of the CO2 fugacity coefficient.

- default: 'Shi92' Shi & Saxena (1992).

- 'Holland91' Holland & Powell (1991).

- 'ideal' Treat CO2 as ideal gas species, fugacity coefficient = 1 at all P.

**y_SO2**: Model for the parameterisation of the SO2 fugacity coefficient.

- default: 'Shi92_Hughes23' Fig.S1 Hughes et al. (2023).

- 'Shi92' Shi & Saxena (1992).

- 'ideal' Treat SO2 as ideal gas species, fugacity coefficient = 1 at all P.

**y_H2S**: Model for the parameterisation of the H2S fugacity coefficient.

- default: 'Shi92_Hughes24' Fig.S1 Hughes et al. (2024).

- 'Shi92' Shi & Saxena (1992).

- 'ideal' Treat H2S as ideal gas species, fugacity coefficient = 1 at all P.

**y_H2**: Model for the parameterisation of the H2 fugacity coefficient.

- default: 'Shaw64' Eq. (4) from Shaw & Wones (1964). [Note: used value of -0.011901 instead of -0.11901 as reported to match their Table 2]

- 'ideal' Treat H2 as ideal gas species, fugacity coefficient = 1 at all P.

**y_O2**: Model for the parameterisation of the O2 fugacity coefficient.

- default: 'Shi92' Shi & Saxena (1992).

- 'ideal' Treat O2 as ideal gas species, fugacity coefficient = 1 at all P.

**y_S2**: Model for the parameterisation of the O2 fugacity coefficient.
        
- default: 'Shi92' Shi & Saxena (1992).

- 'ideal' Treat S2 as ideal gas species, fugacity coefficient = 1 at all P.

**y_CO**: Model for the parameterisation of the CO fugacity coefficient.

- default: 'Shi92' Shi & Saxena (1992).

- 'ideal' Treat CO as ideal gas species, fugacity coefficient = 1 at all P.        

**y_CH4**: Model for the parameterisation of the CH4 fugacity coefficient.

- default: 'Shi92' Shi & Saxena (1992).

- 'ideal' Treat CH4 as ideal gas species, fugacity coefficient = 1 at all P.    

**y_H2O** : Model for the parameterisation of the H2O fugacity coefficient.

- default: 'Holland91' Holland & Powell (1991).

- 'ideal' Treat H2O as ideal gas species, fugacity coefficient = 1 at all P.    
    
**y_OCS**: Model for the parameterisation of the OCS fugacity coefficient.

- default: 'Shi92' Shi & Saxena (1992).

- 'ideal' Treat OCS as ideal gas species, fugacity coefficient = 1 at all P.            

**y_X**: Model for the parameterisation of the X fugacity coefficient.

- default: 'ideal' Treat X as ideal gas species, fugacity coefficient = 1 at all P.

- Only one option available currently, included for future development.