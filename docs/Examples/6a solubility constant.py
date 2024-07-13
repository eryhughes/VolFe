### CURRENT SOLUBILITY CONSTANT FOR H2S 

###################################
### solubility constant for H2S ###
###################################
def C_H2S(run,PT,melt_wf,setup,species,models): # C_H2S = wmH2S/fH2S (ppm H2S, fH2S bar)
    model = models.loc["hydrogen sulfide","option"]
    if model == "basalt":
        K = 10.23 # fitted to basalt data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015 ChemGeol 418:104–116
    elif model == "basaltic andesite":
        K = 6.82 # fitted to basaltic andesite data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015 ChemGeol 418:104–116 
    return K


### ADDING "my new amazing parameterisation" option

###################################
### solubility constant for H2S ###
###################################
def C_H2S(run,PT,melt_wf,setup,species,models): # C_H2S = wmH2S/fH2S (ppm H2S, fH2S bar)
    model = models.loc["hydrogen sulfide","option"]
    if model == "basalt":
        K = 10.23 # fitted to basalt data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015 ChemGeol 418:104–116
    elif model == "basaltic andesite":
        K = 6.82 # fitted to basaltic andesite data from Moune+ 2009 CMP 157:691–707 and Lesne+ 2015 ChemGeol 418:104–116 
    elif model == “my new amazing parameterisation”:
        P, T_K = PT['P'], PT['T']+273.15 # bars, K # extracts pressure in bars and temperature in ‘C from PT and convert temperature to K.
        tot, Si, Ti, Al, FeT, Fe2, Fe3, Mg, Mn, Ca, Na, K, P_, H, C = mg.melt_cation_proportion(run,melt_wf,setup,species,"no","no") # calculates cation mole fractions from the melt composition without volatiles in the melt and does not consider iron speciation (i.e., all Fe is FeO).
        A, B, C = 2.5, 6.8, -9. # sets the constants A, B, and C
        C0 = A*Si + B*Mg + C*K # calculates the compositional term from the cation mole fractions
        DV = 10. # cm3/mol # set volume change term
        DH = -13. # kJ/mol # sets enthalpy change term
        P0 = 1. # bar # sets reference pressure 
        T0 = 1400 # K # sets reference temperature
        R = 83.15 # gives value for R the gas constant
        B = -((DV/(R*T_K))*(P-P0)) + (DH/R)*((1.0/T0) - (1.0/T_K)) # calculates PT term
        K = C0*gp.exp(B) # calculates solubility constant in ppmw/kbar
        K = K/1000. # converts fH2S in kbars to fH2S in bars
    return K

