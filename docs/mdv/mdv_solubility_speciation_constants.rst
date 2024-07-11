===================================================================================
Solubility constants
===================================================================================

Below are the different models for the solubility constants currently in VolFe, which can also be viewed in the :doc:`solubility constants notebook <SolSpecConst>`.
This is how VolFe incorporates the effect of melt composition on the solubility and speciation of volatile species in the melt.

### IN PROGRESS ###

**carbon dioxide**: Model for the parameterisation of the CO2T solubility constant.

- default: 'MORB_Dixon95' Bullet (5) of summary from Dixon et al. (1995).

- 'Basalt_Dixon97' Eq. (7) from Dixon et al. (1997).

- 'NorthArchBasalt_Dixon97' Eq. (8) from Dixon et al. (1997).

- 'Basalt_Lesne11' Eq. (25,26) from Lesne et al. (2011).

- 'VesuviusAlkaliBasalt_Lesne11' VES-9 in Table 4 from Lesne et al. (2011).

- 'EtnaAlkaliBasalt_Lesne11' ETN-1 in Table 4 from Lesne et al. (2011).

- 'StromboliAlkaliBasalt_Lense11' PST-9 in Table 4 from Lesne et al. (2011).

- 'SunsetCraterAlkaliBasalt_Allison19' Sunset Crater in Table 4 from Allison et al. (2019).

- 'SVFVBasalticAndesite_Allison19' SVFV in Table 4 from Allison et al. (2019).

- 'ErebusPhonotephrite_Allison19' Erebus in Table 4 from Allison et al. (2019).

- 'VesuviusPhonotephrite_Allison19' Vesuvius in Table 4 from Allison et al. (2019).

- 'EtnaTrachybasalt_Allison19' Etna in Table 4 from Allison et al. (2019).

- 'StromboliAlkaliBasalt_Allison19' Stromboli in Table 4 from Allison et al. (2019).

- 'Basanite_Holloway94' Basanite in Table 5 from Holloway and Blank (1994).

- 'Leucitite_Thibault94' Leucitite from Thibault & Holloway (1994).

- 'TholeiiteBasalt_Allison22' N72 basalt in Table 2 from Allison et al. (2022).

- 'Rhyolite_Blank93' Fig.2 caption from Blank et al. (1993).

    ### water: Model for the parameterisation for the H2O solubility constant.
        default: 'Basalt_Hughes24' Hughes et al. (2024) AmMin 109(3):422-438 based on data compiliation from Allison et al. (2022) CMP 177(3):40 doi:10.1007/s00410-022-01903-y
    
    hydrogen: Model for the parameterisation of the H2 solubility constant.
        default: 'Basalt_Hughes24' Basalt H2 in Table S4 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
        Other options:
        'Andesite_Hughes24' Andesite H2 in Table S4 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
    
    sulfide: Model for the parameterisation for the *S2- solubility constant (all calibrated over wide range of silicate melt compositions).
        default: 'ONeill21dil' Eq. (10.34) inc. H2O dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        Other options:
        'ONeill21' Eq. (10.34) ex. H2O dilution from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        'ONeill21hyd' (hydrous) Eq. (10.34, 10.49) from O'Neill (2021) in "Magma Redox Geochemistry" doi:10.1002/9781119473206.ch10
        'Boulliung23eq6' Eq. (6) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9
        'Boulliung23eq7' Eq. (7) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9 
    
    sulfate: Model for the parameterisation of the S6+ solubility constant (all calibrated over wide range of silicate melt compositions).
        default: 'ONeill22dil' Eq. (12a) inc. H2O dilution from O'Neill & Mavrogenes (2022) GCA 334:368-382 10.1016/j.gca.2022.06.020
        Other options:
        'ONeill22' Eq. (12a) without H2O dilution from O'Neill & Mavrogenes (2022) GCA 334:368-382 doi:10.1016/j.gca.2022.06.020
        'Boulliung22nP' (no P-dependence) Eq. (5) from Boulliung & Wood (2023) GCA 343:420 doi:10.1016/j.gca.2022.11.025
        'Boulliung22wP' (inc. P-dependece) Eq. (5) from Boulliung & Wood (2023) GCA 343:420 doi:10.1016/j.gca.2022.11.025 and Eq. (8) for P from Boulliung & Wood (2022) GCA 336:150-164 doi:10.1016/j.gca.2022.08.032
        'Boulliung23eq9' Eq. (9) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9
        'Boulliung23eq11' Eq. (11) from Boulliung & Wood (2023) CMP 178:56 doi:10.1007/s00410-023-02033-9
    
    hydrogen sulfide: Model for the parameterisation for the H2S solubility constant.
        default 'Basalt_Hughes24' Fig.S6 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
        Other options:
        'BasalticAndesite_Hughes24' Fig.S6 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
    
    methane: Model for the parameterisation of the CH4 solubility constant.
        default: 'Basalt_Ardia13' Eq. (7a) from Ardia et al. (2013) GCA 114:52-71 doi:10.1016/j.gca.2013.03.028
        Only one option available currently, included for future development.

    carbon monoxide: Model for the parameterisation of the CO solubility constant.
        default: 'Basalt_Hughes24' CO in Table S4 from Hughes et al. (2024) AmMin 109(3):422-438 doi:10.2138/am-2023-8739
        Only one option available currently, included for future development.

    species X solubility: Model for the parameterisation of the X solubility constant. 
        default: 'Ar_Basalt_HughesIP' Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        Other options:
        Ar_Rhyolite_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        Ne_Basalt_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
        Ne_Rhyolite_HughesIP: Hughes et al. (in prep) based on data from Iacono-Marziano et al. (2010) Chemical Geology 279(3–4):145-157
    
    Cspeccomp: Model for the parameterisation of the speciation constant for CO2mol and CO32- in the melt.
        default: 'Basalt' Assume all oxidised carbon in the melt is present as carbonate ions.
        Other options:
        Andesite_Botcharnikov06: Eq. (8) from Botcharnikov et al. (2006) Chem.Geol. 229(1-3)125-143 doi:10.1016/j.chemgeo.2006.01.016
        Dacite_Botcharnikov06: Botcharnikov et al. (2006) Chem.Geol. 229(1-3)125-143 doi:10.1016/j.chemgeo.2006.01.016
        Rhyolite: Assume all oxidised carbon in the melt is present as molecular CO2.
