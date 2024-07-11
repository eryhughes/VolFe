===================================================================================
Equilibrium constants
===================================================================================

Below are the different models for the equilibrium constants currently in VolFe, which can also be viewed in the :doc:`equilibrium constants notebook <EqConst>`. 
The calibration range for these parametrisations is not listed in the original references but is suitable for magmatic systems.
Currently there is only one option for each equilibrium constant (except H2S), but more could be added in the future. 
The model type is in **bold** (i.e., what option will effect this variable) and the options of models to use are in *italics* (default options are indicated).

    KHOg: Model for the parameterisation of the equilibiurm constant for H2 + 0.5O2 = H2O.
        default: 'Ohmoto97' Reaction (d) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Only one option available currently, included for future development.

    KHOSg: Model for the parameterisation of the equilibiurm constant for 0.5S2 + H2O = H2S + 0.5O2.
        default: 'Ohmoto97' Reaction (h) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Other options:
        'no H2S' Stops H2S forming in the vapor (K = 0).
    
    KOSg: Model for the parameterisation of the equilibiurm constant for 0.5S2 + O2 = SO2.
        default: 'Ohmoto97' Reaction (f) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Other options:
        'no SO2' Stops SO2 forming in the vapor (K = 0). As a by-product, OCS will also stop forming.

    KOSg2: Model for the parameterisation of the equilibiurm constant for 0.5S2 + 1.5O2 = SO3.
        default: 'ONeill2' Eq. (6b) from O'Neill & Mavrogenes (2022) GCA 334:368-382 doi:10.1016/j.gca.2022.06.020
        Only one option available currently, included for future development.

    KOCg: Model for the parameterisation of the equilibiurm constant for CO + 0.5O2 = CO2.
        default: 'Ohmoto97' Reaction (c) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Only one option available currently, included for future development. 

    KCOHg: Model for the parameterisation of the equilibiurm constant for CH4 + 2O2 = CO2 + 2H2O.
        default: 'Ohmoto97' Reaction (e) in Table 1 from Ohmoto & Kerrick (1977) AmJSci 277:1013-1044
        Other options:
        'no CH4' Stops CH4 forming in the vapor (K = large number).

    KOCSg: Model for the parameterisation of the equilibiurm constant for OCS.
        default: 'Moussallam19' Eq. (8) for 2CO2 + OCS ⇄ 3CO + SO2 in Moussallam et al. (2019) EPSL 520:260-267 doi:10.1016/j.epsl.2019.05.036 for 
        Other options:
        'no OCS' Stops OCS forming in the vapor (K = large number).  

    KCOs: Model for the parameterisation of the equilibiurm constant for Cgrahite + O2 = CO2.
        default: 'Holloway92' Eq. (3) KI in Holloway et al. (1992) EuropeanJ.Mineralogy 4(1):105-114.
        Only one option available currently, included for future development.

    carbonylsulfide: Reaction equilibrium KOCSg is for. 
        default: 'COS' 2CO2 + OCS ⇄ 3CO + SO2
        Only one option available currently, included for future development.


**KHOg** 

Equilibrium constant for homogeneous vapor equilibrium reaction 0.5H2 + O2 = H2O

- *Ohmoto97*: Reaction (d) in Table 1 of Ohmoto & Kerrick (1997) [default]


**KCOg** 

Equilibrium constant for homogeneous vapor equilibrium reaction CO + 0.5O2 = CO2

- *Ohmoto97*: Reaction (c) in Table 1 of Ohmoto & Kerrick (1997) [default]


**KOSg** 

Equilibrium constant for homogeneous vapor equilibrium reaction 0.5S2 + O2 = SO2

- *Ohmoto97*: Reaction (f) in Table 1 of Ohmoto & Kerrick (1997) [default]


**KCOHg** 

Equilibrium constant for homogeneous vapor equilibrium reaction CH4 + 2O2 = CO2 + 2H2O

- *Ohmoto97*: Reaction (e) in Table 1 of Ohmoto & Kerrick (1997) [default]


**KHOSg**

Equilibrium constant for homogeneous vapor equilibrium reaction 0.5S2 + H2O = H2S + 0.5O2

- *Ohmoto97*: Reaction (h) in Table 1 of Ohmoto & Kerrick (1997) [default]

- *noH2S*: *K* = 0, which prevents H2S forming in the vapor


**KOCSg**

Equilibrium constant for homogeneous vapor equilibrium reaction 2CO2 + OCS = 3CO + SO2

- *Moussallam19*: Eq. (8) in Moussallam et al. (2019) [default]


**KOSg2**

Equilibrium constant for homogeneous vapor equilibrium reaction 0.5S2 + 1.5O2 = SO3

- *ONeill22*: Eq (6b) in O’Neill and Mavrogenes (2022) [default]


**KCOs**

Equilibrium constant for heterogeneous solid-vapor equilibrium reaction Cgraphite + O2 = CO2

- *Holloway92*: Eq. (3) KI in Holloway et al. (1992) [default]