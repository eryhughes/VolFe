===================================================================================
Equilibrium constants
===================================================================================

Below are the different models for the equilibrium constants currently in VolFe, which can also be viewed in the :doc:`equilibrium constants notebook <EqConst>`. 
Currently there is only one option for each equilibrium constant (except H2S), but more could be added in the future. 
The name of the option in VolFe is in **Bold**.

0.5H2 + O2 = H2O
-----

**Ohmoto97:** Reaction (d) in Table 1 of Ohmoto & Kerrick (1997)

.. math:: 12510/T - 0.979log10T + 0.483


CO + 0.5O2 = CO2
-------

**Ohmoto97:** Reaction (c) in Table 1 of Ohmoto & Kerrick (1997)

.. math:: 14751/T - 4.535

0.5S2 + O2 = SO2
-------

**Ohmoto97:** Reaction (f) in Table 1 of Ohmoto & Kerrick (1997)

.. math:: 18929/T - 3.783

CH4 + 2O2 = CO2 + 2H2O
--------

**Ohmoto97:** Reaction (e) in Table 1 of Ohmoto & Kerrick (1997)

.. math:: 12510/T - 0.979log10T + 0.483


0.5S2 + H2O = H2S + 0.5O2
--------

**Ohmoto97:** Reaction (h) in Table 1 of Ohmoto & Kerrick (1997)

**noH2S:** *K* = 0, which prevents H2S forming in the vapor


2CO2 + OCS = 3CO + SO2
------

**Moussallam19:** Eq. (8) in Moussallam et al. (2019)


0.5S2 + 1.5O2 = SO3
----

**ONeill22:** Eq (6b) in O’Neill and Mavrogenes (2022)


Cgraphite + O2 = CO2
----

**Holloway92:** Eq. (3) KI in Holloway et al. (1992)