===================================================================================
Fugacity coefficients
===================================================================================

Below are the different models for the fugacity coefficients currently in VolFe, which can also be viewed in the :doc:`fugacity coefficients notebook <FugCoeff>`. 
The model type is in **bold** (i.e., what option will effect this variable) and the options of models to use are in *italics* (default options are indicated). 
VolFe treats the vapor as an ideal mixture of non-ideal gases; hence, these parameterisations depend no *P* and *T* but *not* the vapor composition.
The range of valid *P* and *T* is stated for each model when known.

Note that all vapor species can be treated as ideal (i.e., their fugacity coeffient = 1 at all *P*) by setting **ideal_gas** to *yes*, which overides the individual options listed below.


**y_O2** 

- *Shi92*: Shi & Saxena (1992), valid between 51–20000 bar and -118–2227 'C [default]

- *ideal*: *y* = 1 at all *P*


**y_H2** 

- *Shaw64*: Shaw & Wones (1964), valid between 1–3040 bar and 0–1000 'C (Note: used value of -0.011901 instead of -0.11901 as reported to match their Table 2) [default] 

- *ideal*: *y* = 1 at all *P*


**y_S2** 

- *Shi92*: Shi & Saxena (1992), valid between 73–20000 bar and -65–2227 'C [default]

- *ideal*: *y* = 1 at all *P*


**y_CO** 

- *Shi92*: Shi & Saxena (1992), valid between 35–20000 bar and -140–2227 'C [default]

- *ideal*: *y* = 1 at all *P*


**y_H2O** 

- *Holloway91*: Holland & Powell (1991), valid between 1–50000 bar and 100–1400 'C  [default]

- *ideal*: *y* = 1 at all *P*


**y_CO2** 

- *Shi92*: Shi & Saxena (1992), valid between 74–20000 bar and 31–2227 'C [default]

- *Holloway91*: Holland & Powell (1991)

- *ideal*: *y* = 1 at all *P*


**y_SO2** 

- *Shi92*: Shi & Saxena (1992), valid between 78–20000 bar and 158–2227 'C

- *Shi92_Hughes24*: Figure S1 from Hughes et al. (2023) based on Shi and Saxena (1992), valid between 78–20000 bar and 158–2227 'C [default]

- *ideal*: *y* = 1 at all *P*


**y_H2S** 

- *Shi92*: Shi & Saxena (1992), valid between 90–20000 bar and 101–2227 'C

- *Shi92_Hughes24*: Figure S1 from Hughes et al. (2024) based on Shi and Saxena (1992), valid between 78–20000 bar and 158–2227 'C [default]

- *ideal*: *y* = 1 at all *P*


**y_CH4** 

- *Shi92*: Shi & Saxena (1992), valid between 46–20000 bar and -82–2227 'C [default]

- *ideal*: *y* = 1 at all *P*


**y_OCS** 

- *Shi92*: Shi & Saxena (1992), valid between 66–20000 bar and 105–2227 'C [default]

- *ideal*: *y* = 1 at all *P*


**y_X** 

- *ideal*: *y* = 1 at all *P* [default]