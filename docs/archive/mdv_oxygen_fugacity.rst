===================================================================================
Oxygen fugacity and Fe3+/FeT
===================================================================================

Below are the different models for the oxygen fugacity buffers and relationships between fO2 and Fe3+/FeT, which can also be viewed in the :doc:`oxygen fugacity notebook <fO2>`. 
The model type is in **bold** (i.e., what option will effect this variable) and the options of models to use are in *italics* (default options are indicated). 
Note that currently the degassing calculation can only be run using **fO2** = *Kress91A*.


**FMQbuffer** 

Model for the parameterisation for the fO2 value of the FMQ buffer.

- *Frost91*: Frost (1991) [default]

- *ONeill87*: O'Neill (1897)


**NNObuffer** 

Model for the parameterisation for the fO2 value of the NNO buffer.

- *Frost91*: Frost (1991) [default]

- Only one option available currently, included for future development.


**fO2**

Model for parameterisation of relationship between fO2 and Fe3+/FeT.       

- *Kress91A*: Eq. (A-5, A-6) in Kress and Carmichael (1991) [default]

- *Kress91*: Eq. (7) in Kress and Carmichael (1991)

- *ONeill18*: Eq. (9a) in O'Neill et al. (2018)

- *Borisov18*: Eq. (4) from Borisov et al. (2018)