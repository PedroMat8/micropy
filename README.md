# pore_sizer
Tool to handle pore size distributions.
It returns PSD and CPD from MIP output. CPD is ordered from smaller to larger
pores.

* Input files:
  * input_mip --> is a two columns txt file from MIP analysis
    * 1st column --> Intrusion pressure [PSI]
    * 2nd column --> Intruded volume [cc]

  * input_cpd --> is a two columns txt file from existing cpd
    * 1st column --> diameters [um]
    * 2nd column --> void ratio or porosity [-]. Assumes it is void ratio

* Input parameters are in class Inputs and InputsGtec:
  * Inputs --> geometrical input parameters
    * intervals --> psd intervals
    * dmax --> distribution max limit [um]
    * dmin --> distributions min limit [um]

  * InputsGtec --> physical input parameters
    * Gs --> specific gravity
    * Ms --> MIP sample dry mass
    * w --> water content
    * teta --> contact angle [degres]
    * surf_tension --> surface tension [N/m]

* Methods:
  * plot_data --> plots psd and cpd
  * plot_mip --> plots psd and cpd against expected e
  * save_output --> saves psd and cpd on txt
  * psd_from_cpd --> computes psd
  * cpd_from_mip --> import cpd from MIP
  * cpd_from_file --> import cpd from txt file
  * cpd_from_array --> create cpd from d and e array
  * cpd.sort_cpd --> order cpd from smaller to larger diameters
  * cpd.reverse_cpd -->         reverse cpd order

* Classes:
  * Inputs --> geometrical input parameters
  * InputsGtec --> physical input parameters
  * PSD --> psd
  * CPD --> cpd

  Created on 29/10/2019
  MIT License - Copyright (c) 2019 Matteo Pedrotti
  @author: PedroMat8
  @email: matteo.pedrotti@strath.ac.uk
