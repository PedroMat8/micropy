## Package --> pore_sizer
Tool to investigate pore space.
## Module --> pore_distribution
It handles cumulative and frequency pore distributions from MIP and/or existing
data.

* ### Input files:
  * **input_mip** --> is a two columns txt file from MIP analysis
    * _1st column_ --> Intrusion pressure [PSI]
    * _2nd column_ --> Intruded volume [cc]

  * **input_cpd** --> is a two columns txt file from existing cpd
    * _1st column_ --> diameters [um]
    * _2nd column_ --> void ratio or porosity [-]. Assumes it is void ratio

* ### Input parameters:
  * **Inputs** --> geometrical input parameters
    * _Inputs.intervals_ --> psd intervals
    * _Inputs.dmax_ --> distribution max limit [um]
    * _Inputs.dmin_ --> distributions min limit [um]

  * **InputsGtec** --> physical input parameters
    * _InputsGtec.Gs_ --> specific gravity
    * _InputsGtec.Ms_ --> MIP sample dry mass
    * _InputsGtec.w_ --> water content
    * _InputsGtec.teta_ --> contact angle [degres]
    * _InputsGtec.surf_tension_ --> surface tension [N/m]

* ### Methods:
  * _plot_data_ --> plots psd and cpd
  * _plot_mip_ --> plots psd and cpd against expected e
  * _save_output_ --> saves psd and cpd on txt
  * _psd_from_cpd_ --> computes psd
  * _cpd_from_mip_ --> import cpd from MIP
  * _cpd_from_file_ --> import cpd from txt file
  * _cpd_from_array_ --> create cpd from d and e array
  * _cpd.sort_cpd_ --> order cpd from smaller to larger diameters
  * _cpd.reverse_cpd_ --> reverse cpd order

* ### Classes:
  * _Inputs_ --> geometrical input parameters
  * _InputsGtec_ --> physical input parameters
  * _PSD_ --> psd
  * _CPD_ --> cpd

  Created on 29/10/2019
  MIT License - Copyright (c) 2019 Matteo Pedrotti

  @author: PedroMat8
  @email: matteo.pedrotti@strath.ac.uk
