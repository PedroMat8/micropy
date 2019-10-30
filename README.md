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
    * 2nd column --> void ratio or porosity [-]
                                (assumes it is void ratio)

* Inputs are in class INPUTS:
  * Gs --> specific gravity
  * Ms --> MIP sample dry mass
  * w --> water content
  * teta --> contact angle [degres]
  * surf_tension --> surface tension [N/m]
  * intervals --> psd intervals
  * dmax --> distribution max limit [um]
  * dmin -->        ** distributions min limit [um]

Methods:
    plot                    -->     plots psd and cpd
    save_output             -->     saves psd and cpd on txt
    psd_from_cpd            -->     computes psd
    cpd_from_mip            -->     import cpd from MIP
    cpd_from_file           -->     import cpd from txt file

Classes:
    Inputs                  -->     input parameters
    PSD                     -->     psd
    CPD                     -->     cpd
