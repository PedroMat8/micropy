# -*- coding: utf-8 -*-

''' Please cite the following DOI:
    DOI:10.5281/zenodo.3525292

    GitHub repository: https://github.com/PedroMat8/micropy

    micropy package installation :
    pip install -i https://test.pypi.org/simple/ micropy'''

''' Created on 29/10/2019
    MIT License - Copyright (c) 2019 Matteo Pedrotti

    @author: PedroMat8
    @email: matteo.pedrotti@strath.ac.uk '''

# Required python packages:
# - numpy       -->     pip install numpy
# - matplotlib  -->     pip install matplotlib
# - sys         -->     pip install os-sys
# - numba       -->     pip install numba

# Keep the input data in a subfolder named "input".
# Input columns are MIP intrusion pressure [PSI] and intruded volume [cc].
# No headers.
##############################################################################
##############################################################################

from micropy import pore_distribution as pore

outputfilename = 'output.txt'
inputfilename = 'input\input_mip.txt'

# Input parameters
Gs = 2.65               # of the clay
Ms = 0.2248             # [g] --> of the specimen
w = 0.47                # [-] --> of the specimen
teta = 147              # [degrees] --> Mercury contact angle
surf_tension = 0.48     # [N/m] --> mercury surface tension

# Distribution
intervals = 35          # Number of intervals for frequency
dmax = 120              # [um] --> Maximum diameter
dmin = 0.004            # [um] --> Minimum diameter

# There should no need to touch anything below this line
inputs = {
        'intervals': intervals,
        'dmax': dmax, 'dmin': dmin
        }
inputs_gtec = {'Gs': Gs, 'Ms': Ms, 'w': w, 'teta': teta,
               'surf_tension': surf_tension}

data = pore.DataElaboration(inputs) # Create a data set
data.cpd_from_mip(inputfilename, inputs_gtec) # Create the cpd

# Calculate the psd
[data.psd.d, data.psd.e] = data.psd.psd_from_cpd(data.cpd.d, data.cpd.e)

# Plot cpd and psd
data.plot_mip(inputs_gtec)

# Save data
data.save_output()
