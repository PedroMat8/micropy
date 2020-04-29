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

# =============================================================================
# Required python packages:
# - numpy       -->     conda install numpy
# - matplotlib  -->     conda install matplotlib
# - numba       -->     conda install numba
#               -->     conda install -n micropy cudatoolkit
#
# Keep the input data in a subfolder named "input".
# Input files columns are pore diameter[um] and cumulative void ratio[-].
# No headers.
# =============================================================================


from micropy import pore_distribution as pore

# =============================================================================
# # Input parameters
# =============================================================================
intervals = 35          # Number of intervals for frequency
dmax = 770              # [um] --> Maximum diameter
dmin = 0.0044           # [um] --> Minimum diameter

# =============================================================================
# # There should no need to touch anything below this line
# =============================================================================
inputs = {
        'intervals': intervals,
        'dmax': dmax, 'dmin': dmin
        }

data = pore.DataElaboration(inputs)
data.get_cpd_from_file('input\input-edo2.txt')
[data.psd.d, data.psd.e] = data.psd.get_psd_from_cpd(data.cpd.d, data.cpd.e)
# data.norm_dist()

# Plot cpd and psd
data.plot_data(data.cpd.d, data.cpd.e, data.psd.d, data.psd.e)

# Save data
data.save_output('output')
