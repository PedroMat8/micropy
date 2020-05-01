# -*- coding: utf-8 -*-
"""
Timing function interpolate in micropy against interpolate in scipy

To run it, the funcrion get_cpd_from_file in micropy needs to:
    return d, e

@author: Matteo
"""

import pore_distribution as dist
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

import timeit

# =============================================================================
# # Prepare data
# =============================================================================
intervals = 60          # Number of intervals for frequency
dmax = 770              # [um] --> Maximum diameter
dmin = 0.0044           # [um] --> Minimum diameter

inputs = {
        'intervals': intervals,
        'dmax': dmax, 'dmin': dmin
        }

data = dist.DataElaboration(inputs)
d, e = data.get_cpd_from_file('..\input\input-edo2.txt')

[d, e] = data.sort_cpd(d, e)

# =============================================================================
# # Funtions to time
# =============================================================================
def micropy(dmin, dmax, d, e, intervals):
    d1, e1 =  data.interpolate_e(dmin, dmax, d, e, intervals)
    return d1, e1

def sci_interp(dmin, dmax, d ,e, intervals):
    delta_log = ((np.log10(dmax)-np.log10(dmin)) / (intervals-1))
    dd = np.log10(d)
    f = interpolate.interp1d(dd,e,kind='cubic')
    d_new = np.arange(np.min(dd), np.max(dd), delta_log)
    e_new = f(d_new)
    return 10**d_new, e_new

# =============================================================================
# # Timeit
# =============================================================================
number = 1000
repeat = 4
time_list = timeit.repeat(
        setup='',
        stmt='[micropy(data.inputs.dmin, data.inputs.dmax, d, e, data.inputs.intervals)]',
        number= number,
        repeat= repeat,
        globals=globals())

time_list2 = timeit.repeat(
        setup='',
        stmt='[sci_interp(data.inputs.dmin, data.inputs.dmax,d,e,data.inputs.intervals)]',
        number= number,
        repeat= repeat,
        globals=globals())

print('interp: ' + str(min(time_list)/number) + ' s')
print('sci_interp: ' + str(min(time_list2)/number) + ' s')

# =============================================================================
# # Plot
# =============================================================================
d1, e1 = micropy(data.inputs.dmin, data.inputs.dmax, d, e, data.inputs.intervals)
d2, e2 = sci_interp(data.inputs.dmin, data.inputs.dmax, d,e,data.inputs.intervals)
plt.semilogx(d,e, label='input')
plt.semilogx(d1,e1, label='micropy')
plt.semilogx(d2,e2, label='scipy')
plt.legend()

