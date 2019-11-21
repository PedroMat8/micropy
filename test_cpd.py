# -*- coding: utf-8 -*-
from micropy import pore_distribution as pore
import matplotlib.pyplot as plt


def main():

#    Gs = 2.65
#    Ms = 0.2248  # [g]
#    w = 0.47  # [-]
#    teta = 147  # [degrees]
#    surf_tension = 0.48  # [N/m]
    intervals = 35
    dmax = 770  # [um]
    dmin = 0.0044  # [um]

    inputs = {
            'intervals': intervals,
            'dmax': dmax, 'dmin': dmin
            }

    data = pore.DataElaboration(inputs)
    data.cpd_from_file('input\input-edo2.txt')
#    data.cpd.sort_cpd()
    [data.psd.d, data.psd.e] = data.psd.psd_from_cpd(data.cpd.d, data.cpd.e)
#    data.plot_data()
#    data.save_output()
    data.norm_dist()

    data.plot_data(data.cpd.d, data.cpd.e, data.psd.d, data.psd.e)
    data.plot_data(data.cpd.d, data.cpd.norm, data.psd.d, data.psd.norm)

#    fig, axs = plt.subplots(2)
#    fig.suptitle('Distributions')
#    axs[0].semilogx(data.psd.d, e_norm)
#    axs[0].semilogx(data.psd.d, data.psd.e, '-r')
#    axs[0].set(xlabel='diameters [um]', ylabel='norm')
#    axs[1].semilogx(data.psd.d, e_frequency)
#    axs[1].set(xlabel='diameters [um]', ylabel='frequency')


if __name__ == "__main__":
    main()
