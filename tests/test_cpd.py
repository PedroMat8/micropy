# -*- coding: utf-8 -*-
from micropy import pore_distribution as pore


def main():

#    Gs = 2.65
#    Ms = 0.2248  # [g]
#    w = 0.47  # [-]
#    teta = 147  # [degrees]
#    surf_tension = 0.48  # [N/m]
    intervals = 40
    dmax = 200  # [um]
    dmin = 0.004  # [um]

    inputs = {
            'intervals': intervals,
            'dmax': dmax, 'dmin': dmin
            }

    data = pore.DataElaboration(inputs)
    data.cpd_from_file('input_cpd.txt')
    data.cpd.sort_cpd()
    data.psd_from_cpd()
    data.plot_data()
    data.save_output()


if __name__ == "__main__":
    main()
