# -*- coding: utf-8 -*-
from micropy import pore_distribution as pore


def main():

    intervals = 40
    dmax = 240  # [um]
    dmin = 0.004  # [um]

    inputs = {'intervals': intervals, 'dmax': dmax, 'dmin': dmin}

    Gs = 2.65
    Ms = 0.2248  # [g]
    w = 0.47  # [-]
    teta = 147  # [degrees]
    surf_tension = 0.48  # [N/m]

    inputs_gtec = {'Gs': Gs, 'Ms': Ms, 'w': w, 'teta': teta,
                   'surf_tension': surf_tension}

    data = pore.DataElaboration(inputs)
    data.cpd_from_mip('input_mip.txt', inputs_gtec)
    data.cpd.sort_cpd()
    data.psd_from_cpd()
    data.plot_mip(inputs_gtec)
    data.save_output()


if __name__ == "__main__":
    main()
