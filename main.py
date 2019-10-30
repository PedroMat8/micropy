# -*- coding: utf-8 -*-
import pore_distribution as pore


def main():

    Gs = 2.65
    Ms = 0.2248  # [g]
    w = 0.47  # [-]
    teta = 147  # [degrees]
    surf_tension = 0.48  # [N/m]
    intervals = 40
    dmax = 200  # [um]
    dmin = 0.004  # [um]

    inputs = {
            'Gs': Gs, 'Ms': Ms, 'w': w, 'teta': teta,
            'surf_tension': surf_tension, 'intervals': intervals,
            'dmax': dmax, 'dmin': dmin}

    data = pore.DataElaboration(inputs)
    data.cpd_from_mip('input/input_mip.txt')
    data.cpd.sort_cpd()
    data.psd_from_cpd()
    data.plot()
    data.save_output()

#    data = pore.DataElaboration(inputs)
#    data.cpd_from_file('input/input_cpd.txt')
#    data.cpd.sort_cpd()
#    data.psd_from_cpd()
#    data.plot()
#    data.save_output()


if __name__ == "__main__":
    main()
