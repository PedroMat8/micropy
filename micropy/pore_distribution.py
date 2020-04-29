"""
pore_distribution module handles cumulative and frequency pore distributions
from MIP and data.

Created on 29/10/2019
MIT License - Copyright (c) 2019 Matteo Pedrotti

@author: PedroMat8
@email: matteo.pedrotti@strath.ac.uk

Please cite the following DOI: 10.5281/zenodo.3524929
GitHub repository: https://github.com/PedroMat8/micropy
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from numba import jit


class DataElaboration():
    """ Stores psd and cpd of the data set"""
    def __init__(self, inputs):
        self.inputs = self.Inputs(inputs)
        self.psd = self.PSD(self.inputs.intervals)
        self.cpd = self.CPD(self.inputs.intervals)


    def plot_mip(self, inputs_gtec):
        """Plot cpd and psd against expected void ratio"""
        saturated = (input('Is e=wGs? [Y or N]: '))
        if saturated in ('YES', 'yes', 'Yes', 'Y', 'y'):
            gtec = self.InputsGtec(inputs_gtec)
            e = gtec.w*gtec.Gs
        else:
            e = [float(input('Input void ratio: '))]

        cpd_d = self.cpd.d
        cpd_e = self.cpd.e
        psd_d = self.psd.d
        psd_e = self.psd.e

# TODO: This part below partially replicates method plot_data"
        fig, axs = plt.subplots(2)
        fig.suptitle('Void ratio comparison')
        axs[0].semilogx(cpd_d, cpd_e)
        axs[0].semilogx(cpd_d, [e] * len(cpd_d))
        axs[0].set(xlabel='diameters [um]', ylabel='void ratio [e]')
        axs[1].semilogx(psd_d, psd_e)
        axs[1].set(xlabel='diameters [um]', ylabel='frequency de/d(logd)')


    @staticmethod
    def plot_data(cpd_d, cpd_e, psd_d, psd_e):
        """Plot cpd and psd"""
        fig, axs = plt.subplots(2)
        fig.suptitle('Distributions')
        axs[0].semilogx(psd_d, psd_e)
        axs[0].set(xlabel='diameters [um]', ylabel='frequency de/d(logd)')
        axs[1].semilogx(cpd_d, cpd_e)
        axs[1].set(xlabel='diameters [um]', ylabel='void ratio [e]')
        return fig, axs

    def save_output(self, name='output'):
        """Save outputs"""
        matrix = np.column_stack((self.cpd.d, self.cpd.e,
                                  np.append(self.psd.d,0),
                                  np.append(self.psd.e,0)))
        header=('diameters_cpd\tvoid_ratio_cpd\tdiameters_psd\tvoid_ratio_psd')
        namefile = name + '.txt'
        np.savetxt(namefile, matrix,
                   header=header, delimiter='\t', fmt='%s')

    @staticmethod
    def interpolate_e(d_min, d_max, d_starting, e_starting, intervals):

        delta_log = ((np.log10(d_max)-np.log10(d_min)) / (intervals-1))
        d_starting_log = np.log10(d_starting)
        e_new = np.empty(intervals)

        n = np.arange(0, intervals)
        ndlog = np.multiply(n ,delta_log)
        d_new_log = np.add(np.log10(d_min), ndlog)
        d_new = 10**d_new_log
        e_new = np.interp(d_new_log, d_starting_log, e_starting)

        return d_new, e_new

    def get_cpd_from_array(self, d, e):
        """Get cpd from array"""
        [d, e] = self.sort_cpd(d, e)
        [self.cpd.d, self.cpd.e] = DataElaboration.interpolate_e(
            self.inputs.dmin, self.inputs.dmax, d, e, self.inputs.intervals)

    def get_cpd_from_mip(self, input_file, inputs_gtec):
        """Get cpd from MIP"""
        gtec = self.InputsGtec(inputs_gtec)
        alf = np.loadtxt(input_file, usecols=(0, 1), skiprows=0)
        p = alf[:, 0]*0.00689475908677536  # [MPa]
        v = alf[:, 1]*1000  # [mm3]
        Vs = gtec.Ms/gtec.Gs*1000
        e = v/Vs
         # 2 for parallel plates and 4 for axial symmetry
        dd = -2*0.48*np.cos(np.radians(147))/p
        if len(dd) < self.inputs.intervals:
            print('Too many intervals. Reduced to ', len(dd))
            self.inputs.intervals = len(dd)

        print('Max available diameter [um]: ', np.max(dd))
        print('Min available diameter [um]: ', np.min(dd))

        if self.inputs.dmax > np.max(dd):
            new_dmax = input('dmax is too large, input new diameter < ' + str(
                    round(np.max(dd))) + ': ')
            self.inputs.dmax = round(float(new_dmax))

        if self.inputs.dmin < np.min(dd):
            new_dmin = input('dmin is too small, input new diameter > ' + str(
                    round(np.min(dd), 4)) + ': ')
            self.inputs.dmin = round(float(new_dmin), 4)

        self.get_cpd_from_array(dd, e)

    def get_cpd_from_file(self, input_file):
        """Get cpd from txt file"""
        alf = np.loadtxt(input_file, usecols=(0, 1), skiprows=0)
        d = alf[:, 0]
        e = alf[:, 1]

        print('Max available diameter [um]: ', np.max(d))
        print('Min available diameter [um]: ', np.min(d))
        if len(d) < self.inputs.intervals:
            print('Too many intervals. Reduced to ', len(d))
            self.inputs.intervals = len(d)

        if self.inputs.dmax > np.max(d):
            new_dmax = input('dmax is too large, input new diameter < ' + str(
                    round(np.max(d))) + ': ')
            self.inputs.dmax = round(float(new_dmax))

        if self.inputs.dmin < np.min(d):
            new_dmin = input('dmin is too small, input new diameter > ' + str(
                    round(np.min(d), 4)) + ': ')
            self.inputs.dmin = round(float(new_dmin), 4)

        self.get_cpd_from_array(d, e)

    @staticmethod
    def sort_cpd(d, e):
        """Sort cpd"""
        if d[0] >= d[-1]:
            if e[0] >= e[-1]:
                return np.sort(d), np.sort(e)
            else:
                return np.sort(d), np.sort(max(e)-e)
        else:
            if e[0] <= e[-1]:
                return np.sort(d), np.sort(e)
            else:
                return np.sort(d), np.sort(max(e)-e)

    class Inputs:
        def __init__(self, inputs):
            self.intervals = inputs['intervals']
            self.dmax = inputs['dmax']
            self.dmin = inputs['dmin']

    class InputsGtec:
        def __init__(self, inputs):
            self.Gs = inputs['Gs']
            self.Ms = inputs['Ms']
            self.w = inputs['w']
            self.teta = inputs['teta']
            self.surf_tension = inputs['surf_tension']

    class CPD():
        def __init__(self, intervals):
            dim = intervals
            self.d = np.empty(dim)
            self.e = np.empty(dim)

    class PSD(CPD):
        def __init__(self, intervals):
            dim = intervals
            self.d = np.empty(dim)
            self.e = np.empty(dim)
            # self.norm = np.empty(dim)

        @staticmethod
        def get_psd_from_cpd(cpd_d, cpd_e):
            """Create psd from cpd"""
            [cpd_d, cpd_e] = DataElaboration.sort_cpd(cpd_d, cpd_e)

            alf = np.size(cpd_d, 0)
            psd_d = np.empty(alf-1)
            psd_e = np.empty(alf-1)

            for i in range(alf-1):
                foo = 10**((np.log10(cpd_d[i]) + np.log10(cpd_d[i+1])) / 2)
                foo1 = (cpd_e[i+1]-cpd_e[i])/(np.log10(cpd_d[i+1]/cpd_d[i]))

                psd_d[i] = foo
                psd_e[i] = foo1
                if foo1 < 0:
                    sys.exit('negative psd')

            psd_e/= np.sum(psd_e)

            return psd_d, psd_e
