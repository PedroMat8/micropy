"""
pore_distribution module handles cumulative and frequency pore distributions
from MIP and data.

Created on 29/10/2019
MIT License - Copyright (c) 2019 Matteo Pedrotti

@author: PedroMat8
@email: matteo.pedrotti@strath.ac.uk
"""

import numpy as np
import matplotlib.pyplot as plt


class DataElaboration:
    """ Returns PSD and CPD from MIP output"""
    def __init__(self, inputs):
        self.psd = self.PSD()
        self.cpd = self.CPD()
        self.inputs = self.Inputs(inputs)

    def plot_mip(self, inputs_gtec):
        """Plot cpd and psd against expected void ratio"""
        saturated = (input('Is e=wGs? :') or 'No')
        if saturated in ('YES', 'yes', 'Yes', 'Y', 'y'):
            gtec = self.InputsGtec(inputs_gtec)
            e = gtec.w*gtec.Gs
        else:
            e = [float(input('Input void ratio: '))]

        axs = self.plot_data()
        axs[1].semilogx(
                self.cpd.d, ([e] * len(self.cpd.d)))

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

    def save_output(self):
        """Savve outputs"""
        np.savetxt(
                'output.txt',
                np.transpose([self.cpd.d, self.cpd.e, self.psd.d, self.psd.e]),
                header=('diameters_cpd\tvoid_ratio_cpd\tdiameters_psd\t' +
                        'void_ratio_psd'), delimiter='\t')

    def cpd_from_array(self, d, e):
        """Get cpd from array"""
        [d, e] = self.sort_cpd(d, e)

        def delta_e(e, d, d_cpd, idx, i):
            value = ((e[idx]-e[idx-1]) / (
                    np.log10(d[idx])-np.log10(d[idx-1])) * (
                            np.log10(d_cpd[i])-np.log10(d[idx-1])))
            return value

        if len(d) < self.inputs.intervals:
            print('Too many intervals. Reduced to ', len(d))
            self.inputs.intervals = len(d)

        self.cpd.d = np.append(self.cpd.d, self.inputs.dmin)
        self.cpd.e = np.append(self.cpd.e, 0)

        delta_log = ((np.log10(self.inputs.dmax)-np.log10(self.inputs.dmin)) /
                     (self.inputs.intervals-1))

        for i in range(1, self.inputs.intervals):
            self.cpd.d = (np.append(self.cpd.d,
                                    10**(np.log10(self.cpd.d[i-1])+delta_log)))
            for diam in d:
                if diam > self.cpd.d[i]:
                    idx, = np.where(d == diam)
                    break

            incr = delta_e(e, d, self.cpd.d, idx, i)
            if incr < 0:
                print('cdp from_array WARNING: void ratio increment negative!')
                break

            self.cpd.e = (np.append(self.cpd.e, (e[idx-1]+incr)))

    def cpd_from_mip(self, input_file, inputs_gtec):
        """Get cpd from MIP"""
        gtec = self.InputsGtec(inputs_gtec)
        alf = np.loadtxt(input_file, usecols=(0, 1), skiprows=0)
        p = alf[:, 0]*0.00689475908677536  # [MPa]
        v = alf[:, 1]*1000  # [mm3]
        Vs = gtec.Ms/gtec.Gs*1000
        e = v/Vs
        dd = -4*gtec.surf_tension*np.cos(np.radians(gtec.teta))/p
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

        self.cpd_from_array(dd, e)

    def cpd_from_file(self, input_file):
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

        self.cpd_from_array(d, e)

    def norm_dist(self):
        """Get normalized distributions"""
        self.cpd.norm_cpd()
        self.psd.norm_psd(self.cpd.d, self.cpd.e)

    @staticmethod
    def sort_cpd(d, e):
        """Sort cpd"""
        if np.argmax(d) == np.argmax(e):
            d = np.sort(d)
            e = np.sort(e)
            return d, e

        elif d[0] > d[-1]:
            d = np.sort(d)
            e = np.sort(max(e)-e)
            return d, e

        elif e[0] > e[-1]:
            e = (max(e)-e)
            return d, e

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

    class PSD:
        def __init__(self):
            self.d = []
            self.e = []
            self.norm = []

        @staticmethod
        def psd_from_cpd(cpd_d, cpd_e):
            """Create psd from cpd"""
            psd_d = []
            psd_e = []

            for i in range(np.size(cpd_d, 0)-1):
                foo = 10**((np.log10(cpd_d[i]) + np.log10(cpd_d[i+1])) / 2)
                foo1 = (cpd_e[i+1]-cpd_e[i])/(np.log10(cpd_d[i+1]/cpd_d[i]))

                if np.argmax(cpd_e) == 0:  # reverse cpd
                    foo1 = -1*foo1

                psd_d = np.append(psd_d, foo)
                psd_e = np.append(psd_e, foo1)

            return psd_d, psd_e

        def norm_psd(self, cpd_d, cpd_e):
            """Normalize psd"""
            cpd_norm = cpd_e / np.max(cpd_e)
            [foo, psd_norm] = DataElaboration.PSD.psd_from_cpd(cpd_d, cpd_norm)
            self.norm = np.append(self.norm, psd_norm)

    class CPD(PSD):
        def reverse_cpd(self):
            """Reverse cpd"""
            [self.d, self.e] = self.sort_cpd(self.d, self.e)
            self.e = np.sort(max(self.e)-self.e)[::-1]
            return self

        def norm_cpd(self):
            """Normalize cpd"""
            cpd_norm = self.e / np.max(self.e)
            self.norm = np.append(self.norm, cpd_norm)
