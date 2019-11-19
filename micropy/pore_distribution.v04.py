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
        print('ciccio')
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

        fig, axs = self.plot_data()
        axs[1].semilogx(
                self.cpd.d, ([e] * len(self.cpd.d)))

    def plot_data(self):
        """Plot cpd and psd"""
        fig, axs = plt.subplots(2)
        fig.suptitle('Distributions')
        axs[0].semilogx(self.psd.d, self.psd.e)
        axs[0].set(xlabel='diameters [um]', ylabel='frequency de/d(logd)')
        axs[1].semilogx(self.cpd.d, self.cpd.e)
        axs[1].set(xlabel='diameters [um]', ylabel='void ratio [e]')
        return fig, axs

    def save_output(self):
        np.savetxt(
                'output.txt',
                np.transpose([self.cpd.d, self.cpd.e, self.psd.d, self.psd.e]),
                header=('diameters_cpd\tvoid_ratio_cpd\tdiameters_psd\t' +
                        'void_ratio_psd'), delimiter='\t')

    def psd_from_cpd(self):
        self.psd.d = np.append(self.psd.d, 0)
        self.psd.e = np.append(self.psd.e, 0)

        for i in range(np.size(self.cpd.d, 0)-1):
            foo = np.exp((np.log(self.cpd.d[i]) +
                          np.log(self.cpd.d[i+1])) / 2)
            foo1 = (self.cpd.e[i+1]-self.cpd.e[i])/(
                            np.log(self.cpd.d[i+1]/self.cpd.d[i]))
            if np.argmax(self.cpd.e) == 0:  # reverse cpd
                foo1 = -1*foo1

            self.psd.d = np.append(self.psd.d, foo)
            self.psd.e = np.append(self.psd.e, foo1)

    def cpd_from_array(self, d, e):

        def delta_e(e, d, d_cpd, idx, i):
            value = (((e[idx]-e[idx-1]) / (np.log(d[idx])-np.log(d[idx-1])) *
                      (np.log(d_cpd[i])-np.log(d[idx]))))
            return value

        if len(d) < self.inputs.intervals:
            print('Too many intervals. Reduced to ', len(d))
            self.inputs.intervals = len(d)

        self.cpd.d = np.append(self.cpd.d, self.inputs.dmin)
        idx = (np.abs(d - self.cpd.d[0])).argmin()
        incr = delta_e(e, d, self.cpd.d, idx, 0)
        self.cpd.e = np.append(self.cpd.e, (e[idx] - incr))

        delta_log = ((np.log(self.inputs.dmax)-np.log(self.inputs.dmin)) /
                     (self.inputs.intervals-1))

        for i in range(1, self.inputs.intervals):
            self.cpd.d = (np.append(self.cpd.d,
                                    np.exp(np.log(self.cpd.d[i-1])+delta_log)))
            idx = (np.abs(d - self.cpd.d[i])).argmin()
            incr = delta_e(e, d, self.cpd.d, idx, i)
            self.cpd.e = (np.append(self.cpd.e, (e[idx]-incr)))

    def cpd_from_mip(self, input_file, inputs_gtec):
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

    class CPD(PSD):
        def sort_cpd(self):
            if np.argmax(self.d) == np.argmax(self.e):
                self.d = np.sort(self.d)
                self.e = np.sort(self.e)
                return self

            elif np.argmax(self.d) == 0:
                self.d = np.sort(self.d)
                self.e = np.sort(max(self.e)-self.e)
                return self

            elif np.argmax(self.e) == 0:
                self.e = np.sort(max(self.e)-self.e)[::-1]
                return self

        def reverse_cpd(self):
            self.sort_cpd()
            self.e = np.sort(max(self.e)-self.e)[::-1]
            return self
