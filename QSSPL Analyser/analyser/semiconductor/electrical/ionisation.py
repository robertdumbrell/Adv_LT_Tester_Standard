#!/usr/local/bin/python
# UTF-8

import numpy as np
import matplotlib.pylab as plt
import os

from semiconductor.helper.helper import HelperFunctions
from . import impurity_ionisation_models as IIm
from semiconductor.material.densityofstates import DOS
import semiconductor.general_functions.carrierfunctions as CF


class Ionisation(HelperFunctions):

    '''
    Depending on a dopant level from a band, and the thermal
    energy available, a dopant is electrical active (donating or
    accepting an electron to the band)  or inactive.
    '''
    _cal_dts = {
        'material': 'Si',
        'temp': 300,
        'author': None,
        'ni_author': None
    }

    author_list = 'ionisation.models'

    def __init__(self, **kwargs):

        # update any values in cal_dts
        # that are passed
        self.calculationdetails = kwargs
        self._init_links()

        # get the address of the authors list
        author_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            self._cal_dts['material'],
            self.author_list)

        # get the models ready
        self._int_model(author_file)

        # initiate the first model
        self.change_model(self._cal_dts['author'])

    def _init_links(self):

        self.Dos = DOS(material=self._cal_dts['material'],
                       temp=self._cal_dts['temp'],
                       author=None
                       )

    def update(self, N_imp, ne, nh, impurity, **kwargs):
        '''
        Calculates the number of ionisied impurities

        Inputs:
            N_imp: (float or numpy array)
                number of the impurity in the material
            ne: (float or numpy array)
                number of electrons
            nh: (float or numpy array)
               number of holes
            impurity: (str)
                the element name of the impurity
            temp: (optional)
                the  temperature in Kelvin to be evaluated
            author: (optional str)
                the author of the impurity model to use

        output:
            the number of ionised impurities
        '''
        self.calculationdetails = kwargs

        # a check to make sure the model hasn't changed
        if 'author' in kwargs.keys():
            self.change_model(self._cal_dts['author'])

        # checks if and get the required density of states model
        if 'dos_author' in self.vals.keys():
            Nc, Nv = self.Dos.update(material=self._cal_dts['material'],
                                     temp=self._cal_dts['temp'],
                                     author=self.vals['dos_author']
                                     )
        else:
            Nc, Nv = 0, 0

        if impurity in self.vals.keys():
            # get the ionisation fraction
            iN_imp = getattr(IIm, self.model)(
                self.vals, N_imp, ne, nh,
                self._cal_dts['temp'], Nc, Nv,
                self.vals[impurity])
            # multiply it by the number of dopants
            iN_imp *= N_imp
        else:
            print('''\nWarning:\n\t'''
                  '''No such impurity, please check your model'''
                  '''and spelling.\n\tReturning zero array\n''')
            iN_imp = np.zeros(np.asarray(N_imp).flatten().shape[0])

        return iN_imp

    def update_dopant_ionisation(self, N_dop, nxc, impurity, **kwargs):
        '''
        This is a special function used to determine the number of
        ionised dopants given a number of excess carriers, and a
        single dopant type.

        inputs:
            N_dop: (float; cm^-3)
                The dopant density
            nxc: (float; cm^-3)
                The excess carrier density
            impurity: (str)
                The name of the dopant used e.g. boron, phosphorous. The
                dopants available depend on the model used

        output:
            N_idop: (float cm^-2)
                The number of ionised dopants
        '''

        self.calculationdetails = kwargs

        # a check to make sure the model hasn't changed
        if 'author' in kwargs.keys():
            self.change_model(self._cal_dts['author'])

        N_idop = N_dop

        if not isinstance(N_idop, np.ndarray):
            N_idop = np.asarray([N_idop])

        if not isinstance(N_dop, np.ndarray):
            N_dop = np.asarray([N_dop])

        if impurity in self.vals.keys():
            # TO DO, change this from just running 10 times to a proper check
            for i in range(10):
                if self.vals['tpe_' + self.vals[impurity]] == 'donor':
                    Nd = np.copy(N_idop)
                    Na = np.zeros(Nd.shape)
                elif self.vals['tpe_' + self.vals[impurity]] == 'acceptor':
                    Na = np.copy(N_idop)
                    Nd = np.zeros(Na.shape)
                else:
                    print('something went wrong in ionisation model')

                ne, nh = CF.get_carriers(
                    Na,
                    Nd,
                    nxc,
                    temp=self._cal_dts['temp'],
                    material=self._cal_dts['material'],
                    ni=self._cal_dts['ni_author'])

                N_idop = self.update(
                    N_dop, ne, nh, impurity)

        else:
            print(r'Not a valid impurity, returning 100% ionisation')

        return N_idop

    def check_models(self):
        '''
        Plots a check of the modeled data against Digitised data from either
        papers or from other implementations of the model.
        '''
        plt.figure('Ionised impurities')

        iN_imp = N_imp = np.logspace(15, 20)

        dn = 1e10
        temp = 300

        for impurity in ['phosphorous', 'arsenic']:

            iN_imp = self.update_dopant_ionisation(N_imp,
                                                   dn,
                                                   impurity,
                                                   temp, author=None)

            if not np.all(iN_imp == 0):
                plt.plot(
                    N_imp, iN_imp / N_imp * 100, label='Altermatt: ' + impurity)

        test_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'Si', 'check data', 'donors.csv')

        data = np.genfromtxt(test_file, delimiter=',', skip_header=1)

        for i in range(0, (data.shape[1] + 2) / 2, 2):
            # print i
            plt.plot(
                data[:, i], data[:, i + 1] * 100, 'r.',
                label='Digitised data')

        plt.semilogx()
        plt.xlabel('Impurity (cm$^{-3}$)')
        plt.ylabel('Fraction of Ionised impurities (%)')

        plt.legend(loc=0)
