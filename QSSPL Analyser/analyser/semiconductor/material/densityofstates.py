#!/usr/local/bin/python
# UTF-8

import numpy as np
import matplotlib.pylab as plt
import sys
import os

from semiconductor.helper.helper import HelperFunctions
import semiconductor.material.densityofstates_models as dos_models
from semiconductor.material.bandgap_intrinsic import IntrinsicBandGap as Egi


class DOS(HelperFunctions):

    '''
    The density of states is a value that determines the
    number of free states for electrons and holes in the conduction
    and valance band
    '''
    _cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'author': None,
        'iEg_author': None
    }
    author_list = 'DOS.models'

    def __init__(self, **kwargs):

        # update any values in cal_dts
        # that are passed
        self.calculationdetails = kwargs

        # get the address of the authors list
        author_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            self._cal_dts['material'],
            self.author_list)

        # get the models ready
        self._int_model(author_file)

        # initiate the first model
        self.change_model(self._cal_dts['author'])

    def update(self, **kwargs):
        '''
        a function to update the density of states

        inputs:
            dictionary with values found in self._cal_dts:
                temperature: (optional)
                             in kelvin
                author:  (optional)
                    the author used.
                    If not provided the last provided author is used
                    If no author has been provided,  Couderc's model is used
        output:
            the density of states of the valance band
            the density of states of the conduction band
        '''
        self.calculationdetails = kwargs

        # a check to make sure the model hasn't changed
        if 'author' in kwargs.keys():
            self.change_model(self._cal_dts['author'])

        if 'ieg_author' in self.vals.keys():

            Eg0 = Egi(
                material=self._cal_dts['material'],
                temp=0,
                author=self.vals['ieg_author'],
            ).update()
            Egratio = Eg0 / Egi(
                material=self._cal_dts['material'],
                temp=self._cal_dts['temp'],
                author=self.vals['ieg_author'],
            ).update()

        else:
            Egratio = None

        self.Nc, self.Nv = getattr(
            dos_models, self.model
        )(
            self.vals,
            temp=self._cal_dts['temp'],
            Egratio=Egratio)

        return self.Nc, self.Nv

    def check_models(self):
        temp = np.logspace(0, np.log10(600))
        num = len(self.available_models())

        fig, ax = plt.subplots(1)
        self.plotting_colours(num, fig, ax, repeats=2)

        for author in self.available_models():
            Nc, Nv = self.update(temp=temp, author=author)
            # print Nc.shape, Nv.shape, temp.shape
            ax.plot(temp, Nc, '--')
            ax.plot(temp, Nv, '.', label=author)

        ax.loglog()
        leg1 = ax.legend(loc=0, title='colour legend')

        Nc, = ax.plot(np.inf, np.inf, 'k--', label='Nc')
        Nv, = ax.plot(np.inf, np.inf, 'k.', label='Nv')

        plt.legend([Nc, Nv], ['Nc', 'Nv'], loc=4, title='Line legend')
        plt.gca().add_artist(leg1)

        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('Density of states (cm$^{-3}$)')
        plt.show()
