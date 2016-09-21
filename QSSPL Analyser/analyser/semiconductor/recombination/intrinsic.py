
import numpy as np
import matplotlib.pylab as plt
import os
import configparser

from semiconductor.helper.helper import HelperFunctions, change_model
from semiconductor.general_functions.carrierfunctions import get_carriers
import semiconductor.recombination.radiative_models as radmdls
import semiconductor.recombination.auger_models as augmdls


class Intrinsic(HelperFunctions):

    _cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'ni_author': None,
        'rad_author': None,
        'aug_author': None,
        'Na': 1,
        'Nd': 1e16,
    }

    def __init__(self, **kwargs):
        # update any values in cal_dts
        # that are passed
        self.calculationdetails = kwargs

        # pass values to models
        self._update_links()

    def _update_links(self):

        self.Radiative = Radiative(
            material=self._cal_dts['material'],
            author=self._cal_dts['rad_author'],
            temp=self._cal_dts['temp'],
            ni_author=self._cal_dts['ni_author'],
            Na=self._cal_dts['Na'],
            Nd=self._cal_dts['Nd'],
        )

        self.Auger = Auger(
            material=self._cal_dts['material'],
            author=self._cal_dts['aug_author'],
            temp=self._cal_dts['temp'],
            ni_author=self._cal_dts['ni_author'],
            Na=self._cal_dts['Na'],
            Nd=self._cal_dts['Nd'],
        )

    def tau(self, nxc, **kwargs):
        '''
        Returns the intrinsic carrier lifetime
        '''
        return 1. / self.itau(nxc, **kwargs)

    def itau(self, nxc, **kwargs):
        '''
        Returns the inverse of the intrinsic carrier lifetime
        '''
        self.calculationdetails = kwargs
        if 'author' in ''.join(kwargs.keys()):
            self._update_links()

        itau = self.Radiative.itau(nxc, **kwargs) +\
            self.Auger.itau(nxc, **kwargs)

        return itau


class Radiative(HelperFunctions):

    author_list = 'radiative.model'

    _cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'author': None,
        'ni_author': None,
        'Na': 1,
        'Nd': 1e16,
        'nxc': 1e12
    }

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

    def tau(self, nxc, **kwargs):
        self.calculationdetails = kwargs
        self.change_model(self._cal_dts['author'])

        ne0, nh0 = get_carriers(
            Na=self._cal_dts['Na'],
            Nd=self._cal_dts['Nd'],
            nxc=self._cal_dts['nxc'],
            ni_author=self._cal_dts['ni_author'],
            temp=self._cal_dts['temp']
        )

        B = self._get_B()

        return getattr(radmdls, self.model)(
            self.vals, nxc, nh0, ne0, B, temp=self._cal_dts['temp']
        )

    def itau(self, nxc, **kwargs):
        return 1. / self.tau(nxc, **kwargs)

    def _get_Blow(self):

        # if there is a model for blow, apply it
        if 'blow_model' in self.vals.keys():
            vals, model = change_model(self.Models, self.vals['blow_vals'])

            B = getattr(radmdls, self.vals['blow_model'])(
                vals=vals, temp=self._cal_dts['temp'])

        # else use the constant value
        else:
            B = self.vals['b']

        return B

    def _get_B(self, **kwargs):
        self.calculationdetails = kwargs
        # if there is a model for blow, apply it
        if 'doping' not in kwargs.keys():
            kwargs['doping'] = abs(self._cal_dts['Na'] -
                                   self._cal_dts['Nd'])

        if 'b_model' in self.vals.keys():
            B = getattr(radmdls, self.vals['b_model'])(
                vals=self.vals, temp=self._cal_dts['temp'],
                nxc=self._cal_dts['nxc'],
                doping=kwargs['doping'],
                Blow=self._get_Blow())
        else:
            B = self._get_Blow()

        return B


class Auger(HelperFunctions):
    author_list = 'auger.model'

    _cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'author': None,
        'ni_author': None,
        'Na': 1,
        'Nd': 1e16,
    }

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

    def tau(self, nxc, **kwargs):
        self.calculationdetails = kwargs

        if 'author' in kwargs.keys():
            self.change_model(self._cal_dts['author'])

        ne0, nh0 = get_carriers(
            Na=self._cal_dts['Na'],
            Nd=self._cal_dts['Nd'],
            nxc=0,
            ni_author=self._cal_dts['ni_author'],
            temp=self._cal_dts['temp']
        )

        return getattr(augmdls, self.model)(
            self.vals, nxc, ne0, nh0, temp=self._cal_dts['temp'])

    def itau(self, nxc, **kwargs):
        return 1. / self.tau(nxc, **kwargs)

    def check(self, author, fig=None, ax=None):
        if ax is None:
            fig, ax = plt.subplots(1)
        self.change_model(author, self.Models)

        func = getattr(augmdls, self.model)

        getattr(augmdls, author + '_check')(self.vals, func, fig, ax)
        ax.set_xlim(left=1e13)
