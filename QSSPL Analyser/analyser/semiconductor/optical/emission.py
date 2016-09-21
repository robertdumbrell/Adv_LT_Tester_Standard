import numpy as np
from matplotlib.pylab import *
import sys
import os
import scipy.constants as const

import semiconductor.material.ni as ni
import semiconductor.optical.opticalproperties as opticalproperties
import semiconductor.optical.absorptance as absorptance
from semiconductor.helper.helper import HelperFunctions

import inspect


class SpontaneousRadiativeEmission(HelperFunctions):

    """
    This class calculates the spectral spontaneous radiative emisison from
    the genralised planks law per wavelength or per photon interval
    Currents it takes in silicons properties by defult_QF_split
    """

    _cal_dts = dict(
        material='Si',
        temp=300.,  # temp in kelvin
        width=0.018,  # width in cm
        ni_author=None,  # author of intrinsic carrier density
        optics_k_author='Green_2008',
        optics_n_author='Green_2008',
        )

    def __init__(self, **kwargs):
        """
        Can provide a specific instance of:
            A material. then it will attempt to look up the optical constants
                intrinsic carrier concentration

            if one of
            optical properties  or  ni module is provided
            These will then be used
        """

        self.calculationdetails = kwargs
        self._update_links()

    def _update_links(self):
        '''
        updates the models used, or takes a passed class
        '''

        self._optics = opticalproperties.TabulatedOpticalProperties(
            material=self._cal_dts['material'],
            temp=self._cal_dts['temp'],
            abs_author=self._cal_dts['optics_k_author'],
            ref_author=self._cal_dts['optics_n_author'])

        self._ni = ni.IntrinsicCarrierDensity(
            material=self._cal_dts['material'],
            author=self._cal_dts['ni_author'],
            temp=self._cal_dts['temp'])

    def blackbody_photon_per_wavelength(self, emn_wavelegnth=None, **kwargs):
        """
        Returns photon emission per wavelength interval per solid angle for a
        black body emission.

        sometimes is is obsrved is with a factor of 4 pi,
        this is when integrated for emission in all direction

        input:
            wavelength: (np array)
                    wavelength in nm
            kwargs: (optional)
                anything returned by the cal_dts function

        Currently:
         1. I divided by 10000, for some reaon? I have temp commented it out,
         2. I have no just multiplied by 1000 for NO reason
        """

        if emn_wavelegnth is None:
            emn_wavelegnth = self._optics.wavelength
        # wavelength to meters
        emn_wavelegnth = emn_wavelegnth * 1e-9

        return 2 * const.c / emn_wavelegnth**4 * 1. / (
            np.exp(const.h * const.c / emn_wavelegnth / const.k /
                   self._cal_dts['temp']) - 1.) * 1000

    def genralised_planks_PerWavelength_Carriers(self, np=1e16, **kwargs):
        """
        generalised planks law.
        Provides emitted photons per wavelength interval
        Uses the format outlined by green

        inputs:
                ni (float, optional):
                    the product of carriers in cm^6
        returns:
                The black body emission spectrum
        """
        if bool(kwargs):
            self.calculationdetails = kwargs
            self._update_links()

        # black body here is per solid angle
        BB = self.blackbody_photon_per_wavelength(temp=self._cal_dts['temp'])

        # The PL spectrum with no QF splitting
        rsp_thermal = (
            BB * self._optics.abs_cof_bb) / self._optics.ref_ind**2

        return rsp_thermal * ((np) / self._ni.update()**2)

    def genralised_planks_PerEnergy(self, QF_split=0.1, **kwargs):
        """
        generalised planks law.
        Provides emitted photons per energy interval
        Uses the traditional form

        inputs:
                QF_split (float, optional):
                    Quasi fermi energly level splitting in eV
        returns:
                The black body emission spectrum
        """

        if bool(kwargs):
            self.calculationdetails = kwargs
            self._update_links()

        QF_split *= const.e

        E = const.h * const.c / (self._optics.wavelength * 1e-9)

        # speed of light in medium
        try:
            c = const.c / self._optics.ref_ind
        except:
            c = const.c

        # Density of state of phtons
        D = E**2 / c**3 * 2**2 * const.pi / const.h**3

        # Note that in Gesikers phd he droped from the denumerator
        # The spectrum with QF splitting
        return self._optics.abs_cof_bb * c * D / (
            np.exp(E / const.k / self._cal_dts['temp']) *
            np.exp(-QF_split / const.k / self._cal_dts['temp']) - 1
        )

    def genralised_planks_PerWavelength(self, **kwargs):
        """
        generalised planks law.
        Provides emitted photons per wavelength interval
        Is just an adjustedment to the energy interval expression
        """

        if bool(kwargs):
            self.calculationdetails = kwargs
            self._update_links()

        # we just need to multip the per energy by the derivative below
        dEdwl = const.h * const.c / (self._optics.wavelength * 1e-9)**2

        # Adjust the values
        return self.genralised_planks_PerEnergy(**kwargs) * dEdwl


class luminescence_emission(HelperFunctions):

    """
    A class that simualted the PL emitted by a device
    i.e. tries to account for reabsorption and reflections
    """

    # Given a deltan V x profile  provides the PL out the system
    # Currently can not adjust for dector

    # Dictionaries
    wafer_optics_dic = {'polished': 'double_side_polished',
                        'textured': 'double_side_lambertian'}
    PL_Dection_side_depth = {'rear': 'escape_rear',
                             'front': 'escape_front'}

    _cal_dts = dict(
        wafer_opitcs='polished',
        detection_side='front',
        material='Si',
        temp=300.,  # temp in kelvin
        width=0.018,  # width in cm
        ni_author=None,  # author of intrinsic carrier density
        optics_k_author='Green_2008',
        optics_n_author='Green_2008',
        nxc=np.ones(10),  # the number of excess carrier with depth
        doping=1e16,  # the doping in cm^-3
        )

    def __init__(self, **kwargs):

        self.calculationdetails = kwargs

        self._index = None

        self._update_x_dist()
        self._update_links()

    def _update_x_dist(self):
        '''
        updates the distance
        '''
        self._x = np.linspace(0, self._cal_dts['width'],
                              self._cal_dts['nxc'].shape[0])

    def _update_links(self):
        '''
        A function where the links to other
        instances's is completely refreshed.
        '''
        self._sre = SpontaneousRadiativeEmission(
            temp=self._cal_dts['temp'],
            optics_k_author=self._cal_dts['optics_k_author'],
            optics_n_author=self._cal_dts['optics_n_author'],
            material=self._cal_dts['material'],
            ni_author=self._cal_dts['ni_author'],
            )

        # I got lasy, so i'm using the previous classes stuff
        self._optics = self._sre._optics

        if self._index is None:
            self._index = self._optics.wavelength > 0
        elif self._index.shape != self._optics.wavelength.shape:
            self._index = self._optics.wavelength > 0

        self._optics.wavelength = self._optics.wavelength[self._index]
        self._optics.abs_cof_bb = self._optics.abs_cof_bb[self._index]
        self._optics.ref_ind = self._optics.ref_ind[self._index]

        self._sre._optics = self._optics

        self._esc = absorptance.EscapeProbability(
            material=self._cal_dts['material'],
            optics_k_author=self._cal_dts['optics_k_author'],
            optics_n_author=self._cal_dts['optics_n_author'],
            temp=self._cal_dts['temp'],
            x=self._x)

        self._esc._optics = self._optics

        self._update_escape()

    def update_carrierdensity(self, deltan, doping=None):
        """
        inputs for carrier density
        If not doping is given, used the doping in self.doping
        """
        if doping is not None:
            self._cal_dts['doping'] = doping

        if deltan.shape != self.x.shape:
            print('number of x-values not equal to delta n values')

        self.np = self._cal_dts['doping'] * deltan

    def limit_wavelegnths(self, wl_min=None, wl_max=None):
        """
        Used for obtained the basic values required for this caculation
        i.e. optical cosntants, ni, an escape fraction
        """

        self._index = self._optics.wavelength > wl_min
        self._index *= self._optics.wavelength < wl_max

    def _update_escape(self):
        """
        Can be used to update the escape fraction, no inputs
        """

        getattr(self._esc, self.wafer_optics_dic[
            self._cal_dts['wafer_opitcs']])()

        self._escapeprob = getattr(
            self._esc,
            self.PL_Dection_side_depth[self._cal_dts['detection_side']])

    def calculate_spectral(self, **kwargs):
        """
        deteries the spectral PL emitted from a sample

        don't think this is correct at the moment
        """
        # ensure inputs are good
        if bool(kwargs):
            self.calculationdetails = kwargs
            self._update_x_dist()
            self._update_links()

        # cacualte the generated PL
        sre = self._sre.genralised_planks_PerWavelength_Carriers(
            self._sre._ni.update()**2)

        # this is the spectral distribution from each point
        # Normalised to deltan = 1, so we can just multi this by deltan
        assert self._cal_dts['nxc'].shape == self._x.shape, (
            "nxc is different length to x spacing")

        Spectral_PL = np.trapz(
            (
                (sre * self._escapeprob).T *
                self._cal_dts['nxc'] *
                self._cal_dts['doping'] /
                self._sre._ni.update()**2),
            self._x,
            axis=1)

        return Spectral_PL

    def calculate_emitted(self, **kwargs):
        """
        multiples the detected PL by an EQE
        currently this does NOTHING
        """
        spectral = self.calculate_spectral(**kwargs)
        return np.trapz(spectral, self._optics.wavelength)
