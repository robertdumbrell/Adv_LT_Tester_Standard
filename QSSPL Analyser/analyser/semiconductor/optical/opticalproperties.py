import numpy as np
import configparser
import os
import sys
import scipy.constants as const
from semiconductor.helper.helper import HelperFunctions


class TabulatedOpticalProperties(HelperFunctions):
    '''
    A class to obtain n, k, and alpha
    '''
    _cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'abs_author': None,
        'ref_author': None,
        'ext_cof': False
    }

    def __init__(self, **kwargs):

        self.calculationdetails = kwargs
        self._update_links()
        self.load()

    def _update_links(self):
        self.tac = TabulatedAbsorptionCoefficient(
            material=self._cal_dts['material'],
            author=self._cal_dts['abs_author'],
            temp=self._cal_dts['temp'],
        )

        try:
            self.tri = TabulatedRefractiveIndex(
                material=self._cal_dts['material'],
                author=self._cal_dts['ref_author'],
                temp=self._cal_dts['temp'],
            )
        except:
            print(sys.exc_info()[0])
            sys.exit('No reactive index data found')

    def load(self, common_range=True, **kwargs):
        self.calculationdetails = kwargs

        if 'author' or 'material' in ''.join(kwargs.keys()):
            self._update_links()

        try:
            self.tri.load(temp=self._cal_dts['temp'])
            self.ref_ind = self.tri.ref_ind
        except:
            pass

        self.tac.load(temp=self._cal_dts['temp'])
        self.abs_cof_bb = self.tac.abs_cof_bb

        if self._cal_dts['ext_cof']:
            self.ext_cof_bb = self.tac.calculate_ext_coef()

        if common_range:
            # set the values
            self.wavelength = self.tac.wavelength
            self.energy = self.tac.energy

            # update the refractive index
            try:
                self.ref_ind = self.tri.ref_ind_at_wls(self.wavelength)
            except:
                pass

            # remove unused variables
            self.wl_abs_cof_bb = None
            self.wl_ref_ind = None
        else:
            # just use as is
            self.wl_abs_cof_bb = self.tac.wl
            self.wl_ref_ind = self.tri.wl

            # remove unused variables
            self.wavelength = None

        pass


class TabulatedAbsorptionCoefficient(HelperFunctions):

    """
    A class containg the optical constants of silicon
    These are temperature dependence.
    """

    _cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'author': None,
    }

    author_file = r'tabulated_absorption_coefficient.const'

    def __init__(self, **kwargs):

        # update any values in cal_dts
        # that are passed
        self.calculationdetails = kwargs

        # get the address of the authors list
        author_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            self._cal_dts['material'],
            self.author_file)

        # get the models ready
        self._int_model(author_file)

        # initiate the first model
        self.change_model(self._cal_dts['author'])

    def load(self, **kwargs):
        """
        Loads alpha and n
        from the provided or from self. the name
        """
        # check to see if its set here out outside of this function
        self.calculationdetails = kwargs

        # if change the model, update it.
        if 'author' in kwargs.keys():
            self.change_model(self._cal_dts['author'])

        # Getting the absorption coefficient from a file
        data = np.genfromtxt(os.path.join(os.path.dirname(__file__),
                                          self._cal_dts['material'],
                                          self.model),
                             names=True, delimiter=',')

        # need something here to get temp dependence
        self.wavelength, self.energy = data[
            'wavelength'], data['energy']

        # if the temp is a float
        if type(self.vals['temp']) is float:
            self.abs_cof_bb = data['alpha']
            if self._cal_dts['temp'] != self.vals['temp']:
                try:
                    self.abs_cof_bb = _temp_power_law(self.abs_cof_bb,
                                                      data['C_ka'],
                                                      self._cal_dts['temp'],
                                                      self.vals['temp'])
                except:
                    print(
                        '''Warning:'''
                        '''\n\tNo tabulated data, or temp cofs for'''
                        ''' {0:.0f} K'''.format(self._cal_dts['temp']) +
                        '''\tfor the author {0}'''.format(
                            self._cal_dts['author']) +
                        '\tusing data for temperature {0:.0f} K.'.format(
                            self.vals['temp'])
                    )

        else:

            # this happens when there are several alpha values, so lets try a
            # specif temp
            name = 'alpha_{0:.0f}K'.format(self._cal_dts['temp'])
            if name in data.dtype.names:
                self.abs_cof_bb = data[name]
            else:
                # if doesn't work just use the stipulated default
                print(
                    'Warning:'
                    '\n\tTabulated data at {0} K does not exist.'.format(
                        self._cal_dts['temp']) +
                    '\n\tfor the author {0}'.format(
                        self._cal_dts['author']) +
                    '\n\tThe value for', self.vals['default_temp'],
                    'K is used'
                )
                name = 'alpha_{0:.0f}K'.format(self.vals['default_temp'])
                self.abs_cof_bb = data[name]

        try:
            # try to get the uncertainty
            self.U = data['U']

        except:
            pass

    def alphaBB_at_wls(self, wavelength):
        return np.interp(wavelength,
                         self.wavelength,
                         self.abs_cof_bb)

    def calculate_ext_coef(self):
        self.ext_cof_bb = self.abs_cof_bb * self.wavelength / 4 / np.pi
        return self.ext_cof_bb


class TabulatedRefractiveIndex(HelperFunctions):

    _cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'author': None,
    }

    author_file = r'tabulated_refractive_index.const'

    def __init__(self, **kwargs):

        # update any values in cal_dts
        # that are passed
        self.calculationdetails = kwargs

        # get the address of the authors list
        author_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            self._cal_dts['material'],
            self.author_file)

        # get the models ready
        self._int_model(author_file)

        # initiate the first model
        self.change_model(self._cal_dts['author'])

    def load(self, **kwargs):
        """
        Loads n
        from the provided or from self. the name
        """

        self.calculationdetails = kwargs

        # a check to make sure the model hasn't changed
        if 'author' in kwargs.keys():
            self.change_model(self._cal_dts['author'])

        # To do
        # need to make a check if there is a temp value
        # if there is use it, if not check if there are temp coefficients
        # there there are use them
        # if nothing return the default temp value.

        # Get n
        data = np.genfromtxt(os.path.join(os.path.dirname(__file__),
                                          self._cal_dts['material'],
                                          self.model),
                             names=True, delimiter=',')

        self.wavelength, self.ref_ind, self.energy = data[
            'wavelength'], data['n'], data['energy']

        if self._cal_dts['temp'] != self.vals['temp']:
            try:
                self.ref_ind = _temp_power_law(self.ref_ind, data['C_n'],
                                               self._cal_dts['temp'],
                                               self.vals['temp'])
            except:
                print(
                    'Temp Warning:'
                    '\tNo tabulated data, or temp cofs for {0:.0f} K'.format(
                        self._cal_dts['temp']) +
                    '\tfor the author {0}'.format(self._cal_dts['author']) +
                    '\tusing data for temperature {0:.0f} K.'.format(
                        self.vals['temp'])
                )

    def ref_ind_at_wls(self, wavelength):
        '''
        returns the refrative index n's
        and the supplied wavelengths

        inputs:
            wavelength (array)
        output:
            n (array)
        '''
        return np.interp(wavelength,
                         self.wavelength,
                         self.ref_ind)


def _temp_power_law(ref_vairable, coef, temp, ref_temp):
    '''
    The coefficient is assumed to be times by 1e-4
    '''
    return ref_vairable * np.power(temp / ref_temp, coef*1e-4*ref_temp)


class ModelledAbsorptionCoefficient(HelperFunctions):

    ''' This purpose of this it to provide access if the absorption
        coefficient have a model.
        This is a work in progress.
    '''
    beta = 0
    gamma = 0

    _cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'author': None,
    }

    author_file = 'modelled_absorption_coefficient.const'

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

    def update_absorptioncoefficients(self, f=None, Input=None):
        '''
        updates the absorption coefficients
        inputs:
            f: can input a frequency range that you want alpha returned for
                if not supplied will use self.f
            Input: (optional) set to wavelength to input wavelength in nm

        returns:
            absorption coefficients in cm^-1
        '''
        if f is None:
            f = self.f
        else:
            if Input == 'wavelength':
                f = self._wavelength2frequency(wavelength)

        self.alpha = getattr(self, self.model)(self.vals, f)

        return self.alpha

    def _EgwithT(self, Eg, T, gamma=None, beta=None):
        if beta is None:
            beta = self.beta

        if gamma is None:
            gamma = self.gamma
        return Eg - beta * T**2 / (T + gamma)

    def _checkf(self, f):
        if f is None:
            f = self.f
        return f

    def _wavelength2frequency(self, Lambda):
        '''
        Takes lambda in nm
        provides f in Hz
        '''
        self.f = const.c / Lambda / 1e-9
        return self.f

    def _alpha_function(self, f, Eth, A, power, T=300.):
        '''
        Generic function to determine alpha
        inputs:
            f is the photon energy in hz
            Eth is the threshold energy in eV
            A is the probability coefficient
            power is the power its applied to

        '''

        # Use provided f or self.f
        f = self._checkf(f)

        # change E with T, assumes that the
        # band narrowing is not a function of Eg
        Eth = self._EgwithT(Eth, T)

        # calculate alpha
        alpha = A * (const.h * f / const.e - Eth)**power

        # make sure we don't have any impossible values
        alpha[const.h * f / const.e - Eth < 0] = 0

        return alpha

    def alpha_p_absorption(self, Eg, Ep, A, T, f=None):
        '''Phonon absorption, with E as power 2'''
        A /= (np.exp(Ep / const.k / T * const.e) - 1)
        return self._alpha_function(f, Eg - Ep, A, 2, T)

    def alpha_p_emission(self, Eg, Ep, A, T, f=None):
        '''Phonon emission, with E as power 2'''
        A /= (1 - np.exp(-Ep / const.k / T * const.e))
        return self._alpha_function(f, Eg + Ep, A, 2, T)

    def alpha_exciton_emission(self, Eg, Ep, Ee, A, T, f=None):
        '''Phonon and exciton emission, E power of 0.5'''
        A /= (1 - np.exp(-Ep / const.k / T * const.e))
        return self._alpha_function(f, Eg + Ep + Ee, A, 0.5, T)

    def alpha_exciton_absorption(self, Eg, Ep, Ee, A, T, f=None):
        '''Phonon and absorption emission, E power of 0.5
        This process may not be possible is high quality semiconductors'''
        A /= (1 - np.exp(-Ep / const.k / T * const.e))
        return self._alpha_function(f, Eg + Ep - Ee, A, 0.5, T)

    def alpha_direct(self, Eg, A, T=300, f=None, power=0.5):
        '''Direct band gap'''
        return self._alpha_function(f, Eg, A, power, T)

    def alpha_indirect(self, Eg, Ephonon, A, T):
        """
        The change in phonon absorption probability between emission
        and absorption is just e^{E_p / kt}, taken from MacFarlane
        """

        alpha_in = np.zeros(self.f.shape)
        for Ep, Aa in zip(Ephonon, A):
            Ae = Aa * np.exp(Ep * const.e / const.k / T)
            alpha_in += self.alpha_p_absorption(Eg, Ep, Aa, T)
            alpha_in += self.alpha_p_emission(Eg, Ep, Ae, T)
        return alpha_in

    def MacFarlane(self, vals, f):

        # Eg = vals['egi']

        Egi_list = [vals[s] for s in vals if 'egi' in s]
        Ep_list = [vals[s] / 1000 for s in vals if 'ep' in s]
        Ap_list = [vals[s.replace('e', 'a')] for s in vals if 'ep' in s]

        Egd_list = [vals[s] for s in vals if 'egd' in s]
        Ad_list = [vals[s.replace('e', 'a')] for s in vals if 'ad' in s]

        # print Ep_list, Ap_list
        alpha = np.zeros(self.f.shape)

        for Eg in Egi_list:
            # print Eg
            alpha += self.alpha_indirect(Eg, Ep_list, Ap_list, 300.)

        for Eg in Egd_list:
            alpha += self.alpha_direct(Eg, Ad_list)

        return alpha

    def Rajkanan(self, vals, f):
        # Based on indirect theory from Elliot

        # print vals
        self.gamma = vals['gamma']
        self.beta = vals['beta']

        Egi_list = [vals[s] for s in vals if 'egi' in s]
        Ep_list = [vals[s] / 1000 for s in vals if 'ep' in s]
        Ap_list = [vals[s.replace('e', 'a')] for s in vals if 'ap' in s]
        C_list = [vals[s.replace('e', 'a')] for s in vals if 'c' in s]

        Egd_list = [vals[s] for s in vals if 'egd' in s]
        Ad_list = [vals[s.replace('e', 'a')] for s in vals if 'ad' in s]

        alpha = np.zeros(self.f.shape)
        T = 300.
        for Eg, Ap in zip(Egi_list, Ap_list):
            for Ep, C in zip(Ep_list, C_list):

                alpha += self.alpha_p_absorption(Eg, Ep, Ap * C, T)
                alpha += self.alpha_p_emission(Eg, Ep, Ap * C, T)

        for Eg in Egd_list:
            alpha += self.alpha_direct(Eg, Ad_list)

        return alpha

    def Bucher(self, vals, f):
        # Based on indirect theory from Elliot

        # print vals
        self.gamma = vals['gamma']
        self.beta = vals['beta']

        Egi_list = [vals[s] for s in vals if 'egi' in s]
        Ep_list = [vals[s] / 1000 for s in vals if 'ep' in s]
        Ap_list = [vals[s.replace('e', 'a')] for s in vals if 'ap' in s]
        C_list = [vals[s.replace('e', 'a')] for s in vals if 'c' in s]

        Egd_list = [vals[s] for s in vals if 'egd' in s]
        Ad_list = [vals[s.replace('e', 'a')] for s in vals if 'ad' in s]

        alpha = np.zeros(self.f.shape)
        T = 300.
        for Eg, Ap in zip(Egi_list, Ap_list):
            for Ep, C in zip(Ep_list, C_list):
                alpha += self.alpha_p_absorption(Eg, Ep, Ap * C, T)
                alpha += self.alpha_p_emission(Eg, Ep, Ap * C, T)

        power = vals['power']

        for Eg in Egd_list:
            alpha += self.alpha_direct(Eg,
                                       Ad_list[0] / const.h / f * const.e, 1.5)

        return alpha
