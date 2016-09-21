
import numpy as np
import sys
import os
import scipy.constants as const
import matplotlib.pylab as plt

from semiconductor.helper.helper import HelperFunctions
from semiconductor.general_functions.carrierfunctions import get_carriers
from semiconductor.material.ni import IntrinsicCarrierDensity as ni
from semiconductor.material.thermal_velocity import ThermalVelocity as Vel_th

sys.path.append(
    os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)))


class SRH(HelperFunctions):

    # ToDo:
    # need to assign which thermal velocity values
    # were used for each value of capture cross section
    # then need to look this value up from the model.
    # currently just have this as an input.

    _cal_dts = {
        'material': 'Si',
        'defect': None,
        'temp': 300.,
        'vth_author': None,
        'Nt': 1e10,  # the number of traps
        'Nd': 0,
        'Na': 1e16,
        'nxc': 1e10,
        'ni_author': None
    }

    vel_th_e = None
    vel_th_h = None

    author_list = 'SRH.defects'

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

        # initiate the a defect
        self._change_model(self._cal_dts['defect'])
        self._update_links()
        self._cal_taun_taup()

    def _update_links(self):

        self.ni = ni(material=self._cal_dts['material'],
                     author=self._cal_dts['ni_author'],
                     temp=self._cal_dts['temp'],
                     )
        self.ni.update()

        self.vel_th_e, self.vel_th_h = Vel_th(
            material=self._cal_dts['material'],
            author=self._cal_dts['vth_author'],
            temp=self._cal_dts['temp']).update()

    def _cal_taun_taup(self):
        '''
        Determines the SRH lifetime values
        from the capture cross sections
        '''

        self.vals['tau_e'] = 1. / self._cal_dts['Nt']\
            / self.vals['sigma_e'] / self.vel_th_e  # s
        self.vals['tau_h'] = 1. / self._cal_dts['Nt']\
            / self.vals['sigma_h'] / self.vel_th_h  # s

    def _change_model(self, defect):
        '''
        A stage to update the thermal velocity model,
        if one way provided with the SRH defect model
        '''

        # load the defect
        self.change_model(self._cal_dts['defect'])
        print(self.vals.keys(), 'insdie')
        # check is a model was provided
        if 'vth_author' in self.vals.keys():
            self._cal_dts['vth_author'] = self.vals['vth_author']

        # get the values from the model
        self._update_links()

    def tau(self, **kwargs):
        '''
        reports the lifetime of the current defect for the given
        excess carrier density.
        '''

        if bool(kwargs):
            self.calculationdetails = kwargs

        if 'defect' in kwargs:
            self._change_model(self._cal_dts['defect'])
            self._cal_taun_taup()

        # if change in a model update the values
        if 'author' in ''.join(kwargs.keys()):
            self._update_links()
            self._cal_taun_taup()

        return self._tau(self._cal_dts['nxc'],
                         self.vals['tau_e'],
                         self.vals['tau_h'],
                         self.vals['et'])

    def itau(self, **kwargs):
        return 1. / self.tau(**kwargs)

    def _tau(self, nxc, tau_e, tau_h, Et):
        """
        The Shockley Read Hall lifetime for a defect level
        in the band gap. This is calculated from the kinetics
        of repairing.
        """

        # the escape from defects
        nh1 = self.ni.ni * \
            np.exp(-Et * const.e / (const.k * self._cal_dts['temp']))

        ne1 = self.ni.ni * \
            np.exp(Et * const.e / (const.k * self._cal_dts['temp']))


        # get the number of carriers
        ne, nh = get_carriers(Na=self._cal_dts['Na'],
                              Nd=self._cal_dts['Nd'],
                              nxc=self._cal_dts['nxc'],
                              temp=self._cal_dts['temp'],
                              material=self._cal_dts['material'],
                              ni=self.ni.ni
                              )

        # calculate the recombination rate
        U = (ne * nh - self.ni.ni**2) / \
            (tau_h * (ne + ne1) + tau_e * (nh + nh1))

        return self._cal_dts['nxc'] / U

    def usr_vals(self, Et=None, sigma_e=None, sigma_h=None,
                 tau_e=None, tau_h=None, Nt=None):
        '''
        a function to provide arbitory values for SRH lifetime.
        inputs:
            Et: (optional float)
                The energy level from the intrinsic level in eV
            sigma_e: (optional float)
                The capture cross section for electrons in seconds
            sigma_h: (optional float)
                The capture cross section for holes in seconds
        '''
        temp = self.vals
        # assign the new values
        self.vals = {
            'et': Et or temp['et'],
            'sigma_e': sigma_e or temp['sigma_e'],
            'sigma_h': sigma_h or temp['sigma_h'],
            'doi':  None,
            'notes': 'User defined value',
            'tau_e': np.array([tau_e]) or temp['tau_e'],
            'tau_h': np.array([tau_h]) or temp['tau_h'],
            }

        # the or above doens't work for zeros
        if Et == 0:
            self.vals['et']  = 0

        if 'vth_author' in self.vals:
            del self.vals['vth_author']

        self._cal_dts['Nt'] = Nt or self._cal_dts['Nt']

        # if tau_n or tau_p's values passed,
        # do dont' caculate them
        if tau_e is None and tau_h is None:
            self._cal_taun_taup()
        else:
            self.vals.update(
                {'sigma_e': None,
                 'sigma_h': None,
                 }
            )

    def _plot_all(self):
        fig, ax = plt.subplots(1, 2, figsize=(16, 6))
        # ax = plt.add_subplot(111)
        counter = 0
        defects = self.AvailableModels()
        deltan = np.logspace(12, 17)
        for defect in defects:
            # ax.plot(np.inf,np.inf,'k-',label = 'Auger')
            self.select_defect(defect)
            self.cal_tau(1e10)
            ax[0].plot(counter, self.Et, 'o', label=defect)
            ax[1].plot(deltan, self.tau(deltan) * 1e6, label=defect)
            counter += 1
        # ax.legend(loc=0)
        # ax.loglog()
        x = np.arange(counter)
        # print defects
        ax[0].set_xticks(x)
        ax[0].set_xticklabels(defects, rotation=90)
        ax[0].set_xlabel('Defect')
        ax[0].set_ylabel('E$_t$ from Ei (eV)')
        ax[1].loglog()
        ax[1].set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax[1].set_ylabel('Lifetime (us)')
