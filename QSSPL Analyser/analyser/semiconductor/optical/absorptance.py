import numpy as np
import matplotlib.pylab as plt
import sys
import os
import scipy.constants as Const
import semiconductor.optical.opticalproperties as opticalproperties
from semiconductor.helper.helper import HelperFunctions


class EscapeProbability(HelperFunctions):
    '''
    A class that contains models for the escape probabilities
    '''

    _cal_dts = dict(
        material='Si',
        temp=300.,  # temp in kelvin
        wafer_opitcs='polished',
        detection_side='front',
        width=0.018,  # width in cm
        ref_front=0.1,
        ref_rear=0.1,
        optics_k_author=None,
        optics_n_author=None,
        )

    def __init__(self,
                 x=None, **kwargs):

        self.calculationdetails = kwargs

        if x is None:
            self.x = np.linspace(0, 0.018, 100)
        else:
            self.x = x

    def _update_links(self):

        self._optics = opticalproperties.TabulatedOpticalProperties(
            material=self._cal_dts['material'],
            temp=self._cal_dts['temp'],
            abs_author=self._cal_dts['optics_k_author'],
            ref_author=self._cal_dts['optics_n_author']
            )

    def width_from_xlegnth(self):

        self._cal_dts['width'] = self.x[-1]

    def double_side_lambertian(self, ax=None, **kwargs):
        """
        This is Rudigers model (2007), though its taken from Schinkes
        2013 paper as he wrote it with some mistakes and in a hard to
        read manner. Its derived from 2 Lambertian surfaces and infinite
        reflections

        This depends on the samples
        Width, refractive index. That seems to be it.
        There is no input for surface reflections
        """
        if kwargs:
            self.calculationdetails = kwargs
            self._update_links()

        self.width_from_xlegnth()

        # self.theta = 60.
        # self.theta_rad = self.theta / 180. * Const.pi
        # self.theta = 0

        # This line is as it is requried in the equation used below
        a = 1. / self._optics.n**2
        b = 1. / \
            (1. - (1 - a)**2 *
             np.exp(-4 * self._optics.abs_cof_bb * self._cal_dts['width']))
        c = (1. - a) * np.exp(-4 * self._cal_dts['width'] *
                              self._optics.abs_cof_bb)

        # Now cac both ways
        xd_W = 2 * self._optics.abs_cof_bb * \
            ((np.ones([self.x.shape[0], self._optics.abs_cof_bb.shape[0]]).T) *
             (self._cal_dts['width'] - self.x)).T

        xd = 2 * self._optics.abs_cof_bb * \
            (np.ones([self.x.shape[0], self._optics.abs_cof_bb.shape[0]]).T *
             self.x).T

        # Calc both methods
        self.escape_front = a * b * (
            np.exp(-xd) + c * np.exp(xd)
        )

        # print np.exp(-xd_W).shape, (c * np.exp(xd_W)).shape

        self.escape_rear = a * b * (
            np.exp(-xd_W) + c * np.exp(xd_W)
        )

        if ax is not None:
            ax.plot(self.x,  (np.exp(-xd)), 'r')
            ax.plot(self.x,  (c * np.exp(xd)), 'b')
            ax.plot(self.x,  self.Escape_front, 'g')

            # ax.plot(self.wavelength_emission, self._optics.abs_cof_bb)
            # ax.plot(self.wavelength_emission, self._optics.n)

    def double_side_polished(self, **kwargs):
        """
        This is taken from Schick1992 paper
        the reflectin values should be provided as a decimal not a fraction
        """

        if kwargs:
            self.calculationdetails = kwargs
            self._update_links()

        self.width_from_xlegnth()
        self.theta = 0  # This is for a polished sample

        xd_rear = self._optics.abs_cof_bb * \
            ((np.ones([self.x.shape[0], self._optics.abs_cof_bb.shape[0]]).T) *
             (self._cal_dts['width'] - self.x)).T

        b_rear = self._cal_dts['ref_front'] * \
            np.exp(-2 * self._optics.abs_cof_bb * self._cal_dts['width'])

        xd_frot = self._optics.abs_cof_bb * \
            (np.ones([self.x.shape[0], self._optics.abs_cof_bb.shape[0]]).T *
             self.x).T
        b_frot = self._cal_dts['ref_rear'] * \
            np.exp(-2 * self._optics.abs_cof_bb * self._cal_dts['width'])

        self.escape_rear = (
            np.exp(-xd_rear) + b_rear * np.exp(xd_rear)
        ) / (1 - b_rear * self._cal_dts['ref_rear'])

        self.escape_front = (
            np.exp(-xd_frot) + b_frot * np.exp(xd_frot)
        ) / (1 - b_frot * self._cal_dts['ref_front'])

        # self._optics.abs_cof_bb*=cos(self.theta*np.pi/180)

    def general_form(self):
        """
        Taken from
        C. Schinke, D. Hinken, J. Schmidt, K. Bothe, and R. Brendel, IEEE
         J. Photovoltaics 3, 1038 (2013).
         in the from
        f(x) =  scale * [A * exp(-x) + B exp(W-x) + C (exp(-x)+D)]

        this paramterisation was also used by
        R. Brendel, M. Hirsch, R. Plieninger, and J.H. Werner, IEEE Trans.
         Electron Devices 43, 1104 (1996).

        Please remember this is a 9 variable fit. it should always give a
         good fit
        """

        # calculated from law of diffraction considering angle of entry
        # theta = gamma  - arcsin(n_air/n_si sin(gamma))
        theta1 = np.pi / 3 * 0
        # the angle light is reflected from the rear into the device
        theta2 = np.pi / 3 * 0
        # the angle light is then reflected from that point
        thetan = np.pi / 3 * 0

        # these are defined as exp{-alpha W / cos(theta_x))}
        T1 = np.exp(-self._optics.abs_cof_bb * self._cal_dts['width'] /
                    np.cos(theta1))
        # print T1
        # experimentally this is taken as the average of T1 and Tn.
        # T2 =  (Lambda * Rbd Tn + (1- Lambda) Rbs T1)/(Lambda * Rbn + (1+
        # Lambda)* Rbs)
        # where:
        # Rbd is the rear surface reflection for lambertian reflected light a
        # Rbs is the rear surface reflectance for specular reflected light
        T2 = np.exp(-self._optics.abs_cof_bb * self._cal_dts['width'] /
                    np.cos(theta2))
        Tn = np.exp(-self._optics.abs_cof_bb * self._cal_dts['width'] /
                    np.cos(thetan))

        Rb1 = 0.0

        # taken from the reflectance measured at short wavelengths and
        # extrapolated, in this way it is only the reflected light
        Rf = 0.

        # this can be related to other things
        # Rf1 = Lambda * Rbd * Tn * Rfn + (1 - Lambda) *Rbs T1 * Rfs /
        # (Lambda * Rbd * Tn + (1-Lambda)+ R bs + T1 )
        # Where Rfs is taken from ray tracing
        Rf1 = 0.
        Rfn = 0.
        #
        Rbn = 0.

        # Lambda and Rb are taken from a fit to the reflectance curve
        # R = 1 - (1-Rf) [1 - T1*Rb1 * T2 * (1-Rf1)- T1*Rb1*T2*Rf1*Rbn *
        # (1-R) * Tn^2]\
        #       1 - Rbn*Rfn * Tn^2

        scale = (1. - Rf)

        A = 1. / np.cos(theta1)
        B = T1 * Rb1 / np.cos(theta2)
        C = T1 * Rb1 * T2 * Rf1 / (1. - Tn * Tn * Rfn * Rbn) / np.cos(thetan)
        D = Tn * Rbn

        # print A, B, C, D

        xd_W = self._optics.abs_cof_bb * \
            ((np.ones([self.x.shape[0], self._optics.abs_cof_bb.shape[0]]).T) *
             (self._cal_dts['width'] - self.x)).T

        xd = self._optics.abs_cof_bb * \
            (np.ones([self.x.shape[0], self._optics.abs_cof_bb.shape[0]]).T *
             self.x).T

        self.Escape_front = scale * (
            A * np.exp(-xd / np.cos(theta1)) +
            B * np.exp(-xd_W / np.cos(theta2)) +
            C * (np.exp(-xd / np.cos(thetan)) +
                 D * np.exp(-xd_W / np.cos(thetan))
                 )
        )


if __name__ == "__main__":
    a = EscapeProbability()

    a.double_side_polished(0, 0)
    plt.plot(a.x,
             a.Escape_front[:, 1] / np.amax(a.Escape_front[:, 1]),
             label='polished')

    a.double_side_lambertian()
    plt.plot(a.x,
             a.Escape_front[:, 1] / np.amax(a.Escape_front[:, 1]), '--',
             label='lambert')

    a.general_form()
    plt.plot(a.x,
             a.Escape_front[:, 1] / np.amax(a.Escape_front[:, 1]), ':',
             label='gen')
    # print a.Escape_front[:, 1]

    plt.legend(loc=0)

    plt.show()
