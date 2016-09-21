#!/usr/local/bin/python
# UTF-8

import numpy as np
from semiconductor.material.ni import IntrinsicCarrierDensity as NI
# import fdint as fd
import scipy.constants as const


def get_carriers(Na, Nd, nxc,
                 temp=300, material='Si', ni_author=None, ni=None):
    '''
    returns the carrier densities given the number of ionised dopants and ni
    and the excess carriers

    input:
    Na = the number of ionised acceptors
    Nd = the number of ionised donors
    nxc = the excess carrier density. In this function assume
                  deltap = deltan
    temp = temperature
    ni: (optional)
        provide  a values so this function doesn't calculate ni

    returns ne, nh

    '''

    # check types
    if not isinstance(Na, np.ndarray):
        Na = np.asarray([Na])
    if not isinstance(Nd, np.ndarray):
        Nd = np.array([Nd])
    if not isinstance(nxc, np.ndarray):
        nxc = np.array([nxc])

    # print('enter', Na.shape, Na.shape, nxc.shape)
    # if ni not provided obtain
    if ni is None:
        ni = NI(material=material).update(author=ni_author, temp=temp)

    # Calculated on the assumption that at thermal equilibrium in the
    # dark n0p0 = ni**2, and that charge neutrality holds. Usually

    # simplified to saying the majority carrier density ~ the doping and min
    # carrier denisty is the number of excess carriers. The below version
    # more accurately incorporates ni though, which is particularly important
    # for temperature dependent measurements.

    maj_car_den = 0.5 * (np.abs(Nd - Na) +
                         np.sqrt((Nd - Na)**2 + 4 * ni**2))

    if not isinstance(maj_car_den, np.ndarray):
        # print('is not instance', maj_car_den)
        maj_car_den = np.array([maj_car_den])
        # print('now is ', maj_car_den)
    if len(maj_car_den.shape) == 0:
        # print('length 0 ', maj_car_den)
        maj_car_den = np.array([maj_car_den])
        # print('after change is ', maj_car_den)

    min_car_den = ni**2 / maj_car_den
    # assign the dark minority carriers
    ne0 = np.copy(min_car_den)
    nh0 = np.copy(min_car_den)

    # check the doping and assign
    # if the number of donars are larger
    index = Na < Nd

    # print(maj_car_den, index, Na, Nd)
    # print(maj_car_den[~index], maj_car_den[index])
    nh0[~index] = maj_car_den[~index]
    ne0[index] = maj_car_den[index]

    # add the number of excess carriers
    assert ne0.shape == nh0.shape

    # make the excess carrier the right shape to add
    nxc_m = nxc * np.ones(nxc.shape[0])

    # add them to the dark carriers
    ne = ne0 + nxc_m
    nh = nh0 + nxc_m
    # print('the adding shapes', ne0.shape, nxc_m.shape)
    # To make sure an array the same shape as the input array is returned
    # Note this function does not play well if  both an array of dopants
    # and excess carrerier densities are passed
    # print(ne0, type(ne0), maj_car_den, type(maj_car_den))
    # print(ne0.shape)

    # if ne0.shape[0] == 1:
    #     ne = ne[0]
    #     nh = nh[0]
    # else:
    #     ne = ne.flatten()
    #     nh = nh.flatten()

    # print('exit', ne.shape, nh.shape)

    return ne, nh


def fermi2carrier_fermi(Ef, ni_author=None, eg_author=None, temp=300,
                        material='Si', Ei=0):
    '''
    This does not work.

    determines the number of carriers from the fermi energy level

    inputs:
        Ef: (array like)
            The Fermi energy level referenced to the intrinsic level

    '''

    # this may have to be the effective ni, i'm not sure
    ni = NI(material=material).update(author=ni_author, temp=temp)
    ni = 2.831801e10
    print(ni)
    Nc = 2.857082e19
    Nv = 2.513669e19

    Eg = -np.log(ni**2. / Nc / Nv) * const.k * temp / const.e

    Ei = (Eg * const.e / 2. + 0.5 * const.k * temp * np.log(Nv / Nc)) / const.e

    dEc = (Ef - (Eg - Ei))
    dEv = -dEc - Eg

    n0i = fd.parabolic(Ef / const.k / temp * const.e) * ni

    n = fd.parabolic(dEc / const.k / temp * const.e) * Nc
    p = fd.parabolic(dEv / const.k / temp * const.e) * Nv

    print(dEc, dEv, Ei, Eg - Ei, )
    print(n0i, Eg)

    print('{0:.2e} {1:.2e}'.format(n * p, ni**2))

    return n, p


def fermi2carrier_boltz(Ef, ni_author=None, eg_author=None, temp=300, material='Si', Ei=0):
    '''
    This does not work.

    determines the number of carriers from the fermi energy level

    inputs:
        Ef: (array like)
            The Fermi energy level referenced to the intrinsic level

    '''

    # this may have to be the effective ni, i'm not sure
    ni = 2.831801e10
    print(ni)
    Nc = 2.857082e19
    Nv = 2.513669e19
    Eg = -np.log(ni**2 / Nc / Nv) * const.k * temp / const.e

    Ei = (Eg * const.e / 2 + 0.5 * const.k * temp * np.log(Nv / Nc)) / const.e

    dEc = (Ef - (Eg - Ei))
    dEv = -dEc - Eg

    n = np.exp(dEc / const.k / temp * const.e) * Nc
    p = np.exp(dEv / const.k / temp * const.e) * Nv

    return n, p


def carrier2fermi_fermi(ne, nh, ni_author=None, eg_author=None, temp=300, material='Si', Ei=0):
    '''
    This does not work.

    determines the number of carriers from the fermi energy level

    inputs:
        Ef: (array like)
            The Fermi energy level referenced to the intrinsic level

    '''

    # this may have to be the effective ni, i'm not sure
    ni = 2.831801e10
    Nc = 2.857082e19
    Nv = 2.513669e19

    # determine the band gap from the input paramters
    Eg = -np.log(ni**2. / Nc / Nv) * const.k * temp / const.e

    # Valance band to Ei
    Ei = (Eg * const.e / 2. + 0.5 * const.k * temp * np.log(Nv / Nc)) / const.e

    dEc = (Ef - (Eg - Ei))
    dEv = -dEc - Eg

    Efe = fd.iparabolic(ne / Nc) - (Eg - Ei)
    Efh = fd.iparabolic(nh / Nv) - Ei

    return Efe, Efh

# fermi2carrier(0.3)
