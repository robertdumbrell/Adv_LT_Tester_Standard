
import numpy as np


def Roosbroeck(vals, nxc, nh0, ne0, B, **kwargs):
    '''
    The classic roosbroeck function
        It simply states that the radiaative recombiation rate is
        proportional to the produce of the carrier, correct to the
        background conccentraion of carriers.
    '''
    if np.all(nxc == 0):
        nxc += 1

    nh = nh0 + nxc
    ne = ne0 + nxc

    R = B * (ne * nh - ne0 * nh0)
    # print R
    return nxc / R


def Roosbroeck_with_screening_B(vals, nxc, doping, temp, Blow):
    """
    This is the roosbroeck model that accounts for many things,
    such as band gap narrowing.
    It needs temperature, nxc, doping and blow to be defined
    """

    bmin = vals['rmax'] + (vals['rmin'] - vals['rmax']) / (
        1. + (temp / vals['r1'])**vals['r2'])
    b1 = (vals['smax'] + (vals['smin'] - vals['smax']) / (
        1. + (temp / vals['s1'])**vals['s2'])) * 2
    b3 = (vals['wmax'] + (vals['wmin'] - vals['wmax']
                          ) / (
        1. + (temp / vals['w1'])**vals['w2'])) * 2

    # print bmin

    B = Blow * (bmin + (vals['bmax'] - bmin) / (
        1. + ((2. * nxc + doping) / b1
              )**vals['b2']
        + ((2. * nxc + doping) / b3)**vals['b4']))

    return B


def Roosbroeck_with_screening(vals, nxc, nh0, ne0, Blow, temp):
    """
    This is the roosbroeck model that accounts for many things
    It needs temperature, nxc, doping and blow to be defined
    """
    B = Roosbroeck_with_screening_B(vals, nxc, np.amax([nh0, ne0]), temp, Blow)
    tau = Roosbroeck(vals, nxc, nh0, ne0, B=B)
    return tau


def cubic_loglog_parm(vals, temp):
    '''
    a cubic fit to tabulated data at different temperatures
    '''

    B = vals['a'] * np.log10(temp)**3 +\
        vals['b'] * np.log10(temp)**2 \
        + vals['c'] * np.log10(temp)**1 \
        + vals['d']
    return 10**B
