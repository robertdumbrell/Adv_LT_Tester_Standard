
import numpy as np
import scipy.constants as const


def Green_1990(vals, temp, Egratio, **kargs):
    """
     This form as described by Green in 10.1063/1.345414.
     inputs:
        vals: (dic)
            the effect mass values
        temp: (float)
            the temperature in kelvin

    outputs:
        vel_th_c: (float)
            the termal velocity for the conduction in cm/s
        vel_th_v: (float)
            the termal velocity for the valance band in cm/s
    """

    # the values relative to the rest mass
    ml = vals['ml'] * const.m_e
    mt = vals['mt'] * Egratio * const.m_e

    delta = np.sqrt((ml-mt)/ml)
    # conduction band effective mass
    mth_c = 4.*ml/(
        1.+np.sqrt(ml/mt) * np.arcsin(delta)/delta)**2

    vel_th_c = np.sqrt(8*const.k * temp / np.pi / mth_c)
    # valance band effective mass, its a 7 order poynomial fit
    mth_v = np.sum(
        [vals['meth_v'+str(i)]*temp**i for i in range(8)]) * const.m_e

    vel_th_v = np.sqrt(8*const.k * temp / np.pi / mth_v)

    # adjust the values from m/s to cm/s and return
    return vel_th_c * 100, vel_th_v * 100


def constants(vals, **kwargs):
    """
    Returns a constant value
    """

    return vals['vth_e'], vals['vth_h']
