
import numpy as np
import matplotlib.pylab as plt
import scipy.constants as const


def EnergyToWavelength(data):
    '''
    appends array with the wavelength data
    data is required to have the headings:
    'energy' and 'alpha'
    '''

    # created the names and type for the named array
    names = list(data.dtype.names + ('wavelength',))
    fmat = ['f4' for i in names]

    datanew = np.zeros((data.shape[0]),
                       dtype={'names': names, 'formats': fmat}
                       )

    # assignes the values
    for name in data.dtype.names:
        datanew[name] = data[name]

     # cacls and appends the wavelength values
    datanew['wavelength'] = const.c * const.h / data['energy'] / const.e * 1e9

    # returns the new array
    return datanew


def wl2nrg(data):
    '''
    appends array with the wavelength data
    data is required to have the headings:
    'energy' and 'alpha'
    '''

    # created the names and type for the named array
    names = list(data.dtype.names + ('energy',))
    fmat = ['f8' for i in names]

    datanew = np.zeros((data.shape[0]),
                       dtype={'names': names, 'formats': fmat}
                       )

    # assignes the values
    for name in data.dtype.names:
        datanew[name] = data[name]

     # cacls and appends the wavelength values
    datanew['energy'] = (const.c * const.h) / data['wavelength'] * 1e9 / const.e

    # returns the new array
    return datanew


def rmfield(a, *fieldnames_to_remove):
    return a[[name for name in a.dtype.names if name not in fieldnames_to_remove]]


def add_wavelength_oren_ergy_to_file(fname):

    data = np.genfromtxt(
        fname, names=True, delimiter=',', filling_values=np.nan)

    data = rmfield(data, 'energy')

    if 'wavelength' in data.dtype.names and 'energy' not in data.dtype.names:
        data = wl2nrg(data)
    elif 'wavelength' not in data.dtype.names and 'energy' in data.dtype.names:
        data = EnergyToWavelength(data)

    # a = tuple(b.split(';'))
    # data.dtype.names = a

    np.savetxt(fname, data, delimiter=',', header=','.join(
        data.dtype.names), comments='', fmt=('%0.4e'), )

# fnames = ['Si_Schinke2014','Si_Green2008','Si_Green1995','Si_Daub1995']
fname = 'Silicon_Nguyen2014test'
add_wavelength_oren_ergy_to_file(fname)

plt.semilogy()
plt.show()
