import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import astropy.units as u

def file_loader(filename, numports=1):
    data = np.genfromtxt(filename,  delimiter=',', skip_header=8)

    print (data.shape)
    frequency = data[:,0] * 1000 * u.MHz
    if numports == 1 :
        Y11 = (data[:,1] + 1j * data[:,2]) / u.Ohm
        return frequency, Y11
    elif numports == 2 :
        Y11 = (data[:,1] + 1j * data[:,2]) / u.Ohm
        Y12 = (data[:,3] + 1j * data[:,4]) / u.Ohm
        Y21 = (data[:,5] + 1j * data[:,6]) / u.Ohm
        Y22 = (data[:,7] + 1j * data[:,8]) / u.Ohm

        return frequency, Y11, Y12, Y21, Y22
    else:
        raise ValueError("Can only use 1 or 2 ports in this module")

# def file_loader(filename, numports=1, delimiter=' '):
#     data = np.genfromtxt(filename,  delimiter=delimiter, skip_header=8)

#     print (data.shape)
#     frequency = data[:,0] * 1000 * u.MHz
#     if numports == 1 :
#         Z11 = (data[:,1] + 1j * data[:,2]) * u.Ohm
#         return frequency, Z11
#     elif numports == 2 :
#         Z11 = (data[:,1] + 1j * data[:,2]) * u.Ohm
#         Z12 = (data[:,3] + 1j * data[:,4]) * u.Ohm
#         Z21 = (data[:,5] + 1j * data[:,6]) * u.Ohm
#         Z22 = (data[:,7] + 1j * data[:,8]) * u.Ohm

#         return frequency, Z11, Z12, Z22
#     else:
#         raise ValueError("Can only use 1 or 2 ports in this module")


def ZtoY(Z11, Z12, Z22):
    detZ = Z11 * Z22 - Z12**2
    Y11 = Z22/detZ
    Y12 = -Z12/detZ
    Y22 = Z11/detZ
    return Y11, Y12, Y22


if __name__=="__main__":
    file_1a = "parasitic_capacitance_1a.csv"
    file_1b = "parasitic_capacitance_1b.csv"
    file_2a = "parasitic_capacitance_2a.csv"
    file_2b = "parasitic_capacitance_2b.csv"

    frequency, Y11_1a = file_loader(file_1a, 1)
    _, Y11_1b = file_loader(file_1b, 1)
    _, Y11_2a, Y12_2a, Y21_2a, Y22_2a = file_loader(file_2a, 2)
    _, Y11_2b, Y12_2b, Y21_2b, Y22_2b = file_loader(file_2b, 2)

    print ("All the files have been loaded")
    # Get the corresponding Y parameters

    # Let's make comparisons first between 1a and 1b.
    C_1a = (np.imag(Y11_1a)/(2*np.pi*frequency)).to('pF')
    C_1b = (np.imag(Y11_1b)/(2*np.pi*frequency)).to('pF')
    C11_2a= (np.imag(Y11_2a)/(2*np.pi*frequency)).to('pF')
    C12_2a= (np.imag(-Y12_2a)/(2*np.pi*frequency)).to('pF')
    C11_2b= (np.imag(Y11_2b)/(2*np.pi*frequency)).to('pF')
    C12_2b= (np.imag(-Y12_2b)/(2*np.pi*frequency)).to('pF')

    C12_est = 1/(1/C_1a + 1/C12_2b)
    # fig, ax = plt.subplots(figsize=(10,10))
    # ax.plot(frequency, C_1a, 'b', label="With Ground Plane")
    # ax.plot(frequency, C_1b, 'r', label="Distant Ground Plane")
    # ax.set_xlabel(r'Frequency [MHz]')
    # ax.set_ylabel(r'C [pF]')
    # ax.grid(which='both')
    # ax.legend(loc='best')
    # ax.yaxis.set_major_formatter(StrMethodFormatter("{x:1.3f}"))
    # ax.axis('tight')

    # Now let's see how well the pi model matches with the data
    reciprocality_2a = Y12_2a - Y21_2a
    reciprocality_2b = Y12_2b - Y21_2b
    # assert(np.all(Y12_2b[1:] == Y21_2b[1:]))
    deltaY_a =  np.imag(Y11_1a) - np.imag(Y11_2a + Y12_2a) 
    deltaC_a = (deltaY_a/(2*np.pi*frequency)).to('pF')

    # fig, ax = plt.subplots(figsize=(10,10))
    # ax.plot(frequency, np.imag(Y11_1a), 'b', label=r"$Y_{11}$ for 1a")
    # ax.plot(frequency, np.imag(Y11_2a + Y12_2a), 'r', label=r"$Y_{11}$ for 2a")
    # ax.set_xlabel(r'Frequency [MHz]')
    # ax.set_ylabel(r'Admittance to Ground Y $[\Omega^{-1}]$')
    # ax.grid(which='both')
    # ax.yaxis.set_major_formatter(StrMethodFormatter("{x:1.3e}"))
    # ax.legend(loc='best')
    # ax.axis('tight')

    # Now let's see how well the pi model matches with the data
    # fig, ax = plt.subplots(figsize=(10,10))
    # ax.plot(frequency, deltaC_a, 'b')
    # ax.set_xlabel(r'Frequency [MHz]')
    # ax.set_ylabel(r'Difference in Capacitance to Ground $\Delta C$ [pF]')
    # ax.grid(which='both')
    # ax.yaxis.set_major_formatter(StrMethodFormatter("{x:1.3f}"))
    # ax.axis('tight')

    # Compare Y_12 for 2a vs 2b using the capacitance
    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot(frequency, np.real(reciprocality_2a), 'b', label=r"$Y_{12}$ ")
    ax.set_xlabel(r'Frequency [MHz]')
    ax.set_ylabel(r'Y $[\Omega^{-1}]$')
    ax.grid(which='both')
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:1.3e}"))
    ax.axis('tight')


    # fig, ax = plt.subplots(figsize=(10,10))
    # ax.plot(frequency, C_12, 'b', label=r"$C_{11}$")
    # ax.set_xlabel(r'Frequency [MHz]')
    # ax.set_ylabel(r'$C_{12}$ [pF]')
    # ax.grid(which='both')
    # ax.set_title("Capacitance to Top Side Ground Plane")
    # ax.yaxis.set_major_formatter(StrMethodFormatter("{x:1.3f}"))
    # ax.axis('tight')

    plt.show()