#! /usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, FuncFormatter

# I'm defining a function that would give me the correct formatting for the
# y-axis labels.
def powersof10(number, *pos):
  exponent = int(np.trunc(np.log10(number)))
  retstr = r"$\texttt{10}^" +  r"{0:1d}$".format(exponent)
  return retstr


filename = "publication.txt"

fs, asd, nep = np.loadtxt(filename, unpack=True, skiprows=1)

photon_neps = [92.10, 74.20, 42.3]
photon_labels = [r'$\nu$ = 270 GHz', r'$\nu$ = 230 GHz', r'$\nu$ = 150 GHz']
photon_fmt = ['r--', 'k--', 'g--']

yticks = np.concatenate((np.arange(10, 100, 10), np.arange(100, 500, 100)))
yticklabels = [powersof10(x) for x in [10, 100]]
fig, ax = plt.subplots(figsize=(10, 10))
ax.loglog(fs, nep)
#ax.yaxis.set_major_formatter(StrMethodFormatter("{x:1.0e}"))
ax.set_xlabel(r'Frequency [Hz]')
ax.set_ylabel(r'NEP [aW/$\sqrt{\texttt{Hz}}$]')
ax.grid(which='both')

ones = np.ones_like(fs)
for i in range(3):
  ax.loglog(fs, ones*photon_neps[i], photon_fmt[i], label=photon_labels[i])
ax.legend(loc='best', fontsize=30, title="Photon NEP")
#ax.set_yticks([10, 100])
#ax.set_yticklabels(yticklabels)
ax.axis([0.1, 500, 10,500])
ax.xaxis.set_major_formatter(StrMethodFormatter("{x:1.1f}"))
ax.yaxis.set_major_formatter(FuncFormatter(powersof10))

plt.savefig('bestNEP.png')
plt.savefig('bestNEP.pdf')
plt.show()
