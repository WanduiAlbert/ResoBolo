#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt


resistance = np.array([3.77, 3.77, 3.48, 3.92, 3.51, 3.45, 3.46, 3.59, 3.48 ])
index = np.arange(len(resistance)) + 1

fig, ax = plt.subplots(figsize=(10,10))
ax.scatter(index + 0.5, resistance)
ax.set_xticks(index + 0.5)
ax.set_xticklabels(['B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'])
ax.set_xlabel(r'Chip')
ax.set_ylabel(r'Resistance [mOhms]')
ax.set_ylim(ymin=3.4, ymax=4.0)
ax.grid(which='both')
ax.vlines(5, 3.4, 4.0, colors='r', linestyles='dashed')

plt.savefig('Dip_probe_results.png')
