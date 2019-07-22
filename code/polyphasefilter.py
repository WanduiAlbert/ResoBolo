#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

P = 4
N = 1 << 10
M = N*P
t = np.arange(M)
x = np.random.randn(M)

filt = np.sinc((t - t[M//2]))

simp_fft = np.fft.fft(x)
