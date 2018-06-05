import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.constants import k, pi

import astropy.units as u


# Have the power spectra in Vrms^2/Hz. Do the final conversion
def vrmstodbm(v):
    return 10*np.log10(v**2 * 20)
def dbmtovrms(P):
    return np.sqrt((10**(P/10))/20)


dBm = u.dB(u.mW)
rate = 1e6
Tmax = 50
N = int(Tmax*rate)
ADC_max = (1<<13)-1
resolution = 1<<14
txgain = 0
ampl = 0.5
dt = Tmax/rate
B = rate/2


# Calculate the noise power expected in each dof of the oscillator
T = 8 #Noise Temperature
noisepower = ((0.5*k*T*u.W).to(dBm)).value
print ("The expected noise power in each dof is %4.3f dBm/Hz"%noisepower)
resistance = 50
sigma = dbmtovrms(noisepower)*B**0.5 # Variance of the white noise spectrum in Vrms
print ("The std of the I/Q fluctuation is %3.2f uVrms"%(sigma*1e6))
R = 5.775e-6 #Vrms/ADC
Vmax = ADC_max * R
Vmin = (ADC_max - resolution)*R

t = np.arange(N)*dt
fr = 320e6
# BWfilter = 1/(1 + 2j*B*t)
BWfilter = np.ones_like(t)

noise_real = np.random.randn(N)
noise_real *= sigma #Rescale the rv by the variance

noise_imag = np.random.randn(N)
noise_imag *= sigma #Rescale the rv by the variance

noise_z = noise_real + 1j*noise_imag

V = (ampl*resolution)*np.ones_like(t, dtype=np.complex128)*R
V *= V
V += noise_z #This should give me a noise time stream like the output of the ADC
V *= BWfilter # Attenuate all fluctuations faster than the bandwidth of the system


# Now I want the power spectra in dBm/Hz
NPERS_DIV = 512
nperseg = N//NPERS_DIV
z_mean = np.mean(V)
z = V * (np.abs(z_mean)/z_mean)
f, S_real = sig.welch(z.real, nperseg=nperseg, fs=rate, detrend='linear', scaling='density')
f, S_imag = sig.welch(z.imag, nperseg=nperseg, fs=rate, detrend='linear', scaling='density')

Sr_dBm = vrmstodbm(S_real**0.5)
Si_dBm = vrmstodbm(S_imag**0.5)

#fig, ax = plt.subplots(figsize=(10,10))
#ax.semilogx(f, Sr_dBm, 'b')
#ax.semilogx(f, Si_dBm, 'r')
#ax.set_xlabel('Frequency [Hz]')
#ax.set_ylabel('Noise Spectrum [dBm/Hz]')
#ax.axis('tight');
#ax.grid(which='both')
#plt.show()

# Let's try and get the phase noise spectrum
phase = np.unwrap(np.angle(V))
f_phase, S_phase = sig.welch(phase, nperseg=nperseg, fs=rate, detrend='linear', scaling='density')
P_carrier = vrmstodbm(np.abs(z_mean))
print ("The carrier power is %3.1f"%P_carrier)

#Now we make a comparison between the spectrum in dBc/Hz and that in rad^2/Hz
# converted to log units
S_phasedB = 10*np.log10(S_phase)

fig, ax = plt.subplots(figsize=(10,10))
ax.semilogx(f_phase, S_phasedB, 'b')
ax.semilogx(f, Si_dBm - P_carrier, 'r')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('Phase Noise Spectrum [dBc/Hz]')
ax.axis('tight');
ax.grid(which='both')
plt.show()
