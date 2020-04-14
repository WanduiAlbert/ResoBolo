
import numpy as np
import matplotlib.pyplot as pl

def S21(f,f0,Qr,Qcmag,Qcphi,A,phi):
    Qc = Qcmag * np.exp(1.0j*Qcphi)
    x = (f-f0)/f0
    z = A * np.exp(1.0j*phi) * (1. - (Qr/Qc) / (1 + 2j*Qr*x))
    return z

### Generate some fake data
def make_fake_data():
    f0 = 300.
    Qr = 10000.
    Qcmag = 13000.
    Qcphi = 0.5
    A = 1.7
    phi = 0.7

    ptrue = (f0,Qr,Qcmag,Qcphi,A,phi)

    df = f0/Qr
    nf = 1000
    f = np.linspace(f0-4*df,f0+6*df,nf)

    zmeas = S21(f,*ptrue)
    namp = 0.02
    zmeas += namp*(np.random.normal(size=nf) + 1j*np.random.normal(size=nf))
    return f,zmeas

### Make a crude approximation to the resonance for initial guesses

def guess_params(f,zmeas):

    pm = 0.5*(zmeas[0] + zmeas[-1])

    # Strip off an estimate of the scale factor to simplify the math
    # use zcorr1 for the rest of the function
    Ag = np.abs(pm)
    phig = np.angle(pm)
    zcorr1 = zmeas / (Ag * np.exp(1.0j*phig))

    # "Kasa" method for circle fitting
    # Find the center of the circle and its diameter
    # This is used to estimate Qr/Qc
    # Circle center is x0,y0
    # radius r0
    x = zcorr1.real
    y = zcorr1.imag
    A = np.vstack([x,y,np.ones(x.size)]).T
    print A.shape
    b = x*x + y*y
    P,resid,rank,s = np.linalg.lstsq(A,b,rcond=None)
    x0 = P[0]/2.
    y0 = P[1]/2.
    r0 = np.sqrt(0.25*(P[0]*P[0] + P[1]*P[1])+P[2])
    z0 = x0 + 1.0j*y0

    # Estimate Qcphig from the angle between the center of the circle and the s21 point far from resonance, which is 1 + 0j after stripping the scale factor
    Qcphig = -np.angle(1. - z0)

    # Estimate the width of the resonance by taking moments of the network analysis
    # This is how we calculate the bandwidth of an antenna
    df = f[1]-f[0]
    z = zcorr1 - np.mean(zcorr1)
    m0 = np.sum(df*np.abs(z))
    m1 = np.sum(f*df*np.abs(z))
    m2 = np.sum(f*f*df*np.abs(z))
    width = np.sqrt((m2/m0) - (m1/m0)**2)

    # Estimate resonant frequency by looking for the diametrically opposed point from S21 far from resonance relative to the center of the circle
    zfg = 2*z0 - 1
    ifg = np.argmin(np.abs(zfg - zcorr1))
    fg = f[ifg]

    '''
    pl.scatter(zcorr1.real,zcorr1.imag)
    ig = np.argmin(np.abs(f-fg))
    #pl.scatter([zfg.real],[zfg.imag],color='k')
    pl.scatter([zcorr1[ifg].real],[zcorr1[ifg].imag],color='k')
    pl.show()
    exit()
    '''

    # The moments calculate sigma, but we want fwhm
    # This is the relationship for a gaussian, which seems to work for the lorentzian resonance too
    fwhm = width / np.sqrt(8*np.log(2))
    Qrg = fg / fwhm

    # diameter of the circle is Qr/Qc
    QrQc = 2*r0
    Qcg = Qrg / QrQc

    print "fg: ",fg
    print "Qrg: ",Qrg
    print "Qcg: ",Qcg
    print "Qcphig: ",Qcphig

    pguess = (fg,Qrg,Qcg,Qcphig,Ag,phig)
    return pguess


f,zmeas = make_fake_data()
pguess = guess_params(f,zmeas)
zguess = S21(f,*pguess)

pl.scatter([0],[0],color='k')
pl.plot(zmeas.real,zmeas.imag)
pl.plot(zguess.real,zguess.imag)
pl.xlabel('Re[S21]')
pl.ylabel('Im[S21]')
pl.title('Initial parameter guess for resonance - no fitting')
pl.gca().set_aspect('equal')
pl.grid()
pl.savefig('demo.png')
pl.show()


