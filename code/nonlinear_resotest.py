
import numpy as np
from scipy import optimize
from math import pi
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pl


MHz = 1e6

def real_of_complex(z):
    ''' flatten n-dim complex vector to 2n-dim real vector for fitting '''
    r = np.hstack((z.real,z.imag))
    return r
def complex_of_real(r):
    assert len(r.shape) == 1
    nt = r.size
    assert nt % 2 == 0
    no = nt//2
    z = r[:no] + 1j*r[no:]
    return z

def model_linear(f,f0,A,phi,D,dQr,dQe_re,dQe_im,a):
    f0 = f0 * 1e6
    cable_phase = np.exp(2.j*pi*(1e-6*D*(f-f0)+phi))
    dQe = dQe_re + 1.j*dQe_im
    x = (f - f0)/f0
    s21 = A*cable_phase*(1. - dQe/(dQr + 2.j*x))
    return s21

def get_roots(y0, a):
    coeffs = [4, -4*y0, 1, -(y0+a)]
    roots = np.roots(coeffs)
    return roots[np.isreal(roots)]

def get_y(y0, a, f):
    low_to_high = np.all(np.diff(f) > 0)
    if low_to_high:
        find_roots = lambda x: np.min(get_roots(x, a))
    else:
        find_roots = lambda x: np.max(get_roots(x, a))
    ys = list(map(find_roots, y0))
    return np.array(ys)

def make_model(dQe_re,dQe_im):
    def f(f,f0,A,phi,D,dQr,a):
        return full_model(f,f0,A,phi,D,dQr,dQe_re,dQe_im,a)
    return f

def model(f,f0,A,phi,D,dQr,dQe_re,dQe_im,a):
    f0 = f0 * 1e6
    cable_phase = np.exp(2.j*pi*(1e-6*D*(f-f0)+phi))
    dQe = dQe_re + 1.j*dQe_im
    x0 = (f - f0)/f0
    y0 = x0/dQr
    k2 = np.sqrt((y0**3/27. + y0/12. + a/8.)**2 - (y0**2/9. - 1/12.)**3, dtype=np.complex128)
    k1 = np.power(a/8. + y0/12. + k2 + y0**3/27., 1./3)
    eps = (-1. + 3**0.5 * 1j)/2.

    #np.seterr(all='raise')
    #print (np.sum(np.abs(k1) == 0.0))
    y1 = y0/3. + (y0**2/9.-1/12.)/k1 + k1
    y2 = y0/3. + (y0**2/9.-1/12.)/eps/k1 + eps*k1
    y3 = y0/3. + (y0**2/9.-1/12.)/eps**2/k1 + eps**2*k1

    y1[np.abs(k1) == 0.0] = y0[np.abs(k1) == 0.0]/3.
    y2[np.abs(k1) == 0.0] = y0[np.abs(k1) == 0.0]/3.
    y3[np.abs(k1) == 0.0] = y0[np.abs(k1) == 0.0]/3.

    # Out of the three roots we need to pick the right branch of the bifurcation
    thresh = 1e-4
    low_to_high = np.all(np.diff(f) > 0)
    if low_to_high:
        y = y2.real
        mask = (np.abs(y2.imag) >= thresh)
        y[mask] = y1.real[mask]
    else:
        y = y1.real
        mask = (np.abs(y1.imag) >= thresh)
        y[mask] = y2.real[mask]

    #if np.sum(k1 == 0.0) > 0 :
    #    y[k1 == 0.0] = y0[k1 == 0.0]/3
    #y = get_y(y0, a, f)
    x = y*dQr
    s21 = A*cable_phase*(1. - (dQe)/(dQr + 2.j*x))
    return s21

def test_fit():
    nt = 1000
    f0 = 300
    A = 1.0
    phi = 0.0
    D = 0.
    dQr = 1./16000
    dQe_re = 1./20000
    dQe_im =  -0.2*dQe_re
    dQe = dQe_re + 1j * dQe_im
    Qe = 1./dQe
    Qi = 1/(dQr - dQe_re)
    all_as = [0.0, 0.5, 0.9, 1.5]
    for a in all_as:
        p0 = (f0,A,phi,D,dQr,dQe_re,dQe_im,a)
        freq = np.linspace(f0-0.1, f0+0.1, nt)*MHz
        freq = freq[::-1]
        z = model(freq, *p0)
        z_lin = model_linear(freq, *p0)
        fmax = f0*MHz*(1-a*dQr)
        print (fmax)
        pt = model(np.array([fmax]), *p0)
        pt_lin = model_linear(np.array([f0*MHz]), *p0)

        pl.figure(figsize=(10,10))
        pl.plot(freq/MHz - f0, np.abs(z), 'k-', label='Non Linear Model a = %1.3f'%a)
        pl.plot(freq/MHz - f0, np.abs(z_lin), 'k--', label='Linear Model')
        pl.vlines(fmax/MHz - f0, 0, 1, colors='r', linestyles='dotted',
                label='Non Linear Model: Max Responsivity')
        pl.vlines(f0 - f0, 0, 1, colors='b', linestyles='dotted',
                label='Linear Model: Max Responsivity')
        pl.xlabel('Relative Frequency (MHz)')
        pl.ylabel('|S21|')
        pl.grid()
        pl.legend(loc='lower right')
        #pl.xlim(-1,1)
        #pl.ylim(-1,1)
        pl.title('Qi=%d Qr=%d Qe_re=%d Qe_im=%d'%(Qi,1./dQr,1./dQe_re, Qe.imag))
        pl.savefig('S21_mag_plot_a%1.2f.png'%a)
        pl.close()

        pl.figure(figsize=(10,10))
        pl.plot(z.real, z.imag, 'k-', label='Non Linear Model a = %1.3f'%a)
        pl.plot(z_lin.real, z_lin.imag, 'k--', label='Linear Model')
        pl.plot(pt.real, pt.imag, 'rd', label='Non Linear Model: Max Responsivity')
        pl.plot(pt_lin.real, pt_lin.imag, 'bs',
                label='Linear Model: Max Responsivity')
        pl.xlabel('I')
        pl.ylabel('Q')
        pl.grid()
        pl.legend(loc='best')
        #pl.xlim(-1,1)
        #pl.ylim(-1,1)
        pl.axis('square')
        pl.title('Qi=%d Qr=%d Qe_re=%d Qe_im=%d'%(Qi,1./dQr,1./dQe_re, Qe.imag))
        pl.savefig('S21_IQ_plot_a%1.2f.png'%a)
        pl.close()

        pl.figure(figsize=(10,10))
        pl.plot(freq/MHz - f0, z.real, 'k-', label='Non Linear Model a = %1.3f'%a)
        pl.plot(freq/MHz - f0, z_lin.real, 'k--', label='Linear Model')
        pl.vlines(fmax/MHz - f0, 0, 1, colors='r', linestyles='dotted',
                label='Non Linear Model: Max Responsivity')
        pl.vlines(f0 - f0, 0, 1, colors='b', linestyles='dotted',
                label='Linear Model: Max Responsivity')
        pl.xlabel('Relative Frequency (MHz)')
        pl.ylabel('Re S21')
        pl.grid()
        pl.legend(loc='lower right')
        #pl.xlim(-1,1)
        #pl.ylim(-1,1)
        pl.title('Qi=%d Qr=%d Qe_re=%d Qe_im=%d'%(Qi,1./dQr,1./dQe_re, Qe.imag))
        pl.savefig('S21_Re_plot_a%1.2f.png'%a)
        pl.close()

        pl.figure(figsize=(10,10))
        pl.plot(freq/MHz - f0, z.imag, 'k-', label='Non Linear Model a = %1.3f'%a)
        pl.plot(freq/MHz - f0, z_lin.imag, 'k--', label='Linear Model')
        pl.vlines(fmax/MHz - f0, -0.6, 0.6, colors='r', linestyles='dotted',
                label='Non Linear Model: Max Responsivity')
        pl.vlines(f0 - f0, -0.6, 0.6, colors='b', linestyles='dotted',
                label='Linear Model: Max Responsivity')
        pl.xlabel('Relative Frequency (MHz)')
        pl.ylabel('Im S21')
        pl.grid()
        pl.legend(loc='lower right')
        #pl.xlim(-1,1)
        #pl.ylim(-1,1)
        pl.title('Qi=%d Qr=%d Qe_re=%d Qe_im=%d'%(Qi,1./dQr,1./dQe_re, Qe.imag))
        pl.savefig('S21_Im_plot_a%1.2f.png'%a)
        pl.close()

if __name__=='__main__':
    test_fit()
