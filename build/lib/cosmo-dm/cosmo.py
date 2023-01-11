"""
Filter functions, growth factor
"""
import scipy
from scipy import misc
from scipy import integrate
import numpy as np



#define window function top hat filter (in Fourier space)
def W_tophat(x):
    return 3*(np.sin(x)-x*np.cos(x))/x**3

def derivative_W_tophat(x):
    return scipy.misc.derivative(W_tophat,x)

#sharp-k filter
def W_sharpk(x):
    return np.heaviside(1-x,0)
def derivative_W_sharpk(x):
    return scipy.misc.derivative(W_sharp,x)


#smooth-k filter
def W_smooth_k(x,param):
    return ((1+x**param.PS.Beta)**-1)
def derivative_W_smooth(x,param):
    Beta = param.PS.Beta
    return -(Beta*x**(Beta-1))/((1+x**Beta)**2)
def W_times_W_prime_smooth(x,param):
    Beta = param.PS.Beta
    return Beta*x**(Beta-1)/((1+x**Beta)**3)


#define the crossing distribution function
def crossing_f_ST(sigm,param):
    """""
    First crossing distribution
    """""
    A = param.PS.A
    q = param.PS.q
    delta_c = param.PS.delta_c
    p = param.PS.p
    return A*np.sqrt(2*q*(delta_c**2/sigm**2)/np.pi)*(1+(sigm**2/(q*delta_c**2))**p)*(np.exp(-q*delta_c**2/(2*sigm**2)))


#define Hubble factor H=H0*E
def E(x,param):
    return np.sqrt(param.cosmo.Om * x**-3 + (1-param.cosmo.Om-param.cosmo.Ol)* x**-2 + param.cosmo.Ol)
def H(z,param):
    """""
    Hubble factor [yr-1] 
    """""
    return param.cosmo.h * 10**(-10) * np.sqrt(param.cosmo.Om**(1+z)**3+(1-param.cosmo.Om*-param.cosmo.Oml)*(1+z)**2+param.cosmo.Ol)

#define D(a) non-normalized
def D_non_normalized(a,param):
    w=integrate.quad(lambda u: 1/(u*E(u,param))**3, 0, a)[0]
    return (5*param.cosmo.Om * E(a,param)/2)*w

#define D normalized
def D(a,param):
    return D_non_normalized(a,param)/D_non_normalized(1,param)


