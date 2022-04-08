import scipy.special as sps
import scipy.optimize as sop
import numpy as np

def gauss(x, A, mu, sigma, y_offset):
    return y_offset + (A/  (np.sqrt(2*np.pi)*sigma)) * np.exp(-(x-mu)**2/(2*sigma**2))

def igauss(x, A, mu, sigma, y_offset):
    return y_offset +  0.5 * A * (
        sps.erf((x+0.5-mu)/(np.sqrt(2)*sigma))
        - sps.erf((x-0.5-mu)/(np.sqrt(2)*sigma)))

def loss_1(params, *args):
    intens, func = args
    n = np.arange(intens.shape[0])
    i = func(n, *params)
    return (i-intens)/np.sqrt(np.abs(i)+1e-8)

def loss_2(params, *args):
    intens, func = args
    n = np.arange(intens.shape[0])
    i = func(n, *params)
    return intens-i

def loss_3(params, *args):
    intens, func = args
    n = np.arange(intens.shape[0])
    i = func(n, *params)
    return (i-intens)/np.sqrt(intens)

def estimate_location(intens, fun, g):

    # check positivity
    intens = intens - min(0, np.min(intens))
    # first guess of parameters
    n = np.arange(intens.shape[0])
    A = np.sum(intens)
    mu = np.sum(intens * n) / np.sum(intens)
    sigma = np.sqrt(np.sum(intens * (n-mu)**2) / np.sum(intens))
    y_min = np.min(intens)
    y_max = np.min(intens)

    params0 = np.array([A, mu, sigma, 0])
    bounds = (np.array([max(0, A-4*np.sqrt(A)), 0, sigma/4, -y_min-1]), 
              np.array([A+4*np.sqrt(A),         len(intens)-1, 4*sigma,  y_max]))
    res = sop.least_squares(
        fun, 
        params0, 
        # bounds=bounds, 
        method='dogbox', 
        args=(intens, g))
    return res.x

