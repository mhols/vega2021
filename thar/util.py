import scipy.special as sps
import scipy.optimize as sop
import numpy as np
import sys

used_params = {}    # for the documentation of used parameters
                    # pass **kwargs to every method
                    # and use get_kwarg(kwargs, 'key', default)


def get_kwarg(kwargs, key, default):
    res  = kwargs.get(key, default)
    used_params[key] = res
    return res

def gauss(x, A, mu, sigma, y_offset):
    return y_offset + (A/  (np.sqrt(2*np.pi)*sigma)) * np.exp(-(x-mu)**2/(2*sigma**2))

def igauss(x, A, mu, sigma, y_offset):
    return y_offset +  0.5 * A * (
        sps.erf((x + 0.5 - mu)/(np.sqrt(2)*sigma))
        - sps.erf((x -0.5 -mu)/(np.sqrt(2)*sigma)))

def cauchy(x, A, mu, sigma, y_offset):
    return y_offset + A*sigma/(sigma**2 + (x-mu)**2)

def loss_1(params, *args):
    intens, func = args
    n = np.arange(intens.shape[0])
    i = func(n, *params)
    return np.abs(i-intens)/np.sqrt(np.abs(i)+1e-8)

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


def mean_histogram_estimator(intens):
    """
    removing baseline followed by histogram mean
    """
    lag = 2 # for continuum estimation
    cont = (np.mean(intens[:lag])+np.mean(intens[-lag:]))/2
    cont = 0 #0.5*(intens[0]+intens[-1])
    intens = intens - max(0.0, np.min(intens))
    n = np.arange(intens.shape[0])
    intens = np.abs(intens)
    return np.sum(n*intens) / np.sum(intens)

def estimate_location(intens, fun, g):

    # backproject on positivity
    intens = intens - max(0, np.min(intens))
    # first guess of parameters
    n = np.arange(intens.shape[0])
    A = np.sum(np.abs(intens))
    Amin = max(0, A-4*np.sqrt(A))
    Amax = A+4*np.sqrt(A)

    mu = len(intens)/2
    mumin = 0
    mumax = len(intens)+1

    sigma = len(intens)
    sigmamax = 4*sigma
    sigmamin = 0.1

    y_offset_min = 0
    y_offset_max = 10 + max(0, np.min(intens))
    y_offset = (y_offset_min+y_offset_max)/2

    params0 = np.array([A, mu, sigma, y_offset])

    bounds = (np.array([Amin, mumin, sigmamin, y_offset_min]),
              np.array([Amax, mumax, sigmamax, y_offset_max]))

    res = sop.least_squares(
        fun, 
        params0, 
        bounds=bounds, 
        #method='dogbox', 
        args=(intens, g))
    return res.x

function_map = {
    'gauss': gauss,
    'igauss': igauss,
    'cauchy': cauchy,
    'loss_1': loss_1,
    'loss_2': loss_2,
    'loss_3': loss_3,
}


def bootstrap_estimate_location(intens, **kwargs):
    g = function_map[kwargs['profile']]
    size = kwargs['n_bootstrap']
    fun = function_map[kwargs['loss_function']]
    intens = intens - min(0, np.min(intens))
    p = intens / np.sum(intens)

    n = int(np.sum(intens))
    intens = np.random.multinomial(n, p, size)
    res = [ estimate_location(inte, fun, g)[1] for inte in intens]

    return np.mean(res), np.std(res)

def sigma_clipping(
    data, 
    epsilon=1, 
    fit=lambda d: 1,
    **kwargs):
    """
    data: list of data records
    fit: maps data into iterable of floats measuring the deviation
    epsilon: cutoff for filtering
    n: number of iterations (<=0 or None means iterate until no change)
    """
    i = 0
    lenold = 0
    n = kwargs.get('n_sigma_clip', 5)
    while i < n:
        I = fit(data)<=epsilon
        print(i, len(data), len(I))
        data = data[I]
        if len(I)==lenold:
            break
        i+=1
        lenold = len(I)
    return data
