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
    A = np.sum(np.abs(intens))
    mu = np.sum(np.abs(intens) * n) / np.sum(np.abs(intens))
    sigma = np.sqrt(np.sum(np.abs(intens) * (n-mu)**2)) / np.sum(np.abs(intens))
    y_min = np.min(intens)
    y_max = np.min(intens)

    params0 = np.array([A, mu, sigma, 0])
    bounds = (np.array([max(0, A-5*np.sqrt(A)), 0, sigma/10, -y_min-1]), 
              np.array([A+4*np.sqrt(A),         len(intens)-1, 4*sigma,  y_max]))
    res = sop.least_squares(
        fun, 
        params0, 
        bounds=bounds, 
        method='dogbox', 
        args=(intens, g))
    return res.x

def bootstrap_estimate_location(intens, fun, g, size=50):
    intens = intens - min(0, np.min(intens))
    p = intens / np.sum(intens)

    n = int(np.sum(intens))
    intens = np.random.multinomial(n, p, size)
    import matplotlib.pyplot as plt
    res = [ estimate_location(inte, fun, g)[1] for inte in intens]
    """plt.figure()
    for inte in intens:
        plt.plot(inte)
    plt.show()

    plt.plot(res)
    plt.show()
    """
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
