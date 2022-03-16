import scipy.special as sps
import scipy.optimize as sop
import numpy as np
import matplotlib.pyplot as plt


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
    return np.sqrt(i)*(1-intens/i)

def loss_2(params, *args):
    intens, func = args
    n = np.arange(intens.shape[0])
    i = func(n, *params)
    return intens-i

def loss_3(params, *args):
    intens, func = args
    n = np.arange(intens.shape[0])
    i = func(n, *params)
    return np.sqrt(intens)*(1-i/intens)

def estimate_location(intens, fun, g):

    # first guess of parameters
    n = np.arange(intens.shape[0])
    A = np.sum(intens)
    mu = np.sum( intens * n) / np.sum(intens)
    sigma = np.sqrt(np.sum(intens * (n-mu)**2) / np.sum(intens))
    y_offset = np.min(intens)

    params0 = np.array([A, mu, sigma, y_offset/2])
    bounds = (np.array([0, mu-3, sigma/3, 0]), np.array([1.5*A, mu+3, 3*sigma, y_offset]))


    res = sop.least_squares(fun, params0, method='dogbox', bounds=bounds, ftol=1e-8, args=(intens, g))
    return res.x

nbins = 7

res1=[]
res2=[]

mus = np.linspace(1.8, 3.2, 50)
for mu in mus:
    y_offset = 0.05
    A = 1
    sigma = 0.9

    #generating data
    n = np.arange(nbins)
    ig = igauss(n, A, mu, sigma, y_offset)
    ig += 0.001*np.random.normal(nbins) * np.sqrt(ig)


    res1.append(estimate_location(ig, loss_2, gauss)[1])
    res2.append(estimate_location(ig, loss_3, igauss)[1])

plt.plot(mus, res1-mus, label='g' )
plt.plot(mus, res2-mus, label='ig')
plt.legend()
plt.show()
