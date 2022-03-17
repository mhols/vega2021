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
    return (i-intens)/np.sqrt(i)

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

def play_1():
    nbins = 9

    res1=[]
    res2=[]
    res3=[]

    mus = np.linspace(1.8, 3.2, 50)
    for mu in mus:
        A = 1000
        y_offset = 0.01*A
        sigma = 0.7

        #generating data
        n = np.arange(nbins)
        ig = igauss(n, A, mu, sigma, y_offset)
        ig += np.random.normal(nbins) * np.sqrt(ig)

        res1.append(estimate_location(ig, loss_2, gauss)[1])
        res2.append(estimate_location(ig, loss_1, igauss)[1])
        res3.append(estimate_location(ig, loss_3, igauss)[1])
    plt.plot(mus, res1-mus, label='g' )
    plt.plot(mus, res2-mus, label='ig')
    plt.plot(mus, res3-mus, label='iig')
    plt.legend()
    plt.show()


def play_2(lo, ga):
    """
    lo: the loss function be used
    ga: the gaussian model to be used
    """
    nbins = 9  # number of bins
    ntrial = 10000  # number of random samples

    # prepare a random collection of parameters
    A = 10000
    mus = np.random.uniform(low=nbins/2-2, high=nbins/2+2, size=ntrial)
    sigmas = np.random.uniform(low=0.5, high=4, size=ntrial)
    y_offsets= np.random.uniform(low=0.001, high=0.1, size=ntrial)*A

    res = []  # to collect the simulation results

    for mu, sigma, y_offset in zip(mus, sigmas, y_offsets):

        #generating data
        n = np.arange(nbins)
        ig = igauss(n, A, mu, sigma, y_offset)
        ig += np.random.normal(size=nbins) * np.sqrt(ig)   #  poor man's Poisson error

        res.append(estimate_location(ig, lo, ga))
    res = np.array(res)
    Ae, mue, sigmae, y_offsete = res.T

    plt.figure()
    plt.title('reduced misfit as function of estimated sigma')
    plt.plot(sigmae, (mue-mus)/sigmae**1.5, 'o')

    plt.show()

play_2(loss_2, gauss)
