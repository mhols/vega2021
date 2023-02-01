import numpy as np
import matplotlib.pyplot as plt
from util import *

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
    nbins =  11 # number of bins
    ntrial = 1000  # number of random samples

    # prepare a random collection of parameters
    As = np.random.uniform(1000, 500000, size=ntrial)
    mus = np.random.uniform(low=nbins/2-1, high=nbins/2+1, size=ntrial)
    sigmas = np.random.uniform(low=1, high=3, size=ntrial)
    y_offsets= np.random.uniform(low=0.001, high=0.1, size=ntrial)*As

    res = []  # to collect the gaussian fit estimates
    ras = []  # to collect the mean as simple estimator

    for A, mu, sigma, y_offset in zip(As, mus, sigmas, y_offsets):

        #generating data
        n = np.arange(nbins)
        ig = igauss(n, A, mu, sigma, y_offset)
        ig += np.random.normal(size=nbins) * np.sqrt(ig)   #  poor man's Poisson error
        ig = np.abs(ig)  # bad backfolding...
        res.append(estimate_location(ig, lo, ga))
        ras.append(np.sum(ig*(n+0.5))/np.sum(ig))
    res = np.array(res)  # transforming into numpy array for better indexing
    ras = np.array(ras)
    Ae, mue, sigmae, y_offset = res.T  # the columns are the estimates of the params


    # this random variable seems to have a distribution which does
    # more or less not depend on the parameters
    reduced_error = Ae**(1/3) * (mue-mus)/sigmae**2

    print (np.mean(reduced_error), np.std(reduced_error))

    plt.figure()
    plt.title('reduced misfit as function of sigma')
    plt.plot(sigmas, reduced_error, 'o')

    plt.figure()
    plt.title('reduced misfit as function of A')
    plt.plot(As, reduced_error, 'o')

    plt.figure()
    plt.plot(mue, ras, 'o')

    plt.show()

play_2(loss_1, igauss)
