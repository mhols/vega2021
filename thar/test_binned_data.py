import scipy.special as sps
import numpy as np
import matplotlib.pyplot as plt


def gauss(x, A, mu, sigma, y_offset):
    return y_offset + (A/  (np.sqrt(2*np.pi)*sigma)) * np.exp(-(x-mu)**2/(2*sigma**2))

def igauss(x, A, mu, sigma, y_offset):
    return y_offset +  0.5 * A * ( \
        sps.erf((x+0.5-mu)/(np.sqrt(2)*sigma)) \
        - sps.erf((x-0.5-mu)/(np.sqrt(2)*sigma)))

nbins = 6
n = np.arange(nbins)

mu = 3.2
y_offset = 0.5
A = 1
sigma = 0.6

g = gauss(n, A, mu, sigma, y_offset) 
ig = igauss(n, A, mu, sigma, y_offset)

plt.plot(g, label='g')
plt.plot(ig, label='ig')
plt.legend()
plt.show()
