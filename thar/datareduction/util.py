import scipy.special as sps
import scipy.optimize as sop
import numpy as np
import sys
from numpy.polynomial import Polynomial
from scipy.interpolate import UnivariateSpline, LSQUnivariateSpline

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
    sigmamin = 0.0001

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
    if size <= 1:
        return estimate_location(intens, fun, g)[1], np.sqrt(sum(intens))
     
    p = intens / np.sum(intens)

    n = int(np.sum(intens))
    intens = np.random.multinomial(n, p, size)
    res = [ estimate_location(inte, fun, g)[1] for inte in intens]

    return np.mean(res), np.std(res)

def background(l, v, nnodes=5, q=0.3, qq=0.8, qqq=0.9):
    N_MAXITER = 100
    
    #s = UnivariateSpline(t, t, s=0)

    I = np.logical_not(np.isnan(v))
    II = np.logical_and(I, v<np.quantile(v[I],0.95))
    I = np.logical_and(II, v>np.quantile(v[I],0.05))

    #l = l[I]
    v = v[I]
    l = l[I] 

    #l = np.arange(len(v))
    t = np.linspace(l[0], l[-1], nnodes+2, endpoint=True)

    I = np.full(len(l), True)
    for i in range(N_MAXITER):
        try:
            p = LSQUnivariateSpline(l[I], v[I], t[1:-1])
        except:
            print('in background ', v, l, t)
            pass  # TODO log the event...
        I1 = p(l)-v <= np.quantile(p(l)-v, qq)
        I2 = p(l)-v >= np.quantile(p(l)-v, q)
        II = np.logical_and(I1, I2)

        s = np.std( p(l[I]) -v[I])
        II = np.abs( p(l) - v) < 3*s 
        if np.alltrue(II==I):
            break
        I = II
        print(i)

    d = []
    delt=t[2]-t[0]
    res = v - p(l)
    for tt, ttt in zip(t[:-4], t[4:]):
        J = np.logical_and(l>=tt, l<ttt)
        d.append( np.quantile(res[J], qqq))
    pp = UnivariateSpline(t[2:-2], d, s=0)
    return p  ####lambda x: pp(x) + p(x)

def pseudo_inverse(x, increasing=True):
    """
    computes the pseudo inversion of x
    """
    if not increasing:
        x = -x

    if len(x) <= 1:
        raise Exception('pseudo_inverse: too few elements')
    
    tmp = [x[0]]
    for a, b in zip(x[:-1], x[1:]):
        tmp.append(max(a, b))
    tmp = np.array(tmp)
    if not increasing:
        tmp = -tmp
    return tmp

    


# The MIT License (MIT)
# Copyright (c) 2016 Vladimir Ignatev
#
# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software 
# is furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
# OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()  # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)


if __name__ == '__main__':

    d = np.loadtxt('voi35.dat')
    l = np.arange(len(d))
    p = background(l,d,nnodes=25,q=0.4, qq=0.6, qqq=0.)

    plot(l, d)
    plot(l, p(l))