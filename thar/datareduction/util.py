import scipy.special as sps
import scipy.optimize as sop
from scipy.stats import wasserstein_distance as wd
import numpy as np
import sys
from units import *
from numpy.polynomial import Polynomial
from scipy.interpolate import interp1d, UnivariateSpline, LSQUnivariateSpline


def gauss(x, A, mu, sigma, y_offset):
    return y_offset + (A/  (np.sqrt(2*np.pi)*sigma)) * np.exp(-(x-mu)**2/(2*sigma**2))

def igauss(x, A, mu, sigma, y_offset):
    xx = np.zeros(x.shape[0]+2)
    xx[1:-1] = x
    xx[0] = x[0]-(x[1]-x[0])
    xx[-1] = x[-1] + (x[-1]-x[-2])

    d = 0.5 * (xx[1:]+xx[:-1])
    return y_offset * (d[1:]-d[:-1]) +  0.5 * A * (
        sps.erf((d[1:]  - mu)/(np.sqrt(2)*sigma))
        - sps.erf((d[:-1] - mu)/(np.sqrt(2)*sigma)))

def cauchy(x, A, mu, sigma, y_offset):
    return y_offset + A*sigma/(sigma**2 + (x-mu)**2)

def loss_1(params, *args):
    intens, func = args
    n = np.arange(intens.shape[0])
    i = func(n, *params)
    return np.abs(i-intens)*np.sqrt(np.abs(i))
    
def loss_2(params, *args):
    intens, func = args
    n = np.arange(intens.shape[0])
    i = func(n, *params)
    return intens-i

def loss_3(params, *args):
    intens, func = args
    n = np.arange(intens.shape[0])
    i = func(n, *params)
    return (i-intens)*np.sqrt(np.where(intens>1, intens, 1))


def loss_4(params, *args):
    lams, intens, func = args
    i = func(lams, *params)
    return (i-intens)*np.sqrt(np.where(intens>1, intens, 1))


def estimate_igauss(lams, intens):
    """
    simple gauss fit without weights...
    """

    # shift by mean....
    mlams = np.mean(lams)
    lams -= mlams

    A0 = 0.5 * (np.max(intens) - np.min(intens)) * (lams[-1]-lams[0]) / (lams[1]-lams[0])
    params0 = [-A0, 0, lams[-1], np.max(intens)/(lams[1]-lams[0])]


    res = sop.least_squares(
        loss_4, 
        params0,
        args = (lams, intens, igauss)

    )

    A, mu, sigma, y_offset = res.x

    return A, mlams + mu, sigma, y_offset
     

function_map = {
    'gauss': gauss,
    'igauss': igauss,
    'cauchy': cauchy,
    'loss_1': loss_1,
    'loss_2': loss_2,
    'loss_3': loss_3,
}

def estimate_location(intens,  absorption=False, **kwargs):

    intens = np.array(intens)

    if absorption:
        intens = -intens

    shift = np.min(intens)

    intens -= shift

    # backproject on positivity by clipping
    # intens = np.where(intens>0, intens, 0)
    # first guess of parameters
    n = np.arange(intens.shape[0])
    A = np.sum(np.abs(intens))
    Amin = max(0, A-4*np.sqrt(A))
    Amax = A+4*np.sqrt(A)

    mu = np.argmax(intens)
    mumin = 0
    mumax = len(intens)-1

    sigma = len(intens)
    sigmamax = 4*sigma
    sigmamin = 1e-6

    y_offset_min = -np.max(intens)
    y_offset_max = np.max(intens)
    y_offset = (y_offset_min+y_offset_max)/2

    params0 = np.array([A, mu, sigma, y_offset])

    bounds = (np.array([Amin, mumin, sigmamin, y_offset_min]),
              np.array([Amax, mumax, sigmamax, y_offset_max]))

    fun = function_map[kwargs['loss_function']]
    g = function_map[kwargs['profile']]
    res = sop.least_squares(
        fun, 
        params0, 
        bounds=bounds, 
        args=(intens, g))
    
    A, mu, sigma, y_offset = res.x

    return  A, mu, sigma, y_offset + shift


def bootstrap_estimate_location(intens, **kwargs):
    if len(intens) <=1:
        raise Exception('intens is empty')
    try:
        isn = np.any(np.isnan(intens))
    except Exception as ex:
        print(intens)
    
    if isn:
        raise Exception('intens not finite...')


    g = function_map[kwargs['profile']]
    size = kwargs['n_bootstrap']
        
    fun = function_map[kwargs['loss_function']]
    intens = intens - min(0, np.min(intens))
    
    if size <= 2:
        
        return estimate_location(intens, **kwargs)


    res = []
    for i in range(size):
        try:
            params = estimate_location(np.random.poisson(intens), **kwargs) 
        except:
            continue
        res.append(params)
    
    res = np.array(res)

    return np.mean(res[:,1]), np.std(res[:,1]), res

def air_to_vac(lam):
	#lambda_vac (see https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion)
	#lambda_air = lambda_vac/refindex
	s = 10**4/(lam/ANGSTROM)
	refindex = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
	return refindex
	
def vac_to_air(lam):
	#lambda_air (see https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion)
	#lambda_vac =  lmabda_air/refindex
	s = 10**4/(lam/ANGSTROM)
	refindex = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2)\
	+ 0.0001599740894897 / (38.92568793293 - s**2)
	return refindex
 
def continuum(l, v, nnodes=10, q=0.3, qq=0.8, qqq=0.9):
    N_MAXITER = 100
    
    #s = UnivariateSpline(t, t, s=0)

    if np.all(np.isnan(v)):
        return v
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

    d = []
    delt=t[2]-t[0]
    res = v - p(l)
    for tt, ttt in zip(t[:-4], t[4:]):
        J = np.logical_and(l>=tt, l<ttt)
        d.append( np.quantile(res[J], qqq))
    pp = UnivariateSpline(t[2:-2], d, s=0)
    return p 

def pseudo_inverse(x):
    """
    computes the pseudo based on the largest 
    interbal on which x is monotoneous
    """
    if len(x) <= 1:
        raise Exception('pseudo_inverse: too few elements')
    
    mc, inc, dec = monotoneous_chunks(x) 
    

def local_maxima(v):
    v = np.array(v)
    vv = np.zeros(len(v)+2)
    vv[1:-1] = v
    d = np.diff(vv)
    lma = np.where(np.logical_and(d[:-1] > 0, d[1:] < 0))[0] 
    lmi = np.where(np.logical_and(d[:-1] < 0, d[1:] > 0))[0] 

    
    tmp = []
    for im in lma:
        try:
            a = np.max(lmi[lmi<im])
        except:
            a = 0
        try:
            b = np.min(lmi[lmi>im])
        except Exception as ex:
            b = len(v)-1

        tmp.append([im, a, b])
    
    tmp = sorted(tmp, key=lambda r: -v[r[0]])
    return np.array(tmp)

def monotoneous_chunks(v):

    lm = local_maxima(v)
    chunks = np.row_stack((lm[:, [1,0]], lm[:, [0,2]]))
    I = np.argsort(-np.abs(chunks[:,0]- chunks[:,1]))
    chunks = chunks[I]
    chunks = chunks[chunks[:,0] < chunks[:,1]]

    growing = v[chunks[:,0]] < v[chunks[:,1]]
    decreasing = v[chunks[:,0]] > v[chunks[:,1]]

    return chunks, growing, decreasing

def sigma_clipping_general_map(fitmachine, clipmachine, I0):
    """
    fitmachine (x, y, I) -> fitted 
    clipmachine(fitted) -> indices retained
    """
    NMAX = 200

    p = fitmachine(I0)
    I = clipmachine(p)
    # I = np.logical_and(I0, I)

    nclip = 1
    for i in range(NMAX):
        p = fitmachine(I)
        #II = np.logical_and( I, clipmachine(p) )
        II = clipmachine(p)
        if np.all(I == II):
            break
    
        I = II
        nclip += 1
    
    return p, I, nclip



def homothetie_wasserstein(x,y,xx,yy, rb0, rb1, ra0, ra1):
    """
    best homothetie based on Wasserstein distance

    """
    params0 = np.array([(rb0+rb1)/2, (ra0 + ra1)/2])
    bounds = (np.array([rb0, ra0]), np.array([rb1, ra1])) 
    res = sop.least_squares(
        lambda ba: wd(x, ba[1]*xx+ba[0], y, yy ),
        x0=params0, 
        bounds=bounds 
    )
    return res.x
   
def clean_nans(x, default=0):
    return np.where(np.isnan(x), default, x)

def homothetie(x,y, xx, yy, rb0, rb1, ra0, ra1):
    """
    returns mapping f: x -> xx so that
    y(f^{-1}xx) = yy(xx)  <==> yy(f(x)) = y(x)
    
    f is searched for in the space of affine mappings
    """
    N =  10*(len(x) + len(xx))
    yx = interp1d(x,y,fill_value='extrapolate')
    yyxx = interp1d(xx, yy, fill_value='extrapolate')

    minx = max(np.min(x), np.min(xx*ra0+rb0))
    maxx = min(np.max(x), np.max(xx*ra1+rb1))

    xxx = np.linspace(minx, maxx, N)
    def L(ba):
        b, a = ba
        return yx(xxx) - yyxx( (xxx-b)/a) 

    params0 = np.array([(rb0+rb1)/2, (ra0 + ra1)/2])
    bounds = (np.array([rb0, ra0]), np.array([rb1, ra1])) 
    res = sop.least_squares(
        L, 
        x0=params0, 
        bounds=bounds 
    )
    return res.x
   

class MonotoneFunction:
    """
    operations on monotone relations x <--> y
    the function are linearly extendet to +/- infinity
    """
    def __init__(self, a, b, f, df):
        """
        on the interval [a, b] the function f should be strictly monoton
        """
        assert a < b, 'a should be smaller than b'
        self.f = f
        self.df = df
        self.a, self.b = a, b
        self.ra, self.rb = min(f(a), f(b)), max(f(a), f(b))

        # the toggle between the direct evaluation or its inverse
        self._direct = True

    @property
    def is_growing(self):
        return self.df(0.5 * (self.a + self.b)) > 0    

    @property
    def domain(self):
        return (self.a, self.b) if self._direct else (self.ra, self.rb)

    @property 
    def range(self):
        return (self.ra, self.rb) if self._direct else (self.a, self.b)


    def deriv(self):
        if not self._direct:
            raise Exception('for inverse function on yet implementet')
        return self._df
    

    def _der_at_lower_upper(self):
        """
        derivative at lower and larger part of range
        """
        return (self.df(self.a), self.df(self.b)) if self.is_growing else (self.df(self.b), self.df(self.a))


    def _inverse(self, y):

        r1, r2 = self.ra, self.rb
        mi, ma = np.min(y), np.max(y)
        d1, d2 = self._der_at_lower_upper()

        
        xmi = self.a + (mi - r1)/d1 if mi < r1 else self.a
        xma = self.b + (ma - r2)/d2 if ma > r2 else self.b

        n = 10001
        x = np.linspace(xmi, xma, n)
        f = self._call(x)

        xx = np.interp(y, f, x) if self.is_growing else np.interp(y, f[::-1], x[::-1])

        # some Newton iterations....
        for i in range(10):
            xx = xx - (self._call(xx) - y) / self._df(xx)
        
        return xx
       
    def _df(self, x):
        a, b = self.a, self.b
        dfa = self.df(a)
        dfb = self.df(b)
        I1, = np.where( x < a)
        I2, = np.where( x > b)
        I3, = np.where( np.logical_and (x >= a, x <= b))
        res = np.zeros_like(x)
        res[I3] = self.df(x[I3])
        res[I1] = dfa
        res[I2] = dfb

        return res

        
    def _call(self, x):
        a, b = self.a, self.b
        dfa = self.df(a)
        dfb = self.df(b)
        I1, = np.where( x < a)
        I2, = np.where( x > b)
        I3, = np.where( np.logical_and (x >= a, x <= b))
        res = np.zeros_like(x)
        res[I3] = self.f(x[I3])
        res[I1] = self.f(a) + dfa * ( x[I1] - a)
        res[I2] = self.f(b) + dfb * ( x[I2] - b)

        return res


    def inverse(self):
        tmp = MonotoneFunction(self.a, self.b, self.f, self.df)
        tmp._direct = not self._direct
        return tmp

    def deriv(self):
        return self._df
        
    def __call__(self, x):
        return self._call(x) if self._direct else self._inverse(x)
   
#--------------
# Matching sequeces
#------------
def matching(v, w, dvw, dwv=None):
    """
    v: iteratble
    w: iterable of same length
    dvw: match function v -> w
    dwc: match function w -> v
    """

    if dwv is None:
        dwv = dvw
    m1 =  dvw(v[:, np.newaxis], w[np.newaxis,:])
    m2 =  dwv(v[:, np.newaxis], w[np.newaxis,:])

    mm = m1 & m2
    mm = mm & (np.count_nonzero(mm, axis=1)==1)[:, np.newaxis] \
            & (np.count_nonzero(mm, axis=0)==1)[np.newaxis, :]
    I, J = np.where( mm )

    return I, J


def extract_snippets(l, v, lmin, lmax):
    ll = []
    vv = []
    for l1, l2 in zip(lmin, lmax):
        J = (l1 <= l) & (l <= l2)
        if sum(J) > 0:
            ll.append(l[J])
            vv.append(v[J])
    return ll, vv




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
    if count >= total:
        sys.stdout.write('\n')


if __name__ == '__main__':

    import matplotlib.pyplot as plt    
    p = np.polynomial.Polynomial([-1,-0.5, -0.3])

    pp = MonotoneFunction(0,4, p, p.deriv(1))
    pi = pp.get_inverse()

    x = np.linspace(-1, 5, 1000)
    y = pp(x)
    #plt.plot(x, pp(x))

    plt.plot(x, pi(pp(x))-x)
    plt.plot(x, pp(pi(x))-x)

    plt.show()