import json
from util import *
import sys
import pandas as pd
import numpy as np
from scipy.optimize import bisect
from numpy.polynomial import Polynomial
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.dpi'] = 200

## put the parameters of the probram here...
default_kwargs = {
    "datafile": "ThAr2D_voie_1_new.json",
    'epsilon_sigma_bootstrap': 2.0,
    'epsilon_sigma_clipp': 8,
    'epsilon_sigma_clipp_2d' : 5,
    'orders': np.arange(23, 59+1),
    'order_n': 6,
    'order_o': 7,
    'crossmax': 10,
    'n_sigma_clip' : 5,
    'file_lambda_list': 'arturo.dat',
}


def prepare_jsons():
    global data1, data2
    with open("./ThAr2D_voie_1.dat.json", "r") as f:
        data = json.load(f)

    res = []
    shit = 0
    for l in data:
        try:
            dp, sigma = bootstrap_estimate_location(
                np.array(l['flux_values_extract']), 
                loss_1, # the loss function loss_2 L2
                gauss   # the flux model
                )
            p = l['pixels_extract'][0] + dp

        except Exception as ex:
            print(l, ex)
            shit += 1
            plt.figure()
            plt.plot(l['flux_values_extract'])
            plt.show()
            continue

        l['new_mean_pixel_pos'] = p
        l['sigma_new_mean'] = sigma
        res.append(l)

    with open('ThAr2D_voie_1_new.json', 'w') as f:
        json.dump(res, f)

    print ('n negative ', shit)

def link_ol(ol):
    return ol/1000

def fit_polynome_direct_order_by_order(data, **kwargs):
    orders = get_orders(data, **kwargs)
    porder_ol = get_kwarg(kwargs, 'porder_ol', 6)

    eps = pd.Series(index=data.index, dtype=float)
    for o in orders:
        I = data['true_order_number']==o
        s = data[I]['sigma_new_mean'].to_numpy()
        l = data[I]['ref_lambda'].to_numpy()
        x = data[I]['new_mean_pixel_pos'].to_numpy()
        p = np.polyfit(link_ol(o*l), x, porder_ol, w=1./s)
        eps[I] = np.abs((np.polyval(p, link_ol(o*l)) - x))/s
    return eps

def fit_2d_polynome(data, **kwargs):
    """
    two d polynomial fit
    """
    P = PolySets()
    P.estimate_polynome(data)
    print(P.coeffs)
    o = data['true_order_number']
    l = data['ref_lambda']
    x = data['new_mean_pixel_pos']
    s = data['sigma_new_mean']

    eps = np.abs((P(o*l, o) - x))/s
    return eps

def load_data1(**kwargs):
    """
    returns data1 with simple quality filter (sigma < epsilon)
    """
    data1 = pd.DataFrame(pd.read_json(kwargs['datafile']))
    epsilon = kwargs['epsilon_sigma_bootstrap']
    data1 = data1[data1['sigma_new_mean'] < epsilon]

    return data1

class PolySets:
    """
    two 2 polynomial interpolation class
    highly NON optimized.....
    """
    def __init__(self, **kwargs):
        self._2dpoly(**kwargs)

    def _2dpoly(self, **kwargs):
        self.n=get_kwarg(kwargs, 'order_n',  6)
        self.o=get_kwarg(kwargs, 'order_o', 7)
        self.crossmax=get_kwarg(kwargs, 'crossmax', 6)

        res=[]
        for i in range (self.n+1):
            for j in range (self.o+1):
                if (i+j <= self.crossmax) | (i==0) | (j==0):
                    res.append((i,j))
        self._monome = res

    def nfunc(self):
        return len(self._monome)

    def func(self, i, ol, o):
        n,m = self._monome[i]
        return ol**n*o**m

    def dxfunc(self,i, ol, o):
        """
        hier nur Ableitung nach lambda (ol), da o bekannt sind
        """
        n,m = self._monome[i]
        return n*ol**(n-1)*o**m

    def estimate_polynome(self, data, **kwargs):

        n_data = data.shape[0]
        n_poly = self.nfunc()

        F = np.zeros((n_data, n_poly))
        i=0
        for k, row in data.iterrows():
            for j in range(n_poly):
                o = row['true_order_number']
                l = row['ref_lambda']
                F[i,j] = self.func(j,link_ol(o*l), o)
            i+=1
        x = data['meas_pixel_position']
        w = 1./ data['sigma_new_mean']
        Fw = F * np.sqrt(w)[:, np.newaxis]
        xw = x * np.sqrt(w)
        self.coeffs = np.linalg.lstsq(Fw, xw)[0]

    def __call__(self, ol, o):
        """
        evaluation at (ol, o)
        """
        return sum (
            [self.coeffs[i]*self.func(i, link_ol(ol), o) for i in range(self.nfunc())]
        )

class CCD2d:
    """
    resolution of ccd imaging
    """
    def __init__(self, data=None, **kwargs):
        self.kwargs = kwargs  # saving for later use
        self.P = PolySets(**kwargs)  # the polynomials
        if data is None:
            data = pd.DataFrame(pd.read_json(kwargs['datafile']))
        epsilon = kwargs.get('epsilon_sigma_bootstrap', 1.0)
        data = data[data['sigma_new_mean'] < epsilon]
        self.data = data

        # sigma clipping
        self.P.estimate_polynome(self.data, **kwargs)

    def sigma_clipped(self, **kwargs):
        self.kwargs.update(kwargs)
        # sigma clipping
        data = sigma_clipping(self.data,
                          epsilon=self.kwargs['epsilon_sigma_clipp'],
                          fit=fit_polynome_direct_order_by_order, **self.kwargs)
        data = sigma_clipping(data,
                          epsilon=self.kwargs['epsilon_sigma_clipp_2d'],
                          fit=fit_2d_polynome, **self.kwargs)
        return CCD2d(data, **self.kwargs)

    def bare_x_to_lambda(self, x, o):
        """
        solution of implicit equation for given order to
        recover lambda from x
        """
        ol = self.o * self.l
        p = np.polyfit(self.x, ol, 1)
        ol0 = np.polyval(p, x)
        i=1
        while (self.P(ol0-i,o)-x)*(self.P(ol0+i,o)-x)>0:
            i+=1
        lam = bisect(
            lambda ol: self.P(ol, o)-x,
            ol0-i,
            ol0+i,
        ) / o

        return lam

    @property
    def x(self):
        return self.data['new_mean_pixel_pos']

    @property
    def sigma(self):
        return self.data['sigma_new_mean']

    @property
    def l(self):
        return self.data['ref_lambda']

    @property
    def o(self):
        return self.data['true_order_number']

    def lambda_to_x(self, l, o):
        return self.P(o*l, o)

    def index_order(self, o):
        return self.data['true_order_number']==o

    def get_orders(self):
        orders = list(set(self.data['true_order_number']))
        orders = get_kwarg(self.kwargs, 'orders', orders)
        return orders

    def polynomes_order_per_order(self):
        n = get_kwarg(self.kwargs, 'n_poly', 6)
        res = {}
        for o in self.get_orders():
            I = self.index_order(o)
            res[o] = np.polyfit(self.l[I]*o, self.x[I], w = 1/self.sigma[I])
        return res

    def rms(self):
        return [ 3e8 * (self.bare_x_to_lambda(self.x[i],self.o[i])/self.l[i]-1)
                for i in self.data_sigma_clipped.index]

    def get_lambda_list(self):
        res = np.zeros(39*7800)
        i = 0
        for o in range(21, 60):
            print('working on order ', o)
            for x in range(7800):
                res[i] = self.bare_x_to_lambda(x, o)
                i +=1
        np.savetxt(self.kwargs['file_lambda_list'], res)

def get_orders(data, **kwargs):
    orders = list(set(data['true_order_number']))
    orders = get_kwarg(kwargs, 'orders', orders)
    return orders


def plot_4(**kwargs):
    data = load_data1(**kwargs)
    data = sigma_clipping(data,
        epsilon=get_kwarg(kwargs, 'epsilon_sigma_clipp', 3),
                          fit=fit_polynome_direct_order_by_order, **kwargs)

    o0 = 45                       # reference order for plotting
    orders = get_orders(data, **kwargs)

    # reference polynomial
    ol = data['true_order_number'] * data['ref_lambda']
    x = data['new_mean_pixel_pos']
    I = data['true_order_number']==o0
    s = data['sigma_new_mean']
    p0 = np.polyfit(link_ol(ol[I]), x[I], 4, w=1./s[I])

    plt.figure()

    for o in orders:
        I = data['true_order_number']==o
        ol_continuum =np.linspace(ol[I].min(), ol[I].max(), 1024) 
        p = np.polyfit(link_ol(ol[I]), x[I], 4, w=1./s[I])
        plt.plot(ol_continuum, 
                 np.polyval(p-p0, link_ol(ol_continuum)))
        plt.plot(ol[I], x[I] - np.polyval(p0, link_ol(ol[I])), 'o', label=o)
    plt.legend()

def plot_5(**kwargs):
    # data = load_data1(**kwargs)

    CCD = CCD2d(**kwargs).sigma_clipped()
    data = CCD.data
    orders = get_orders(data, **kwargs)

    l = data['ref_lambda']
    ol = l * data['true_order_number']
    x = data['new_mean_pixel_pos']


    fig = plt.figure()
    ax = fig.add_subplot(111)

    for o in orders:
        I = data['true_order_number']==o
        l_continuum =np.linspace(ol.min(), ol.max(), 1024)
        line, = ax.plot(l_continuum/o, CCD.lambda_to_x(l_continuum/o, o))
        ax.plot(l[I], x[I], 'o', label=o, color = line.get_color())
    plt.legend()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    p0 = np.polyfit(ol,x,4)

    for o in orders:
        I = data['true_order_number']==o
        l_continuum =np.linspace(l[I].min(), l[I].max(), 1024)
        line, = ax.plot(o*l_continuum, CCD.lambda_to_x(l_continuum, o)- np.polyval(p0, o*l_continuum))
        plt.plot(o*l[I], x[I]-np.polyval(p0, o*l[I]), 'o', label=o, color=line.get_color())
    plt.legend()

    o = 25
    lam =CCD.bare_x_to_lambda(4000, o)
    print( 'lam...', lam, CCD.lambda_to_x(lam, o))
    print( CCD.rms() )

def pilote_1(**kwargs):
    CCD = CCD2d(**kwargs).sigma_clipped()
    CCD.get_lambda_list()
    
if __name__=='__main__':
    prepare_jsons()
    sys.exit(0)
    """plot_1()
    #plot_2(int(sys.argv[1]))
    plt.show()
    """

    #coeff, res = estimate_polynome()
    # print (np.std(res))
    # print (np.mean(res))
    kwargs = default_kwargs
    kwargs.update({
        "datafile": "ThAr2D_voie_1_new.json",
        'epsilon_sigma_bootstrap': 2.0,
        'epsilon_sigma_clipp': 8,
        'epsilon_sigma_clipp_2d' : 5,
        'orders': np.arange(26, 54),
        'order_n': 6,
        'order_o': 7,
        'crossmax': 7,
        'n_sigma_clip' : 5,
    })
    # plot_4(**kwargs)
    plot_5(**kwargs)
    # pilote_1(**kwargs)

    plt.show()

