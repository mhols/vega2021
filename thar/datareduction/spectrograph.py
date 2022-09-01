import json
from util import *
import sys
import pandas as pd
import numpy as np
import scipy as sp
from scipy.optimize import bisect
from scipy.interpolate import interp1d
from numpy.polynomial import Polynomial
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
matplotlib.rcParams['figure.dpi'] = 200


class CCD2d:
    """
    #  resolution of ccd imaging
    """
    def __init__(self, data=None, **kwargs):
        self.kwargs = kwargs  # saving for later use
        if data is None:
            data = pd.DataFrame(pd.read_json(kwargs['datafile']))
        self.data = data

        if self.kwargs['bootstrap_data']:
            self.bootstrap_data()

        # basic outlier removal / quality filter
        self.outlier_removal_1()
        self.polynomial_fit = {}

    def generate_fixed_effects(self):
        Nlo = self.kwargs['Nlo']
        No = len(self.get_orders())
        N = max(Nlo, No)
        ol_range = ((self.o * self.l).min(), (self.o * self.l).max())
        o_range = (self.get_orders().min(), self.get_orders().max())
        I = np.identity(No + 1)
        Cheb_o = [np.polynomial.chebyshev.Chebyshev(c, window=[-1,1], domain=o_range) for c in I]
        lo_range = (self.get_orders().min(), self.get_orders().max())
        I = np.identity(Nlo + 1)
        Cheb_o = [np.polynomial.chebyshev.Chebyshev(c, window=[-1,1], domain=lo_range) for c in I]

    def bootstrap_data(self):
        res = []
        new_mean_pixel_pos = []
        sigma_new_mean = []
        shit = 0
        print('\nbootstrapping data...\n')
        for i, line in self.data.iterrows():
            position = np.nan
            sigma = np.nan
            try:
                delta_p, sigma = bootstrap_estimate_location(
                    np.array(line['flux_values_extract']),
                    **self.kwargs
                    )
                position = line['pixels_extract'][0] + delta_p

            except Exception as ex:
                print(i, line, ex)
                shit += 1
            new_mean_pixel_pos.append(position)
            sigma_new_mean.append(sigma)
            progress(i,len(self.data), status='')
        self.data['new_mean_pixel_pos'] = pd.Series(new_mean_pixel_pos)
        self.data['sigma_new_mean'] = pd.Series(sigma_new_mean)

        if self.kwargs['save_bootstraped_data']:
            self.data.to_json(self.kwargs['bootstraped_file'])

        ## print ('n negative ', shit)


    def color_of_order(self, o):
        cmap = cm.get_cmap(self.kwargs['palette_order'])
        orders = self.get_orders()
        colors = cmap(np.linspace(0,1, max(orders) - min(orders) +1))
        return colors[o-min(orders)]

    def outlier_removal_1(self):
        """
        removing outliers based on bootstrap uncertainty
        of histotram expectation
        """
        epsilon = self.kwargs['epsilon_sigma_bootstrap']
        self.data = self.data[self.data['sigma_new_mean'] < epsilon].copy().reset_index(drop=True)


    def eval_polynomial(self, ol, o):
        return np.array([ self.polynomial_fit[oo](ooll) for oo, ooll in zip(o, ol)])

    def sigma_clipping(self):
        # sigma clipping

        epsilon=self.kwargs['epsilon_sigma_clipp']
        thd = epsilon
        for i in range(self.kwargs['n_sigma_clip']):
            self.fit_polynomial_order_by_order()
            res = self.eval_polynomial(self.o*self.l, self.o) - self.x

            if self.kwargs['clipp_method'] == 'pixerror':
                thd = 1 * epsilon
            elif self.kwargs['clipp_method'] == 'rel_std':
                thd = self.sigma * epsilon
            else:
                thd = epsilon * np.std(res)
            I = np.abs(res)<=thd
            if np.all(I):
                break
            self.data = self.data[I].copy().reset_index(drop=True)

        epsilon=self.kwargs['epsilon_sigma_clipp_2d']
        for i in range(self.kwargs['n_sigma_clip']):
            self.fit_2d_polynomial()
            res = self.eval_polynomial(self.ol, self.o) - self.x
            if self.kwargs['clipp_method'] == 'pixerror':
                thd = 1 * epsilon
            elif self.kwargs['clipp_method'] == 'rel_std':
                thd = self.sigma * epsilon
            else:
                thd = epsilon * np.std(res)
            I = np.abs(res)<=thd
            if np.all(I):
                break
            self.data = self.data[I].copy().reset_index(drop=True)

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

    def lambda_at_x_o(self, x, o):
        nn = 10001
        olams = np.linspace(self.ol.min()/1.05, self.ol.max()*1.05, nn, endpoint=True)
        p = interp1d(
            self.polynomial_fit[o] (olams),
            olams)
        return p(x)/o

    def lambda_map_at_o(self, o):
        nn = 10001
        p = self.polynomial_fit[o]
        olams = np.linspace(self.ol.min(), self.ol.max(), nn, endpoint=True)
        der = p.deriv()
        I, = np.where(der(olams) > 0)
        II, = np.where(I[1:]-I[:-1] > 1)
        if len(II)==0:
            olmin, olmax = olams[I[0]], olams[I[-1]]
        else:
            III = np.argmax(II[1:]-II[:-1])
            olmin = olams[I[II[III]]]
            olmax = olams[I[II[III-1]]]

        olams = np.linspace(olmin, olmax, nn, endpoint=True)
        p = interp1d(self.polynomial_fit[o](olams), olams/o)
        return p

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

    @property
    def ol(self):
        return self.o * self.l

    def index_order(self, o):
        return self.data['true_order_number']==o

    def get_orders(self):
        orders = list(set(self.data['true_order_number']))
        orders = self.kwargs.get('orders', orders)
        return orders

    def get_report_orders(self):
        report_only_nonempty = self.kwargs.get('report_only_nonempty', False)
        if report_only_nonempty:
            return self.get_orders()
        return self.kwargs.get('report_orders', self.get_orders())

    def fit_global_polynomial(self, **kwargs):
        n = kwargs.get('order_ol', self.kwargs['order_ol'])
        ol = self.o * self.l
        domain = [ol.min(), ol.max()]
        p = None
        try:
            p = np.polynomial.chebyshev.Chebyshev.fit(ol, self.x,
                                                    domain=domain,
                                                    deg=n,
                                                    w = 1/self.sigma)

        except Exception as ex:
            print ('problem fitting global polynomial\n', ex)
        self.polynomial_fit = { o: p for o in self.get_report_orders()}

        return p

    def fit_polynomial_order_by_order(self, **kwargs):
        n = kwargs.get('order_ol', self.kwargs['order_ol'])
        res = {}
        ol = self.o * self.l
        domain = [ol.min(), ol.max()]
        orders = self.get_orders()
        for o in self.get_report_orders():
            if not o in orders:
                res[o] = None
                continue
            I = self.index_order(o)
            try:
                res[o] = np.polynomial.chebyshev.Chebyshev.fit(
                    self.l[I]*o, self.x[I],
                    domain=domain,
                    deg=n,
                    w = 1/self.sigma[I]
                )
            except Exception as ex:
                print('------------\nAt order {} the following error occured\n-------\n'.format(o),ex)
        self.polynomial_fit = res
        return res

    def fit_2d_polynomial(self):
        """
        # computes a polynomial proxy
        """

        olrange = [self.ol.min(), self.ol.max()]
        orange = [self.o.min(), self.o.max()]

        Nol = self.kwargs['order_ol']+1
        No  = self.kwargs['order_o'] +1

        nolko = [(n, k) for n in range(Nol) for k in range(No)]

        Tshebol = [np.polynomial.chebyshev.Chebyshev.basis(
                    window=[-1,1],
                    domain=olrange,
                    deg=n) for n in range(Nol)]

        Tshebo = [np.polynomial.chebyshev.Chebyshev.basis(
                    window=[-1,1],
                    domain=orange,
                    deg=n) for n in range(No)]

        G = np.zeros((self.ndata(), len(nolko)))
        for i in range(self.ndata()):
            for j in range(len(nolko)):
                nol, no = nolko[j]
                G[i,j] = Tshebol[nol](self.o[i]*self.l[i]) * Tshebo[no](self.o[i])

        sigma_min = self.kwargs['sigma_min']
        G = np.array((1./ np.sqrt(self.sigma**2+sigma_min**2)))[:, np.newaxis] * G
        coeffs = np.linalg.lstsq(G, (1./self.sigma) * self.x , rcond=None)[0]

        self.polynomial_fit = {
            o : sum([ coeffs[i] * Tshebol[nol] * Tshebo[no](o)
                for i, (nol, no) in enumerate(nolko)])
                     for o in self.get_report_orders()
        }

        return self.polynomial_fit

    def ndata(self):
        return len(self.o)

    def get_lambda_list(self):
        """
        2Dmap for all orders as returned by self.get_report_orders()
        saves to kwargs['file_lambda_list']
        """
        res = np.zeros(len(self.get_report_orders())*7800)
        i = 0
        x = np.arange(7800)
        for i, o in enumerate(self.get_report_orders()):
            p = self.lambda_map_at_o(o);
            for j, xx in enumerate(x):
                try:
                    res[i*7800+j] = p(xx)
                except:
                    res[i*7800+j] =  np.NaN
        np.savetxt(self.kwargs['file_lambda_list'], res)
        return res

    def get_lambda_map(self):
        tmp = self.get_lambda_list()
        return {o: tmp[i*7800:(i+1)*7800] for i, o in enumerate(self.get_report_orders())}


    def get_rms_per_spectrum(self):
        res = pd.Series(index=self.data.index)
        for o in self.get_orders():
            I = self.index_order(o)
            res.iloc[I] = 3e5 * (self.lambda_at_x_o(self.x[I], o)/self.l-1).abs()
        np.savetxt(self.kwargs['rms_report_file'], res)

class FP:
    """
    Fabry-Perrot modeling
    """
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.fpdata = json.load(open(kwargs['fp_file'],'r'))
        self.ccd = kwargs.get('ccd2d', CCD2d(**kwargs))
        self.ccd.sigma_clipping()

        self.data = pd.concat([ pd.DataFrame({
            'true_order_number': int(r['true_order_number']),
            'wave3': pd.Series(r['wave3']),
            'flux3': pd.Series(r['flux3'])
            }) for r in self.fpdata if len(r)>0]
        )

    def flux_at_order(self, o):
        return self.data[self.data['true_order_number']==o]['flux3'].values

    def lam_at_order(self, o):
        p = self.ccd.lambda_map_at_o(o)
        return p(self.data[self.data['true_order_number']==o]['wave3'].values)

    def extract_fp_peaks_at_order(self, o):
        f = self.flux_at_order(o)
        tf = pd.DataFrame({'flux':f})
        #f = tf['f'].array
        tf['locmax'] = False
        tf['locmax'] = np.logical_and((tf['flux'].diff() > 0).shift(), (tf['flux'].diff()<0))
        tf['locmax'] = tf['locmax'].shift(-1, fill_value=False)
        r = self.kwargs['fp_window']
        # eliminate points that are too close to the border
        tf['locmax'].iloc[:r] = False
        tf['locmax'].iloc[-r:] = False
        tf['locmax'][tf['flux']<self.kwargs['cutoff_fp_flux']] = False
        I = tf.index[tf['locmax']]
        snipp = tf['flux'].array[(I[:, np.newaxis] + np.arange(-r,r+1)).T]
        snapp = np.arange(-r,r+1)
        # local polynomial fit
        p = np.polyfit(snapp, snipp, deg=2)
        # shift of max position
        delta = -p[1]/(2*p[0])
        # possition of local max in pixel
        pixpos = I+delta
        tf['peak'] = np.NaN
        tf['peak'].iloc[I] = self.ccd.lambda_map_at_o(o)(pixpos)
        tf['peakf'] = 1./tf['peak']
        tf['lam']=self.ccd.lambda_map_at_o(o)(tf.index)
        tf['freq']=1./tf['lam']

        self.tf = tf

        return tf


if __name__=='__main__':

    kwargs = {
    "datafile": "NEO_20220219_173048_th2_voie1.json",
    'bootstrap_data': True,
    'n_bootstrap': 100,                             # number of bootstrap experiments
    'profile': 'gauss',
    'loss_function': 'loss_1',                     # weighted L2-loss
    'save_bootstraped_data': True,
    'bootstraped_file': 'ThAr2D_voie_1_new.json',
    'epsilon_sigma_bootstrap': 3,         # locations with larger uncertainty are removed
    'epsilon_sigma_clipp': 3,
    'epsilon_sigma_clipp_2d' : 2,
    'clipp_method' : 'rel_std',             # 'rel_std' or 'pixerror' or 'est_std'
    'orders': np.arange(23, 59+1),          # orders to use, set None if taken from data
    'palette_order': 'gist_rainbow',
    'order_ol': 6,
    'order_o': 7,
    'n_sigma_clip' : 100,
    'sigma_min' : 0.,                    # minimial sigma to avoid overfitting
    'file_lambda_list': 'arturo.dat',
    }

    data = CCD2d(**kwargs)
#    data.get_lambda_list()
    sys.exit(0)

