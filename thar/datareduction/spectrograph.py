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
## matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm

matplotlib.rcParams['figure.dpi'] = 200

from settings import *

##  local units
M = 1.0         # Meter
S = 1.0         # Second
PIXEL = 1       # Pixel unit
HZ = 1/S        # Herz
KILO = 1000     # Kilo
KM = KILO * M   # Kilometer

C_LIGHT = 300000 * KM / S    # speed of ligth
##  -------------


def polyfit_2d(ol, o, x, w, **kwargs):
    olrange = kwargs['olrange']
    orange = kwargs['orange']
    Nol = kwargs['order_ol']+1
    No  = kwargs['order_o'] +1

    nolko = [(n, k) for n in range(Nol) for k in range(No)]

    Tshebol = [np.polynomial.chebyshev.Chebyshev.basis(
                    window=[-1,1],
                    domain=olrange,
                    deg=n) for n in range(Nol)]

    Tshebo = [np.polynomial.chebyshev.Chebyshev.basis(
                    window=[-1,1],
                    domain=orange,
                    deg=n) for n in range(No)]
    
    ndata = len(x)
    G = np.zeros((ndata, len(nolko)))
    sigma_min = kwargs['sigma_min']
    
    s = np.sqrt(s**2 + sigma_min**2)

    for i in range(ndata):
        for j, (nol, no) in enumerate(nolko):
            G[i,j] = Tshebol[nol](ol[i]) * Tshebo[no](o[i])

    G = (1./ s)[:, np.newaxis] * G
    coeffs = np.linalg.lstsq(G, w * x )[0]

    polynomial_fit = {
        o : sum([ coeffs[i] * Tshebo[no](o) * Tshebol[nol] 
            for i, (nol, no) in enumerate(nolko)])
                for o in orange 
    }
    return polynomial_fit





class CCD2d:
    """
    class for wavemap computation (i.e. 2D polynomial and friends)
    """
    def __init__(self, data=None, **kwargs):
        self.kwargs = kwargs  # saving for later use
        if data is None:
            data = pd.DataFrame(pd.read_json(kwargs['datafile']))
        self.data = data

        self.data['selected'] = True # we are using only this subset
        total_flux = np.array([sum(flux) for flux in self.data['flux_values_extract']])
        self.data['total_flux'] = total_flux

        if self.kwargs.get('bootstrap_data', False):
            self.bootstrap_data()

        # basic outlier removal / quality filter
        self.outlier_removal_1()
        self.polynomial_fit = {
            o: np.polynomial.Polynomial([0]) 
            for o in self.all_order() 
        }

    def generate_fixed_effects(self):
        """
        for later use when replacing the 2d polynomial with 2d splines
        """
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
        """
        bootstrap estimate of location uncertainty
        """
        res = []
        new_mean_pixel_pos = []
        sigma_new_mean = []
        shit = 0
        print('\nbootstrapping data...\n')   ## TODO: better logging
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
                print("problem at ", i, line, ex)  ## TODO: better logging
                shit += 1
            new_mean_pixel_pos.append(position)
            sigma_new_mean.append(sigma)
            progress(i,len(self.data), status='')
        self.data['new_mean_pixel_pos'] = pd.Series(new_mean_pixel_pos)
        self.data['sigma_new_mean'] = pd.Series(sigma_new_mean)

        if self.kwargs.get('save_bootstraped_data', False):
            try:
                self.data.to_json(self.kwargs['bootstraped_file'])
            except:
                raise Exception('bootstraped_file has not not been specified')
        
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

    def eval_dpolyomial_dl(self, ol, o):
        return np.array([ oo*self.polynomial_fit[oo].deriv()(ooll) for oo, ooll in zip(o, ol)])
    
    def sigma_clipping(self):
        # sigma clipping

        # order by order fitting first
        epsilon = self.kwargs['epsilon_sigma_clipp']
        thd = epsilon 
        o = self.data['true_order_number']
        l = self.data['ref_lambda']
        x = self.data['new_mean_pixel_pos']
        sigma = self.data['sigma_new_mean']

        for i in range(self.kwargs['n_sigma_clipp']):

            self.fit_polynomial_order_by_order()
            res = self.eval_polynomial(o*l, o) - x
            dxdl = self.eval_dpolyomial_dl(o*l, o)

            # how to define the clipping
            if self.kwargs['clipp_method'] == 'pixerror':   # absolute error in pixel
                thd = 1 * epsilon
            elif self.kwargs['clipp_method'] == 'rel_std':  # relative error in std of each pixel
                thd = sigma * epsilon
            elif self.kwargs['clipp_method'] == 'vrad':     # radial velocity uncertainty
                thd = l * dxdl * epsilon / C_LIGHT
            else:
                raise Exception ('no such method') 
                # thd = epsilon * np.std(res)
            # clipping 
            I = np.abs(res)<=thd

            if np.all(I==self.data['selected']):
                break
            self.data.loc[:, 'selected'] = I[:]

        epsilon=self.kwargs['epsilon_sigma_clipp_2d']


        for i in range(self.kwargs['n_sigma_clipp_2d']):
            self.fit_2d_polynomial()
            res = self.eval_polynomial(o*l, o) - x
            dxdl = self.eval_dpolyomial_dl(o*l, o)

            if self.kwargs['clipp_method'] == 'pixerror':
                thd = 1 * epsilon
            elif self.kwargs['clipp_method'] == 'rel_std':
                thd = sigma * epsilon
            elif self.kwargs['clipp_method'] == 'vrad':     # radial velocity uncertainty
                thd = l * dxdl * epsilon / C_LIGHT
            else:
                thd = epsilon * np.std(res)
            I = np.abs(res)<=thd
            if np.all(I==self.data['selected']):
                break
            if len(I) < 10:
                epsilon *= 1.5
            self.data.loc[:, 'selected'] = I

        if self.kwargs.get('save_sigma_clipped_data', False):
            self.data.to_json(self.kwargs['sigma_clipped_file'])

    def _delta_lam_over_lam_factor(self, l, o):
        """

        dl / l = (dl / dx) * x / l * dx = (dx / dl)^{-1} * x / l * dx
        at o we have x = p(ol) => dx/dl = p'(ol) * o
        returns (dx / dl)^{-1} * x / l 
        """
        p = self.polynomial_fit[o]
        dp = p.deriv()
        x = p(o*l)
        factor = (o * dp(o*l))**(-1) * x / l
        return factor


    def delta_radial_velocity_model(self):
        res = {}
        for o in self.all_order():
            I = self.index_order(o)
            p = self.polynomial_fit[o]
            dp = p.deriv()
            dxdl = o*dp(o*self.l[I])
            dx = p(o*self.l[I]) - self.x[I]
            factor = dxdl**(-1) / self.l[I]
            tmp = C_LIGHT * factor * dx / (M/S) 
            res[o] = np.sqrt(np.mean(tmp)**2+np.var(tmp))
        return res

    def delta_l_over_l(self):
        lmaps = self.lambda_maps()
        return np.array([(lmaps[o](x)/l-1) for o, l,x in zip(self.o, self.l, self.x)])

    def delta_pixel(self):
        return self.eval_polynomial(self.ol, self.o) - self.x

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

    def lambda_maps(self):
        return {o: self.lambda_map_at_o(o) for o in self.all_order()}

    @property
    def x(self):
        return self.data['new_mean_pixel_pos'][self.data['selected']]

    @property
    def sigma(self):
        return self.data['sigma_new_mean'][self.data['selected']]

    @property
    def l(self):
        return self.data['ref_lambda'][self.data['selected']]

    @property
    def o(self):
        return self.data['true_order_number'][self.data['selected']]

    @property
    def ol(self):
        return self.o * self.l

    @property
    def total_flux(self):
        return self.data['total_flux'][self.data['selected']]

    @property
    def index(self):
        return self.data.index[self.data['selected']]

    def index_order(self, o):
        return self.data['true_order_number'][self.data['selected']]==o

    def index_order_unselected(self, o):
        return self.data['true_order_number']==o

    def all_order(self):
        return self.kwargs.get('ORDERS', ORDERS)

    def get_orders(self):
        orders = list(set(self.data['true_order_number'][self.data['selected']]))
        orders = self.kwargs.get('orders', orders)
        return orders

    def get_report_orders(self):
        report_only_nonempty = self.kwargs.get('report_only_nonempty', False)
        if report_only_nonempty:
            return self.get_orders()
        return self.kwargs.get('report_orders', self.get_orders())

    def get_global_polynomial(self):
        n = self.kwargs['order_ol'] 
        ol = self.ol
        domain = [ol.min(), ol.max()]
        p = None
        try:
            p = np.polynomial.chebyshev.Chebyshev.fit(ol, self.x,
                                                    domain=domain,
                                                    deg=n,
                                                    w=self.fit_weight())
        except:
            raise Exception ('problem fitting global polynomial\n')
        return p

    def fit_global_polynomial(self):
        p = self.get_global_polynomial()
        self.polynomial_fit = {
            o: p # np.polynomial.Polynomial([0]) 
            for o in self.all_order()
        }
        return p

    def len_total_data(self):
        return len(self.data['true_order_number'])

    def fit_weight(self):

        match self.kwargs.get('fitweight', 'sigma'):
            case 'equal':
                weight = 1
            case 'sigma':
                weight = 1/self.sigma
            case 'vrad':
                weight = np.ones(self.ndata)
                for o in self.all_order():
                    I = self.index_order(o)
                    weight[I] = 1/np.abs((self.ol[I]*self.polynomial_fit[o].deriv(1)(self.ol[I])))
        return weight      

    def fit_polynomial_order_by_order(self):
        n = kwargs['order_ol']
        weight = self.fit_weight()
        res = {}
        ol = self.o * self.l
        domain = [ol.min(), ol.max()]
        orders = self.all_order()

        pglobal = self.fit_global_polynomial()

        for o in orders:
            I = self.index_order(o)
            try:
                res[o] = np.polynomial.chebyshev.Chebyshev.fit(
                    self.l[I]*o, self.x[I],
                    domain=domain,
                    deg=n,
                    w = weight[I]
                )
            except Exception as ex:
                res[o] = pglobal
                
                #print('------------\nAt order {} the following error occured\n-------\n'.format(o),ex)
        
        self.polynomial_fit = res
        return res

    def fit_2d_polynomial(self, **kwargs):
        """
        # computes a polynomial proxy
        """

        self.kwargs.update(kwargs)
        o = self.data['true_order_number']
        l = self.data['ref_lambda']

        olrange = [(o*l).min(), (o*l).max()]
        orange = [o.min(), o.max()]
        self.polynomial_fit = polyfit_2d(
            self.ol.to_numpy(),
            self.o.to_numpy(),
            self.x.to_numpy(),
            self.fit_weight().to_numpy(),
            olrange=olrange,
            orange=orange,
            **self.kwargs
        )
        return self.polynomial_fit
        """
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

        G = np.zeros((self.ndata, len(nolko)))
        o = self.o.to_numpy()
        l = self.l.to_numpy()
        x = self.x.to_numpy()
        s = self.sigma.to_numpy()
        s = np.sqrt(s**2 + self.kwargs['sigma_min']**2)

        for i in range(self.ndata):
            for j, (nol, no) in enumerate(nolko):
                G[i,j] = Tshebol[nol](o[i]*l[i]) * Tshebo[no](o[i])

        G = (1./ s)[:, np.newaxis] * G
        coeffs = np.linalg.lstsq(G, (1./s) * x )[0]

        ## compute restriction to each order of 2d polynomial
        self.polynomial_fit = {
            o : sum([ coeffs[i] * Tshebo[no](o) * Tshebol[nol] 
                for i, (nol, no) in enumerate(nolko)])
                     for o in self.all_order() 
        }

        return self.polynomial_fit
        """
    @property
    def ndata(self):
        return len(self.o)

    def get_lambda_list(self):
        """
        2Dmap for all orders as returned by self.get_report_orders()
        saves to kwargs['file_lambda_list']
        """
        res = np.zeros(len(self.get_report_orders())*NROWS)
        i = 0
        x = np.arange(NROWS)
        for i, o in enumerate(self.get_report_orders()):
            p = self.lambda_map_at_o(o);
            for j, xx in enumerate(x):
                try:
                    res[i*NROWS+j] = p(xx)
                except:
                    res[i*NROWS+j] =  np.NaN
        np.savetxt(self.kwargs['file_lambda_list'], res)
        return res

    def get_lambda_map(self):
        tmp = self.get_lambda_list()
        return {o: tmp[i*NROWS:(i+1)*NROWS] for i, o in enumerate(self.get_report_orders())}


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
    'n_sigma_clipp' : 100,
    'n_sigma_clipp_2d' : 100,
    'sigma_min' : 0.,                    # minimial sigma to avoid overfitting
    'file_lambda_list': 'arturo.dat',
    }

    data = CCD2d(**kwargs)
#    data.get_lambda_list()
    sys.exit(0)
