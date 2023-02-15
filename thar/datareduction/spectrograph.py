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


"""
def polyfit_2d(ol, o, x, w, **kwargs):
    olrange = kwargs['olrange']
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
    
    w = np.sqrt(w**2 + sigma_min**2)

    for i in range(ndata):
        for j, (nol, no) in enumerate(nolko):
            G[i,j] = Tshebol[nol](ol[i]) * Tshebo[no](o[i])

    G = (1./ w)[:, np.newaxis] * G
    coeffs = np.linalg.lstsq(G, w * x )[0]

    polynomial_fit = {
        o : sum([ coeffs[i] * Tshebo[no](o) * Tshebol[nol] 
            for i, (nol, no) in enumerate(nolko)])
                for o in orange 
    }
    return polynomial_fit
"""

def to_list(x):
    if hasattr(x, '__iter__'):
        return x
    else:
        return [x]

def inverse_map(p):   
    nn = 10001
    mi, ma = p.domain
    x = np.linspace(mi, ma, nn, endpoint=True)
    der = p.deriv(1)
    I, = np.where(der(x) > 0)
    II, = np.where(I[1:]-I[:-1] > 1)
    if len(II)==0:
        xmin, xmax = x[I[0]], x[I[-1]]
    else:
        III = np.argmax(II[1:]-II[:-1])
        xmin = x[I[II[III]]]
        xmax = x[I[II[III-1]]]

    x = np.linspace(xmin, xmax, nn, endpoint=True)
    y = p(x)
    mi, ma = y.min(), y.max()
    return interpolate_extend(y, x, mi, ma)
       

class interpolate_extend:
    
    def __init__(self, x, y, mi, ma):
        self._mi, self._ma = np.min(x), np.max(x)
        self._interp = interp1d(x,y, fill_value="extrapolate")

    @property    
    def domain(self):
        return np.array([self._mi, self._ma])

    def inverse(self):
        return inverse_map(self)

    def __call__(self, x):
        return np.where((x-self._mi)*(self._ma-x)>0, self._interp(x), 0)


class CCD2d:
    """
    class for wavemap computation (i.e. 2D polynomial and friends)
    """
    def __init__(self, data=None, **kwargs):
        self.kwargs = kwargs  # saving for later use
        if data is None:
            try:
                data = pd.DataFrame(pd.read_json(kwargs['datafile']))
            except:
                raise Exception('could not read datafile. Does it exist ?')
        
        self.data = data

        self.data['selected'] = True # we are using only this subset
        total_flux = np.array([sum(flux) for flux in self.data['flux_values_extract']])
        self.data['total_flux'] = total_flux

        if self.kwargs.get('bootstrap_data', False):
            self.bootstrap_data()

        # basic outlier removal / quality filter
        self._outlier_removal()
        self.fit_global_polynomial()

        self._map_1D_x_ol_o = None
        self._map_1D_ol_x_o = None
        
        self._map_2D_x_ol_o = None
        self._map_2D_ol_x_o = None
        
        
    def generate_fixed_effects(self):
        """
        for later use when replacing the 2d polynomial with 2d splines
        """
        Nlo = self.kwargs['Nlo']
        No = len(self.all_order())
        N = max(Nlo, No)
        ol_range = ((self.o * self.l).min(), (self.o * self.l).max())
        o_range = (self.all_order().min(), self.all_order().max())
        I = np.identity(No + 1)
        Cheb_o = [np.polynomial.chebyshev.Chebyshev(c, window=[-1,1], domain=o_range) for c in I]
        lo_range = (self.all_order().min(), self.all_order().max())
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
        self._bootstraped = True

    def color_of_order(self, o):
        cmap = cm.get_cmap(self.kwargs['palette_order'])
        orders = self.all_order()
        colors = cmap(np.linspace(0,1, max(orders) - min(orders) +1))
        return colors[o-min(orders)]

    def _outlier_removal(self):
        """
        removing outliers based on bootstrap uncertainty
        of histotram expectation
        the removed data is dropped from the .data
        """
        epsilon = self.kwargs['epsilon_sigma_bootstrap']
        self.data = self.data[self.data['sigma_new_mean'] < epsilon].copy().reset_index(drop=True)

    def clear_index(self):
        """
        makes all data selected again
        """
        self.data['selected'] = True

    """
    def eval_polynomial(self, ol, o):
        ol = to_list(ol)
        o = to_list(o)
        return np.array([ self.polynomial_fit[oo](ooll) for oo, ooll in zip(o, ol)])

    def eval_dpolynomial_dl(self, ol, o):
        ol = to_list(ol)
        o = to_list(o)
        return np.array([ oo*self.polynomial_fit[oo].deriv()(ooll) for oo, ooll in zip(o, ol)])
    """    
    
    def sigma_clipping(self):
        # sigma clipping for x=P_o(ol) (i.e. 1D) and x = P(ol, o) (i.e. 2D)

        

        def _get_threshold(epsilon, fit_now):
            d_fit_now = {o: p.deriv(1) for o, p in fit_now.items() }
            dxdl = self._o * np.abs(self._eval_order_by_order_full(d_fit_now, self._ol))

            # how to define the clipping
            if self.kwargs['clipp_method'] == 'pixerror':   # absolute error in pixel
                thd = 1 * epsilon
            elif self.kwargs['clipp_method'] == 'rel_std':  # relative error in std of each pixel
                thd = self._sigma * epsilon
            elif self.kwargs['clipp_method'] == 'vrad':     # radial velocity uncertainty
                thd = self._l * dxdl * epsilon / C_LIGHT
            else:
                raise Exception ('no such method') 
                # thd = epsilon * np.std(res)
            return thd


        # order by order fitting first
        epsilon = self.kwargs['epsilon_sigma_clipp']
        thd = epsilon 
        self.data['selected'] = True

        for i in range(self.kwargs['n_sigma_clipp']):
            fit_now = self._fit_polynomial_order_by_order() # fit on selected data set
            res = self._eval_order_by_order_full(fit_now, self._ol) - self._x # give all points a chance

            # clipping 
            thd = _get_threshold(epsilon, fit_now)
            I = np.abs(res)<=thd

            if np.all(I==self.data['selected']):
                print('stable clipping after {} iterations'.format(i))
                break
            self.data.loc[:, 'selected'] = I[:]

        epsilon=self.kwargs['epsilon_sigma_clipp_2d']

        for i in range(self.kwargs['n_sigma_clipp_2d']):
            fit_now = self._fit_2d_polynomial()
            res = self._eval_order_by_order_full(fit_now, self._ol) - self._x
            thd = _get_threshold(epsilon, fit_now)
            I = np.abs(res)<=thd

            if np.all(I==self.data['selected']):
                print('stable clipping after {} iterations'.format(i))
                break

            self.data.loc[:, 'selected'] = I[:]

        self._map_1D_x_ol_o = self._fit_polynomial_order_by_order()
        self._map_2D_x_ol_o = fit_now
        self._map_1D_ol_x_o = {o: inverse_map(p) for o, p in  self._map_1D_x_ol_o.items()}
        self._map_2D_ol_x_o = {o: inverse_map(p) for o, p in  self._map_2D_x_ol_o.items()}

        self._report_fitresult()

    def sigma_clipping_l_x_o(self):
        # sigma clipping for l=P_o(x)/o (i.e. 1D) and l = P(ol, o)/o (i.e. 2D)
        
        def _get_threshold(epsilon):
            dldx = np.abs(self.eval_dpolynomial_dl(self._o*self._l, self._o))
            # how to define the clipping
            if self.kwargs['clipp_method'] == 'pixerror':   # absolute error in pixel
                thd = 1 * epsilon
            elif self.kwargs['clipp_method'] == 'rel_std':  # relative error in std of each pixel
                thd = self._sigma * epsilon
            elif self.kwargs['clipp_method'] == 'vrad':     # radial velocity uncertainty
                thd = self._l * dxdl * epsilon / C_LIGHT
            else:
                raise Exception ('no such method') 
                # thd = epsilon * np.std(res)
            return thd


        # order by order fitting first
        epsilon = self.kwargs['epsilon_sigma_clipp']
        thd = epsilon 
        self.data['selected'] = True
        o = self.o
        l = self.l
        x = self.x

        for i in range(self.kwargs['n_sigma_clipp']):

            self.fit_polynomial_order_by_order()
            res = self.eval_polynomial(self._o*self._l, self._o) - self._x

           # clipping 
            thd = _get_threshold(epsilon)
            I = np.abs(res)<=thd

            if np.all(I==self.data['selected']):
                print('stable clipping after {} iterations'.format(i))
                break
            self.data.loc[:, 'selected'] = I[:]

        epsilon=self.kwargs['epsilon_sigma_clipp_2d']

        for i in range(self.kwargs['n_sigma_clipp_2d']):
            self.fit_2d_polynomial()
            res = self.eval_polynomial(self._o*self._l, self._o) - self._x
            thd = _get_threshold(epsilon)
            I = np.abs(res)<=thd

            if np.all(I==self.data['selected']):
                print('stable clipping after {} iterations'.format(i))
                break

            self.data.loc[:, 'selected'] = I[:]

        self._report_fitresult()

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

    def _inverse_map_at_o(self, p, o):
        nn = 10001
        # p = self.polynomial_fit[o]
        olams = np.linspace(self._ol.min(), self._ol.max(), nn, endpoint=True)
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
        xx = self.polynomial_fit[o](olams)
        mi, ma = xx.min(), xx.max()
        return interpolate_extend(xx, olams/o, mi, ma)
       

    def _compute_lambda_maps(self, dxolo):
        """
        from a solution x = P(ol, o) compute the corresponding solution
        for l = P(x, o)
        """
        self._lambda_map = {o: self.lambda_map_at_o(p, o) for p, o in listitems(dxolo)} 

    @property
    def lambda_map(self):
        return self._lambda_map


    @property
    def _x(self):
       return self.data['new_mean_pixel_pos']

    @property
    def _l(self):
        return self.data['ref_lambda'] 

    @property
    def _o(self):
        return self.data['true_order_number']

    @property
    def _ol(self):
        return self._o * self._l

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

    def get_orders_in_data(self):
        orders = list(set(self.data['true_order_number'][self.data['selected']]))
        orders = self.kwargs.get('orders', orders)
        return orders

    """
    def get_report_orders(self):
        report_only_nonempty = self.kwargs.get('report_only_nonempty', False)
        if report_only_nonempty:
            return self.get_orders()
        return self.kwargs.get('report_orders', self.get_orders())
    """

    def get_global_polynomial(self, full=False):
        """
        fits a global polynomial on all data
        :full: if True take all data else only the selected data (default False) 
        """
        n = self.kwargs['order_ol']
        ol = self.ol if not full else self._o * self._l
        domain = [ol.min(), ol.max()]
        p = None
        w = self.fit_weight() if not full else np.ones(self.ndata_total)
        x = self.x if not full else self.data['new_mean_pixel_pos']
        try:
            p = np.polynomial.chebyshev.Chebyshev.fit(ol, x,
                                                    domain=domain,
                                                    deg=n,
                                                    w=w)
        except Exception as ex:
            print(ex)
            raise Exception ('problem fitting global polynomial\n')
        return p

    def fit_global_polynomial(self):
        """
        make frequmap a global polynomial
        """
        p = self.get_global_polynomial(full=False)
        self.polynomial_fit = {
            o: p # np.polynomial.Polynomial([0]) 
            for o in self.all_order()
        }
        return p

    @property
    def len_total_data(self):
        return len(self.data['true_order_number'])

    def fit_weight(self):  ## TODO: change name to fit_weight_x_ol_o
        tmp =self.kwargs.get('fitweight', 'sigma') 
        if (tmp == 'equal'):
            weight = 1
        elif (tmp == 'sigma'):
            weight = 1/self.sigma
        elif (tmp == 'vrad'):
            weight = np.ones(self.ndata)
            for o in self.all_order():
                I = self.index_order(o)
                weight[I] = \
                1/(C_LIGHT * np.abs((self.ol[I]*self.eval_dpolynomial_dl(self.ol[I], o))))
        else:
            raise Exception('no such fitweight' +tmp )
        return weight      

    def _eval_order_by_order_full(self, dofmaps, xol):
        """
        evaluates a dictionary of maps at its arguments
        """
        tmp = np.zeros(self.ndata_total)
        for o in self.all_order():
            I = self.index_order_unselected(o)
            tmp[I] = dofmaps[o](xol[I])
        return tmp


    def _fit_polynomial_order_by_order(self):
        """
        x = P_o(ol)
        """
        n = kwargs['order_ol']
        weight = self.fit_weight()
        res = {}
        ol = self.o * self.l
        domain = [ol.min(), ol.max()]
        orders = self.all_order()

        pglobal = self.get_global_polynomial()

        for o in orders:
            I = self.index_order(o)
            try:
                if sum(I)<n+1:
                    raise Exception('not enough points for fit in order')
                res[o] = np.polynomial.chebyshev.Chebyshev.fit(
                    self.l[I]*o, self.x[I],
                    domain=domain,
                    deg=n,
                    w = weight[I]
                )
            except Exception as ex:
                res[o] = pglobal
                print('------------\nAt order {} the following error occured\n-------\n'.format(o),ex)
 
        # self.polynomial_fit = res
        return res

    def _fit_2d_polynomial(self):
        """
        # computes a polynomial proxy
        """
        Nol = self.kwargs['order_ol']+1
        No  = self.kwargs['order_o'] +1

        o = self.o
        l = self.l

        olrange = [(o*l).min(), (o*l).max()]
        orange = [o.min(), o.max()]
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
        s = self.fit_weight()
        s = np.sqrt(s**2 + self.kwargs['sigma_min']**2)

        for i in range(self.ndata):
            for j, (nol, no) in enumerate(nolko):
                G[i,j] = Tshebol[nol](o[i]*l[i]) * Tshebo[no](o[i])

        G = (1./ s)[:, np.newaxis] * G
        coeffs = np.linalg.lstsq(G, (1./s) * x, rcond=None)[0]

        ## compute restriction to each order of 2d polynomial
        res = {
            o : sum([ coeffs[i] * Tshebo[no](o) * Tshebol[nol] 
                for i, (nol, no) in enumerate(nolko)])
                     for o in self.all_order() 
        }
        return res


    def _report_fitresult(self):
        """
        updates the self.data field with the fit results
        """
        self.data["x_1D_ol_o"] = self._eval_order_by_order_full(
            self._map_1D_x_ol_o,
            self._ol)
        
        self.data["x_2D_ol_o"] = self._eval_order_by_order_full(
            self._map_2D_x_ol_o,
            self._ol)

        self.data["l_1D_x_o"] = self._eval_order_by_order_full(
            self._map_1D_ol_x_o,
            self._x
        ) / self._o

        self.data["l_2D_x_o"] = self._eval_order_by_order_full(
            self._map_2D_ol_x_o,
            self._x
        ) / self._o
       
            
    @property
    def ndata(self):
        return len(self.o)

    @property
    def ndata_total(self):
        return len(self._o)

    def ndata_in_order(self, o):
        return sum(self.index_order(o))
    
    def get_lambda_list(self):
        """
        2Dmap for all orders as returned by self.get_report_orders()
        saves to kwargs['file_lambda_list']
        """
        res = np.zeros(len(self.all_order())*NROWS)
        i = 0
        x = np.arange(NROWS)
        for i, o in enumerate(self.all_order()):
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
        return {o: tmp[i*NROWS:(i+1)*NROWS] for i, o in enumerate(self.all_order())}


    def get_rms_per_spectrum(self):
        res = pd.Series(index=self.data.index)
        for o in self.all_order():
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
