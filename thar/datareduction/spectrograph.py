import json
from util import *
import sys
import os
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
import cvxopt as cx
from cvxopt.blas import dot
from cvxopt.solvers import qp

##  local units
from units import *
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
    """
    mi, ma = p.domain
    x = np.linspace(mi, ma, nn, endpoint=True)
    y = p(x)
    yy = pseudo_inverse(y)
    return interpolate_extend(yy, x)
    """
    mi, ma = p.domain
    x = np.linspace(mi, ma, nn, endpoint=True)
    der = p.deriv(1)
    ## locate connected region where the functions monoton and hence has an inverse.
    I, = np.where(der(x) > 0)
    II, = np.where(I[1:]-I[:-1] > 1)
    if len(II)<=1:
        xmin, xmax = x[I[0]], x[I[-1]]
    else:
        III = np.argmax(II[1:]-II[:-1])
        xmin = x[I[II[III]]]
        xmax = x[I[II[III-1]]]

    x = np.linspace(xmin, xmax, nn, endpoint=True)
    y = p(x)
    mi, ma = y.min(), y.max()
    
    return interpolate_extend(y, x)
       
def derivate(p):
    return {o: q.deriv(1) for o, q in p.items()}


class interpolate_extend:
    
    def __init__(self, x, y):
        self._x, self._y = x, y
        self._mi, self._ma = np.min(x), np.max(x)
        self._interp = interp1d(x,y, fill_value="extrapolate")

    @property    
    def domain(self):
        return np.array([self._mi, self._ma])

    def inverse(self):
        return interpolate_extend(self._y, self._x)

    def __call__(self, x):
        return np.where((x-self._mi)*(self._ma-x)>0, self._interp(x), 0)


class CCD2d:
    """
    class for wavemap computation (i.e. 2D polynomial and friends)
    """
    def __init__(self, data, **kwargs):
        self.kwargs = kwargs  # saving for later use
       
        self._data =  data    # reference data
        self._mdata = None   # matching data
        self.kwargs = kwargs  # saving for later use
        if data is None:
            try:
                data = pd.DataFrame(pd.read_json(kwargs['datafile']))
            except:
                raise Exception('could not read datafile. Does it exist ?')
        
        self._data = data
        self.ORDERS = self.kwargs['ORDERS']

        self._data.loc[:, 'selected'] = True # we are using only this subset
        # self._mdata['selected'] = True # matching data
        total_flux = np.array([sum(flux-np.min(flux)) for flux in self._data['bare_voie']])
        self._data.loc[:, 'total_flux'] = total_flux

        #if self.kwargs.get('bootstrap_data', 'True')=='True':
        # self.bootstrap_data()   # TODO move to snippets or extract

        self._map_1D_x_ol_o = None
        self._map_1D_ol_x_o = None
        self._map_2D_x_ol_o = None
        self._map_2D_ol_x_o = None
        
        # basic outlier removal / quality filter
        print('before outlier removal', self.ndata)
        self._outlier_removal()
        print('after outlier removal left ', self.ndata)
        self.sigma_clipping()
        
    @property
    def NROWS(self):
        return self.kwargs['NROWS']

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

    # def bootstrap_data(self):
    #     """
    #     bootstrap estimate of location uncertainty
    #     """
    #     res = []
    #     new_mean_pixel_pos = []
    #     sigma_new_mean = []
    #     shit = 0
    #     print('\n-------------\nbootstrapping\n')
    #     for i, line in self._data.iterrows():
    #         position = np.nan
    #         sigma = np.nan
    #         try:
    #             delta_p, sigma = bootstrap_estimate_location(
    #                 np.array(line['flux_values_extract']),
    #                 **self.kwargs
    #                 )
    #             position = line['pixels_extract'][0] + delta_p

    #         except Exception as ex:
    #             print("problem at ", i, line, ex)  ## TODO: better logging
    #             shit += 1
    #         new_mean_pixel_pos.append(position)
    #         sigma_new_mean.append(sigma)
    #         progress(i,len(self._data), status='')
    #     self._data['new_mean_pixel_pos'] = pd.Series(new_mean_pixel_pos)
    #     self._data['sigma_new_mean'] = pd.Series(sigma_new_mean)
    #     self._bootstraped = True

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
        self._data = self._data[self._sigma < epsilon].copy().reset_index(drop=True)

    def clear_index(self):
        """
        makes all data selected again
        """
        self._data['selected'] = True
    
    @property
    def ndata(self):
        return len(self.o)

    @property
    def ndata_total(self):
        return len(self._o)

    def ndata_in_order(self, o):
        return sum(self.index_order(o))
  
    @property
    def ndata_per_order(self, full=False):
        tmp = {}
        for o in self.all_order():
            I = self.index_order_unselected(o) if full else self.index_order(o)
            tmp[o] = sum(I)
        return tmp

    @property
    def noverlap(self):
        return len(self.mo)

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
        self._fit_global_polynomial(full=True)  ## just to have a start
        epsilon = self.kwargs['epsilon_sigma_clipp']
        thd = epsilon 
        self._data['selected'] = True   ## start with all data
        fit_now = self._fit_polynomial_order_by_order(weight=1./self._sigma) # fit on selected data set

        for i in range(self.kwargs['n_sigma_clipp']):
            res = self._eval_order_by_order_full(fit_now, self._ol) - self._x # give all points a chance
            res = np.abs(res)
            # res = np.abs(res) / self._sigma

            # clipping 
            #thd = _get_threshold(epsilon, fit_now)
            thd = np.quantile(res,0.8)
            I = np.abs(res)<=thd

            if np.all(I==self._data['selected']):
                print('\nstable 1D clipping after {} iterations'.format(i))
                break
            self._data.loc[:, 'selected'] = I[:]
            fit_now = self._fit_polynomial_order_by_order(weight=self._fit_weight(fit_now)) # fit on selected data set

        epsilon=self.kwargs['epsilon_sigma_clipp_2d']

        fit_now2 = self._fit_2d_polynomial(weight=self._fit_weight(fit_now))
        for i in range(self.kwargs['n_sigma_clipp_2d']):
            
            res = self._eval_order_by_order_full(fit_now2, self._ol) - self._x
            res = np.abs(res)
            # res /= self._sigma
            #thd = _get_threshold(epsilon, fit_now2)
            thd = np.quantile(res,0.8)
            
            I = np.abs(res)<=thd

            if np.all(I==self._data['selected']):
                print('stable 2D clipping after {} iterations'.format(i))
                break

            self._data.loc[:, 'selected'] = I[:]
            fit_now2 = self._fit_2d_polynomial(weight=self._fit_weight(fit_now2))

        self._map_1D_x_ol_o = fit_now #self._fit_polynomial_order_by_order(weight=self._fit_weight(fit_now2)) 
        self._map_2D_x_ol_o = fit_now2
        self._map_1D_ol_x_o = {o: inverse_map(p) for o, p in  self._map_1D_x_ol_o.items()}
        self._map_2D_ol_x_o = {o: inverse_map(p) for o, p in  self._map_2D_x_ol_o.items()}

        n = np.arange(self.NROWS)
        self._final_map_l_x_o = {o: 
                interp1d( n, self._map_2D_ol_x_o[o](n)/o, fill_value='extrapolate')
                  for o in self.ORDERS}

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
        self._data['selected'] = True
        o = self.o
        l = self.l
        x = self.x

        for i in range(self.kwargs['n_sigma_clipp']):

            self.fit_polynomial_order_by_order()
            res = self.eval_polynomial(self._o*self._l, self._o) - self._x

           # clipping 
            thd = _get_threshold(epsilon)
            I = np.abs(res)<=thd

            if np.all(I==self._data['selected']):
                print('stable clipping after {} iterations'.format(i))
                break
            self._data.loc[:, 'selected'] = I[:]

        epsilon=self.kwargs['epsilon_sigma_clipp_2d']

        for i in range(self.kwargs['n_sigma_clipp_2d']):
            self.fit_2d_polynomial()
            res = self.eval_polynomial(self._o*self._l, self._o) - self._x
            thd = _get_threshold(epsilon)
            I = np.abs(res)<=thd

            if np.all(I==self._data['selected']):
                print('stable clipping after {} iterations'.format(i))
                break

            self._data.loc[:, 'selected'] = I[:]

        self._report_fitresult()

    def delta_lam_over_lam_factor(self, l, o):
        """
        dl / l = (dl / dx)  dx / l = ((dx / dl)^{-1} / l) * dx
        at o we have x = p(ol) => dx/dl = p'(ol) * o
        returns (dx / dl)^{-1} * x / l 
        """
        p = self.polynomial_fit[o]
        dp = p.deriv(1)
        factor = (o * dp(o*l))**(-1) / l
        return self.kwargs['C_LIGHT'] * factor

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

    @property
    def _x(self):
       return self._data['pixel_mean']

    @property
    def _l(self):
        return self._data['ref_lambda'] 

    @property
    def _o(self):
        return self._data['true_order_number']

    @property
    def _ol(self):
        return self._o * self._l

    @property
    def _sigma(self):
        return self._data['pixel_std']    

    @property
    def _mx(self):
       return self._mdata['pixel_mean_o']

    @property
    def _mxx(self):
       return self._mdata['pixel_mean_oo']

    @property
    def _ml(self):
        return self._mdata['ref_lambda'] 

    @property
    def _mo(self):
        return self._mdata['true_order_number_o']

    @property
    def _moo(self):
        return self._mdata['true_order_number_oo']

    @property
    def _msigmao(self):
        return self._data['pixel_std_o']    

    @property
    def _msigmaoo(self):
        return self._data['pixel_std_oo']    

    @property
    def x(self):
        return self._data['pixel_mean'][self._data['selected']]
        
    @property
    def sigma(self):
        return self._data['pixel_std'][self._data['selected']]

    @property
    def l(self):
        return self._data['ref_lambda'][self._data['selected']]

    @property
    def o(self):
        return self._data['true_order_number'][self._data['selected']]

    @property
    def ol(self):
        return self.o * self.l

    @property
    def mo(self):
        return self._mo[self._mdata['selected']]

    @property
    def mx(self):
        return self._mx[self._mdata['selected']]

    @property
    def moo(self):
        return self._moo[self._mdata['selected']]

    @property
    def mxx(self):
        return self._mxx[self._mdata['selected']]

    @property
    def total_flux(self):
        return self._data['total_flux'][self._data['selected']]

    @property
    def data(self):
        return self._data[self._data['selected']]

    @property
    def mdata(self):
        return self._mdata[self._mdata['selected']]

    def index_order(self, o):
        return self._data['true_order_number'][self._data['selected']]==o

    def index_order_unselected(self, o):
        return self._data['true_order_number']==o

    def all_order(self):
        return self.ORDERS

    def get_orders_in_data(self):
        orders = list(set(self._data['true_order_number'][self._data['selected']]))
        orders = self.kwargs.get('orders', orders)
        return orders

    def get_global_polynomial(self, full=False):
        """
        fits a global polynomial on all data
        :full: if True take all data else only the selected data (default False) 
        the weight is based on sigma as obtained from the bootstrap
        """
        n = self.kwargs['order_ol']
        ol = self.ol if not full else self._o * self._l
        domain = [ol.min(), ol.max()]
        p = None
        w = 1./self.sigma if not full else 1./self._sigma
        x = self.x if not full else self._x
        try:
            p = np.polynomial.chebyshev.Chebyshev.fit(ol, x,
                                                    domain=domain,
                                                    deg=n,
                                                    w=w)
        except Exception as ex:
            print(ex)
            raise Exception ('problem fitting global polynomial\n')
        return p

    def _fit_global_polynomial(self, full=False):
        """
        make frequmap a global polynomial
        """
        p = self.get_global_polynomial(full=full)
        self.polynomial_fit = {
            o: p 
            for o in self.all_order()
        }

    def _fit_weight(self, p):  ## TODO: change name to fit_weight_x_ol_o
        tmp =self.kwargs.get('fitweight', 'sigma') 
        if (tmp == 'equal'):
            weight = np.ones(self.ndata) 
        elif (tmp == 'sigma'):
            weight = 1/self.sigma
        elif (tmp == 'vrad'):
            weight = np.ones(self.ndata)
            dp = self._eval_order_by_order_full(derivate(p), self._ol)
            oldp = self._ol * dp
            oldp = oldp[self._data['selected']]
            for o in self.all_order():
                I = self.index_order(o)
                weight[I] = C_LIGHT / (np.abs(oldp[I]))
        elif (tmp == 'vradsigma'):
            weight = np.ones(self.ndata)
            dp = self._eval_order_by_order_full(derivate(p), self._ol)
            oldps = self._sigma * self._ol * dp
            oldps = oldps[self._data['selected']]
            for o in self.all_order():
                I = self.index_order(o)
                weight[I] = C_LIGHT / (np.abs(oldps[I]))
        elif (tmp == "flux"):
            weight = 1./np.sqrt(self.data['total_flux'])
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

    def get_polynomial_fit_of_order(self, o, weight=None, full=False):
        I = self.index_order_unselected(o) if full else self.index_order(o)
        if weight is None:
            weight = 1./self._sigma if full else 1./self.sigma
        w = weight[I]
        x = self._x[I]  if full else self.x[I]
        l = self._l[I]  if full else self.l[I]
        n = self.kwargs['order_ol']
        return np.polynomial.chebyshev.Chebyshev.fit(
                    o*l, x,
                    domain= (np.min(o*l)*0.99, np.max(o*l)*1.01),
                    deg=n,
                    w = w
                )



    def _fit_polynomial_order_by_order(self, weight=None):
        """
        x = P_o(ol)
        """
        if weight is None:
            weight = 1./self.sigma
        ENLARGEMENT_FACTOR = 1.01
        n = self.kwargs['order_ol']
        res = {}
        ol = self.o * self.l
        domain = [ol.min()/ENLARGEMENT_FACTOR, ol.max()*ENLARGEMENT_FACTOR]
        orders = self.all_order()

        pglobal = self.get_global_polynomial()

        for o in orders:
            I = self.index_order(o)
            try:
                res[o] = np.polynomial.chebyshev.Chebyshev.fit(
                    self.ol[I], self.x[I],
                    domain=domain,
                    deg=min(n, max(0, len(I)-1)),
                    w = weight[I]
                )
            except Exception as ex:
                res[o] = pglobal
                # print('------------\nAt order {} the following error occured\n-------\n'.format(o),ex)
 
        return res

    def _fit_2d_polynomial(self, weight):
        """
        # computes a polynomial proxy
        """
        ENLARGEMENT_FACTOR = 1.1

        Nol = self.kwargs['order_ol']+1
        No  = self.kwargs['order_o'] +1

        o = self.o.to_numpy()
        l = self.l.to_numpy()
        x = self.x.to_numpy()
        s = weight

        #mo = self.mo.to_numpy()
        #moo = self.moo.to_numpy()
        #mx = self.mx.to_numpy()
        #mxx = self.mxx.to_numpy()        

        olrange = [(o*l).min() / ENLARGEMENT_FACTOR, (o*l).max() * ENLARGEMENT_FACTOR]
        orange = [o.min(), o.max()]

        #mxrange = [-200, self.kwargs['NROWS']+200]


        nolko = [(n, k) for n in range(Nol) for k in range(No)]

        ## preparing basis functions for P(ol,o)
        Tshebol = [np.polynomial.chebyshev.Chebyshev.basis(
                    window=[-1,1],
                    domain=olrange,
                    deg=n) for n in range(Nol)]

        Tshebo = [np.polynomial.chebyshev.Chebyshev.basis(
                    window=[-1,1],
                    domain=orange,
                    deg=n) for n in range(No)]

        G = np.zeros((self.ndata, len(nolko)))

        s = np.sqrt(s**2 + self.kwargs['sigma_min']**2)

        for i in range(self.ndata):
            for j, (nol, no) in enumerate(nolko):
                G[i,j] = Tshebol[nol](o[i]*l[i]) * Tshebo[no](o[i])

        G = (1./ s)[:, np.newaxis] * G
        xs = (1./s) * x
        ### preparing matrix for pixel map
        """
        TshebolF = [np.polynomial.chebyshev.Chebyshev.basis(
                    window=[-1,1],
                    domain=mxrange,
                    deg=n) for n in range(2)] 

        TsheboF = [np.polynomial.chebyshev.Chebyshev.basis(
                    window=[-1,1],
                    domain=orange,
                    deg=n) for n in range(10)]
        
        enumbasis = [ (n, m) for n in range(2) for m in range(10)]        
        basis = [ lambda x, o: TshebolF[n](x) * TsheboF[m](o) for n,m in enumbasis]

        H = np.zeros((self.ndata, len(enumbasis)))
        for j, b in enumerate(basis):
            H[:,j] = -b(x, o)
        
        L = np.zeros((self.noverlap, len(enumbasis)))
        for j, b in enumerate(basis):
            L[:, j] = b(mx, mo) - b(mxx, moo)

        o0 = int(np.median(self.ORDERS)) ## central order
        I0 = self.mdata['true_order_number_o'] == o0

        A = np.zeros((np.count_nonzero(I0), len(enumbasis)))
        for j, b in enumerate(basis):
            A[:,j] = b(mx[I0], o0)
        A = cx.matrix(A)

        Xi = np.row_stack((np.column_stack( [G, -H] ), 
                          np.column_stack( [np.zeros((self.noverlap, len(nolko))), L ])))
        
        #### TODO BELOW
        #sig = np.column_stack((1./s, np.full(len(enumbasis), 1e4)))
        #Xi = sig[:,None] * Xi
        #Xi = np.dot(Xi.T, Xi)
        #Xi = cx.matrix(Xi)

        #nab = len(nolko) + len(enumbasis)
        #coeffs = qt( Xi, cx.matrix(0., (1,nab)),   ## quadratic cost function
        #    cx.matrix(0., (1, nab)), cx.matrix(0, (1,1)), # lienear inequality
        #    A, b)['x']                        ## linear constraint

        """


        #coeffs = np.linalg.lstsq(G, xs, rcond=None)[0]
        #G = np.array((1./ s))[:, np.newaxis] * G
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
        updates the self._data field with the fit results
        """
        self._data["x_1D_ol_o"] = self._eval_order_by_order_full(
            self._map_1D_x_ol_o,
            self._ol)
        
        self._data["x_2D_ol_o"] = self._eval_order_by_order_full(
            self._map_2D_x_ol_o,
            self._ol)

        self._data["l_1D_x_o"] = self._eval_order_by_order_full(
            self._map_1D_ol_x_o,
            self._x
        ) / self._o

        self._data["l_2D_x_o"] = self._eval_order_by_order_full(
            self._map_2D_ol_x_o,
            self._x
        ) / self._o

        self._data["dvrad_1D"] = C_LIGHT * ((self._data['l_1D_x_o']/self._l) -1)/(M/S)
        self._data["dvrad_2D"] = C_LIGHT * ((self._data['l_2D_x_o']/self._l) -1)/(M/S)

    def quality(self):
        totalrms2D=np.sqrt(np.mean(self.data["dvrad_2D"]**2))
        totalrms1D=np.sqrt(np.mean(self.data["dvrad_1D"]**2))
        order_rms2D={}
        order_rms1D={}
        nlines={}
        res = {}
        for o in self.all_order():
            I=self.index_order(o)
            res[o]= (
                    np.sum(I), 
                    np.sqrt(np.mean(self.data[I]["dvrad_2D"]**2)),
                    np.sqrt(np.mean(self.data[I]["dvrad_1D"]**2)) 
            )
            nlines[o]=np.sum(I)
            order_rms2D[o]=np.sqrt(np.mean(self.data[I]["dvrad_2D"]**2))
            order_rms1D[o]=np.sqrt(np.mean(self.data[I]["dvrad_1D"]**2))
        return nlines, totalrms2D, order_rms2D, totalrms1D, order_rms1D
   
    def get_lambda_list(self):
        """
        returns lambda_mapping of 2D polynomial
        """
        res = np.zeros(len(self.all_order())*self.NROWS)
        i = 0
        x = np.arange(self.NROWS)
        for i, o in enumerate(self.all_order()):
            ip = self._map_2D_ol_x_o[o]
            res[i*self.NROWS:(i+1)*self.NROWS] = ip(x)/o 
        return np.array(res)

    def get_lambda_map(self):
        tmp = self.get_lambda_list()
        return {o: tmp[i*self.NROWS:(i+1)*self.NROWS] for i, o in enumerate(self.all_order())}


class FP:
    """
    Fabry-Perrot modeling
    """
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.fpdata = json.load(open(kwargs['fp_file'],'r'))
        self.ccd = kwargs.get('ccd2d', CCD2d(**kwargs))
        self.ccd.sigma_clipping()

        self._data = pd.concat([ pd.DataFrame({
            'true_order_number': int(r['true_order_number']),
            'wave3': pd.Series(r['wave3']),
            'flux3': pd.Series(r['flux3'])
            }) for r in self.fpdata if len(r)>0]
        )

    def flux_at_order(self, o):
        return self._data[self._data['true_order_number']==o]['flux3'].values

    def lam_at_order(self, o):
        p = self.ccd.lambda_map_at_o(o)
        return p(self._data[self._data['true_order_number']==o]['wave3'].values)

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

    """
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
    """
    import snippets
    snip=snippets.Snippets(voie=1, tharfits=os.path.join(DATADIR,'NEO_20220903_191404_th0.fits'))
    data = CCD2d(snip.snippets, **kwargs)
