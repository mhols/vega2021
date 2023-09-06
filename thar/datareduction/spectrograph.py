import json
from util import *
import sys
from collections import UserDict
import os
import pandas as pd
import numpy as np
import scipy as sp
from scipy.optimize import bisect
from scipy.interpolate import interp1d
from numpy.polynomial import Polynomial
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import cvxopt as cx
from cvxopt.blas import dot
from cvxopt.solvers import qp
from astropy.utils.decorators import lazyproperty

import traceback

##  local units
from units import *
##  -------------


def to_list(x):
    if hasattr(x, '__iter__'):
        return x
    else:
        return [x]

def inverse_map(p, domain=None):   
    nn = 10001
    """
    mi, ma = p.domain
    x = np.linspace(mi, ma, nn, endpoint=True)
    y = p(x)
    yy = pseudo_inverse(y)
    return interpolate_extend(yy, x)
    """
    if domain is None:
        domain = p.domain
    mi, ma = domain
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
       
def deriv(p):
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


class DictOfMapsPerOrder(UserDict):
    """
    class of maps per order
    """

    def __init__(self, ccd, data):
        super(DictOfMapsPerOrder, self).__init__(data)
        
        self.ccd = ccd
        # self._dict_of_maps =  dict_of_maps

    @property
    def index(self):
        return self.ccd._data.index

    @lazyproperty
    def deriv(self):
        return DictOfMapsPerOrder(
            self.ccd,
            {o: p.deriv() for o, p in self.items()}
        )    

    # def _inverse_map(self, p, factor = 1.0, I=None, domain=None):   
    #     nn = 100001
    #     """
    #     mi, ma = p.domain
    #     x = np.linspace(mi, ma, nn, endpoint=True)
    #     y = p(x)
    #     yy = pseudo_inverse(y)
    #     return interpolate_extend(yy, x)
    #     """
    #     if p is None:
    #         return None
    #     if domain is None:
    #         mi, ma = p.domain
    #     else:
    #         pass
    #     mi, ma = domain
        
    #     x = np.linspace(mi, ma, nn, endpoint=True)
    #     y = p(x)
    #     n = np.arange(self.ccd.NROWS)
    #     a, b = I
    #     J = (a <= y & y <= b)

    #     try:
    #         c, grow, decr = monotoneous_chunks(y[J])
    #     except:
    #         return None
    #     a, b = c[0]
    #     return interpolate_extend(y[J][a:b+1], factor * x[J][a:b+1])
    #     #return interpolate_extend(y, factor * x)
    
    def inverse_map(self):   
        return DictOfMapsPerOrder(
            self.ccd,
            { o: p.inverse() for o, p in self.items()},
        )

    def inverse_o_map(self): 
        res = {}
        for o, p in self.items():
            q = p.f.copy()
            q.domain = q.domain / o
            qq = MonotoneFunction(p.a/o, p.b/o, q, q.deriv(1))
            res[o] = qq.inverse()
        return DictOfMapsPerOrder(self.ccd, res)


    def __call__(self, xol, o):
        """
        evaluates a dictionary of maps at its arguments
        xol is dataframe or series
        its index is alligned with the index of extract._data
        This allows to determine, which map to apply to which part
        """
        tmp = pd.Series(np.nan, index=xol.index)
        # np.zeros(self.ndata_total)

        for oo in self.ccd.all_order():
            I = (o==oo)
            try: 
                tmp.loc[I] = self[oo](xol.loc[I].values)
            except:
                tmp.loc[I] = np.nan
        return tmp

class CCD2d:
    """
    class for wavemap computation (i.e. 2D polynomial and friends)
    """
    def __init__(self, extractor, data, **kwargs):
        self.kwargs = kwargs        # saving for later use
        self.extractor = extractor  # back reference to central extractor objects
        self._data =  data.copy()   # data coming from snippets
        self._data.loc[self._data.index, 'selected'] = True # we are using only this subset
        self._mdata = None          # matching data
        self.ORDERS = self.kwargs['ORDERS']     # convenince field
        self.NROWS = self.kwargs['NROWS']

        total_flux = np.array([sum(flux-np.min(flux)) for flux in self._data['bare_voie']])
        self._data.loc[self._data.index, 'total_flux'] = 1.0* total_flux
        
        self._map_1D_x_ol_o = None
        self._map_1D_ol_x_o = None
        self._map_2D_x_ol_o = None
        self._map_2D_ol_x_o = None
        self._final_map_l_x_o = None

       
        # all fits are done in this range for ol
        self._olmin = self._ol.min()
        self._olmax = self._ol.max()


        p = self.get_global_polynomial()
        self._global_polynomial =  DictOfMapsPerOrder(self, {
                o: MonotoneFunction( 
                    self._olmin, self._olmax,
                    p, p.deriv(1))  
                for o in self.ORDERS
            })
        
        self._set_global_polynomial()

        self._final_map_x_ol_o = self._global_polynomial # shall be updated later

        if kwargs.get('USE_1D_POLYNOMIAL', 'False') == 'True':
            self.sigma_clipping_polynomial_order_by_order()
        else:
            self.sigma_clipping_2D_polynomial()
        
    def color_of_order(self, o):
        cmap = cm.get_cmap(self.kwargs['palette_order'])
        orders = self.all_order()
        colors = cmap(np.linspace(0,1, max(orders) - min(orders) +1))
        return colors[o-min(orders)]


    def _set_global_polynomial(self):
        """
        set a global polynomial order by order 
        """

        # # as default start polynome use a global polynomial
        p = self.get_global_polynomial()

        self._global_polynomial =  DictOfMapsPerOrder(self, {
                o: MonotoneFunction( 
                    self._olmin, self._olmax,
                    p, p.deriv(1))  
                for o in self.ORDERS
            })

    #===================================    
    # clipping procedures
    #===================================
    def _clippings(self, p):
        """
        Actual approximation is p.
        the residuals are computed by self._mismatch(p) 
        they are then used to  define 
        valid sets of retained points. The following methods are available and
        can be specified with 'CLIPMETHOD'. The corresponding threshold is in 
        'CLIPSHRESHOLD'

        'maxvrad':  the error is measured in vrad
        'quantile': das alpha quantile wird behalten (can be set with CLIPSHRESHOLD) 
        'quantile_order_by_order': the same order by order

        p: map(ol, l) -> x

        returns: boolean array of data rows that fullfill accepted criteria

        """
        clip = self.kwargs.get('CLIPMETHOD', 'quantile')
        alpha = self.kwargs.get('CLIPTHRESHOLD', 0.7)
        maxvrad = self.kwargs.get('CLIP_MAX_VRAD', 1000 * M/S)
        clip_quantity = self.kwargs.get('CLIP_QUANTITY', 'deltavr')
        
        # compute clip_quantiy for the approximating mapping p
        r = self._mismatches(clip_quantity, p)
        absr = r.abs()

        # global clipping
        def clip_max_vrad():
            deltavr = self._mismatches('deltavr', p)
            return deltavr.abs() < maxvrad
        I = clip_max_vrad()
        
        def clip_quantile():
            return np.logical_and(I, absr < np.quantile(absr, alpha))

        def clip_quantile_order_by_order():
            res = pd.Series(False, self._data.index)
            for o in self.ORDERS:
                II = self.index_order_unselected(o)
                try:
                    res.loc[II] = np.logical_and(I, absr.loc[II] < np.quantile(absr.loc[II], alpha))
                except:
                    print('in order ', o, 'fall back to default')
                    res.loc[II] = False
            return res

        def clip_est_std():
            return np.logical_and(I, absr < alpha * r.std())

        def clip_threshold():
            return np.logical_and(I, absr < alpha)
        
        def clip_noclip():
            return absr > -1.0 ## always true.....

        # the collection of our possible clipping functions
        # this avoids a length if statement
        clippings = {
            'quantile': clip_quantile,
            'quantile_order_by_order': clip_quantile_order_by_order,
            'noclip': clip_noclip,
            'max_vrad': clip_max_vrad,
            'threshold': clip_threshold,
            'est_std': clip_est_std,
        }

        # evaluate the requested clipping function
        return clippings[clip]()

    # mismaches
    def _mismatches(self, mismatch, p):
        """
        The wavemap p us used to compute the residuums measured with different weights

        return: DataFrame of mismatches
        """

        # polynomial approximation of x
        px = p(self._ol, self._o)
        r = self._x - px # the 'pixel residuum'

        mismatches = {
            'deltax':   lambda: r,
            'deltal':   lambda: r * self._dldx, 
            'deltavr':  lambda: r * C_LIGHT * self._dldx / self._l,
            'chi2':     lambda: r / self._sigma,
            'fitweightchi2': lambda: r / self._fit_weight()
        }

        # evaluate lambda expression in local context
        return mismatches[mismatch]()


    def _fit_weight(self): 
        """
        Returns the weight for the fit in x.
        """
        
        tmp =self.kwargs.get('FITWEIGHT', 'sigma') 


        if (tmp == 'equal'):
            weight = pd.Series(1.0, index=self._data.index) 
        
        elif (tmp == 'sigma'):
            weight = 1/self._sigma
        
        elif (tmp == 'vrad'):
            weight = C_LIGHT * self._dldx / self._l

        elif (tmp == 'intens_weighted_vrad'):
            weight = C_LIGHT * self._dldx * np.sqrt(self._data['total_flux']) / self._l

        elif (tmp == 'sigma_weighted_vrad'):
            weight = C_LIGHT * self._dldx / (self._sigma * self._l)

        elif (tmp == "flux"):
            weight = np.sqrt(self._data['total_flux'])

        else:
            raise Exception('no such fitweight' +tmp )

        return weight      

    def sigma_clipping_polynomial_order_by_order(self):

        w = self._fit_weight()

        fitmachine = lambda I: self._fit_polynomial_order_by_order(
            weight=w[I], I=I
        )

        clipmachine = self._clippings     

        p, I, nclip = sigma_clipping_general_map(
            fitmachine,
            clipmachine,
            I0 = self._data.index   # start from all points
        )

        self._data.loc[:, 'selected_order_by_order'] = False
        self._data.loc[I, 'selected_order_by_order'] = True

        print('stable clipping after ', nclip, 'operations')
         

        self._map_1D_x_ol_o = p
        self._map_1D_ol_x_o = p.inverse_map()
        self._final_map_x_ol_o = p
        self._final_map_l_x_o = p.inverse_o_map()

        self._data["x_1D_ol_o"] = self._map_1D_x_ol_o(self._ol, self._o)
        self._data["l_1D_x_o"] = self._map_1D_ol_x_o(self._x, self._o) / self._o
        self._data["dvrad_1D"] = C_LIGHT * (self._data['x_1D_ol_o'] - self._x) * self._dldx /self._l/(M/S)

        return p

    def sigma_clipping_2D_polynomial(self):
        w = self._fit_weight()

        fitmachine = lambda I: self._fit_2d_polynomial(
            weight=w[I], I=I
        )

        clipmachine = self._clippings

        p, I, nclip = sigma_clipping_general_map(
            fitmachine,
            clipmachine,
            I0 = self._data.index
        )


        self._data.loc[:, 'selected'] = False
        self._data.loc[I, 'selected'] = True

        print('stable clipping after ', nclip, 'operations')
         
        self._map_2D_x_ol_o = p
        self._map_2D_ol_x_o = p.inverse_map()
        self._final_map_x_ol_o = p
        self._final_map_l_x_o = p.inverse_o_map()
        self._data["x_2D_ol_o"] = self._map_2D_x_ol_o(self._ol, self._o)
        self._data["l_2D_x_o"] = self._map_2D_ol_x_o(self._x, self._o) / self._o
        self._data["dvrad_2D"] = C_LIGHT * (self._data['x_2D_ol_o'] - self._x) * self._dldx /self._l/(M/S)
    
    @property
    def ndata(self):
        return len(self.o)

    @property
    def _ndata(self):
        return len(self._o)

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
    def dll(self):
        lmaps = self._final_map_l_x_o
        return 1 - lmaps(self.x, self.o)/self.l

    @property
    def _dll(self):
        lmaps = self._final_map_l_x_o
        return 1 - lmaps(self._x, self._o)/self._l
    
    @property
    def dxdl(self):
        return self.o * self._global_polynomial.deriv(self.ol, self.o)

    @property
    def dldx(self):
        return 1./ self.dxdl

    @property
    def _dxdl(self):
        return self._o * self._global_polynomial.deriv(self._ol, self._o)
    
    @property
    def _dldx(self):
        return 1. / self._dxdl 
   
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

    @property
    def selected(self):
        return self._data['selected']
    
    def index_order(self, o):
        return self._data['true_order_number'][self._data['selected']]==o

    def index_order_unselected(self, o):
        """
        boolean arry true on order o
        """
        return self._data['true_order_number']==o

    def all_order(self):
        return self.ORDERS

    #def get_orders_in_data(self):
    #    orders = list(set(self._data['true_order_number'][self._data['selected']]))
    #    orders = self.kwargs.get('orders', orders)
    #    return orders

    def get_global_polynomial(self, full=False):
        """
        fits a global polynomial on all data
        :full: if True take all data else only the selected data (default False) 
        the weight is based on sigma as obtained from the bootstrap
        """
        n = self.kwargs['order_ol']
        ol = self.ol if not full else self._o * self._l
        domain = [self._olmin, self._olmax]
        p = None
        #w = 1./self.sigma if not full else 1./self._sigma
        x = self.x if not full else self._x
        try:
            p = np.polynomial.chebyshev.Chebyshev.fit(ol, x,
                                                    domain=domain,
                                                    deg=n)
        except Exception as ex:
            print(ex)
            raise Exception ('problem fitting global polynomial\n')
        return p

    def _fit_polynomial_order_by_order(self, weight=None, I = None, deg=None):
        """
        x = P_o(ol)
        """
        if I is None:
            I = self._data.index

        if weight is None:
            weight = 1./self._sigma.loc[I]

        if deg is None:
            deg = self.kwargs['order_ol']

        ENLARGEMENT_FACTOR = 0.01

        # order of polyunomial

        res = {}
        oo = self._o.loc[I]
        ol = (self._ol).loc[I]
        x = self._x.loc[I]

        ols = np.linspace(self._ol.min(), self._ol.max(), 10001)

        domain = [ol.min()*(1-ENLARGEMENT_FACTOR), 
                  ol.max()*(1+ENLARGEMENT_FACTOR)]


        pglobal = self.get_global_polynomial()

        for o in self.ORDERS:
            Io = (oo==o)  # boolean series True on order o
            try:
                p = np.polynomial.Polynomial.fit(
                    ol[Io], x[Io],
                    domain=domain,
                    deg=min(deg, max(0, sum(Io)-1)),
                    w = weight[Io]
                )
                res[o] = MonotoneFunction( *self.extractor.olambda_range_voie1(o), p, p.deriv(1))
            except Exception as ex:
                print(traceback.format_exc())
                print('could not fit polynomial in order ', o, ex)
                res[o] = MonotoneFunction( *self.extractor.olambda_range_voie1(o), pglobal, pglobal.deriv(1))
        
        return DictOfMapsPerOrder(self, res)
       
    def _get_basis(self):
        Nol = self.kwargs['order_ol']+1
        No  = self.kwargs['order_o'] +1



    def _fit_2d_polynomial(self, weight=None, I=None):
        """
        # computes a polynomial proxy
        """
        if I is None:
            I = self._data.index
        if weight is None:
            weight = 1./ self._sigma[I]

        ENLARGEMENT_FACTOR = 1

        Nol = self.kwargs['order_ol']+1
        No  = self.kwargs['order_o'] +1

        o = self._o[I].to_numpy()
        l = self._l[I].to_numpy()
        x = self._x[I].to_numpy()
        s = weight.to_numpy()
        xg = self._global_polynomial(self._ol[I], self._o[I]).to_numpy()

        #mo = self.mo.to_numpy()
        #moo = self.moo.to_numpy()
        #mx = self.mx.to_numpy()
        #mxx = self.mxx.to_numpy()        

        olrange = [self._olmin, self._olmax]
        orange = [self._o.min(), self._o.max()]

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

        G = np.zeros((len(o), len(nolko)))

        # TODO: minimal sigma...
        # if self.kwargs.get('USE_SIGMA_MIN', 'False') == 'True':
        #    s = np.sqrt(s**2 + self.kwargs['sigma_min']**2)

        for i in range(len(o)):
            for j, (nol, no) in enumerate(nolko):
                G[i,j] = Tshebol[nol](o[i]*l[i]) * Tshebo[no](o[i])

        G = (1./ s)[:, np.newaxis] * G
        self._G = G

        xs = (1./s) * (x - xg) # remove ....
        ### preparing matrix for pixel map
        coeffs = np.linalg.lstsq(G, xs, rcond=None)[0] 

        ## compute restriction to each order of 2d polynomial
        res = {
            o : sum([ coeffs[i] * Tshebo[no](o) * Tshebol[nol]
                for i, (nol, no) in enumerate(nolko)])
                     for o in self.all_order() 
        }
        # restore
        res = { o: p + self._global_polynomial[o].f for o, p in res.items()}
        
        res = {o: MonotoneFunction(*self.extractor.olambda_range_voie1(o), p, p.deriv(1)) for o, p in res.items()}

        return DictOfMapsPerOrder(self, res)

    def quality(self):
        J = self._data['selected_order_by_order']
        totalrms2D=np.sqrt(np.mean(self.data["dvrad_2D"]**2))
        totalrms1D=np.sqrt(np.mean(self._data.loc[J, "dvrad_1D"]**2))
        order_rms2D=[]
        order_rms1D=[]
        nlines=[]
        nlines_1D=[]
        totlines = []
        for o in self.all_order():
            I=self.index_order(o)
            JJ = self.index_order_unselected(o)
            nlines.append(np.sum(I))
            totlines.append(np.sum(JJ))
            nlines_1D.append(np.sum(JJ*J))
            order_rms2D.append(np.sqrt(np.mean(self.data.loc[I, "dvrad_2D"]**2)))
            order_rms1D.append(np.sqrt(np.mean(self._data.loc[J*JJ, "dvrad_1D"]**2)))
        return pd.DataFrame({'true_order_number': self.all_order(), 
                             'totlines': totlines,
                             'nlines_1D':  nlines_1D,
                             'nlines_2D': nlines, 
                             'rms_1D': order_rms1D,
                             'rms_2D': order_rms2D}),\
                            totalrms1D, totalrms2D
   

class FP:
    """
    Fabry-Perrot modeling
    """
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.fpdata = json.load(open(kwargs['fp_file'],'r'))
        self.ccd = kwargs.get('ccd2d', CCD2d(**kwargs))
        # self.ccd.sigma_clipping()

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
