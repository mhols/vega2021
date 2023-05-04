import util
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import regal
import snippets
import spectrograph
from plottextractmixin import PlotExtractMixin
import astropy.io.fits as pyfits
import os
import sys
from scipy.ndimage import minimum_filter1d, generic_filter
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import convolve2d
from numpy.polynomial import Polynomial
import pandas as pd
from astropy.utils.decorators import lazyproperty

### all constants are in settings.py
# SETTINGSMODULE = os.environ.get('SETTINGSMODULE', 'settings')

#try:
#    exec('from %s import *'%(SETTINGSMODULE,))
#except:
#    raise Exception('could not import { SETTINGSMODULE }')
#
#print('SETTINGSMODULE=', SETTINGSMODULE)

##### utility functions


def removeCross(image, **kwargs):

    NCROSS = kwargs['NCROSS']
    NROWSBLOCK = kwargs['NROWSBLOCK']
    NCOLSBLOCK = kwargs['NCOLSBLOCK']

    img1=np.zeros(image.shape)
    img1[:NCROSS,:] = image[NROWSBLOCK:NROWSBLOCK+NCROSS,:]
    img1[NCROSS:NCROSS+NROWSBLOCK,:NCROSS] = image[:NROWSBLOCK,NCOLSBLOCK:NCOLSBLOCK+NCROSS]
    img1[NCROSS+NROWSBLOCK:,:NCROSS] = image[NCROSS+NROWSBLOCK:,NCOLSBLOCK:NCOLSBLOCK+NCROSS]
    img1[NCROSS:NROWSBLOCK+NCROSS,NCROSS:NCOLSBLOCK+NCROSS] = image[:NROWSBLOCK,:NCOLSBLOCK]    #UL
    img1[NCROSS+NROWSBLOCK:,NCROSS:NCOLSBLOCK+NCROSS] = image[NROWSBLOCK+NCROSS:,:NCOLSBLOCK]   #DL
    img1[NCROSS:NROWSBLOCK+NCROSS,NCROSS+NCOLSBLOCK:] = image[:NROWSBLOCK,NCROSS+NCOLSBLOCK:]   #UR
    img1[NCROSS+NROWSBLOCK:,NCROSS+NCOLSBLOCK:] = image[NCROSS+NROWSBLOCK:,NCROSS+NCOLSBLOCK:]  #DR
    return img1

def gauss(x, A, mu, sigma, y_offset):
# A, sigma, mu, y_offset = p
   return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y_offset

def gaussfit(xx,yy):

    indymin=np.argmin(yy)
    indymax=np.argmax(yy)
    indxmin=np.argmin(xx)
    indxmax=np.argmax(xx)
    A=-(yy[indymax]-yy[indymin])
    sigma=0.5*(xx[indxmax]-xx[indxmin])
    mu=xx[indymin]
    y_offset=yy[indymax]
    #            gy=gauss(reflambda[indrefextract],A,mu,sigma,y_offset)
    init_vals = [A, mu,sigma,y_offset]  # for [amp, cen, wid]
    #   print ("init_vals gaussfit:", init_vals)
    gopt, covar = curve_fit(gauss,xx,yy,p0=init_vals)
    #        print(k,gopt[1],gopt[0],gopt[2],gopt[3])
    #Return co-effs for fit and covariance

    """
    plt.plot(xx,gauss(xx,gopt[0],gopt[1],gopt[2],gopt[3]))
    """

    #[gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
    refposi=gopt[1]
    #    print("refposi", refposi)
    return refposi

def secondpoly(xx,yy):

    p=Polynomial.fit(xx,yy,2)
#   model=np.polyval(p,xx)
#   f(x)=ax^2+bx+c, f'(x)= 2ax+b, f'(x)=0 => x=-b/2a
    c=p.convert().coef
    # here coefficients are c[0]*x^2+c[1]*x+c[2]
    refposi= -c[1]/(2*c[2])
    
#   plt.plot(xx,model)
    
    return refposi, c[2]

def load_image_from_fits(fitsfile, **kwargs):
    # REMOVECROSS=kwargs['REMOVECROSS']==int(True)
    a=pyfits.open(fitsfile)
    # print('fitsfile ', fitsfile)
    image = np.clip(a[0].data,-100,65535)
    a.close()
    if has_cross(fitsfile):
        return removeCross(image, **kwargs)
    else:
        return image
        
def header_info_from_fits(fitsfile, keyword):
    a=pyfits.open(fitsfile)
    try:
        key = a[0].header[keyword]
    except:
        raise Exception('keyword not defined')
    a.close()
    return key

def gettimestamp(fitsfile):
        return header_info_from_fits(fitsfile, 'DATE_JUL')

def is_flatfield(fitsfile):
    res = False
    try:
        res = header_info_from_fits(fitsfile, 'OBJECT') == 'Flat'
    except:
        pass
    return res

def is_bias(fitsfile):
    try:
        return header_info_from_fits(fitsfile, 'OBJECT') == 'Bias'
    except:
        return False
    
def is_thorium(fitsfile):
    try:
        return header_info_from_fits(fitsfile, 'OBJECT') == 'Thorium' 
    except:
        return str(fitsfile).endswith('th0.fits')

def is_star(fitsfile, name):
    if fitsfile.endswith('st0.fits'):
        return True
    try:
        res = header_info_from_fits(fitsfile, 'OBJECT') == name
    except:
        pass
    return res

def has_cross(fitsfile):
    return int(header_info_from_fits(fitsfile, 'PORTMODE')) == 0

#liste = listallfits('/Users/boehm/Desktop/extract/Vega_2022TBL')
def listallfits(dirname):
    filelist = []
    for f in os.listdir(dirname):
        if f.endswith('.fits'):
            filelist.append(os.path.join(dirname,f))
    return filelist

def getallflatfits(dirname,texp):
    for f in listallfits(dirname):
        if is_flatfield(f) and header_info_from_fits(f, 'EXPTIME')==texp:
            yield f

def getallbiasfits(dirname):
    for f in listallfits(dirname):
        if is_bias(f):
            yield f

def getallthoriumfits(dirname):
    for f in listallfits(dirname):
        if is_thorium(f) and str(f).endswith('_th0.fits'):
            yield f

def getallstarfits(dirname, name=''):
    for f in listallfits(dirname):
        if is_star(f, name):
            yield f

def meanfits(*fitsfiles, **kwargs):
    su = 0.0
    for f in fitsfiles:
        su += load_image_from_fits(f, **kwargs)
    return su / len(fitsfiles)
    

def masterbias(dirname, **kwargs):

    bia = meanfits(*getallbiasfits(dirname), **kwargs)
    return bia


def masterflat(dirname, **kwargs):
    HIGHEXP = kwargs["HIGHEXP"]
    LOWEXP = kwargs["LOWEXP"]
    CUTORDER = kwargs['CUTORDER']
    NROWS = kwargs['NROWS']

    high = meanfits(*getallflatfits(dirname,HIGHEXP), **kwargs)
    low = meanfits(*getallflatfits(dirname,LOWEXP), **kwargs)

    pblue = beamfit(high, CUTORDER)
    pred =  beamfit(low, CUTORDER-1)

    masklow = np.zeros(high.shape)
    for i in range(NROWS):
        masklow[i,0:int(np.round(0.5*(pblue(i)+pred(i))))]=1
    mf=masklow*low + (1.-masklow)*high

    return mf


def index_along_offset(y, x, delta=[0]):
    x = np.round(x).astype(int)
    delta = np.asarray(delta).astype(int)
    A = y.astype(int)
    B = x[:,np.newaxis] + delta
    return A[:, np.newaxis], B

def mask_along_offset(y, x, delta=[0], **kwargs):
    """
    returns a mask along and offset
    """
    NROWS = kwargs['NROWS']
    NCOLS = kwargs['NCOLS']

    I = np.full((NROWS, NCOLS), False, dtype=bool)
    A, B = index_along_offset(y, x, delta)
    I[A, B] = True
    return I

def extract_along_offset(image, y, x, delta=[0]):
    A,B = index_along_offset(y, x, delta)
    return image[A, B]


def _followorder(image,xstart,ystart,up=True, **kwargs):
    """
    x are column indices
    y are row indices
    orders are "aligned" with columns
    """
   
    delta = kwargs['ABSORPTIONHALFW']
    NROWS = kwargs['NROWS']
    NCROSS = kwargs['NCROSS']
    JUMP = kwargs['JUMP']

    if up:
        rows=range(ystart,NROWS)
        IC = range(0, ystart) # complement
    else:
        rows=range(ystart-1,NCROSS-1,-1)
        IC = range(ystart, NROWS) # complement
    res=np.zeros(NROWS,dtype = float)
    # vals=np.zeros(NROWS)

    x=xstart
    
    for row in rows:
        #extract central position betwwen two beams, and this for each row
        positions=np.arange(int(round(x))-delta, int(round(x))+delta+1)
        values=image[row,positions]
        
        #centralposi, c2 =secondpoly(positions,values)
        #centralposi = x + MEMORY_POSITION * (centralposi-x)

        i = np.argmin(values)
        #ooops = False
        #try:
        #    centralposi, c = secondpoly(positions[i-1:i+2],values[i-1:i+2]) 
        #except:
        #    ooops = True
        centralposi = positions[i]        
        if abs(centralposi - x) <= JUMP:
            res[row] = centralposi
            x = centralposi
        else:
            res[row]=x # only update position but remain for search at same

    # smooth
    res[IC] = 0
    return res
    
def followorder(image,xstart,ystart, **kwargs):
    # SMOOTHWIDTH_BEAM =  kwargs['SMOOTHWIDTH_BEAM']
    res = _followorder(image,xstart,ystart,up=True) + _followorder(image, xstart,ystart,up=False)
    res = generic_filter(res, 
        lambda x: Polynomial.fit(np.arange(len(x)), x, deg=4)(len(x)/2), size=SMOOTHWIDTH_BEAM
    )

    return res
    
def get_lambda(order, orders, **kwargs):
    NROWS=kwargs['NROWS']
    ORDERS = orders
    hobolambda = np.loadtxt(kwargs['LAMBDAFILE'])
    mult = order - min(ORDERS)
    tmp = np.zeros(NROWS)
    lamb = hobolambda[mult*NROWS:(mult+1)*NROWS]
    return lamb

class BeamOrder:
    """
    the beam along an order
    """
    DELTA_BEAM_INTERP = 10
    DEGREE_BEAM_INTERP = 4

    def __init__(self, order, masterflat, centralposition, final=False, **kwargs):
        self._masterflat = masterflat
        self.kwargs = kwargs
        self._order = order
        self.NROWS = kwargs['NROWS']
        self.CENTRALROW = kwargs['CENTRALROW']
        self.NROWS = kwargs['NROWS']
        self.CENTRALPOSITION = centralposition #self.kwargs.get('CENTRALPOSITION')
        self.NCROSS = kwargs['NCROSS']
        self.SHIFT_MASK_VOIE1 = kwargs['SHIFT_MASK_VOIE1']
        self.SHIFT_MASK_VOIE2 = kwargs['SHIFT_MASK_VOIE2']
        self.FLUX_LIMIT = kwargs['FLUX_LIMIT']
        self.BLAZE_RANGE = kwargs['BLAZE_RANGE']
        self.DEGREEBLAZE = kwargs['DEGREEBLAZE']

        self._x = np.arange(self.NROWS)
        self._y = _followorder(masterflat, self.CENTRALPOSITION[order], self.CENTRALROW,up=True, **kwargs) +\
                  _followorder(masterflat, self.CENTRALPOSITION[order], self.CENTRALROW,up=False, **kwargs)

        self._xx = self._x
        self._yy = self._y

        I = np.logical_not(np.logical_or(np.isnan(self._x), np.isnan(self._y)))
        
        ### first approximation
        self._xx = self._x[I]
        self._yy = self._y[I]

        self._lower, self._upper = 0, self.NROWS
        self._evaluator = interp1d(self._xx, self._yy, fill_value='extrapolate') 
        #lambda x: np.interp(x, self._xx, self._yy)
        self._x = np.arange(self.NROWS)
        self._y = self(self._x)

        if not final:
            return


        ### go along with this proxy
        self._lower, self._upper = self._beam_limits()

        try:
            self._evaluator = Polynomial.fit(
                np.arange(self._lower, self._upper),
                self(np.arange(self._lower, self._upper)),            
                deg=4)
        except:
            print('propblem', self._order)
        
        self._y = self(self._x)

    def logging(self, message):
        print('BeamOrder: '+message)
    
    @property
    def x(self):
        return self._x
    
    @property
    def y(self):
        return self._y
    
    @property
    def I(self):
        """
        all indices (usually from NCROSS to NROWS)
        half open interval [lower, upper )
        """
        return np.arange(self._lower, self._upper) 

    @property
    def pixel_range(self):
        return np.array([self._lower, self._upper])

    @property
    def II(self):
        """
        logical mask of I
        """
        tmp = np.full(self.NROWS, False)
        tmp[self.I] = True
        return tmp

    def mask_central_voie12(self):
        x=np.arange(self.NCROSS+1, self.NROWS)
        return mask_along_offset(x, self(x), [0], **self.kwargs)

    def mask_voie1(self):
        x=np.arange(self.NCROSS+1, self.NROWS)
        return mask_along_offset(x, self(x), self.SHIFT_MASK_VOIE1, **self.kwargs)

    def mask_voie2(self):
        x=np.arange(self.NCROSS+1, self.NROWS)
        return mask_along_offset(x, self(x), self.SHIFT_MASK_VOIE2, **self.kwargs)

    def beam_sum_voie1(self, image):
        #addition of 0.5 to take into account true start and end of pixel -0.5:+0.5
        y = 0.5+ self(np.arange(self.NROWS))   ### TODO: check limits
        voie_part_1 = (np.ceil(y)-y) * np.sum(image*self.mask_central_voie12(), axis=1)
        
        return voie_part_1 + np.sum(self.mask_voie1() * image, axis=1)
        

    def beam_sum_voie2(self, image):
        y = 0.5+self(np.arange(self.NROWS))   ### TODO: check limits
        voie_part_2 = (y -np.floor(y)) * np.sum(image*self.mask_central_voie12(), axis=1)

        return voie_part_2 + np.sum(self.mask_voie2() * image, axis=1)
        

    def _beam_limits(self):
        mask = mask_along_offset(self._xx, self._yy, self.BLAZE_RANGE, **self.kwargs)
        # mask = self.mask_blaze()
        tmp = np.sum(mask * self._masterflat, axis=1)

        # TODO tmp = medianfilter(tmp)
        tmp=generic_filter(tmp,lambda x:np.median(x),size=10)
        j = np.argmax(tmp)
        for i in range(j,self.NROWS):
             if tmp[i] < self.FLUX_LIMIT:
                break
        upper = i
        j = np.argmax(tmp)
        for i in range(j,self.NCROSS-1,-1):
             if tmp[i] < self.FLUX_LIMIT:
                break
        lower = i
        if upper==lower:
            # artificial treatment of limits for bad beams
            # TODO introduce good beams
            lower = self.CENTRALROW - 3
            upper = self.CENTRALROW + 3
            self.bad_beam = True
        else:
            self.bad_beam = False
        return lower, upper
               

    def mask_blaze(self):
        return mask_along_offset(np.arange(self.NROWS), self(np.arange(self.NROWS)), self.BLAZE_RANGE, **self.kwargs)

    def beam_sum_blaze(self, image=None):
        if image is None:
            try:
                image = self._extractor.masterflat
            except:
                raise Exception('no extractor in beam')
        return np.sum(self.mask_blaze() * image, axis=1)

    def blaze(self, image):
        """
        returns: polynomial proxy for blaze
        """
        return Polynomial.fit(self.I, self.beam_sum_blaze(image)[self.I], deg=self.DEGREEBLAZE) # TODO check if polynom is good proxy

    @property
    def blaze_along_beam(self):
        return self.blaze()(self.I)

    def __call__(self, x):
        return np.where(np.logical_and(x>=self._lower, x<self._upper), 
            self._evaluator(x), 0)


#def beamfit(image, order, extractor, **kwargs):
#    CENTRALPOSITION = extractor.CENTRALPOSITION
#    CENTRALROW = kwargs['CENTRALROW']
#    NROWS = kwargs['NROWS']

#    y = followorder(image, CENTRALPOSITION[order], CENTRALROW)
#    return BeamOrder(order, np.arange(NROWS), y, extractor=extractor, **kwargs)


#### factory methods for Extractor objects
store = regal.Store()    # we only have one global store

def get_ext(f_thar, level="level_1", **kwargs):
    classmap = {
        'level_1': Extractor_level_1,
        'level_2': Extractor_level_2,
        'level_3': Extractor
    }
    try:
        myext = store.get(f_thar) # retrieval without kwargs is reload of existing
        print('retrieving Extract from ', f_thar)
        return myext
    except Exception as ex:
        print('could not retrieve. ', f_thar, 'Reason: ', ex, 'generating a new one\n----\n')
    try:
        myext = classmap[level](f_thar, **kwargs)
        store.store(f_thar, myext)
        return myext
    except Exception as ex:
        raise Exception('could not create Extract. Reason: ', ex)


class Extractor_level_1:
    """
    A class for the data reduction pipeline

    """

    def __init__(self, fitsfile, **kwargs):
        """
        if fitsfile specified class takes it
        the following keyword arguments are recognized
        imagedata: np.array of image data
        flatfield: np.array of flatfield image
        background: np.array of background image
        dir: path to folder containing all fits files

        fitsfile must be a Thorium!!

        """

        assert is_thorium(fitsfile),\
            'you need to specify a thorium file to generate an extractor'
        self.kwargs = kwargs
        
        self._tharfits = fitsfile        
        self._fitsfile = fitsfile

        self._logindent = 0 
        self.logging('Creation of Level 1 Extract')
        
        self.CENTRALROW = kwargs['CENTRALROW']
        self.NROWS = kwargs['NROWS']

        # all properties can be lasily evaluated TODO: use decorator
        self._masterflat = None
        self._masterbias = None
        self._beams = None
        self._Blaze = None
        self._voie1 = None
        self._voie2 = None
        self._voie2 = None
        self._bare_image = None
        self._image = None 
        self._pix_to_lambda_map_2 = None
        self._pix_to_lambda_map_3 = None
        self._fitsfile_set_for_reduction = False
        self.end_logging()
 
    def update_kwargs(self, **kwargs):
        self.kwargs.update(kwargs)
   

    def logging(self, message):
        """
        write a log string with info
        """
        pref = self._logindent * ' '
        print('\n' + pref + 'Extractor, SETTING_ID: ' + self.SETTINGS_ID + ',\n' + pref + 'ThArg: ', 
              self._tharfits + '\n' + pref +message + '.....')
        self._logindent += 2

    def end_logging(self):

        self._logindent = max(0, self._logindent-2)
        print(self._logindent*' ' + '......done')

    def message(self, message):
        pref = self._logindent * ' '
        print (pref + message + '\n')
    
    
    def get_kwarg(self, key, default):
        if key in self.kwargs:
            return self.kwargs[key]
        else:
            self.kwargs[key] =  default
            return default


    @property
    def SETTINGS_ID(self):
        return self.kwargs.get('SETTING_ID', 'WARNING: no setting id, you should specify one in the settings module')    

    def _restart(self):
        # local quantity
        self._image = None
        self._bare_image = None
        self._voie1 = None
        self._voie2 = None

    @property
    def fitsfile(self):
        """
        convenince method for fitsfile
        """
        if self._fitsfile is None:
            raise Exception('no fitsfile available')
        else:
            return self._fitsfile

    def set_fitsfile(self, fitsfile):
        self._restart()
        self._fitsfile = fitsfile
        self._fitsfile_set_for_reduction = True

    def loadimage(self):
        return load_image_from_fits(self.fitsfile, **self.kwargs)

    @property
    def ORDERS(self):
        return self.kwargs['ORDERS']

    @property
    def PREFIX(self):
        return self.kwargs.get('PREFIX', 'HOBO_') #TODO recompute 

    @property
    def bare_image(self):
        if self._bare_image is None:
            self._bare_image = load_image_from_fits(self.fitsfile, **self.kwargs)
        return self._bare_image

    @property
    def image(self):
        """
        bias removed
        background removed
        """
        if self._image is None:
            image = self.bare_image 
            self._image = image  - self.masterbias 
            self._image -= self._estimate_background(self._image)
        return self._image

    @property
    def masterbias(self):
        if self._masterbias is None:
            try:
                self.logging('computing masterbias')
                self._masterbias = masterbias(self.DATADIR, **self.kwargs)
                self.end_logging()
            except:
                self.end_logging()
                raise Exception('please specify DATADIR or MASTERBIAS')
        return self._masterbias


    @lazyproperty
    def _pix_to_lambda_map_level1(self):
        """
        The preliminary lambda map
        based on a datafile containing the lambdas
        """
        return { o: interp1d(np.arange(self.NROWS), \
                                 get_lambda(o, self.ORDERS, **self.kwargs), 
                                 fill_value='extrapolate') 
                for o in self.ORDERS }
        

    @property
    def DATADIR(self):
        try:
            datadir = self.kwargs['DATADIR']
        except:
            pass
        try:
            datadir = os.path.dirname(self._fitsfile)
        except Exception as ex:
            raise Exception('DIRNAME, reason :', ex)
        return datadir

    @property
    def pix_to_lambda_map_voie1(self):
        """
        is dictionary of mappings fractionary pixel -> wavelength
        """
        return self._pix_to_lambda_map_level1

    @property
    def pix_to_lambda_map_voie2(self):
        #self.message('pix_to_lambda_map_voie1 on level 2')
        return self._pix_to_lambda_map_level1

    @property
    def pix_to_lambda_map_voie3(self):
        return self._pix_to_lambda_map_level1
    
    @property
    def pix_to_lambda_map_voie(self):
        return {
            1: self.pix_to_lambda_map_voie1, 
            2: self.pix_to_lambda_map_voie2, 
            3: None
        }

    @property
    def lambdas_per_order_voie1(self):
        n = np.arange(self.NROWS)
        return {o: self.pix_to_lambda_map_voie1[o](n) for o in self.ORDERS}

    @property
    def lambdas_per_order_voie2(self):
        n = np.arange(self.NROWS)
        return {o: self.pix_to_lambda_map_voie2[o](n) for o in self.ORDERS}

    @property
    def lambdas_per_order_voie3(self):
        n = np.arange(self.NROWS)
        return {o: self.pix_to_lambda_map_voie3[o](n) for o in self.ORDERS}

    @property
    def lambdas_per_order_voie(self):
        return {1: self.lambdas_per_order_voie1, 
                2: self.lambdas_per_order_voie2,
                3: None
                }
    

    def lam_to_o(self, lam):
        """the orders of lambda """
        return [ o for o in self.ORDERS if 
                self.pix_to_lambda_map_voie1[o](0) <= lam and lam <= self.pix_to_lambda_map_1[o](self.NROWS-1)]

    def get_lambda_intens1(self, o):
        I = self.beams[o].I
        return self.lambdas_per_order_voie1[o], self.voie1[o], I

    def get_lambda_intens2(self, o):
        I = self.beams[o].I
        return self.lambdas_per_order_voie2[o], self.voie2[o], I

    def get_lambda_intens3(self, o):
        return None

    def get_lambda_intens(self, voie, o):
        tmp = {
            1: self.get_lambda_intens1,
            2: self.get_lambda_intens2,
            3: None
        }
        return tmp[voie](o)

    def get_lambda_list_voie1(self):
        return np.array([self.lambdas_per_order_voie1[o] for o in self.ORDERS]).ravel()

    def get_lambda_list_voie2(self):
        return np.array([self.lambdas_per_order_voie2[o] for o in self.ORDERS]).ravel()

    def get_lambda_list_voie3(self):
        return np.array([self.lambdas_per_order_voie3[o] for o in self.ORDERS]).ravel()

    def get_lambda_list_voie(self, voie):
        if voie==1:
            return self.get_lambda_list_voie1()
        elif voie==2:
            return self.get_lambda_list_voie2()
        else:
            raise Exception('no such voie')
    
    def lambda_range_voie1(self, o):
        return self.pix_to_lambda_map_voie1[o]([self.I[o][0], self.I[o][-1]])    

    def olambda_range_voie1(self, o):
        l1, l2 = self.lambda_range_voie1(o)
        return o * l1, o * l2

    def _compute_voie1et2(self):
        choice =  self.kwargs['VOIE_METHOD']
        if choice == 'SUM_DIVIDE':
            self._voie1 = {
                    o: self.beams[o].beam_sum_voie1(self.image) \
                    / self.beams[o].beam_sum_voie1(self.masterflat) for o in self.ORDERS
                }
            self._voie2 = {
                o: self.beams[o].beam_sum_voie2(self.image) \
                / self.beams[o].beam_sum_voie2(self.masterflat) for o in self.ORDERS
            }
        elif choice == 'DIVIDE_SUM':
            image = self.image / (self.masterflat + 1)
        
            self._voie1 = {
                o: self.beams[o].beam_sum_voie1(image) for o in self.self.ORDERS
            }

            self._voie2 = {
                o: self.beams[o].beam_sum_voie2(image) for o in self.ORDERS
            }
        elif choice == 'DIVIDE_BLAZE':
            self._voie1 = {
                o: self.beams[o].beam_sum_voie1(self.image) \
                    / self.blaze(o)(np.arange(self.NROWS)) \
                    for o in self.ORDERS
            }
            self._voie2 = {
                o: self.beams[o].beam_sum_voie2(self.image) \
                    / self.blaze(o)(np.arange(self.NROWS)) \
                    for o in self.ORDERS
            }

        elif choice == "SUM_DIVIDE_CENTRALROW":
            self._voie1 = {
                    o: self.beams[o].beam_sum_voie1(self.image) \
                    * self.blaze(o)(self.CENTRALROW) \
                    / self.beams[o].beam_sum_voie1(self.masterflat) for o in self.ORDERS
                }
            self._voie2 = {
                o: self.beams[o].beam_sum_voie2(self.image) \
                * self.blaze(o)(self.CENTRALROW) \
                / self.beams[o].beam_sum_voie2(self.masterflat) for o in self.ORDERS
            }
        else:
            raise Exception('no such method ' + choice)
        return
        

    @property
    def voie1(self):
        if self._voie1 is None:
            self.logging('computing voie1 and voie2')
            self._compute_voie1et2()
            self.end_logging()
        return self._voie1
    
    @property
    def voie2(self):
        if self._voie2 is None:
            self.logging('computing voie1 and voie')
            self._compute_voie1et2()
            self.end_logging()
        return self._voie2

    #---------------------------
    # Continuum
    #----------------------------
    def continuum_voie1(self, o, nnodes, q, qq, qqq):
        v = self.voie1[o]
        l = np.arange(len(v))
        p = util.continuum(l,v,nnodes,q,qq,qqq)
        return p(l) 

    def continuum_voie2(self, o, nnodes, q, qq, qqq):
        v = self.voie2[o]
        l = np.arange(len(v))
        p = util.continuum(l,v,nnodes,q,qq,qqq)
        return p(l) 

    def voie1_all(self, nnodes=10, q=0.3, qq=0.7, qqq=0.95):
        res = []
        for o in self.ORDERS:
            res.extend(self.voie1[o]/self.continuum_voie1(o, nnodes, q, qq, qqq))
        return np.array(res)

    def voie2_all(self, nnodes=10, q=0.3, qq=0.7, qqq=0.95):
        res = []
        for o in self.ORDERS:
            res.extend(self.voie2[o]/self.continuum_voie2(o, nnodes, q, qq, qqq))
        return np.array(res)

    #-------------------------
    # the unmodified voies as sum over the bemas
    #-------------------------
    def _bare_voie1(self, o):
        v = self.beams[o].beam_sum_voie1(self.image)
        v[np.isnan(v)] = 0
        return v

    def _bare_voie2(self, o):
        v = self.beams[o].beam_sum_voie2(self.image)
        v[np.isnan(v)] = 0
        return v

    def _bare_voie3(self, o):
        return [0] # TODO

    @property
    def bare_voie1(self):
        return {o: self._bare_voie1(o) for o in self.ORDERS}  

    @property    
    def bare_voie2(self):
        return {o: self._bare_voie2(o) for o in self.ORDERS}

    @property    
    def bare_voie3(self):
        return {o: self._bare_voie3(o) for o in self.ORDERS}
  
    @property
    def voie(self):
        return {1: self.voie1, 2: self.voie2, 3: None}

    @property
    def bare_voie(self):
        return {1: self.bare_voie1, 2: self.bare_voie2, 3: None}
    
    #-----------------------
    # ranges of 'good' part of voie
    #-----------------------
    @lazyproperty
    def I(self):
        """
        The index set of the good part for each order
        """
        return {o: self.beams[o].I for o in self.ORDERS}

    @lazyproperty
    def Ic(self):
        """
        The logical mask for each order of the 'bad' part (i.e. the complement of good)
        """
        tmp = {}
        for o in self.ORDERS:
            res = np.full(self.NROWS, True)
            res[self.I[o]] = False
            tmp[o] = res
        return tmp

    @lazyproperty
    def Il(self):
        """
        The logical mask for each order of the 'good' part
        """
        tmp = {}
        for o in self.ORDERS:
            res = np.full(self.NROWS, False)
            res[self.I[o]] = True
            tmp[o] = res
        return tmp
    
    @lazyproperty
    def Ibounds(self):
        """
        the pixel ranges of each order
        """
        return {o: (self.I[o][0], self.I[o][-1]) for o, p in self.ORDERS.items()}


    @property
    def n(self):
        return np.arange(self.NROWS)


    #--------------------------------
    # the master flat 
    #--------------------------------

    @property
    def masterflat_high(self):
        return meanfits(*getallflatfits(self.DATADIR, self.kwargs['HIGHEXP']), **self.kwargs)

    @property
    def masterflat_low(self):
        return meanfits(*getallflatfits(self.DATADIR, self.kwargs['LOWEXP']), **self.kwargs)


    def _compute_masterflat(self):
        """
        Bias removed flat
        """
        CUTORDER = self.kwargs['CUTORDER']

        high = self.masterflat_high
        low = self.masterflat_low

        high -= self.masterbias
        low -= self.masterbias

        high -= self._estimate_background(high)
        low -= self._estimate_background(low)

        pb = BeamOrder(CUTORDER, high, self.CENTRALPOSITION, **self.kwargs)    #beamfit(high, CUTORDER)
        pr = BeamOrder(CUTORDER-1, low,  self.CENTRALPOSITION, **self.kwargs)    #beamfit(low, CUTORDER-1)
        
        masklow = np.zeros(high.shape)

        for i in range(self.NROWS):
            I = np.arange( int( round(0.5 * ( pb(i) + pr(i)))))
            masklow[i,I]=1

        return masklow*low + (1.-masklow)*high


    @property
    def CENTRALPOSITION(self):
        """
        extraction of orders and their positions
        """
        try:
            return self.kwargs['CENTRALPOSITION']
        except:
            raise Exception('For a level 1 Extract CENTRALPOSITION must be specified')

    @property
    def masterflat(self):
        if self._masterflat is None:
            self.logging('computing masterflat')
            try:
                DIR = self.DATADIR
            except Exception as ex:
                self.end_logging()
                raise Exception('no masterflat and no DATADIR specified...'+ str(ex))
            self._masterflat = self._compute_masterflat()
            self.end_logging()
        return self._masterflat

    def headerinfo(self,keyword):
        """
        convenince function to extract header info
        """
        return header_info_from_fits(self.fitsfile, keyword)

    def blaze(self,order):
        """
        OFFSETRED: towards the RED side of the prism dispersion, "left" on the image
        OFFSEBLUE: towards the BLUE side of the prism dispersion, "right on the image
                    Careful= OFFSETBLUE can be 2*OFFSETBLUE in case of a 3rd beam!
        """
        return self.beams[order].blaze(self.masterflat)
 

    @property
    def beams(self):
        if self._beams is None:
            self.logging('computing beams')
            self._beams = { o: BeamOrder(o, self.masterflat, self.CENTRALPOSITION, final=True, **self.kwargs) for o in self.ORDERS}
            self.end_logging()
        return self._beams

    # def _compute_flat_voie1et2(self):
    #     self.logging('computing flat_voie1 and flat_voie2')
    #     self._flat_voie1 = {
    #         o: self.beams[o].beam_sum_voie1(self.masterflat) for o in self.ORDERS
    #     }
    #     self._flat_voie2 = {
    #         o: self.beams[o].beam_sum_voie2(self.masterflat) for o in self.ORDERS
    #     }
    #     self.end_logging()

    # @property
    # def flat_voie1(self):
    #     if self._flat_voie1 is None:
    #         self._compute_flat_voie1et2()
    #     return self._flat_voie1


    # @property
    # def flat_voie2(self):
    #     if self._flat_voie2 is None:
    #         self._compute_flat_voie1et2()
    #     return self._flat_voie2

    @property
    def non_normalized_intens_1(self):
        res = []
        for o in self.ORDERS:
            inte = self.bare_voie1[o]
            I =  self.beams[o].II
            inte[np.logical_not(I)] = 0.0
            res.extend(inte)
        return np.array(res)

    @property
    def non_normalized_intens_2(self):
        res = []
        for o in self.ORDERS:
            inte = self.bare_voie2[o]
            I = self.beams[o].II
            inte[np.logical_not(I)] = 0
            res.extend(inte)
        return np.array(res)

    @property
    def non_normalized_intens_3(self):
        return self.non_normalized_intens_2

    
    def get_snippets_voie_order(self, voie, o, lmin, lmax):
        """
        returns a list of snippets in good range [lmin, lmax]
        voie: voie
        o: order
        lmin: iterable of lower bounds
        lmax: iterable of upper bounds
        """
        l, v, I = self.get_lambda_intens(voie, o)
        l = l[I]
        v = v[I]
        return util.extract_snippets(l, v, lmin, lmax)

    
    @property
    def noise_1(self):
        return np.sqrt(self.non_normalized_intens_1)

    @property
    def noise_2(self):
        return np.sqrt(self.non_normalized_intens_2)

    @property
    def noise_3(self):
        return np.sqrt(self.non_normalized_intens_3)


    def _estimate_background(self, image):
        MIN_FILTER_SIZE = 50
        SMOOTHING_FILTER_SIZE = 50
        res = np.zeros(image.shape)
        F = np.ones(SMOOTHING_FILTER_SIZE)/SMOOTHING_FILTER_SIZE
        F = np.convolve(F,F)
        for i, r in enumerate(image):
            d = minimum_filter1d(r, size=MIN_FILTER_SIZE)
            res[i,:] = np.convolve(F,d, 'same')
        return res

    @property
    def background(self):
        return self._estimate_background(self.bare_image)

    def cutoff_mask(self, *args):
        CUTOFF_PERCENTILE = args[0]
        tmp = self.background
        v = np.quantile(tmp.ravel(), CUTOFF_PERCENTILE/100)
        return tmp > v

    #------------------------
    # STORE manipulation
    #---------------------------
    def save_to_store(self):
        store.store(self._tharfits, self)


    #-------------------
    # export results to fits
    #-------------------

    def save_fits(self):
        # collecting data
        order = []
        for o in self.ORDERS:
            for i in range(self.NROWS):
                order.append(o)
      
        mask = []
        for o in self.ORDERS:
            for i in range(self.NROWS):
                if i in self.beams[o].I:
                    mask.append(1)
                else:
                    mask.append(0)
        
        # print(len(order), len(mask))
        fitstable = pyfits.BinTableHDU.from_columns ([
        pyfits.Column(
            name="true_order_number",
            format='I',
            array=order
        ),
        pyfits.Column(
            name='flux_mask',
            format = 'I',
            array = mask
        ),
        pyfits.Column(
            name="wavelength_1",
            format="D",
            array=self.get_lambda_list_voie1()
            ),
        pyfits.Column(
            name="wavelength_2",
            format="D",
            array=self.get_lambda_list_voie2()
        ),
        pyfits.Column(
            name="wavelength_3",
            format="D",
            array=np.zeros(len(self.ORDERS)*self.NROWS)
        ),
        pyfits.Column(
            name='flux_1',
            format='E',
            array=self.voie1_all()
        ),
        pyfits.Column(
            name='flux_2',
            format='E',
            array=self.voie2_all()
        ),
        pyfits.Column(
            name='flux_3',
            format='E',
            array=self.voie2_all() ### TODO voie3
        ),
        pyfits.Column(
            name='noise_1',
            format='E',
            array=self.noise_1
        ),
        pyfits.Column(
            name='noise_2',
            format='E',
            array=self.noise_2
        ),
        pyfits.Column(
            name='noise_3',
            format='E',
            array=self.noise_3
        )])

        # fits.Column(
        #     name="blaze_1",
        #     format='E',
        #     array=myext.blaze_1
        # ),
        # fits.Column(
        #     name="blaze_2",
        #     format='E',
        #     array=myext.blaze_2
        # ),
        # fits.Column(
        #     name="blaze_3",
        #     format='E',
        #     array=myext.blaze_3
        # ),
        # fits.Column(
        #     name="continuum_1_1",
        #     format='E',
        #     array=myext.continuum_1_1
        # ),
        # fits.Column(
        #     name="continuum_1_2",
        #     format='E',
        #     array=myext.continuum_1_2
        # ),
        # fits.Column(
        #     name="continuum_1_3",
        #     format='E',
        #     array=myext.continuum_1_3
        # )

        page1 = pyfits.PrimaryHDU()

        PREFIX = self.PREFIX
        filename = os.path.basename(self.fitsfile)
        postfix = ''
        if is_thorium(self.fitsfile):
            postfix = 'th1.fits'
        else:
            postfix = 'st1.fits'
        filename = PREFIX + filename[:-8]+postfix

        RESULTDIR = self.kwargs['RESULTDIR']
        if not os.path.exists(RESULTDIR):
            os.mkdir(RESULTDIR)

        with pyfits.open(self.fitsfile) as sfi:
            page1.header = sfi[0].header

        page1.header.append((PREFIX+"LEVEL", 1))
        #for key in self.kwargs['HEADER_ITEMS']:
        #    page1.header.append((PREFIX+key, globals()[key]))
        # TODO find a way to compute HEADER_ITEMS ....
        newfits = pyfits.HDUList(
            [ page1, fitstable ]
        )

        newfits.writeto(self.result_path, overwrite=True)

    @property
    def result_path(self):
        PREFIX = self.PREFIX
        filename = os.path.basename(self.fitsfile)
        RESULTDIR = self.kwargs['RESULTDIR']
        filename = PREFIX + filename[:-8]+'st1.fits'

        return os.path.join(RESULTDIR, filename)

    def __del__(self):
       pass


######
#  Level 2 lambdamaps are obrained by matchin the level1 reference
####

class Extractor_level_2(Extractor_level_1):
    """
    On this level, the orders can be found by matching a reference
    Extract. This has to be passed in the kwargs
    'REFFITSFILE': fitsfile of ThAr used to compute orders by matching
    'REFKWARGS': kwargs for the creation of level_1 reffits
    """

    def __init__(self, fitsfile, **kwargs):

        super(Extractor_level_2, self).__init__(fitsfile, **kwargs)
        self.logging('Creation of Level 2 Extract')
        self.end_logging()


    @lazyproperty
    def reference_extract(self):
        try:
            extract = get_ext(self.kwargs['REFFITSFILE'])
            return extract
        except:
            self.logging('Could not retrieve reference extractor. Please create one')
            self.end_logging()

    @lazyproperty
    def ORDERS(self):
        return self.reference_extract.ORDERS

    @lazyproperty
    def CENTRALPOSITION(self):
        """
        extraction of orders and their positions
        """
        tmpext = self.reference_extract
        reflow = tmpext.masterflat_low[tmpext.CENTRALROW,:]
        mylow = self.masterflat_low[self.CENTRALROW,:]
        
        Nr = reflow.shape[0]
        Nm = mylow.shape[0]

        nr = np.arange(Nr)
        nm = np.arange(Nm)

        d = np.argmax(np.correlate(reflow, mylow, 'full'))-len(reflow)+1
        b, a = util.homothetie(nr, reflow, nm, mylow, d-5, d+5, 0.9, 1.1)
        res = {o: int(np.round((i-b)/a)) for o, i in tmpext.CENTRALPOSITION.items()}
        return res    

   
    def _pix_to_lambda_map_level2(self, voie):
        """
        The preliminary lambda map
        Either you specify a file containing the lambdas
        or you try to obtain by homothetic mapping from some reference
        """

        self.logging('using a homothetic mapping to get lambda map in level2')
        REFORDER = 33  ## the order used to estimate the homothetie
        
        # TODO use good range self.I
        rv1 = self.reference_extract.voie[voie][REFORDER]
        rv1 = np.where(np.isnan(rv1), 0, rv1)

        # setting myself on first fitsfile (TODO STORE)
        mv1 = self.voie[voie][REFORDER]
        mv1 = np.where(np.isnan(mv1), 0, mv1)
 
        n = np.arange(self.NROWS)
        d = np.argmax(np.correlate(rv1, mv1, 'full')) - self.NROWS + 1

        # TODO: check gridsearch fro homathetie 
        b, a = util.homothetie(n, rv1, n, mv1, d-1, d+1, 0.99999, 1.00001) 

        tmp = {}

        n = np.arange(self.NROWS)

        for o in self.ORDERS:
            tmp[o] = interp1d(n, self.reference_extract.pix_to_lambda_map_voie[voie][o](n+b),
                        fill_value='extrapolate')
        
        self.end_logging()
        return tmp

    @lazyproperty
    def pix_to_lambda_map_voie1(self):
        self.logging('pix_to_lambda_map_voie1 on level2')
        self.end_logging()
        return self._pix_to_lambda_map_level2(1)
    
    @lazyproperty
    def pix_to_lambda_map_voie2(self):
        self.logging('pix_to_lambda_map_voie2 on level2')
        self.end_logging()
        return self._pix_to_lambda_map_level2(2)
    
    @lazyproperty
    def pix_to_lambda_map_voie3(self):
        return self._pix_to_lambda_map_level2(3)
    
    
#####
#  Final level 
#####
class Extractor(PlotExtractMixin, Extractor_level_2):
    """
    A class for the data reduction pipeline
    """

    def __init__(self, fitsfile, **kwargs):
        Extractor_level_2.__init__(self, fitsfile, **kwargs)
        self.snippets_voie1.prepare_snippets()
        self.snippets_voie2.prepare_snippets()
        """
        for i in range(3):
            self.update()
        self.save_to_store()
        self.message(f'Saved to store, you may retrieve it with \n\
                       get_ext({self.fitsfile}')
        """

    @lazyproperty
    def snippets_voie1(self):
        snip = snippets.Snippets(voie=1, extractor=self)
        return snip

    @lazyproperty
    def snippets_voie2(self):
        snip = snippets.Snippets(voie=2, extractor=self)
        return snip


    def snippets_voie(self):
        return {1: self.snippets_voie1, 2: self.snippets_voie2, 3: None}
    
    @lazyproperty
    def ccd_voie1(self):
        self.kwargs["ORDERS"] = self.ORDERS
        return spectrograph.CCD2d(self, self.snippets_voie1.sn, **self.kwargs)
    
    @lazyproperty
    def ccd_voie2(self):
        self.kwargs["ORDERS"] = self.ORDERS
        return spectrograph.CCD2d(self, self.snippets_voie2.sn,  **self.kwargs)

    @property
    def ccd_voie(self):
        return {1: self.ccd_voie1, 2: self.ccd_voie2, 3: None}

    def update_snippets(self):
        #del self.snippets_voie1
        #del self.snippets_voie2

        self.snippets_voie1.update_snippets()
        self.snippets_voie2.update_snippets()

        self.snippets_voie1.sn
        self.snippets_voie2.sn
    
    def update_ccds(self):
        del self.ccd_voie1
        del self.ccd_voie2

        self.ccd_voie1
        self.ccd_voie2

    def update_lambdamap(self):

        self.pix_to_lambda_map_voie1 = self.ccd_voie1._final_map_l_x_o
        self.pix_to_lambda_map_voie2 = self.ccd_voie2._final_map_l_x_o
         
    def update(self):
        self.update_snippets()
        self.update_ccds()
        self.update_lambdamap()


    def clear_ccds(self):
        del self.ccd_voie1
        del self.ccd_voie2
        del self.pix_to_lambda_map_voie1
        del self.pix_to_lambda_map_voie2
        self.logging('lambda map back to level 2')
        self.end_logging()
                    

    def __del__(self):
        del self.snippets_voie1
        del self.snippets_voie2
        del self.reference_extract

"""
TODO: kwargs fill in globally
TODO: index ranges uniformize to [a, b) i.e. half-open
TODO: interorder-background contains strange negative values and perturbates the continuum
"""


def setup_reference():
    from settings_reference import kwargs
    myext = Extractor_level_1(kwargs['REFFITSFILE'], **kwargs)
    # myext.voie # TODO make in constructor
    myext.save_to_store()

    # make better lambda map
    myext = Extractor(kwargs['REFFITSFILE'], **kwargs)
    myext.update()
    myext.update()
    myext.update()
    myext.save_to_store()



if __name__ == '__main__':
    from settings_reference import kwargs as kwargs_reference
    from settingspegasus import kwargs as kwargs_pegasus
    from settingsmoon import kwargs as kwargs_moon
    from settings_vega import kwargs as kwargs_vega
    from units import *
    fitsfile_reference = os.path.join(kwargs_reference['BASEDIR'], 'vega_reference/NEO_20220903_191404_th0.fits')
    fitsfile_pegasus = os.path.join(kwargs_reference['BASEDIR'], '51Peg_raw/2020_0917/NEO_20200917_173122_th0.fits')
    fitsfile_vega = os.path.join(kwargs_reference['BASEDIR'], 'vega_raw/NEO_20220903_191404_th0.fits')
    fitsfile_moon = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230406_190457_th0.fits')

    #myext = Extractor(fitsfile_moon, **kwargs_moon)

    #myext = get_ext(fitsfile_moon)