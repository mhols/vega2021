import matplotlib.pyplot as plt
import numpy as np
import pickle
import astropy.io.fits as pyfits
import os
import sys
from scipy.ndimage import minimum_filter1d, generic_filter
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline, bisplrep, bisplev
from scipy.signal import convolve2d
from numpy.polynomial import Polynomial
import pandas as pd

NROWS = 4208        #number of rows
NCOLS = 4196        #number of columns
NCROSS = 100        #number of rows or columns in central cross
NROWSBLOCK = 2054   #number of rows in individual blocks
NCOLSBLOCK = 2048   #number of cols in individual blocks
HIGHEXP = 60
LOWEXP = 15
CUTORDER = 35   #means that cutting flats is between 34 and 35
ABSORPTIONHALFW = 6 # central region beteween orders
JUMP = 2.     # allowed jump for beam extraction
SMOOTHWIDTH_BEAM =   101 # width of local polynomial fit
BACKGROUNDHW = 5                  
VOIE1WIDTH =   18                 #right of separator
VOIE2WIDTH =   18                 #left (redwards) of separator
VOIE3WIDTH =   16
FLUX_LIMIT =   500                # below, the beam extraction is discarded

SHIFT_MASK_VOIE1 = range(1, VOIE2WIDTH + 2)
SHIFT_MASK_VOIE2 = range(-VOIE1WIDTH-1, 0)        # realtive indices of extraction mask
OFFSETRED=16
OFFSETBLUE=16
MEMORY_POSITION = 0.7  # memory of AR1 process for line following
BLAZE_RANGE = range(-OFFSETRED, OFFSETBLUE+1)                     # range for the blase function
DEGREEBLAZE = 7                                                   # polynomial degree for blaze funtion fit
"""
voie3 must be redetermined with new separation beam along minima between 1 and 3
"""

DIRNAME ='./datafiles'
DATADIR=DIRNAME
LAMBDAFILE = './reffiles/artlambda.dat'

"""
Neo-Narval fits files contain nli ligns and  ncol columns. Red orders are low pixel indices, blue ones correspond to high indices. Low pixel indices correspond to lower wavelength within one given order, wavelength is increasing towards higher column indices.

    image[500,:] is row at rowindex = 500
        here it cuts all the orders
    self.image[:,300] is column at columnindex = 300
        it cuts the image parallel to the orders

    Centralrow = image[Centralrowindex,:] = row cutting orders at blaze max

CENTRALPOSITION: corresponds to the lowest value pixel position between two             beams of a given order, the order number beeing the true                order number. Coordinates are image values of ncol x nli of             the fits files. [trueordernumber, columnposition]

OFFSETLIG: offset in line (this zone contains noise information not used                here)
OFFSETCOL: same for column
"""

CENTRALROW = 2161
CENTRALPOSITION = {
    21: 824,
    22: 866,
    23: 907,
    24: 948,
    25: 990,
    26:1032,
    27:1074,
    28:1117,
    29:1161,
    30:1206,
    31:1252,
    32:1299,
    33:1347,
    34:1398,
    35:1449,
    36:1501,
    37:1555,
    38:1611,
    39:1669,
    40:1728,
    41:1789,
    42:1852,
    43:1918,
    44:1985,
    45:2054,
    46:2126,
    47:2200,
    48:2276,
    49:2355,
    50:2436,
    51:2520,
    52:2606,
    53:2697,
    54:2789,
    55:2885,
    56:2984,
    57:3086,
    58:3192,
    59:3302}

# ORDERS = list(CENTRALPOSITION.keys())
ORDERS = range(21, 58)
REPORT_ORDERS = range(21, 60)

# utility functions


def removeCross(image):

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

def load_image_from_fits(fitsfile):
    a=pyfits.open(fitsfile)
    image = np.clip(a[0].data,-100,65535)
    a.close()
    return removeCross(image)

def header_info_from_fits(fitsfile, keyword):
    a=pyfits.open(fitsfile)
    try:
        key = a[0].header[keyword]
    except:
        raise Exception('keyword not defined')
    a.close()
    return key

def is_flatfield(fitsfile):
    return header_info_from_fits(fitsfile, 'OBJECT') == 'Flat'

def is_bias(fitsfile):
    return header_info_from_fits(fitsfile, 'OBJECT') == 'Bias'

def is_thorium(fitsfile):
    return header_info_from_fits(fitsfile, 'OBJECT') == 'Thorium'

#liste = listallfits('/Users/boehm/Desktop/extract/Vega_2022TBL')
def listallfits(dirname=DIRNAME):
    filelist = []
    for f in os.listdir(dirname):
        if f.endswith('.fits'):
            filelist.append(os.path.join(dirname,f))
    return filelist

def getallflatfits(dirname=DIRNAME,texp=HIGHEXP):
    for f in listallfits(dirname):
        if is_flatfield(f) and header_info_from_fits(f, 'EXPTIME')==texp:
            yield f

def getallbiasfits(dirname=DIRNAME):
    for f in listallfits(dirname):
        if is_bias(f):
            yield f

def getallthoriumfits(dirname=DIRNAME):
    for f in listallfits(dirname):
        if is_thorium(f):
            yield f

def meanfits(*fitsfiles):
    sum = 0.0
    for f in fitsfiles:
        sum += load_image_from_fits(f)
    return sum / len(fitsfiles)
    

def masterbias(dirname=DIRNAME):

    bia = meanfits(*getallbiasfits(dirname))
    return bia


def masterflat(dirname=DIRNAME):
    high = meanfits(*getallflatfits(dirname,HIGHEXP))
    low = meanfits(*getallflatfits(dirname,LOWEXP))

    pblue = beamfit(high, CUTORDER)
    pred =  beamfit(low, CUTORDER-1)

    masklow = np.zeros(high.shape)
    for i in range(NROWS):
        masklow[i,0:int(np.round(0.5*(pblue(i)+pred(i))))]=1
    masterflat=masklow*low + (1.-masklow)*high

    return masterflat


def index_along_offset(y, x, delta=[0]):
    x = np.round(x).astype(int)
    delta = np.asarray(delta).astype(int)
    A = y.astype(int)
    B = x[:,np.newaxis] + delta
    return A[:, np.newaxis], B

def mask_along_offset(y, x, delta=[0]):
    """
    returns a mask along and offset
    """
    I = np.full((NROWS, NCOLS), False, dtype=bool)
    A, B = index_along_offset(y, x, delta)
    I[A, B] = True
    return I

def extract_along_offset(image, y, x, delta=[0]):
    A,B = index_along_offset(y, x, delta)
    return image[A, B]


def _followorder(image,xstart,ystart,up=True):
    """
    x are column indices
    y are row indices
    orders are "aligned" with columns
    """
    delta = ABSORPTIONHALFW
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
    
    
def followorder(image,xstart,ystart):
    res = _followorder(image,xstart,ystart,up=True) + _followorder(image, xstart,ystart,up=False)
    res = generic_filter(res, 
        lambda x: Polynomial.fit(np.arange(len(x)), x, deg=4)(len(x)/2), size=SMOOTHWIDTH_BEAM
    )

    return res
    
def get_lambda(order):
    artlambda = loadtxt(LAMBDAFILE)
    selectedorder = o
    mult = selectedorder - 21
    lamb = artlambda[mult*NROWS:(mult+1)*NROWS]
    return lamb

class BeamOrder:
    """
    the beam along an order
    """
    DELTA_BEAM_INTERP = 10
    DEGREE_BEAM_INTERP = 4

    def __init__(self, order, masterflat, final=False):
        self._masterflat = masterflat
        self._order = order

        self._x = np.arange(NROWS)
        self._y = _followorder(masterflat, CENTRALPOSITION[order], CENTRALROW,up=True) +\
                  _followorder(masterflat, CENTRALPOSITION[order], CENTRALROW,up=False)

        I = np.logical_not(np.logical_or(np.isnan(self._x), np.isnan(self._y)))
        
        ### first approximation
        self._xx = self._x[I]
        self._yy = self._y[I]

        self._lower, self._upper = 0, NROWS
        self._evaluator = lambda x: np.interp(x, self._xx, self._yy)
        self._x = np.arange(NROWS)
        self._y = self(self._x)

        if not final:
            return


        ### go along with this proxy
        self._lower, self._upper = self._beam_limits()

        self._evaluator = Polynomial.fit(
            np.arange(self._lower, self._upper),
            self(np.arange(self._lower, self._upper)),            
            deg=4)
        
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
        all indices (usually from NCROSS to RROWS)
        """
        return np.arange(self._lower, self._upper) 
 
    def mask_central_voie12(self):
        x=np.arange(NCROSS+1, NROWS)
        return mask_along_offset(x, self(x), [0])

    def mask_voie1(self):
        x=np.arange(NCROSS+1, NROWS)
        return mask_along_offset(x, self(x), SHIFT_MASK_VOIE1)

    def mask_voie2(self):
        x=np.arange(NCROSS+1, NROWS)
        return mask_along_offset(x, self(x), SHIFT_MASK_VOIE2)

    def beam_sum_voie1(self, image):
        #addition of 0.5 to take into account true start and end of pixel -0.5:+0.5
        y = 0.5+ self(np.arange(NROWS))   ### TODO: check limits
        voie_part_1 = (np.ceil(y)-y) * np.sum(image*self.mask_central_voie12(), axis=1)
        
        return voie_part_1 + np.sum(self.mask_voie1() * image, axis=1)
        

    def beam_sum_voie2(self, image):
        y = 0.5+self(np.arange(NROWS))   ### TODO: check limits
        voie_part_2 = (y -np.floor(y)) * np.sum(image*self.mask_central_voie12(), axis=1)

        return voie_part_2 + np.sum(self.mask_voie2() * image, axis=1)
        

    def _beam_limits(self):
        mask = self.mask_blaze()
        tmp = np.sum(mask * self._masterflat, axis=1)

        # TODO tmp = medianfilter(tmp)

        j = np.argmax(tmp)
        for i in range(j,NROWS):
             if tmp[i] < FLUX_LIMIT:
                break
        upper = i
        j = np.argmax(tmp)
        for i in range(j,NCROSS-1,-1):
             if tmp[i] < FLUX_LIMIT:
                break
        lower = i
        return lower, upper
               

    def mask_blaze(self):
        return mask_along_offset(np.arange(NROWS), self(np.arange(NROWS)), BLAZE_RANGE)

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
        return Polynomial.fit(self.I, self.beam_sum_blaze(image)[self.I], deg=DEGREEBLAZE) # TODO check if polynom is good proxy

    @property
    def blaze_along_beam(self):
        return self.blaze()(self.I)

    def __call__(self, x):
        return np.where(np.logical_and(x>=self._lower, x<self._upper), 
            self._evaluator(x), 0)


def beamfit(image, order, extractor=None):
    y = followorder(image, CENTRALPOSITION[order], CENTRALROW)
    return BeamOrder(order, np.arange(NROWS), y, extractor=extractor)

class Extractor:
    """
    A class for the data reduction pipeline

    """

    def __init__(self,fitsfile=None, **kwargs):
        """
        if fitsfile specified class takes it
        the following keyword arguments are recognized
        imagedata: np.array of image data
        flatfield: np.array of flatfield image
        background: np.array of background image
        dir: path to folder containing all fits files

        """
        self._restart()
        self.kwargs = kwargs
        self._fitsfile = fitsfile

        # all properties can be lasily evaluated

        # masterflat and masterbias can be computed or passed
        self._masterflat = kwargs.get('masterflat', None)
        self._masterbias = kwargs.get('masterbias', None)
        self._beams = None
        self._order_bounds = None
        self._Blaze = None


    def logging(self, message):
        print('Extractor: ' + message)

    def _restart(self):
        # local quantity
        self._image = None
        self._bare_image = None
        self._voie1 = None
        self._voie2 = None
        self._flat_voie1 = None
        self._flat_voie2 = None
        self._bias_voie1 = None
        self._bias_voie2 = None

    def get_fitsfile(self):
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

    fitsfile = property(get_fitsfile, set_fitsfile)

    def loadimage(self):
        return load_image_from_fits(self.fitsfile)

    @property
    def bare_image(self):
        if self._bare_image is None:
            self._bare_image = load_image_from_fits(self.fitsfile)
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
                self._masterbias = masterbias(self.kwargs['DATADIR'])
            except:
                raise Exception('please specify DATADIR or MASTERBIAS')
        return self._masterbias


    def bias_voie1(self, o):
        return self.beams[o].beam_sum_voie1(self.masterbias) 

    def bias_voie2(self, o):
        return self.beams[o].beam_sum_voie2(self.masterbias) 

    
    def _compute_voie1et2(self):
        choice =  self.kwargs.get('VOIE_METHOD', 'SUM_DIVIDE') 
        if choice == 'SUM_DIVIDE':
            self._voie1 = {
                    o: self.beams[o].beam_sum_voie1(self.image) \
                    / self.beams[o].beam_sum_voie1(self.masterflat) for o in ORDERS
                }
            self._voie2 = {
                o: self.beams[o].beam_sum_voie2(self.image) \
                / self.beams[o].beam_sum_voie2(self.masterflat) for o in ORDERS
            }
        elif choice == 'DIVIDE_SUM':
            image = self.image / (self.masterflat + 1)
        
            self._voie1 = {
                o: self.beams[o].beam_sum_voie1(image) for o in ORDERS
            }

            self._voie2 = {
                o: self.beams[o].beam_sum_voie2(image) for o in ORDERS
            }
        elif choice == 'DIVIDE_BLAZE':
            self._voie1 = {
                o: self.beams[o].beam_sum_voie1(self.image) \
                    / self.blaze(o)(np.arange(NROWS)) \
                    for o in ORDERS
            }
            self._voie2 = {
                o: self.beams[o].beam_sum_voie2(self.image) \
                    / self.blaze(o)(np.arange(NROWS)) \
                    for o in ORDERS
            }

        elif choice == "SUM_DIVIDE_CENTRALROW":
            self._voie1 = {
                    o: self.beams[o].beam_sum_voie1(self.image) \
                    * self.blaze(o)(CENTRALROW) \
                    / self.beams[o].beam_sum_voie1(self.masterflat) for o in ORDERS
                }
            self._voie2 = {
                o: self.beams[o].beam_sum_voie2(self.image) \
                * self.blaze(o)(CENTRALROW) \
                / self.beams[o].beam_sum_voie2(self.masterflat) for o in ORDERS
            }
        else:
            raise Exception('no such method ' + choice)
        return
        
 

    @property
    def voie1(self):
        if self._voie1 is None:
            self.logging('computing voie1 and voie2')
            self._compute_voie1et2()
        return self._voie1

    @property
    def voie2(self):
        if self._voie2 is None:
            self.logging('computing voie1 and voie')
            self._compute_voie1et2()
        return self._voie2

    
    def _compute_masterflat(self, dirname=DIRNAME):
        """
        Bias removed flat
        """
        high = meanfits(*getallflatfits(dirname,HIGHEXP))
        low = meanfits(*getallflatfits(dirname,LOWEXP))

        high -= self.masterbias
        low -= self.masterbias

        high -= self._estimate_background(high)
        low -= self._estimate_background(low)

        pb = BeamOrder(CUTORDER, high)    #beamfit(high, CUTORDER)
        pr = BeamOrder(CUTORDER-1, low)    #beamfit(low, CUTORDER-1)
        
        masklow = np.zeros(high.shape)

        for i in range(NROWS):
            I = np.arange( int( round(0.5 * ( pb(i) + pr(i)))))
            masklow[i,I]=1

        return masklow*low + (1.-masklow)*high

    @property
    def masterflat(self):
        if self._masterflat is None:
            self.logging('computing masterflat')
            try:
                DIR = self.kwargs['DATADIR']
            except:
                raise Exception('no masterflat and no DIRNAME specified...')
            self._masterflat = self._compute_masterflat(DIR)
        return self._masterflat

    def headerinfo(self,keyword):
        """
        convenince function to extract header info
        """
        return header_info_from_fits(self.fitsfile, keyword)

    def extimage(self):
        return self.image[NCROSS:,NCROSS:]

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
            self._beams = { o: BeamOrder(o, self.masterflat, final=True) for o in ORDERS}
        return self._beams


    @property
    def flat_voie1(self):
        if self._flat_voie1 is None:
            self.logging('computing flat_voie1 and flat_voie2')
            self._flat_voie1 = {
                o: self.beams[o].beam_sum_voie1(self.masterflat) for o in ORDERS
            }
            self._flat_voie2 = {
                o: self.beams[o].beam_sum_voie2(self.masterflat) for o in ORDERS
            }
        return self._flat_voie1


    @property
    def flat_voie2(self):
        if self._flat_voie2 is None:
            self._flat_voie1 = {
                o: self.beams[o].beam_sum_voie1(self.masterflat) for o in ORDERS
            }
            self._flat_voie2 = {
                o: self.beams[o].beam_sum_voie2(self.masterflat) for o in ORDERS
            }
        return self._flat_voie2

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

    def save_pixel_res(self, outfile):
        os = []
        pix = []
        v1 = []
        v2 = []

        for o in ORDERS:
            for i in range(NROWS):
                os.append(o)
                pix.append(i)
                v1.append(self.voie1[o][i])
                v2.append(self.voie2[o][i])
    
        tmp = {
            'true_order_number' : os,
            'pixel': pix,
            'voie1': v1,
            'voie2': v2
        }
        try:
            with open(outfile, 'w') as f:
                pd.DataFrame(tmp).to_csv(f)
        except:
            pass

        return tmp


"""
TODO: kwargs fill in globally
TODO: index ranges uniformize to [a, b) i.e. half-open
TODO: interorder-background contains strange negative values and perturbates the continuum
"""
