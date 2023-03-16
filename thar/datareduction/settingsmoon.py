"""
global parameters for the datareduction pipeline
"""
import os
import numpy as np
from dotenv import load_dotenv
load_dotenv()

comment = {}  ## this comments will be included into the .fits document


RECOMPUTE_2D_POLYNOMIAL = (os.environ.get('RECOMPUTE_2D_POLYNOMIAL', 'True') == 'True')

STARNAME = 'Thorium';   comment["STARNAME"] = "Name of object, used to select the starfiles"

## ----- directory layout
BASEDIR = os.path.join(os.path.dirname(__file__), '../')  # put jour actual base path here
#DATADIR = os.path.abspath(os.path.join(BASEDIR, 'datafiles'));
DATADIR = os.path.abspath(os.path.join(BASEDIR, 'lune_raw'));
REFFILES = os.path.abspath(os.path.join(BASEDIR, 'reffiles')); comment['REFFILES'] = "Referenzfiles"

## ----- extractor constants
NROWS = 4208;           comment['NROWS'] = "number of rows"
NCOLS = 4196;           comment['NCOLS'] = "number of columns"
NCROSS = 100;           comment['NCROSS'] = "number of rows/columns in central cross"
NROWSBLOCK = 2054   #number of rows in individual blocks
NCOLSBLOCK = 2048   #number of cols in individual blocks
REMOVECROSS = int(False)
HIGHEXP = 15
LOWEXP = 4
CUTORDER = 35   #means that cutting flats is between 34 and 35
ABSORPTIONHALFW = 6 # central region beteween orders
JUMP = 4.;              comment["JUMP"] = "allowed jump for beam extraction"
SMOOTHWIDTH_BEAM = 101; comment["SMOOTHWIDTH_BEAM"] = "width of local polynomial fit"
BACKGROUNDHW = 5                  
VOIE1WIDTH = 18                 #right of separator
VOIE2WIDTH = 18                 #left (redwards) of separator
VOIE3WIDTH = 16
FLUX_LIMIT = 500                # (usual: 500) below, the beam extraction is discarded

SHIFT_MASK_VOIE1 = list(range(1, VOIE2WIDTH + 2))
SHIFT_MASK_VOIE2 = list(range(-VOIE1WIDTH-1, 0))        # realtive indices of extraction mask
OFFSETRED=16
OFFSETBLUE=16
MEMORY_POSITION = 0.7  # memory of AR1 process for line following
BLAZE_RANGE = list(range(-OFFSETRED, OFFSETBLUE+1))                     # range for the blase function
DEGREEBLAZE = 7                                                   # polynomial degree for blaze funtion fit
CLUM = 3e5

## ----- unities 
M = 1.0         # Meter
S = 1.0         # Second
PIXEL = 1       # Pixel unit
TIMES = 1       # number unit (i.e. 50*TIMES )
HZ = 1/S        # Herz
KILO = 1000     # Kilo
KHZ = KILO * HZ # Kiloherz
KM = KILO * M   # Kilometer
C_LIGHT = 300000 * KM / S    # speed of ligth
ADU = 1.0       # astronical data unit ????

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
CENTRALPOSITION = {  ### TODO: move to reffiles...
    21: 849,
    22: 891,
    23: 932,
    24: 974,
    25:1015,
    26:1057,
    27:1099,
    28:1143,
    29:1187,
    30:1232,
    31:1278,
    32:1325,
    33:1374,
    34:1422,
    35:1475,
    36:1527,
    37:1582,
    38:1638,
    39:1696,
    40:1755,
    41:1816,
    42:1880,
    43:1945,
    44:2012,
    45:2081,
    46:2153,
    47:2227,
    48:2304,
    49:2382,
    50:2464,
    51:2548,
    52:2635,
    53:2725,
    54:2817,
    55:2914,
    56:3013,
    57:3115,}
    # 58:3192,} TODO make snipets not to break ...
    # 59:3302}  TODO make robust for empty beams

ORDERS = list(CENTRALPOSITION.keys())
#LAMBDAFILE = os.path.join(REFFILES, 'artlambda2254correct.dat')  # TODO to be removed....

LAMBDAFILE = os.environ.get('LAMBDAFILE', os.path.join(REFFILES, 'hobo.dat') ) # TODO to be removed....
OFFSET_LAMBDA=int(os.environ.get('OFFSET_LAMBDA', 0))
SAVE_LAMS = os.environ.get('SAVE_LAMS', 'False') == 'True'

## ------ snippet constants
REF_SPECTRUM = os.path.join(REFFILES, 'thar_spec_MM201006.dat')
REF_ATLASLINES = os.path.join(REFFILES, 'thar_UVES_MM090311.dat')
EXCLUSION = os.path.join(REFFILES, 'excluded.dat')

SEUIL = 0.2 * ADU   # seuil en ADU 
SEUILR = 800.

VRANGE = 9.0 * KM/S      # vrange in km/s

VOIE_METHOD = 'SUM_DIVIDE_CENTRALROW'   # defines flux_123 in .fits

## ------ spectrograph paramter
voie_method = VOIE_METHOD
datadir = DATADIR
n_bootstrap = 3                        # number of bootstrap experiments
profile = 'gauss'                      # fit profile for bootstrap estimate of centroid
loss_function = 'loss_1'               # weighted L2-loss for bootstrap estimate of centroid
epsilon_sigma_bootstrap = 3*PIXEL      # locations with larger pixel uncertainty are removed
epsilon_sigma_clipp = 200*M/S          # sigma clip for 1d polynomial
epsilon_sigma_clipp_2d = 200*M/S      # sigma clip for 2d polynomial
clipp_method = 'vrad'                 # 'rel_std' or 'pixerror' or 'est_std'
n_sigma_clipp =  100*TIMES                   # maximal number of sigma clips
n_sigma_clipp_2d = 100*TIMES                # maximal number of sigma clips
fitweight = 'flux'                     # flux (total flux of line) weight based on pixel uncertainty of snippet centroid 
sigma_min = 5                         # minimial sigma to avoid overfitting
palette_order = 'gist_rainbow'         # palette of orders
order_ol = 7                           # order polynomial in ol
order_o = 5                            # order polynomial in o

kwargs = { k: v for k, v in globals().items() if '_'!=k[0] and 
    (type(v) is str or type(v) is float or type(v) is int or type(v) is list or type(v) is dict)
}

## the following parameters are included into the fits files header
PREFIX = 'HOBO_'
HEADER_ITEMS = [k for k in comment.keys() if k[:2]!='__']
