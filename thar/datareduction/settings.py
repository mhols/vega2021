"""
global parameters for the datareduction pipeline
"""
import os
import numpy as np
from dotenv import load_dotenv
load_dotenv()

comment = {}  ## this comments will be included into the .fits document


RECOMPUTE_2D_POLYNOMIAL = (os.environ.get('RECOMPUTE_2D_POLYNOMIAL', 'False') == 'True')

STARNAME = 'Vega';      comment["STARNAME"] = "Name of object, used to select the starfiles"  
# directory layout

BASEDIR = os.path.join(os.path.dirname(__file__), '../')  # put jour actual base path here

DATADIR = os.path.abspath(os.path.join(BASEDIR, 'datafiles')); 
REFFILES = os.path.abspath(os.path.join(BASEDIR, 'reffiles')); comment['REFFILES'] = "Referenzfiles"
## TMPFILES = os.path.abspath(os.path.join(BASEDIR, 'tmpfiles'))
## RESFILES = os.path.abspath(os.path.join(BASEDIR, 'mappings'))


## ----- extractor constants
NROWS = 4208;           comment['NROWS'] = "number of rows"
NCOLS = 4196;           comment['NCOLS'] = "number of columns"
NCROSS = 100;           comment['NCROSS'] = "number of rows/columns in central cross"
NROWSBLOCK = 2054   #number of rows in individual blocks
NCOLSBLOCK = 2048   #number of cols in individual blocks
HIGHEXP = 60
LOWEXP = 15
CUTORDER = 35   #means that cutting flats is between 34 and 35
ABSORPTIONHALFW = 6 # central region beteween orders
JUMP = 2.;              comment["JUMP"] = "allowed jump for beam extraction"
SMOOTHWIDTH_BEAM = 101; comment["SMOOTHWIDTH_BEAM"] = "width of local polynomial fit"
BACKGROUNDHW = 5                  
VOIE1WIDTH = 18                 #right of separator
VOIE2WIDTH = 18                 #left (redwards) of separator
VOIE3WIDTH = 16
FLUX_LIMIT = 500                # below, the beam extraction is discarded

SHIFT_MASK_VOIE1 = range(1, VOIE2WIDTH + 2)
SHIFT_MASK_VOIE2 = range(-VOIE1WIDTH-1, 0)        # realtive indices of extraction mask
OFFSETRED=16
OFFSETBLUE=16
MEMORY_POSITION = 0.7  # memory of AR1 process for line following
BLAZE_RANGE = range(-OFFSETRED, OFFSETBLUE+1)                     # range for the blase function
DEGREEBLAZE = 7                                                   # polynomial degree for blaze funtion fit
CLUM = 3e5

## unities 
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
    57:3086,}
    # 58:3192,} TODO make snipets not to break ...
    # 59:3302}  TODO make robust for empty beams

ORDERS = list(CENTRALPOSITION.keys())
LAMBDAFILE = os.path.join(REFFILES, 'artlambda2254correct.dat')  # TODO to be removed....



## ----- snippet constants
REF_SPECTRUM = os.path.join(REFFILES, 'thar_spec_MM201006.dat')
REF_ATLASLINES = os.path.join(REFFILES, 'thar_UVES_MM090311.dat')
EXCLUSION = os.path.join(REFFILES, 'excluded.dat')

SEUIL = 0.2 * ADU   # seuil en ADU 
SEUILR = 800.
VRANGE = 9.0 * KM/S      # vrange in km/s
VOIE_METHOD = 'SUM_DIVIDE_CENTRALROW'

## ----- spectrograph paramter  TODO: uppercase
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

kwargs = {k: v for k, v in globals().items if k[:1]!='_'}

# kwargs = {
#     'datadir': DATADIR,
#     'voie_method': VOIE_METHOD,
#     # ----------------

#     'n_bootstrap': 3,                        # number of bootstrap experiments
#     'profile': 'gauss',                      # fit profile for bootstrap estimate of centroid
#     'loss_function': 'loss_1',               # weighted L2-loss for bootstrap estimate of centroid
#     'epsilon_sigma_bootstrap': 3*PIXEL,      # locations with larger pixel uncertainty are removed
#     'epsilon_sigma_clipp': 200*M/S,          # sigma clip for 1d polynomial
#     'epsilon_sigma_clipp_2d' : 200*M/S,      # sigma clip for 2d polynomial
#     'clipp_method' : 'vrad',                 # 'rel_std' or 'pixerror' or 'est_std'
#     'n_sigma_clipp' :    100*TIMES,                   # maximal number of sigma clips
#     'n_sigma_clipp_2d' : 100*TIMES,                # maximal number of sigma clips
#     'fitweight': 'flux',                     # flux (total flux of line) weight based on pixel uncertainty of snippet centroid 
#     'sigma_min' : 5,                         # minimial sigma to avoid overfitting
#     'palette_order': 'gist_rainbow',         # palette of orders
#     'order_ol': 7,                           # order polynomial in ol
#     'order_o': 5,                            # order polynomial in o
# }

## the following parameters are included into the fits files header
PREFIX = 'HOBO_'
HEADER_ITEMS = [k for k in comment.keys() if k[:2]!='__']
