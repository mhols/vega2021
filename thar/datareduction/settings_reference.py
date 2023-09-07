"""
global parameters for the datareduction pipeline
"""
import os
import numpy as np
from dotenv import load_dotenv

from units import *

load_dotenv()

comment = {}  ## this comments will be included into the .fits document

SETTING_ID = 'vega_reference'

## ----- directory layout
BASEDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))  # put jour actual base path here
DATADIR = os.path.abspath(os.path.join(BASEDIR, 'vega_reference')); 
REFFILES = os.path.abspath(os.path.join(BASEDIR, 'reffiles')); comment['REFFILES'] = "Referenzfiles"

REFFITSFILE = os.path.abspath(os.path.join(DATADIR, 'NEO_20220903_191404_th0.fits'))  # our refernce thorium
RESULTDIR = os.path.abspath(os.path.join(BASEDIR, 'resultfiles'))

## ----- extractor constants
NROWS = 4208;           comment['NROWS'] = "number of rows"
NCOLS = 4196;           comment['NCOLS'] = "number of columns"
NCROSS = 100;           comment['NCROSS'] = "number of rows/columns in central cross"
REMOVECROSS = int(True)
NROWSBLOCK = 2054   #number of rows in individual blocks
NCOLSBLOCK = 2048   #number of cols in individual blocks
#old moon
##HIGHEXP = 15
##LOWEXP = 4
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

SHIFT_MASK_VOIE1 = list(range(1, VOIE2WIDTH + 2))
SHIFT_MASK_VOIE2 = list(range(-VOIE1WIDTH-1, 0))        # realtive indices of extraction mask
OFFSETRED=16
OFFSETBLUE=16
MEMORY_POSITION = 0.7  # memory of AR1 process for line following
BLAZE_RANGE = list(range(-OFFSETRED, OFFSETBLUE+1))                     # range for the blase function
DEGREEBLAZE = 7                                                   # polynomial degree for blaze funtion fit

CENTRALROW = 2161 # arbitrarily selected central row
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
    57:3086,}
    # 58:3192,} TODO make snipets not to break ...
    # 59:3302}  TODO make robust for empty beams

ORDERS = list(CENTRALPOSITION.keys())

LAMBDAFILE = os.path.join(REFFILES, 'NEXTRA_base_lambda_map.txt')
# LAMBDAFILE = os.path.join(REFFILES, 'HOBO_NEO_20220903_191404_th1_wave1.txt')

# snippets extraction
REF_SPECTRUM = os.path.join(REFFILES, 'thar_spec_MM201006.dat')

#REF_ATLASLINES = os.path.join(REFFILES, 'thar_UVES_MM090311.dat')
REF_ATLASLINES_REEDMAN = os.path.join(REFFILES, 'Redman_table6.dat')
REF_ATLASLINES_UVES = os.path.join(REFFILES, 'thar_UVES_MM090311.dat')
REF_ATLASLINES_CLICKED = os.path.join(REFFILES, 'thar_clicked_uves.csv')
ATLAS_FOR_SNIPPETS = "UVES" #'CLICKED'   # choose from 'UVES', 'REEDMAN', 'CLICKED'
WAVEMAP_IN_VACUUM_AIR = "VACUUM" # 'VACUUM' # or AIR

EXCLUSION = os.path.join(REFFILES, 'excluded.dat')

SEUIL = 0.2 * ADU   # seuil en ADU 
SEUILR = 800.
VRANGE = 6.0 * KM/S      # vrange in km/s
VOIE_METHOD = 'SUM_DIVIDE_CENTRALROW'   # defines flux_123 in .fits

## ------ spectrograph paramter
#voie_method = VOIE_METHOD
n_bootstrap =  0                       # number of bootstrap experiments
profile = 'igauss'                     # fit profile for bootstrap estimate of centroid
loss_function = 'loss_3'               # weighted L2-loss for bootstrap estimate of centroid
CLIPMETHOD = 'threshold'
CLIP_QUANTITY = 'deltavr'
CLIPTHRESHOLD = 300 * M / S
CLIP_MAX_VRAD = 400 * M / S

FITWEIGHT = 'flux'                   # weighted empirical risk
USE_SIGMA_MIN = 'True'                # do not use a minimal sigma in fitting
sigma_min = 0.5 * M / S               # minimial sigma to avoid overfitting

palette_order = 'gist_rainbow'         # palette of orders
order_ol = 5                           # order polynomial in ol
order_o = 7                            # order polynomial in o

kwargs = { k: v for k, v in globals().items() if '_'!=k[0] and
    (type(v) is str or type(v) is float or type(v) is int or type(v) is list or type(v) is dict)
}

## the following parameters are included into the fits files header
PREFIX = 'HOBO_'
HEADER_ITEMS = [k for k in comment.keys() if k[:2]!='__']
