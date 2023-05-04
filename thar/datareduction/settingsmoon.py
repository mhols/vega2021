"""
global parameters for the datareduction pipeline
"""
import os
import numpy as np
from dotenv import load_dotenv
from units import *
from settings_reference import *

load_dotenv()


SETTING_ID = "moon"


## ----- directory layout
BASEDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))  # put jour actual base path here
#DATADIR = os.path.abspath(os.path.join(BASEDIR, 'datafiles'));
DATADIR = os.path.abspath(os.path.join(BASEDIR, 'lune_raw'));
REFFILES = os.path.abspath(os.path.join(BASEDIR, 'reffiles')); comment['REFFILES'] = "Referenzfiles"

REFFITSFILE = os.path.join(BASEDIR, 'vega_reference/NEO_20220903_191404_th0.fits')


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
n_bootstrap = 0                        # number of bootstrap experiments
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
