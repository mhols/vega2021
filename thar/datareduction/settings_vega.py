"""
global parameters for the datareduction pipeline
"""
import os
import numpy as np
from dotenv import load_dotenv
load_dotenv()

from settings_reference import *
from units import *


SETTING_ID = 'vega'

## ----- directory layout
#BASEDIR = os.path.join(os.path.dirname(__file__), '../')  # put jour actual base path here
#DATADIR = os.path.abspath(os.path.join(BASEDIR, 'datafiles')); 
#REFFILES = os.path.abspath(os.path.join(BASEDIR, 'reffiles')); comment['REFFILES'] = "Referenzfiles"

#HIGHEXP = 60
#LOWEXP = 15
#CUTORDER = 35   #means that cutting flats is between 34 and 35
#ABSORPTIONHALFW = 6 # central region beteween orders
#JUMP = 2.;              comment["JUMP"] = "allowed jump for beam extraction"
#SMOOTHWIDTH_BEAM = 101; comment["SMOOTHWIDTH_BEAM"] = "width of local polynomial fit"
#BACKGROUNDHW = 5                  
#VOIE1WIDTH = 18                 #right of separator
#VOIE2WIDTH = 18                 #left (redwards) of separator
#VOIE3WIDTH = 16
#FLUX_LIMIT = 500                # below, the beam extraction is discarded
#
#SHIFT_MASK_VOIE1 = list(range(1, VOIE2WIDTH + 2))
#SHIFT_MASK_VOIE2 = list(range(-VOIE1WIDTH-1, 0))        # realtive indices of extraction mask
#OFFSETRED=16
#OFFSETBLUE=16
#BLAZE_RANGE = list(range(-OFFSETRED, OFFSETBLUE+1))                     # range for the blase function
#DEGREEBLAZE = 7                                                   # polynomial degree for blaze funtion fit

## ------ snippet constants
REF_SPECTRUM = os.path.join(REFFILES, 'thar_spec_MM201006.dat')
REF_ATLASLINES = os.path.join(REFFILES, 'thar_UVES_MM090311.dat')
EXCLUSION = os.path.join(REFFILES, 'excluded.dat')

SEUIL = 0.2 * ADU   # seuil en ADU 
SEUILR = 800.
VRANGE = 9.0 * KM/S      # vrange in km/s
VOIE_METHOD = 'SUM_DIVIDE_CENTRALROW'   # defines flux_123 in .fits

## ------ spectrograph paramter
#voie_method = VOIE_METHOD
n_bootstrap =  10                       # number of bootstrap experiments
profile = 'gauss'                      # fit profile for bootstrap estimate of centroid
loss_function = 'loss_1'               # weighted L2-loss for bootstrap estimate of centroid
sigma_min = 0.001                          # minimial sigma to avoid overfitting
palette_order = 'gist_rainbow'         # palette of orders
order_ol = 5                           # order polynomial in ol
order_o = 7                            # order polynomial in o

kwargs = { 
    k: v for k, v in globals().items() if 
        '_'!=k[0] and 
        'kwargs' != k and 
        (
        type(v) is str or 
        type(v) is float or 
        type(v) is int or 
        type(v) is list or 
        type(v) is dict 
        )
}

HEADER_ITEMS = [k for k in comment.keys() if k[:2]!='__']
