"""
Global parameters for the datareduction pipeline
Normally you should not change these values in this file directly
Rather create a new setting_xxxx.py module in which you import

from nextra.setting_reference.py import kwargs as kwargs_xxxx

You then edit the keys of this dictionary for you specific project, e.g.

kwargs_xxxx['SETTING_ID'] = 'Moon'

Alternatively you may copy this module as a template and adapt it
"""
import os
import numpy as np
from dotenv import load_dotenv

from nextra.units import *

load_dotenv()

comment = {}  ## this comments will be included into the .fits document

class SettingsReference:
    #:
    SETTING_ID = 'vega_reference'
    """
    The setting ID for logging
    """
    #--------------------------------------------
    #    reference settings for Extractor_level_1
    #--------------------------------------------

    BASEDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
    REFFILES = os.path.abspath(os.path.join(BASEDIR, 'assets/reffiles'))
    REFFITSFILE = os.path.abspath(os.path.join(REFFILES, 'refthar/NEO_20220903_191404_th0.fits'))



    #----------------------------------
    #   geometric constans of ccd
    #---------------------------------
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
    USE_PICKED_CENTRAL_POSITIONS = True    # if True use CENTRALPOSITIONS othterwise try matching 
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

    ORDERS = list(range(21,58))


    #----------------------------------------
    #    wave resolution for Extractor_level_1
    #---------------------------------------
    LAMBDAFILE = os.path.join(REFFILES, 'NEXTRA_base_lambda_map.txt')

#REF_ATLASLINES = os.path.join(REFFILES, 'thar_UVES_MM090311.dat')
REF_ATLASLINES_REEDMAN = os.path.join(REFFILES, 'Redman_table6.dat')
REF_ATLASLINES_UVES = os.path.join(REFFILES, 'thar_UVES_MM090311.dat')
REF_ATLASLINES_CLICKED = os.path.join(REFFILES, 'thar_clicked_uves.csv')
REF_ATLASLINES_IVAN = os.path.join(REFFILES, 'Calibration_Ivan.pkl')
ATLAS_FOR_SNIPPETS = "IVAN" #'CLICKED'   # choose from 'UVES', 'REEDMAN', 'CLICKED'
WAVEMAP_IN_VACUUM_AIR = "AIR" #"VACUUM" # 'VACUUM' # or AIR



    # snippets extraction
    #REF_SPECTRUM = os.path.join(REFFILES, 'thar_spec_MM201006.dat')

    ESTIMATE_BACKGROUND = "BACKGROUND_1D"

    #----------------------------------------
    #  atlas lines available
    #-------------------------------------------
    REF_ATLASLINES_REEDMAN = os.path.join(REFFILES, 'Redman_table6.dat')
    REF_ATLASLINES_UVES = os.path.join(REFFILES, 'thar_UVES_MM090311.dat')
    REF_ATLASLINES_CLICKED = os.path.join(REFFILES, 'thar_clicked_uves.csv')

    # choice of catalog to use
    # TODO: use enum type
    ATLAS_FOR_SNIPPETS = "CLICKED" #'CLICKED'   # choose from 'UVES', 'REEDMAN', 'CLICKED'

    # zones not used (beware of Argon lines)
    EXCLUSION = os.path.join(REFFILES, 'excluded.dat')


    WAVEMAP_IN_VACUUM_AIR = "VACUUM" # 'VACUUM' # or AIR

    SEUIL = 0.2 * ADU   # seuil en ADU
    SEUILR = 800.
    VRANGE = 6.0 * KM/S      # vrange in km/s

    SNIPPETS_PIXEL_DELTA = 2 * PIXEL    # pixel interval around pixel mean for matching catalog
    VOIE_METHOD = 'SUM_DIVIDE_CENTRALROW'   # defines flux_123 in .fits

    # ------ spectrograph paramter

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
    order_ol = 5                           # order polynomial in ol
    order_o = 7                            # order polynomial in o


    #----------------------------
    #   plotting settings
    #----------------------------------
    palette_order = 'gist_rainbow'         # palette of orders




def get_kwargs():
    return { k: v for k, v in SettingsReference.__dict__.items() if not k.startswith('_') }

## the following parameters are included into the fits files header
PREFIX = 'NEXTRA_'
HEADER_ITEMS = [k for k in comment.keys() if k[:2]!='__']
