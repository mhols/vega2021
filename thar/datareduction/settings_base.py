"""
global parameters for the datareduction pipeline
"""
import os
import numpy as np
from dotenv import load_dotenv
load_dotenv()

from units import *

comment = {}  ## this comments will be included into the .fits document

SETTING_ID = ''
PREFIX = 'HOBO_'

## ----- extractor constants
NROWS = 4208;           comment['NROWS'] = "number of rows"
NCOLS = 4196;           comment['NCOLS'] = "number of columns"
NCROSS = 100;           comment['NCROSS'] = "number of rows/columns in central cross"
REMOVECROSS = int(True)
NROWSBLOCK = 2054   #number of rows in individual blocks
NCOLSBLOCK = 2048   #number of cols in individual blocks
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

CENTRALROW = 2161

## ------ snippet constants
def get_kwargs():
    return { k: v for k, v in globals().items() if '_'!=k[0] and
                (type(v) is str or type(v) is float or type(v) is int or type(v) is list or type(v) is dict)
            }

## the following parameters are included into the fits files header
