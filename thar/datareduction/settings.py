"""
global parameters for the datareduction pipeline
"""
import os
import numpy as np

# directory layout

BASEDIR = os.path.join(os.path.dirname(__file__), '../')  # put jour actual base path here

DATADIR = os.path.abspath(os.path.join(BASEDIR, 'datafiles'))
REFFILES = os.path.abspath(os.path.join(BASEDIR, 'reffiles'))
TMPFILES = os.path.abspath(os.path.join(BASEDIR, 'tmpfiles'))
RESFILES = os.path.abspath(os.path.join(BASEDIR, 'mappings'))

# processing parameters
kwargs = {
    'report_orders': np.arange(21, 61),     # wavemap shall be computed for these orders
    'n_bootstrap': 100,                     # number of bootstrap experiments
    'profile': 'gauss',                     #
    'loss_function': 'loss_1',              # weighted L2-loss
    'epsilon_sigma_bootstrap': 4,           # locations with larger uncertainty are removed
    'epsilon_sigma_clipp': 4,               # sigma clip for 1d polynomial
    'epsilon_sigma_clipp_2d' : 2,           # sigma clip for 2d polynomial
    'clipp_method' : 'pixerror',            # 'rel_std' or 'pixerror' or 'est_std'
    'n_sigma_clipp' : 100,                   # maximal number of sigma clips
    'n_sigma_clipp_2d' : 100,                   # maximal number of sigma clips
    'sigma_min' : 0.01,                       # minimial sigma to avoid overfitting
    'palette_order': 'gist_rainbow',        # palette of orders
    'order_ol': 7,                          # order polynomial in ol
    'order_o': 0,                           # order polynomial in o
}

# parameters may be added or changed using kwargs.update('param': value)

