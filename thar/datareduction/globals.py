"""
global parameters for the datareduction pipeline
"""
import os

# directory layout

BASEDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../')  # put jour actual base path here
DATADIR = os.path.join(BASEDIR, './datafiles')
REFFILES = os.path.join(BASEDIR, './reffiles')
TMPFILES = os.path.join(BASEDIR, './tmpfiles')

# processing parameters
kwargs = {
    'n_bootstrap': 100,                     # number of bootstrap experiments
    'profile': 'gauss',                     #
    'loss_function': 'loss_1',              # weighted L2-loss
    'save_bootstraped_data': True,          #
    'epsilon_sigma_bootstrap': 3,           # locations with larger uncertainty are removed
    'epsilon_sigma_clipp': 3,               # sigma clip for 1d polynomial
    'epsilon_sigma_clipp_2d' : 2,           # sigma clip for 2d polynomial
    'clipp_method' : 'rel_std',             # 'rel_std' or 'pixerror' or 'est_std'
    'n_sigma_clip' : 100,                   # maximal number of sigma clips
    'sigma_min' : 0.,                       # minimial sigma to avoid overfitting
    'orders': None,                         # orders to use, set None if taken from data
    'palette_order': 'gist_rainbow',        # palette of orders
    'order_ol': 6,                          # order polynomial in ol
    'order_o': 7,                           # order polynomial in o
}








}
