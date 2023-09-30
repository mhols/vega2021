"""
Running vega from scratch
"""
import unittest
from nextra import *
import os
import pathlib

# recompute the reference spectrum
setup_reference()



class TestVega(unittest.TestCase):

    
    def test_luna(self):
        # adjusting parameters
        kwargs = settings_reference.get_kwargs()
        kwargs['SETTING_ID'] = 'Moon'
        kwargs['HIGHEXP'] = 15
        kwargs['LOWEXP'] = 4
        kwargs['USE_PICKED_CENTRAL_POSITIONS'] = False

        f_luna = pathlib.Path(__file__)

        f_luna = f_luna.absolute().parents[3] / 'lune_raw/NEO_20200202_182642_st0.fits'

        store.rm(f_luna)

        OK = False
        try:
            reduce_star(f_luna, **kwargs)
            OK = True
        except Exception as ex:
            OK = False
        
        self.assertTrue(OK)
        
    def test_Pegasus(self):
        # adjusting parameters
        kwargs = settings_reference.get_kwargs()
        kwargs['SETTING_ID'] = 'Pegasus'
        kwargs['ORDERS'] = list(range(21,57))   # order 56 not present
        kwargs['USE_PICKED_CENTRAL_POSITIONS'] = False



        f = pathlib.Path(__file__)

        f = f.absolute().parents[3] / '51Peg_raw/2020_0913/NEO_20200913_233957_st0.fits'

        store.rm(f)

        OK = False
        try:
            reduce_star(f, **kwargs)
            OK = True
        except Exception as ex:
            OK = False
        
        self.assertTrue(OK)
    

    def test_HD(self):
        # adjusting parameters
        kwargs = settings_reference.get_kwargs()
        kwargs['SETTING_ID'] = 'HD'
        kwargs['USE_PICKED_CENTRAL_POSITIONS'] = False


        f_hd = pathlib.Path(__file__)

        f_hd = f_hd.absolute().parents[3] / 'HD185144_raw/15apr22/NEO_20220416_011157_st0.fits'

        store.rm(f_hd)

        OK = False
        try:
            reduce_star(f_hd, **kwargs)
            OK = True
        except Exception as ex:
            OK = False
        
        self.assertTrue(OK)
 
if __name__ == '__main__':
    unittest.main()

