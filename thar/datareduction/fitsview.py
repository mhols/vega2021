from astropy.io import pyfits
from numpy import *
import extract

class FitsView:

    def __init__(self, fitsfile):
        self.fitsfile = fitsfile
    
    @property
    def is_star(self):
        return extract.is_star(self.fitsfile)
    
    @property
    def is_thorium(self):
        return extract.is_thorium(self.fitsfile)
    
    
