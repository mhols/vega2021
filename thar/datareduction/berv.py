
#verifie flux sur 3eme voie
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

from astropy.constants import c


loc = EarthLocation(42.9333,0.1333,2869.4)   #Pic du Midi
#sc = SkyCoord.from_name("Gamma Equ")
sc = SkyCoord.from_name("beta CrB")

s = '/Users/boehm/Desktop/vega2021/thar/06oct21_Gam_Equ/NEO_20211006_214255_st0.fits'


hdulist=pyfits.open(s)
hdr= hdulist[0].header

jd = hdr.get('DATE_JUL')
dateobs = hdr['DATE-FIT']       #'2021-10-06T21:42:55.408' / Data Observation
exptime = hdr['EXPTIME']        #Integration time(sec)
stokesnum = hdr['STOKESEQ']     #Number of the seq in the serie (1..4)

t = Time(jd, format='jd', scale='utc')
#t.jd ou t.mjd

vcorr = sc.radial_velocity_correction(kind='barycentric', obstime=t, location=loc)
vcorr = vcorr.value

bjd = t.tcb           #barycentric Coordinate Time
bjd = bjd.value


#rv = rv + vcorr + rv * vcorr / c
print("barycentric julian date (bjd):  ",bjd, "barycentric velocity correction (m/s):  ", vcorr)
