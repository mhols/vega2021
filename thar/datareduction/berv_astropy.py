
#verifie flux sur 3eme voie
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.constants import c

#if query unsuccessful: "NameResolveError"

loc = EarthLocation(0.1333,42.9333,2869.4)   #Pic du Midi
#targ = SkyCoord.from_name("Gamma Equ")
#sc = SkyCoord.from_name("beta CrB")
targ = SkyCoord.from_name("AD Leo")

#s = '/Users/boehm/Desktop/vega2021/thar/06oct21_Gam_Equ/NEO_20211006_214255_st0.fits'
s =  '/Users/boehm/Desktop/vega2021/thar/06apr23_ADLeo/NEO_20230406_215211_st0.fits'

hdulist=pyfits.open(s)
hdr= hdulist[0].header

jd = hdr.get('DATE_JUL')
dateobs = hdr['DATE-FIT']       #'2021-10-06T21:42:55.408' / Data Observation
exptime = hdr['EXPTIME']        #Integration time(sec)
stokesnum = hdr['STOKESEQ']     #Number of the seq in the serie (1..4)

exptimed = exptime/(60*60*24)    #exposure time in days


#    here add half the exposure time to jd
t = Time(jd, format='jd', scale='utc')

#t.jd ou t.mjd

"""
from astropy.utils.iers import conf
conf.auto_max_age = None
"""

vcorr = targ.radial_velocity_correction(kind='barycentric', obstime=t, location=loc)
vcorr = vcorr.value##vcorrt = vcorrt.value
ltt_bary = t.light_travel_time(targ, location=loc)

bjd = t.tdb + ltt_bary          #barycentric Coordinate Time
bjd = bjd.value


#rv = rv + vcorr + rv * vcorr / c
print("barycentric julian date (bjd):  ",bjd, "barycentric velocity correction (m/s):  ", vcorr)

