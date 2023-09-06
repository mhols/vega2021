
#verifie flux sur 3eme voie
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.constants import c

"""
from astropy.utils.iers import conf as iers_conf
iers_conf.iers_auto_url = 'https://astroconda.org/aux/astropy_mirror/iers_a_1/finals2000A.all'   # outdated see: https://github.com/mzechmeister/serval/issues/47#issuecomment-997350597
iers_conf.iers_auto_url = 'https://datacenter.iers.org/data/9/finals2000A.all'
iers_conf.auto_max_age = None
"""

import barycorrpy

#see also http://astroutils.astronomy.ohio-state.edu/exofast/barycorr.html

#Pic du Midi
obs_long = 0.1333
obs_lat = 42.9333
obs_alt = 2869.4

loc = EarthLocation.from_geodetic(obs_long, obs_lat, obs_alt)




s = '/Users/boehm/Desktop/vega2021/thar/06oct21_Gam_Equ/NEO_20211006_214255_st0.fits'


hdulist=pyfits.open(s)
hdr= hdulist[0].header

jd = hdr.get('DATE_JUL')
dateobs = hdr['DATE-FIT']       #'2021-10-06T21:42:55.408' / Data Observation
exptime = hdr['EXPTIME']        #Integration time(sec)
stokesnum = hdr['STOKESEQ']     #Number of the seq in the serie (1..4)
exptimed = exptime/(60*60*24)   #exposure time in days
object = hdr['OBJECT']          #name of object

sc = SkyCoord.from_name(object)
ra = sc.ra.value
dec =  sc.dec.value
targ = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')

#    here add half the exposure time to jd
#t = Time(jd + 0.5 * exptimed, format='jd', scale='utc')
t = Time(jd, format='jd', scale='utc')
loc = EarthLocation.from_geodetic(obs_long, obs_lat, obs_alt)


"""
   jd_utc: Julian date (UTC)
   ra: RA (J2000) [deg]
   dec: Dec (J2000) [deg]
   obsname: Observatory name (overrides coordinates if set)
   lat: Observatory latitude  [deg]
   lon: Observatory longitude (E) [+/-360 deg]
   elevation: Observatory elevation [m]
   pmra: Proper motion (RA*cos(Dec)) [mas/yr]
   pmdec: Proper motion (Dec) [mas/yr]
   parallax: Parallax [mas]
   rv: Radial velocity (within 100 km/s) [m/s]
   zmeas: Measured redshift
   epoch: Epoch (default 2448348.56250, J2000)
   tbase: Baseline subtracted from times (default 0.0)
"""

kwargs = dict(ra=ra, dec=dec, epoch=2000.0, pmra=0.0, pmdec=0.0, rv=0.0,
              lat=obs_lat, longi=obs_long, alt=obs_alt)
bervres = barycorrpy.get_BC_vel(jd, zmeas=0.0, **kwargs)
bjdres = barycorrpy.utc_tdb.JDUTC_to_BJDTDB(jd, **kwargs)
ltt_bary = t.light_travel_time(targ, location=loc)

berv = bervres[0][0]
bjd = bjdres[0][0]
#works also:     bjdastropy=t.tdb.value+ltt_bary.value

print(bjd,berv)


#return berv, bjd

