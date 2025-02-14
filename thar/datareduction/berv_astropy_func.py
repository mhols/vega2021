
#verifie flux sur 3eme voie
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.constants import c
from astroquery.simbad import Simbad

#if query unsuccessful: "NameResolveError"
def berv(jd,object,obs_long, obs_lat, obs_alt):
    targ = SkyCoord.from_name(object)
    loc = EarthLocation(obs_long, obs_lat, obs_alt)
    #    here add half the exposure time to jd
    t = Time(jd, format='jd', scale='utc')

    #t.jd ou t.mjd

    from astropy.utils.iers import conf as iers_conf
    iers_conf.iers_auto_url = 'https://datacenter.iers.org/data/9/finals2000A.all'
    iers_conf.auto_max_age = 14.  #or none

    """
    from astropy.utils.iers import conf
    #conf.auto_max_age = None
    """

    vcorr = targ.radial_velocity_correction(kind='barycentric', obstime=t, location=loc)
    vcorr = vcorr.value##vcorrt = vcorrt.value
    ltt_bary = t.light_travel_time(targ, location=loc)

    bjd = t.tdb + ltt_bary          #barycentric Coordinate Time
    bjd = bjd.value

    # accessing the radial velocity values of a star
    # utils.py in barycorrpy uses:

    customSimbad = Simbad()
    customSimbad.add_votable_fields('rv_value')
    obj = customSimbad.query_object(object)
    rv = obj['RV_VALUE'][0] * 1000
    
    
    return bjd,vcorr,rv



#Pic du Midi
obs_long = 0.1333
obs_lat = 42.9333
obs_alt = 2869.4


#sc = SkyCoord.from_name("beta CrB")
#targ = SkyCoord.from_name("AD Leo")

#s = '/Users/boehm/Desktop/vega2021/thar/06oct21_Gam_Equ/NEO_20211006_214255_st0.fits'
#s =  '/Users/boehm/Desktop/vega2021/thar/06oct21_Gam_Equ/NEO_20211006_214656_st0.fits'
s = '/Users/boehm/Desktop/vega2021/thar/06apr23_ADLeo/NEO_20230406_221242_st0.fits'

hdulist=pyfits.open(s)
hdr= hdulist[0].header

jd = hdr.get('DATE_JUL')
dateobs = hdr['DATE-FIT']       #'2021-10-06T21:42:55.408' / Data Observation
exptime = hdr['EXPTIME']        #Integration time(sec)
stokesnum = hdr['STOKESEQ']     #Number of the seq in the serie (1..4)
object = "AD Leo"
"""
# another way to calculate julian date
t = Time(dateobs)
t.jd
"""

exptimed = exptime/(60*60*24)    #exposure time in days
jdmid = jd + exptimed/2.

bjd, vcorr, rv = berv(jdmid,object,obs_long, obs_lat, obs_alt)
        
print("barycentric julian date (bjd):  ",bjd, "barycentric velocity correction (m/s):  ", vcorr)


clum = 3e8
rv_new = rv + vcorr + rv * vcorr / clum
print("rv:",rv, "rv_new:", rv_new)
#rv = rv + vcorr * (1 + rv/c)


# delta_lambda/lambda =  delta_v/c
# lambda_new = lambda + delta_lambda
#            = lambda * (1 + delta_v/c)
# with delta_v = rv_new - rv


