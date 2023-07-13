from extract import *
from settings_reference import kwargs as kwargs_reference
import pickle
import util
from units import *
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.constants import c
import barycorrpy


#
#fitsname = '/Users/boehm/Desktop/vega2021/thar/06apr23_ADLeo/NEO_20230406_215211_st0.fits'
    
def photometry(fitsname):
    s = fitsname
    #s ='/Users/boehm/Desktop/vega2021/thar/06oct21_Gam_Equ/NEO_20211006_214054_st0.fits'
    hdulist=pyfits.open(s)
    hdr= hdulist[0].header
    jd = hdr.get('DATE_JUL')
    exptime = hdr['EXPTIME']        #Integration time(sec)
    exptimed = exptime/(60*60*24)    #exposure time in days



    counts=hdulist[1].data['count']
    time=hdulist[1].data['julien']

    fractime=time-time[0]

    btime=0.
    sumcount=0.

    for i in range(fractime.size):
        btime=btime+counts[i]*fractime[i]
        sumcount=sumcount+counts[i]
        
    btime=btime/sumcount

    jd_firstmoment = time[0] + btime

    print("total count time (s): ", (fractime[-1]-fractime[0])*24.*3600., "total exposure time (s):", exptime)
    print("half count time (s) :", (fractime[-1]-fractime[0])*24.*3600.*0.5)
    print("center of exposure (s): ", (btime-fractime[0])*24.*3600., "half exposure time: ",exptime/2.)
    print("time_beg (jd): ", time[0], "jd_firstmoment (jd): ", jd_firstmoment, "time_end (jd): ",time[-1])
    print("center of exposure (jd): ", time[0]+(time[-1]-time[0])/2.)

    return jd_firstmoment

### should not be needed.....
"""
def photometry2(fitsfile):

    #s ='/Users/boehm/Desktop/vega2021/thar/06oct21_Gam_Equ/NEO_20211006_214054_st0.fits'
    hdulist=pyfits.open(fitsfile)
    hdr= hdulist[0].header
    jd = hdr.get('DATE_JUL')
    exptime = hdr['EXPTIME']        #Integration time(sec)
    exptimed = exptime/(60*60*24)    #exposure time in days



    counts=hdulist[1].data['count']
    time=hdulist[1].data['julien']

    fractime=time-time[0]

    btime=0.
    sumcount=0.

    for i in range(fractime.size):
        btime=btime+counts[i]*fractime[i]
        sumcount=sumcount+counts[i]
        
    btime=btime/sumcount

    jdfirstmoment = time[0] + btime
    # self.jdfirstmoment = jdfirstmoment

    print("total count time (s): ", (fractime[-1]-fractime[0])*24.*3600., "total exposure time (s):", exptime)
    print("half count time (s) :", (fractime[-1]-fractime[0])*24.*3600.*0.5)
    print("center of exposure (s): ", (btime-fractime[0])*24.*3600., "half exposure time: ",exptime/2.)
    print("time_beg (jd): ", time[0], "jd_firstmoment (jd): ", jdfirstmoment, "time_end (jd): ",time[-1])
    print("center of exposure (jd): ", time[0]+(time[-1]-time[0])/2.)

    return jdfirstmoment
"""

def barycorr(fitsfile, obsname='TBL', method='barycorrpy', julbase='juld', pmra=0., pmdec=0., parallax=0., rv=0., zmeas=0., epoch=2451545.0, tbase=0.):
        
    """
    obsname =   here TBL if nothing else is specified
    method  =   barycorrpy or astropy
    julbase =   juld:  means normal jd_utc, half of extimed is to be added
                jdcentroid: the true first moment of the flux is used using photometry2
                
    ra: RA (J2000) [deg]
    dec: Dec (J2000) [deg]
    
    pmra: Proper motion (RA*cos(Dec)) [mas/yr]
    pmdec: Proper motion (Dec) [mas/yr]
    parallax: Parallax [mas]
    rv: Radial velocity (within 100 km/s) [m/s]
    zmeas: Measured redshift
    epoch: Epoch (default 2448348.56250, J2000)
    tbase: Baseline subtracted from times (default 0.0)
    """
    s = fitsfile
    
    # TODO use methods allready implemented in extract  DRY DRY!

    hdulist=pyfits.open(s)
    hdr= hdulist[0].header
    
    jd_utc = hdr.get('DATE_JUL')
    exptime = hdr['EXPTIME']        #Integration time(sec)
    exptimed = exptime/(60*60*24)    #exposure time in days
    starname = hdr.get('OBJECT')
    ra = hdr.get('RA')
    dec = hdr.get('DEC')
    
    
    # TODO: create a class / setting entry observatory
    # here are prvi
    if obsname=='TBL':
        obs_lat = 42.9333
        obs_lon = 0.1333
        obs_elevation = 2869.4

    loc = EarthLocation(obs_lon,obs_lat,obs_elevation)
  
  
    try:
        targ = SkyCoord.from_name(starname)
        ra = targ.ra.value
        dec = targ.dec.value
        
    except:
        #some stellar fitsfiles do have ra,dec = -999.0 ie not provided
        #if offline and -999.0 the program cannot work!!
        #in that case nothing works anymore
        #if ra and dec are in deg.fraction then icrs, otherwise xx:xx:xx fk4
        targ = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='fk4')
        ra = targ.ra.value
        dec = targ.dec.value

    #    here add half the exposure time to jd
    
    if julbase=='juld':
        jd =  jd_utc + 0.5 * exptimed
        t = Time(jd, format='jd', scale='utc')
    else:
        jd = photometry(fitsfile) # self.jdfirstmoment
        t = Time(jd, format='jd', scale='utc')
        

    if method=='astropy':

        vcorr = targ.radial_velocity_correction(kind='barycentric', obstime=t, location=loc)
        vcorr = vcorr.value
        ltt_bary = t.light_travel_time(targ, location=loc)

        bjd = t.tdb + ltt_bary          #barycentric Coordinate Time
        bjd = bjd.value


        #rv = rv + vcorr + rv * vcorr / c
        
        rvmeas = zmeas * c.value
        rv = rvmeas + vcorr * (1 + rvmeas/c.value)
        berv = vcorr
        
        #self.berv = berv
        #self.bjdtdb = bjd
    
        print("barycentric julian date (bjd):  ",bjd)
        print("rvmeas: ", rvmeas, "barycentric velocity correction - berv (m/s):  ", berv, "rv (m/s): ", rv)

    if method=='barycorrpy':
        kwargs = dict(ra=ra, dec=dec, epoch=2000.0, pmra=0.0, pmdec=0.0, rv=0.0,
              lat=obs_lat, longi=obs_lon, alt=obs_elevation)
             
        
        bervres = barycorrpy.get_BC_vel(jd, zmeas=0.0, **kwargs)
        bjdres = barycorrpy.utc_tdb.JDUTC_to_BJDTDB(jd, **kwargs)

        berv = bervres[0][0]
        bjd = bjdres[0][0]
        
        #self.berv = berv
        #self.bjdtdb = bjd
        #c = np.array(c)
        rvmeas = zmeas * c.value
        vcorr = berv
        rv = rvmeas + vcorr * (1 + rvmeas/c.value)
    
        print("barycentric julian date (bjd):  ",bjd)
        print("rvmeas: ", rvmeas, "barycentric velocity correction - berv (m/s):  ", berv, "rv (m/s): ", rv)

    return berv, bjd
    





