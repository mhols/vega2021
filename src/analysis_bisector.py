'''
Created on May 22, 2014

@author: hols
'''

from spectralutil import *
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectral

PROGDIR = os.path.dirname(__file__)
DATADIR = os.path.join(PROGDIR, '../data')

#DATAFILE = os.path.join(DATADIR, 'filematrix.dat')
DATAFILE = os.path.join(DATADIR, 'filematrix_nn.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_takeda.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_strong_s.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_strong_os.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_superstrong.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_superstrong_notnorm.dat')
#DATAFILE = os.path.join(DATADIR,'filematrix_453.dat')

time, velocity, intensity, list_time, list_inte = load_data(DATAFILE)

# overall scaling
#intensity = 1 - (intensity-1)/(intensity.min()-1)

dmin, dmax = 0.2, 0.7 # bounds for bisector analysis
nn = 50                # levels for bisector

d = np.linspace(dmin, dmax, nn)

#allbisec = np.zeros((len(time), nn))
#for i in xrange(len(time)):
#    allbisec[i, :] = bisector_v(velocity, intensity[i], d)
        




def doit00():
    """
    plotting the prfiles
    """
    plt.figure()
    for i in xrange(len(time)):
        plt.plot(velocity, intensity[i])

def doit0():
    """
        plots all the bisectors
    """
    
    print "plotting all the bisectors"
    nn = 50
    tmp = np.zeros((len(time), nn))
    
    d = np.linspace(dmin, dmax, nn)
    for i in xrange(len(time)):
        tmp[i,:] = bisector_v(velocity, intensity[i], d)
        plt.plot(tmp[i,:], d, '.k')
    
    plt.plot(tmp.mean(axis=0), d, '-r', linewidth=5)


def doit1():
    """
    Lomb scargal analyse der bisektoren
    """
    
    print "Lomb Scargel analyse of the bisectors"
    print "--->start"
    
    
    nn = 50
    d = np.linspace(dmin, dmax, nn)
    
    tmp = np.zeros((len(time), nn))
    for i in xrange(len(time)):
        tmp[i, :] = bisector_v(velocity, intensity[i], d)
        
    m = tmp.mean(axis=0)
    
    nf = 1024
    freqs = np.linspace(2*np.pi*0.001, 2*np.pi*2.5, nf)
    
    res = np.empty((nf, nn))
    
    print "--->computing the specra"
    for i in xrange(nn):
        spk = spectral.lombscargle(time, tmp[:,i]-m[i], freqs)
        res[:,i] = np.abs(spk)
        
    plt.figure()
    plt.title("LS of bisector")
    res/=res.max()
    #res = np.log(10*res + 1)
    plt.contourf(d, freqs/(2*np.pi), res)

def doit2():
    """
        average bisector and variability
    """
    upper = (0.5, 0.7)
    lower = (0.15, 0.3)
    
    nn = 20
    tmp = np.zeros(len(time))
    
    du = np.linspace(upper[0], upper[1], nn)
    dl = np.linspace(lower[0], lower[1], nn)
    
    for i in xrange(len(time)):
        bu = bisector_v(velocity, intensity[i], du)
        bl = bisector_v(velocity, intensity[i], dl)
        tmp[i] = np.median(bu) - np.median(bl)
    
    nf=1024
    freqs = np.linspace(2*np.pi*0.001, 2*np.pi*15, nf)
    
    np.savetxt('vspan.dat', np.column_stack((time, tmp)), fmt=['%.5f', '%.7f'])
    plt.figure()
    plt.plot(freqs/(2*np.pi), np.abs(spectral.lombscargle(time, tmp-tmp.mean(), freqs)))
    

def _doit3():
    """
    plotting vspam as function of time
    """
    upper = (0.6, 0.8)
    lower = (0.15, 0.3)
    
    nn = 20
    tmp = np.zeros(len(time))
    
    du = np.linspace(upper[0], upper[1], nn)
    dl = np.linspace(lower[0], lower[1], nn)
    
    for i in xrange(len(time)):
        bu = bisector_v(velocity, intensity[i], du)
        bl = bisector_v(velocity, intensity[i], dl)
        tmp[i] = np.median(bu) - np.median(bl)
    

    return tmp

def doit3():
    tmp = _doit3()
    np.savetxt('vspantime.dat', np.column_stack((time,tmp)))
    plt.figure()
    plt.title('')
    plt.plot(time, tmp)

def doit4():
    """
    relative area
    """
    res = np.zeros(len(time))
    for i in xrange(len(time)):
        area = (velocity.max()-velocity.min()) * (1-intensity.min())
        res[i] = (1-intensity[i,:]).sum() / area
    
    plt.figure()
    plt.title('relative area')
    plt.plot(time, res)
    
    nf=1024
    freqs = np.linspace(2*np.pi*0.001, 2*np.pi*15, nf)
    
    plt.figure()
    plt.title('LS analysis of relative area')
    plt.plot(freqs/(2*np.pi), np.abs(spectral.lombscargle(time, res-res.mean(), freqs)))
    
def doit5():
    """
    relative area with individual night polynomail
    trend remoed
    """
    res = np.zeros(len(time))
    for i in xrange(len(time)):
        area = (velocity.max()-velocity.min()) * (1-intensity.min())
        res[i] = (1-intensity[i,:]).sum() / area
    
    res = remove_night_trend(time, res)
    
    plt.figure()
    plt.title('relative area')
    plt.plot(time, res)
    
    nf=1024
    freqs = np.linspace(2*np.pi*0.001, 2*np.pi*15, nf)
    
    plt.figure()
    plt.title('LS analysis of relative area')
    plt.plot(freqs/(2*np.pi), np.abs(spectral.lombscargle(time, res-res.mean(), freqs)))
    
def doit6():
    """
    Window function of observations
    """
    nf=1024
    freqs = np.linspace(2*np.pi*0.001, 2*np.pi*15, nf)
    
    tmp = np.ones(len(time))
    plt.figure()
    plt.title("LS of observation window")
    plt.plot(freqs/(2*np.pi), np.abs(spectral.lombscargle(time, tmp, freqs)))
    
def doit7():
    """
    plotting nightly trends
    """
    n=0
    for time,inte in zip (list_time, list_inte):
        plt.figure()
        n+=1
        plt.title("night %d"%n)
        print time.shape
        print inte.shape
        rest, p, pp = remove_trend(time, inte)
        d,m = pp.shape
        offset = np.linspace(0,0.001, m)
        plt.plot(time, pp-pp[0]+offset)
        plt.ylim(-0.002,0.002)
    
def doit8():
    """
    pixelwise LS analysis with nightly trend removed
    """
    inte = remove_night_trend(time, intensity)
    nt, nv = inte.shape
    nf = 1024 # number of LS frequencies
    freqs = np.linspace(2*np.pi*0.001, 2*np.pi*15, nf)
    res = np.empty( (nf, nv))

    for i in xrange(nv):
        print "working on velocity bin nr ", i
        spk = spectral.lombscargle(time, inte[:,i], freqs)
        res[:,i] = np.abs(spk)
    
    plt.figure()
    res/=res.max()
    #res = np.log(10*res + 1)
    plt.contourf(velocity, freqs/(2*np.pi), res)


def _doit9():
    relative_depth = 0.6
    depth = 1- (1-intensity.min())*relative_depth
    
    nt,nv = intensity.shape
    
    res = np.zeros(nt)
    for i in xrange(nt):
        depth = 1- (1-intensity[i,:].min())*relative_depth
        print depth, intensity[i,:].min()
        I = np.where(intensity[i,:]<= depth)[0]
        iii = intensity[i,I]
        iii = iii.max()-iii
        res[i] = np.sum(iii*velocity[I]) / np.sum(iii)    
    return res

def doit9():
    """
    radial velocity based on mean
    """
    res = _doit9()
    np.savetxt('lsvrad.dat', np.column_stack((time, res)), fmt=['%.5f', '%.7f'])

    plt.figure()
    plt.plot(time, res)

def doit10():
    res = _doit9()
    res -= res.mean()
    nf = 1024 # number of LS frequencies
    freqs = np.linspace(2*np.pi*0.001, 2*np.pi*15, nf)
    spk = spectral.lombscargle(time, res, freqs)
    
    plt.figure()
    plt.plot( freqs/(2*np.pi), np.abs(spk))
    
def doit11():
    vr = _doit9()
    vs = _doit3()
    plt.figure()
    plt.plot(vr,vs, 'o')
    
if __name__ == "__main__":
    doit00() # the spectra
    #doit0()  # the bisectors
    #doit1()   # LS Analyse of bisectors
    #doit2()   #LS Analyse of vspan
    #doit3()   # vspam analyse
    #doit4()
    #doit5()
    #doit6() # window function od obs
    #doit7() # plotting the nightly polynomial trend
    #doit8() # LS analysis with nightly trend for each bin removed
    #doit9() # radial velocity based on mean
    #doit10() # LS abakysus if resutl of 9
    #doit11()
    plt.show()  