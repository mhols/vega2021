'''
Created on Jul 9, 2014

@author: hols
'''

from spectralutil import *
from dycos.lmm import *
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectral

PROGDIR = os.path.dirname(__file__)
DATADIR = os.path.join(PROGDIR, '../data')

#DATAFILE = os.path.join(DATADIR, 'filematrix.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_nn.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_takeda.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_strong_s.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_strong_os.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_superstrong.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_superstrong_notnorm.dat')
#DATAFILE = os.path.join(DATADIR,'filematrix_453.dat')
DATAFILE = os.path.join(DATADIR, 'filematrix_shifted.dat')

time, velocity, intensity, list_time, list_inte = load_data(DATAFILE)
ntimes = time.size

tim = time-time[0]
vsp = vspan( upper=(0.6,0.8), lower=(0.15,0.3), time=tim, velocity=velocity, intensity=intensity)
vrd = vrad( relative_depth=0.5, velocity=velocity, intensity=intensity)

def Fsystem(t, freq, nharm=1):
    """
    returns a system matrix for a given frequency freq
    and times t
    """
    n = t.size
    F = np.zeros((n, 2*nharm))
    for i in xrange(nharm):
        F[:,2*i  ] = np.cos ( 2 * (i+1) * np.pi * freq * t)
        F[:,2*i+1] = np.sin ( 2 * (i+1) * np.pi * freq * t)
    
    return F

def doit0():
    tim=time-time[0] #np.sin(2*np.pi*0.7*tim)#1+np.sin(2*np.pi*1.3*tim)+
    #tim = np.linspace(0,50,1024)
    fmin, fmax = 0.1, 5
    freqs = np.linspace(fmin, fmax, 128)
    res = []
    for f in freqs:
        F = np.column_stack((np.ones((len(tim),1)), Fsystem(tim, f, 3)))
        calc = LM(F, vsp)
        res.append(calc.getLogLHMarginalFixedMarginalSigma())
    
    res=np.asarray(res)
    res = res-res.max()
    plt.figure()
    plt.plot(freqs, res)
    
def doit1():
    tim=time-time[0] #np.sin(2*np.pi*0.7*tim)#1+np.sin(2*np.pi*1.3*tim)+
    #tim = np.linspace(0,50,1024)
    fmin, fmax = 0.1, 5
    freqs = np.linspace(fmin, fmax, 128)
    res = []
    for f in freqs:
        F = np.column_stack((np.ones(len(tim)), Fsystem(tim, f, 2)))
        calc = LM(F, vrd)
        res.append(calc.getLogLHMarginalFixedMarginalSigma())
    
    res=np.asarray(res)
    res = res-res.max()
    plt.figure()
    plt.plot(freqs, res)
    

def doit2():
    tim=time-time[0] #np.sin(2*np.pi*0.7*tim)#1+np.sin(2*np.pi*1.3*tim)+
    #tim = np.linspace(0,50,1024)
    fmin, fmax = 0.1, 3
    freqs = np.linspace(fmin, fmax, 128)
    
    plt.figure()
    plt.plot(freqs, np.abs(spectral.lombscargle(tim, vsp, 2*np.pi*freqs)))
    
def doit3():
    tim=time-time[0] #np.sin(2*np.pi*0.7*tim)#1+np.sin(2*np.pi*1.3*tim)+
    #tim = np.linspace(0,50,1024)
    fmin, fmax = 0.1, 3
    freqs = np.linspace(fmin, fmax, 128)
    res = []
    for f in freqs:
        F = np.column_stack((np.ones(len(tim)), Fsystem(tim, f, 2)))
        calc = LM(F, vsp)
        res.append(calc.getLogLHMarginalFixedMarginalSigma())
    
    res=np.asarray(res)
    res = res-res.max()
    plt.figure()
    plt.plot(freqs, np.exp(res))


def doit4():
    """
    Likelihood ratio test
    """
    tim=time-time[0] #np.sin(2*np.pi*0.7*tim)#1+np.sin(2*np.pi*1.3*tim)+
    #tim = np.linspace(0,50,1024)
    F =  np.ones((len(tim),1))
    calc = LM(F, vsp)
    lh1 = calc.getLogLHMarginalFixedMarginalSigma()
    F = F = np.column_stack((np.ones(len(tim)), Fsystem(tim, 1.45, 1)))
    calc = LM(F, vsp)
    lh2 = calc.getLogLHMarginalFixedMarginalSigma()
    
    print "log(lh2/lh1) ", lh2-lh1
    
def doit5():
    plt.figure()
    plt.plot(vrd,vsp, 'o')

if __name__=="__main__":
    
    doit0()
    #doit1()
    #doit2()
    #doit3()
    #doit4() $likelihood ratio
    doit5()
    plt.show()