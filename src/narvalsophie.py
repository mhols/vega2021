from scipy.interpolate.fitpack2 import UnivariateSpline
from matplotlib import gridspec
import spectralutil as sp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys


#### some pictures for exploration...

def binned_spectrum(*specdat, **kwargs):
    nn = len(specdat)
    plt.figure('period = '+str(kwargs.get('period', sp.VEGAPERIOD)), figsize=(8*nn, 16))
    for i, spec in enumerate(specdat):
        plt.subplot(100+10*nn+i+1)
        F = np.array([1,2,1]); F = F / np.sum(F)
        quantity = spec.filtered_intensity(F)
        quantity = 1-quantity
        quantity = quantity / np.sum(quantity, axis=1)[:, np.newaxis]
        quantity = quantity - np.median(quantity, axis=0)[np.newaxis, :]
        binnedspec, bins  = sp.spectrum_matrix(
            time=spec.time-spec.time[0],
            nphase=256,
            quantity=quantity,
            **kwargs
        )
        # plt.figure(figsize=(5,10))
        plt.title(spec.name)
        plt.imshow(np.sign(binnedspec)*np.abs(binnedspec)**0.5,
                   aspect='auto', origin='lower',      cmap='gist_gray',
                   extent=[
                        np.min(spec.velocity),
                        np.max(spec.velocity),
                        0, 1]
         )
        plt.axvline(x=spec.v0)
        plt.axvline(x=spec.v0+22)
        plt.axvline(x=spec.v0+37)
        plt.axvline(x=spec.v0-22)
        plt.axvline(x=spec.v0-37)
     
def radial_velocity(*specdat, relative_depth=0.1):
    ax = plt.figure()
    colors = ['r', 'g', 'b', 'k', 'm', 'y']
    name = ''
    i = 0
    for spdat in specdat:
        col = colors[i]
        name += spdat.name
        plt.plot(np.mod(spdat.time, sp.VEGAPERIOD), 
                 spdat.rv_mean(relative_depth), 
                    'o', color=col, label=spdat.name,
                    markersize=10)
        i+=1
    plt.legend()

    plt.title(name)


def radial_velocity_correlation(*specdat, relative_depth=1.0):
    ax = plt.figure('radial velocity correlation')
    colors = ['r', 'g', 'b', 'k', 'm', 'y']
    name = ''
    i = 0
    for spdat in specdat:
        col = colors[i]
        name += spdat.name
        plt.plot(np.mod(spdat.time, sp.VEGAPERIOD), 
                spdat.rv_corr(relative_depth), 
                    'o', color=col, label=spdat.name,
                    markersize=10)
        i+=1
    plt.legend()

    plt.title(name)

def radial_velocity_bisector(*specdat, atdepth=0.5, **kwargs):
    plt.figure('radial velocity bisector {}'.format(kwargs.get('period', '')))
    colors = ['r', 'g', 'b', 'k', 'm', 'y']
    name = ''
    i = 0
    for spdat in specdat:
        col = colors[i]
        name += spdat.name
        plt.plot(np.mod(spdat.time, kwargs.get('period', sp.VEGAPERIOD)), 
                spdat.rv_bis(atdepth), 
                    'o', color=col, label=spdat.name,
                    markersize=10)
        i+=1
    plt.legend()

    plt.title(name+str(kwargs.get('period','')))

def lomb_scargel(*specs, atdepth=0.3):
    nn = len(specs)
    cpdmin, cpdmax = 0.2, 20
    oms = 2*np.pi * np.linspace(cpdmin, cpdmax, 1024)
    plt.figure(figsize=(18, nn*4))
    for i, spec in enumerate(specs):
        plt.subplot(nn*100 + 10 + i+1)
        plt.title(spec.name)
        plt.plot(oms/(2*np.pi), spec.lomb_skagel_vr_bis(oms, atdepth))

def lomb_scargel_vspan_old(*specs, uu=0.2, ul=0.3, lu=0.5, ll=0.6):
    nn = len(specs)
    cpdmin, cpdmax = 0.2, 20
    oms = 2*np.pi * np.linspace(cpdmin, cpdmax, 1024)
    plt.figure(figsize=(18, nn*4))
    for i, spec in enumerate(specs):
        plt.subplot(nn*100 + 10 + i+1)
        plt.title(spec.name+' vspan')
        plt.plot(oms/(2*np.pi), spec.lomb_scargel_vspan_old(oms, uu, ul, lu, ll))


def lomb_scargel_vspan(*specs, depth=0.9, uu=0.2, ul=0.3, lu=0.5, ll=0.6):
    nn = len(specs)
    cpdmin, cpdmax = 0.2, 20
    oms = 2*np.pi * np.linspace(cpdmin, cpdmax, 1024)
    plt.figure(figsize=(18, nn*4))
    for i, spec in enumerate(specs):
        plt.subplot(nn*100 + 10 + i+1)
        plt.title(spec.name+' vspan')
        data = spec.lomb_scargel_vspan(oms, uu, ul, lu, ll)
        plt.plot(oms/(2*np.pi), data)


def bisector_test(*specs):
    for spec in specs:
        bise = spec.bisector_borders()
        atdepth = np.linspace (0.01, 0.9, 64)
        min_intens = np.min(spec.intensity, axis=1)
        inte = spec.intensity
        v = spec.velocity
        plt.figure()
        plt.title(spec.name)
        for (i, mini) in enumerate(min_intens):
            at = mini + atdepth * (1-mini)  # linear interpolation 
            l, r = bise[i]
            plt.plot((l(at)+r(at))/2, at)
            plt.plot(v, inte[i])


def intens_bubble(*specs):
    for spec in specs:
        intens = spec.intensity
        #v = spec.velocity
        plt.figure()
        plt.title(spec.name)
        for inte in intens:
            mini = np.min(inte)
            plt.plot((inte-mini)/(1-mini))



"""
def bisector_time():
    name = "bisector_time"
    plt.figure()
    gs=gridspec.GridSpec(1, 2, width_ratios=[3, 1])
    plt.subplot(gs[0])
    plt.axis([-14,-12.75,0.,1.])
    plt.title('bisector variations')
    plt.xlabel('velocity (km/s)')
    plt.ylabel('Profile depth')
    n, d = bisector.shape
    di = 0.5 * (self.depth[1] - self.depth[0])
    for i in range(int(n/100)):
        u = np.random.uniform(-di, di, d)
        plt.plot(self.bisector[i*100, :], 1-self.depth + u, ',b')
#        for i in range(n):
#            u = np.random.uniform(-di, di, d)
#            plt.plot(self.bisector[i, :], 1-self.depth + u, ',b')

    m = self.bisector.mean(axis=0)
    s = self.bisector.std(axis=0)
    plt.plot(m, 1-self.depth, '-k')
    plt.plot(m + 1.96 * s, 1-self.depth, '-g')
    plt.plot(m - 1.96 * s, 1-self.depth, '-g')

    plt.subplot(gs[1])
    plt.title('std of bisector')
    plt.axis([0, 0.15, 0, 1])
    plt.plot(s, 1-self.depth, '-g')
    plt.yticks([])
    plt.xticks([])
    plt.xticks([0,0.1],['0','0.1'])
    plt.savefig(name + self.format)
"""

if __name__ == '__main__':
    matplotlib.rcParams.update({'font.size': 22})

    sophie2018 = sp.SpectralAnalyser('sophie_reduced.json')
    sophie2012 = sp.SpectralAnalyser('sophie12_reduced.json')
    narval = sp.SpectralAnalyser("narval_reduced.json")

    #for f in np.linspace(0.9, 1.1, 21):
    #    binned_spectrum(sophie2012, sophie2018, period=f*sp.VEGAPERIOD)

    #lomb_scargel_vspan(sophie2012, sophie2018, narval, uu=0.5, ul=0.66, lu=0.75, ll=0.9)
    #lomb_scargel_vspan_old(sophie2012, sophie2018, narval, uu=0.5, ul=0.34, lu=0.25, ll=0.1)
    #radial_velocity(narval, sophie2018, sophie2012, relative_depth=0.8)
    #radial_velocity_correlation(narval, sophie2018, sophie2012, relative_depth=0.8)
    #radial_velocity_bisector(narval, sophie2018, sophie2012, atdepth=0.5)


    # bisector_test(sophie2012, sophie2018, narval)
    intens_bubble(sophie2012, sophie2018, narval)
    plt.show()
