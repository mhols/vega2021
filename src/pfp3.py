'''
cration jan 24 2025
@author: hols and boehm
'''

from matplotlib import gridspec
import os
from scipy.signal import lombscargle
from scipy.stats import gaussian_kde
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from spectralutil_old import *


PROGDIR = os.path.dirname(__file__)
DATADIR = os.path.join(PROGDIR, '../data')
DATADIR = os.environ.get('DATADIR', DATADIR)

DATAFILE7 = os.path.join(DATADIR, 'Vega_cfht.dat')
DATAFILE8 = os.path.join(DATADIR, 'Vega_9500.40.03-10.dat')  # lines between 0.3 and 1.0 depth, medium
DATAFILE9 = os.path.join(DATADIR, 'Vega_tbl10.dat')  # lines between 0.3 and 1.0 depth, medium
DATAFILE10 = os.path.join(DATADIR, 'Vega_2018_0310.dat')  # Sophie lines between 0.3 and 1.0 depth, medium

DATAFILE11 = os.path.join(DATADIR, 'Vega_Narval_2018_031.dat')  # Narval lines between 0.3 and 1.0 depth, medium
DATAFILE12 = os.path.join(DATADIR, 'Vega_2024.dat')
#DATAFILE12 = os.path.join(DATADIR, 'Vega_2023_maskvega_folsom.clean.old.clean')
#DATAFILE12 = os.path.join(DATADIR, 'Vega_2023_maskvega_folsom.clean')
#DATAFILE11 = os.path.join(DATADIR, 'test.dat')  # lines between 0.3 and 1.0 depth, medium
DATAFILE13 = os.path.join(DATADIR, 'Vega_2018_maskvega_folsom.clean_07.clean')
DATAFILE14 = os.path.join(DATADIR, 'Vega_2018_maskvega_folsom.clean_1.7.clean')
os.path.join(DATADIR, 'Vega_2018_0310.dat')

ames = ["full", "single"]

ranges = {}
ranges[DATAFILE7] = None
ranges[DATAFILE8] = (72, 112)
ranges[DATAFILE9] = None
ranges[DATAFILE10] = (72, 112)
ranges[DATAFILE11] = (72, 112)
ranges[DATAFILE12] = (61, 104)
ranges[DATAFILE13] = None
ranges[DATAFILE14] = None


vranges = {}
vranges[DATAFILE7] = (-60.0, 40.0)
vranges[DATAFILE8] = None
vranges[DATAFILE9] = (-60.0, 40.0)
vranges[DATAFILE10] = (-60.0, 40.0)
vranges[DATAFILE11] = (-60.0, 40.0)
vranges[DATAFILE12] = (-60.0, 40.0)
vranges[DATAFILE13] = (-60.0, 40.0)
vranges[DATAFILE14] = (-60.0, 40.0)

DATAFILE = DATAFILE11 # DATAFILE8


class Experiment:
    def __init__(self, **kwargs):
        self.kwargs = kwargs



VEGA_2018_SOPHIE = Experiment(
    DATAFILE = os.path.join(DATADIR, 'Vega_2018_0310.dat'), # Sophie lines between 0.3 and 1.0 depth, medium
    vrange = [-60.0, 40.0],
    vrad = -13.4,  
    prot = 0.678, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5','s6', 's123456','s123','s456'], 
        [[0], [1], [2], [3], [4], [5], [0,1,2,3,4,5], [0,1,2],[3,4,5]])
)



experiment = VEGA_2018_SOPHIE.kwargs

class Pictures(object):

    def __init__(self, *args):

        # parameters

        self.normalize = False  # normalise spectra to interval (0,1) ?
        self.r_depth = 0.5  # depth for
#        self.upper, self.lower = (0.35, 0.5), (0.15, 0.25)  # limits for vspan
        self.upper, self.lower = (0.35, 0.5), (0.1, 0.25)  # limits for vspan
#        self.extension = (0.15, 0.3) # limits of bisector area for median calculation vrad_bis

        self.rotperiod =  0.678 # 0.66149 #0.678# #0.678 #rotation period

        self.extension = (0.15, 0.3) # limits of bisector area for median calculation vrad_bis
        self.d0, self.d1 = (0.1, 0.9)  # limits for bisector
        self.cpd = [24.0 / 12.5 ]  # cycles per day, where to plot a vertical lline (FRot
        self.cpdl = 1. / (0.678 + 0.036)  # cycles per day, low, Alina
        self.cpdh = 1. / (0.678 - 0.029)  # cycles per day, high, Alina
        self.cpd2l = 1. / 0.678 + 1. / (0.678 + 0.036)
        self.cpd2h = 1. / 0.678 + 1. / (0.678 - 0.029)
        self.cpd3l = 2. / 0.678 + 1. / (0.678 + 0.036)
        self.cpd3h = 2. / 0.678 + 1. / (0.678 - 0.029)
        self.cpd4l = 3. / 0.678 + 1. / (0.678 + 0.036)
        self.cpd4h = 3. / 0.678 + 1. / (0.678 - 0.029)
        self.cpd5l = 4. / 0.678 + 1. / (0.678 + 0.036)
        self.cpd5h = 4. / 0.678 + 1. / (0.678 - 0.029)
        self.cpd6l = 5. / 0.678 + 1. / (0.678 + 0.036)
        self.cpd6h = 5. / 0.678 + 1. / (0.678 - 0.029)
        self.cpd7l = 6. / 0.678 + 1. / (0.678 + 0.036)
        self.cpd7h = 6. / 0.678 + 1. / (0.678 - 0.029)
        self.cpd8l = 7. / 0.678 + 1. / (0.678 + 0.036)
        self.cpd8h = 7. / 0.678 + 1. / (0.678 - 0.029)
        self.cpd9l = 8. / 0.678 + 1. / (0.678 + 0.036)
        self.cpd9h = 8. / 0.678 + 1. / (0.678 - 0.029)
        self.cpd10l = 9. / 0.678 + 1. / (0.678 + 0.036)
        self.cpd10h = 9. / 0.678 + 1. / (0.678 - 0.029)

        self.cpdnew = 1.78

        self.cpdb = 1.6062959  # cycles per day, Butkovskaya
        self.nfreq = 1024  # number of frquencies for spectral analysis
        self.min_cpd, self.max_cpd = 0.2, 50.  # limits of frequency analysis (cycles per day)
        self.format = ".pdf"  # all pictures as .pdf files

        # laoding data and compute some quantities

        self.comment = "description of plots\n"
        self.analyzer = SpectralAnalyser(experiment['DATAFILE'],
                                         experiment['normalize'],
                                         vrange=experiment['vrange'],
                                         noiselevel=experiment['noiselevel'],
                                        )
        ## comodity fields
        self.time = self.analyzer.time
        self.nspec = len(self.time)
        self.inte = self.analyzer.intensity
        self.velocity = self.analyzer.velocity
        self.cycles_per_day = np.linspace (self.min_cpd, self.max_cpd, self.nfreq)  # the frequencies
        self.freq = 2 * np.pi * self.cycles_per_day

        # quantities that shall be computed "on demand"
        # this is realized through @properites

        self._intens_mean = None
        self._vrad_eqwidth = None
        self._vrad_mean = None
        self._vrad_corr = None
        self._vspan = None
        self._vrad_bis = None
        self._vrad_skew = None
        self._vrad_std = None
        self._ls_all3 = None
        self._eqwidth = None
        self._ls_vrad_skew = None
        self._ls_vrad_mean = None
        self._ls_vrad_bis = None
        self._ls_vrad_corr = None
        self._ls_vspan = None
        self._ls_eqwidth = None
        self._bisector = None
        self._depth = None
        self._window = None

    def saveData(self, filename, data):
        """
        saves the computed data in files
        """
        np.savetxt(filename, np.column_stack(data), fmt='%12.8f')

    @property
    def vrad_eqwidth(self):
        if None == self._vrad_mean:
            self._vrad_eqwidht = self.analyzer.rv_eqwidth(relative_depth=self.r_depth)
        return self._vrad_mean

    @property
    def vrad_mean(self):
        if None == self._vrad_mean:
            self._vrad_mean = self.analyzer.rv_mean(relative_depth=self.r_depth)
        return self._vrad_mean

    @property
    def intens_mean(self):
#       if None == self._intens_mean:
#            self._intens_mean = self.analyzer.meanIntensity
        self._intens_mean = self.analyzer.meanIntensity
        return self._intens_mean

    @property
    def vrad_corr(self):
        if None is self._vrad_corr:
            self._vrad_corr = self.analyzer.rv_corr(relative_depth=self.r_depth)
        return self._vrad_corr

    @property
    def vspan(self):
        if self._vspan is None:
            self._vspan = self.analyzer.vspan(upper=self.upper, lower=self.lower)
        return self._vspan

    @property
    def vrad_bis(self):
        if None is self._vrad_bis:
            self._vrad_bis = self.analyzer.vrad_bis(extension=self.extension)
        return self._vrad_bis

    @property
    def vrad_skew(self):
        if None is self._vrad_skew:
            self._vrad_skew = self.analyzer.rv_skew(self.r_depth)
        return self._vrad_skew

    @property
    def vrad_std(self):
        if None is self._vrad_std:
            self._vrad_std = self.analyzer.rv_std(self.r_depth)
        return self._vrad_std

    @property
    def ls_vrad_skew(self):
        if None is self._ls_vrad_skew:
            mn = np.mean(self.vrad_skew)
            self._ls_vrad_skew = lombscargle(self.time, self.vrad_skew - mn, self.freq)
        return self._ls_vrad_skew

    @property
    def ls_vrad_mean(self):
        if None is self._ls_vrad_mean:
            mn = np.mean(self.vrad_mean)
            self._ls_vrad_mean = lombscargle(self.time, self.vrad_mean - mn, self.freq)
        return self._ls_vrad_mean

    @property
    def ls_vrad_corr(self):
        if None is self._ls_vrad_corr:
            mn = np.mean(self.vrad_corr)
            self._ls_vrad_corr = lombscargle(self.time, self.vrad_corr - mn, self.freq)
        return self._ls_vrad_corr

    @property
    def ls_vspan(self):
        if None is self._ls_vspan:
            mn = np.mean(self.vspan)
            self._ls_vspan = lombscargle(self.time, self.vspan - mn, self.freq)
        return self._ls_vspan

    @property
    def ls_vrad_bis(self):
        if None is self._ls_vrad_bis:
            mn = np.mean(self.vrad_bis)
            self._ls_vrad_bis = lombscargle(self.time, self.vrad_bis - mn, self.freq)
        return self._ls_vrad_bis

    @property
    def bisector(self):
        if None is self._bisector:
            self._bisector, self._depth = self.analyzer.bisector (upper=self.d0, lower=self.d1)
        return self._bisector

    @property
    def depth(self):
        if None is self._depth:
            self._bisector, self._depth = self.analyzer.bisector (upper=self.d0, lower=self.d1)
        return self._depth

    @property
    def window(self):
        if None is self._window:
            self._window = lombscargle(self.time, np.ones(len(self.time)), self.freq)
        return self._window

    @property
    def eqwidth(self):
        if None is self._eqwidth:
            self._eqwidth = self.analyzer.eqwidth()
        return self._eqwidth

    @property
    def ls_eqwidth(self):
        if None is self._ls_eqwidth:
            mn = np.mean(self.eqwidth)
            self._ls_eqwidth = lombscargle(self.time, self.eqwidth - mn, self.freq)
        return self._ls_eqwidth


    def vrad_mean_vrad_corr(self):
        name = "vrad_mean_vrad_corr"  # name of file

        plt.figure()
#        plt.title(name)
        plt.title(' ')
        plt.plot(self.vrad_mean, self.vrad_corr, 'o')
        plt.xlabel('vrad (correlation of profiles)  (m/s)')
        plt.ylabel('vrad (first moment) (m/s)')
        plt.savefig(name + self.format)

        self.comment += "\n----------------\n"
        self.comment += name + ":\n"
        self.comment += ""

    def ts_vrad(self):
        name = "ts_vrad_mean_corr"

        plt.figure()
        plt.title(name)
        plt.title('vrad_mean')
        plt.plot(self.analyzer.time, self.vrad_mean, '.')
        plt.plot(self.analyzer.time, self.vrad_corr, 'o')
        plt.ylabel('vrad')

        plt.savefig(name + self.format)

    def intens(self):
        name = "intens"

        plt.figure()
        mn = self.intens_mean
        mn = mn-1.
        mn /= mn.min()
        mn = 1.-mn
        plt.title('intensity profile')
        plt.axis([-60.,40.,0.,1.])
        plt.plot(self.velocity, mn)
        plt.xlabel('vrad (first moment) [km/s]')
        plt.ylabel('scaled intensity')
        plt.savefig(name + self.format)

    def intens_all(self):
        name = 'intens_all'

        plt.figure()

        mnint = self.intens_mean
        mnint = mnint - 1.
        mnint /= mnint.min()
        mnint = 1. - mnint
        for n in range(self.nspec):
            mn = self.inte[n, :]
            mn = mn-1.
            mn /= mn.min()
            mn = 1.-mn
            plt.axis([-60.,40.,-0.1,0.1])
            plt.plot(self.velocity, mn - mnint)
        plt.title('intensity profile')
        plt.xlabel('vrad (first moment) [km/s]')
        plt.ylabel('scaled intensity')
        plt.savefig(name + self.format)


    def bisector_time(self):
        name = "bisector_time"
        plt.figure()
        gs=gridspec.GridSpec(1, 2, width_ratios=[3, 1])
        plt.subplot(gs[0])
        plt.axis([-14,-12.75,0.,1.])
        plt.title('bisector variations')
        plt.xlabel('velocity (km/s)')
        plt.ylabel('Profile depth')
        n, d = self.bisector.shape
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

    def bisector_width(self):
        name = "bisector_width"

        plt.figure()
        plt.title('bisector standard deviation')
        n, d = self.bisector.shape
        s = self.bisector.std(axis=0)
        plt.plot(s, -self.depth, '-b', linewidth=4)
        plt.ylabel('Profile depth')
        plt.xlabel('std of bisector (km/s)')
        plt.ylim(-1, 0)
        plt.xlim(0, np.max(s) * 1.1)

        plt.savefig(name + self.format)

    def vrad_mean_vspan(self):
        name = "vrad_mean_vspan"

        x = self.vrad_mean
        y = self.vspan
        xy = np.vstack([x,y])
        z =gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

#        fig, ax = plt.subplots()
#        ax.scatter(x, y, c=z, s=20)
        colours = np.zeros( (len(z),3) )
        norm = Normalize( vmin=z.min(), vmax=z.max() )
        colours = [cm.ScalarMappable( norm=norm, cmap='rainbow').to_rgba( val ) for val in z]

        plt.figure()
        plt.title(name)
        plt.scatter( x, y, color=colours )
#


#        plt.axis([-13.06,-12.96,0.1,0.6])
        plt.title('Correlation vspan -- radial velocity')
#        plt.plot(self.vrad_mean, self.vspan, 'o')
        plt.xlabel('radial velocity (first moment) (km/s)')
        plt.ylabel('vspan (km/s)')
        plt.xlim(-13.1, -13.01)
        plt.ylim(0.2, 0.6)
        plt.savefig(name + self.format)


    def vrad_corr_vspan(self):
        name = "vrad_corr_vspan"


#        print(self.vrad_corr)
        x = self.vrad_corr
        y = self.vspan
        xy = np.vstack([x,y])
        z =gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

#        fig, ax = plt.subplots()
#        ax.scatter(x, y, c=z, s=20)
        colours = np.zeros( (len(z),3) )
        norm = Normalize( vmin=z.min(), vmax=z.max() )
        colours = [cm.ScalarMappable( norm=norm, cmap='rainbow').to_rgba( val ) for val in z]

        plt.figure()
        plt.title(name)
        plt.scatter( x, y, color=colours )
##        plt.plot(self.vrad_corr, self.vspan, 'o', color='#ED7F10', ms=12, alpha=0.7)
        plt.xlabel('vrad (m/s)')
        plt.ylabel('vspan (m/s)')
        plt.title('vspan as a function of radial velocity')
        plt.xlim(-0.17, 0.17)
        plt.ylim(0.15, 0.55)

        plt.savefig(name + self.format)

    def vrad_mean_skew(self):
        name = "vrad_mean_skew"
        plt.figure()
        plt.title(name)
        plt.plot(self.vrad_mean, self.vrad_skew, 'o')
        plt.savefig(name + self.format)

    def vrad_mean_std(self):
        name = "vrad_mean_std"
        plt.figure()
        plt.title(name)
        plt.plot(self.vrad_mean, self.vrad_std, 'o')
        plt.savefig(name + self.format)


    def vrad_corr_skew(self):
        name = "vrad_corr_skew"
        plt.figure()
        plt.title(name)
        plt.plot(self.vrad_corr, self.vrad_skew, 'o')
        plt.savefig(name + self.format)

    def skew_vspan(self):
        name = "skew_vspan"
        plt.figure()
        plt.title(name)
        plt.plot(self.vrad_skew, self.vspan, 'o')
        plt.savefig(name + self.format)

    def ts_skew(self):
        name = "ts_skew"
        plt.figure()
        plt.title(name)
        plt.plot(self.analyzer.time, self.vrad_skew, '-')
        plt.savefig(name + self.format)

    def ts_vspan(self):
        name = "ts_vspan"
        plt.figure()
        plt.title(name)
        plt.plot(self.analyzer.time, self.vspan, '-')
        plt.savefig(name + self.format)

    def ts_eqwidth(self):
        name = "ts_eqwidth"
        plt.figure()
        plt.title(name)
        plt.plot(self.analyzer.time, self.eqwidth, '-')
        plt.savefig(name + self.format)

    def _plt_ls(self, amp, box=[], bars=True):
        ll = [self.cpdl, self.cpd2l, self.cpd3l, self.cpd4l, self.cpd5l, self.cpd6l, self.cpd7l, self.cpd8l, self.cpd9l, self.cpd10l]
        hh = [self.cpdh, self.cpd2h, self.cpd3h, self.cpd4h, self.cpd5h, self.cpd6h, self.cpd7h, self.cpd8h, self.cpd9h, self.cpd10h]

        if bars:
            for l, h in zip (ll, hh):
                plt.gca().add_patch(plt.Rectangle((l,0), h-l, 0.3, fc='0.85', color='0.85'))
                plt.gca().add_patch(plt.Rectangle((l-1,0), h-l, 0.15, fc='0.7', color='0.7'))
                plt.gca().add_patch(plt.Rectangle((l+1,0), h-l, 0.15, fc='0.7', color='0.7'))
                plt.gca().add_patch(plt.Rectangle((l-2,0), h-l, 0.02, fc='0.5', color='0.5'))
                plt.gca().add_patch(plt.Rectangle((l+2,0), h-l, 0.02, fc='0.5', color='0.5'))

            for n in range(1,5):
                plt.vlines([n*self.cpdnew],0,0.3)
                plt.vlines([n*self.cpdnew-1],0,0.15)
                plt.vlines([n*self.cpdnew+1],0,0.15)
                plt.vlines([n*self.cpdnew -2 ],0,0.02)
                plt.vlines([n*self.cpdnew +2],0,0.02)

        plt.xlim(box[0],box[1])
        plt.ylim(box[2],box[3])
        plt.plot(self.cycles_per_day, amp, '-k',linewidth=2, color = 'r')
        plt.minorticks_on()

    def _plt_ls_vr(self, data):
        amp = np.abs(data)
#        amp /= np.max(amp)
        ll = [self.cpdl, self.cpd2l, self.cpd3l, self.cpd4l, self.cpd5l, self.cpd6l, self.cpd7l, self.cpd8l, self.cpd9l]
        hh = [self.cpdh, self.cpd2h, self.cpd3h, self.cpd4h, self.cpd5h, self.cpd6h, self.cpd7h, self.cpd8h, self.cpd9h]
        for l, h in zip (ll, hh):
            plt.gca().add_patch(plt.Rectangle((l,0), h-l, 10., fc='0.85', color='0.85'))
            plt.gca().add_patch(plt.Rectangle((l-1,0), h-l, 5., fc='0.7', color='0.7'))
            plt.gca().add_patch(plt.Rectangle((l+1,0), h-l, 5., fc='0.7', color='0.7'))
            plt.gca().add_patch(plt.Rectangle((l-2,0), h-l, 1., fc='0.5', color='0.5'))
            plt.gca().add_patch(plt.Rectangle((l+2,0), h-l, 1., fc='0.5', color='0.5'))

        for n in range(1,5):
            plt.vlines([n*self.cpdnew],0,10.)
            plt.vlines([n*self.cpdnew-1],0,5.)
            plt.vlines([n*self.cpdnew+1],0,5.)
            plt.vlines([n*self.cpdnew -2 ],0,1.)
            plt.vlines([n*self.cpdnew +2],0,1.)

        maxa = max(amp)
        plt.axis([0,15.,0,maxa])
        plt.plot(self.cycles_per_day, amp, '-k',linewidth=2, color = 'r')
        plt.minorticks_on()

    def _plt_ls_all(self, data, data2, data3, box=None):
        amp = np.abs(data)
        amp2 = np.abs(data2)
        amp3 = np.abs(data3)
        amp /= np.max(amp)
        amp2 /= np.max(amp2)
        amp3 /= np.max(amp3)

        ll = [self.cpdl, self.cpd2l, self.cpd3l, self.cpd4l, self.cpd5l, self.cpd6l, self.cpd7l, self.cpd8l, self.cpd9l]
        hh = [self.cpdh, self.cpd2h, self.cpd3h, self.cpd4h, self.cpd5h, self.cpd6h, self.cpd7h, self.cpd8h, self.cpd9h]
        for l, h in zip (ll, hh):
            plt.gca().add_patch(plt.Rectangle((l,0), h-l, 10., fc='0.85', color='0.85'))
            plt.gca().add_patch(plt.Rectangle((l-1,0), h-l, 5., fc='0.7', color='0.7'))
            plt.gca().add_patch(plt.Rectangle((l+1,0), h-l, 5., fc='0.7', color='0.7'))
            plt.gca().add_patch(plt.Rectangle((l-2,0), h-l, 1., fc='0.5', color='0.5'))
            plt.gca().add_patch(plt.Rectangle((l+2,0), h-l, 1., fc='0.5', color='0.5'))

        for n in range(1,5):
            plt.vlines([n*self.cpdnew],0,10.)
            plt.vlines([n*self.cpdnew-1],0,5.)
            plt.vlines([n*self.cpdnew+1],0,5.)
            plt.vlines([n*self.cpdnew -2 ],0,1.)
            plt.vlines([n*self.cpdnew +2],0,1.)

        maxa = max(amp)
        plt.axis([0,15.,0,maxa])
        plt.plot(self.cycles_per_day, amp, '-k',linewidth=2, color = 'r')
        plt.plot(self.cycles_per_day, amp2, '-k',linewidth=2, color = 'b')
        plt.plot(self.cycles_per_day, amp3, '-k',linewidth=2, color = 'g')

        plt.minorticks_on()


    def _plt_ls_vb(self, data):
        amp = np.abs(data)
#        amp /= np.max(amp)
        ll = [self.cpdl, self.cpd2l, self.cpd3l, self.cpd4l, self.cpd5l, self.cpd6l, self.cpd7l, self.cpd8l, self.cpd9l]
        hh = [self.cpdh, self.cpd2h, self.cpd3h, self.cpd4h, self.cpd5h, self.cpd6h, self.cpd7h, self.cpd8h, self.cpd9h]
        for l, h in zip (ll, hh):
            plt.gca().add_patch(plt.Rectangle((l,0), h-l, 0.15, fc='0.85', color='0.85'))
            plt.gca().add_patch(plt.Rectangle((l-1,0), h-l, 0.075, fc='0.7', color='0.7'))
            plt.gca().add_patch(plt.Rectangle((l+1,0), h-l, 0.075, fc='0.7', color='0.7'))
            plt.gca().add_patch(plt.Rectangle((l-2,0), h-l, 0.01, fc='0.5', color='0.5'))
            plt.gca().add_patch(plt.Rectangle((l+2,0), h-l, 0.01, fc='0.5', color='0.5'))

        for n in range(1,5):
            plt.vlines([n*self.cpdnew],0,0.15)
            plt.vlines([n*self.cpdnew-1],0,0.075)
            plt.vlines([n*self.cpdnew+1],0,0.075)
            plt.vlines([n*self.cpdnew -2 ],0,0.01)
            plt.vlines([n*self.cpdnew +2],0,0.01)

        maxa = max(amp)
        plt.axis([0,15.,0,maxa])
        plt.plot(self.cycles_per_day, amp, '-k',linewidth=2, color = 'r')
        plt.minorticks_on()




    def _plt_ls_pure(self,data):
        amp = np.abs(data)
        amp /= np.max(amp)

        maxa = max(amp)
#tata
        plt.axis([0,15.,0,maxa])
        plt.plot(self.cycles_per_day, amp, '-k',linewidth=2, color = 'r')
        plt.minorticks_on()



    def ls_spec_vrad_skew(self):
        name = 'ls_spec_vrad_skew'
        plt.figure(figsize=(10,4))
        plt.title(name)
        self._plt_ls(self.ls_vrad_skew)
        plt.savefig(name + self.format)

    def ls_spec_vrad_mean(self):
        name = 'ls_spec_vrad_mean'

        data = np.abs(self.ls_vrad_mean)
        data /= np.max(data)
        plt.figure(figsize=(10,4))
        plt.title('Lomb Scargle periodogram of vrad (first moment)')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')
#        self._plt_ls_vr(1000*self.ls_vrad_mean)

        maxa = np.max(data)
        box = [0,15,0,1]
        self._plt_ls(data, box=box)
        plt.savefig(name + self.format)

        plt.figure(figsize=(10,4))
        plt.title('Lomb Scargle periodogram of vrad (first moment)')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')
#        self._plt_ls_vr(1000*self.ls_vrad_mean)

        box = [15,50,0,0.15]
        self._plt_ls(data,box=box, bars=False)
        plt.savefig(name + "_2"+self.format)

    def ls_spec_vrad_bis(self):
        name = 'ls_spec_vrad_bis'

        data = np.abs(self.ls_vrad_bis)
        data /= np.max(data)
        plt.figure(figsize=(10,4))
        plt.title('Lomb Scargle periodogram of vrad (bisector)')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')
#        self._plt_ls_vr(1000*self.ls_vrad_mean)

        maxa = np.max(data)
        box = [0,15,0,1]
        self._plt_ls(data, box=box)
        plt.savefig(name + self.format)

        plt.figure(figsize=(10,4))
        plt.title('Lomb Scargle periodogram of vrad (bisector)')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')
#        self._plt_ls_vr(1000*self.ls_vrad_mean)

        box = [15,50,0,0.15]
        self._plt_ls(data,box=box, bars=False)
        plt.savefig(name + "_2"+self.format)

    def ls_spec_eqwidth(self):
        name = 'ls_spec_eqwidth'

        data = np.abs(self.ls_eqwidth)
        data /= np.max(data)
        plt.figure(figsize=(10,4))
        plt.title('Lomb Scargle periodogram of eqwidth')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')

        maxa = np.max(data)
        box = [0,15,0,1]
        self._plt_ls(data, box=box)
        plt.savefig(name + self.format)

        plt.figure(figsize=(10,4))
        plt.title('Lomb Scargle periodogram of eqwidth')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')
        box = [15,50,0,0.15]
        self._plt_ls(data,box=box, bars=False)
        plt.savefig(name + "_2"+self.format)

    def ls_spec_all3(self):
        name = 'ls_spec_all3'

        data = np.abs(self.ls_vspan)
        data2 = np.abs(self.ls_vrad_mean)
        data3 = np.abs(self.ls_vrad_bis)
        data /= np.max(data)
        plt.figure(figsize=(10,4))
        plt.title('3 Lomb Scargle periodograms')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')

        maxa = np.max(data)
        box = [0,15,0,1]
        self._plt_ls_all(data,data2,data3, box=box)
        plt.savefig(name + self.format)

        plt.figure(figsize=(10,4))
        plt.title('3 Lomb Scargle periodograms')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')

        box = [15,50,0,0.15]
        self._plt_ls(data,box=box, bars=False)
        plt.savefig(name + "_2"+self.format)



    def ls_spec_vrad_corr(self):
        name = 'ls_spec_vrad_corr'
        plt.figure()
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative amplitude')
        plt.title(name)
        self._plt_ls(self.ls_vrad_corr)
        plt.savefig(name + self.format)

    def ls_window(self):
        name = 'window_function'
        plt.figure(figsize=(6,4))
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')
        plt.title('window function')
        self._plt_ls_pure(self.window)
        plt.axis([0.,7.,0.,1.])
        plt.savefig(name + self.format)

    def ls_spec_vspan(self):
        name = 'ls_spec_vspan'

        data = np.abs(self.ls_vspan)
        data /= np.max(data)
        plt.figure(figsize=(10,4))
        plt.title('Lomb Scargle periodogram of vspan')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')
#        self._plt_ls_vr(1000*self.ls_vrad_mean)

        maxa = np.max(data)
        box = [0,15,0,1]
        self._plt_ls(data, box=box)
#[left, bottom, width, height]
        a = plt.axes([.65, .6, .2, .2], facecolor='white')
        plt.title('window function')
        wf = self.window
        wf /= wf.max()
        plt.plot(self.cycles_per_day, wf, 'b-')
        plt.xlim(0,7)
        plt.ylim(0,1)

        plt.savefig(name + self.format)

        plt.figure(figsize=(10,4))
        plt.title('Lomb Scargle periodogram of vspan')
        plt.xlabel('Frequency (c/d)')
        plt.ylabel('Relative power spectral density')
#        self._plt_ls_vr(1000*self.ls_vrad_mean)

        box = [15,50,0,1.0]
        self._plt_ls(data,box=box, bars=False)
        plt.savefig(name + "_2"+self.format)

    def _bayes_freq(self, value, name):
        name = "bayes_freq_" + name
        plt.figure()
        plt.title(name)
        plt.plot(self.cycles_per_day, self.analyzer.PosteriorFreq(self.freq, value, nharm=1))
        plt.vlines(self.cpd, 0, 1)

        plt.savefig(name + self.format)

    def bayes_freq_vrad_mean(self):
        self._bayes_freq(self.vrad_mean, name='mean')

    def moving_peaks(self):
        name = 'moving_peaks'
        plt.figure(figsize=(10,4))
        val = 1-self.inte
        val /= val.sum(axis=1)[:,np.newaxis]
        val -= np.median(val,axis=0)
        #val /= val.max(axis=1)[:,np.newaxis]
        #Sophie 2018:2.458331340824170038e+06
        #Narval 2018:2.458332337050000206e+06
#        tdiff=2.458332337050000206e+06-2.458331340824170038e+06
        tmp, bins = self.analyzer.spectrum_matrix(self.time, val, period=self.rotperiod, nphase=32)

        plt.imshow(np.abs(tmp[::-1,:])**0.5, interpolation='none', cmap=plt.cm.Blues)
        plt.plot(tmp)
        plt.savefig(name + self.format)

    def ew_noise_corr(self):
        name = 'moving_peaks_eq_width'
        #return #TODO REMOVE
        vul = 1.-self.inte
        eqwidth = self.analyzer.eqwidth()
        signois = self.analyzer.meansignoise()

        for na, nightlist in zip(['s1','s2','s3','s4','s5', 's12345','s123','s45'], [[0], [1], [2], [3], [4], [0,1,2,3,4], [0,1,2],[3,4]]):
            VV =[]
            TT =[]
            plt.title('night'+na)
            for night in nightlist:
                print("night = ", night)
                I = self.analyzer.list_index[night]
                print(night, len(I))
                TT += list(self.time[I])
                pp = np.poly1d(np.polyfit(signois[I], eqwidth[I], deg=3))

                plt.figure()
                plt.plot(signois[I],eqwidth[I],'o')
                maxs = np.max(signois[I])
                mins = np.min(signois[I])
                newx = np.linspace(mins,maxs,1000)
                plt.plot(newx,pp(newx))
                plt.show()



    def moving_peaks_time(self):
        name = 'moving_peaks_time'
        
        vul = 1.-self.inte
        eqwidth = self.analyzer.eqwidth()
        signois = self.analyzer.meansignoise()


#        for na, nightlist in zip(['s:1', 's:2', 's:3', 's:123456'], [[0], [1], [2], [0,1,2,3,4,5]]):

        minmax=[]
        

        for na, nightlist in zip(['s1','s2','s3','s4','s5','s6','s7', 's1234567','s123','s4567'], [[0], [1], [2], [3], [4], [5], [6], [0,1,2,3,4,5,6], [0,1,2],[3,4,5,6]]):
        #for na, nightlist in zip(['s1','s2','s3','s4'], [[0], [1], [2], [3]]):
        #for na, nightlist in zip(['s1'], [[0]]):
        #for na, nightlist in zip(['s:12','s:3','s:123'], [[0,1],[2],[0,1,2]]):
            VV =[]
            TT =[]

            mintime=1e8
            maxtime=0
            plt.figure(figsize=(6,10))
            plt.title('night'+na)
            for night in nightlist:
                print("night = ", night)
                I = self.analyzer.list_index[night]
                print(night, len(I))
                TT += list(self.time[I])

                mintime=min(np.min(self.time[I]), mintime)
                maxtime=max(np.max(self.time[I]), maxtime)

                pp = np.poly1d(np.polyfit(signois[I], eqwidth[I], deg=3))
                fac = pp(signois[I])[:,np.newaxis]
                temp=vul[I]/fac
                temp-=np.median(temp,axis=0)
                VV += list(temp)

            TT = [tt - mintime for tt in TT]
                
            val = np.row_stack(VV)

            nphase=128
            tmp, bins, mask = self.analyzer.spectrum_matrix_full(
                TT, val, period=maxtime-mintime, nphase=nphase, method=np.median)
            #print("minmax:",np.min(tmp),np.max(tmp))
            minmax.append([np.min(tmp),np.max(tmp)])
            
            #plt.contourf(self.velocity, bins, tmp, cmap=plt.cm.Reds)


            #tmp = np.sign(tmp) * np.abs(tmp)**0.7

            #ax = plt.axes([self.velocity[0],self.velocity[-1],0,1], frameon=False)
            #ax.set_axis_off()
            #ax.set_xlim(self.velocity[0],self.velocity[-1])
            #ax.set_ylim(0,1)
            v0=-13.01

            #plt.imshow(tmp, cmap=plt.cm.gray_r, aspect='auto',interpolation='none', origin='lower',extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1])
            plt.imshow(tmp, cmap=plt.cm.gray_r, aspect='auto',interpolation='None', origin='lower',extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],vmin=-0.001,vmax=0.0015)
            #plt.plot(tmp)
            plt.xticks([])
            rv = np.array([-20, -10, 0, 10, 20 ])
            st = [ str(i) for i in rv]
            #plt.xticks(rv, st)
            #plt.yticks([])
            yt = [0,0.2,0.4,0.6,0.8, 1]
            #plt.yticks(yt,[str(t) for t in yt])
            #plt.vlines([-22, 0, 22], 0,1)
            #plt.vlines([-35, 0, 37], 0,1,linestyles='dashed')
            #plt.hlines([1./2], -33,32)
            plt.xlim(self.velocity[0]-v0, self.velocity[-1]-v0)
            plt.xlabel('velocity [km/s]')
            plt.ylabel('phase fraction of period]')

            plt.savefig(name + na +self.format)
            
        print("final minmax:",np.min(np.array(minmax)),np.max(np.array(minmax)))

    def moving_peaks_simple(self):
        nbins = 128 #256
        sa = self.analyzer
        time = sa.time

        meanprofile = np.median(sa.intensity, axis=0)

        F = np.zeros((sa.nvelocity, 2))
        F[:,0] = 1
        F[:,1] = meanprofile

        diff = np.zeros((sa.nobs, sa.nvelocity))

        for i in range(sa.nobs):
            pp = np.linalg.lstsq(F, sa.intensity[i], rcond=None)[0]
            diff[i,:] = sa.intensity[i] - (pp[0] + pp[1] * meanprofile)

        print(pp)


        # diff = sa.intensity - np.median(sa.intensity, axis=0)

        phasebins = np.linspace(0, self.rotperiod, nbins+1)
        #I = np.digitize(np.mod(sa.time, self.rotperiod), phasebins)
        
        res = np.zeros((nbins, sa.nvelocity))

        tt = np.mod(sa.time, self.rotperiod)
        for i, p1p2 in enumerate(zip(phasebins[0:-1], phasebins[1:])):
            p1, p2 = p1p2
            I = (p1 <= tt) & (tt < p2)
            if sum(I) >=4:
                res[i,:] = np.mean(diff[I], axis=0)
        

        minmax = np.max(np.abs(res))
        print(minmax)

        res /= minmax
        res = np.sign(res) * np.abs(res)**0.75

        plt.imshow(
            res, 
            cmap=plt.cm.gray_r, 
            aspect='auto',
            interpolation='None', 
            origin='lower',
            #extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],
            vmin=-0.8, vmax=0.8)








    def moving_peaks_signoise(self):
        name = 'moving_peaks_eq_width'
        #return #TODO REMOVE
        vul = 1.-self.inte
        eqwidth = self.analyzer.eqwidth()
        signois = self.analyzer.meansignoise()


#        for na, nightlist in zip(['s:1', 's:2', 's:3', 's:123456'], [[0], [1], [2], [0,1,2,3,4,5]]):

        minmax=[]
        

        #for na, nightlist in zip (['s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s123456','s123','s456','s789','all'], [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9],[0,1,2,3,4,5], [0,1,2],[3,4,5],[6,7,8], [0,1,2,3,4,5,6,7,8,9]]):
        #for na, nightlist in zip(['s1','s2','s3','s4'], [[0], [1], [2], [3]]):
        
        #2024:
        #for na, nightlist in zip(['s1','s2','s3','s4','s5', 's12345','s123','s45'], [[0], [1], [2], [3], [4], [0,1,2,3,4], [0,1,2],[3,4]]):
        
        #2018 TBL:
        #for na, nightlist in zip(['s1','s2','s3','s4','s5','s6','s7', 's1234567','s123','s4567'], [[0], [1], [2], [3], [4], [5], [6], [0,1,2,3,4,5,6], [0,1,2],[3,4,5,6]]):
        
        #2018 OHP:
        for na, nightlist in experiment['nights']:

            VV =[]
            TT =[]
            plt.figure(figsize=(6,10))
            plt.title('night'+na)
            for night in nightlist:
                print("night = ", night)
                I = self.analyzer.list_index[night]
                print(night, len(I))
                TT += list(self.time[I])
                pp = np.poly1d(np.polyfit(signois[I], eqwidth[I], deg=3))
                fac = pp(signois[I])[:,np.newaxis]
                temp=vul[I]/fac
                temp-=np.median(temp,axis=0)
                VV += list(temp)
                
            val = np.row_stack(VV)
            val -= np.median(val,axis=0)

            nphase=128
            tmp, bins, mask = self.analyzer.spectrum_matrix_full(TT, val, period=self.rotperiod, nphase=nphase, method=np.median)
            #print("minmax:",np.min(tmp),np.max(tmp))
            minmax.append([np.min(tmp),np.max(tmp)])
            
            #plt.contourf(self.velocity, bins, tmp, cmap=plt.cm.Reds)


            #tmp = np.sign(tmp) * np.abs(tmp)**0.7

            #ax = plt.axes([self.velocity[0],self.velocity[-1],0,1], frameon=False)
            #ax.set_axis_off()
            #ax.set_xlim(self.velocity[0],self.velocity[-1])
            #ax.set_ylim(0,1)
            v0=-13.01

            #plt.imshow(tmp, cmap=plt.cm.gray_r, aspect='auto',interpolation='none', origin='lower',extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1])
            plt.imshow(tmp, cmap=plt.cm.gray_r, aspect='auto',interpolation='bicubic', origin='lower',extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],vmin=-0.0003,vmax=0.0003)
            #plt.plot(tmp)
            plt.xticks([])
            rv = np.array([-20, -10, 0, 10, 20 ])
            st = [ str(i) for i in rv]
            plt.xticks(rv, st)
            plt.yticks([])
            yt = [0,0.2,0.4,0.6,0.8, 1]
            plt.yticks(yt,[str(t) for t in yt])
            plt.vlines([-22, 0, 22], 0,1)
            plt.vlines([-35, 0, 37], 0,1,linestyles='dashed')
            plt.hlines([1./2], -33,32)
            plt.xlim(self.velocity[0]-v0, self.velocity[-1]-v0)
            plt.xlabel('velocity [km/s]')
            plt.ylabel('phase fraction of period]')

            plt.savefig(name + na +self.format)
            
        print("final minmax:",np.min(np.array(minmax)),np.max(np.array(minmax)))

    def ew_noise_corr(self):
        name = 'moving_peaks_eq_width'
        #return #TODO REMOVE
        vul = 1.-self.inte
        eqwidth = self.analyzer.eqwidth()
        signois = self.analyzer.meansignoise()


#       for na, nightlist in zip(['s1','s2','s3','s4','s5','s6'], [[0], [1], [2], [3], [4], [5]]):
#        for na, nightlist in zip(['s1','s2','s3','s4','s5','s6','s7'], [[0], [1], [2], [3], [4], [5], [6]]):


        for na, nightlist in zip(['s1','s2','s3','s4','s5','s6','s7'], [[0], [1], [2], [3], [4], [5], [6]]):

            VV =[]
            TT =[]

            for night in nightlist:
                print("night = ", night)
                try:
                    I = self.analyzer.list_index[night]
                    print(night, len(I))
                except:
                    continue
                TT += list(self.time[I])
                pp = np.poly1d(np.polyfit(signois[I], eqwidth[I], deg=3))

                plt.figure()
                plt.title('night'+na)
                plt.plot(signois[I],eqwidth[I],'o')
                maxs = np.max(signois[I])
                mins = np.min(signois[I])
                newx = np.linspace(mins,maxs,1000)
                plt.plot(newx,pp(newx))


    def signoise_eqwidth(self):
        bins = [0,0.4,0.5, 0.55, 0.6,0.7]
        clrs = ['-r', '-g', '-b', '-k','-c']
        tt = np.digitize(np.mod(self.time,1), bins)

        for j in range(5):
            name = 'signoise_eqwidth' + str(j)
            plt.figure(figsize=(10,4))
            plt.title('night '+str(j+1))
            for i in range(1,6):
                II = self.analyzer.list_index[j]
                I = np.where(tt[II]==i)[0]
                plt.plot(self.analyzer.eqwidth()[II][I], self.analyzer.meansignoise()[II][I], clrs[i-1])
            plt.savefig(name + self.format)

    def estrotentropy(self):
        """
        estimating the rotation using local generalized entroypy
        """

        def ecart(x, axis):
            #return np.std(x,axis=axis)
            return np.percentile(x,95, axis=axis)-np.percentile(x, 5, axis=axis)

        plt.figure()
        name = 'spreadphasemap'
        nrot = 256 # number of rotperiods
        p0, p1 = 0.55, 0.8
        pr = np.linspace(p0,p1,nrot)

        vul = 1.-self.inte
        eqwidth = self.analyzer.eqwidth()
        signois = self.analyzer.meansignoise()

        nnoise = 5
        res = np.zeros(nrot)
        noiseres = np.zeros((nnoise, nrot))
        for i,p in enumerate(pr):
            VV =[]
            TT =[]
            for night in [0,1,2,3,4]:
                I = self.analyzer.list_index[night]
                TT += list(self.time[I])
                pp = np.poly1d(np.polyfit(signois[I], eqwidth[I], deg=3))
                fac = pp(signois[I])[:,np.newaxis]
                VV += list(vul[I]/fac)
            val = np.row_stack(VV)
            val -= np.median(val,axis=0)
            tmp, bins, mask = self.analyzer.spectrum_matrix_full(TT, val,
                            period=p, nphase=128, method=ecart)
            res[i]=np.mean(tmp[mask])
            for j in range(nnoise):
                tmp, bins, mask = self.analyzer.spectrum_matrix_full(
                         np.random.permutation(TT), val,
                            period=p, nphase=128, method=ecart)
                noiseres[j,i]=np.mean(tmp[mask])
        plt.title('spread of velocity phase map')
        plt.plot(pr,res)
        for j in range(nnoise):
            plt.plot(pr, noiseres[j,:], '.g')
        #adding random shuffles

        plt.xlabel('rotation period [days]')
        plt.ylabel('quantile spread in arbitrary units')
        plt.yticks([])
        plt.vlines(self.rotperiod, [0], [1], colors=['red'])
        plt.ylim(0.6*res.max(), 1.2*res.max())
        plt.yticks([0.9*res.max(), 1.0*res.max()],['0.9','1.0'])
        plt.xlim(p0,p1)
        plt.savefig(name+self.format)

if __name__ == '__main__':
    myPics = Pictures(DATAFILE)
    #myPics.vrad_mean_vrad_corr()
    #myPics.ts_vrad()
    #myPics.ts_eqwidth()
#    myPics.intens()
#    myPics.intens_all()
#    myPics.vrad_mean_vspan()
#    myPics.vrad_corr_vspan()
#    myPics.vrad_mean_skew()
#    myPics.vrad_mean_std()
#    myPics.vrad_corr_skew()
#    myPics.skew_vspan()
#    myPics.ts_skew()
#    myPics.ts_vspan()
#    myPics.ls_spec_vrad_skew()
#    myPics.ls_spec_vrad_mean()
#    myPics.ls_spec_all3()
#    myPics.ls_spec_vrad_corr()
#    myPics.ls_spec_vrad_bis()
     # myPics.ls_spec_vspan()
#    myPics.ls_spec_eqwidth()
###    myPics.bisector_time()
#    myPics.bisector_width()
#    myPics.ls_window()
#    alldata = [self.time, self.inte, self.vrad_mean, self.vrad_corr, self.vspan, self.vrad_skew, self.vrad_std]
 #   myPics.bayes_freq_vrad_mean()
    # myPics.moving_peaks_signoise()
    myPics.moving_peaks_simple()
   # myPics.moving_peaks_time()
###    myPics.estrotentropy()

    #myPics.saveData("time_vrad_mean.dat", [6142.+myPics.time, myPics.vrad_mean])  # first column
    #myPics.saveData("time_vrad_corr.dat", [myPics.time, myPics.vrad_corr])  # first column
    #myPics.saveData("time_vrad_bis.dat", [6142.+myPics.time, myPics.vrad_bis])  # first column
    #myPics.saveData("time_vspan.dat", [6142.+myPics.time, myPics.vspan])  # first column
    #myPics.saveData("time_skew.dat", [myPics.time, myPics.vrad_skew])  # first column




    # myPics.ew_noise_corr()

    plt.show()
