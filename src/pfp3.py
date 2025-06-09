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
from mpl_toolkits.basemap import Basemap
from spectralutil_old import *

radial_velocity = -13.9 #km/s
rotation_vsini = 21.6 #km/s
full_prof_add = 10 #km/s

# rotation period as of alina et al
# rotper = 0.678 #d upper_errorbar = 0.0042, lower_errorbar = - 0.0042

# rotation period as of 2012 and 2018 FELIX calculation
#rotper = 0.6705 #d error bar 0.0019

# rotation period FELIX for 2012 data
# rotper = 0.6679 # error bar 0.0016

# rotation period FELIX for 2018 data
# rotper = 0.6730 # error bar 0.0001

# rotation period fixed for comparision with PPetit
rotper = 0.678


prof_beg = radial_velocity-rotation_vsini #km/s
prof_end = radial_velocity+rotation_vsini #km/s
full_prof_beg =  prof_beg-full_prof_add #km/s
full_prof_end =  prof_end+full_prof_add#km/s


PROGDIR = os.path.dirname(__file__)
DATADIR = os.path.join(PROGDIR, '../data')
DATADIR = os.environ.get('DATADIR', DATADIR)

DATAFILE7 = os.path.join(DATADIR, 'duofilematrix.dat')
DATAFILE8 = os.path.join(DATADIR, 'Vega_9500.40.03-10.dat')  # lines between 0.3 and 1.0 depth, medium
DATAFILE9 = os.path.join(DATADIR, 'Vega_tbl10.dat')  # lines between 0.3 and 1.0 depth, medium
#DATAFILE10 = os.path.join(DATADIR, 'Vega_2018_0310.dat')  # Sophie lines between 0.3 and 1.0 depth, medium
DATAFILE10 = os.path.join(DATADIR, 'Vega_2018_SOPHIE_maskvega_folsom.clean.dat')
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


# OLD DATA REDUCTION Sophie 2018 Sophie_pipeline LSD (JFD) mask_0310
VEGA_2018_SOPHIE_LSDJFD = Experiment(
    DATAFILE = os.path.join(DATADIR, 'Vega_2018_0310.dat'), # Sophie lines between 0.3 and 1.0 depth, medium
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days Petit et al 2022
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5','s6', 's12','s34','s56'],
        [[0], [1], [2], [3], [4], [5], [0,1], [2,3],[4,5]]),
    lamfilter = -0.465
)

# Sophie 2012 Sophie_pipeline LSDpy 2.0km/s mask_folsom
VEGA_2012_SOPHIE_LSDPY = Experiment(
    DATAFILE = os.path.join(DATADIR, 'Vega_2012_SOPHIE_maskvega_folsom.clean.1.7'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5', 's12','s34','2012 SOPHIE (1,2,3,4,5)'],
        [[0], [1], [2], [3], [4], [0,1], [2,3],[0,1,2,3,4]]),
    lamfilter = 1000
    )


# Sophie 2012 Sophie_pipeline LSD Donati
VEGA_2012_SOPHIE_LSD = Experiment(
    DATAFILE = os.path.join(DATADIR, 'Vega_Sophie_2012_9500.40.03-10.dat'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5', 's12','s34','2012 SOPHIE (1,2,3,4)'],
        [[0], [1], [2], [3], [4], [0,1], [2,3],[0,1,2,3]]),
        #['s12'],
        #[[0,1]]),
    lamfilter = 1000
    )


# Sophie 2018 Sophie_pipeline LSDpy 2.0km/s mask_folsom
VEGA_2018_SOPHIE_LSDPY = Experiment(
DATAFILE = os.path.join(DATADIR, 'Vega_2018_SOPHIE_maskvega_folsom.clean.dat.1.7'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5', 's12','s34', '2018 SOPHIE (1,2,3,4,5,6)'],
        [[0], [1], [2], [3], [4], [0,1], [2,3], [0,1,2,3,4,5]]),
    lamfilter = -0.194
    )


# OLD DATA REDUCTION Narval 2018 LE LSD (JFD) mask_0310
VEGA_2018_NARVAL_LE = Experiment(
    DATAFILE = os.path.join(DATADIR, 'Vega_Narval_2018_031.dat'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,  
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5','s6', 's7', '2018 NARVAL (1,2,3,4,5,6,7)','s123','s456'],
        [[0], [1], [2], [3], [4], [5], [6], [0,1,2,3,4,5,6], [0,1,2],[3,4,5]]),
    lamfilter = -0.425
    )

# Narval 2018 NEXTRA LSDpy 2.Okm/s mask_folsom
VEGA_2018_NARVAL_NEXTRA = Experiment(
DATAFILE = os.path.join(DATADIR, 'Vega_2018_maskvega_folsom.clean'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5', 's12','s34', '2018 NARVAL (1,2,3,4,5,6)'],
        [[0], [1], [2], [3], [4], [0,1], [2,3], [0,1,2,3,4,5]]),
    lamfilter = 1000 #-0.174
    )


# Narval 2018 NEXTRA LSDpy 1.7km/s mask_folsom, Quantile 1
VEGA_2018_NARVAL_NEXTRA_KM17 = Experiment(
    DATAFILE = os.path.join(DATADIR, 'Vega_2018_maskvega_folsom.clean_1.7.clean'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,  
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5','s6', 's7', 's12','s34','s56','2018 NARVAL (1,2,3,4,5,6)'],
        [[0], [1], [2], [3], [4], [5], [6], [0,1], [2,3],[4,5],[0,1,2,3,4,5]]),
    lamfilter = 1000
    )



# Neo-Narval 2023 NEXTRA LSDpy 2.0km/s mask_folsom
VEGA_2023_NEONARVAL_NEXTRA = Experiment(
    DATAFILE = os.path.join(DATADIR, 'Vega_2023_maskvega_folsom.clean'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = list(zip(
        ['NEO-NARVAL 2023 (4,6,7,10)'],
        [[3,5,6,9]])),
        #['s1','s2','s3','s4','s5','s6', 's7', 's8', 's9', 's10', 's12','s34','s56','s78','s910','s46710'],
        #[[0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [0,1], [2,3], [4,5] , [6,7], [8,9],[3,5,6,9]])),
    lamfilter = -0.163
    )

# Neo-Narval 2023 NEXTRA LSDpy 2.0km/s mask_folsom
VEGA_2023_NEONARVAL_NEXTRA_SELECT = Experiment(
    DATAFILE = os.path.join(DATADIR, 'Vega_2023_maskvega_folsom.clean_31_6_7_12.clean'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4'],
        [[0], [1], [2], [3]]),
    lamfilter = 1000
    )


# Neo-Narval 2023 NEXTRA LSDpy 2.0km/s mask_folsom, quantile 0.9 selection in filematrix
VEGA_2023_NEONARVAL_NEXTRA_Q09 = Experiment(
    DATAFILE = os.path.join(DATADIR,'Vega_2023_maskvega_folsom.clean_0.9.dat'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5','s6', 's7', 's8', 's9', 's10', 's12','s34','s56','s78','s910','s46710'],
        [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [0,1], [2,3], [4,5] , [6,7], [8,9],[3,5,6,9]]),
    lamfilter = 1000
    )

# Neo-Narval 2024 NEXTRA LSDpy 2.0km/s mask_folsom
VEGA_2024_NEONARVAL_NEXTRA = Experiment(
    DATAFILE = os.path.join(DATADIR, 'Vega_2024.dat'),
    #DATAFILE = os.path.join(DATADIR, 'Vega_2024_sel.dat'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5', 's12','s34', 'NEO-NARVAL 2024 (1,2,3,4,5)'],
        [[0], [1], [2], [3], [4], [0,1], [2,3], [0,1,2,3,4]]),
    lamfilter = 1000
    )

# Neo-Narval 2018 SOPHIE & NEXTRA LSDpy 1.7 km/s mask_folsom
VEGA_2018_SOPHIE_NARVAL_DUOFILE = Experiment(
    DATAFILE = os.path.join(DATADIR, 'duofilematrix.dat'),
    vrange = [-60.0, 40.0],
    vrad = -13.4,
    prot = rotper, #days
    noiselevel = 1.3,
    normalize = True,
    nights = zip(
        ['s1','s2','s3','s4','s5','s6', 's7', 's12','s34','s56', '2018 NARVAL and SOPHIE s1234567'],
        [[0], [1], [2], [3], [4], [5], [6], [0,1], [2,3], [4,5], [0,1,2,3,4,5,6]]),
    lamfilter = 1000
    )


#experiment = VEGA_2012_SOPHIE_LSDPY.kwargs

#experiment = VEGA_2018_SOPHIE_LSDJFD.kwargs
experiment = VEGA_2018_SOPHIE_LSDPY.kwargs

#experiment = VEGA_2018_NARVAL_LE.kwargs
#experiment = VEGA_2018_NARVAL_NEXTRA.kwargs
#experiment = VEGA_2018_NARVAL_NEXTRA_KM17.kwargs

#experiment = VEGA_2018_SOPHIE_NARVAL_DUOFILE.kwargs

#experiment = VEGA_2023_NEONARVAL_NEXTRA.kwargs
#experiment = VEGA_2023_NEONARVAL_NEXTRA_Q09.kwargs
#experiment = VEGA_2023_NEONARVAL_NEXTRA_SELECT.kwargs

#experiment = VEGA_2024_NEONARVAL_NEXTRA.kwargs



class Pictures(object):

    def __init__(self, *args):

        # parameters

        self.normalize = False  # normalise spectra to interval (0,1) ?
        self.r_depth = 0.5  # depth for
#        self.upper, self.lower = (0.35, 0.5), (0.15, 0.25)  # limits for vspan
        self.upper, self.lower = (0.35, 0.5), (0.1, 0.25)  # limits for vspan
#        self.extension = (0.15, 0.3) # limits of bisector area for median calculation vrad_bis
        """
        self.upper_errorbar = 0.036
        self.lower_errorbar = - 0.029
        self.rotperiod =  0.68 # 0.66149 #self.rotperiod# #self.rotperiod #rotation period
        """
        
        self.upper_errorbar = 0.0001
        self.lower_errorbar = - 0.0001
        self.rotperiod =  0.678 # 0.66149 #self.rotperiod# #self.rotperiod #rotation period
        ##self.rotperiod =  0.678
        self.extension = (0.15, 0.3) # limits of bisector area for median calculation vrad_bis
        self.d0, self.d1 = (0.1, 0.9)  # limits for bisector
        self.cpd = [24.0 / 12.5 ]  # cycles per day, where to plot a vertical lline (FRot
        
        self.cpdl = 1. / (self.rotperiod + self.upper_errorbar)  # cycles per day, low, Alina
        self.cpdh = 1. / (self.rotperiod + self.lower_errorbar)  # cycles per day, high, Alina
        self.cpd2l = 1. / self.rotperiod + 1. / (self.rotperiod + self.upper_errorbar)
        self.cpd2h = 1. / self.rotperiod + 1. / (self.rotperiod + self.lower_errorbar)
        self.cpd3l = 2. / self.rotperiod + 1. / (self.rotperiod + self.upper_errorbar)
        self.cpd3h = 2. / self.rotperiod + 1. / (self.rotperiod + self.lower_errorbar)
        self.cpd4l = 3. / self.rotperiod + 1. / (self.rotperiod + self.upper_errorbar)
        self.cpd4h = 3. / self.rotperiod + 1. / (self.rotperiod + self.lower_errorbar)
        self.cpd5l = 4. / self.rotperiod + 1. / (self.rotperiod + self.upper_errorbar)
        self.cpd5h = 4. / self.rotperiod + 1. / (self.rotperiod + self.lower_errorbar)
        self.cpd6l = 5. / self.rotperiod + 1. / (self.rotperiod + self.upper_errorbar)
        self.cpd6h = 5. / self.rotperiod + 1. / (self.rotperiod + self.lower_errorbar)
        self.cpd7l = 6. / self.rotperiod + 1. / (self.rotperiod + self.upper_errorbar)
        self.cpd7h = 6. / self.rotperiod + 1. / (self.rotperiod + self.lower_errorbar)
        self.cpd8l = 7. / self.rotperiod + 1. / (self.rotperiod + self.upper_errorbar)
        self.cpd8h = 7. / self.rotperiod + 1. / (self.rotperiod + self.lower_errorbar)
        self.cpd9l = 8. / self.rotperiod + 1. / (self.rotperiod + self.upper_errorbar)
        self.cpd9h = 8. / self.rotperiod + 1. / (self.rotperiod + self.lower_errorbar)
        self.cpd10l = 9. / self.rotperiod + 1. / (self.rotperiod + self.upper_errorbar)
        self.cpd10h = 9. / self.rotperiod + 1. / (self.rotperiod + self.lower_errorbar)

        self.cpdnew = 1.78

        #self.cpdb = 1.6062959  # cycles per day, Butkovskaya
        self.nfreq = 1024  # number of frquencies for analysis
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
        if self._vrad_mean is None:
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
            filename="lsvspan_timeseries"
            np.savetxt(filename,np.column_stack([self.time,self.vspan-mn]))
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
        #plt.xlim(-13.1, -13.01)
        #plt.ylim(0.2, 0.6)
        plt.savefig(name + self.format)
        

        name="vspan_time"
        plt.figure()
        plt.title(name)
        plt.title('vspan')
        plt.plot(self.time, self.vspan, 'o')
        plt.xlabel('time (BJD)')
        plt.ylabel('vspan (km/s)')
        plt.savefig(name + self.format)
        np.savetxt(name, np.column_stack((self.time,self.vspan)), fmt='%12.8f')
        
        name="vrad_mean_time"
        plt.figure()
        plt.title(name)
        plt.title('vrad_mean')
        plt.plot(self.time, self.vrad_mean, 'o')
        plt.xlabel('time (BJD)')
        plt.ylabel('vrad_mean (km/s)')
        plt.savefig(name + self.format)
        np.savetxt(name, np.column_stack((self.time,self.vrad_mean)), fmt='%12.8f')



    def vrad_corr_vspan(self):
        name = "vrad_corr_vspan"


#        print(self.vrad_corr)
        x = self.vrad_corr
        y = self.vspan
        print(self.time.shape, y.shape, x.shape)
        xy = np.vstack([x,y])
        z = np.mod(self.time / self.rotperiod, 1) #gaussian_kde(xy)(xy)
        I = ((z>0.15) & (z < 0.35)) | ((z >0.65) & ( z < 0.85))
        z = np.where(I, 0, 1)

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
        #plt.xlim(-0.17, 0.17)
        #plt.ylim(0.15, 0.55)

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
        ll = [self.cpdl, self.cpd2l, self.cpd3l, self.cpd4l, self.cpd5l]
        hh = [self.cpdh, self.cpd2h, self.cpd3h, self.cpd4h, self.cpd5h]

        if bars:
            for l, h in zip (ll, hh):
                #plt.gca().add_patch(plt.Rectangle((l,0), h-l, 0.3, fc='0.85', color='0.85'))
                #plt.gca().add_patch(plt.Rectangle((l,0), h-l, 0.3, fc='0.3', color='0.3'))
                #plt.gca().add_patch(plt.Rectangle((l-1,0), h-l, 0.15, fc='0.5', color='0.5'))
                #plt.gca().add_patch(plt.Rectangle((l+1,0), h-l, 0.15, fc='0.5', color='0.5'))
                #plt.gca().add_patch(plt.Rectangle((l-2,0), h-l, 0.02, fc='0.7', color='0.7'))
                #plt.gca().add_patch(plt.Rectangle((l+2,0), h-l, 0.02, fc='0.7', color='0.7'))
                plt.gca().add_patch(plt.Rectangle((l,0), h-l, 0.3, fc='orangered', color='orangered'))
                plt.gca().add_patch(plt.Rectangle((l,0), h-l, 0.3, fc='orangered', color='orangered'))
                plt.gca().add_patch(plt.Rectangle((l-1,0), h-l, 0.15, fc='darkorange', color='darkorange'))
                plt.gca().add_patch(plt.Rectangle((l+1,0), h-l, 0.15, fc='darkorange', color='darkorange'))
                plt.gca().add_patch(plt.Rectangle((l-2,0), h-l, 0.02, fc='wheat', color='wheat'))
                plt.gca().add_patch(plt.Rectangle((l+2,0), h-l, 0.02, fc='wheat', color='wheat'))
            
            """
            for n in range(1,5):
                plt.vlines([n*self.cpdnew],0,0.3)
                plt.vlines([n*self.cpdnew-1],0,0.15)
                plt.vlines([n*self.cpdnew+1],0,0.15)
                plt.vlines([n*self.cpdnew -2 ],0,0.02)
                plt.vlines([n*self.cpdnew +2],0,0.02)
            """

        maxa = max(amp)
        plt.axis([0,15.,0,maxa])
        plt.plot(self.cycles_per_day, amp, '-k',linewidth=2, color = 'dimgrey')
        plt.minorticks_on()

    def _plt_ls_vr(self, data):
        amp = np.abs(data)
#        amp /= np.max(amp)
        ll = [self.cpdl, self.cpd2l, self.cpd3l, self.cpd4l, self.cpd5l]
        hh = [self.cpdh, self.cpd2h, self.cpd3h, self.cpd4h, self.cpd5h]
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
            v0=-13.9

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

    def moving_peaks_simple_time(self):
        nbins = 512
        sa = self.analyzer
        time = sa.time

        meanprofile = np.median(sa.intensity, axis=0)

        F = np.zeros((sa.nvelocity, 2))
        F[:,0] = 1
        F[:,1] = meanprofile

        diff = np.zeros((sa.nobs, sa.nvelocity))

        lams = np.zeros(sa.nobs)
        offsets = np.zeros(sa.nobs)

        for i in range(sa.nobs):
            pp = np.linalg.lstsq(F, sa.intensity[i], rcond=None)[0]
            diff[i,:] = sa.intensity[i] - (pp[0] + pp[1] * meanprofile)
            lams[i] = pp[1]
            offsets[i] = pp[0]


        plt.figure()
        plt.title('lambda')
        plt.plot( time, lams, '.')


        print(pp)


        time = sa.time[ lams > 0.98 ]

        # diff = sa.intensity - np.median(sa.intensity, axis=0)

        phasebins = np.linspace(0, self.rotperiod, nbins+1)
        #I = np.digitize(np.mod(sa.time, self.rotperiod), phasebins)


        tt = np.mod(sa.time, self.rotperiod)

        for na, nightlist in zip(['s1','s2','s3','s4','s5'], [[0], [1], [2], [3], [4], [5] ]):

            plt.figure(figsize=(6,10))
            plt.title('night'+na)

            res = np.zeros((nbins, sa.nvelocity))

            for night in nightlist:
                I = self.analyzer.list_index[night]

                ttt = time - time[I][0]
                for i, p1p2 in enumerate(zip(phasebins[0:-1], phasebins[1:])):
                    p1, p2 = p1p2
                    II = (p1 <= ttt[I]) & (ttt[I] < p2)
                    if sum(II) >=1:
                        res[i,:] = np.mean(diff[I][II], axis=0)


            minmax = np.max(np.abs(res))
            print(minmax)

            res /= minmax
            res = np.sign(res) * np.abs(res)**0.75

            v0=-13.01

            plt.imshow(
                res,
                cmap=plt.cm.gray_r,
                aspect='auto',
                interpolation='None',
                origin='lower',
                extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],
                vmin=-0.8, vmax=0.8)

    def moving_peaks_simple_per_night_play(self):
        nbins = 128
        sa = self.analyzer
        time = sa.time

        meanprofile = np.median(sa.intensity, axis=0)
        meanp = (1-meanprofile) / np.sqrt(np.sum((1-meanprofile)**2))

        F = np.zeros((sa.nvelocity, 3))     #linear model
        F[:,0] = 1          #constant
        F[:,1] = meanp      # mean profile
        d = meanp[1:]-meanp[:-1]    #first derivative
        F[:-1,2] = 0.5 * d
        F[1:, 2] += 0.5 * d


        diff = np.zeros((sa.nobs, sa.nvelocity))

        lams = np.zeros(sa.nobs)
        offsets = np.zeros(sa.nobs)
        translat = np.zeros(sa.nobs)

        for i in range(sa.nobs):
            pp = np.linalg.lstsq(F, sa.intensity[i], rcond=None)[0]
            alpha = np.sum((1-sa.intensity[i]) * meanp)
            #diff[i,:] = sa.intensity[i] - (1-alpha*meanp)
            #pp[2]=0     #translation has no significant impact
            #pp[1]=0    scaling has major impact!
            #pp[0]=0    has impact
            diff[i,:] = sa.intensity[i] - np.dot(F, pp) #difference between obs and 3 component adjutstment
            #diff[i,:] = sa.intensity[i] - meanprofile
            lams[i] = pp[1]             #coefficient of mean profile
            offsets[i] = pp[0]          #shift in intensity (continuum normalization error)
            translat[i] = pp[2]         #shift in velocity

        #plt.figure()
        #plt.plot(meanp)

        #plt.figure()
        #plt.plot(meanp[1:] - meanp[:-1])

        #plt.figure()
        #plt.title('lambda')
        #plt.plot( time, lams, '.')

        #plt.figure()
        #plt.title('offset')
        #plt.plot( time, offsets, '.')

        #plt.figure()
        #plt.title('translat')
        #plt.plot( time, translat, '.')

        #lamfilter = lams < experiment['lamfilter']      #throw out too extreme profile changes

        #diff = sa.intensity - np.median(sa.intensity, axis=0)




        #phasebins = np.linspace(0, self.rotperiod, nbins+1)
        #I = np.digitize(np.mod(sa.time, self.rotperiod), phasebins)


        #tt = np.mod(time, self.rotperiod)       #time modulo rotationperiod
        
        for na, nightlist in experiment['nights']:          #for all nights

            plt.figure(figsize=(6,10))
            plt.title('nights '+na)

            res = np.zeros((nbins, sa.nvelocity))

            ## collecting data for list of nights
            TT = []
            val = []

            for night in nightlist:
    
                I = self.analyzer.list_index[night]         #indices showing belonging to a night

                dd = sa.intensity[I]
                H = np.median(dd, axis=0)#-1
                #H = np.convolve([0.25,0.5,0.25], H)
                #H = np.convolve([0.25,0.5,0.25], H)
                #H = H[2:-2]+1
                #plt.figure()
                # plt.plot(H)
                dd -= H

                TT += list(time[I])
                val += list(dd)

            val = np.row_stack(val)

            res, bins, mask = self.analyzer.spectrum_matrix_full(
                    TT, val, period=self.rotperiod, nphase=nbins, method=np.median)
            """
                for i, p1p2 in enumerate(zip(phasebins[0:-1], phasebins[1:])):
                    p1, p2 = p1p2
                    II = (p1 <= tt[I]) & (tt[I] < p2) #& lamfilter[I]
                    if sum(II) >=1:
                        res[i,:] = np.median(dd[II], axis=0)
            """

            #minmax = np.max(np.abs(res))
            #print(minmax)

            #res /= minmax
            #res = np.sign(res) * np.abs(res)**0.75


            np.savetxt('moving_simple_per_night.txt', res)
            np.savetxt('velocitybins.txt', self.velocity)

            v0=-13.9
            plt.imshow(-res, cmap=plt.cm.gray_r,
                       aspect='auto',interpolation='bicubic',
                       origin='lower',
                       extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],
                       vmin=-0.0002,vmax=0.0002
                       )
            plt.savefig('nights'+na + self.format)

            """
            plt.imshow(
                res,
                #cmap=plt.cm.bwr,  #plt.cm.gray_r,
                cmap=plt.cm.gray_r,
                aspect='auto',
                interpolation='bicubic',
                origin='lower',
                extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],
                vmin=-0.8, vmax=0.8)
            """
            plt.xticks([])
            rv = np.array([-20, -10, 0, 10, 20 ])
            st = [ str(i) for i in rv]
            plt.xticks(rv, st)
            plt.yticks([])
            yt = [0,0.2,0.4,0.6,0.8, 1]
            plt.yticks(yt,[str(t) for t in yt])
            prof_beg_0 = -rotation_vsini
            prof_end_0 =  rotation_vsini
            full_prof_beg_0 = prof_beg_0-full_prof_add
            full_prof_end_0 = prof_end_0+full_prof_add

            plt.vlines([prof_beg_0,0,prof_end_0],0,1)
            plt.vlines([full_prof_beg_0,0,full_prof_end_0],0,1,linestyles="dashed")

            #plt.vlines([-22, 0, 22], 0,1)
            #plt.vlines([-33, 0, 33], 0,1,linestyles='dashed')
            plt.hlines([1./2],  full_prof_beg_0,full_prof_end_0)
            plt.xlim(self.velocity[0]-v0, self.velocity[-1]-v0)
            plt.xlabel('velocity [km/s]')
            plt.ylabel('phase fraction of period]')

    def moving_peaks_simple_per_night(self):
        nbins = 128
        sa = self.analyzer
        time = sa.time

        meanprofile = np.median(sa.intensity, axis=0)
        meanp = (1-meanprofile) / np.sqrt(np.sum((1-meanprofile)**2))

        F = np.zeros((sa.nvelocity, 3))     #linear model
        F[:,0] = 1          #constant
        F[:,1] = meanp      # mean profile
        d = meanp[1:]-meanp[:-1]    #first derivative
        F[:-1,2] = 0.5 * d
        F[1:, 2] += 0.5 * d


        diff = np.zeros((sa.nobs, sa.nvelocity))

        lams = np.zeros(sa.nobs)
        offsets = np.zeros(sa.nobs)
        translat = np.zeros(sa.nobs)

        for i in range(sa.nobs):
            pp = np.linalg.lstsq(F, sa.intensity[i], rcond=None)[0]
            alpha = np.sum((1-sa.intensity[i]) * meanp)
            #diff[i,:] = sa.intensity[i] - (1-alpha*meanp)
            #pp[2]=0     #translation has no significant impact
            #pp[1]=0    scaling has major impact!
            #pp[0]=0    has impact
            diff[i,:] = sa.intensity[i] - np.dot(F, pp) #difference between obs and 3 component adjutstment
            #diff[i,:] = sa.intensity[i] - meanprofile
            lams[i] = pp[1]             #coefficient of mean profile
            offsets[i] = pp[0]          #shift in intensity (continuum normalization error)
            translat[i] = pp[2]         #shift in velocity

        plt.figure()
        plt.plot(meanp)

        plt.figure()
        plt.plot(meanp[1:] - meanp[:-1])

        plt.figure()
        plt.title('lambda')
        plt.plot( time, lams, '.')

        plt.figure()
        plt.title('offset')
        plt.plot( time, offsets, '.')

        plt.figure()
        plt.title('translat')
        plt.plot( time, translat, '.')

        lamfilter = lams < experiment['lamfilter']      #throw out too extreme profile changes

        #diff = sa.intensity - np.median(sa.intensity, axis=0)




        phasebins = np.linspace(0, self.rotperiod, nbins+1)
        #I = np.digitize(np.mod(sa.time, self.rotperiod), phasebins)


        tt = np.mod(time, self.rotperiod)       #time modulo rotationperiod
        
        for na, nightlist in experiment['nights']:          #for all nights

            plt.figure(figsize=(6,10))
            plt.title('nights '+na)

            res = np.zeros((nbins, sa.nvelocity))

            ## collecting data for list of nights
            TT = []
            val = []

            for night in nightlist:
    
                I = self.analyzer.list_index[night]         #indices showing belonging to a night

                dd = diff[I]
                dd -= np.median(dd, axis=0)

                TT += list(time[I])
                val += list(dd)

            val = np.row_stack(val)

            res, bins, mask = self.analyzer.spectrum_matrix_full(
                    TT, val, period=self.rotperiod, nphase=nbins, method=np.median)
            """
                for i, p1p2 in enumerate(zip(phasebins[0:-1], phasebins[1:])):
                    p1, p2 = p1p2
                    II = (p1 <= tt[I]) & (tt[I] < p2) #& lamfilter[I]
                    if sum(II) >=1:
                        res[i,:] = np.median(dd[II], axis=0)
            """

            #minmax = np.max(np.abs(res))
            #print(minmax)

            #res /= minmax
            #res = np.sign(res) * np.abs(res)**0.75


            np.savetxt('moving_simple_per_night.txt', res)
            np.savetxt('velocitybins.txt', self.velocity)

            v0=-13.9
            plt.imshow(-res, cmap=plt.cm.gray_r,
                       aspect='auto',interpolation='bicubic',
                       origin='lower',extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],
                       vmin=-0.0002,vmax=0.0002)
            plt.savefig('nights'+na + self.format)

            """
            plt.imshow(
                res,
                #cmap=plt.cm.bwr,  #plt.cm.gray_r,
                cmap=plt.cm.gray_r,
                aspect='auto',
                interpolation='bicubic',
                origin='lower',
                extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],
                vmin=-0.8, vmax=0.8)
            """
            plt.xticks([])
            rv = np.array([-20, -10, 0, 10, 20 ])
            st = [ str(i) for i in rv]
            plt.xticks(rv, st)
            plt.yticks([])
            yt = [0,0.2,0.4,0.6,0.8, 1]
            plt.yticks(yt,[str(t) for t in yt])
            prof_beg_0 = -rotation_vsini
            prof_end_0 =  rotation_vsini
            full_prof_beg_0 = prof_beg_0-full_prof_add
            full_prof_end_0 = prof_end_0+full_prof_add

            plt.vlines([prof_beg_0,0,prof_end_0],0,1)
            plt.vlines([full_prof_beg_0,0,full_prof_end_0],0,1,linestyles="dashed")

            #plt.vlines([-22, 0, 22], 0,1)
            #plt.vlines([-33, 0, 33], 0,1,linestyles='dashed')
            plt.hlines([1./2],  full_prof_beg_0,full_prof_end_0)
            plt.xlim(self.velocity[0]-v0, self.velocity[-1]-v0)
            plt.xlabel('velocity [km/s]')
            plt.ylabel('phase fraction of period]')
            plt.savefig('nights'+na + self.format)


    def moving_peaks_simple(self):
        nbins = 128
        sa = self.analyzer
        time = sa.time

        meanprofile = np.median(sa.intensity, axis=0)
        meanp = (1-meanprofile) / np.sqrt(np.sum((1-meanprofile)**2))

        F = np.zeros((sa.nvelocity, 2))
        F[:,0] = 1
        F[:,1] = meanp

        diff = np.zeros((sa.nobs, sa.nvelocity))

        lams = np.zeros(sa.nobs)
        offsets = np.zeros(sa.nobs)

        for i in range(sa.nobs):
            pp = np.linalg.lstsq(F, sa.intensity[i], rcond=None)[0]
            alpha = np.sum((1-sa.intensity[i]) * meanp)
            #diff[i,:] = sa.intensity[i] - (1-alpha*meanp)
            diff[i,:] = sa.intensity[i] - (pp[0] + pp[1] * meanp)
            #diff[i,:] = sa.intensity[i] - meanprofile
            lams[i] = pp[1]
            offsets[i] = pp[0]

        plt.figure()
        plt.plot(meanp)


        plt.figure()
        plt.title('lambda')
        plt.plot( time, lams, '.')

        plt.figure()
        plt.title('offset')
        plt.plot( time, offsets, '.')

        #lamfilter = lams < -0.465
        # lamfilter = lams > 0.425
        lamfilter = lams > -1000

        # diff = sa.intensity - np.median(sa.intensity, axis=0)

        phasebins = np.linspace(0, self.rotperiod, nbins+1)
        #I = np.digitize(np.mod(sa.time, self.rotperiod), phasebins)


        tt = np.mod(time, self.rotperiod)

        for na, nightlist in experiment['nights']:

            plt.figure(figsize=(6,10))
            plt.title('night'+na)

            res = np.zeros((nbins, sa.nvelocity))

            for night in nightlist:
                I = self.analyzer.list_index[night]

                for i, p1p2 in enumerate(zip(phasebins[0:-1], phasebins[1:])):
                    p1, p2 = p1p2
                    II = (p1 <= tt[I]) & (tt[I] < p2) & lamfilter[I]
                    if sum(II) >=1:
                        res[i,:] = np.mean(diff[I][II], axis=0)


            minmax = np.max(np.abs(res))
            print(minmax)

            res /= minmax
            res = np.sign(res) * np.abs(res)**0.75

            v0=-13.01

            plt.imshow(
                res,
                #cmap=plt.cm.bwr,  #plt.cm.gray_r,
                cmap=plt.cm.gray_r,
                aspect='auto',
                interpolation='None',
                origin='lower',
                extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],
                vmin=-0.8, vmax=0.8)








    def moving_peaks_signoise(self):
        name = 'moving_peaks_eq_width'
        #return #TODO REMOVE
        vul = 1.-self.inte
        eqwidth = self.analyzer.eqwidth()
        signois = self.analyzer.meansignoise()


#        for na, nightlist in zip(['s:1', 's:2', 's:3', 's:123456'], [[0], [1], [2], [0,1,2,3,4,5]]):

        minmax=[]

        for na, nightlist in experiment['nights']:

            VV =[]
            TT =[]
            plt.figure(figsize=(6,10))
            plt.title('signoise night'+na)
            for night in nightlist:
                I = self.analyzer.list_index[night]
                TT += list(self.time[I])
                pp = np.poly1d(np.polyfit(signois[I], eqwidth[I], deg=3))
                fac = pp(signois[I])[:,np.newaxis]
                temp =vul[I]/fac
                temp -=np.median(temp,axis=0)
                VV += list(temp)

            val = np.row_stack(VV)
            val -= np.median(val,axis=0)

            nphase=128
            tmp, bins, mask = self.analyzer.spectrum_matrix_full(TT, val, period=self.rotperiod, nphase=nphase, method=np.median)
            #print("minmax:",np.min(tmp),np.max(tmp))
            minmax.append([np.min(tmp),np.max(tmp)])

            v0=-13.01

            #plt.imshow(tmp, cmap=plt.cm.gray_r, aspect='auto',interpolation='none', origin='lower',extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1])
            plt.imshow(tmp, cmap=plt.cm.gray_r, aspect='auto',interpolation='bicubic', origin='lower',extent=[self.velocity[0]-v0, self.velocity[-1]-v0,0,1],vmin=-0.0002,vmax=0.0002)
            #vmin=-0.0003,vmax=0.0003
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

        #print("final minmax:",np.min(np.array(minmax)),np.max(np.array(minmax)))

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

    
    def spot_density(self):
        nbins = 64
        theta, phi, res = self.analyzer.spot_density(nbins=nbins)
        map = Basemap(projection='moll', lon_0=0)
        plt.figure()
        #map.drawcoastlines()
        map.imshow(res[:,::-1].T, interpolation='bicubic')

        map.drawmeridians(np.linspace(-180, 180, 12))
        map.drawparallels(np.linspace(-90, 90, 6))
        #map.plot(180*phi/np.pi, 128*[180*self.kwargs['angle']/np.pi], '-b', latlon=True)
        map.plot(180*phi/np.pi, 128*[0], '-r', latlon=True)
        #map.plot(180*phi/np.pi, 128*[-180*self.kwargs['angle']/np.pi], '-b', latlon=True)
        plt.show()

 
        plt.savefig('spot_density.pdf')

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
    myPics.vrad_mean_vspan()
    myPics.vrad_corr_vspan()
#    myPics.vrad_mean_vspan()
#    myPics.vrad_corr_vspan()
#    myPics.vrad_mean_skew()
#    myPics.vrad_mean_std()
#    myPics.vrad_corr_skew()
#    myPics.skew_vspan()
#    myPics.ts_skew()
#    myPics.ts_vspan()
##   myPics.ls_spec_vrad_skew()
##    myPics.ls_spec_vrad_mean()
#    myPics.ls_spec_all3()
##    myPics.ls_spec_vrad_corr()
##    myPics.ls_spec_vrad_bis()
##    myPics.ls_spec_vspan()
#    myPics.ls_spec_eqwidth()
#    myPics.bisector_time()
#    myPics.bisector_width()
#    myPics.ls_window(-rw-rw-r-- 1 hols hols   22112 Jun  4 14:59  spectralutil_old.py

#    alldata = [self.time, self.inte, self.vrad_mean, self.vrad_corr, self.vspan, self.vrad_skew, self.vrad_std]
#   myPics.bayes_freq_vrad_mean()
<<<<<<< HEAD
    myPics.moving_peaks_simple_per_night()
    # myPics.moving_peaks_simple_per_night_play()
=======
>>>>>>> 524f648700e93cd37ef94be20a5241850df2dd71

    myPics.moving_peaks_simple_per_night()
    #myPics.moving_peaks_simple_per_night_play()


##    myPics.spot_density()

    #myPics.moving_peaks_signoise()

    # myPics.moving_peaks_simple()
    # myPics.moving_peaks_simple_time()
    # myPics.moving_peaks_time()
    ###    myPics.estrotentropy()

    #myPics.saveData("time_vrad_mean.dat", [6142.+myPics.time, myPics.vrad_mean])  # first column
    #myPics.saveData("time_vrad_corr.dat", [myPics.time, myPics.vrad_corr])  # first column
    #myPics.saveData("time_vrad_bis.dat", [6142.+myPics.time, myPics.vrad_bis])  # first column
    #myPics.saveData("time_vspan.dat", [6142.+myPics.time, myPics.vspan])  # first column
    #myPics.saveData("time_skew.dat", [myPics.time, myPics.vrad_skew])  # first column




    # myPics.ew_noise_corr()

    plt.show()
