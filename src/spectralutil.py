'''
Created on May 22, 2014

@author: hols
'''
from scipy import interpolate
from scipy.signal import spectral
import matplotlib.pyplot as plt
import numpy as np
import json
import os

from scipy.interpolate.fitpack2 import UnivariateSpline
from scipy.interpolate import interp1d

VEGAPERIOD = 0.678 # 0.98739 * 0.678  # period in days

#for full profile 60,124
#for short profile 72,112
def load_data(DATAFILE, nval, rangei, vrange, noiselevel):
    """
    loads a filematrix type data file

    nval       : number of bins for the spectrum
    rangei     : (a,b) taking indicies from a to b (included)
    vrange     : alternatively prescribing velocity range
    noiselevel : for outlier removal
    """

    coltime = 0  # colum of time values
    colspec = 1
    colval = colspec + nval
    colvul = colval + nval

    data = np.loadtxt(DATAFILE)
    time = data[:, coltime:colspec].ravel()
    velocity = data[0, colspec:colval]  # velocities of bins
    intens = data[:, colval:colvul]  # intensities
    signoise = data[:, colvul:colvul+nval]
    meani    = intens.mean(axis=0)  # mean intensity
    diff     = intens - meani  # fluctuation around mean

    # selecting velocity bins
    if vrange is not None:          # wenn vrange verschieden ist von none
        I1, I2 = np.where(velocity >= vrange[0])[0], np.where(velocity<=vrange[1])[0]
        rangei = (I1[0], I2[-1])    # -1 ist letztes element des vektors
    rangel   = np.arange(rangei[0], rangei[1])  # range of spectral line

    stdi = diff.ravel().std()  # variance of f TODO better wirh std ?       #macht vektor flach, dann varianz
    #ueber langen vektor (ntimes*nval), dh varianz vom gesamten Bild
    t0 = time[0]  # first rime

    # outlier removal
    I = np.where(diff.std(axis=1) < noiselevel * stdi)[0]
    plt.plot(np.std(diff[I, :], axis=1))

    plt.show()

    print("reducing from ", time.shape , "to ", len(I), " spectral lines")

    time = time[I]
    velocity = velocity[rangel]

    intens = intens[I, :]
    intens = intens[:, rangel]

    signoise = signoise[I,:]
    signoise = signoise[:,rangel]

    # recomputing mean and diff
    meani = intens.mean(axis=0)
    diff = intens - meani

    #grouping into nights
    tt = time - int(time[0])

    listOfNightsIndices = []
    listOfNightsTimes = []
    listOfNightsIntens = []

    for i in range(100):
        I = np.where( (tt-i) * (i+1-tt) > 0)[0]     # quadratische Ungleichung waehlt Tag aus
                                # zwischen Tag i und Tag i+1
                                # TODO: ok fuer Europa, fuer Chile kann Nacht
                                # ueber integersprung gehen...
                                # Klasse Observatorium einrichten!!
        if len(I) > 0:
            listOfNightsIndices.append(I)
            listOfNightsTimes.append( time[I] )
            listOfNightsIntens.append ( intens[I] )

    return time, velocity, intens, signoise, listOfNightsTimes, listOfNightsIntens, listOfNightsIndices

def spectrum_matrix(time, quantity, **kwargs):
    """
    binns the spec into nphase bins
    """
    period = kwargs.get('period', VEGAPERIOD)
    nphase = kwargs.get('nphase', 128)

    bins = np.linspace(0, period, nphase + 1)
    res = np.zeros((nphase, quantity.shape[1]))
    nn = np.zeros(nphase)
    ii = np.digitize(np.mod(time, period), bins)
    for i, s in zip(ii, quantity):
        res[i-1, :] += s
        nn[i-1] += 1

    nn = np.where(nn>0,nn,1)
    res[:,:] = res / nn[:,np.newaxis]
    return res[:, 1:-1], bins[:-1]


def bisector(velocity, intensity, depth):
    """
    computes the bisector of the
    absorption (!!) profile velocity, intensity
    at the point depth
    """

    n = intensity.size # number of bins

    i = np.min(np.where(intensity < depth))
    j = np.max(np.where(intensity < depth))

    if i<=0 or j>=n-1:
        raise Exception('wrong depth')

    x = velocity[i - 1] +  (depth - intensity[i - 1]) * (velocity[i] - velocity[i - 1]) / (intensity[i] - intensity[i - 1])
    y = velocity[j + 1] + (depth - intensity[j + 1]) * (velocity[j] - velocity[j + 1]) / (intensity[j] - intensity[j + 1])
    return (x + y) / 2.
   
def bisector_v(velocity, intensity, depth):
    """
    the bisector for the values in depth
    """
    res = [bisector(velocity, intensity, d) for d in depth]
    return res

def remove_night_trend(time, value):
    """
    removing individual trends for each night
    """
    time -= int(time[0])
    nnight = int(time[-1] - time[0]) + 1
    t0 = int(time[-1])
    for n in range(nnight):
        I = np.where((time - n) * (n + 1 - time) >= 0.)
        value[I], p, pp = remove_trend(time[I], value[I])

    return value

def remove_trend(time, value):
    """
    removing individual trends for each night
    """
    p = np.polyfit(time, value, 2)

    print("shap of p ", p.shape)
    n, m = value.shape
    pp = np.zeros(value.shape)
    for i in range(m):
        pp[:, i] = np.polyval(p[:, i], time)
    value -= pp

    return value, p, pp


def vspan(upper, lower, time, velocity, intensity):
    """
    compute vspan as function of time
    """
    nn = 20
    tmp = np.zeros(len(time))


    du = np.linspace(upper[0], upper[1], nn)
    dl = np.linspace(lower[0], lower[1], nn)

    for i in range(len(time)):
        minval = intensity[i].min()

        bu = bisector_v(velocity, intensity[i], minval + du * (1 - minval))
        bl = bisector_v(velocity, intensity[i], minval + dl * (1 - minval))
        tmp[i] = np.median(bu) - np.median(bl)


    return tmp


def vrad_bis(extension, time, velocity, intensity):
    """
    calculates vrad as the median of part of bisector
    """
    nn = 20
    tmp = np.zeros(len(time))


    dl = np.linspace(extension[0], extension[1], nn)

    for i in range(len(time)):
        minval = intensity[i].min()

        bl = bisector_v(velocity, intensity[i], minval + dl * (1 - minval))
        tmp[i] = np.median(bl)

    return tmp


def vrad(relative_depth, velocity, intensity):

    nt, nv = intensity.shape

    res = np.zeros(nt)
    for i in xrange(nt):
        minval = intensity[i, :].min()
        depth = minval + (1 - minval) * relative_depth
        I = np.where(intensity[i, :] <= depth)[0]
        iii = intensity[i, I]
        iii = iii.max() - iii
        res[i] = np.sum(iii * velocity[I]) / np.sum(iii)
    return res

def vrad_corr(velocity, intensity, usemean=True, depth=0.5):
    """
    correlation based vrad
    """

    base = intensity[0]
    if usemean:
        base = intensity.mean(axis=0)

    minval = base.min()
    I = np.where(base >= minval + depth * (1 - minval))[0]
    # base = base[I]
    # intensity=intensity[I]
    n, d = intensity.shape
    # d = len(I)
    mask = np.correlate(np.ones(d), np.ones(d), mode='same')

    res = np.zeros(n)
    for r, inte in enumerate (intensity):
        tmp = np.correlate(1 - inte, 1 - base, mode='same')
        # tmp /= mask

        # plt.plot(velocity, tmp)

        # plt.show()

        i = np.argmax(tmp)
        poly = np.polyfit(velocity[i - 3:i + 4], tmp[i - 3:i + 4], deg=2)

        res[r] = -0.5 * poly[1] / poly[0]


    return res

    # (base, intensity[i], mode=)

def estimate_translat_scaling(intensity):
    """
    removing translational movements and scaling by a simple least square fit
    """
    n, d = intensity.shape
    median = np.mean(intensity, axis=0)
    v = np.linspace(0, 1, median.size)
    tck = interpolate.splrep(v, median, s=2)
    dermed = interpolate.splev(v, tck, der=1)

    F = np.column_stack((median, -dermed))

    print(F.shape)

    res = np.linalg.lstsq(F, intensity.T)
    scaling, translat = res[0]

    movement = (np.dot(F, res[0])).T
    delta = intensity - movement
    at = movement - median
    return scaling, translat, median, at, delta

def vrad_translat(intensity):
    scaling, translat, median, at, delta = estimate_translat_scaling(intensity)
    return translat


def FILEMATRIX_to_JSON(
        DATAFILE,     # datafile containing the filematrix
        nval=201,     # number of velocity bins
    ):
    data = np.loadtxt(DATAFILE)
    time = data[:, 0].ravel()
    print(data.shape, time.shape, time[0].shape)
    velocity = data[0, 1:nval+1]  # velocities of bins
    intensity = data[:, nval+1:2*nval+1]  # intensities
    errors = data[:, 2*nval+1:]
    name, ext = os.path.splitext(DATAFILE)

    I = np.argsort(time)
    res = {
        'name': name + '.json',
        'description': "a little info",
        'noiselevel': 0.7,
        'nvals': nval,
        'range': [72, 201],
        'time': time[I].tolist(),
        'velocity': velocity.tolist(),
        'intensity': intensity[I].tolist(),
        'errors': errors[I].tolist()
    }
    with open(os.path.join(os.path.dirname(DATAFILE), name+".json"), 'w') as outfile:
        json.dump(res, outfile, indent=2)

    print (time[I].shape, intensity[I].shape, errors[I].shape)

class SpectralAnalyser:
    """
    a class for the modeling of spectral lines a la Boehm
    """

    def __init__(self, jsonfile):

        res = self._load_json(jsonfile)

        # reading the file matrix
        self._name = res['name']
        self._time = np.array(res['time'])
        self._velocity = np.array(res['velocity'])
        self._intensity = np.array(res['intensity'])  # bar intensity before selection
        self._errors = np.array(res['errors'])

        # initializing indices for data reduction
        self.usedindex = res.get('usedindex', np.full(self._intensity.shape[0], True))
        self.vrange = res.get('vrange', np.full(self._intensity.shape[1], True))

    @property
    def name(self):
        return self._name

    @property
    def time(self):
        """
        return: time stamps of used data
        """
        return self._time[self.usedindex]

    @property
    def velocity(self):
        """
        return: used velocities
        """
        return self._velocity[self.vrange]

    def reverse_intensity_m(self):
        """
        substracts continuum based on continuum = 1
        suffix _m stands for modifies data
        """
        self._intensity = 1 - self._intensity

    @property
    def v0(self):
        """
        location of minimum spectrum
        """
        return self.velocity[np.argmin(self.mean_intensity)]

    def normalize_flux_m(self):
        tmp = 1 - self._intensity
        tmp = tmp / np.sum(tmp, axis=1)[:, np.newaxis]
        self._intensity = 1-tmp

    def outlier_removal(self, noiselevel, **kwargs):
        # outlier_removal based on usedindex set

        intensity = self._intensity[:, self.vrange]
        mean = self.mean_intensity
        tmp = intensity - mean[np.newaxis, :]  # fluctuation around mean
        I, = np.where(tmp.std(axis=1) >= noiselevel)
        self.usedindex = np.full(self._time.shape[0], True)
        self.usedindex[I] = False

        print("reducing from ", self._time.shape , "to ",
              np.sum(self.usedindex), " spectral lines")

    def filtered_intensity(self, F):
        intensity = self.intensity
        n = (F.shape[0]-1)//2
        for i in range(intensity.shape[0]):
            intensity[i] = np.convolve(F, intensity[i], 'full')[n:-n]
        return intensity

    def get_intensity(self):
        """
        return: used intensity
        TODO: fade out and replace by intensity
        """
        tmp = self._intensity[self.usedindex]
        return tmp[:, self.vrange]

    @property
    def intensity(self):
        tmp = self._intensity[self.usedindex]
        return tmp[:, self.vrange]

    @property
    def errors(self):
        tmp = self._errors[self.usedindex]
        return tmp[:, self.vrange]

    @property
    def median_intensity(self, **kwargs):
        return np.median(self.intensity, axis=0)

    @property
    def mean_intensity(self, **kwargs):
        return np.mean(self.intensity, axis=0)

    @property
    def min_intensity(self):
        return np.min(self.intensity)

    @property
    def deltavbin(self):
        return self.velocity[1]-self.velocity[0]

    @property
    def nobs(self):
        return np.sum(self.usedindex)


    def std_over_time(self):
        """
        return: time serries of deviation (std) of actual
        spectrum from global mean
        """
        mean = self.mean_intensity
        return np.std(self._intensity[:, self.vrange] - mean[np.newaxis, :], axis=1)

    def std_normalized_over_time(self):
       tmp = 1 - self._intensity[:, self.vrange]
       tmp = tmp / tmp.sum(axis=1)[:,np.newaxis]
       mean = tmp.mean(axis=0)
       return np.std(tmp - mean[np.newaxis, :], axis=1)

    def list_of_new_night_indices(self, **kwargs):
        # grouping into nights
        # delta: float fraction of day used to determine observational gap

        delta = kwargs.get('delta', 0.5)

        dt = self._time[1:] - self._time[:-1]
        I, = np.where(dt > delta)  # gap larger than delta

        return I

    def number_of_nights(self, **kwargs):
        return len(self.list_of_new_night_indices(**kwargs))+1

    def indices_of_night(self, n, **kwargs):
        I = self.list_of_new_night_indices(**kwargs)
        if n == 0:
            return np.arange(I[0]+1)
        if n == len(I):
            return np.arange(I[n-1]+1, len(self._time))
        return np.arange(I[n-1]+1, I[n]+1)

    def mask_of_night(self, n, **kwargs):
        """
        out: ndarray of boolean which is true only in night n
        """
        tmp = np.full(self._time.shape[0], False)
        tmp[self.indices_of_night(n, **kwargs)] =  True
        return tmp

    def used_indices_of_night(self, n, **kwargs):
        """
        out: ndarray of indices of night n
        """
        I = self.mask_of_night(n, **kwargs)
        return np.arange(self._time.shape[0])[I * self.usedindex]

    def indices_of_used_nights(self, **kwargs):
        """
        out: list of ndarrays of inidices of used (i.e. downsampled) data
        """
        return [self.used_indices_of_night(n) for n in range(self.number_of_nights())]

    def remove_indices(self, I):
        """
        sets usedindices to false on I
        """
        self.usedindex[I] = False

    def time_interval(self, tmin, tmax):
        I, = np.where((self._time-tmin)*(self._time-tmax) <= 0)
        return I

    def remove_time_interval_indices(self, tmin, tmax):
        """
        indices with respect to original
        time selected time stamps
        """
        self.remove_indices(self.time_interval(tmin, tmax))

    def velocity_range(self, vrange):
        """
        setting the velocity range
        vrange: list of indices
        """
        self.vrange[:] = False
        self.vrange[vrange] = True

    def export_json(self, **kwargs):
        res = {
            'name': kwargs.get('outfile', 'vega.json'),
            'description': "reduction from {}".format(self._name),
            # 'nvals': int(sum(self.vrange)),
            # 'range': [0, int(sum(self.vrange))-1],
            'time': self._time.tolist(),
            'velocity': self._velocity.tolist(),
            'intensity': self._intensity.tolist(),
            'errors': self._errors.tolist(),
            'usedindex': self.usedindex.tolist(),
            'vrange': self.vrange.tolist(),
        }
        with open(kwargs.get('outfile', 
            os.path.join(os.path.dirname(__file__), 'vega.json')), 'w') as outfile:
            json.dump(res, outfile, indent=2)

    def _load_json(self, fname):
        with open(fname, 'r') as infile:
            res = json.load(infile)
        return res

    def spectrum_smooth(self, factor=None):
        """
        cubic spline interpolation of spectrum
        """
        intens = self.intensity
        res = [ 
            UnivariateSpline(self.velocity, intens[i], s=factor)(self.velocity)
                for i in range(intens.shape[0])
                ]
        return np.array(res)

    def rv_mean(self, relative_depth=1):
        """
        mean position of spectrum
        """

        # computing depth for truncation
        # relative_depth = 0 at the bottom
        # relative_depth = 1 truncate at 1
        # note intensity is between mv and 1 (<-- continuum)
        mv = self.mean_intensity.min()
        depth = mv + (1 - mv) * relative_depth
        I, = np.where(self.mean_intensity >= depth)

        vv = self.velocity
        res = 1 - self.intensity
        res[:, I] = 0
        v = np.sum(res * vv[np.newaxis, :], axis=1)
        v /= np.sum(res, axis=1)
        return v


    def rv_corr(self, relative_depth=1, use='mean'):
        """
        correlation based vrad
        """

        base = self.intensity[0]
        if 'mean' == use:
            base = self.mean_intensity

        mv = self.mean_intensity.min()
        depth = mv + (1 - mv) * relative_depth
        I = np.where(self.mean_intensity <= depth)[0]

        # working only with a reduced part of the spectrum
        inte = self.intensity[:, I]
        # m = np.ones(I.size)
        # mask = np.correlate(m, m, mode='full')

        # computing factor to map correlation based 
        # shifts to vrads...
        eps = 0.001  # small number to shift 
        v0 = self.velocity - 0.5 * self.deltavbin  # evaluation at center
        v1 = v0 - eps * self.deltavbin  # shifted centers
        tmp = interp1d(v0, base)
        ii0 = tmp(v0[I])
        ii1 = tmp(v1[I])
        tmp = np.correlate(1 - ii1, 1 - ii0, mode='full')

        n = 3  # offset
        i0 = np.argmax(tmp)
        poly = np.polyfit(np.arange(-n, n + 1), tmp[i0 - n:i0 + n + 1], deg=2)
        pos = -0.5 * poly[1] / poly[0]  # position of maximum of polynomial

        factor = self.deltavbin * eps / pos

        # now doing the same for the spectral data
        res = np.zeros(self.nobs)
        for r, i in enumerate (inte):
            tmp = np.correlate(1 - i, 1 - base[I], mode='full')
            i = np.argmax(tmp)
            poly = np.polyfit(np.arange(-n, n + 1), tmp[i - n:i + n + 1], deg=2)
            res[r] = -0.5 * poly[1] / poly[0] + (i - i0)

        return factor * res

    def bisector_borders(self):
        #mv = self.intensity.min(axis=0)
        #depth = mv + (1 - mv) * depth
        #I = np.where(self.mean_intensity <= depth)
        intens = self.intensity
        v = self.velocity

        res = []
        for inte in intens:
            i0 = np.argmin(inte)
            """
            in this simple form there are multiple values on the intensity axis
            this is due to the fact that towards the edges the left ant right part of 
            the intensity are not
            monotonic functions of the the velocity..
            TODO: add diagnostic tool for the minimal depth
            """
            left = interp1d(inte[:i0+1][::-1], v[:i0+1][::-1], assume_sorted=False)
            right = interp1d(inte[i0:], v[i0:], assume_sorted=False)

            res.append((left, right))

        return res


    def rv_bis(self, atdepth=0.5, **kwargs):
        """
        at: 0 -> at minimum  1-> at 1
        """
        bise = self.bisector_borders()
        min_intens = np.min(self.intensity, axis=1)

        at = min_intens + atdepth * (1-min_intens)  # linear interpolation between min_intens and 1
    
        res = ((r(a) + l(a))/2 for (l,r), a in zip(bise, at))
        return np.fromiter(res, dtype=float)

    def vspan(self, uu, ul, lu, ll, nn=32):
        """
        vspan based on upper interval [uu, ul] and lower interval [lu, ll]
        of relative depth...
        """
        ui = np.linspace(uu, ul, nn)
        li = np.linspace(lu, ll, nn)

        vu = np.row_stack([ self.rv_bis(atdepth=d) for d in ui ])
        vl = np.row_stack([ self.rv_bis(atdepth=d) for d in li ])
        data = np.median(vu, axis=0) - np.median(vl, axis=0)
        np.savetxt('vspan_{}.dat'.format(self.name), np.column_stack((self.time, data)))
       
        return data

    def vspan_old(self, upper, lower):
        """
        plotting vspan as function of time
        """
        nn = 20
        tmp = np.zeros(len(self.time))


        du = np.linspace(upper[0], upper[1], nn)
        dl = np.linspace(lower[0], lower[1], nn)

        for t, i in enumerate(self.intensity):
            minval = i.min()

            bu = bisector_v(self.velocity, i, minval + du * (1 - minval))
            bl = bisector_v(self.velocity, i, minval + dl * (1 - minval))
            tmp[t] = np.median(bu) - np.median(bl)

        return tmp


    def lomb_scargel_vspan(self, freqs, uu=0.2, ul=0.4, lu=0.6, ll=0.8):
        rv = self.vspan(uu, ul, lu, ll)
        return spectral.lombscargle(self.time, rv - np.mean(rv), freqs)


    def lomb_scargel_vspan_old(self, freqs, uu, ul, lu, ll):
        rv = self.vspan_old((uu, ul), (lu, ll))
        return spectral.lombscargle(self.time, rv - np.mean(rv), freqs)



    def lomb_skagel_vr_bis(self, freqs, atdepth=0.5):
        rv = self.rv_bis(atdepth)
        return spectral.lombscargle(self.time, rv - np.mean(rv), freqs)

    def done_data_selection(self):
        self._intensity = self.intensity
        self._velocity = self.velocity
        self._errors = self.errors
        self._time = self.time
        self.usedindex = np.full(self._intensity.shape[0], True)
        self.vrange = np.full(self._intensity.shape[1], True)

        return self


class tobemodified(SpectralAnalyser):

    def _Fsystem(self, freq, nharm=1):
        """
        returns a system matrix for a given frequency freq
        and times t
        """
        t = self.time

        n = t.size
        F = np.zeros((n, 2 * nharm))
        for i in range(nharm):
            F[:, 2 * i  ] = np.cos (2 * (i + 1) * np.pi * freq * t)
            F[:, 2 * i + 1] = np.sin (2 * (i + 1) * np.pi * freq * t)

        return F

    def _FindividualConst(self):
        F = np.zeros((self.nobs, 5))
        for i in range(5):
            F[self.list_index[i], i] = 1
        return F

    def rv_eqwidth(self, relative_depth):
        """
        equivalent widht of profile
        """
        mv = self.meanIntensity.min()
        depth = mv + (1 - mv) * relative_depth
        I = np.where(self.meanIntensity > depth)[0]

        vv = self.velocity
        res = 1 - self.intensity
        res[:, I] = 0
        v = np.sum(vv[np.newaxis, :], axis=1)
        v /= np.sum(res, axis=1)
        return v

    def eqwidth(self):
        return np.sum(1-self.intensity,axis=1)

    def meansignoise(self):
        return np.mean(1./self.signoise, axis=1)

    def rv_mean(self, relative_depth):
        """
        mean position of spectrum
        """
        mv = self.meanIntensity.min()
        depth = mv + (1 - mv) * relative_depth
        I = np.where(self.meanIntensity > depth)[0]

        vv = self.velocity
        res = 1 - self.intensity
        res[:, I] = 0
        v = np.sum(res * vv[np.newaxis, :], axis=1)
        v /= np.sum(res, axis=1)
        return v

    def rv_std(self, relative_depth):
        """
        stdv  of spectrum
        """
        res = 1 - self.intensity
        mv = self.meanIntensity.min()
        depth = mv + (1 - mv) * relative_depth
        I = np.where(self.meanIntensity > depth)[0]

        res[:, I] = 0
        vv = self.velocity
        ss = np.sum(res, axis=1)
        mn = np.sum(res * vv[np.newaxis, :], axis=1) / ss

        v = np.sum(res * (vv[np.newaxis, :] - mn[:, np.newaxis]) ** 2 , axis=1) / ss
        v = np.sqrt(v)
        return  v

    def rv_skew(self, relative_depth):
        """
        skewness  of spectrum
        """
        res = 1 - self.intensity
        mv = self.meanIntensity.min()
        depth = mv + (1 - mv) * relative_depth
        I = np.where(self.meanIntensity > depth)[0]


        res[:, I] = 0
        vv = self.velocity
        ss = np.sum(res, axis=1)
        mn = np.sum(res * vv[np.newaxis, :], axis=1) / ss
        vr = np.sum(res * (vv[np.newaxis, :] - mn[:, np.newaxis]) ** 2  , axis=1) / ss
        sk = np.sum(res * (vv[np.newaxis, :] - mn[:, np.newaxis]) ** 3  , axis=1) / ss
        sk /= vr ** 1.5
        return  sk

    def rv_corr(self, relative_depth, use='mean'):
        """
        correlation based vrad
        """

        if 'first' == use:
            base = self.intensity[0]
        elif 'mean' == use:
            base = self.meanIntensity
        else:
            pass

        mv = self.meanIntensity.min()
        depth = mv + (1 - mv) * relative_depth
        I = np.where(self.meanIntensity <= depth)[0]

        inte = self.intensity[:, I]
        m = np.ones(I.size)
        mask = np.correlate(m, m, mode='full')

        # computing factor
        eps = 0.001
        v0 = self.velocity - 0.5 * self.deltavbin
        v1 = v0 - eps * self.deltavbin
        ii0 = self.mean_spectrum_interp3(v0)[I]
        ii1 = self.mean_spectrum_interp3(v1)[I]
        tmp = np.correlate(1 - ii1, 1 - ii0, mode='full')

        n = 3  # offset
        i0 = np.argmax(tmp)
        print(i0)
        # poly = np.polyfit(np.arange(-n, n+1), tmp[i0-n:i0+n+1], deg=2)
        poly = np.polyfit(np.arange(-n, n + 1), tmp[i0 - n:i0 + n + 1], deg=2)
        pos = -0.5 * poly[1] / poly[0]

        print(pos, poly)
        factor = self.deltavbin * eps / pos
        print(factor, self.deltavbin)
        res = np.zeros(self.nobs)
        for r, i in enumerate (inte):
            # i = ii
            tmp = np.correlate(1 - i, 1 - base[I], mode='full')
            # tmp /= mask
            # plt.plot(self.velocity, tmp)

        # plt.show()

            i = np.argmax(tmp)
            poly = np.polyfit(np.arange(-n, n + 1), tmp[i - n:i + n + 1], deg=2)

            res[r] = -0.5 * poly[1] / poly[0] + (i - i0)

        res -= np.mean(res)
        return factor * res

    def vspan(self, upper, lower):
        """
        plotting vspan as function of time
        """
        nn = 20
        tmp = np.zeros(len(self.time))


        du = np.linspace(upper[0], upper[1], nn)
        dl = np.linspace(lower[0], lower[1], nn)

        for t, i in enumerate(self.intensity):
            minval = i.min()

            bu = bisector_v(self.velocity, i, minval + du * (1 - minval))
            bl = bisector_v(self.velocity, i, minval + dl * (1 - minval))
            tmp[t] = np.median(bu) - np.median(bl)

        return tmp

    def mesurementfunktional_vrad(self, upper, lower):
        """
        comptuing the measuremnt functional of vrad
        """
        eps = 0.01
        nn = 20
        tmp = np.zeros(len(self.velocity))


        du = np.linspace(upper[0], upper[1], nn)
        dl = np.linspace(lower[0], lower[1], nn)

        intens = self.meanIntensity
        minval = intens.min()
        bu = bisector_v(self.velocity, intens, minval + du * (1 - minval))
        bl = bisector_v(self.velocity, intens, minval + dl * (1 - minval))

        phi0 = np.mean(bu) - np.mean(bl)

        for v in range( len(self.velocity)):
            i = +intens
            i[v]+= eps
            minval = i.min()

            bu = bisector_v(self.velocity, i, minval + du * (1 - minval))
            bl = bisector_v(self.velocity, i, minval + dl * (1 - minval))
            tmp[v] = ((np.mean(bu) - np.mean(bl)) - phi0) / eps

        print("done", tmp)
        return tmp

    def vrad_bis(self, extension):
        """
        calculates vrad as the median of part of bisector
        """
        nn = 20
        tmp = np.zeros(len(self.time))


        dl = np.linspace(extension[0], extension[1], nn)

        for t, i in enumerate(self.intensity):
            minval = i.min()

            bl = bisector_v(self.velocity, i, minval + dl * (1 - minval))
            tmp[t] = np.median(bl)

        return tmp

    def bisector(self, upper, lower, nn=128):
        tmp = np.zeros((len(self.time), nn))
        d = np.linspace(upper, lower, nn)

        for t, i in enumerate(self.intensity):
                minval = i.min()

                b = bisector_v(self.velocity, i, 1 - d * (1 - minval))
                tmp[t] = b

        return tmp, d



#       def PosteriorFreq(self, fre, dat, nharm=1):
#           """
#           posterior for a sinusoidal wave
#           and individual nightly constant
#           """
#
#           res = np.zeros(fre.size)
#           for i, f in enumerate (fre):
#               # F = np.column_stack( (self._FindividualConst(), self._Fsystem(f, nharm)))
#               F = np.column_stack((np.ones(len(self.time)), self._Fsystem(f, nharm)))
#               calc = LM(F, dat)
#               res[i] = -calc.getLogR()  # calc.getLogDetF()#calc.getLogLHMarginalFixedMarginalSigma()
#               # res -= res.max()
#           return res  # np.exp(res)

    def spectrum_matrix_full(self, time, spec, period=1, nphase=128, method=None):
        """
        binns the spec into nphase bins
        """
        bins = np.linspace(0, period, nphase + 1)
        res = np.zeros((nphase, self.nvelocity))
        ii = np.digitize(np.mod(time, period), bins)
        mask = np.zeros( (nphase, self.nvelocity), dtype = 'bool')
        mask[:,:] = False
        for i, s in zip(range(1,nphase+1), spec):
            I = np.where( ii == i)[0]
            if len(I)>0:
                mask[i-1,:] = True
                res[i-1, :] = method(spec[I,:], axis=0)
        return res, bins[:-1], mask

if __name__ == '__main__':
    narval = SpectralAnalyser('../data/Vega_Narval_2018_031.json')
    sophie = SpectralAnalyser('../data/Vega_2018_0310.json')

    # --------- narval ------------
    dataset = narval

    dataset.outlier_removal(0.003)
    F = np.array([1,2,3,4,3,2,1]); F = F / np.sum(F)
    quantity = dataset.filtered_intensity(F)
    quantity = 1-quantity
    quantity = quantity / np.sum(quantity, axis=1)[:, np.newaxis]
    quantity = quantity - np.mean(quantity, axis=0)[np.newaxis, :]
    binnedspec, bins = spectrum_matrix(
        time=dataset.time,
        nphase=128,
        quantity=quantity
    )
    plt.figure(figsize=(10,4))
    plt.title('narval')
    plt.imshow(np.sign(binnedspec)*np.abs(binnedspec)**0.5)

    #------ sophie ----------

    dataset = sophie

    I = dataset.used_indices_of_night(0)
    dataset.outlier_removal(0.003)
    F = np.array([1,2,3,4,3,2,1]); F = F / np.sum(F)
    quantity = dataset.filtered_intensity(F)
    quantity = 1-quantity
    quantity = quantity / np.sum(quantity, axis=1)[:, np.newaxis]
    quantity = quantity - np.mean(quantity, axis=0)[np.newaxis, :]
    binnedspec, bins = spectrum_matrix(
        time=dataset.time,
        nphase=128,
        quantity=quantity
    )
    plt.figure(figsize=(10,4))
    plt.title('sophie')
    plt.imshow(np.sign(binnedspec)*np.abs(binnedspec)**0.5)
    plt.show()
