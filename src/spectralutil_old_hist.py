'''
Created on May 22, 2014

@author: hols
'''
#from dycos.lmm import LM
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import json
import os

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
    data = data[data[:,0].argsort()]
    # data = tmp[:100]
    time = data[:, coltime:colspec].ravel()
    print(data.shape, time.shape, time[0].shape)

    #time = time-int(time[0])
    time = time-2458331.

    velocity = data[0, colspec:colval]  # velocities of bins
    intens = data[:, colval:colvul]  # intensities
    tmp = +intens
    for i in range(intens.shape[0]):
        intens[i] = np.convolve(intens[i], np.array([1,2,3,4,3,2,1])/16, 'full')[3:-3]
#        intens[i] = np.convolve(intens[i], np.array([1,2,3,2,1])/9, 'full')[2:-2]
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
#    plt.plot(np.std(diff[I, :], axis=1))

#    plt.show()

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


def bisector(velocity, intensity, depth):
    """
    computes the bisector of the
    absorption (!!) profile velocity, intensity
    at the point depth
    """

    n = intensity.size # number of bins

    i = np.min(np.where(intensity < depth))
    j = np.max(np.where(intensity < depth))

    if (i > 0) and (j < n - 1):
        x = velocity[i - 1] + (depth - intensity[i - 1]) * (velocity[i] - velocity[i - 1]) / (intensity[i] - intensity[i - 1])
        y = velocity[j + 1] + (depth - intensity[j + 1]) * (velocity[j] - velocity[j + 1]) / (intensity[j] - intensity[j + 1])
        return (x + y) / 2.
    else:
        return -10000

def bisector_v(velocity, intensity, depth):
    """
    the bisector for the values in depth
    """
    res = np.zeros(len(depth))
    for i in range(len(depth)):
        res[i] = bisector(velocity, intensity, depth[i])
    return res


def remove_night_trend(time, value):
    """
    removing individual trends for each night
    """
    time -= int(time[0])
    nnight = int(time[-1] - time[0]) + 1
    t0 = int(time[-1])
    for n in xrange(nnight):
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
    for i in xrange(m):
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

    for i in xrange(len(time)):
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

    for i in xrange(len(time)):
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
    intensity = data[:, nval+1:2*nval+2]  # intensities
    errors = data[:, 2*nval+2:]
    name, ext = os.path.splitext(DATAFILE)
    res = {
        'name': name + '.json',
        'description': "a little info",
        'noiselevel': 0.7,
        'nvals': nval,
        'range': [72, 201],
        'time': time.tolist(),
        'velocity': velocity.tolist(),
        'intensity': intensity.tolist(),
        'errors': errors.tolist()
    }
    with open(os.path.join(os.path.dirname(DATAFILE), name+".json"), 'w') as outfile:
        json.dump(res, outfile, indent=2)





class SpectralAnalyser:
    """
    a class for the modeling of spectral lines a la Boehm
    """


    def __init__(self,
                DATAFILE,     # datafile containing the filematrix
                nval,     # number of velocity bins
                normalise,    # shall the spectra be normlized True or False
                rangei,       # intensity ranges
                vrange,       # velocity ranges
                noiselevel    # variance of noise TODO CHECK std ?
        ):

        # reading the file matrix

        self.time, self.velocity, self.intensity, self.signoise, \
        self.list_time, self.list_inte, self.list_index = \
        load_data(DATAFILE, nval, rangei, vrange, noiselevel)

        # normalize the spectrum if required TODO
        if normalise:
            fac = 1. / np.max(self.intensity)
            self.intensity = 1 - fac * (1 - self.intensity)

        self.nobs, self.nvelocity = self.intensity.shape
        self.intshape = (self.nobs, self.nvelocity)

        self.meanIntensity = np.mean(self.intensity, axis=0)

        self.variation = self.intensity - self.meanIntensity[np.newaxis, :]
        self.deltavbin = self.velocity[2] - self.velocity[1]

        self._tck = interpolate.splrep(self.velocity, self.meanIntensity, s=0)

    def export_json(self):
        res = {
            'name': 'VEGA_384.json',
            'description': "This data file is based on Boehm",

            'time': self.time.tolist()
        }
        with open("data.json", 'w') as outfile:
            json.dump(res, outfile)

    def mean_spectrum_interp3(self, v):
        """
        cubic spline interpolation of spectrum
        """
        return interpolate.splev(v, self._tck, der=0)

    def mean_spectrum_interp3(self, v):
        """
        cubic spline interpolation of spectrum
        """
        return interpolate.splev(v, self._tck, der=0)


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

    def spectrum_matrix(self, time, spec, period=1, nphase=128):
        """
    binns the spec into nphase bins
    """
        bins = np.linspace(0, period, nphase + 1)
        res = np.zeros((nphase, self.nvelocity))
        nn = np.zeros(nphase)
        ii = np.digitize(np.mod(time, period), bins)
        for i, s in zip(ii, spec):
            #i = i%nphase
            res[i-1, :] += s
            nn[i-1] +=1.
        nn = np.where(nn>0,nn,1)
        res[:,:] = nn[:,np.newaxis]#res/nn[:,np.newaxis]
        #res[:,:] = nn[:,np.newaxis]
        return res, bins[:-1]

    def spectrum_matrix_full(self, time, spec, period=1, nphase=128, method=None):
        """
    binns the spec into nphase bins
    """
        #print("time",time)
        bins = np.linspace(0, period, nphase + 1)
        #print("bins",bins)
        #print("nphase",nphase)
        res = np.zeros((nphase, self.nvelocity))
        ii = np.digitize(np.mod(time, period), bins)
        #print("ii:",ii)
        mask = np.zeros( (nphase, self.nvelocity), dtype = 'bool')
        mask[:,:] = False
        for i, s in zip(ii, spec):
            #print("ii,i",ii,i)
            I = np.where( ii == i)[0]
            if len(I)>0:
                mask[i-1,:] = True
                res[i-1, :] = method(spec[I,:], axis=0)
        return res, bins[:-1], mask

if __name__ == '__main__':
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import spectral

    PROGDIR = os.path.dirname(__file__)
    DATADIR = os.path.join(PROGDIR, '../data')

    DATAFILE = os.path.join(DATADIR, 'Vega_2018_0310.dat')
    # DATAFILE = os.path.join(DATADIR, 'filematrix.dat')


    myanal = SpectralAnalyser( DATAFILE,     # datafile containing the filematrix
        nval=201,     # number of velocity bins
        normalise=False,    # shall the spectra be normlized True or False
        rangei=(72, 112),       # intensity ranges
        vrange=(-60, 40),       # velocity ranges
        noiselevel=0.5    # variance of noise TODO CHECK std ?
    )
    myanal.export_json()

