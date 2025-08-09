import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys
from scipy import sparse
#from util import *

def TPXYZ(theta, phi):
    return np.vstack(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta)
        ]
    )

def XYZTP(x,y,z):
    r = np.sqrt( x**2 + y**2 + z**2)
    return np.vstack(
        [
            np.arccos(z),
            np.arctan2(y,x)
        ]
    )

def YZTP(y,z):
    x = np.sqrt(1-y**2 -z**2)
    return XYZTP(x,y,z)

def rotation_matrix_Z(om):
    co = np.cos(om)
    so = np.sin(om)
    return np.array(
        [[co, -so, 0],[so, co, 0],[0,0,1]]
    )

def rotation_matrix_X(om):
    co = np.cos(om)
    so = np.sin(om)
    return np.array(
        [[1, 0, 0],[0, co, -so],[0,so,co]]
    )

def rotation_matrix_Y(om):
    co = np.cos(om)
    so = np.sin(om)
    return np.array(
        [[co, 0, so],[0, 1, 0],[-so, 0, co]]
    )

def rot_Y_Z(omx, omz, p):
    
    v = np.dot(rotation_matrix_Z(omz),p)
    return np.dot( rotation_matrix_Y(omx), v)

class Spot:

    def __init__(self, **kwargs):
        """
        'period' : rotation period
        'angle'  : orientation of star (in rad) 0 is pole on  pi/2 equator on
        'nvelocity' : number of velocitybins
        'vbinmin' : minimal v
        'vbinmax' : maximal v
        'velocitybins' : array of velocity values
        'nphase' : number of phases
        'ntheta' : number of thetas on grid on star
        'nphi'   : number of longitudes on grid on star
        """
        self.kwargs = kwargs
        self.velocitybins = self._vbins() # limits of velcity bins (one more than actual bins)
        self.velocitybinscenters = 0.5 * (self.velocitybins[1:] + self.velocitybins[:-1])

    def angle(self):
        return self.kwargs['angle']

    def omega(self):
        return np.array([np.cos(self.angle()), 0, np.sin(self.angle())])

    def _vbins(self):

        try:
            bins = self.kwargs['velocitybins']
            d = (np.max(bins) - np.min(bins)) / (len(bins)-1)
            bins = np.min(bins)  + d * (np.arange(len(bins)+1) - 0.5)
        except:
            d = (self.kwargs['vbinmax'] - self.kwargs['vbinmin']) / (self.kwargs['nvelocity'] -1 )
            bins = self.kwargs['vbinmin'] + d * (np.arange(self.kwargs['nvelocity']+1) - 0.5)
        

        return bins
    
    @property
    def nvelocity(self):
        return self.velocitybinscenters.shape[0]

    def empty_phasemap(self):
        return np.zeros( (self.kwargs['nphase'], self.nvelocity))


    def orbit_on_star_tp(self, theta, phi, t):
        """
        orbit starting at theta phi on start
        """

        phase = 2 * np.pi * np.mod(t, self.kwargs['period']) / self.kwargs['period']
        return theta*np.ones_like(phase), phi + phase
    
    def orbit_on_star_XYZ_Obs(self, theta, phi, t):
        th, phi = self.orbit_on_star_tp(theta, phi, t)
        xyz = TPXYZ(th, phi)
        R = rotation_matrix_Y(np.pi/2 - self.kwargs['angle'])
        return np.dot(R, xyz)
    
    def velocity(self, y):
        return self.kwargs['v0'] + y * self.kwargs['vmax']
    
    def y(self, v):
        return (v-self.kwargs['v0']) / self.kwargs['vmax']
    
    def times(self):
        return np.linspace(0, self.kwargs['period'], self.kwargs['nphase']+1)[:-1]

    def phis(self):
        return np.linspace(0, 2*np.pi, self.kwargs['nphi']+1)[:-1]
    
    def thetas(self):
        angle = self.kwargs['angle']
        if not self.kwargs.get('full', True):
            return np.linspace(0, np.pi/2 + angle, self.kwargs['ntheta']+1)[1:]
        else:
            return np.linspace(0, np.pi, self.kwargs['ntheta']+2)[1:-1]

    def phasemap_of_point_spot(self, theta, phi):
        """
        image of delta surface spot (the directional cos is taken into account)
        """
        t = self.times()
        
        o = self.orbit_on_star_XYZ_Obs(theta, phi, t)
        v = self.velocity(o[1, :]) # y coordinate is velocity
        visi = np.where(o[0,:]>0, o[0, :], 0)

        res = self.empty_phasemap()
        i = np.digitize(v, self.velocitybins)-1
        
        for k in range(res.shape[0]):
            res[k, i[k]] +=  visi[k]

        return res
    
    def theta_phi_grid(self):
        if not hasattr(self, '_theta_phi_grid'):
            self._theta_phi_grid = np.meshgrid(
                #np.linspace(0, np.pi, self.kwargs['ntheta']+2, endpoint=True)[1:-1], 
                #np.linspace(2*np.pi, 0., self.kwargs['nphi'], endpoint=False)
                self.thetas(), 
                self.phis()
        )
        return self._theta_phi_grid
    
    def matrix_star_wavemap(self):

        if not hasattr(self, '_matrix_star_wavemap'):
            th, ph = self.theta_phi_grid()
            th = th.ravel()
            ph = ph.ravel()
            ti = self.times()

            res = np.zeros((self.kwargs['nphase'] * self.nvelocity, th.shape[0]))

            for i, (tth, pph) in enumerate(zip(th, ph)):
                res[:, i ] = self.phasemap_of_point_spot(tth, pph).ravel() * np.sin(tth)

            self._matrix_star_wavemap = res

        return self._matrix_star_wavemap
    
    def matrix_wavemap_star(self):
        #
        # return self.matrix_star_wavemap().T
        if not hasattr(self, '_matrix_wavemap_star'):
            th, ph = self.theta_phi_grid()
            th = th.ravel()
            ph = ph.ravel()
            ti = self.times()

            res = np.zeros((th.shape[0], self.kwargs['nphase'] * self.nvelocity ) )

            for i, (tth, pph) in enumerate(zip(th, ph)):
                res[i, : ] = np.where( self.phasemap_of_point_spot(tth, pph).ravel() * np.sin(tth) > 0, 1, 0)

            b = np.dot( self.matrix_star_wavemap(), np.ones(th.shape[0]))
            b = np.dot( res, b)

            res = res * (b / (b**2 + 0.1 * np.max(b.ravel()**2)))[:, None]
            #bias = self._bias(res)

            #for i in range(res.shape[0]):
            #    j = i % self.kwargs['ntheta']
            #    res[i, :] /= 1 #bias[j]

            self._matrix_wavemap_star = res

        return self._matrix_wavemap_star
    
   
    def plot_on_sphere(self, value):

        assert self.kwargs.get('full', True), 'plot_on_sphere needs full grid'

        map = Basemap(projection='moll', lon_0=0)
        if len(value.shape)==1:
            value = np.reshape(value, (self.kwargs['nphi'], self.kwargs['ntheta']))

        #dark spot on star = emission relative to average profile, generates bright
        #S-shaped signature. Therefore here we have to plot -value!!
        
        map.imshow(-1e8*value[:, ::-1].T, interpolation='none',vmin=-0.0003,vmax=0.0003)
        #map.imshow(-1e8*value[:, ::-1].T, interpolation='none')
        
        map.drawmeridians(np.linspace(-180, 180, 13))
        map.drawparallels(np.linspace(-90, 90 , 7))

        map.drawparallels([0], color = "red", linewidth=2)
        map.drawparallels([-self.kwargs['angle'] * 180 / np.pi, self.kwargs['angle']* 180 / np.pi], color = "blue", linewidth=2)
        #map.plot(180*phi./np.pi, 128*[0], '-r', latlon=True)
        #plt.show()
        #plt.figure()
        plt.imshow(value[:,28:68], interpolation='none')
        #plt.imshow(value, interpolation='none')
        plt.show()
        plt.plot(np.sum(value[:,28:68],axis=1))
        plt.show()
        np.savetxt("courbe_activity18.txt",np.sum(value[:,28:68],axis=1))

    
    def plot_visible_on_cylinder(self, value):
        assert not self.kwargs.get('full', True), 'plot_visible_on_sphere needs partial grid'

        if len(value.shape)==1:
            value = np.reshape(value, (self.kwargs['nphi'], self.kwargs['ntheta']))
        angle = self.kwargs['angle'] * 180 / np.pi
        plt.imshow(value.T, extent=(0, 360, -angle, 90))
        plt.yticks([])
        plt.yticks( [-angle, 0,angle, 30, 60, 90])
        plt.hlines([0, angle], 0, 360)
        plt.xticks()
        plt.xticks(np.arange(0,360,30))


    def plot_phasemap(self, phasemap):

        if len(phasemap.shape)==1:
            phasemap = self.reshape(phasemap)
        vbins = self.velocitybins
        v0 = self.kwargs['v0']
        plt.xlim(vbins[0], vbins[-1])
        plt.ylim(0, 1)
        plt.imshow(phasemap, aspect='auto',interpolation='none', origin='lower',
                extent=[vbins[0], vbins[-1],0,1])
        plt.xticks([])
        rv = np.array([-20, -10, 0, 10,])
        st = [ str(i) for i in rv]
        plt.xticks(rv+v0, st)
        plt.yticks([])
        yt = [0,0.2,0.4,0.6,0.8, 1]
        plt.yticks(yt,[str(t) for t in yt])
        plt.vlines([-22+v0, v0, 22+v0], 0,1)
        plt.vlines([-37+v0, v0, 37+v0], 0,1,linestyles='dashed')
        plt.hlines([1./2], -33,32)
        plt.xlabel('velocity [km/s]')
        plt.ylabel('phase fraction of period]')

    def extended_spot(self, theta, phi, angle):
        t, p = self.theta_phi_grid()

        center = TPXYZ(theta, phi)
        xyz = TPXYZ(t.ravel(), p.ravel())

        return np.where( np.sum(center * xyz, axis=0) < np.cos(angle), np.random(0,5e-5),  1./70000)
    
    def svd(self):
        return np.linalg.svd(self.matrix_star_wavemap(), full_matrices=False, compute_uv=False, hermitian=False)
    
    def phasemap_from_star(self, values):
        th = self.theta_phi_grid()[0]
        th = th.ravel()

        m  = np.dot(self.matrix_star_wavemap(), values.ravel())
        return s.reshape(m)
    
    def star_from_phasemap(self, wavemap, sigma=1.e-5, niter=1):
        F = self.matrix_wavemap_star()
        w = wavemap.ravel()
        th = self.theta_phi_grid()[0].ravel()
        tmp = sigma * np.dot(F, w)

        return tmp 
    
    def reshape(self, flatfilematrix):
        return np.reshape(+flatfilematrix, (self.kwargs['nphase'], self.nvelocity  ))

#vbins = np.loadtxt('velocitybins.txt')
movingpeaks = -np.loadtxt('moving_simple_per_night.txt')
#v0=-13.01
#GRAD = np.pi/180

s = Spot(
    period=0.7, 
    angle=7*np.pi/180,
    nvelocity=59,
    vbinmin= -58.9,
    vbinmax=38.0,
    #velocitybins=vbins,
    nphase=128, 
    v0=0, #-13.01, 
    vmax=22,
    ntheta=128,
    nphi = 128
    )
 
if __name__=="__main__":
    
    import matplotlib.pyplot as plt

    #vbins = np.loadtxt('velocitybins.txt')
    #movingpeaks = -np.loadtxt('moving_simple_per_night.txt')
    #v0=-13.01
    GRAD = np.pi/180

    s = Spot(
        period=0.7, 
        angle=7*np.pi/180,
        nvelocity=58,
        vbinmin= -58.9,
        vbinmax=38.0,
        full=False,
        #velocitybins=vbins,
        nphase=128, 
        v0=-13.01, 
        vmax=22,
        ntheta=70,
        nphi = 128
        )
 

    spots = s.extended_spot(70*GRAD, -10*GRAD, 10*GRAD)  + s.extended_spot(20*GRAD, -120*GRAD,10*GRAD) + s.extended_spot(50*GRAD, -270*GRAD, 10*GRAD)
    #spots =  s.extended_spot(80*GRAD, -10*GRAD, 10*GRAD) + s.extended_spot(80*GRAD, -30*GRAD, 10*GRAD)
    #spots =  s.extended_spot(80*GRAD, -10*GRAD, 10*GRAD)

    phasemap = s.phasemap_from_star(spots)
    print( s.matrix_star_wavemap().shape )
    #plt.figure()

    #plt.imshow(s.reshape(s.matrix_star_wavemap()[:, 678]))
    ##s.plot_phasemap(s.reshape(s.matrix_star_wavemap()[:, 678]))

    

    
    #plt.figure()

    onepixel = s.phasemap_of_point_spot(70*GRAD, 0)
    #s.plot_phasemap(onepixel)



    plt.figure()
    plt.title("spot on map")
    #s.plot_on_sphere( spots )
    s.plot_visible_on_cylinder( spots )
    
    plt.figure()
    plt.title("dynamical spectrum")
    s.plot_phasemap( s.reshape(phasemap) )
    
    plt.figure()
    plt.title("reconstruction")
    #s.plot_on_sphere(s.star_from_phasemap(phasemap, niter=1))
    s.plot_visible_on_cylinder(s.star_from_phasemap(phasemap))
    #plt.figure()
    #s.plot_on_sphere(s.star_from_phasemap(movingpeaks))

    #plt.figure()
    
    """
    s.plot_on_sphere(s.star_from_phasemap(movingpeaks))

    plt.figure()
    s.plot_on_sphere(np.sum(s.matrix_wavemap_star()**2, axis=1))
    s.plot_visible_on_cylinder(np.sum(s.matrix_wavemap_star()**2, axis=1))
    
    plt.figure()
    pm =  s.phasemap_of_point_spot(45*GRAD, 0) \
        + s.phasemap_of_point_spot(70*GRAD, 100*GRAD) \
        + s.phasemap_of_point_spot(10*GRAD, -10*GRAD)
    s.plot_on_sphere(s.star_from_phasemap(pm))
    """

    #s.plot_on_sphere(s.star_from_phasemap(pm)) 
    #s.plot_visible_on_cylinder(s.star_from_phasemap(phasemap))
    
    #plt.figure()
    #s.plot_visible_on_cylinder(s.star_from_phasemap(movingpeaks))

    plt.show()
