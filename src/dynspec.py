'''
Created on Jun 13, 2013

@author: hols,boehm boehm ist natuerlich der Beste


coordinae system convention:

observation line is along X-axis [1,0,0]^T
right handed system

'''


from math import *
import numpy as np
import matplotlib.pyplot as plt
#from scipy import weave

def rotation_matrix(axis, theta, mat = None):
    if mat == None:
        mat = np.eye(3,3)

    x = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2])
    a = cos(theta / 2.0);
    b = -(axis[0] / x) * sin(theta / 2.0);
    c = -(axis[1] / x) * sin(theta / 2.0);
    d = -(axis[2] / x) * sin(theta / 2.0);

    mat[0,0] = a*a + b*b - c*c - d*d;
    mat[0,1] = 2 * (b*c - a*d);
    mat[0,2] = 2 * (b*d + a*c);

    mat[1,0] = 2*(b*c+a*d);
    mat[1,1] = a*a+c*c-b*b-d*d;
    mat[1,2] = 2*(c*d-a*b);

    mat[2,0] = 2*(b*d-a*c);
    mat[2,1] = 2*(c*d+a*b);
    mat[2,2] = a*a+d*d-b*b-c*c;
    support = "#include <math.h>"
    code = """
        double x = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
        double a = cos(theta / 2.0);
        double b = -(axis[0] / x) * sin(theta / 2.0);
        double c = -(axis[1] / x) * sin(theta / 2.0);
        double d = -(axis[2] / x) * sin(theta / 2.0);

        mat[0] = a*a + b*b - c*c - d*d;
        mat[1] = 2 * (b*c - a*d);
        mat[2] = 2 * (b*d + a*c);

        mat[3*1 + 0] = 2*(b*c+a*d);
        mat[3*1 + 1] = a*a+c*c-b*b-d*d;
        mat[3*1 + 2] = 2*(c*d-a*b);

        mat[3*2 + 0] = 2*(b*d-a*c);
        mat[3*2 + 1] = 2*(c*d+a*b);
        mat[3*2 + 2] = a*a+d*d-b*b-c*c;
    """

    #weave.inline(code, ['axis', 'theta', 'mat'], support_code = support, libraries = ['m'])

    return mat

GRAD = pi/180
KM = 1000

class Star(object):
    """
    a class modling the star with respect to
    doppler imaging
    """
    
    def __init__(self, obsdir, ey, ez):
        self._obsdir = obsdir
        self._ey = ey
        self._ez = ez
    
    
    def PositionsOfConstantVrel(self,v,nn):
        """
        For given v returns n positions on the star
        along the line of constant v
        """
        pass
    
    def DifferentialIntensity(self, r):
        """
        the 
        """
        pass
    
    def IntensityVelocity(self, direction, offset, time):
        """
        the intensity of light and the velocity (doppler shift)
        along the line of sight offset + l * dirction, l \in R
        at time t
        """
        pass
    
    def Rmax(self):
        return 1.
    
    def SpectralLine(self, time, nn=517, ns=35):
        """
        The spectral line at time time
        nn: sample density over star disk
        direction : observation direction towards the star
        ey : orthogonal direction
        ey : orthogonal direction
        """
        epsilon = 0.00
        rm = self.Rmax()
        dy = np.linspace(-rm , rm, nn, endpoint=True)[:, np.newaxis, np.newaxis]+np.random.uniform(-epsilon, epsilon, (nn,nn))[:,:,np.newaxis]
        dz = np.linspace(-rm , rm, nn, endpoint=True)[np.newaxis, :, np.newaxis]+np.random.uniform(-epsilon, epsilon, (nn,nn))[:,:,np.newaxis]
        
        dd = self._ey[np.newaxis, np.newaxis, :] *dy + self._ez[np.newaxis, np.newaxis, :] * dz
        
        intens, vel = self.IntensityVelocity(self._obsdir, dd, time)
        hist, bin_edge = np.histogram(vel, bins=ns, weights=intens, density=True)
        return hist
        

class RotatingNonPulsatingStar(Star):
    """
    a rotating star 
    """
    def __init__(self, Omega, obsdir, R=1):
        self._R = R
        self._Om = np.linalg.norm(Omega, 2)
        self._OmegaUnit = Omega / self._Om
        self._obsdir = np.asanyarray(obsdir)
        self._obsdir /= np.linalg.norm(self._obsdir,2)
        self._ey = np.cross(self._OmegaUnit, self._obsdir); self._ey /= np.linalg.norm(self._ey,2)
        self._ez = np.cross(self._ey, self._obsdir); self._ez /= np.linalg.norm(self._ez,2)
        self._v = np.array([ Omega[1]*obsdir[2]-Omega[2]*obsdir[1], 
                             Omega[2]*obsdir[0]-Omega[0]*obsdir[2],
                             Omega[0]*obsdir[1]-Omega[1]*obsdir[0]])
    
    def contrast(self, position):
        """
        the contrast at position at time = 0 (in the non-rotating frame)
        """
        pass

    def vmax(self):
        return sin(np.dot(self._OmegaUnit, self._obsdir))*self._Om*self._R

    def SpectralLine2(self, time, nn=128):
        d = 2*self.vmax() / nn
        vv=np.linspace(-self.vmax()+d/2, self.vmax()-d/2, nn, endpoint=True)
        
        # integrate the intesity over a line 
        pass
    
    def _contrastPoint(self, direction, offset):
        """
        returns the visible point
        we suppose that direction are unit vecotrs and orthogonal to the offset
        """
        assert (abs(direction*offset).sum())<1e-6, "not orthogonal...."
        print offset.shape, offset.size, len(offset)
        offset = offset.reshape((offset.size/3,3))
        print offset.shape
        r2 = (offset*offset).sum(axis=1)
        I = np.where(r2 <= 1)
        print r2[I].shape, offset[I].shape
        point = -np.sqrt(1.0 - r2[I])[:,np.newaxis] * direction
        point  += offset[I]
        return point
    
    
    def IntensityVelocity(self, direction, offset, time):
        """
        for a rotating star this is realizes with the help
        of the rotation of the observational lines
        """
        r = rotation_matrix(self._OmegaUnit, self._Om * time)
        point = self._contrastPoint(np.dot(direction, r.T), np.dot(offset, r.T))
        intens = self.contrast(point)
        print point.shape, intens.shape
        v = (self._v*point).sum(axis=1)
        print "v.shape", v.shape
        return intens, v

    def Rmax(self):
        return 1
    
class Vega(RotatingNonPulsatingStar):
    
    def __init__(self):
        RotatingNonPulsatingStar.__init__(self, Omega=[0,0,1.0], obsdir=[-1.0,0,0])
        
    def contrast(self, r):
        d = np.array([1.,1,1.]); d /= np.linalg.norm(d,2)
        print "d", d
        return 1-(1-(r*d).sum(axis=1))**5

def CUniform(r):
    """
    contrast function at position r
    """
    print r.shape
    return np.ones(r.shape[:-1])



class CBlackSpot(object):
    
    def __init__(self, pos, size):
        self.pos = np.asarray(pos)
        self.size = cos(size)
    
    def __call__(self, r):
        return np.where( np.dot(r,self.pos)>self.size, 0., 1.)



class ContrastSolidRotation(object):
    
    def __init__(self, igrad, Omega, R, C):
        self.i = igrad
        self.Omega = Omega
        self.R = R
        self.C = C
        """
        i: inclination angle (degrees) between line of sight and rotation axis (xz plane)
        Omega: angular velocity
        
        C: function object to compute contrast function in standard frame (rotation axis along z)
        """
        self.rotAxis = np.array([cos(self.i),0,-sin(self.i)])
        self.ey = np.array([0,1,0])
        self.ex = np.array([sin(self.i),0,cos(self.i)])

        
    def vmax(self):
        return self.Omega * self.R * sin(self.i)
    
    def eval(self, t, r):
        
        ax = np.dot(r, self.ex)
        ay = np.dot(r, self.ey)
        az = np.dot(r, self.rotAxis)
        
        om = self.Omega*t
        s = sin(om)
        c = cos(om)
        
        aax = ax * c - ay * s
        aay = ax * s + ay * c
        
        rrot = np.outer(aax,self.ex) + np.outer(aay,self.ey) + np.outer(az,self.rotAxis)
        
        return self.C(rrot)


def testRotation():
    CC = ContrastSolidRotation(igrad = 0, Omega = 1, R=1, C = CUniform(0.2))
    
    
class DynamicSpectrum(object):
    
    def __init__(self, Contrast):
        self.Contrast = Contrast
        
        # setting up the grid for numerical integration
        
        nv = 512
        x = np.linspace(-1+0.5/nv,1-0.5/nv,nv)
        y = np.linspace(-1+0.5/nv,1-0.5/nv,nv)

        
        xv, yv = np.meshgrid(x,y)
        x = xv.ravel()
        y = yv.ravel()
        
        a = 0.1
        xv = cos(a)*x - sin(a)*y
        yv = sin(a)*x + cos(a)*y
        self.I = np.where(xv**2+yv**2 < 1.)
        self.xv=xv[self.I]
        self.yv=yv[self.I]
        
    
        self.rrr = np.zeros((self.xv.size,3))
        self.rrr[:,2] = self.xv
        self.rrr[:,1] = self.yv
        self.rrr[:,0] = np.sqrt(self.rrr[:,0]**2+self.rrr[:,1]**2)
        
        self.vrr = self.xv
        print "vr shape", self.vrr.shape
        self.iii = np.argsort(self.vrr)

        self.ev = np.linspace(-1,1,128)

        # sampling the contrast function over the sphere
        
        
        #plt.figure()
        #plt.plot(self.rrr[:,0], self.rrr[:,1], '.')
    
    def eval(self, t):
        """
        computes velocity profile at time t
        for contrast Object Contrast
        """
        
        ccc = self.Contrast.eval(t, self.rrr)
        ccc = ccc[self.iii]
        vvv = self.vrr[self.iii]
        
        tmp = np.interp(self.ev, vvv, ccc.cumsum(), left=0)
        return tmp[1:]-tmp[:-1]
        
        
        #res = np.zeros(self.ev.size-1)
        #for i in xrange(res.size-1):
        #    I = np.logical_and(vvv > self.ev[i],  vvv<=self.ev[i+1])
        #    res[i] = np.sum(ccc[I])
        #return res

if __name__=='__main__':
    """
    some test program
    """
    
    vega = Vega()
    
    plt.figure()
    for t in np.linspace(0, 6, 1):
        spec = vega.SpectralLine(t)
        plt.plot(spec)
    plt.show()
    
#     T = 0.7
#     a = 20*GRAD
#     C = CBlackSpot( pos = [sin(a), 0, cos(a)], size = 20*GRAD)
#     #C = CUniform
#     vega = ContrastSolidRotation(igrad = 90*GRAD, Omega = 2*pi / T, R = 5000, C = C)
#     dyns = DynamicSpectrum (vega)
#     
#     plt.figure()
#     
#     for t in np.linspace(0.1, 0.1+T, 11)[:-1]:
#         v = dyns.eval(t)
#     
#         #plt.plot(v, '-o')
#         plt.plot(-v)
#     plt.show()
