'''
Created on Jan 2, 2014

@author: hols
'''

import numpy as np
from math import sin, cos, pi

GRAD = pi/180

def Rotation(Omega, Om, p):
    co = np.cos(Om)
    so = np.sin(Om)
    
    #ps = p.shape
    #p=p.reshape((p.size/3, 3))
    print (np.outer(Omega, Omega))
    print (Omega*p).shape
    res = co*p + np.dot(p, np.outer(Omega, Omega)) * (1.-co)
    res += np.cross(Omega, p, axis = -1) * so
    
    return res #co*p + (Omega*p).sum(axis=-1) * (1.-co) * Omega + so * np.cross(Omega, p, axis = -1)
    

def normalize(p):
    p = np.asarray(p)
    return p / np.linalg.norm(p,2)

class RotStar(object):
    """
    A class for the modelling of a rotation star without oscillations
    """
    
    def __init__(self, Omega, Om, obsline):
        self._Om = Om
        self._Omega = normalize(Omega)
        self._obsline = normalize(obsline)
        
        self._v = normalize(np.cross(self._Omega, self._obsline))
        self._u = np.cross(self._obsline, self._v)
        self._g = np.cross(self._v, self._Omega)
        self._nn = 512
        self._ni = 1024
    
    def Intensity(self, p):
        """
        The intensity at point p. North pole is self._Omega
        """
        pass
    
    def CosPhi(self, p):
        return (p * self._g).sum(axis=-1)
    
    def SinPhi(self, p):
        return (p * self._v).sum(axis=-1)
    
    def Phi(self, p):
        return np.arctan2(self.SinPhi(p), self.CosPhi(p))
    
    def CosTheta(self, p):
        return (p * self._Omega).sum(axis=-1)
    
    def Theta(self, p):
        return np.arccos(self.CosTheta(p))
    
    def spectrum(self, time = 0, nn=None, ni = None):
        """
        the spectrum sampled at nn points in velocity and nt points in time
        the integrals are performed using ni nodes
        """
        
        if None != nn:
            self._nn = nn
        else:
            nn = self._nn
            self._ni = ni
        a = np.linspace(-1+0.5/nn,1-0.5/nn, nn) # the velocities
        
        b = np.sqrt(1-a**2) # radius of obseration circle at velocity a
        
        theta = np.linspace(0.5/ni, np.pi - 0.5/ni, ni)
        ct = np.cos(theta)
        st = np.sin(theta)
        ep = a[:,np.newaxis,np.newaxis] * self._v + b[:,np.newaxis,np.newaxis] \
            *(st[np.newaxis,:,np.newaxis]*self._obsline \
            + ct[np.newaxis,:,np.newaxis]*self._u) 
        ep = Rotation(self._Omega, -time * self._Om, ep)
        
        res = b*(self.Intensity(ep)*st[np.newaxis,:]).sum(axis = 1)
        return res
        
class Vega(RotStar):
    """
    Vega as rotating star
    """
    
    def __init__(self):
        RotStar.__init__(self, Omega = [0,0,1], Om = 1, obsline = [sin(7*GRAD),0,cos(7*GRAD)])
        
    def Intensity(self, p):

        phi = self.Phi(p)
        the = self.Theta(p)
        ii = 1+ (1 + 0.0 * np.cos(phi))*np.exp(-((the - pi/4.0)/0.2)**2 )*(the-pi/4.0)
        
        return ii



if __name__=='__main__':
    import matplotlib.pyplot as plt
    vega = Vega()
    
    res = np.zeros((10,512))
    for i, time in enumerate (np.linspace(0, 2*np.pi, 10, endpoint=False)):
        res[i,:] = vega.spectrum(time = time)

    plt.figure()
    for i in range(10):
        plt.plot(res[i,:]- np.mean(res, axis=0))
    plt.show()
    
