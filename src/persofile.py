import spectralutil as sp
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

sophie = sp.SpectralAnalyser('sophie_reduced.json')

velocity = sophie.velocity()
intens = sophie.intensity()
isp = UnivariateSpline(velocity, intens[3], s=0.00001)
vv = sp.np.linspace(velocity[0], velocity[-1], 1024)
plt.plot(velocity, intens[3], 'ro')
plt.plot(vv, isp(vv))
plt.show()
