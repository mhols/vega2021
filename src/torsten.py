from scipy.interpolate.fitpack2 import UnivariateSpline
import spectralutil as sp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys




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
        
def ls_spec_vspan(spec, uu=0.2, ul=0.3, lu=0.5, ll=0.6):
    name = 'ls_spec_vspan'
    cpdmin, cpdmax = 0.2, 20
    oms = 2*np.pi * np.linspace(cpdmin, cpdmax, 1024)
        
    data = np.abs(spec.lomb_scargel_vspan(oms, uu, ul, lu, ll))
    data /= np.max(data)
    plt.figure(figsize=(10,4))
    plt.title('Lomb Scargle periodogram of vspan')
    plt.xlabel('Frequency (c/d)')
    plt.ylabel('Relative power spectral density')
    plt.plot(oms/(2*np.pi), data)
    maxa = np.max(data)
    box = [0,15,0,1]
    #self._plt_ls(data, box=box)
        
    #a = plt.axes([.4, .4, .4, .4], axisbg='white')
    #plt.title('window function')
    #wf = self.window
    #wf /= wf.max()
    #plt.plot(self.cycles_per_day, wf, 'b-')
    #plt.xlim(0,7)
    #plt.ylim(0,1)
        
    plt.savefig(name + '.pdf')

    #plt.figure(figsize=(10,4))
    #plt.title('Lomb Scargle periodogram of vspan')
    #plt.xlabel('Frequency (c/d)')
    #plt.ylabel('Relative power spectral density')
#        self._plt_ls_vr(1000*self.ls_vrad_mean)

    #box = [15,50,0,0.04]
    #self._plt_ls(data,box=box, bars=False)
    #plt.savefig(name + "_2"+self.format)



if __name__ == '__main__':
    matplotlib.rcParams.update({'font.size': 10})
    # this need to to run only once
    # after the first round comment it out (TODO in own preparation module)
    # preparing_sophie2012()
    # preparing_sophie2018()
    # preparing_narval2018()

    # work_from_selected_data()
    
    sophie2018 = sp.SpectralAnalyser('sophie_reduced.json')
    sophie2012 = sp.SpectralAnalyser('sophie12_reduced.json')
    narval = sp.SpectralAnalyser("narval_reduced.json")

   # for f in np.linspace(0.9999, 1.0001, 21):
    #    binned_spectrum(sophie2012, sophie2018, period=f*sp.VEGAPERIOD)

    # lomb_scargel(narval, sophie2018, sophie2012, atdepth=0.6)
    # lomb_scargel_vspan(sophie2012, sophie2018, narval, depth=0.99, uu=0.5, ul=0.34, lu=0.25, ll=0.1)


    # lomb_scargel_vspan_old(sophie2012, sophie2018, narval, uu=0.5, ul=0.34, lu=0.25, ll=0.1)
    # radial_velocity(narval, sophie2018, sophie2012, relative_depth=0.8)
    # radial_velocity_correlation(narval, sophie2018, sophie2012, relative_depth=0.8)
    # radial_velocity_bisector(narval, sophie2018, sophie2012, depth=0.9, atdepth=0.5)
    
    ls_spec_vspan(sophie2018)

    plt.show()
