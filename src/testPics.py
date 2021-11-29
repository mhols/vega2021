'''
Created on Oct 14, 2014

@author: hols
'''

from picturesForPaper import *

myPics = Pictures(DATAFILE)
#myPics.vrad_mean_vrad_corr()
#    myPics.ts_vrad()
#myPics.intens()
#myPics.vrad_mean_vspan()
#    myPics.vrad_corr_vspan()
#myPics.vrad_mean_skew()
#    myPics.vrad_mean_std()
#myPics.vrad_corr_skew()
#myPics.skew_vspan()
#    myPics.ts_skew()
#myPics.ts_vspan()
#myPics.ls_spec_vrad_skew()
#myPics.ls_spec_vrad_mean()
#myPics.ls_spec_vrad_bis()

#myPics.ls_spec_vrad_corr()
#myPics.ls_spec_vspan()
#myPics.bisector_time()
#myPics.bisector_width()
#myPics.ls_depth_spec_vspan()
#myPics.ls_window()
#alldata = [self.time, self.inte, self.vrad_mean, self.vrad_corr, self.vspan, self.vrad_skew, self.vrad_std]
#myPics.bayes_freq_vrad_mean()
#myPics.saveData("time_vrad_mean.dat", [myPics.time, myPics.vrad_mean])  # first column
#myPics.saveData("time_vrad_corr.dat", [myPics.time, myPics.vrad_corr])  # first column
#myPics.saveData("time_vspan.dat", [myPics.time, myPics.vspan])  # first column
#myPics.saveData("time_skew.dat", [myPics.time, myPics.vrad_skew])  # first column
#myPics.moving_peaks()
#myPics.signoise_eqwidth()
#myPics.moving_peaks_signoise()
#myPics.estrotentropy()
#myPics.mesurementfunktional_vrad()
myPics.measurementfunctions_vrad_perfect()
plt.show()