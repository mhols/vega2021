from scipy.interpolate.fitpack2 import UnivariateSpline
import spectralutil as sp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys


def preparing_narval2018(plotting=False):
    # loading narval
    narval = sp.SpectralAnalyser('../data/Vega_Narval_2018_031.json')
    nnights = narval.number_of_nights()
    
    # downsampling of data
    eps = 0.001
    narval.velocity_range(np.arange(50, 130))
    for i in range(10):
        narval.outlier_removal(eps)
    
    # keeping night 5 and  entirely
    narval.usedindex[narval.indices_of_night(5)] = True
    narval.usedindex[narval.indices_of_night(6)] = True
    
   
    if plotting:
        # presenting the selected data
        nights = [narval.indices_of_night(i) for i in range(nnights)]
        used_nights = [narval.used_indices_of_night(i) for i in range(nnights)]
        time = narval._time
        stds = narval.std_normalized_over_time()  # std_over_time()
        for i, (I, II) in enumerate(zip(nights, used_nights)):
            plt.figure()
            plt.title('night narval ' + str(i))
            # t0 = time[I][0]
            plt.plot(I, stds[I])
            plt.plot(II, stds[II], 'o')
            plt.ylim(0.*eps, 1.5*eps)
    
    # exporting reduced data to json
    narval.done_data_selection().export_json(outfile='narval_reduced.json') #  exporting reduced data to json)
 
def preparing_sophie2012(plotting=False):
    # loading sophie 2012
    sophie12 = sp.SpectralAnalyser('../data/Vega_Sophie_2012_9500.40.03-10.json')
    nnights = sophie12.number_of_nights()
    
    # downsampling of data
    sophie12.velocity_range(np.arange(50, 130))
    eps = 0.0015
    for i in range(10):
        sophie12.outlier_removal(eps)
    
    # removing by hand
    sophie12.usedindex[[2021, 2037, 2038, 1726, 1774, 1800, 1413, 1414, 1593, 534, 84, 445 ]] = False
    sophie12.usedindex[range(1788, 1800)] = False
    sophie12.usedindex[range(1326, 1377)] = False
    sophie12.usedindex[range(2490, 2534)] = False
    sophie12.usedindex[range(1543, 1587)] = False
    sophie12.usedindex[range(857, 1016)] = False
    sophie12.usedindex[range(13, 21)] = False
    sophie12.usedindex[range(156, 171)] = False
    sophie12.usedindex[range(2400, 2490)] = False
    sophie12.usedindex[range(1946, 2020)] = False
    sophie12.usedindex[range(1002, 1021)] = False
    
    
    if plotting:
        # presenting the selected data
        nights = [sophie12.indices_of_night(i) for i in range(nnights)]
        used_nights = [sophie12.used_indices_of_night(i) for i in range(nnights)]
        time = sophie12._time
        stds = sophie12.std_over_time()
        for i, (I, II) in enumerate(zip(nights, used_nights)):
            plt.figure()
            plt.title('night sophie 2012' + str(i))
            plt.plot(I, stds[I])
            plt.plot(II, stds[II], 'o')
            plt.ylim(0.*eps, 1.5*eps)
    
    # exporting reduced data as json
    sophie12.done_data_selection().export_json(outfile='sophie12_reduced.json') #  exporting reduced data to json)
    #
    

def preparing_sophie2018(plotting=False):
    # loading sophie
    sophie = sp.SpectralAnalyser('../data/Vega_2018_0310.json')
    nnights = sophie.number_of_nights()
    
    # downsampling of data
    sophie.velocity_range(np.arange(50, 130))
    eps = 0.00138
    for i in range(10):
        sophie.outlier_removal(eps)
    eps = 0.00125
    for i in range(10):
        sophie.outlier_removal(eps)
    # removing by hand
    sophie.usedindex[[289, 686, 2112]] = False
    sophie.usedindex[range(986, 1037)] = False
    
    if plotting:
        # presenting the selected data
        nights = [sophie.indices_of_night(i) for i in range(nnights)]
        used_nights = [sophie.used_indices_of_night(i) for i in range(nnights)]
        time = sophie._time
        stds = sophie.std_normalized_over_time()  # std_over_time()
        for i, (I, II) in enumerate(zip(nights, used_nights)):
            plt.figure()
            plt.title('night sophie' + str(i))
            t0 = time[I][0]
            plt.plot(I, stds[I])
            plt.plot(II, stds[II], 'o')
            plt.ylim(0.*eps, 1.5*eps)
    
    # exporting reduced data as json
    sophie.done_data_selection().export_json(outfile='sophie_reduced.json') #  exporting reduced data to json)

# reduce the data

preparing_sophie2018()
preparing_sophie2012()
preparing_narval2018()

