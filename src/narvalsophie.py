from scipy.interpolate.fitpack2 import UnivariateSpline
import spectralutil as sp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def working_from_rawdata():
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
    
    # exporting reduced data to json
    narval.export_json(outfile='narval_reduced.json') #  exporting reduced data to json)
    # presenting the selected data
    nights = [narval.indices_of_night(i) for i in range(nnights)]
    used_nights = [narval.used_indices_of_night(i) for i in range(nnights)]
    time = narval._time
    stds = narval.std_normalized_over_time()  # std_over_time()
    for i, (I, II) in enumerate(zip(nights, used_nights)):
        plt.figure()
        plt.title('night narval' + str(i))
        t0 = time[I][0]
        plt.plot(I, stds[I])
        plt.plot(II, stds[II], 'o')
        plt.ylim(0.*eps, 1.5*eps)
    
    F = np.array([1,2,3,4,3,2,1]); F = F / np.sum(F)
    quantity = narval.filtered_intensity(F)
    quantity = 1-quantity
    quantity = quantity / np.sum(quantity, axis=1)[:, np.newaxis]
    quantity = quantity - np.median(quantity, axis=0)[np.newaxis, :]
    binnedspec_narval, bins = sp.spectrum_matrix(
        time=narval.time(),
        nphase=256,
        quantity=quantity
    )
    
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
    
    # exporting reduced data as json
    sophie.export_json(outfile='sophie_reduced.json') #  exporting reduced data to json)
    
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
    
    F = np.array([1,2,1]); F = F / np.sum(F)
    quantity = sophie.filtered_intensity(F)
    quantity = 1-quantity
    quantity = quantity / np.sum(quantity, axis=1)[:, np.newaxis]
    quantity = quantity - np.median(quantity, axis=0)[np.newaxis, :]
    binnedspec_sophie, bins  = sp.spectrum_matrix(
        time=sophie.time(),
        nphase=256,
        quantity=quantity
    )
    
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
    
    # exporting reduced data as json
    sophie12.export_json(outfile='sophie12_reduced.json') #  exporting reduced data to json)
    #
    
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
    
    F = np.array([1,2,1]); F = F / np.sum(F)
    quantity = sophie12.filtered_intensity(F)
    quantity = 1-quantity
    quantity = quantity / np.sum(quantity, axis=1)[:, np.newaxis]
    quantity = quantity - np.median(quantity, axis=0)[np.newaxis, :]
    binnedspec_sophie2012, bins  = sp.spectrum_matrix(
        time=sophie12.time(),
        nphase=256,
        quantity=quantity
    )
    
    
    
    
    plt.figure(figsize=(10,10))
    plt.title('sophie narval')
    plt.subplot(121)
    plt.title('sophie 2018')
    plt.imshow(np.sign(binnedspec_sophie)*np.abs(binnedspec_sophie)**0.5,
               aspect='auto', origin='lower', cmap='gist_gray',
                extent=[
                    np.min(sophie.velocity()),
                    np.max(sophie.velocity()),
                    0, 1]
    )
    plt.axvline(x=sophie.v0())
    plt.axvline(x=sophie.v0()+22)
    plt.axvline(x=sophie.v0()+37)
    plt.axvline(x=sophie.v0()-22)
    plt.axvline(x=sophie.v0()-37)
    
    plt.subplot(122)
    plt.title('narval 2018')
    plt.imshow(np.sign(binnedspec_narval)*np.abs(binnedspec_narval)**0.5,
               aspect='auto', origin='lower', cmap='gist_gray',
               extent=[
                    np.min(narval.velocity()),
                    np.max(narval.velocity()),
                    0, 1]
    )
    plt.axvline(x=narval.v0())
    plt.axvline(x=narval.v0()+22)
    plt.axvline(x=narval.v0()+37)
    plt.axvline(x=narval.v0()-22)
    plt.axvline(x=narval.v0()-37)
    
    plt.figure(figsize=(10,10))
    plt.title('sophie2012 sophie2018')
    plt.subplot(121)
    plt.title('sophie 2018')
    plt.imshow(np.sign(binnedspec_sophie)*np.abs(binnedspec_sophie)**0.5,
               aspect='auto', origin='lower',      cmap='gist_gray',
               extent=[
                    np.min(sophie12.velocity()),
                    np.max(sophie12.velocity()),
                    0, 1]
     )
    plt.axvline(x=sophie.v0())
    plt.axvline(x=sophie.v0()+22)
    plt.axvline(x=sophie.v0()+37)
    plt.axvline(x=sophie.v0()-22)
    plt.axvline(x=sophie.v0()-37)
    plt.subplot(122)
    plt.title('sophie 2012')
    plt.imshow(np.sign(binnedspec_sophie2012)*np.abs(binnedspec_sophie2012)**0.5,
               cmap='gist_gray', aspect='auto', origin='lower',
               extent=[
                    np.min(sophie12.velocity()),
                    np.max(sophie12.velocity()),
                    0, 1]
               )
    
    plt.axvline(x=sophie12.v0())
    plt.axvline(x=sophie12.v0()+22)
    plt.axvline(x=sophie12.v0()+37)
    plt.axvline(x=sophie12.v0()-22)
    plt.axvline(x=sophie12.v0()-37)
    plt.show()

def work_from_selected_data():
    sophie2018 = sp.SpectralAnalyser('sophie_reduced.json')
    sophie2012 = sp.SpectralAnalyser('sophie12_reduced.json')
    narval = sp.SpectralAnalyser("narval_reduced.json")

    specdata = sophie2018

    colors = ['r', 'g', 'b', 'k', 'm', 'y']
    

    for factor in np.linspace(0.99, 1.05, 15):
        plt.figure()
        F = np.array([1,2,1]); F = F / np.sum(F)
        quantity = specdata.filtered_intensity(F)
        quantity = quantity[specdata.used_indices_of_night(1)]
        quantity = 1-quantity
        quantity = quantity / np.sum(quantity, axis=1)[:, np.newaxis]
        quantity = quantity - np.median(quantity, axis=0)[np.newaxis, :]
        binnedspec_sophie, bins  = sp.spectrum_matrix(
            time=sophie2018.time(),
            nphase=256,
            quantity=quantity,
            period = factor * sp.VEGAPERIOD
        )
        plt.title('sophie 2018')
        plt.imshow(np.sign(binnedspec_sophie)*np.abs(binnedspec_sophie)**0.5,
                   cmap='gist_gray', aspect='auto', origin='lower',
                   extent=[
                        np.min(sophie2018.velocity()),
                        np.max(sophie2018.velocity()),
                        0, 1]
                   )
        for ofs in [0, 22, 37, -22, -37]: 
            plt.axvline(x=specdata.v0()+ofs)
           
    """
    for factor in np.linspace(0.95, 1.1, 20):
        plt.figure()
        for n, I  in enumerate(sophie2018.indices_of_used_nights()):
            t = sophie2018.time()[I]
            d = sophie2018.intensity()[I, 41]
            v = UnivariateSpline(t, d-d.mean(), s=0.0001)(t)
    
            period = factor * sp.VEGAPERIOD
    
            plt.plot(
                    np.mod(t, period), 
                    d-np.mean(d), 'o', label='night '+str(n), color=colors[n])
            plt.plot(np.mod(t, period), v, '-', linewidth=5, color=colors[n])
        plt.legend()
    """
    plt.show()


def radial_velocity(*specdat):
    ax = plt.figure()
    name = ''
    for spdat in specdat:
        name += spdat.name()
        for I in spdat.indices_of_used_nights():
            plt.plot(np.mod(spdat.time()[I], sp.VEGAPERIOD), spdat.rv_mean()[I], 'o')

    plt.title(name)

if __name__ == '__main__':
    matplotlib.rcParams.update({'font.size': 22})
    # working_from_rawdata()

    sophie2018 = sp.SpectralAnalyser('sophie_reduced.json')
    sophie2012 = sp.SpectralAnalyser('sophie12_reduced.json')
    narval = sp.SpectralAnalyser("narval_reduced.json")

    # work_from_selected_data()


    radial_velocity(sophie2012)
    radial_velocity(sophie2018)
    radial_velocity(narval, sophie2018)
    plt.show()