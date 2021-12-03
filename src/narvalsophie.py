import spectralutil as sp
import matplotlib.pyplot as plt
import numpy as np

# loading narval
narval = sp.SpectralAnalyser('../data/Vega_Narval_2018_031.json')
nnights = narval.number_of_nights()

# downsampling of data
eps = 0.001
narval.velocity_range(np.arange(50, 130))
for i in range(10):
    narval.outlier_removal(eps)

# keeping night 5 entirely
narval.usedindex[narval.indices_of_night(5)] = True

# presenting the selected data
nights = [narval.indices_of_night(i) for i in range(nnights)]
used_nights = [narval.used_indices_of_night(i) for i in range(nnights)]
time = narval._time
stds = narval.std_over_time()
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

# presenting the selected data
nights = [sophie.indices_of_night(i) for i in range(nnights)]
used_nights = [sophie.used_indices_of_night(i) for i in range(nnights)]
time = sophie._time
stds = sophie.std_over_time()
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
#eps = 0.00138
eps = 0.0016

for i in range(10):
    sophie12.outlier_removal(eps)
#eps = 0.00125
eps = 0.0015
for i in range(10):
    sophie12.outlier_removal(eps)
# removing by hand
sophie12.usedindex[[2037, 2038, 1726, 1774, 1800, 1413, 1414, 534, 84, 445 ]] = False
sophie12.usedindex[range(1788, 1800)] = False
sophie12.usedindex[range(1326, 1377)] = False
sophie12.usedindex[range(2490, 2534)] = False
sophie12.usedindex[range(1543, 1587)] = False
sophie12.usedindex[range(857, 1016)] = False
sophie12.usedindex[range(13, 21)] = False
sophie12.usedindex[range(156, 171)] = False
sophie12.usedindex[range(2400, 2490)] = False
sophie12.usedindex[range(1946, 2020)] = False

# presenting the selected data
nights = [sophie12.indices_of_night(i) for i in range(nnights)]
used_nights = [sophie12.used_indices_of_night(i) for i in range(nnights)]
time = sophie12._time
stds = sophie12.std_over_time()
for i, (I, II) in enumerate(zip(nights, used_nights)):
    plt.figure()
    plt.title('night sophie 2012' + str(i))
    t0 = time[I][0]
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




plt.figure(figsize=(10,4))
plt.title('sophie narval')
plt.subplot(121)
plt.title('sophie 2018')
plt.imshow(np.sign(binnedspec_sophie)*np.abs(binnedspec_sophie)**0.5, cmap='gist_gray')
plt.subplot(122)
plt.title('narval 2018')
plt.imshow(np.sign(binnedspec_narval)*np.abs(binnedspec_narval)**0.5, cmap='gist_gray')

plt.figure(figsize=(10,4))
plt.title('sophie2012 sophie2018')
plt.subplot(121)
plt.title('sophie 2018')
plt.imshow(np.sign(binnedspec_sophie)*np.abs(binnedspec_sophie)**0.5, cmap='gist_gray')
plt.subplot(122)
plt.title('sophie 2012')
plt.imshow(np.sign(binnedspec_sophie2012)*np.abs(binnedspec_sophie2012)**0.5, cmap='gist_gray')



plt.show()
