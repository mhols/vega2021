from master import *
from spectrograph import *
from matplotlib.pyplot import *

import pickle

with open('ccds.pickle', 'wb') as f:
    pickle.dump(ccds, f)

with open('ccds.pickle', 'rb') as f:
    ccds = pickle.load(f)

voie1 = 0
voie2 = 1


c = ccds['/Users/boehm/Desktop/vega2021/thar/datafiles/NEO_20220903_191404_th0.fits']['ccd'][voie1]

figure(figsize=(10, 6))
for o in c.all_order():
    I = c.index_order(o)
    plot(c.ol[I], c.x[I], 'o', color=c.color_of_order(o))


figure(figsize=(10,6))
p = c.fit_global_polynomial()
for o in c.all_order():
    I = c.index_order(o)
    plot(c.ol[I], c.x[I]-p(c.ol[I]), '-', 
    color=c.color_of_order(o))
    plot(c.ol[I], c.x[I]-p(c.ol[I]), '.', 
    color=c.color_of_order(o))


figure(figsize=(10,6))
c.fit_polynomial_order_by_order()
for o in c.all_order():
    p = c.polynomial_fit[o]
    I = c.index_order(o)
    plot(c.ol[I], c.x[I]-p(c.ol[I]), 'o', 
    color=c.color_of_order(o), label=str(o))

figure(figsize=(10,6))
c.fit_2d_polynomial()
for o in c.all_order():
    p = c.polynomial_fit[o]
    I = c.index_order(o)
    plot(c.ol[I], c.x[I]-p(c.ol[I]), 'o', 
    color=c.color_of_order(o), label=str(o))




figure()
for o in c.all_order():
    plot([o], [o], 'o', color=c.color_of_order(o))
