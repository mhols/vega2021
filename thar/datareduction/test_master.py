from master import *
from spectrograph import *
from settings import *
import matplotlib.pyplot as plt
import pickle

if RECOMPUTE_2D_POLYNOMIAL:
    with open('ccds.pickle', 'wb') as f:
        pickle.dump(ccds, f)
else:
    with open('ccds.pickle', 'rb') as f:
        ccds = pickle.load(f)
        

voie1 = 0
voie2 = 1


c = ccds[0]['ccd'][voie1]
c.kwargs = kwargs

def plot1():

    plt.figure(figsize=(10, 6))
    for o in c.all_order():
        I = c.index_order(o)
        plt.plot(c.ol[I], c.x[I], 'o', color=c.color_of_order(o))
    plt.show()

def plot2():
    plt.figure(figsize=(10,6))
    p = c.fit_global_polynomial()
    for o in c.all_order():
        I = c.index_order(o)
        plt.plot(c.ol[I], c.x[I]-p(c.ol[I]), '-', 
        color=c.color_of_order(o))
        plt.plot(c.ol[I], c.x[I]-p(c.ol[I]), '.', 
        color=c.color_of_order(o))
    plt.show()

def plot3():
    plt.figure(figsize=(10,6))
    c.fit_polynomial_order_by_order()
    for o in c.all_order():
        p = c.polynomial_fit[o]
        I = c.index_order(o)
        plt.plot(c.ol[I], c.x[I]-p(c.ol[I]), 'o', 
        color=c.color_of_order(o), label=str(o))
    plt.show()
    
def plot4():
    plt.figure(figsize=(10,6))
    c.fit_2d_polynomial()
    for o in c.all_order():
        p = c.polynomial_fit[o]
        I = c.index_order(o)
        plt.plot(c.ol[I], c.x[I]-p(c.ol[I]), 'o', 
        color=c.color_of_order(o), label=str(o))
    plt.show()

def plot5():
    plt.figure()
    olrange = np.linspace(c.ol.min(), c.ol.max(), 2048)
    pglobal = c.get_global_polynomial(full=True)
    for o in c.all_order():
        p = c._map_1D_x_ol_o[o]
        I = c.index_order(o)
        plt.plot(c.ol[I], c.x[I]-pglobal(c.ol[I]), '.', 
        color=c.color_of_order(o))
        plt.plot(olrange, (p(olrange)-pglobal(olrange)), '-', color=c.color_of_order(o))
    plt.ylim([-40, 60])
    plt.show()

def plot6():
    plt.figure()
    olrange = np.linspace(c.ol.min(), c.ol.max(), 2048)
    pglobal = c.get_global_polynomial(full=True)
    for o in c.all_order():
        p = c._map_2D_x_ol_o[o]
        I = c.index_order(o)
        plt.plot(c.ol[I], c.x[I]-pglobal(c.ol[I]), '.', 
        color=c.color_of_order(o))
        plt.plot(olrange, (p(olrange)-pglobal(olrange)), '-', color=c.color_of_order(o))
    plt.ylim([-40, 60])
    plt.show()

def plot7():
    plt.figure()
    xrange = np.arange(NROWS)
    pglobal = c.get_global_polynomial(full=True)
    for o in c.all_order():
        l = c._map_1D_l_x_o[o]
        I = c.index_order(o)
        plt.plot(c.ol[I], c.x[I]-pglobal(c.ol[I]), '.', 
        color=c.color_of_order(o))
        plt.plot(l(xrange)*o, ((olrange)-pglobal(olrange)), '-', color=c.color_of_order(o))
    plt.ylim([-40, 60])
    plt.show()

def palette():

    plt.figure()
    for o in c.all_order():
        plt.plot([o], [o], 'o', color=c.color_of_order(o))
