'''
Created on Jun 15, 2013

@author: hols
'''
import sys
import os

import numpy as np
import matplotlib.pyplot as plt

PROGDIR = os.path.dirname(__file__)
DATADIR = os.path.join(PROGDIR, '../data')
DATAFILE = os.path.join(DATADIR, 'filematrix.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_nn.dat')


nval = 201 # number of data points

coltime = 0
colspec = 1
colval  = colspec + nval
colvul  = colval + nval

data = np.loadtxt(DATAFILE)
time = data[:,coltime]
velocity = data[:, coltime:colval]
intens = data[:, colval:colvul]


meani =intens.mean(axis=0)

def doit0():
    plt.figure()
    plt.title('mean value')
    plt.plot(meani)


vari = (intens-meani).ravel().var()

print vari

t0 = time[0]


diff = intens-meani


# outlier removal
I = np.where( diff.var(axis=1) < 0.05*vari)

print "used data: ", len(I[0]), " out of ", len(time)
print I

data = data[I]
time = data[:,coltime]
velocity = data[:, coltime:colval]
intens = data[:, colval:colvul]
meani = intens.mean(axis=0)
diff = intens-meani

rangel = np.arange(50,130) # range of spectral line

def doit1():
    plt.figure()
    plt.title('difference to mean value')
    for d in diff:
        plt.plot(d)

def doit2():
    plt.figure()
    plt.title('difference over time')
    plt.plot(time, (diff**2).sum(axis=1))
    plt.semilogy()


def doit3():
    plt.figure()
    plt.title('distribution of singular values')
    U, D, Vt = np.linalg.svd(intens[:,rangel], compute_uv = True, full_matrices=True)
    plt.plot(D, 'o')
    plt.semilogy()
    
    print U.shape, D.shape, Vt.shape, time.shape
    
    for i in range(4):
        plt.figure()
        plt.title('EOF Nr. %d' %i)
        plt.plot(rangel, Vt[i,:])
        
        plt.figure()
        plt.title('evolution of component %d' %i)
        plt.plot(time, U[:,i])

def doit4():
    maxamp = 30 / np.abs(intens-meani).std()
    for n in range(5): # five nights...
        plt.figure()
        plt.title("night #%d" %(n+1))
        for t, l in zip(time, (intens-meani)[:, rangel]):
            plt.plot(rangel, maxamp*l+10000*(t-t0), '-')
            plt.ylim(10000*n,10000*n+3000) #TODO die Anfaenge der naechte besser verwalten...


def doit5():
    
    ml = data[:, colval+72:colval+87]
    ll = data[:, colspec+72:colspec+87]
    
    mr = data[:, colval+94 :colval+109]
    lr = data[:, colspec+94:colspec+109]
    
    ml = ml / ml.sum(axis=1)[:,np.newaxis]
    mr = mr / mr.sum(axis=1)[:,np.newaxis]
    
    plt.figure()
    for l, r, u, v in zip(ml, mr, ll, lr):
        plt.plot(u,l, '.')
        plt.plot(v, r, '.')
    
    mll = (ll*ml).sum(axis=1)
    mlr = (lr*mr).sum(axis=1)
    
    plt.figure()
    plt.plot(time, mll-mll.ravel().mean(), '.')
    plt.plot(time, mlr-mlr.ravel().mean(), '.')
    
    
    
    
    plt.figure()
    plt.plot((ml*ll).sum(axis=1), (mr*lr).sum(axis=1), 'o')
    

    plt.figure()
    inte = intens[:,50:125]
    inte /= inte.sum(axis=1)[:,np.newaxis]
    velo = (velocity[:,50:125]*inte).sum(axis=1)
    plt.plot(time, velo, '.')


    
    plt.figure()
    pl = np.polyfit(ll[0,:],ml.T,deg=1)
    pr = np.polyfit(lr[0,:],mr.T,deg=1)
    
    
    print pl.shape
    
    plt.plot( pl[0,:], pr[0,:], 'o')
    
    
    plt.figure()
    plt.plot( time, pl[0,:] , '.b')
    plt.plot( time, pr[0,:] , '.r')
    plt.plot( time, pr[0,:]+pl[0,:] , '.k')
    
    
    plt.figure()
    plt.plot( time, (1-pl[1,:])/pl[0,:] , '.b')
    plt.plot( time, (1-pr[1,:])/pr[0,:] , '.r')
    
    
    plt.figure()
    plt.plot(time, np.median(intens[:,85:95], axis=1))


if __name__=='__main__':
    #doit0()
    #doit1()
    #doit2()
    doit3()
    doit4()
    #doit5()
    plt.show()
