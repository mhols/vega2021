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

#DATAFILE = os.path.join(DATADIR, 'filematrix.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_nn.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_takeda.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_strong_s.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_strong_os.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_superstrong.dat')
#DATAFILE = os.path.join(DATADIR, 'filematrix_superstrong_notnorm.dat')
DATAFILE = os.path.join(DATADIR,'filematrix_453.dat')




nval      = 201                # number of data points
rangel    = np.arange(67,113)  # range of spectral line


coltime = 0 # colum of time values
colspec = 1
colval  = colspec + nval
colvul  = colval + nval

data = np.loadtxt(DATAFILE)
time = data[:,coltime:colspec].ravel()

velocity = data[0, colspec:colval]      # velocities of bins
intens   = data[:, colval:colvul]       # intensities
meani    = intens.mean(axis=0)          # mean intensity
diff     = intens-meani                 # fluctuation around mean
vari     = diff.ravel().var()           # variance of f
t0       = time[0]                      # first rime

# outlier removal
#I = np.arange(len(data[:,0]))
I = np.where( diff.var(axis=1) < 0.05*vari)[0]

print I
print len(I)
#subdata = data[I,rangel]
time = time[I]
velocity = velocity[rangel]
intens = intens[I,:]
intens = intens[:, rangel]

# recomputing mean and diff
meani = intens.mean(axis=0)
diff  = intens-meani

tt=time-time[0]

Inight5 = np.where( (tt-4)*(5-tt) > 0 )
time5 = time[Inight5]
inte5 = intens[Inight5,:]

Inight4 = np.where( (tt-3)*(4-tt) > 0 )
time4 = time[Inight4]
inte4 = intens[Inight4,:]

Inight3 = np.where( (tt-2)*(3-tt) > 0 )
time3 = time[Inight3]
inte3 = intens[Inight3,:]

Inight2 = np.where( (tt-1)*(2-tt) > 0 )
time2 = time[Inight2]
inte2 = intens[Inight2,:]

Inight1 = np.where( (tt-0)*(1-tt) > 0 )
time1 = time[Inight1]
inte1 = intens[Inight1,:]


list_time = [time1, time2, time3, time4, time5]
list_inte = [inte1, inte2, inte3, inte4, inte5]


night1_intens = 5

U, D, Vt = np.linalg.svd(intens, compute_uv = True, full_matrices=True)

print U.shape, D.shape, Vt.shape

def medianfilter(s,n):
    res = np.zeros(len(s)-n)
    for i in xrange(len(s)-n):
        res[i] = np.median(s[i:i+n])
    return res

def doit0():
    """
    The mean value of the spectral line
    """
    plt.figure()
    plt.title('mean value')
    plt.plot(velocity, meani)


def doit01():
    F=intens
    FtF = np.dot(F.T,F)
    v,e = np.linalg.eig(FtF)
    for i in range(4):
        plt.figure()
        plt.title("egeinf nr %d"%i)
        plt.plot(velocity, e[:,i])
        #plt.plot(velocity[rangel], Vt[i,:])

# subselection

#I = np.where( (velocity[0,:]+50)*(velocity[0,:]-20) < 0)

#print I
#velocity = velocity[:,I[0]]
#intens = intens[:,I[0]]



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
    plt.plot(D, 'o')
    plt.semilogy()
    
    plt.figure()
    plt.title('mean value')
    plt.plot(velocity, meani)
    
    for i in range(4):
        plt.figure()
        plt.title('EOF Nr. %d' %i)
        plt.plot(velocity, Vt[i,:])
        
        plt.figure()
        plt.title('evolution of component %d' %i)
        plt.plot(time, U[:,i])



def doit31(m):
    """
    median filter for time evolution
    """
    for i in range(4):
        
        tmp = medianfilter(U[:, i], m)
    
        plt.figure()
        plt.title('median filtered evolution of component %d' %i)
        plt.plot(time[:-m], tmp)


def doit4():
    maxamp = 30 / np.abs(intens-meani).std()
    for n in range(5): # five nights...
        plt.figure()
        plt.title("night #%d" %(n+1))
        for t, l in zip(time, (intens-meani)):
            plt.plot(velocity, maxamp*l+10000*(t-t0), '-')
            plt.ylim(10000*n,10000*n+3000) #TODO die Anfaenge der naechte besser verwalten...


def doit41(n,m):
    tmp = np.dot(np.dot(U[:,n:m],np.diag(D[n:m])),Vt[n:m,:])
    maxamp = 30 / np.abs(tmp).std()
    for n in range(5): # five nights...
        plt.figure()
        plt.title("night #%d" %(n+1))
        for t, l in zip(time, tmp):
            plt.plot(velocity, maxamp*l+10000*(t-t0), '-')
            plt.ylim(10000*n,10000*n+3000) #TODO die Anfaenge der naechte besser verwalten...
        plt.xlim( velocity[0], velocity[-1])

def doit5():
    
    ml = data[:, colval+72:colval+87]
    ll = data[:, colspec+72:colspec+87]
    
    mr = data[:, colval+94 :colval+109]
    lr = data[:, colspec+94:colspec+109]
    
    ml = ml / ml.sum(axis=1)[:,np.newaxis]
    mr = mr / mr.sum(axis=1)[:,np.newaxis]
    
    plt.figure()
    for l, r, u, v in zip(ml, mr, ll, lr):
        #plt.plot(u,l, '.')
        plt.plot(v, r, '.')
    
    #mll = (ll*ml).sum(axis=1)
    #mlr = (lr*mr).sum(axis=1)
    
    #plt.figure()
    #plt.plot(time, mll-mll.ravel().mean(), '.')
    #plt.plot(time, mlr-mlr.ravel().mean(), '.')
    
    
    
    
    #plt.figure()
    #plt.plot((ml*ll).sum(axis=1), (mr*lr).sum(axis=1), 'o')
    

    #plt.figure()
    #inte = intens[:,50:125]
    #inte /= inte.sum(axis=1)[:,np.newaxis]
    #velo = (velocity[:,50:125]*inte).sum(axis=1)
    #plt.plot(time, velo, '.')


    
    #plt.figure()
    #pl = np.polyfit(ll[0,:],ml.T,deg=1)
    #pr = np.polyfit(lr[0,:],mr.T,deg=1)
    
    
    #print pl.shape
    
    #plt.plot( pl[0,:], pr[0,:], 'o')
    
    
    #plt.figure()
    #plt.plot( time, pl[0,:] , '.b')
    #plt.plot( time, pr[0,:] , '.r')
    #plt.plot( time, pr[0,:]+pl[0,:] , '.k')
    
    
    #plt.figure()
    #plt.plot( time, (1-pl[1,:])/pl[0,:] , '.b')
    #plt.plot( time, (1-pr[1,:])/pr[0,:] , '.r')
    
    
    #plt.figure()
    #plt.plot(time, np.median(intens[:,85:95], axis=1))


from scipy.signal import spectral
vmin, vmax = -35, -25
I = np.where((velocity-vmin)*(vmax-velocity)>= 0)[0]
freqs = np.linspace(2*np.pi*0.1, 2*np.pi*50, 1024)

def doit6():
    
    plt.figure()
    for k in I:
        plt.figure()
        plt.plot(time, diff[:,k])
    
    res = np.zeros(len(freqs))
    plt.figure()
    
    #np.savetxt('sproutfort.dat', np.column_stack((time,diff[:,I[len(I)/2]])))
    wsp = spectral.lombscargle(time, np.ones(len(time)), freqs)
    wsp = abs(wsp)
    wsp /= wsp.max()
    for k in I:
        
        inten = diff[:,k]
        inten -= inten.mean()
        spk = spectral.lombscargle(time, inten, freqs)
        res += np.abs(spk)
        plt.plot(freqs/(2*np.pi), np.abs(spk))
        
        

    plt.figure()
    plt.plot(freqs/(2*np.pi), res)
    plt.plot(freqs/(2*np.pi), wsp*res.max(), '-k', linewidth=3)
    
def doit7():
    plt.figure()
    plt.plot(time, np.max(diff[:,I], axis=1))
    plt.plot(time, np.min(diff[:,I], axis=1))
    plt.plot(time, np.mean(diff[:,I], axis=1))
    
def doit8(n,m):
    """
    same as doit7 with eofs
    """
    dd = np.dot(np.dot(U[:,n:m],np.diag(D[n:m])),Vt[n:m,:])

    plt.figure()
    plt.plot(time, np.max(dd[:,I], axis=1))
    plt.plot(time, np.min(dd[:,I], axis=1))
    plt.plot(time, np.mean(dd[:,I], axis=1))
    
def doit9():
    plt.figure()
    plt.plot(time, intens.max(axis=1))
    plt.plot(time, intens.min(axis=1))
    plt.plot(time, intens.mean(axis=1))

def doit10():
    nt, nv = diff.shape
    nf = 1024
    freqs = np.linspace(2*np.pi*0.001, 2*np.pi*4, nf)
    res = np.empty( (nf, nv))
    
    dd = diff
    for i in xrange(nv):
        print i
        spk = spectral.lombscargle(time, dd[:,i], freqs)
        res[:,i] = np.abs(spk)
    
    plt.figure()
    res/=res.max()
    res = np.log(10*res + 1)
    plt.contourf(velocity, freqs/(2*np.pi), res)

def doit11():
    wsp = spectral.lombscargle(time, U[:,2], freqs)
    wsp = abs(wsp)
    plt.figure()
    plt.plot(freqs/(2*np.pi), wsp)

def doit12():
    plt.figure()
    print intens.shape
    plt.plot(time, intens[:,20])


def doit13():
    """
    mean profile for each night
    """
    for tt, ii in zip(list_time, list_inte ):
        d, nt, nv = ii.shape
        ii = ii.reshape((nt,nv))
        
        plt.figure()
        m = ii.mean(axis=0)
        
        for l in ii:
            plt.plot(velocity, l-m)
    

def doit14():
    """
    plotting the lomb scargle for nights individually
    """
    for tt, ii in zip(list_time, list_inte ):
        d, nt, nv = ii.shape
        ii = ii.reshape((nt,nv))
        nf = 1024
        freqs = np.linspace(2*np.pi*0.001, 2*np.pi*20, nf)
        res = np.empty( (nf, nv))
        m = ii.mean(axis=0)
        
        for i in xrange(nv):
            print i
            spk = spectral.lombscargle(tt, ii[:,i]-m[i], freqs)
            res[:,i] = np.abs(spk)
    
        plt.figure()
        res/=res.max()
        res = np.log(10*res + 1)
        plt.contourf(velocity, freqs/(2*np.pi), res)
    
    
if __name__=='__main__':
    doit0()
    #doit01()
    #doit1()
    #doit2()
    #doit3()
    #doit31(20)
    #doit4()
    #doit41(1,10)
    #doit5()
    #doit6()
    #doit7()
    #doit8(0,6)
    #doit9()
    #doit10() # Lomb skargel
    #doit11()
    #doit12()
    #doit13()
    #doit14()
    plt.savefig('spectogram.pdf')
    plt.show()
