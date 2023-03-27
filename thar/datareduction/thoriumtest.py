#compatible python 3, il a ete converti avec 2to3 et modification pickle
#ce code doit comparer deux spectres:
# 1) ThAr de reference UVES, calibre en longueur d'onde au mieux
#       voir http://astronomy.swin.edu.au/~mmurphy/thar/index.html
#  ici comparaison non pas au spectre même, mais à la liste de raies correspondante!
# 2) ThAr reduit par le DRS de Neo-Narval
#le but est de quantifier la precision du polynome de dispersion du DRS, ordre par ordre,
#region par region, le mieux serait raie par raie
#il faut lire les deux spectres, puis faire une fenetre glissante qu'on crosscorrelle avec la ref

import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import pickle
import numpy as np
import math
import string
import sys
import scipy
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from scipy import exp
from scipy.signal import savgol_filter
from scipy.signal import correlate
from scipy.optimize import curve_fit
from settings import *
from extract import *

#

myext = Extractor(**kwargs)

#ThAr moon
a=pyfits.open('../lune_res/HOBO_NEO_20200202_173811_th1.fits')
wave1=a[1].data['wavelength_1']
intens1=a[1].data['flux_1']
wave2=a[1].data['wavelength_2']
intens2=a[1].data['flux_2']
ordernumber=a[1].data['true_order_number']

#Atlas de reference
atlasline = []
with open(REF_ATLASLINES, 'r') as f:
    alines = f.readlines()
atlasline = np.array([float(l.split()[1]) for l in alines])

def gauss(x, A, mu, sigma, y_offset):
# A, sigma, mu, y_offset = p
   return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y_offset


    
def plot_xcorr(x, y):
    "Plot cross-correlation (full) between two signals."
    N = max(len(x), len(y))
    n = min(len(x), len(y))

    if N == len(y):
        lags = np.arange(-N + 1, n)
    else:
        lags = np.arange(-n + 1, N)
    c = correlate(x / np.std(x), y / np.std(y), 'full')

    #plt.plot(lags, c / n)
    
    return lags,c,n
    
    
montrer = 'montremoi1'
#montrer = 'non'
    


 
 
 
 
 

out1 = open("polycorr1.dat","w")
out2 = open("polycorr2.dat","w")
out12 = open("polycorr12.dat","w")
outmean = open("meansigma.dat","w")



# orderlim va de 0 a 152100, la 41eme valeur est a 0. Il y a donc 39 ordres en tout
# in range 38,39 veut dire que tout ce qui est en dessous de 39, a partir de 0, donc 38
# uniquement

meanall1 = zeros (orderlim.size, dtype =  float)
sigmaall1 = zeros (orderlim.size, dtype = float)
meanall1r = zeros (orderlim.size, dtype =  float)
sigmaall1r = zeros (orderlim.size, dtype = float)

meanall2 = zeros (orderlim.size, dtype =  float)
sigmaall2 = zeros (orderlim.size, dtype = float)
meanall2r = zeros (orderlim.size, dtype =  float)
sigmaall2r = zeros (orderlim.size, dtype = float)

meanall12 = zeros (orderlim.size, dtype =  float)
sigmaall12 = zeros (orderlim.size, dtype = float)
meanall12r = zeros (orderlim.size, dtype =  float)
sigmaall12r = zeros (orderlim.size, dtype = float)

meanall1ms = zeros (orderlim.size, dtype =  float)
sigmaall1ms = zeros (orderlim.size, dtype =  float)

meanall2ms = zeros (orderlim.size, dtype =  float)
sigmaall2ms = zeros (orderlim.size, dtype =  float)

meanall12ms = zeros (orderlim.size, dtype =  float)
sigmaall12ms = zeros (orderlim.size, dtype =  float)

meanresol1 = zeros (orderlim.size, dtype =  float)
sigmaresol1 = zeros (orderlim.size, dtype =  float)

meanresol2 = zeros (orderlim.size, dtype =  float)
sigmaresol2 = zeros (orderlim.size, dtype =  float)

medianresol1 = zeros (orderlim.size, dtype =  float)
medianresol2 = zeros (orderlim.size, dtype =  float)

ordernum = zeros (orderlim.size, dtype = float)
centrallam = np.zeros (orderlim.size, dtype = float)
tordernum = np.zeros (orderlim.size, dtype = float)

for j in range (2,37):

    ordernum[j] = j
# vrai numero d'ordre
    tordernum[j] = ordernum[j]+20
    print("j: ",j, " ordernum[j]", ordernum[j], " true ordernum[j]: ",tordernum[j])
#ici on ne selectionne que les raies de l'atlas dans l'ordre j, index.size est ce nombre de raies du catalogue UVES
    indexx=where((wave1[orderlim[j]]<atlasline) & (atlasline < wave1[orderlim[j+1]-1]))
    atlasext=atlasline[indexx]
    wext1=wave1[orderlim[j]:orderlim[j+1]-1]
    iext1=intens1[orderlim[j]:orderlim[j+1]-1]
    wext2=wave2[orderlim[j]:orderlim[j+1]-1]
    iext2=intens2[orderlim[j]:orderlim[j+1]-1]
    wext3=wave3[orderlim[j]:orderlim[j+1]-1]
    iext3=intens3[orderlim[j]:orderlim[j+1]-1]
    centrallam[j]=(wave1[orderlim[j+1]-1]+wave1[orderlim[j]])/2.
  

#
#  ICI on pourrait sortir les valeurs min et max de thararturo pour chaque ordre!
    
    
    
# seuil en ADU pour qu'on puisse selectionner une raie de l'atlas. Il faut qu'elle depasse de ce seuil le zero. Ensuite, il faut donner une zone +/- autour de la raie du catalogue, pour detecter la valeur maxi du flux tu ThAr obs.
# le vrange est donc cette largeur en km/s

    seuil = 2000.
    vrange = 9.
# on ne veut choisir que les raies de l'atlas qui se retrouvent dans la zone atlasext[k] +/- range, et qui ont un flux max au dessus du seuil. Ca reduit la liste. En plus, il faut le faire sur les deux voies et appliquer le mini.
    
    latlasext=atlasext*(1.-vrange/clum)
    ratlasext=atlasext*(1.+vrange/clum)
    
    numlines=int(np.array(indexx).size)
#    numlinesa=int(np.array(indexxa).size)
    maxi = zeros (numlines, dtype =  float)
    maxir = zeros (numlines, dtype =  float)
    intline1 = zeros (numlines, dtype =  float)
    intline2 = zeros (numlines, dtype =  float)
    
# selectionner que les raies du ThAr observes au dessus du seuil.
    for k in range (numlines):
        indext1 =  where((wext1<ratlasext[k]) & (wext1>latlasext[k]))
        w1=wext1[indext1]
        i1=iext1[indext1]
    
        indext2 =  where((wext2<ratlasext[k]) & (wext2>latlasext[k]))
        w2=wext1[indext2]
        i2=iext2[indext2]
        
        indextr =  where((refwave<ratlasext[k]) & (refwave>latlasext[k]))
        wr=refwave[indextr]
        ir=refintens[indextr]
  
        maxi[k]=min(max(i1),max(i2))
        maxir[k] = max(ir)
        
        
    indseuil= where((maxi > seuil) & (maxir > seuil))
    #number of selected reference lines
    nslines=int(np.array(indseuil).size)
    satlasext=atlasext[indseuil]
    slatlasext=latlasext[indseuil]
    sratlasext=ratlasext[indseuil]
    
# ici boucle sur toutes les raies selectionnees, correlation sur une petite zone autour de chaque raie entre l'intensite des deux voies individuellement et l'intensite du spectre de ref.

    corr1 = zeros (nslines, dtype =  float)
    corr2 = zeros (nslines, dtype =  float)
    corr12 = zeros (nslines, dtype =  float)
    corr1ms = zeros (nslines, dtype =  float)
    corr2ms = zeros (nslines, dtype =  float)
    corr12ms = zeros (nslines, dtype =  float)
    resol1 = zeros (nslines, dtype =  float)
    resol2 = zeros (nslines, dtype =  float)
    gaussposi1 = zeros (nslines, dtype =  float)
    gaussposi2 = zeros (nslines, dtype =  float)
    diff1 = zeros (nslines, dtype =  float)
    diff2 = zeros (nslines, dtype =  float)
    diff1ms= zeros (nslines, dtype =  float)
    diff2ms = zeros (nslines, dtype =  float)


    amp1 = zeros (nslines, dtype =  float)
    amp2 = zeros (nslines, dtype =  float)
    
    for k in range (nslines):
# attention, si on test pour seulement 2 raies, le programme se plante plus tard!
#    for k in range (5,6):


# je sais que la raie en question est la raie satlasext[k].
# si je trouve maintenant le pixel (sur 4208) avec sa fraction correspondant
# au centroide fitte, je peux faire directement le P(x) =  lambda sans
# passer par artlambda.dat et seconde correction
# aussi inclure le test sur les zones interdites (exclude.dat)








        indextr = where((refwave <sratlasext[k]) & (refwave>slatlasext[k]))
        wr=refwave[indextr]
        ir=refintens[indextr]
       
        indext1 =  where((wext1<max(wr)) & (wext1>min(wr)))
        w1=wext1[indext1]
        i1=iext1[indext1]
    
        indext2 =  where((wext2<max(wr)) & (wext2>min(wr)))
        w2=wext1[indext2]
        i2=iext2[indext2]
        

        
##        plt.plot(w1,i1,'r')
##        plt.plot(w2,i2,'y')
##        plt.plot(wr,ir,'b')
##        plt.show()
#
#        # on interpole ici le spectre de reference, et on l'evalue aux points de W1 et w2
#        cubic_interp= interp1d(wr, ir, kind='cubic')
#        cubic_results1 = cubic_interp(w1)
#        cubic_results2 = cubic_interp(w2)
#
        
##*************************** 1 versus 2 ***************************
#        #correlation entre raie voie 1 et voie 2 ThAr_obs
#        vals=plot_xcorr(i1,i2)
#        val=np.array(vals,dtype=object)
#        lags=val[0]
#        c=val[1]
#        n=val[2]
#
#        mu=0.
#        A=max(c/n)
#        sigma=0.2
#        y_offset=0.
#
#        try:
#            init_vals = [A, mu,sigma,y_offset]
#            gopt, covar = curve_fit(gauss,lags,c/n,p0=init_vals)
##        plt.plot(lags,gauss(lags,gopt[0],gopt[1],gopt[2],gopt[3]))
##        plt.vlines(gopt[1],0.,5.,'b')
##        plt.show()
#
#
#
#        #shift in km/s
##        corr12[k]=gopt[1]/satlasext[k]*clum*(max(w1)-min(w1))/int(w1.size)
#        #ici shift en Angstroms
#        #posi_courbe2 = posi_courbe1+shift
#        #raie_voie1 + corr12[k] = raie_voie2
#            corr12[k]=-gopt[1]*(max(w1)-min(w1))/int(w1.size)
#            corr12ms[k]=-gopt[1]*(max(w1)-min(w1))/(int(w1.size)-1)/satlasext[k]*clum*1000.
#            print(j,k,(max(w1)-min(w1))/(int(w1.size)-1)/satlasext[k]*clum*1000.)
##            print("shift A_12:", corr12[k])
#
#        #ici il faut correler ir,i1 et ir,i2; les dimensions doivent etre les memes, ainsi que les grilles de lambda
#        #il faut donc interpoler ir(wr) sur la grille de w1 et/ou w2 un sur l'autre par un cubic spline?
#        except (RuntimeError,TypeError):
#            print('Error on 12',k)
#
#
# #*************************** 1 versus ref ***************************
#
#
#        #correlation entre raie voie 1 et voie ref ThAr_obs
#        #puisqu'on ne correlle que les intensites, le resultat est en pas, cad taille du
#        #pas moyen en Angstrom.
#        vals=plot_xcorr(i1,cubic_results1)
#        val1=np.array(vals,dtype=object)
#        lags=val1[0]
#        c=val1[1]
#        n=val1[2]
#
#        mu=0.
#        A=max(c/n)
#        sigma=0.2
#        y_offset=0.
#        try:
#            init_vals = [A, mu,sigma,y_offset]
#            print('init_vals_1ref:', init_vals)
#            gopt, covar = curve_fit(gauss,lags,c/n,p0=init_vals)
#
#
#
##        [gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
#
#            if (montrer == 'montremoir'):
#                plt.figure()
#                plt.plot(lags,gauss(lags,gopt[0],gopt[1],gopt[2],gopt[3]))
#                plt.vlines(gopt[1],0.,5.,'b')
#                plt.show()
#                plt.close()
#
#
#        #shift in km/s
##        corr12[k]=gopt[1]/satlasext[k]*clum*(max(w1)-min(w1))/int(w1.size)
#        #ici shift en Angstroms
#        #raie_voie1 + corr1[k] = raie_ref
#            corr1[k]=-gopt[1]*(max(w1)-min(w1))/int(w1.size)
#            corr1ms[k]=-gopt[1]*(max(w1)-min(w1))/(int(w1.size)-1)/satlasext[k]*clum*1000.
##            print("shift A1:", corr1[k])
#        except (RuntimeError,TypeError):
#            print('Error on 1',k)
#
##############################resolution voie 1#############################################
            
        xxx=w1
        yyy=i1
        print("satlasext[k]", satlasext[k])
        print("xxx:", xxx)
        print("yyy:", yyy)
        print("mean(xxx):",mean(xxx))
        indrefextractmin=np.argmin(yyy)
        indrefextractmax=np.argmax(yyy)
        indrefmin=np.argmin(xxx)
        indrefmax=np.argmax(xxx)
        rA=(yyy[indrefextractmax]-yyy[indrefextractmin])
        rsigma=0.2*(xxx[indrefmax]-xxx[indrefmin])
#        rmu=xx[indrefextractmin]
        rmu=0.
        ry_offset=yyy[indrefextractmin]
#        print("ry_offset: ",y_offset)
#            gy=gauss(reflambda[indrefextract],A,mu,sigma,y_offset)
        rinit_vals = [rA, rmu,rsigma,ry_offset]  # for [amp, cen, wid]
        print ("rinit_vals fit EWref:", rinit_vals)
        errorbar=np.sqrt(yyy)
        try:
            rgopt, rcovar = curve_fit(gauss,xxx-mean(xxx),yyy,p0=rinit_vals,sigma=errorbar)
                
            
#            print("rmu:", gopt[1]," rA:",gopt[0]," rSigma:",gopt[2]," rY_offset:",gopt[3])
        #Return co-effs for fit and covariance
#            plt.figure()
#            plt.plot(xxx,yyy)
#            plt.plot(xxx,gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]),"r")
#            plt.show()
            FWHM =  rgopt[2] * 2 * sqrt(2*np.log(2))
            resol = mean(xxx)/FWHM
            resol1[k]=resol
            amp1[k]=rgopt[0]
            gaussposi1[k] = rgopt[1]+mean(xxx)
            print("gaussposi1[k]:", gaussposi1[k])
# difference satlasext - gaussposi
            diff1[k] = satlasext[k]-gaussposi1[k]
            diff1ms[k] = diff1[k]/satlasext[k]*clum*1000.
            print("Resolution:", resol)
            print("gauss:",gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]))
            print("gopt",rgopt)
        except (RuntimeError,TypeError):
            print("Ca n'a pas marche")
#[gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
# EW= A*(sqrt(2pi) sigma) si le continu est a 1!! Donc renormalisation AVANT
# la division par gopt[3] est une normalisation locale




############################################################################
#
#
##*************************** 2 versus ref ***************************
#
#        #correlation entre raie voie 2 et voie ref ThAr_obs
#        vals=plot_xcorr(i2,cubic_results1)
#        val2=np.array(vals,dtype=object)
#        lags=val2[0]
#        c=val2[1]
#        n=val2[2]
#
#        mu=0.
#        A=max(c/n)
#        sigma=0.2
#        y_offset=0.
#
#        try:
#            init_vals = [A, mu,sigma,y_offset]
#            gopt, covar = curve_fit(gauss,lags,c/n,p0=init_vals)
##        plt.plot(lags,gauss(lags,gopt[0],gopt[1],gopt[2],gopt[3]))
##        plt.vlines(gopt[1],0.,5.,'b')
##        plt.show()
#
#
#
#        #shift in km/s
##        corr12[k]=gopt[1]/satlasext[k]*clum*(max(w1)-min(w1))/int(w1.size)
#        #ici shift en Angstroms
#        #raie_voie2 + corr2[k] = raie_ref
#            corr2[k]=-gopt[1]*(max(w2)-min(w2))/int(w2.size)
#            corr2ms[k]=-gopt[1]*(max(w2)-min(w2))/(int(w2.size)-1)/satlasext[k]*clum*1000.
##            print("shift A2:", corr1[k])
#        except (RuntimeError,TypeError):
#            print('Error on 2',k)
#
##############################resolution voie 2#############################################
            
        xxx=w2
        yyy=i2
        indrefextractmin=np.argmin(yyy)
        indrefextractmax=np.argmax(yyy)
        indrefmin=np.argmin(xxx)
        indrefmax=np.argmax(xxx)
        rA=(yyy[indrefextractmax]-yyy[indrefextractmin])
        rsigma=0.2*(xxx[indrefmax]-xxx[indrefmin])
#        rmu=xx[indrefextractmin]
        rmu=0.
        ry_offset=yyy[indrefextractmin]
#        print("ry_offset: ",y_offset)
#            gy=gauss(reflambda[indrefextract],A,mu,sigma,y_offset)
        rinit_vals = [rA, rmu,rsigma,ry_offset]  # for [amp, cen, wid]
        print ("rinit_vals fit EWref:", rinit_vals)
        errorbar=np.sqrt(yyy)
        try:
            rgopt, rcovar = curve_fit(gauss,xxx-mean(xxx),yyy,p0=rinit_vals,sigma=errorbar)
                
            
#            print("rmu:", gopt[1]," rA:",gopt[0]," rSigma:",gopt[2]," rY_offset:",gopt[3])
        #Return co-effs for fit and covariance
#            plt.figure()
#            plt.plot(w1,i1)
#            plt.plot(xx,gauss(xx-mean(xx),gopt[0],gopt[1],gopt[2],gopt[3]))
#            plt.show()
            FWHM =  rgopt[2] * 2 * sqrt(2*np.log(2))
            resol = mean(xxx)/FWHM
            resol2[k]=resol
            amp2[k]=rgopt[0]
            gaussposi2[k] = rgopt[1]+mean(xxx)
            print("gaussposi1[k]:", gaussposi2[k])
# difference satlasext - gaussposi
            diff2[k] = satlasext[k]-gaussposi2[k]
            diff2ms[k] = diff2[k]/satlasext[k]*clum*1000.


#            print("Resolution 2:", resol)
        except (RuntimeError,TypeError):
            print("Ca n'a pas marche")
#[gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
# EW= A*(sqrt(2pi) sigma) si le continu est a 1!! Donc renormalisation AVANT
# la division par gopt[3] est une normalisation locale




###########################################################################
        
#*************************** sigma clipping ***************************

    
    #   sigmaclipping 3 fois
    satlasext1=satlasext
    satlasext2=satlasext
    
    satlasext1ms=satlasext
    satlasext2ms=satlasext

#*************************** 1 versus ref ***************************

    try:
        fitposi1=satlasext1-mean(satlasext1)
        fitposi1ms=fitposi1
        for ll in range (0,5):
            p1 = np.poly1d(np.polyfit(fitposi1, diff1, 5))
            sigma=np.std(p1(fitposi1)-diff1)
            index = where(abs(p1(fitposi1)-diff1) < 3*sigma)
            fitposi1=fitposi1[index]
            diff1=diff1[index]
            satlasext1=satlasext1[index]
        sigma1fit=sigma
        sigma1=np.std(diff1)
        mean1fit = np.mean(p1(fitposi1)-diff1)
        mean1=np.mean(diff1)
        sigmaall1[j]=sigma1
        meanall1[j]=mean1
    except (LinAlgError,TypeError):
        print('LinAlgError on ',j)
    
    try:
        for ll in range (0,10):
            p1ms = np.poly1d(np.polyfit(fitposi1ms, diff1ms, 5))
            sigma=np.std(p1ms(fitposi1ms)-diff1ms)
            indexms = where(abs(p1ms(fitposi1ms)-diff1ms) < 3*sigma)
            fitposi1ms=fitposi1ms[indexms]
            diff1ms=diff1ms[indexms]
            satlasext1ms=satlasext1ms[indexms]
        sigma1fitms = sigma
        mean1fitms = np.mean(p1ms(fitposi1ms)-diff1ms)
        sigma1ms = np.std(diff1ms)
        mean1ms = np.mean(diff1ms)
    except (LinAlgError,TypeError):
        print('LinAlgError on ',j)
    meanall1ms[j]=np.mean(diff1ms)
    sigmaall1ms[j]=np.std(diff1ms)
   

#    print('order',j)
#    print(diffms)
#    print(corr1ms)
#    print(corr2ms)
   
#    out.write("%15.7f  %15.7f  %15.7f \n" % (j,orderwave[j] ,p))
#    out.write("%15.7f  %15.7f  %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (j,orderwave[j],p[0],p[1],p[2],p[3],p[4],p[5]))
###    out1.write("%15.7f    %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (j,mean(satlasext1), p1[0],p1[1],p1[2],p1[3]))


#*************************** 2 versus ref ***************************

    try:
        fitposi2=satlasext2-mean(satlasext2)
        fitposi2ms=fitposi2
        for ll in range (0,10):
            p2 = np.poly1d(np.polyfit(fitposi2, diff2, 5))
            sigma=np.std(p2(fitposi2)-diff2)
            index = where(abs(p2(fitposi2)-diff2) < 3*sigma)
            fitposi2=fitposi2[index]
            diff2=diff2[index]
            satlasext2=satlasext2[index]
        sigma2fit=sigma
        sigma2=np.std(diff2)
        mean2fit = np.mean(p2(fitposi2)-diff2)
        mean2=np.mean(diff2)
        sigmaall2[j]=sigma2
        meanall2[j]=mean2
    except (LinAlgError,TypeError):
        print('LinAlgError on ',j)
    
    try:
        for ll in range (0,5):
            p2ms = np.poly1d(np.polyfit(fitposi2ms, diff2ms, 5))
            sigma=np.std(p2ms(fitposi2ms)-diff2ms)
            indexms = where(abs(p2ms(fitposi2ms)-diff2ms) < 3*sigma)
            fitposi2ms=fitposi2ms[indexms]
            diff2ms=diff2ms[indexms]
            satlasext2ms=satlasext2ms[indexms]
        sigma2fitms = sigma
        mean2fitms = np.mean(p2ms(fitposi2ms)-diff2ms)
        sigma2ms = np.std(diff2ms)
        mean2ms = np.mean(diff2ms)
    except (LinAlgError,TypeError):
        print('LinAlgError on ',j)
    meanall2ms[j]=np.mean(diff2ms)
    sigmaall2ms[j]=np.std(diff2ms)
   

#    print('order',j)
#    print(diffms)
#    print(corr1ms)
#    print(corr2ms)
   
#    out.write("%15.7f  %15.7f  %15.7f \n" % (j,orderwave[j] ,p))
#    out.write("%15.7f  %15.7f  %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (j,orderwave[j],p[0],p[1],p[2],p[3],p[4],p[5]))
###    out2.write("%15.7f    %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (j,mean(satlasext2), p2[0],p2[1],p2[2],p2[3]))

#
#*************************** 1 versus 2 ***************************
#
#
#    try:
#        fitposi12=satlasext12-mean(satlasext12)
#        fitposi12ms=fitposi12
#        for ll in range (0,5):
#            p12 = np.poly1d(np.polyfit(fitposi12, corr12, 3))
#            sigma=np.std(p12(fitposi12)-corr12)
#            index = where(abs(p12(fitposi12)-corr12) < 3*sigma)
#            fitposi12=fitposi12[index]
#            corr12=corr12[index]
#            satlasext12=satlasext12[index]
#            slatlasext12=slatlasext12[index]
#            sratlasext12=sratlasext12[index]
#        sigma12 = sigma
#        mean12 = np.mean(p12(fitposi12)-corr12)
#        sigma12r = np.std(corr12)
#        mean12r = np.mean(corr12)
#        sigmaall12[j]=sigma12
#        meanall12[j]=mean12
#        sigmaall12r[j]=sigma12r
#        meanall12r[j]=mean12r
#    except (LinAlgError,TypeError):
#        print('LinAlgError on ',j)
#
#
#    try:
#        for ll in range (0,5):
#            p12ms = np.poly1d(np.polyfit(fitposi12ms, corr12ms, 3))
#            sigma=np.std(p12ms(fitposi12ms)-corr12ms)
#            indexms = where(abs(p12ms(fitposi12ms)-corr12ms) < 3*sigma)
#            fitposi12ms=fitposi12ms[indexms]
#            corr12ms=corr12ms[indexms]
#            satlasext12ms=satlasext12ms[indexms]
#        sigma12ms = sigma
#        mean12ms = np.mean(p12ms(fitposi12ms)-corr12ms)
#        sigma12rms = np.std(corr12ms)
#        mean12rms = np.mean(corr12ms)
#    except (LinAlgError,TypeError):
#        print('LinAlgError on ',j)
#    meanall12ms[j]=np.mean(corr12ms)
#    sigmaall12ms[j]=np.std(corr12ms)
#
#
#
    try:
        fitresol1=satlasext-mean(satlasext)
        satlasextresol1=satlasext
        sigma1r = 0.
        sigma2r = 0.
        for ll in range (0,5):
            presol1 = np.poly1d(np.polyfit(fitresol1, resol1, 3))
            sigma1r=np.std(presol1(fitresol1)-resol1)
            indexresol1 = where(abs(presol1(fitresol1)-resol1) < 2.5*sigma1r)
            fitresol1=fitresol1[indexresol1]
            resol1=resol1[indexresol1]
            amp1=amp1[indexresol1]
            gaussposi1=gaussposi1[indexresol1]
            
            satlasextresol1=satlasextresol1[indexresol1]


    except (LinAlgError,TypeError):
        print('LinAlgError on ',j)
    sigmaresol1[j] = sigma1r
    meanresol1[j] = np.mean(resol1)
    medianresol1[j] = np.median(resol1)

    try:
        fitresol2=satlasext-mean(satlasext)
        satlasextresol2=satlasext
        for ll in range (0,5):
            presol2 = np.poly1d(np.polyfit(fitresol2, resol2, 3))
            sigma2r=np.std(presol2(fitresol2)-resol2)
            indexresol2 = where(abs(presol2(fitresol2)-resol2) < 2.5*sigma2r)
            fitresol2=fitresol2[indexresol2]
            resol2=resol2[indexresol2]
            amp2=amp2[indexresol2]
            gaussposi2=gaussposi2[indexresol2]
            satlasextresol2=satlasextresol2[indexresol2]

    except (LinAlgError,TypeError):
        print('LinAlgError on ',j)
    sigmaresol2[j] = sigma2r
    meanresol2[j] = np.mean(resol2)
    medianresol2[j] = np.median(resol2)
#
##    out.write("%15.7f  %15.7f  %15.7f \n" % (j,orderwave[j] ,p))
##    out.write("%15.7f  %15.7f  %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (j,orderwave[j],p[0],p[1],p[2],p[3],p[4],p[5]))
#    out12.write("%15.7f    %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (j,mean(satlasext12), p12[0],p12[1],p12[2],p12[3]))
#
#    outmean.write("%15.7f    %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (j,mean(satlasext12), mean1rms,sigma1rms,mean1ms,sigma1ms))

    plt.figure(figsize=(16,6))
#    print(j,wave1[orderlim[j]],wave1[orderlim[j+1]-1])
    plt.xlim(wext1[0],wext1[-1])
    plt.plot(wext1,iext1,'r')
    tit= "order:" + str(j)
    plt.title(tit)
    plt.plot(refwave,refintens*6.,'b')
 #   plt.plot(wave2D[0:304200],intens1,"y")
 #   plt.plot(wave2[orderlim[j]:orderlim[j+1]-1],intens2[orderlim[j]:orderlim[j+1]-1],'g')
    
    for ll in range (numlines):
        plt.vlines(atlasext[ll],30000.,50000.)
    for ll in range (satlasext1.size):
        plt.vlines(satlasext1[ll],20000.,50000.,'y')
#        plt.vlines(slatlasext1[ll],0.,3000.,'b')
#        plt.vlines(sratlasext1[ll],0.,3000.,'r')
#    for kk in range (numlinesa):
#        plt.vlines(thararturoext[kk],10000.,20000.,'b')
    if (montrer == 'montremoi1'):
        plt.show()
    plt.close()

    
#    plt.figure(figsize=(16,6))
#    plt.plot(satlasext12,corr12,'o')
#    plt.plot(wext1,p12(wext1-mean(satlasext12)))
#    # shift pour chaque raie, decalage du shift_k, moyenne des sigma par rapport a la
#    # courbe ajustee, moyenne de cette distance
#    tit= " 1 + shift_k = 2, order:" + str(j)+" mean:" +str(mean12r) +" sigma:" +str(sigma12r)
#    nomsave = "shift_12_order_" +str(j)
#    plt.title(tit)
#    xlabel('Angstrom')
#    ylabel('Angstrom')
#    plt.savefig(nomsave,dpi=300)
#    if (montrer == 'montremoiA'):
#        plt.show()
#    plt.close()
#
    plt.figure(figsize=(16,6))
    plt.plot(amp1,resol1,'bo')
    plt.plot(amp2,resol2,'ro')
    #
    tit= " resolution (amplitude), order:" + str(j)
    nomsave = "./results_hobo/resolution_ampli_" +str(j)
    plt.title(tit)
    xlabel('amplitude')
    ylabel('Resolution')
    plt.savefig(nomsave,dpi=300)
    if (montrer == 'montremoi1'):
        plt.show()
    plt.close()
    
    
    plt.figure(figsize=(16,6))
    plt.plot(gaussposi1,resol1,'bo')
    plt.plot(gaussposi2,resol2,'ro')
    #
    tit= " resolution (amplitude), order:" + str(j)
    nomsave = "./results_hobo/resolution_pix_" +str(j)
    plt.title(tit)
    xlabel('lambda')
    ylabel('Resolution')
    plt.savefig(nomsave,dpi=300)
    if (montrer == 'montremoi1'):
        plt.show()
    plt.close()
#    plt.figure(figsize=(16,6))
#    plt.plot(wext1,p12(wext1-mean(satlasext12)))
#    # shift pour chaque raie, decalage du shift_k, moyenne des sigma par rapport a la
#    # courbe ajustee, moyenne de cette distance
#    tit= " 1 + shift_k = 2, order:" + str(j)+" mean:" +str(mean12r) +" sigma:" +str(sigma12r)
#    nomsave = "shift_12_order_" +str(j)
#    plt.title(tit)
#    xlabel('Angstrom')
#    ylabel('Angstrom')
#    plt.savefig(nomsave,dpi=300)
#    if (montrer == 'montremoiA'):
#        plt.show()
#    plt.close()

    
    plt.figure(figsize=(16,6))
    plt.plot(satlasext1,diff1,'o')
    plt.plot(wext1,p1(wext1-mean(satlasext1)))
    #plt.plot(satlasext1,corr1-p1(satlasext1-mean(satlasext1)),'o')
    tit= " diff(atlas - posi1), order:" + str(j) +" mean:" +str(mean1)+" sigma:" +str(sigma1)
    nomsave = "./results_hobo/shift_1ref_order_" +str(j)
    plt.title(tit)
    xlabel('Angstrom')
    ylabel('Angstrom')
    plt.savefig(nomsave,dpi=300)
    if (montrer == 'montremoi1'):
        plt.show()
    plt.close()

    
    plt.figure(figsize=(16,6))
    plt.plot(satlasext2,diff2,'o')
    plt.plot(wext2,p2(wext2-mean(satlasext2)))
    tit= "diff(atlas - posi2), order:" + str(j) +" mean:" +str(mean2) +" sigma:" +str(sigma2)
    nomsave = "./results_hobo/shift_2ref_order_" +str(j)
    plt.title(tit)
    xlabel('Angstrom')
    ylabel('Angstrom')
    plt.savefig(nomsave,dpi=300)
    if (montrer == 'montremoi1'):
        plt.show()
    plt.close()
    
#    plt.figure(figsize=(16,6))
#    plt.plot(satlasext12ms,diff12ms,'o')
#    plt.plot(wext1,p12ms(wext1-mean(satlasext12ms)))
#    # shift pour chaque raie, decalage du shift_k, moyenne des sigma par rapport a la
#    # courbe ajustee, moyenne de cette distance
#    tit= " 1 + shift_k = 2, order:" + str(j) +" mean:" +str(mean12ms)+" sigma:" +str(sigma12ms)
#    nomsave = "ms_shift_12_order_" +str(j)
#    plt.title(tit)
#    xlabel('Angstrom')
#    ylabel('m/s')
#    plt.savefig(nomsave,dpi=300)
#    if (montrer == 'montremoi2'):
#        plt.show()
#    plt.close()
#
#
    plt.figure(figsize=(16,6))
    plt.plot(satlasext1ms,diff1ms,'o')
    plt.plot(wext1,p1ms(wext1-mean(satlasext1ms)))
    tit= "diff(atlas - posi1), order:" + str(j) +" mean:" +str(mean1ms)+" sigma:" +str(sigma1ms)+"/"+str(sigma1fitms)
    nomsave = "./results_hobo/ms_shift_1ref_order_" +str(j)
    plt.title(tit)
    xlabel('Angstrom')
    ylabel('m/s')
    plt.savefig(nomsave,dpi=300)
    if (montrer == 'montremoi1'):
        plt.show()
    plt.close()
#    plt.show()
    
    plt.figure(figsize=(16,6))
    plt.plot(satlasext2ms,diff2ms,'o')
    plt.plot(wext2,p2ms(wext2-mean(satlasext2ms)))
    tit= "diff(atlas - posi2), order:" + str(j) +" mean:" +str(mean2ms)+" sigma:" +str(sigma2ms)+"/"+str(sigma2fitms)
    nomsave = "./results_hobo/ms_shift_2ref_order_" +str(j)
    plt.title(tit)
    xlabel('Angstrom')
    ylabel('m/s')
    plt.savefig(nomsave,dpi=300)
    if (montrer == 'montremoi1'):
        plt.show()
    plt.close()
    
    plt.figure(figsize=(16,6))
    plt.plot(satlasextresol1,resol1,'o')
    tit= "resolution voie 1, order:" + str(j) +" mean:" +str(meanresol1[j])+" sigma:" +str(sigmaresol1[j])
    nomsave = "./results_hobo/resol1_order_" +str(j)
    plt.title(tit)
    xlabel('Angstrom')
    ylabel('Resolution')
    plt.savefig(nomsave,dpi=300)
    if (montrer == 'montremoi1'):
        plt.show()
    plt.close()
    
    plt.figure(figsize=(16,6))
    plt.plot(satlasextresol2,resol2,'o')
    tit= "resolution voie 2, order:" + str(j) +" mean:" +str(meanresol1[j])+" sigma:" +str(sigmaresol1[j])
    nomsave = "./results_hobo/resol2_order_" +str(j)
    plt.title(tit)
    xlabel('Angstrom')
    ylabel('Resolution')
    plt.savefig(nomsave,dpi=300)
    if (montrer == 'montremoi1'):
        plt.show()
    plt.close()
    
#plt.show()

#plt.figure(figsize=(8,10))
#plt.subplot(3,2,1)
#tit= "mean shift 1 versus 2 before correction (o) and after (x)"
#plt.title(tit,fontsize=5)
#ylabel('Angstrom')
#plt.plot(ordernum,meanall12r,'o')
#plt.plot(ordernum,meanall12,'x')
#
#plt.subplot(3,2,2)
#tit= "standard deviation shift 1 versus 2 before correction (o) and after (x)"
#plt.title(tit,fontsize=5)
#plt.plot(ordernum,sigmaall12r,'o')
#plt.plot(ordernum,sigmaall12,'x')
#
#plt.subplot(3,2,3)
#tit= "mean shift 1 versus ref before correction (o) and after (x)"
#plt.title(tit, fontsize=5)
#ylabel('Angstrom')
#plt.plot(ordernum,meanall1r,'o')
#plt.plot(ordernum,meanall1,'x')
#
#plt.subplot(3,2,4)
#tit= "standard deviation shift 1 versus ref before correction (o) and after (x)"
#plt.title(tit,fontsize=5)
#plt.plot(ordernum,sigmaall1r,'o')
#plt.plot(ordernum,sigmaall1,'x')
#
#plt.subplot(3,2,5)
#tit= "mean shift 2 versus ref before correction (o) and after (x)"
#plt.title(tit,fontsize=5)
#xlabel('ordernum')
#ylabel('Angstrom')
#plt.plot(ordernum,meanall2r,'o')
#plt.plot(ordernum,meanall2,'x')
#
#plt.subplot(3,2,6)
#tit= "standard deviation shift 2 versus ref before correction (o) and after (x)"
#plt.title(tit,fontsize=5)
#xlabel('ordernum')
#ylabel('Angstrom')
#plt.plot(ordernum,sigmaall2r,'o')
#plt.plot(ordernum,sigmaall2,'x')
#nomsave = "DRS_calibration_polynome.jpg"
#plt.savefig(nomsave,dpi=300)
#plt.show()


## plot ici en m/s
#
#plt.figure(figsize=(8,10))
#plt.subplot(3,2,1)
#tit= "mean shift 1 versus 2 before correction (o)"
#plt.title(tit,fontsize=5)
#ylabel('m/s')
#plt.plot(ordernum,meanall12ms,'o')
#
#
#plt.subplot(3,2,2)
#tit= "standard deviation shift 1 versus 2"
#plt.title(tit,fontsize=5)
#plt.plot(ordernum,sigmaall12ms,'o')
#
#
#plt.subplot(3,2,3)
#tit= "mean shift 1 versus ref"
#plt.title(tit, fontsize=5)
#ylabel('m/s')
#plt.plot(ordernum,meanall1ms,'o')
#
#plt.subplot(3,2,4)
#tit= "standard deviation shift 1 versus ref"
#plt.title(tit,fontsize=5)
#plt.plot(ordernum,sigmaall1ms,'o')
#
#
#plt.subplot(3,2,5)
#tit= "mean shift 2 versus ref"
#plt.title(tit,fontsize=5)
#xlabel('ordernum')
#ylabel('m/s')
#plt.plot(ordernum,meanall2ms,'o')
#
#plt.subplot(3,2,6)
#tit= "standard deviation shift 2 versus ref"
#plt.title(tit,fontsize=5)
#xlabel('ordernum')
#plt.plot(ordernum,sigmaall2ms,'o')
#
#nomsave = "DRS_calibration_polynome_ms_without_correction.jpg"
#plt.savefig(nomsave,dpi=300)
#plt.show()

ymin = 30000.
ymax = 80000.
xmin = 1.
xmax = 39.


plt.figure(figsize=(16,6))
tit="Resolution voie 1(b) et voie 2(r) mediane: ("+ str(int(median(meanresol1))) + "/" +str(int(median(meanresol2))) +")"
xlabel('ordernum')
ylabel('Resolution')
plt.title(tit)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.plot(ordernum,meanresol1,'bo')
plt.plot(ordernum,medianresol1,'bx')
plt.plot(ordernum,meanresol2,'ro')
plt.plot(ordernum,medianresol2,'rx')
nomsave = "./results_hobo/DRS_resolution.jpg"
plt.savefig(nomsave,dpi=300)
plt.show()
plt.close()

out1.close()
out2.close()
out12.close()
outmean.close()


