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
from pylab import *
import matplotlib.pyplot as plt
import math
import string
import sys
import scipy
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from scipy.special import erf
from numpy import exp
from scipy.signal import savgol_filter
from scipy.signal import correlate
from scipy.optimize import curve_fit
#
#
f=open('ThAr_Line_Cal_v2.pkl','rb')
b=pickle.load(f,encoding='bytes')
# liste de raies utilises par Arturo
thararturo = np.array(b[b'Lines Atlas'])
clum=3e5

#lecture du fichier ascii du spectre ThAr de reference (UVES)
#ce spectre contient une premiere colonne lambda, une seconde intensite,...
#213358 lignes

refname = 'thar_spec_MM201006.dat'
#thar
#213358 4
#3033.0867480159  103.99224526863826    0.36633407377209    0.36633407377209     97.132933537024

#nombre de lignes du catalogue de raies
num_lines=len(open(refname,'r').readlines())
npt=num_lines-2
linecount = 0
refwave = zeros (npt, dtype =  float)
refintens = zeros (npt, dtype =  float)


fr = open(refname, 'r')
linefr = fr.readline()
linefr = fr.readline()
for l in range (0,npt):
        linefr = fr.readline()
        sfr = linefr.split()
        refwave[l] = float(sfr[0])
        refintens[l] = float(sfr[1])
fr.close


#lecture du fichier ascii catalogue de raies thar_UVES_MM090311.dat
#33322.4046  3000.109256  0.450 Ar II  W
#33319.9704  3000.328442 -0.007 XX 0   P
atlasname = 'thar_UVES_MM090311.dat'
num_alines=len(open(atlasname,'r').readlines())
npt=num_alines
atlasline = zeros (npt, dtype =  float)
ar = open(atlasname, 'r')
for l in range (0,npt):
        linear = ar.readline()
        sfr = linear.split()
        atlasline[l] = float(sfr[1])
fr.close



#lecture du fichier ThAr reduit avec le DRS de Neo-Narval
#a=pyfits.open('NEO_20200202_173811_th2.fits')
#a=pyfits.open('NEO_20200202_173811_thB.fits')

#urversion
#a=pyfits.open('NEO_20200901_184654_th2.fits')
#####a=pyfits.open('NEO_20200901_184654_thX.fits')
a=pyfits.open('NEO_20220219_173048_th2.fits')
wave1=a[1].data['Wavelength1']
intens1=a[1].data['Beam1']
wave2=a[1].data['Wavelength2']
intens2=a[1].data['Beam2']

orderlim=a[2].data['Orderlimit']

def gauss(x, A, mu, sigma, y_offset):
# A, sigma, mu, y_offset = p
   return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y_offset

"""
def _gauss(x, A, mu, sigma, y_offset):
    return y_offset + 0.5 * np.sqrt(np.pi) * (erf((x+0.5-mu)/(np.sqrt(2)*sigma)) - erf((x-0.5-mu)/(np.sqrt(2)*sigma)))
"""

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
    
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


outtableau1 = open("ThAr2D_voie1.dat","w")
outtableau2 = open("ThAr2D_voie2.dat","w")

#montrer = 'montremoi1'

montrer = 'non'

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
centrallam = zeros (orderlim.size, dtype = float)
tordernum = zeros (orderlim.size, dtype = float)

#for j in range (10,11):
#for j in range (31,orderlim.size-2):
for j in range (2,orderlim.size-2):

    ordernum[j] = j
# vrai numero d'ordre
    tordernum[j] = ordernum[j]+21
    print("j: ",j, " ordernum[j]", ordernum[j], " true ordernum[j]: ",tordernum[j])
#ici on ne selectionne que les raies de l'atlas dans l'ordre j, index.size est ce nombre de raies du catalogue UVES
    indexx=where((wave1[orderlim[j]]<atlasline) & (atlasline < wave1[orderlim[j+1]-1]))
    atlasext=atlasline[indexx]
    wext1=wave1[orderlim[j]:orderlim[j+1]-1]
    iext1=intens1[orderlim[j]:orderlim[j+1]-1]
    wext2=wave2[orderlim[j]:orderlim[j+1]-1]
    iext2=intens2[orderlim[j]:orderlim[j+1]-1]
    indexxa=where((wave1[orderlim[j]]<thararturo) & (thararturo < wave1[orderlim[j+1]-1]))
    thararturoext=thararturo[indexxa]
    centrallam[j]=(wave1[orderlim[j+1]-1]+wave1[orderlim[j]])/2.
#
#  ICI on pourrait sortir les valeurs min et max de thararturo pour chaque ordre!
#
#orderlim[j+1]-orderlim[j] vaut 7800, il y a donc 7800 pixels (=Arturos) par ordre
    
    
# seuil en ADU pour qu'on puisse selectionner une raie de l'atlas. Il faut qu'elle depasse de ce seuil le zero. Ensuite, il faut donner une zone +/- autour de la raie du catalogue, pour detecter la valeur maxi du flux tu ThAr obs.
# le vrange est donc cette largeur en km/s
    seuil = 100.
    vrange = 9.
# on ne veut choisir que les raies de l'atlas qui se retrouvent dans la zone atlasext[k] +/- vrange, et qui ont un flux max au dessus du seuil. Ca reduit la liste. En plus, il faut le faire sur les deux voies et appliquer le mini.
    
    latlasext=atlasext*(1.-vrange/clum)
    ratlasext=atlasext*(1.+vrange/clum)
    
    numlines=int(np.array(indexx).size)
    numlinesa=int(np.array(indexxa).size)
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
    pixelposi1 = zeros (nslines, dtype =  float)
    pixelposi2 = zeros (nslines, dtype =  float)
    erreurposi1 = zeros (nslines, dtype =  float)
    erreurposi2 = zeros (nslines, dtype =  float)

    diff1 = zeros (nslines, dtype =  float)
    diff2 = zeros (nslines, dtype =  float)
    diff1ms= zeros (nslines, dtype =  float)
    diff2ms = zeros (nslines, dtype =  float)


    amp1 = zeros (nslines, dtype =  float)
    amp2 = zeros (nslines, dtype =  float)
    
    for k in range (nslines):
# attention, si on teste pour seulement 2 raies, le programme se plante plus tard!
#    for k in range (5,6):

# on incremente le compteur general de raies spectrales de 1
# ici on pretend que voies 1 et 2 sont suffisamment proches pour que le nombre
# de raies est identique
        print("linecount: ",linecount)
        linecount += 1
        indextr = where((refwave <sratlasext[k]) & (refwave>slatlasext[k]))
        wr=refwave[indextr]
        ir=refintens[indextr]
       
        indext1 =  where((wext1<max(wr)) & (wext1>min(wr)))
        w1=wext1[indext1]
        i1=iext1[indext1]
        pix1=np.array(indext1)
        pixel1=pix1[0]
    
        indext2 =  where((wext2<max(wr)) & (wext2>min(wr)))
        w2=wext1[indext2]
        i2=iext2[indext2]
        pix2=np.array(indext2)
        pixel2=pix2[0]


##############################resolution voie 1#############################################
            
        xxx=w1
        yyy=i1
        print("satlasext[k]", satlasext[k])
#        print("xxx:", xxx)
#        print("yyy:", yyy)
#        print("pixel1:",pixel1)
#        print("mean(xxx):",mean(xxx))
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
        rsigmap=0.2*(np.argmax(pixel1)-np.argmin(pixel1))
        rinit_valsp = [rA, rmu, rsigmap, ry_offset]
        print ("rinit_vals fit EWref:", rinit_vals)
        try:
            rgopt, rcovar = curve_fit(gauss,xxx-mean(xxx),yyy,p0=rinit_vals)
            rgoptp, rcovarp = curve_fit(gauss,pixel1-mean(pixel1),yyy,p0=rinit_valsp)
# la matrice de covariance rcovarp (en pixels) contient sur sa diagonale les erreurs associes
# aux differents parametres. L'erreur sur la position du fit, donc relatif a rgoptp[1]
# s'exprime alors
            errs=np.diag(rcovarp)
            erreurposi1[k] = errs[1]
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
            pixelposi1[k] = rgoptp[1]+mean(pixel1)
            print("gaussposi1[k]:", gaussposi1[k])
            print("pixelposi1[k]:", pixelposi1[k])

# ici on ecrit le tableau de sortie (i, lam_i, true_ordre_i, posi_i,err_posi_i)
            outtableau1.write("%4d %20.10e %2d %20.10e %20.10e" % (linecount, gaussposi1[k], tordernum[j], pixelposi1[k], erreurposi1[k])+str(pixel1)+str(yyy)+"\n")
            outtableau1.write(pixel1,yyy)
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


##############################resolution voie 2#############################################
            
        xxx=w2
        yyy=i2
        print("satlasext[k]", satlasext[k])
#        print("xxx:", xxx)
#        print("yyy:", yyy)
#        print("pixel2:",pixel2)
#        print("mean(xxx):",mean(xxx))
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
        rsigmap=0.2*(np.argmax(pixel2)-np.argmin(pixel2))
        rinit_valsp = [rA, rmu, rsigmap, ry_offset]
        print ("rinit_vals fit EWref:", rinit_vals)
        try:
            rgopt, rcovar = curve_fit(gauss,xxx-mean(xxx),yyy,p0=rinit_vals)
            rgoptp, rcovarp = curve_fit(gauss,pixel2-mean(pixel2),yyy,p0=rinit_valsp)
            errs=np.diag(rcovarp)
            erreurposi2[k] = errs[1]

            
#            print("rmu:", gopt[1]," rA:",gopt[0]," rSigma:",gopt[2]," rY_offset:",gopt[3])
        #Return co-effs for fit and covariance
#            plt.figure()
#            plt.plot(xxx,yyy)
#            plt.plot(xxx,gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]),"r")
#            plt.show()
            FWHM =  rgopt[2] * 2 * sqrt(2*np.log(2))
            resol = mean(xxx)/FWHM
            resol2[k]=resol
            amp2[k]=rgopt[0]
            gaussposi2[k] = rgopt[1]+mean(xxx)
            pixelposi2[k] = rgoptp[1]+mean(pixel2)
            print("k: ",k)
            print("gaussposi1[k]:", gaussposi2[k])
            print("pixelposi1[k]:", pixelposi2[k])
            # ici on ecrit le tableau de sortie (i, lam_i, true_ordre_i, posi_i,err_posi_i)
            outtableau2.write("%15.7f %20.10e %20.10e %20.10e %20.10e\n" % (linecount, gaussposi2[k], tordernum[j], pixelposi2[k], erreurposi2[k])+str(pixel2)+str(yyy)+"\n")
# difference satlasext - gaussposi
            diff2[k] = satlasext[k]-gaussposi2[k]
            diff2ms[k] = diff2[k]/satlasext[k]*clum*1000.
#            print("Resolution:", resol)
#            print("gauss:",gauss(xxx-mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]))
#            print("gopt",rgopt)
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
        for ll in range (0,10):
            p1 = np.poly1d(np.polyfit(fitposi1, diff1, 3))
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
        for ll in range (0,5):
            p1ms = np.poly1d(np.polyfit(fitposi1ms, diff1ms, 3))
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

#    out1.write("%15.7f    %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (j,mean(satlasext1), p1[0],p1[1],p1[2],p1[3]))


#*************************** 2 versus ref ***************************

    try:
        fitposi2=satlasext2-mean(satlasext2)
        fitposi2ms=fitposi2
        for ll in range (0,10):
            p2 = np.poly1d(np.polyfit(fitposi2, diff2, 3))
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
            p2ms = np.poly1d(np.polyfit(fitposi2ms, diff2ms, 3))
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
   
#    out2.write("%15.7f    %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (j,mean(satlasext2), p2[0],p2[1],p2[2],p2[3]))


    try:
        fitresol1=satlasext-mean(satlasext)
        satlasextresol1=satlasext
        for ll in range (0,5):
            presol1 = np.poly1d(np.polyfit(fitresol1, resol1, 3))
            sigma1r=np.std(presol1(fitresol1)-resol1)
            indexresol1 = where(abs(presol1(fitresol1)-resol1) < 2.5*sigma1r)
            fitresol1=fitresol1[indexresol1]
            resol1=resol1[indexresol1]
            amp1=amp1[indexresol1]
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
    print(j,wave1[orderlim[j]],wave1[orderlim[j+1]-1])
    plt.xlim(wave1[orderlim[j]],wave1[orderlim[j+1]-1])
    plt.plot(wave1[orderlim[j]:orderlim[j+1]-1],intens1[orderlim[j]:orderlim[j+1]-1],'r')
    tit= "order:" + str(j)
    plt.title(tit)
    plt.plot(refwave,refintens,'b')
    plt.plot(wave2[orderlim[j]:orderlim[j+1]-1],intens2[orderlim[j]:orderlim[j+1]-1],'g')
    
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

"""
index=where(tordernum>0.)
lam=centrallam[index]
ord=tordernum[index]

plt.figure()
plt.plot(lam*ord,'b')
plt.plot(lam*(ord+1))
plt.plot(lam*(ord-1))
plt.plot(lam*(ord+2))
plt.plot(lam*(ord-2))

plt.title("true ordernum identification")
plt.show()

plt.figure()
plt.plot(ord,lam*ord,'b')
plt.plot(ord,lam*(ord+1))
plt.plot(ord,lam*(ord-1))
plt.plot(ord,lam*(ord+2))
plt.plot(ord,lam*(ord-2))

plt.title("true ordernum identification")
plt.show()
"""

outtableau1.close()
outtableau2.close()

