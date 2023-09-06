#compatible python 3, il a ete converti avec 2to3 et modification pickle
#ici on compare les resultats du DRS avec celui de LE, en gardant les limites d'ordre bcp plus petits de LE,
#sinon on penalise la qualite du DRS.

import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import pickle
import numpy
from pylab import *
import matplotlib.pyplot as plt
import math
import string
import sys
import scipy
from scipy.interpolate import UnivariateSpline
from scipy import exp
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
#



resolution=30000.
clum=3e8
eadu=2.0
totalsnoises=[]
totalEWs=[]
totalEWrefs=[]
totalline=[]
totalposi=[]
#on decale le spectre de la lune de delta_v (ici 1089m/s)
#Vesta:
wave2Dname = 'NEO_20200202_173811_th2.wave2D.voie_1.dat'
wave2D = np.loadtxt(wave2Dname)


#delta_v=(8085.18-8084.50)/8085.*clum
#delta_v=0.
#     350.016      26.0000     0.607000      2.87600      2.25000       1
fname = 'sun_mask.5750.45.p00.correct'
#nombre de lignes du catalogue de raies
num_lines=len(open(fname,'r').readlines())
npt=num_lines

wave = zeros (npt, dtype =  float)
prof = zeros (npt, dtype =  float)
fr = open(fname, 'r')
for l in range (0,npt):
        linefr = fr.readline()
        sfr = linefr.split()
        wave[l] = float(sfr[0])*10.
        prof[l] = float(sfr[2])
fr.close

#lecture du fichier de raies telluriques

telfile = '/Users/boehm/Desktop/vega2021/thar/reffiles/telluric_lines/telluric_mask_atlas_short.dat'

TEL = np.loadtxt(telfile)
nlines=TEL.size
telwvl = zeros (nlines, dtype =  float)
telflag = zeros (nlines, dtype =  float)

telwvl = TEL[:,0]
telflag = TEL[:,1]
    


oname = 'orderlim_LE.dat'
orderlines=len(open(oname,'r').readlines())


onum = zeros(orderlines, dtype = float)
obeg = zeros(orderlines, dtype = float)
oend = zeros(orderlines, dtype = float)
osize = zeros(orderlines, dtype = float)
orderbeg = zeros(orderlines, dtype = float)
orderend = zeros(orderlines, dtype = float)

f = open(oname, 'r')
for j in range (0, orderlines):
    line = f.readline()
    s = line.split()
    onum[j] = float(s[0])
    obeg[j] = float(s[1])
    oend[j] = float(s[2])
    osize[j] = float(s[3])
f.close()



#LE n'a pas de raies de ref dans son ordre le plus bleu, celui ci s'arretant a 3820A, or
#les raies les plus bleues identifiees de orderlim.size-2 sont pour DRS_NEO dans la partie rouge
#de cet ordre au dela de 3820A. On decide alors pour la comparaison de limiter à orderlim.size-3
#DRS = 'LE'
DRS = 'DRS_Neo'
montrer = 'montremoi'
#montrer = 'non'

if (DRS == 'LE'):
    delta_v=(7090.39-7090.91)/7090.*clum
    lefile='12jan12_norm_moon_031.s'
    l=open(lefile,'r')
    line = l.readline()
    line = l.readline()
    s = line.split()
    nlines=int(s[0])

    wvl1 = zeros (nlines, dtype =  float)
    int1 = zeros (nlines, dtype =  float)
    noise1 = zeros (nlines, dtype =  float)

    for j in range (0, nlines):
        line = l.readline()
        s = line.split()
        wvl1[j] = 10*float(s[0])
        int1[j] = float(s[1])
        noise1[j] = float(s[2])
    l.close
    wvl1=wvl1*(1+delta_v/clum)


    orderl=[wvl1[0]]

    for j in range (1,nlines-1):
        if ((wvl1[j]<wvl1[j-1]) or (wvl1[j] > wvl1[j-1]+0.5)):
            orderl.append(wvl1[j-1])

    orderl.append(wvl1[nlines-1])
    orderli=np.array(orderl)
    osize=orderli.size

    orderlim = zeros(osize, dtype = int)
    sorderlim = orderlim

    for k in range (0,osize):
        index = where(wvl1 == orderli[k])
#    print(index)
        orderlim[k]=index[0][0]

else:


# delta_v=(6609.13-6609.26)/6609.*clum
    delta_v=0.
#    a=pyfits.open('NEO_20200202_182606_stX.fits')
    nnfile = 'moon_int_Neo-Narval.s'
    
    l=open(nnfile,'r')
    line = l.readline()
    line = l.readline()
    s = line.split()
    nlines=int(s[0])

    wvl1 = zeros (nlines, dtype =  float)
    int1 = zeros (nlines, dtype =  float)
    noise1 = zeros (nlines, dtype =  float)

    for j in range (0, nlines):
        line = l.readline()
        s = line.split()
        wvl1[j] = 10*float(s[0])
        int1[j] = float(s[1])
        noise1[j] = float(s[2])
    l.close
    wvl1=wvl1*(1+delta_v/clum)


    orderl=[wvl1[0]]
    
    
    
    

    for j in range (1,nlines-1):
        if ((wvl1[j]<wvl1[j-1]) or (wvl1[j] > wvl1[j-1]+0.5)):
            orderl.append(wvl1[j-1])

    orderl.append(wvl1[nlines-1])
    orderli=np.array(orderl)
    osize=orderli.size

    orderlim = zeros(osize, dtype = int)
    sorderlim = orderlim

    for k in range (0,osize):
        index = where(wvl1 == orderli[k])
#    print(index)
        orderlim[k]=index[0][0]
         

f=open('SolarSpectrum.pkl3','rb')
b=pickle.load(f)
#f=open('SolarAtlas_lines.pkl','rb')
f=open('Selectedlines_fit.pkl3','rb')
c=pickle.load(f)

tamp=[]
tamp=list(set(c['sellines']))
lines=sort(np.array(tamp))
print("number of selected lines: ", lines.size)
#lines = numpy.array(c['sellines'])
llim = zeros(lines.size, dtype =  float)
rlim = zeros(lines.size, dtype =  float)
equw = zeros(lines.size, dtype =  float)


reflambda = numpy.array(b['Wavelength'])
refintens = numpy.array(b['Intensity'])
isoline = numpy.zeros(lines.size)
wavel = zeros (lines.size, dtype =  float)
waver = zeros (lines.size, dtype =  float)

def gauss(x, A, mu, sigma, y_offset):
# A, sigma, mu, y_offset = p
   return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y_offset

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

for k in range (1,lines.size-1):
    equw[k]=1.
for k in range (1,lines.size-1):
#    print (rlim[k-1],llim[k],rlim[k],llim[k+1])
######   deltalam=1.5*lines[k]/resolution
    deltalam=2.*lines[k]/resolution
#    if ((rlim[k-1] < lines[k]-deltalam) & (lines[k]+deltalam < llim[k+1])):
#        isoline[k]=1.
    isoline[k]=1.
isol=where(isoline == 1)
isolines=lines[isol]
print(isolines)
#leftline=llim[isol]
#rightline=rlim[isol]
leftline=isolines*(1-1./resolution)
rightline=isolines*(1+1./resolution)
Rleftline=isolines*(1-1./resolution)
Rrightline=isolines*(1+1./resolution)
equws=equw[isol]

#ici definition des mesures d'erreur concernant EW_mes versus EW_th
#distance  = 1/2 * sqrt((EW - EWref)**2+(EWref - EW)**2)
#          = 1/2 * sqrt(2*(EW-EWref)**2)
#signe: Si EW > EWref alors +
#       Si EW < EWref alors -
#       Si EW = EWref alors diagonale!
taille=orderlim.size-2
mean_EW_dist = zeros (taille, dtype =  float)
median_EW_dist = zeros (taille, dtype =  float)
sigma_EW_dist = zeros (taille, dtype =  float)
mean_EW_dist_A = zeros (taille, dtype =  float)
median_EW_dist_A = zeros (taille, dtype =  float)
sigma_EW_dist_A = zeros (taille, dtype =  float)
mean_EW_dist_n = zeros (taille, dtype =  float)
median_EW_dist_n = zeros (taille, dtype =  float)
sigma_EW_dist_n = zeros (taille, dtype =  float)
mean_EW_dist_A_n = zeros (taille, dtype =  float)
median_EW_dist_A_n = zeros (taille, dtype =  float)
sigma_EW_dist_A_n = zeros (taille, dtype =  float)
sigmaposi = zeros (taille, dtype =  float)
ordernumber = zeros (taille, dtype =  float)
orderwave = zeros (taille, dtype =  float)

# orderlim va de 0 a 152100, la 41eme valeur est a 0. Il y a donc 39 ordres en tout
# in range 38,39 veut dire que tout ce qui est en dessous de 39, a partir de 0, donc 38
# uniquement

ormin=2
ormax=orderlim.size-20

#ormin=15
#ormax=16

#ormin=orderlim.size-12
###ormax=orderlim.size-8

sormin=ormin
sormax=ormax
#ormin=29
#ormax=33



print("ormin: ",ormin," ormax: ",ormax)
print("trueormin: ",ormin+21, " trueormax: ", ormax+21)
#ormax=orderlim.size-2
#ormin=2,ormax=orderlim.size-2)
"""
if (DRS == 'LE'):
    print("Libre Esprit")
    ormin=orderlim.size-sormax
    ormax=orderlim.size-sormin-1
    print("ormin: ",ormin," ormax: ",ormax)
    

if (DRS == 'LE'):
    if (ormax == 38):
        ormax=ormax-1
"""
ormin=orderlim.size-sormax
ormax=orderlim.size-sormin-1
print("ormin: ",ormin," ormax: ",ormax)

if (ormax == 38):
    ormax=ormax-1
    
sorderlim = orderlim
for j in range (ormin,ormax):

 #   plt.figure
    plt.figure(figsize=(16,6))
   
    
    
    
    print(j,wvl1[orderlim[j]],wvl1[orderlim[j+1]-1])
    plt.xlim(wvl1[orderlim[j]],wvl1[orderlim[j+1]-1])
    plt.plot(wvl1[orderlim[j]:orderlim[j+1]-1],int1[orderlim[j]:orderlim[j+1]-1])
    plt.vlines(obeg[j],0.1, 1.0, linewidth=10, alpha=0.3)
    plt.vlines(oend[j],0.1, 1.0, linewidth=10, alpha=0.3)
    plt.plot(b['Wavelength'],b['Intensity'],'r')

    indextel=where((wvl1[orderlim[j]]<telwvl) & (telwvl < wvl1[orderlim[j+1]-1]))
    telw = telwvl[indextel]
    telf = telflag[indextel]
    
    for ll in range (0,telw.size-1):
        if telf[ll] == 0:
            plt.vlines(telw[ll],0.1, 0.2, "r")
        else:
            plt.vlines(telw[ll],0.1, 0.6, "b")
        
    
#    plt.plot(reflambda,smooth(refintens,20))
##    ysavgol = savgol_filter(int1[orderlim[j]:orderlim[j+1]-1], 5, 2) # window size 5, polynomial order 3
##    plt.plot(wvl1[orderlim[j]:orderlim[j+1]-1],ysavgol,label='Savgol Filter')
#    spl=UnivariateSpline(wvl1[orderlim[j]:orderlim[j+1]-1],int1[orderlim[j]:orderlim[j+1]-1], np.gradient(int1[orderlim[j]:orderlim[j+1]-1]), k=5)
#   spl.set_smoothing_factor(0.0001)
#   plt.plot(wvl1[orderlim[j]:orderlim[j+1]-1],spl(wvl1[orderlim[j]:orderlim[j+1]-1]),label='Smooth Fct 3')
#    spl.set_smoothing_factor(10)
#    plt.plot(wvl1[orderlim[j]:orderlim[j+1]-1],spl(wvl1[orderlim[j]:orderlim[j+1]-1]),label='Smooth Fct 10')
    plt.legend(loc='lower left')

#    max_idx = np.argmax(spl(np.arange(int1[orderlim[j]:orderlim[j+1]-1])))
#    plt.vlines(max_idx, -5, 9, linewidth=5, alpha=0.3)

    
    
    
    
#ici on ne selectionne que les raies isolees dans l'ordre j, index.size est ce nombre de raies du catalogue .pkl

#orderbeg et end definissent le domaine spectral en longueur d'onde commun entre DRS et LE
#pour l'ordre j donne (attention, il y a inversion d'ordres entre DRS et LE.

    orderbeg[j] = max(obeg[j],wvl1[orderlim[j]+1])
    orderend[j] = min(oend[j],wvl1[orderlim[j+1]-1])
#    print(j, 'obeg ',obeg[j],' oend ',oend[j], wvl1[orderlim[ormax-j]+1], wvl1[orderlim[ormax-j+1]-1])
    
# ATTENTION ICI PROBLEME D'INVERSION DU SENS DES ORDRES...PB A REGLER POUR LE OBEG ET OEND
    
    
    
#ici il faut trouver l'indice correspondant à la longueur d'onde la plus proche dans le DRS ou LE
    
    extsize=wvl1[orderlim[j]+1:orderlim[j+1]-1].size
    barg = zeros (extsize, dtype =  float)
    earg = zeros (extsize, dtype =  float)

# barg determine la longueur d'onde la plus proche de orderbeg, idem pour earg
    for k in range (orderlim[j]+1,orderlim[j+1]-1):
        barg[k-(orderlim[j]+1)] = abs(wvl1[k]-orderbeg[j])
        earg[k-(orderlim[j]+1)] = abs(wvl1[k]-orderend[j])
        
    obarg = np.argmin(barg)
    oearg = np.argmin(earg)
 
    obarg = obarg + orderlim[j]+1
    oearg = oearg + orderlim[j]+1
    
    """
    if (DRS == 'LE'):
        obarg = orderlim[j]
        oearg = orderlim[j+1]
    """
    obarg = orderlim[j]
    oearg = orderlim[j+1]
 
 
 
 
        
#
#   index=where((wvl1[orderlim[j]+1]<isolines) & (isolines < wvl1[orderlim[j+1]-1]))
    index=where((wvl1[obarg]<isolines) & (isolines < wvl1[oearg]))
#en plus selection sur flux dans continu
    clevel=0.1
    
#    wext=wvl1[orderlim[j]+1:orderlim[j+1]-1]
#    wint=int1[orderlim[j]+1:orderlim[j+1]-1]
#    wnoise=noise1[orderlim[j]+1:orderlim[j+1]-1]

#    wcont=cont1[orderlim[j]:orderlim[j+1]-1]
    
# attention, wvl1 est decroissant, pour cela orderend:orderbeg
    wext=wvl1[obarg:oearg]
    wint=int1[obarg:oearg]
    wnoise=noise1[obarg:oearg]

#    wcont=cont1[oearg:obarg]
#
#    werr=error1[orderlim[j]:orderlim[j+1]-1]
    
#    wcontmax=max(wcont)
#    cindex=where(
    
    orderline=isolines[index]
    leftlines=leftline[index]
    rightlines=rightline[index]
    Rleftlines=Rleftline[index]
    Rrightlines=Rrightline[index]
    eq=equws[index]
    init_vals=[]
    
    R = zeros (orderline.size, dtype =  float)
    EWref = zeros (orderline.size, dtype =  float)
    refposi = zeros (orderline.size, dtype =  float)
    posi = zeros (orderline.size, dtype =  float)
    deltaposi = zeros (orderline.size, dtype =  float)
    EW = zeros (orderline.size, dtype =  float)
    snoise = zeros (orderline.size, dtype = float)
    EWint = zeros (orderline.size, dtype =  float)
    erreur = zeros (orderline.size, dtype =  float)
    EWreferr = zeros (orderline.size, dtype =  float)
    EWerr = zeros (orderline.size, dtype =  float)

#    EW_dist_lin = zeros (orderline.size, dtype =  float)
    print(orderline.size)
    for k in range (orderline.size):
#        print("line number:", k, orderline[k])
        plt.axvline(x=orderline[k], ymin=0., ymax=1000.,c='y')
        plt.axvline(x=leftlines[k],ymin=0.7,ymax=2,c='b')
        plt.axvline(x=rightlines[k],ymin=0.7,ymax=2,c='r')
        plt.vlines(Rleftlines[k],0.9, 1.0, linewidth=5, alpha=0.3)
        plt.vlines(Rrightlines[k],0.9, 1.0, linewidth=5, alpha=0.3)
        
# calcul d'un fit gaussien pour chaque raie
# -> attention, eventuellement normalisation locale necessaire entre left- et rightlines
        
#travail sur le spectre de reference
        try:
            indrefextract=where((reflambda >= leftlines[k]) & (reflambda <= rightlines[k]))
            xx=reflambda[indrefextract]
            yy=refintens[indrefextract]
            indrefextractmin=np.argmin(yy)
            indrefextractmax=np.argmax(yy)
            indrefmin=np.argmin(xx)
            indrefmax=np.argmax(xx)
            A=-(yy[indrefextractmax]-yy[indrefextractmin])
            sigma=0.2*(xx[indrefmax]-xx[indrefmin])
            mu=xx[indrefextractmin]
            y_offset=yy[indrefextractmax]
#            gy=gauss(reflambda[indrefextract],A,mu,sigma,y_offset)
            init_vals = [A, mu,sigma,y_offset]  # for [amp, cen, wid]
#            print ("init_vals fit EWref:", init_vals)
            try:
                gopt, covar = curve_fit(gauss,xx,yy,p0=init_vals)
#        print(k,gopt[1],gopt[0],gopt[2],gopt[3])
        #Return co-effs for fit and covariance
                plot(xx,gauss(xx,gopt[0],gopt[1],gopt[2],gopt[3]))
        #[gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
        # EW= A*(sqrt(2pi) sigma) si le continu est a 1!! Donc renormalisation AVANT

        # la division par gopt[3] est une normalisation locale
#                print("para fit EWref: ",gopt)
                EWref[k] = -gopt[0]*sqrt(2*pi)*gopt[2]/gopt[3]
                errref = np.sqrt(np.diag(covar))
                perr = np.sqrt(np.diag(covar))
                refposi[k]=gopt[1]
#                print('perr ref: ', perr)
                EWreferr[k] = abs(perr[0]*sqrt(2*pi)*perr[2])
            except:
                EWref[k] = 1000.
                refposi[k] = 1000.
                pass
        except (ValueError):
            print('value error ici')
            pass
            
            
#        print("on est ici")
        windextract=where((wext>=Rleftlines[k]) & (wext<=Rrightlines[k]))
        
        
###### ici on a un pb pour windextract!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
#        yphat = savgol_filter(py, 101, 3) # window size 51, polynomial order 3
#        wx=wvl1[windextract]
#        wy=int1[windextract]
#        iq = np.argsort(wx)
#        x=wx[iq]
#        y=wy[iq]
         
         
#  ici on a selectionne dans le spectre observe une zone correspondant a la raie, mais a +/- lam/resolution
        wx=wext[windextract]
        wy=wint[windextract]
        wn=wnoise[windextract]
        windextractmin=np.argmin(wy)
        windextractmax=np.argmax(wy)
        windmin=np.argmin(wx)
        windmax=np.argmax(wx)
        A=-(wy[windextractmax]-wy[windextractmin])
#        sigma=0.2*(y[windmax]-y[windmin])
        mu=wx[windextractmin]
        sigma=mu/resolution
        y_offset=wy[windextractmax]
        #        gy=gauss(reflambda[indrefextract],A,mu,sigma,y_offset)
        xxx=linspace(wx[windmin],wx[windmax],100)
        init_vals = [A, mu,sigma,y_offset]  # for [amp, cen, wid]
#        print ("init_vals fit EW:", init_vals)
        try:
            gopt, covar = curve_fit(gauss,wx,wy,p0=init_vals,maxfev=1000)
            errs = np.sqrt(np.diag(covar))
#Return co-effs for fit and covariance
#            plot(wx,gauss(wx,gopt[0],gopt[1],gopt[2],gopt[3]),'y')
            plot(xxx,gauss(xxx,gopt[0],gopt[1],gopt[2],gopt[3]),'y')
#            print(k,gopt[1],gopt[0],gopt[2],gopt[3])
#            print(wx)
#            print(wy)
#            print("para fit EW: ",gopt)
            #gopt[2]=sigma est souvent negatif ici pour EW...va comprendre...
#            EW[k] = -gopt[0]*sqrt(2*pi)*gopt[2]/gopt[3]
            EW[k] = -gopt[0]*sqrt(2*pi)*abs(gopt[2])/gopt[3]
            posi[k] = gopt[1]
            snoise[k] = sqrt(np.max(wn)*eadu)
#           gy=gauss(xxx,gopt[0],gopt[1],gopt[2],gopt[3])
#            gindmin=np.argmin(gy)
            R[k]=abs(gopt[1]/(2*gopt[2]*sqrt(2*log(2))))
#            print(k,EW[k],EWref[k],gopt[1],gopt[0],gopt[2],gopt[3])
            #gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset
            gy=gauss(wx,gopt[0],gopt[1],gopt[2],gopt[3])
            erreur[k]=sqrt(sum((wy-gy)**2))/wx.size
            print(gopt[1], erreur[k])
            perr = np.sqrt(np.diag(covar))
#            print('perr obs: ',perr)
            EWerr[k] = abs(perr[0]*sqrt(2*pi)*perr[2])
        except (RuntimeError,TypeError,ValueError):
#            print('Error on',k)
            EW[k] = 1000.
            posi[k]=0.
            R[k]=0.
            erreur[k]=1000.
            pass
            
        def normalisation(x,y):
            indmin=np.argmin(x)
            indmax=np.argmax(x)
            pente=(y[indmax]-y[indmin])/(x[indmax]-x[indmin])
            droite=y[indmin]+pente*(x-x[indmin])
            return (droite)

        droite = normalisation(wx,wy)
 #       plot(wx,droite)
            
 #   indexw=where((wvl1[orderlim[j]]<wave) & (wave < wvl1[orderlim[j+1]-1]))
    indexw=where((wvl1[obarg]<wave) & (wave < wvl1[oearg]))
    profex=prof[indexw]
    waveex=wave[indexw]
     
     
    for l in range (0,waveex.size):
        plt.vlines(waveex[l],1.-profex[l], 1.0, linewidth=5, alpha=0.3)
    tit= DRS + " spectre - order " + str(j)
    plt.title(tit)
    if (montrer == 'montremoi'):
        plt.show()
    
    indice = (erreur != 1000.)
    sig  = std(erreur[indice])
    ind = (erreur < 3 *sig)
    
#   $$$$$$$$$$$$$$$$$$$$$$$$$ on a calcule un refposi, posi et on a prépare un deltaposi [k]
#   maintenant il faut voir par rapport a numok ci dessous
#   le but est de pouvoir plotter deltaposi en fonction de refposi, calculer le sigma de la distribution de deltaposi
#   (et le mean) et puis de l integrer dans un rapport final pour detecter les erreurs sur la calib. Un moment donne il
#   faudra aussi regarder les deux premiers ordres...

    EWrefmax=max(EWref[ind])
    eqmax=max(eq[ind])
#    maxi=max(eqmax,EWrefmax)
    maxi=EWrefmax
    xmin=0.
    xmax=maxi
#    xmax=1.
    ymin=0.
    ymax=maxi
#    ymin=-10.
#    ymax=10.
    xxp=linspace(0,maxi,100)
    yyp=xxp
    
    inds = ((abs(EW) < 2*maxi) & (abs(EW) < 1.) & (abs(EWref) < 1.))
    
    # ici calcul des distances de chaque EW a la diagonale
    numok = erreur[inds].size
    EW_dist_lin = zeros (numok, dtype =  float)
    EW_dist_lin_A = zeros (numok, dtype =  float)
    EW_dist_lin_n = zeros (numok, dtype =  float)
    EW_dist_lin_A_n = zeros (numok, dtype =  float)
    eqs = eq[inds]
    EWrefs = EWref[inds]
    EWs = EW[inds]
    snoises = snoise[inds]
    
    posi=posi[inds]
    refposi=refposi[inds]

    deltaposi=((posi-refposi)/posi)*(clum/1000.)

    sigmaposi[j]=np.std(deltaposi)
    orderwave[j]=(max(posi)+min(posi))/2.
    print("j",j,sigmaposi[j])
    
    for ll in range (0,numok-1):
#        print EWs[ll],EWrefs[ll], EWs[ll]-EWrefs[ll]
#        if (abs(EWs[ll]) >= EWrefs[ll]):
#            EW_dist_lin[ll] = 0.5 * sqrt(2*(abs(EWs[ll])-EWrefs[ll])**2)
#            EW_dist_lin_n[ll] = 0.5 * sqrt(2*(abs(EWs[ll])-EWrefs[ll])**2)/abs(EWrefs[ll])
#        else:
#            EW_dist_lin[ll] = - 0.5 * sqrt(2*(abs(EWs[ll])-EWrefs[ll])**2)
#            EW_dist_lin_n[ll] = - 0.5 * sqrt(2*(abs(EWs[ll])-EWrefs[ll])**2)/abs(EWrefs[ll])
        EW_dist_lin[ll] =  0.5 * sqrt(2*(abs(EWs[ll])-EWrefs[ll])**2)
        EW_dist_lin_n[ll] =  0.5 * sqrt(2*(abs(EWs[ll])-EWrefs[ll])**2)/abs(EWrefs[ll])

        if (abs(EWs[ll]) >= eqs[ll]):
            EW_dist_lin_A[ll] = 0.5 * sqrt(2*(abs(EWs[ll])-eqs[ll])**2)
            EW_dist_lin_A_n[ll] = 0.5 * sqrt(2*(abs(EWs[ll])-eqs[ll])**2)/abs(EWrefs[ll])
        else:
            EW_dist_lin_A[ll] = - 0.5 * sqrt(2*(abs(EWs[ll])-eqs[ll])**2)
            EW_dist_lin_A_n[ll] = - 0.5 * sqrt(2*(abs(EWs[ll])-eqs[ll])**2)/abs(EWrefs[ll])
   
    mean_EW_dist[j] = mean(EW_dist_lin)
    median_EW_dist[j] = median(EW_dist_lin)
    sigma_EW_dist[j]= std(EW_dist_lin)
    mean_EW_dist_A[j] = mean(EW_dist_lin_A)
    median_EW_dist_A[j] = median(EW_dist_lin_A)
    sigma_EW_dist_A[j]= std(EW_dist_lin_A)
    
    mean_EW_dist_n[j] = mean(EW_dist_lin_n)
    median_EW_dist_n[j] = median(EW_dist_lin_n)
    sigma_EW_dist_n[j]= std(EW_dist_lin_n)
    mean_EW_dist_A_n[j] = mean(EW_dist_lin_A_n)
    median_EW_dist_A_n[j] = median(EW_dist_lin_A_n)
    sigma_EW_dist_A_n[j]= std(EW_dist_lin_A_n)
    
    
    plt.figure(figsize=(8,10))
    #    plt.figure
    plt.subplot(3,2,1)
    tit= DRS + " Posi diff obs-ref (ref) - order "+ str(j)
    plt.title(tit)
    ylabel('km/s')
    ymin=-2.
    ymax=2.
    plt.ylim(ymin,ymax)
    plot(refposi,deltaposi,'o', c='r')
    
#    plt.subplot(3,2,2)
#    tit= "EW_dist_lin (EWrefs)- order "+ str(j)
#    plt.title(tit)
#    plot(EWrefs,EW_dist_lin,'o')
    
    plt.subplot(3,2,2)
    tit= "Intensity - order "+ str(j)
    plt.title(tit)
    plot(wext,1./(wnoise**2))

    plt.subplot(3,2,3)
    tit= "EW_dist_lin_n (lambda)- order "+ str(j)
    plt.title(tit)
    orl = orderline[inds]
    plot(orl,EW_dist_lin_n,'o')
    
    plt.subplot(3,2,4)
    tit= "EW_dist_lin_n (EWref)- order "+ str(j)
    plt.title(tit)
    plot(EWrefs,EW_dist_lin_n,'o')

    plt.subplot(3,2,5)
    tit= "EW_obs / EWref (S/N) - order "+ str(j)
    plt.title(tit)
    plot(1./snoises, EWs/EWrefs,'o')

    xmin=0.
    xmax=0.20
    ymin=0.
    ymax=0.20
#    plt.subplot(3,2,(5,6))
    plt.subplot(3,2,6)
    tit= "EW_obs(EW_ref_Tor) - order " + str(j)
    plt.title(tit)
    xlabel('EW_ref')
    ylabel('EW')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.plot(EWref[ind],abs(EW[ind]),'.')
    plt.vlines(EWref[ind],abs(EW[ind])-EWerr[ind],abs(EW[ind])+EWerr[ind], linewidth=2, alpha=0.2)
    plot(xxp,yyp,'-')

    nomsave= DRS + "_summary_order_" + str(j) +".jpg"
    plt.savefig(nomsave,dpi=300)
    if (montrer == 'montremoi'):
        plt.show()
    plt.close()
    
    totalsnoises=totalsnoises+list(snoises)
    totalEWs=totalEWs+list(EWs)
    totalEWrefs=totalEWrefs+list(EWrefs)
    totalline=totalline+list(orl)
        
    
           
            
        #travail sur le spectre observe
#        indextract=where((wvl1 >= leftlines[k]) & (wvl1 <= rightlines[k]))
#        indextractmin=np.argmin(int1[indextract])
        
        
#        A=-(1.-fyhat[findex_min])
#        sigma=0.2*(x[indmax]-x[indmin])
#        mu=fx[findex_min]
#       y_offset=yhat[indmin]

#      gy=gauss(fx,A,mu,sigma,y_offset)


#     init_vals = [A, mu,sigma,y_offset]  # for [amp, cen, wid]
#     gopt, covar = curve_fit(gauss,fx,fy,p0=init_vals) #Return co-effs for fit and covariance

montrer = 'montremoi'

indi = np.arange(mean_EW_dist.size)

plt.close()


plt.figure(figsize=(8,10))
plt.subplot(3,2,1)
tit= DRS + " mean_EW_dist"
plt.title(tit)
plot(indi,mean_EW_dist,'o')
    
plt.subplot(3,2,2)
tit= "median_EW_dist"
plt.title(tit)
plot(indi,median_EW_dist,'o')
    
plt.subplot(3,2,3)
tit= "sigma_EW_dist"
plt.title(tit)
plot(indi,sigma_EW_dist,'o')

plt.subplot(3,2,4)
tit= "mean_EW_dist_n"
plt.title(tit)
plot(indi,mean_EW_dist_n,'o')

plt.subplot(3,2,5)
tit= "median_EW_dist_n"
plt.title(tit)
plot(indi,median_EW_dist_n,'o')
plt.vlines(indi,median_EW_dist_n+sigma_EW_dist,median_EW_dist_n-sigma_EW_dist, linewidth=2, alpha=0.2)

plt.subplot(3,2,6)
tit= "sigma_EW_dist_n"
plt.title(tit)
plot(indi,sigma_EW_dist_n,'o')

nomsave= DRS +"_overall summary.jpg"
plt.savefig(nomsave,dpi=300)
if (montrer == 'montremoi'):
    plt.show()
plt.close()



totsn=np.array(totalsnoises)
totEWs=np.array(totalEWs)
totEWrefs=np.array(totalEWrefs)
totline=np.array(totalline)

plt.figure(figsize=(16,6))
tit= DRS + " EWs/EWrefs (S/N)"
plt.title(tit)
plot(1./totsn, totEWs/totEWrefs,'o')
plt.axhline(1.,c="r")
nomsave= DRS + "_overall_EWratio_SN.jpg"
plt.savefig(nomsave,dpi=300)
if (montrer == 'montremoi'):
    plt.show()
plt.close()

plt.figure(figsize=(16,6))
valeur=np.std(totEWs/totEWrefs)
mvaleur=median(totEWs/totEWrefs)
tit= DRS + " EWs/EWrefs (lambda) - median/std dev: " + str(mvaleur) +"/" +str(valeur)
plt.title(tit)
plot(totline, totEWs/totEWrefs,'o')
plt.axhline(1.,c="r")
nomsave= DRS + "_overall_EWratio_lambda.jpg"
plt.savefig(nomsave,dpi=300)
if (montrer == 'montremoi'):
    plt.show()
plt.close()

ymin=-1.
ymax=3.
plt.figure(figsize=(16,6))
index = where (np.abs(totEWs/totEWrefs-1) < 2.)
valeur=np.std(totEWs[index]/totEWrefs[index])
mvaleur=median(totEWs[index]/totEWrefs[index])
tit= DRS + " EWs/EWrefs (lambda) - median/std dev: " + str(mvaleur) +"/" +str(valeur)
plt.title(tit)
plt.ylim(ymin,ymax)
plot(totline, totEWs/totEWrefs,'o')
plt.axhline(1.,c="r")
nomsave= DRS + "_overall_EWratio_lambda_zoom.jpg"
plt.savefig(nomsave,dpi=300)
if (montrer == 'montremoi'):
    plt.show()
plt.close()

plt.figure(figsize=(16,6))
tit= DRS + " EWs/EWrefs (EWrefs)"
plt.title(tit)
plot(totEWrefs, totEWs/totEWrefs,'o')
plt.axhline(1.,c='r')
nomsave= DRS + "_overall_EWratio_EW.jpg"
plt.savefig(nomsave,dpi=300)
if (montrer == 'montremoi'):
    plt.show()
plt.close()

plt.figure(figsize=(16,6))
tit= DRS + " EWs/EWrefs (EWrefs) zoom"
plt.title(tit)
plt.ylim(ymin,ymax)
plot(totEWrefs, totEWs/totEWrefs,'o')
plt.axhline(1.,c='r')
nomsave= DRS + "_overall_EWratio_EW_zoom.jpg"
plt.savefig(nomsave,dpi=300)
if (montrer == 'montremoi'):
    plt.show()
plt.close()



plt.figure(figsize=(16,6))
valeur = median(sigmaposi[2:39])
"""
if (DRS == 'LE'):
    valeur=median(sigmaposi[2:37])
"""
valeur=median(sigmaposi[2:37])
tit= DRS + " sigmaposi (order_centralwavelength) - mediane:" + str(valeur)
plt.title(tit)
ylabel('km/s')
xmin=3700.
xmax=10000.
plt.xlim(xmin,xmax)
plot(orderwave, sigmaposi,'o')
nomsave= DRS +"_sigmaposi.jpg"
plt.savefig(nomsave,dpi=300)
if (montrer == 'montremoi'):
    plt.show()
plt.close()

ymin=0.
ymax=3.
plt.figure(figsize=(16,6))
tit= DRS + " sigmaposi (order_centralwavelength) zoom"
plt.title(tit)
ylabel('km/s')
xmin=3700.
xmax=10000.
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plot(orderwave, sigmaposi,'o')
nomsave= DRS +"_sigmaposi_zoom.jpg"
plt.savefig(nomsave,dpi=300)
if (montrer == 'montremoi'):
    plt.show()
plt.close()

maxi=0.3
xxp=linspace(0,maxi,100)
yyp=xxp

plt.figure(figsize=(16,6))
tit= DRS + " EWs (EWrefs)"
plt.title(tit)
plot(totEWrefs, totEWs,'o')
plot(xxp,yyp,'-')
nomsave= DRS + "_overall_EWs_EWrefs.jpg"
plt.savefig(nomsave,dpi=300)
if (montrer == 'montremoi'):
    plt.show()
plt.close()

