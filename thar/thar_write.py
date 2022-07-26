#ce code genere fichiers .json contenant les snippets utilises pour la calibration 2D
#ce code doit comparer deux spectres:
# 1) ThAr de reference UVES, calibre en longueur d'onde au mieux
#       voir http://astronomy.swin.edu.au/~mmurphy/thar/index.html
#  ici comparaison non pas au spectre même, mais à la liste de raies correspondante!
# 2) ThAr reduit par le DRS de Neo-Narval
#


import astropy.io.fits as pyfits
import numpy as np
import sys
import os
import json
import scipy
from scipy.optimize import curve_fit
#


clum=3e5

def thar_write(refname,atlasname,fitsfile,outfile1,outfile2,outfile3):
    """
    refname: fichier ascii du spectre ThAr de reference (UVES)
    atlasname: fichier ascii des raies retenues ThAr de reference
    fitsfile: fichier input Neo-Narval ThAr (th2)
    outfile1 (2,3): output json file
    """
    #lecture du fichier ascii du spectre ThAr de reference (UVES)
    #ce spectre contient une premiere colonne lambda, une seconde intensite,...
    #213358 lignes

    #refname = 'thar_spec_MM201006.dat'
    #thar
    #213358 4
    #3033.0867480159  103.99224526863826    0.36633407377209    0.36633407377209     97.132933537024

    #nombre de lignes du catalogue de raies
    num_lines=len(open(refname,'r').readlines())
    npt=num_lines-2
    linecount = 0
    refwave = np.zeros (npt, dtype =  float)
    refintens = np.zeros (npt, dtype =  float)


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
    #atlasname = 'thar_UVES_MM090311.dat'
    num_alines=len(open(atlasname,'r').readlines())
    npt=num_alines
    atlasline = np.zeros (npt, dtype =  float)
    ar = open(atlasname, 'r')
    for l in range (0,npt):
            linear = ar.readline()
            sfr = linear.split()
            atlasline[l] = float(sfr[1])
    fr.close

    a=pyfits.open(fitsfile)
    wave1=a[1].data['Wavelength1']
    intens1=a[1].data['Beam1']
    wave2=a[1].data['Wavelength2']
    intens2=a[1].data['Beam2']
    wave3=a[1].data['Wavelength3']
    intens3=a[1].data['Beam3']

    orderlim=a[2].data['Orderlimit']

    def gauss(x, A, mu, sigma, y_offset):
    # A, sigma, mu, y_offset = p
       return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y_offset


    # orderlim va de 0 a 304200, 304200/7800= 39. la 41eme valeur est a 0. Il y a donc 39 ordres en tout
    # in range 38,39 veut dire que tout ce qui est en dessous de 39, a partir de 0, donc 38
    # uniquement

    ordernum = np.zeros (orderlim.size, dtype = float)
    centrallam = np.zeros (orderlim.size, dtype = float)
    tordernum = np.zeros (orderlim.size, dtype = float)

    res1=[]
    res2=[]
    res3=[]


    #for j in range (19,20):
    #for j in range (31,orderlim.size-2):
    for j in range (2,orderlim.size-2):

        ordernum[j] = j
    # vrai numero d'ordre
        tordernum[j] = ordernum[j]+21
    #    print("j: ",j, " ordernum[j]", ordernum[j], " true ordernum[j]: ",tordernum[j])
    #ici on ne selectionne que les raies de l'atlas dans l'ordre j, index.size est ce nombre de raies du catalogue UVES
        indexx=np.where((wave1[orderlim[j]]<atlasline) & (atlasline < wave1[orderlim[j+1]-1]))
        atlasext=atlasline[indexx]
        wext1=wave1[orderlim[j]:orderlim[j+1]]
        iext1=intens1[orderlim[j]:orderlim[j+1]]
        wext2=wave2[orderlim[j]:orderlim[j+1]]
        iext2=intens2[orderlim[j]:orderlim[j+1]]
        wext3=wave3[orderlim[j]:orderlim[j+1]]
        iext3=intens3[orderlim[j]:orderlim[j+1]]
        centrallam[j]=(wave1[orderlim[j+1]-1]+wave1[orderlim[j]])/2.
    #
    #orderlim[j+1]-orderlim[j] vaut 7800, il y a donc 7800 pixels (="Arturos") par ordre
    #ces pixels ne sont pas les pixels du CCD mais correspondent a un sousechantillonnage afin d arriver a la superresolution
        
        
    # seuil en ADU pour qu'on puisse selectionner une raie de l'atlas. Il faut qu'elle depasse de ce seuil le zero. Ensuite, il faut donner une zone +/- autour de la raie du catalogue, pour detecter la valeur maxi du flux tu ThAr obs.
    # le vrange est donc cette largeur en km/s
        seuil = 100.
        vrange = 9.
    # on ne veut choisir que les raies de l'atlas qui se retrouvent dans la zone atlasext[k] +/- vrange, et qui ont un flux max au dessus du seuil. Ca reduit la liste. En plus, il faut le faire sur les deux voies et appliquer le mini.
        
        latlasext=atlasext*(1.-vrange/clum)
        ratlasext=atlasext*(1.+vrange/clum)
        
        numlines=int(np.array(indexx).size)
        maxi = np.zeros (numlines, dtype =  float)
        maxir = np.zeros (numlines, dtype =  float)
        intline1 = np.zeros (numlines, dtype =  float)
        intline2 = np.zeros (numlines, dtype =  float)
        
    # selectionner que les raies du ThAr observes au dessus du seuil.
        for k in range (numlines):
            indext1 =  np.where((wext1<ratlasext[k]) & (wext1>latlasext[k]))
            w1=wext1[indext1]
            i1=iext1[indext1]
        
            indext2 =  np.where((wext2<ratlasext[k]) & (wext2>latlasext[k]))
            w2=wext2[indext2]
            i2=iext2[indext2]
            
            indext3 =  np.where((wext3<ratlasext[k]) & (wext3>latlasext[k]))
            w3=wext3[indext3]
            i3=iext3[indext3]
            
            indextr =  np.where((refwave<ratlasext[k]) & (refwave>latlasext[k]))
            wr=refwave[indextr]
            ir=refintens[indextr]
      
      # attention: j'avais ici modifie en maxi[k]=min(max(i1),max(i2),max(i3)), mais sans flux sur la voie 3
      # ca enleve toutes les raies!!
            maxi[k]=min(max(i1),max(i2))
            maxir[k] = max(ir)
            
            
        indseuil= np.where((maxi > seuil) & (maxir > seuil))
        #number of selected reference lines
        nslines=int(np.array(indseuil).size)
        satlasext=atlasext[indseuil]
        slatlasext=latlasext[indseuil]
        sratlasext=ratlasext[indseuil]
        
    # ici boucle sur toutes les raies selectionnees, correlation sur une petite zone autour de chaque raie entre l'intensite des deux voies individuellement et l'intensite du spectre de ref.
     
        resol1 = np.zeros (nslines, dtype =  float)
        resol2 = np.zeros (nslines, dtype =  float)
        resol3 = np.zeros (nslines, dtype =  float)
        gaussposi1 = np.zeros (nslines, dtype =  float)
        gaussposi2 = np.zeros (nslines, dtype =  float)
        gaussposi3 = np.zeros (nslines, dtype =  float)
        pixelposi1 = np.zeros (nslines, dtype =  float)
        pixelposi2 = np.zeros (nslines, dtype =  float)
        pixelposi3 = np.zeros (nslines, dtype =  float)
        erreurposi1 = np.zeros (nslines, dtype =  float)
        erreurposi2 = np.zeros (nslines, dtype =  float)
        erreurposi3 = np.zeros (nslines, dtype =  float)
        amp1 = np.zeros (nslines, dtype =  float)
        amp2 = np.zeros (nslines, dtype =  float)
        amp3 = np.zeros (nslines, dtype =  float)
        
        for k in range (nslines):
    # attention, si on teste pour seulement 2 raies, le programme se plante plus tard!
    #    for k in range (5,6):

    # on incremente le compteur general de raies spectrales de 1
    # ici on pretend que voies 1 et 2 sont suffisamment proches pour que le nombre
    # de raies est identique
    #        print("linecount: ",linecount)
            linecount += 1
            indextr = np.where((refwave <sratlasext[k]) & (refwave>slatlasext[k]))
            wr=refwave[indextr]
            ir=refintens[indextr]
           
            indext1 =  np.where((wext1<max(wr)) & (wext1>min(wr)))
            w1=wext1[indext1]
            i1=iext1[indext1]
            pix1=np.array(indext1)
            pixel1=pix1[0]
        
            indext2 =  np.where((wext2<max(wr)) & (wext2>min(wr)))
            w2=wext1[indext2]
            i2=iext2[indext2]
            pix2=np.array(indext2)
            pixel2=pix2[0]

            indext3 =  np.where((wext3<max(wr)) & (wext3>min(wr)))
            w3=wext3[indext3]
            i3=iext3[indext3]
            pix3=np.array(indext3)
            pixel3=pix3[0]

    ##############################resolution voie 1#############################################
                
            xxx=w1
            yyy=i1
    #        print("satlasext[k]", satlasext[k])
    #        print("xxx:", xxx)
    #        print("yyy:", yyy)
    #        print("pixel1:",pixel1)
    #        print("np.mean(xxx):",np.mean(xxx))
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
    #        print ("rinit_vals fit EWref:", rinit_vals)
            try:
                rgopt, rcovar = curve_fit(gauss,xxx-np.mean(xxx),yyy,p0=rinit_vals)
                rgoptp, rcovarp = curve_fit(gauss,pixel1-np.mean(pixel1),yyy,p0=rinit_valsp)
    # la matrice de covariance rcovarp (en pixels) contient sur sa diagonale les erreurs associes
    # aux differents parametres. L'erreur sur la position du fit, donc relatif a rgoptp[1]
    # s'exprime alors
                errs=np.diag(rcovarp)
                erreurposi1[k] = errs[1]
    #            print("rmu:", gopt[1]," rA:",gopt[0]," rSigma:",gopt[2]," rY_offset:",gopt[3])
            #Return co-effs for fit and covariance
    #            plt.figure()
    #            plt.plot(xxx,yyy)
    #            plt.plot(xxx,gauss(xxx-np.mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]),"r")
    #            plt.show()
                FWHM =  rgopt[2] * 2 * np.sqrt(2*np.log(2))
                resol = np.mean(xxx)/FWHM
                resol1[k]=resol
                amp1[k]=rgopt[0]
                gaussposi1[k] = rgopt[1]+np.mean(xxx)
                pixelposi1[k] = rgoptp[1]+np.mean(pixel1)
    #            print("gaussposi1[k]:", gaussposi1[k])
    #            print("pixelposi1[k]:", pixelposi1[k])
                

    # ici on ecrit le tableau de sortie (i, lam_i, true_ordre_i, posi_i,err_posi_i)
    #            outtableau1.write("%4d %20.10e %2d %20.10e %20.10e" % (linecount, gaussposi1[k], tordernum[j], pixelposi1[k], erreurposi1[k])+str(pixel1)+str(yyy)+"\n")
                res1.append({
                    'ThAr_line_number': linecount,
                    'ref_lambda': satlasext[k],
                    'meas_lambda': gaussposi1[k],
                    'true_order_number': tordernum[j],
                    'meas_pixel_position': pixelposi1[k],
                    'fit_error_pixel_position': erreurposi1[k],
                    'pixels_extract': pixel1.tolist(),
                    'flux_values_extract': yyy.tolist()
                })
                """
                print("Resolution:", resol)
                print("gauss:",gauss(xxx-np.mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]))
                print("gopt",rgopt)
                """
#            except (RuntimeError,TypeError):
            except Exception as ex:
                print("Ca n'a pas marche voie 1", ex, np.mean(xxx))
    #[gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
    # EW= A*(np.sqrt(2pi) sigma) si le continu est a 1!! Donc renormalisation AVANT
    # la division par gopt[3] est une normalisation locale


    ##############################resolution voie 2#############################################
                
            xxx=w2
            yyy=i2
    #        print("satlasext[k]", satlasext[k])
    #        print("xxx:", xxx)
    #        print("yyy:", yyy)
    #        print("pixel2:",pixel2)
    #        print("np.mean(xxx):",np.mean(xxx))
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
    #        print ("rinit_vals fit EWref:", rinit_vals)
            try:
                rgopt, rcovar = curve_fit(gauss,xxx-np.mean(xxx),yyy,p0=rinit_vals)
                rgoptp, rcovarp = curve_fit(gauss,pixel2-np.mean(pixel2),yyy,p0=rinit_valsp)
                errs=np.diag(rcovarp)
                erreurposi2[k] = errs[1]

                
    #            print("rmu:", gopt[1]," rA:",gopt[0]," rSigma:",gopt[2]," rY_offset:",gopt[3])
            #Return co-effs for fit and covariance
    #            plt.figure()
    #            plt.plot(xxx,yyy)
    #            plt.plot(xxx,gauss(xxx-np.mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]),"r")
    #            plt.show()
                FWHM =  rgopt[2] * 2 * np.sqrt(2*np.log(2))
                resol = np.mean(xxx)/FWHM
                resol2[k]=resol
                amp2[k]=rgopt[0]
                gaussposi2[k] = rgopt[1]+np.mean(xxx)
                pixelposi2[k] = rgoptp[1]+np.mean(pixel2)
    #            print("k: ",k)
    #            print("gaussposi2[k]:", gaussposi2[k])
    #            print("pixelposi2[k]:", pixelposi2[k])
                # ici on ecrit le tableau de sortie (i, lam_i, true_ordre_i, posi_i,err_posi_i)
    #        outtableau2.write("%15.7f %20.10e %20.10e %20.10e %20.10e\n" % (linecount, gaussposi2[k], tordernum[j], pixelposi2[k], erreurposi2[k])+str(pixel2)+str(yyy)+"\n")
                res2.append({
                    'ThAr_line_number': linecount,
                    'ref_lambda': satlasext[k],
                    'meas_lambda': gaussposi2[k],
                    'true_order_number': tordernum[j],
                    'meas_pixel_position': pixelposi2[k],
                    'fit_error_pixel_position': erreurposi2[k],
                    'pixels_extract': pixel2.tolist(),
                    'flux_values_extract': yyy.tolist()
                })
                
            except Exception as ex:
                print("Ca n'a pas marche voie 2", ex, np.mean(xxx))
    #[gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
    # EW= A*(np.sqrt(2pi) sigma) si le continu est a 1!! Donc renormalisation AVANT
    # la division par gopt[3] est une normalisation locale


    ##############################resolution voie 3#############################################
                
            xxx=w3
            yyy=i3
    #       print("satlasext[k]", satlasext[k])
    #        print("xxx:", xxx)
    #        print("yyy:", yyy)
    #        print("pixel2:",pixel2)
    #        print("np.mean(xxx):",np.mean(xxx))
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
            rsigmap=0.2*(np.argmax(pixel3)-np.argmin(pixel3))
            rinit_valsp = [rA, rmu, rsigmap, ry_offset]
    #         print ("rinit_vals fit EWref:", rinit_vals)
            try:
                rgopt, rcovar = curve_fit(gauss,xxx-np.mean(xxx),yyy,p0=rinit_vals)
                rgoptp, rcovarp = curve_fit(gauss,pixel3-np.mean(pixel3),yyy,p0=rinit_valsp)
                errs=np.diag(rcovarp)
                erreurposi3[k] = errs[1]

                
    #            print("rmu:", gopt[1]," rA:",gopt[0]," rSigma:",gopt[2]," rY_offset:",gopt[3])
            #Return co-effs for fit and covariance
    #            plt.figure()
    #            plt.plot(xxx,yyy)
    #            plt.plot(xxx,gauss(xxx-np.mean(xxx),rgopt[0],rgopt[1],rgopt[2],rgopt[3]),"r")
    #            plt.show()
                FWHM =  rgopt[2] * 2 * np.sqrt(2*np.log(2))
                resol = np.mean(xxx)/FWHM
                resol3[k]=resol
                amp3[k]=rgopt[0]
                gaussposi3[k] = rgopt[1]+np.mean(xxx)
                pixelposi3[k] = rgoptp[1]+np.mean(pixel2)
    #             print("k: ",k)
    #             print("gaussposi3[k]:", gaussposi3[k])
    #             print("pixelposi3[k]:", pixelposi3[k])
                # ici on ecrit le tableau de sortie (i, lam_i, true_ordre_i, posi_i,err_posi_i)
    #        outtableau2.write("%15.7f %20.10e %20.10e %20.10e %20.10e\n" % (linecount, gaussposi2[k], tordernum[j], pixelposi2[k], erreurposi2[k])+str(pixel2)+str(yyy)+"\n")
                res3.append({
                    'ThAr_line_number': linecount,
                    'ref_lambda': satlasext[k],
                    'meas_lambda': gaussposi3[k],
                    'true_order_number': tordernum[j],
                    'meas_pixel_position': pixelposi3[k],
                    'fit_error_pixel_position': erreurposi3[k],
                    'pixels_extract': pixel3.tolist(),
                    'flux_values_extract': yyy.tolist()
                })
            except Exception as ex:
                print("Ca n'a pas marche voie 3", ex, np.mean(xxx))
    #[gopt[0]=A,   gopt[1]= mu, gopt[2]=sigma, gopt[3] = y_offset]
    # EW= A*(np.sqrt(2pi) sigma) si le continu est a 1!! Donc renormalisation AVANT
    # la division par gopt[3] est une normalisation locale

    with open(outfile1, 'w') as outfile:
        json.dump(res1, outfile, indent=2)
    with open(outfile2, 'w') as outfile:
        json.dump(res2, outfile, indent=2)
    with open(outfile3, 'w') as outfile:
        json.dump(res3, outfile, indent=2)

if __name__ == "__main__":
    refname = "./reffiles/thar_spec_MM201006.dat"
    atlasname = "./reffiles/thar_UVES_MM090311.dat"
    fitsfile = "./datafiles/NEO_20220219_173048_th2.fits"
    outfile1 = "NEO_20220219_173048_th2_voie1.json"
    outfile2 = "NEO_20220219_173048_th2_voie2.json"
    outfile3 = "NEO_20220219_173048_th2_voie3.json"
    
    thar_write(refname,atlasname,fitsfile,outfile1,outfile2,outfile3)
