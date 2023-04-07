#ce code genere les snippets utilises pour la calibration 2D
#ce code doit comparer deux spectres:
# 1) ThAr de reference UVES, calibre en longueur d'onde au mieux
#       voir http://astronomy.swin.edu.au/~mmurphy/thar/index.html
#  ici comparaison non pas au spectre même, mais à la liste de raies correspondante!
# 2) ThAr reduit par hobo sur Neo-Narval
# au depart on s'est servi de artlambda.dat (genere par arthols.py) de la même nuit,
# afin d'avoir une grille en lambda approximative
#


import astropy.io.fits as pyfits
from astropy.utils.decorators import lazyproperty
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
import json
import scipy
import util
from scipy.optimize import curve_fit
import extract

### all constants are in settings.py

# try:
#     settingsmodule = os.environ['SETTINGSMODULE']
# except:
#     raise Exception('wrong SETTINGSMODULE')

# try:
#     exec('from %s import *'%(settingsmodule,))
# except:
#     raise Exception('could not import {settingsmodule}')

# with open(REF_SPECTRUM, 'r') as f:
#     lines = f.readlines()

# num_lines=len(lines)
# npt=num_lines-2
# linecount = 0
# refwave = np.zeros (npt, dtype =  float)
# refintens = np.zeros (npt, dtype =  float)

# for l in range (0,npt):
#     linefr = lines[2+l]
#     sfr = linefr.split()
#     refwave[l] = float(sfr[0])
#     refintens[l] = float(sfr[1])


# def _snippets(extractor,nvoie,order):

#     """
#     #lecture du fichier ascii du spectre ThAr de reference (UVES)
#     #ce spectre contient une premiere colonne lambda, une seconde intensite,...
#     #213358 lignes

#     #refname = 'thar_spec_MM201006.dat'
#     #thar
#     #213358 4
#     #3033.0867480159  103.99224526863826    0.36633407377209    0.36633407377209     97.132933537024
#     """
#     #nombre de lignes du catalogue de raies
#     #lecture du fichier ascii catalogue de raies thar_UVES_MM090311.dat
#     #33322.4046  3000.109256  0.450 Ar II  W
#     #33319.9704  3000.328442 -0.007 XX 0   P
#     #atlasname = 'thar_UVES_MM090311.dat'
#     atlasline = []
#     with open(REF_ATLASLINES, 'r') as f:
#         alines = f.readlines()

#     # extract information...
#     atlasline = np.array([float(l.split()[1]) for l in alines])
#     #num_alines=len(alines)
#     #npt=num_alines
#     #atlasline = np.zeros (npt, dtype =  float)
#     #for l in range (0,npt):
#     #    linear = alines[l]
#     #    sfr = linear.split()
#     #    atlasline[l] = float(sfr[1])

#     # atlasline

#     o = order
#     if nvoie == 1:
#         lam, _flux, I = extractor.get_lambda_intens1(o)
#         bare_voie = extractor.bare_voie1(o)
#     elif nvoie == 2:
#         lam, _flux, I = extractor.get_lambda_intens2(o)
#         bare_voie = extractor.bare_voie2(o)
#     elif nvoie == 3:
#         lam, _flux, I = extractor.get_lambda_intens3(o)
#         bare_voie = extractor.bare_voie3(o)
#     else:
#         raise Exception('no such voie')

#     flux = np.zeros(NROWS)
#     flux[I] = _flux[I]

#     res = []


#     #####

#     #selectionner les raies de reference dans l'intervalle spectral
#     #lam[-1] deniere valeur de lambda

#     minlambda = np.min(lam[I])
#     maxlambda = np.max(lam[I])
#     indexx=np.where((minlambda < atlasline) & (atlasline < maxlambda))
#     atlasext=atlasline[indexx]


#     # ICI ON DOIT LIRE exlude.dat et filter les atlasline d'office
#     # cad ne pas prendre les raies tombant dans des endroits interdits
#     exclusion = np.loadtxt(EXCLUSION)

#     I, = np.where(exclusion[:,0] == o)

#     goodlines = []
#     for l in atlasext:

#         ##for i in I:
#         ##    if l >= exclusion[i,1] and l <= exclusion[i,2]:
#         ##       break

#         goodlines.append(l)


#     atlasext=np.array(goodlines)
#     #(atlasext)


#     """
#     plt.figure(figsize=(16,6))
#     plt.plot(lam,flux,"b")
#     plt.plot(refwave,-refintens/np.max(refintens),"r")
# #    plt.show()
#     """

#     """
#     for ll in atlasext:
#         plt.vlines(ll,0.,20000.,'y')
#     """

#     # on ne veut choisir que les raies de l'atlas qui se retrouvent dans la zone atlasext[k] +/- vrange, et qui ont un flux max au dessus du seuil. Ca reduit la liste.

#     latlasext=atlasext*(1.-VRANGE/C_LIGHT)
#     ratlasext=atlasext*(1.+VRANGE/C_LIGHT)

#     numlines=atlasext.shape[0]
#     maxi = np.zeros (numlines, dtype =  float)
#     maxir = np.zeros (numlines, dtype =  float)

#     #ici on selectionne des snippets autour de chaque raie du catalogue
#     #de reference. Il faut que la grille en longueur d'onde soit des le
#     #depart suffisamment bonne pour que la raie observee tombe dans le snippet


#     snip = []
#     for l,c,r in zip(latlasext,atlasext,ratlasext):
#         indext,  =  np.where((lam > l) & (lam < r))
#         wave=lam[indext]
#         inte=flux[indext]
#         bare_inte = bare_voie[indext]

#         """
#         indextr, = np.where((refwave > l) & (refwave < r))
#         waver = refwave[indextr]
#         inter = refintens[indextr]
#         print(c,wave,inte)
#         """


#     # selectionner que les raies du ThAr observes au dessus du seuil.
#     # pour chaque raie k on determine le maximum de flux
#         distmax = 0   ## TODO: make global constant
#         goodsnippet = True
#    ###     goodsnippet = goodsnippet and (np.max(inte) - np.min(inte)) >= SEUIL
#    #     goodsnippet = goodsnippet and (np.max(inter) - np.min(inter)) >= SEUILR
#         try:
#             goodsnippet = goodsnippet \
#                 and (np.argmax(inte)>= distmax) \
#                 and (np.argmax(inte) <= inte.shape[0]-distmax)
#         except:
#             print('order ', o, 'problem line')
#             continue
#         if goodsnippet:
#             snip.append({
#                 "true_order_number": o,
#                 "ref_lambda": c ,
#                 "pixels_extract": indext,
#                 "pixel_mean" : np.mean(indext),
#                 "wave": wave,
#                 "reduced_flux_values_extract": inte,
#                 "flux_values_extract" : bare_inte,
#             })
#             """
#             plt.vlines(c,-10.,10.,'y')
#             plt.plot(wave,inte,"y")

#     plt.show()
#     """

#     return pd.DataFrame(snip)

# def snippets(extractor,nvoie,orders):
#     snipets = []

#     for o in np.array(orders):
#         tmp = _snippets(extractor, nvoie, o)
#         snipets.append(tmp)
#         print('snippets: order {} nr of snippets {}'.format(o, len(tmp)) )
#     return pd.concat(snipets, ignore_index=True, axis=0)


class Snippets:

    def __init__(self, voie, extractor=None):
        self.voie = voie
        self._atlasline = None
        self._snippets = None
        #self._lam=None
        #self._I=None
        #self._flux=None
        #self._bare_voie=None
        #self._snippets=None
        #self._lines_voie=None

        self.extractor = extractor

        self.NROWS = self.extractor.NROWS
        self.ORDERS = self.extractor.ORDERS

    """
    def _prepare(self):
        self._lam={}
        self._I={}
        self._flux={}
        self._bare_voie={}
        # self._lines_voie=

        for o in self.ORDERS:

            if self.voie == 1:
                lam, _flux, I = self.extractor.get_lambda_intens1(o)
                bare_voie = self.extractor.ba_voie1[o]
            elif self.voie == 2:
                lam, _flux, I = self.extractor.get_lambda_intens2(o)
                bare_voie = self.extractor.ba_voie2[o]
            elif self.voie == 3:
                lam, _flux, I = self.extractor.get_lambda_intens3(o)
                bare_voie = self.extractor.ba_voie3[o]
            else:
                raise Exception('no such voie')


            flux = np.zeros(self.NROWS)
            flux[I] = _flux[I]

            self._lam[o] = lam
            self._I[o] = I
            self._flux[o] = flux
            self._bare_voie[o] = bare_voie
    """
    @property
    def kwargs(self):
        return self.extractor.kwargs

    @property
    def atlasline(self):
        if not self._atlasline is None:
            return self._atlasline

        with open(self.kwargs['REF_ATLASLINES'], 'r') as f:
            alines = f.readlines()

        # extract information...
        atlaslines =  [float(l.split()[1]) for l in alines]
        # self._atlasline = atlaslines
        self._atlasline = pd.DataFrame({'atlas_line': pd.Series(atlaslines)})
        return self._atlasline

    def lambda_range(self, o):
        return  self.pix_to_lambda_map[o](self.extractor.beams[o].pixel_range)

    def atlasext(self, o):
        """
        the extracted lines of atlas
        """
        minlambda, maxlambda = self.lambda_range(o)

        indexx = (minlambda < self.atlasline['atlas_line']) & (self.atlasline['atlas_line'] < maxlambda)
        tmp = self.atlasline[indexx]

        exclusion = np.loadtxt(self.kwargs['EXCLUSION'])

        I, = np.where(exclusion[:,0] == o)
        exc = exclusion[I]

        j = pd.Series(True, index=tmp.index)
        for e in exc:
            j = j & ( (tmp['atlas_line'] < e[1]) | (tmp['atlas_line'] > e[2]) )
        return tmp[j]


    def _lines(self, o):
        NMAX_LINES = 50 # maximal number of lines to extract
        v = self.bare_voie(o)
        lm = util.local_maxima(v)
        bs = []
        for i,xab in enumerate ( lm[:min(NMAX_LINES, len(lm))]):
            x, a, b = xab
            if a==b:
                continue
            s = v[np.arange(a,b+1)]
            if len(s)==0:
                continue
            try:
                A, mu, sigma, offset = util.estimate_location(s, **self.kwargs)
            except Exception as ex:
                print(s, ex)
                continue
            bs.append( {
                'true_order_number': o, 
                'index_pixel_snippet': i, 
                'posmax': x,
                'left': a, 
                'right': b, 
                'pixel_A': A,
                'pixel_mu': mu,
                'pixel_sigma': sigma,
                'pixel_offset': offset,
                'pixel_mean': a + mu,
                'pixel_std' : 1./ np.sqrt(np.sum(s)),  
                'pixel_sum_intens': np.sum(s),
                'pixel_max_intens': v[x],
                'bare_voie': s,
                'bootstraped': False
            })
        return pd.DataFrame(bs)

    def bare_voie(self, o):
        if self.voie == 1:
            return self.extractor.bare_voie1[o]
        elif self.voie == 2:
            return self.extractor.bare_voie2[o]
        else:
            raise Exception('line 3 not implemented')      

    #@lazyproperty
    #def lines_voie(self):
    #    res ={o: self._lines(self.ba_voie(o)) for o in self.ORDERS }
    #    return res



    @property
    def pix_to_lambda_map(self):
        if self.voie == 1:
            return self.extractor.pix_to_lambda_map_voie1
        elif self.voie == 2:
            return self.extractor.pix_to_lambda_map_voie2
        else:
            raise Exception('voie 3 not yet implemented')

    def _snippet(self, o):

        atlasext = self.atlasext(o)
        bare_voie = self.bare_voie(o)
        lines_voie = self._snippets.loc[self._snippets['true_order_number'] == o]

        lams = self.pix_to_lambda_map[o](lines_voie['pixel_mean'].values)

        for ic, c in atlasext.itertuples():


            l = (1.-self.kwargs['VRANGE']/self.kwargs['C_LIGHT']) * c
            r = (1.+self.kwargs['VRANGE']/self.kwargs['C_LIGHT']) * c

            I, = np.where( np.logical_and ( l <= lams, lams <= r))
            goodlam = len(I) == 1 ## exacly one peak in interval
            goodlam = goodlam and True ## SEUIL in extract TODO

            if not goodlam:
                continue

            i = I[0]  # index of single good snippet belonging to l

            a = lines_voie['left'].iloc[i]
            b = lines_voie['right'].iloc[i]
            idx = lines_voie.index[i]

            intens = bare_voie[a:b+1]

            # bootstrapping
            if self.kwargs['n_bootstrap'] > 1:
                mu, s, sample = util.bootstrap_estimate_location(intens, **self.kwargs)
            else:
                mu = util.estimate_location(intens, **self.kwargs)[1]
                s = 1/np.sqrt(sum(intens))
            self._snippets.loc[idx,'goodsnippet'] = True
            self._snippets.loc[idx,'ref_lambda'] = c
            self._snippets.loc[idx,'pixel_mean'] = a+mu
            self._snippets.loc[idx,'pixel_std'] = s
            self._snippets.loc[idx,'catalog_index'] = int(ic)
    
    @property 
    def snippets(self):
        if not self._snippets is None:
            return self._snippets
        self.extractor.logging('computing snippets for voie '+str(self.voie))
        
        # generating all possible snippets
        tmp = []
        for i, o in enumerate(self.ORDERS):
            util.progress(i, len(self.ORDERS))
            tmp.append(self._lines(o))
        tmp = pd.concat(tmp, ignore_index=True, axis=0)
        self._snippets = tmp

        self._snippets['goodsnippet'] = False
        self._snippets['ref_lambda'] = 0.0
        self._snippets['catalog_index'] = 0

        for i, o in enumerate(self.ORDERS):
            util.progress(i, len(self.ORDERS))
            self._snippet(o)
        return self._snippets

    def index_good_order(self, o):
        """
        indices of good snippets of order o
        """
        return self.snippets['true_order_number']==o & self.good_snippet

    

    @lazyproperty
    def overlapping_snippets(self):
        tmp = []
        for o, oo in zip(self.ORDERS[:-1], self.ORDERS[1:]):

            Io = self.snippets['true_order_number'] == o
            Ioo = self.snippets['true_order_number'] == oo
            pp1 = self.snippets.loc[Io, 'pixel_mean']
            pp2 = self.snippets.loc[Ioo,'pixel_mean']


            p1 = self.pix_to_lambda_map[o](pp1)
            p2 = self.pix_to_lambda_map[oo](pp2)

            # TODO add global constant
            f = 0.2*self.kwargs['VRANGE']/self.kwargs['C_LIGHT']

            d = abs(p1[:, np.newaxis]-p2[np.newaxis,:])

            mm = (d < f*p1[:,np.newaxis])  &  (d < f * p2[np.newaxis,:])
            mm = mm & (np.count_nonzero(mm, axis=1)==1)[:, np.newaxis] & (np.count_nonzero(mm, axis=0)==1)[np.newaxis, :]
            I, J = np.where( mm )

            tmp.append(pd.DataFrame({
                'true_order_number_o': o,
                'true_order_number_oo': oo,
                'index_line_o':  I,
                'index_line_oo': J,
                #'pixel_mean_o' :pp1.values[I],
                #'pixel_mean_oo': pp2.values[J],
                #'pixel_std_o': self.snippets[Io]['pixel_std'].values[I],
                #'pixel_std_oo': self.snippets[Ioo]['pixel_std'].values[J]
                }
            ))

        return pd.concat(tmp).reset_index(drop=True)

    ## convenience methods
    
    @property
    def good_snippet(self):
        return self.snippets['goodsnippet'] == True
    
    @property
    def _x(self):
        return self.snippets['pixel_mean']
    
    @property
    def x(self):
        return self._x[self.good_snippet]
    
    @property
    def _o(self):
        return self.snippets['true_order_number']
    
    @property
    def o(self):
        return self._o[self.good_snippet]

    @property
    def _l(self):
        return self.snippets['ref_lambda']
    
    @property
    def l(self):
        return self._l[self.good_snippet]
      
    @property
    def _sigma(self):
        return self.snippets['pixel_std']
    
    @property
    def sigma(self):
        return self._sigma[self.good_snippet]

    @property
    def _sn(self):
        return self.snippets

    @property
    def sn(self):
        return self.snippets[self.good_snippet]
    
    @property
    def ci(self):
        return self.snippets[self.good_snippet]['catalog_index']    


    @property
    def ngood(self):
        return sum(self.good_snippet)

if __name__ == "__main__":
    import os
    from settingspegasus import kwargs
    import extract
    myext = extract.get_ext(
        '/home/hols/vega2021/thar/51Peg_raw/2020_1029/NEO_20201029_172538_th0.fits',
        **kwargs)
    myext.update_kwargs(**kwargs)
    snip = Snippets(voie=1, extractor=myext)
