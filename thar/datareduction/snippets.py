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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
import json
import scipy
from scipy.optimize import curve_fit
import extract

### all constants are in settings.py
settingsmodule = os.environ.get('SETTINGSMODULE', 'settings')

try:
    exec('from %s import *'%(settingsmodule,))
except:
    raise Exception('could not import {settingsmodule}')

with open(REF_SPECTRUM, 'r') as f:
    lines = f.readlines()

num_lines=len(lines)
npt=num_lines-2
linecount = 0
refwave = np.zeros (npt, dtype =  float)
refintens = np.zeros (npt, dtype =  float)

for l in range (0,npt):
    linefr = lines[2+l]
    sfr = linefr.split()
    refwave[l] = float(sfr[0])
    refintens[l] = float(sfr[1])


def _snippets(extractor,nvoie,order):

    """
    #lecture du fichier ascii du spectre ThAr de reference (UVES)
    #ce spectre contient une premiere colonne lambda, une seconde intensite,...
    #213358 lignes

    #refname = 'thar_spec_MM201006.dat'
    #thar
    #213358 4
    #3033.0867480159  103.99224526863826    0.36633407377209    0.36633407377209     97.132933537024
    """
    #nombre de lignes du catalogue de raies
    #lecture du fichier ascii catalogue de raies thar_UVES_MM090311.dat
    #33322.4046  3000.109256  0.450 Ar II  W
    #33319.9704  3000.328442 -0.007 XX 0   P
    #atlasname = 'thar_UVES_MM090311.dat'
    atlasline = []
    with open(REF_ATLASLINES, 'r') as f:
        alines = f.readlines()

    # extract information...
    atlasline = np.array([float(l.split()[1]) for l in alines])
    #num_alines=len(alines)
    #npt=num_alines
    #atlasline = np.zeros (npt, dtype =  float)
    #for l in range (0,npt):
    #    linear = alines[l]
    #    sfr = linear.split()
    #    atlasline[l] = float(sfr[1])

    # atlasline

    o = order
    if nvoie == 1:
        lam, _flux, I = extractor.get_lambda_intens1(o)
        bare_voie = extractor.bare_voie1(o)
    elif nvoie == 2:
        lam, _flux, I = extractor.get_lambda_intens2(o)
        bare_voie = extractor.bare_voie2(o)
    elif nvoie == 3:
        lam, _flux, I = extractor.get_lambda_intens3(o)
        bare_voie = extractor.bare_voie3(o)
    else:
        raise Exception('no such voie')

    flux = np.zeros(NROWS)
    flux[I] = _flux[I]

    res = []


    #####

    #selectionner les raies de reference dans l'intervalle spectral
    #lam[-1] deniere valeur de lambda

    minlambda = np.min(lam[I])
    maxlambda = np.max(lam[I])
    indexx=np.where((minlambda < atlasline) & (atlasline < maxlambda))
    atlasext=atlasline[indexx]


    # ICI ON DOIT LIRE exlude.dat et filter les atlasline d'office
    # cad ne pas prendre les raies tombant dans des endroits interdits
    exclusion = np.loadtxt(EXCLUSION)

    I, = np.where(exclusion[:,0] == o)

    goodlines = []
    for l in atlasext:

        ##for i in I:
        ##    if l >= exclusion[i,1] and l <= exclusion[i,2]:
        ##       break

        goodlines.append(l)


    atlasext=np.array(goodlines)
    #(atlasext)


    """
    plt.figure(figsize=(16,6))
    plt.plot(lam,flux,"b")
    plt.plot(refwave,-refintens/np.max(refintens),"r")
#    plt.show()
    """

    """
    for ll in atlasext:
        plt.vlines(ll,0.,20000.,'y')
    """

    # on ne veut choisir que les raies de l'atlas qui se retrouvent dans la zone atlasext[k] +/- vrange, et qui ont un flux max au dessus du seuil. Ca reduit la liste.

    latlasext=atlasext*(1.-VRANGE/C_LIGHT)
    ratlasext=atlasext*(1.+VRANGE/C_LIGHT)

    numlines=atlasext.shape[0]
    maxi = np.zeros (numlines, dtype =  float)
    maxir = np.zeros (numlines, dtype =  float)

    #ici on selectionne des snippets autour de chaque raie du catalogue
    #de reference. Il faut que la grille en longueur d'onde soit des le
    #depart suffisamment bonne pour que la raie observee tombe dans le snippet


    snip = []
    for l,c,r in zip(latlasext,atlasext,ratlasext):
        indext,  =  np.where((lam > l) & (lam < r))
        wave=lam[indext]
        inte=flux[indext]
        bare_inte = bare_voie[indext]

        """
        indextr, = np.where((refwave > l) & (refwave < r))
        waver = refwave[indextr]
        inter = refintens[indextr]
        print(c,wave,inte)
        """


    # selectionner que les raies du ThAr observes au dessus du seuil.
    # pour chaque raie k on determine le maximum de flux
        distmax = 0   ## TODO: make global constant
        goodsnippet = True
   ###     goodsnippet = goodsnippet and (np.max(inte) - np.min(inte)) >= SEUIL
   #     goodsnippet = goodsnippet and (np.max(inter) - np.min(inter)) >= SEUILR
        try:
            goodsnippet = goodsnippet \
                and (np.argmax(inte)>= distmax) \
                and (np.argmax(inte) <= inte.shape[0]-distmax)
        except:
            print('order ', o, 'problem line')
            continue
        if goodsnippet:
            snip.append({
                "true_order_number": o,
                "ref_lambda": c ,
                "pixels_extract": indext,
                "pixel_mean" : np.mean(indext),
                "wave": wave,
                "reduced_flux_values_extract": inte,
                "flux_values_extract" : bare_inte,
            })
            """
            plt.vlines(c,-10.,10.,'y')
            plt.plot(wave,inte,"y")

    plt.show()
    """

    return pd.DataFrame(snip)

def snippets(extractor,nvoie,orders):
    snipets = []

    for o in np.array(orders):
        tmp = _snippets(extractor, nvoie, o)
        snipets.append(tmp)
        print('snippets: order {} nr of snippets {}'.format(o, len(tmp)) )
    return pd.concat(snipets, ignore_index=True, axis=0)


class Snippets:

    def __init__(self, voie, tharfits, extractor=None, **kwargs):
        self.kwargs = kwargs
        self.voie = voie
        self.tharfits = tharfits
        self._atlasline = None
        self._lam=None
        self._I=None
        self._flux=None
        self._bare_voie=None
        self._snippets=None

        self.ORDERS = kwargs['ORDERS']
        self.NROWS = kwargs['NROWS']
        self.extractor = extractor if not extractor is None else extract.Extractor(**kwargs)
        self.extractor.set_fitsfile(self.tharfits)

    def _prepare(self):
        self._lam={}
        self._I={}
        self._flux={}
        self._bare_voie={}

        for o in self.kwargs.get('ORDERS', self.ORDERS):

            if self.voie == 1:
                lam, _flux, I = self.extractor.get_lambda_intens1(o)
                bare_voie = self.extractor.bare_voie1(o)
            elif self.voie == 2:
                lam, _flux, I = self.extractor.get_lambda_intens2(o)
                bare_voie = self.extractor.bare_voie2(o)
            elif self.voie == 3:
                lam, _flux, I = self.extractor.get_lambda_intens3(o)
                bare_voie = self.extractor.bare_voie3(o)
            else:
                raise Exception('no such voie')


            flux = np.zeros(self.NROWS)
            flux[I] = _flux[I]

            self._lam[o] = lam
            self._I[o] = I
            self._flux[o] = flux
            self._bare_voie[o] = bare_voie

    @property
    def atlasline(self):
        if not self._atlasline is None:
            return self._atlasline

        with open(self.kwargs['REF_ATLASLINES'], 'r') as f:
            alines = f.readlines()


        # extract information...
        atlaslines =  np.array([float(l.split()[1]) for l in alines])
        self._atlasline = atlaslines
        #self._data = pd.DataFrame({'atlas_lines': atlaslines, 'in_atlas_sni'})
        return self._atlasline

    def atlasext(self, o):
        """
        the extracted lines of atlas
        """
        minlambda = np.min(self.lam[o][self.I[o]])   ## min lambda of good lams in order
        maxlambda = np.max(self.lam[o][self.I[o]])   ## max lambda of good lams in order

        indexx, = np.where((minlambda < self.atlasline) & (self.atlasline < maxlambda))
        tmp = self.atlasline[indexx]

        #self._data['in_atlas_snippet'] \
        #    = (minlambda < self._data['atlas_lines']) & (maxlambda > self._data['atlas_lines'])

        ## use exclusions
        exclusion = np.loadtxt(self.kwargs['EXCLUSION'])
        I, = np.where(exclusion[:,0] == o)
        exc = exclusion[I]

        #self._data['excluded'] = \
        #    [np.all( np.logical_or( l < exc[:,1], l > exc[:,2]))
        #     for l in self._data['atlas_lines']]

        # goodlines = tmp
        #goodlines = []
        #for i in I:
        #    for l in tmp:
        #        if l >= exclusion[i,1] and l <= exclusion[i,2]:
        #            continue
        #        goodlines.append(l)

        ## excluding lines in exc
        goodlines = [l for l in tmp
            if np.all( np.logical_or( l < exc[:,1], l > exc[:,2]))]
        return np.array(goodlines)

    @property
    def flux(self):
        if self._flux is None:
            self._prepare()
        return self._flux

    @property
    def lam(self):
        if self._lam is None:
            self._prepare()
        return self._lam

    @property
    def bare_voie(self):
        if self._bare_voie is None:
            self._prepare()
        return self._bare_voie

    @property
    def I(self):
        if self._I is None:
            self._prepare()
        return self._I

    def _snippet(self, o):

        atlasext = self.atlasext(o)
        latlasext=atlasext*(1.-self.kwargs['VRANGE']/self.kwargs['C_LIGHT'])
        ratlasext=atlasext*(1.+self.kwargs['VRANGE']/self.kwargs['C_LIGHT'])

        numlines=atlasext.shape[0]
        #maxi = np.zeros (numlines, dtype =  float)
        #maxir = np.zeros (numlines, dtype =  float)

        #ici on selectionne des snippets autour de chaque raie du catalogue
        #de reference. Il faut que la grille en longueur d'onde soit des le
        #depart suffisamment bonne pour que la raie observee tombe dans le snippet
        lam = self.lam[o]
        flux = self.flux[o]
        bare_voie = self.bare_voie[o]
        snip = []
        prob = [] # problematic snippets
        for l,c,r in zip(latlasext,atlasext,ratlasext):
            indext,  =  np.where((lam >= l) & (lam <= r))
            wave=lam[indext]
            inte=flux[indext]
            bare_inte = bare_voie[indext]

            # selectionner que les raies du ThAr observes au dessus du seuil.
            # pour chaque raie k on determine le maximum de flux
            distmax = 1   ## TODO: make global constant
            goodsnippet = True
            # goodsnippet = goodsnippet and (np.max(inte) - np.min(inte)) >= SEUIL
            # goodsnippet = goodsnippet and (np.max(inter) - np.min(inter)) >= SEUILR
            try:
                goodsnippet = goodsnippet \
                    and (np.argmax(inte)>= distmax) \
                    and (np.argmax(inte) <= inte.shape[0]-distmax-1)
            except:
                goodsnippet = False
            sni ={
                "true_order_number": o,
                "ref_lambda": c ,
                "pixels_extract": indext,
                "pixel_mean" : np.mean(indext),
                "wave": wave,
                "reduced_flux_values_extract": inte,
                "flux_values_extract" : bare_inte,
            }
            if goodsnippet:
                snip.append(sni)
            else:
                prob.append(sni)
        return snip, prob

    @property
    def snippets(self):
        if self._snippets is None:
            tmp = [ pd.DataFrame(self._snippet(o)[0])
                for o in self.kwargs['ORDERS']
            ]
            self._snippets = pd.concat(tmp, ignore_index=True, axis=0)
        return self._snippets

if __name__ == "__main__":
    import os
    snip = Snippets(voie=1, tharfits=os.path.join(DATADIR,'NEO_20200202_173811_th0.fits'), **kwargs)
