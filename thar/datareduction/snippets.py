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
from units import *

class Snippets:

    def __init__(self, voie, extractor):
        self.voie = voie
        self._atlasline = None
        self._snippets = None
        self.extractor = extractor

        self.NROWS = self.extractor.NROWS
        self.ORDERS = self.extractor.ORDERS

    @property
    def kwargs(self):
        return self.extractor.kwargs


    # TODO implement specific atlas routines
    @lazyproperty
    def atlasline_uves(self):

        # print('using atlas UVES')
        
        with open(self.kwargs['REF_ATLASLINES'], 'r') as f:
            alines = f.readlines()

        # extract information...
        atlaslines =  np.array([float(l.split()[1]) for l in alines])
        # self._atlasline = atlaslines
        l = (1.-self.kwargs['VRANGE']/self.kwargs['C_LIGHT'])
        r = (1.+self.kwargs['VRANGE']/self.kwargs['C_LIGHT'])

        atlasline = pd.DataFrame(
            {'ref_lambda': pd.Series(atlaslines),
             'lower': pd.Series( l * atlaslines ),
             'upper': pd.Series( r * atlaslines)}
        )
    
        return atlasline

    @lazyproperty 
    def atlasline_redman(self):
        """
        usage of Redman catalogue 
        """

        # print('using atlas REDMAN')

        f = os.path.join(self.kwargs['REFFILES'], 'Redman_table6.dat')

        d = pd.read_fwf(f, names=[i for i in range(1,14)],infer_nrows=10000)
        # extract information...

        # atlaslines =  np.array([float(l.split()[1]) for l in alines])
        # self._atlasline = atlaslines
        l = (1.-self.kwargs['VRANGE']/self.kwargs['C_LIGHT'])
        r = (1.+self.kwargs['VRANGE']/self.kwargs['C_LIGHT'])

        tmp = pd.DataFrame(
            {'ref_lambda': d[3] / util.vac_to_air(d[3] * ANGSTROM),
             'lower':  l * d[3] / util.vac_to_air(d[3] * ANGSTROM),
             'upper':  r * d[3] / util.vac_to_air(d[3] * ANGSTROM),
             'uncertainty': d[4] / util.vac_to_air(d[3] * ANGSTROM),
             'relative_intensity': d[7]}
        )
    
        return tmp
    
    @property
    def atlasline(self):
        atlas_dict = {
            'REDMAN': self.atlasline_redman,
            'UVES': self.atlasline_uves
        }
        key = self.kwargs.get('ATLAS_FOR_SNIPPETS', 'UVES')
        
        # print('using now atlas', key)
        return atlas_dict[key]


    def lambda_range(self, o):
        return  self.extractor.pix_to_lambda_map_voie[self.voie][o](
            self.extractor.beams[o].pixel_range
        )

    def atlasext(self, o):
        """
        the extracted lines of atlas
        """
        minlambda, maxlambda = self.lambda_range(o)

        indexx = (minlambda < self.atlasline['ref_lambda']) & (self.atlasline['ref_lambda'] < maxlambda)
        tmp = self.atlasline[indexx]

        exclusion = np.loadtxt(self.kwargs['EXCLUSION'])

        I, = np.where(exclusion[:,0] == o)
        exc = exclusion[I]

        j = pd.Series(True, index=tmp.index)
        for e in exc:
            j = j & ( (tmp['ref_lambda'] < e[1]) | (tmp['ref_lambda'] > e[2]) )
        return tmp[j]


    def _potential_snippets(self, o):
        """
        collect all information that describes the snippets
        at order o
        this feature vector is later used for the matching
        with the catalog
        """

        NMAX_LINES = self.kwargs.get('MAX_LINES', 50) # maximal number of lines to extract

        # the signal used to define the snippets
        # TODO wrong name and use additional voices... singal to noise
    
        v = self.bare_voie(o)

        # hierarchical decomposition of signal into chuncs
        lm = util.local_maxima(v)

        # collect all information about the chuncs
        bs = []
        for i,xab in enumerate ( lm[:min(NMAX_LINES, len(lm))]):
            x, a, b = xab  # the position of max and the limits
            if b-a<=1:
                continue

            s = v[np.arange(a,b+1)]
            if len(s)==0:
                continue

            #try:
            #    A, mu, sigma, offset = util.estimate_location(s, **self.kwargs)
            #except Exception as ex:
            #    continue
            
            bs.append( {
                'true_order_number': o,
                'index_pixel_snippet': i, 
                'posmax': x,
                'left': a, 
                'right': b, 
                'pixel_A': 0.0, #A,
                'pixel_mu': 0.0, #mu,
                'pixel_sigma': 1.0, #sigma,
                'pixel_offset': 0.0, #offset,
                'pixel_mean': x, #a + mu,   # TODO change name to pixel_lage
                'pixel_std' : 1./ np.sqrt(np.sum(s)),  
                'pixel_sum_intens': np.sum(s),  # TODO A???
                'pixel_max_intens': v[x],
                'pixel_range' : np.arange(a, b+1),
                'bare_voie': s,      # TODO include true bare voie
                'bootstraped': False
            })
        bs = pd.DataFrame(bs)
        # map the central positions of the snippets to the lambda axis
        try:
            bs['est_lambda'] = self.extractor.pix_to_lambda_map_voie[self.voie][o](
                bs.loc[:, 'pixel_mean']
            )
        except:
            print('problem in order ', o)
            pass
        return bs

    def bare_voie(self, o):
        """
        convinience method. TODO: put into Extractor
        """
        l, v, I = self.extractor.get_lambda_intens(self.voie, o)
        
        res = np.zeros(self.NROWS)
        res[I] = v[I]
        return util.clean_nans(res)
    
    def clean_voie(self, o):
        pass

    def _match_snippets(self, o):
        """
        central part of algorithm. Here we construct tha match bewteen
        snippets and catalog lines
        """

        # the dataframe of the cleaned catalog at order o
        atlasext = self.atlasext(o)
        
        # select potential snippet lines at order o
        pot_snip = self.snippets.loc[self.snippets['true_order_number'] == o]

        # put your matchin logic here....
        # possibly use a matching table but for the moment we
        # have implemented a one way search starting from the catalog lines
        # and finding the matching snippets

        # take 50 strongest in each order

        n = len(atlasext)

        NCATAL = self.kwargs.get('NCATAL', 50)

        if NCATAL < n:
            try:
                atlasext = atlasext[atlasext.relative_intensity >= 
                        np.quantile(atlasext.relative_intensity, 1 - NCATAL / n)].copy()
            except:
                pass


        matchings = []


        # take  

        for ic in atlasext.index:

            row = atlasext.loc[ic]
            l = row['lower']
            c = row['ref_lambda']
            r = row['upper']

            I, = np.where( np.logical_and ( 
                l <= pot_snip['est_lambda'].values, pot_snip['est_lambda'].values <= r)
            )

            goodlam = len(I) == 1 ## exacly one peak in interval
            goodlam = goodlam and True ## SEUIL in extract TODO

            if not goodlam:
                continue

            # Gaussian fit for good snippets
            

            i = I[0]  # index of single good snippet belonging to l
            idx = pot_snip.index[i]

            if not self._snippets.loc[idx, 'gaussfit']:
                s = self._snippets.loc[idx,'bare_voie'] 
                try:
                    A, mu, sigma, offset = util.estimate_location(s, **self.kwargs)
                except Exception as ex:
                    continue
                self._snippets.loc[idx, 'pixel_A'] = A,
                self._snippets.loc[idx, 'pixel_mu'] = mu,
                self._snippets.loc[idx, 'pixel_sigma'] = sigma,
                self._snippets.loc[idx, 'pixel_offset'] = offset,
                self._snippets.loc[idx, 'gaussfit'] = True
                self._snippets.loc[idx, 'pixel_mean'] = self._snippets.loc[idx, 'left'] + mu
             # else, we have a matchi

            # keep track of matchings for later use 

            
            self._snippets.loc[idx,'goodsnippet'] = True
            self._snippets.loc[idx,'ref_lambda'] = c
            self._snippets.loc[idx,'catalog_index'] = int(ic)
           
            # bootstrapping
            if self.kwargs['n_bootstrap'] > 1 and not \
                self._snippets.loc[idx, 'bootstraped']:
                intens = self._snippets.loc[idx, 'bare_voie']
                mu, s, sample = util.bootstrap_estimate_location(intens, **self.kwargs)

                self._snippets.loc[idx,'pixel_mean'] = pot_snip.loc[idx, 'left']+mu
                self._snippets.loc[idx,'pixel_std'] = s
                self._snippets.loc[idx, 'bootstraped'] = True
            else:
                pass

            self._snippets.loc[idx, 'est_lambda'] = \
                self.extractor.pix_to_lambda_map_voie[self.voie][o](self._snippets.loc[idx, 'pixel_mean'])
            matchings.append([ic, idx])


    @property 
    def snippets(self):
        if not self._snippets is None:
            return self._snippets
        self.extractor.logging('computing snippets for voie '+str(self.voie))
        
        # generating all possible snippets
        tmp = []
        for i, o in enumerate(self.ORDERS):
            util.progress(i, len(self.ORDERS))
            tmp.append(self._potential_snippets(o))
        tmp = pd.concat(tmp, ignore_index=True, axis=0)

        self.extractor.end_logging()

        self._snippets = tmp
        self._snippets['goodsnippet'] = False
        self._snippets['ref_lambda'] = -1.0
        self._snippets['catalog_index'] = -1
        self._snippets['selected'] = True
        self._snippets['total_flux'] = 0.0
        self._snippets['bootstraped'] = False
        self._snippets['gaussfit'] = False

        return self.update_snippets()
        
    def filter_snippets_max_amplitude(self,o):
        alpha = self.kwargs.get("FILTER_AMPLITUDE_QUANTILE", 1.)
        # filter says true or false for each snippe
        res = pd.Series(False,index=self._snippets.index)
        II = self._snippets["true_order_number"] == o
        sn = self._snippets[II]
        A_crit = np.quantile(sn["pixel_A"],alpha)
        III = sn["pixel_A"] > A_crit
        res.loc[III.index] = True
        return res
        
    def filter_snippets_width(self,o):
        alpha = self.kwargs.get("FILTER_WIDTH_QUANTILE", 1.)
        # filter says true or false for each snippe
        res = pd.Series(False,index=self._snippets.index)
        II = self._snippets["true_order_number"] == o
        sn = self._snippets[II]
        sigma_crit = np.quantile(sn["pixel_sigma"],alpha)
        
        #vorsicht in m/s
        
        III = sn["pixel_A"] > sigma_crit
        res.loc[III.index] = True
        return res
               
    def update_snippets(self):

        self.extractor.logging('matching snippets for voie '+str(self.voie))
        self.snippets.loc[:,'goodsnippet'] = False

        for i, o in enumerate(self.ORDERS):
            util.progress(i, len(self.ORDERS))
            self._match_snippets(o)
        self.extractor.end_logging()
        return self._snippets


    @lazyproperty
    def overlapping_snippets(self):
        tmp = []
        for o, oo in zip(self.ORDERS[:-1], self.ORDERS[1:]):

            Io = self.snippets['true_order_number'] == o
            Ioo = self.snippets['true_order_number'] == oo
            pp1 = self.snippets.loc[Io, 'pixel_mean']
            pp2 = self.snippets.loc[Ioo,'pixel_mean']


            p1 = self.extractor.pix_to_lambda_map_voie[self.voie][o](pp1)
            p2 = self.extractor.pix_to_lambda_map_voie[self.voie][oo](pp2)

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
                }
            ))

        return pd.concat(tmp).reset_index(drop=True)

    ## convenience methods

    def index_good_order(self, o):
        """
        indices of good snippets of order o
        """
        return self.snippets['true_order_number']==o & self.good_snippet

   
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

