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
import pickle
import sys
import os
import json
import scipy

import util
import units

from scipy.optimize import curve_fit
import nextra.extract as extract
from nextra.units import *

class Atlas:
    """
    Default Atlas class. Should be overwritten for specific catalogs (see below)

    :param file_of_lines: path to file containing list of wavelength. defaults to None
    :param **kwargs: list of options. Typically provided by some setting module. defaults to None
    """
    def __init__(self, snippets, file_of_lines=None, **kwargs):
        self.kwargs = kwargs
        self.snippets = snippets
        self.extractor = snippets.extractor
        if not file_of_lines is None:
            self._lambda_vacuum = np.loadtxt(file_of_lines)
        self._lambda_air = None
       
    def lambda_ref(self, VACUUM_AIR='VACUUM'):
        if VACUUM_AIR == 'VACUUM':
            return self._lambda_vacuum
        elif VACUUM_AIR == 'AIR':
            return self._lambda_air
        else:
            raise Exception('wrong choice '+VACUUM_AIR)

    def lmin(self, VACUUM_AIR='VACUUM'):
        """
        default implementation for catalog uncertainty range

        l = 1.-self.kwargs['VRANGE']/units.C_LIGHT
        return l * self.lambda_ref(VACUUM_AIR)

    def lmax(self, VACUUM_AIR='VACUUM'):
        r = 1.+self.kwargs['VRANGE']/units.C_LIGHT

        return r * self.lambda_ref(VACUUM_AIR)
    
    def atlasext(self, o):
        """
        the extracted lines of atlas
        """
        minlambda, maxlambda = self.extractor.lambda_range_voie(self.snippets.voie, o)
        VACUUM_AIR = self.kwargs['WAVEMAP_IN_VACUUM_AIR']

        I = (minlambda < self.lambda_ref(VACUUM_AIR)) & (self.lambda_ref(VACUUM_AIR) < maxlambda)
        
        l = self.lambda_ref(VACUUM_AIR)

        exclusion = np.loadtxt(self.kwargs['EXCLUSION'])

        J = exclusion[:,0] == o
        exc = exclusion[J]

        #pd.Series(True, index=tmp.index)

        # TODO put excluded vacuum in file....
        for e in exc:
            if VACUUM_AIR == 'AIR':
                I = I & ( (l < e[1]) | (l > e[2]) )
            else:
                I = I & ( (l < e[1]*util.air_to_vac(e[1])) | (l > e[2]*util.air_to_vac(e[2])) )
        return I


class Atlas_UVES(Atlas):

    def __init__(self, snippets, **kwargs):
        super(Atlas_UVES, self).__init__(snippets, **kwargs)
        uves = pd.read_csv(kwargs['REF_ATLASLINES_UVES'],  
                 header=None, 
                     sep=r"\s+", 
                     usecols=[0,1,2,3,4,5], 
                     engine='python')

        
        self._lambda_vacuum = 1e8 / uves.loc[:, 0].to_numpy()
        self._lambda_air = uves.loc[:,1].to_numpy()
        
class Atlas_IVAN(Atlas):
    def __init__(self, snippets, **kwargs):
        super(Atlas_IVAN, self).__init__(snippets, **kwargs)
        ivan = open(kwargs['REF_ATLASLINES_IVAN'],'rb')
        c = pickle.load(ivan)
        wave = np.array(c['Selected lines'])
        ivan.close()
        
        self.kwargs = kwargs
        self._lambda_air = wave
        self._lambda_vacuum = util.air_to_vac(wave)*wave
        

class Atlas_Redman(Atlas):
    pass

class Atlas_Clicked(Atlas):
    def __init__(self, snippets, **kwargs):
        super(Atlas_Clicked, self).__init__(snippets, **kwargs)
        clicked= pd.read_csv(
            kwargs['REF_ATLASLINES_CLICKED'],  
                header=0, 
                sep=r",", 
                engine='python'
            )
        I = clicked['click_selected']
        
        uves = pd.read_csv(kwargs['REF_ATLASLINES_UVES'],  
                 header=None, 
                     sep=r"\s+", 
                     usecols=[0,1,2,3,4,5], 
                     engine='python')

        lambda_vacuum = 1e8 / uves.loc[:, 0]
        lambda_air = uves.loc[:,1]

 
        self._lambda_vacuum = lambda_vacuum[I].values
        self._lambda_air = lambda_air[I].values


class Snippets:

    def __init__(self, voie, extractor):
        extractor.logging(f"initializing snippets for voie {voie}")
        
        self.voie = voie
        self.extractor = extractor

        self.NROWS = self.extractor.NROWS
        self.ORDERS = self.extractor.ORDERS
        self.kwargs = extractor.kwargs

        self.snippets = self._prepare_snippets()
        # put a link to the snippets on the extractor (Mixin like behavior ...)
        setattr(self.extractor, "snippets_voie" + str(voie), self.snippets)

        # produce a first guess of snippet's lambda based on local maxima
        self._compute_est_lambda()

        # filter, match, and compute gaussfits eventually bootstrap for matching snippets
        self.update_snippets()
        self.snippets.loc[:,  'selected'] = self.snippets.loc[:, 'goodsnippet']


        extractor.end_logging()

        

    def _prepare_snippets(self):
        
        # generating all possible snippets
        tmp = []
        for o in self.ORDERS:
            tmp.append(self._potential_snippets(o))
        tmp =  pd.concat(tmp, ignore_index=True, axis=0)

        return tmp

    def _potential_snippets(self, o):
        """
        collect all information that describes the snippets
        at order o
        this feature vector is later used for the matching
        with the catalog
        """


        # the signal used to define the snippets
    
        v = self.bare_voie(o)

        # hierarchical decomposition of signal into chuncs
        # the snippets are in decreasing order with respect to value of maximum
        lm = util.local_maxima(v)

        # collect all information about the chuncs
        bs = []


        for i,xab in enumerate (lm):
            
            x, a, b = xab  # the position of max and the limits
            s = v[a:b+1]
            if (b-a<=2):
                continue
            if (np.sum(s) <= 0):
                continue
            if (np.any(s<=0)):
                continue            
                   
            bs.append( {
                'true_order_number': o,
                'index_pixel_snippet': i, 
                'posmax': x,
                'a': a,
                'b': b,
                'left': float(x)-self.kwargs['SNIPPETS_PIXEL_DELTA'],    # match interval [left, right] 
                'right': float(x)+self.kwargs['SNIPPETS_PIXEL_DELTA'],
                'pixel_A': 0.0, #A,
                'pixel_mu': float(x-a), #mu,
                'pixel_sigma': 1.0, #sigma,
                'pixel_offset': 0.0, #offset,
                'pixel_mean': float(x), 
                'pixel_std' : 1./ np.sqrt(np.sum(np.abs(s))),  
                'pixel_sum_intens': np.sum(s),  # TODO A???
                'pixel_max_intens': v[x]-np.min(v),
                'total_flux': np.sum(s - np.min(s)),
                'pixel_range' : np.arange(a, b+1),
                'bare_voie': s,      # TODO include true bare voie
                'bootstraped': False,
                'selected': False,
                'goodsnippet': False,
                'usablesnippet': False,
                'gaussfit': False,
                'ref_lambda': 0.0,
                'lambda_left': 0.0,
                'lambda_right': 0.0,
            })
        bs = pd.DataFrame(bs)

        return bs

    def update_snippets(self):
        self._filter_snippets()
        self._match_snippets()
        self._update_snippet_parameters()
        self._compute_est_lambda()


    @property
    def _atlas(self):
        atlas_dict = {
            # 'REDMAN': self.atlasline_redman,
            'UVES': Atlas_UVES(self, **self.kwargs),
            'CLICKED': Atlas_Clicked(self, **self.kwargs),
            'IVAN': Atlas_IVAN(self, **self.kwargs)
        }
        return atlas_dict

    @lazyproperty
    def atlas(self):
        return self._atlas[self.kwargs['ATLAS_FOR_SNIPPETS']]

    def lambda_range(self, o):
        return  self.extractor.lambda_range_voie(self.voie, o)

  
    def find_snippet(self, lam):
        """
        :returns: dataframe of snippets containing wavelength lam
        """
        df = self.snippets
        I = (df['lambda_left'] <= lam) & (lam <= df['lambda_right'])
        return df[I]
    
   
    def bare_voie(self, o):
        """
        convinience method. TODO: put into Extractor
        """
        return self.extractor.bare_voie[self.voie][o]

    
    def _update_snippet_parameters(self):
        """
        For "good" snippets that have not been treated for Gaussian fit
        compute the relevant paramters. If bootstrap is required do it now.
        """


        self.extractor.logging('updating snippet parameters for voie '+ str(self.voie))

        index = self.snippets[
            np.logical_and(
                np.logical_not (self.snippets['gaussfit']),
                self.snippets['goodsnippet']
            )].index

        for i, idx in enumerate(index):

            util.progress(i, len(index))

            s = self.snippets.loc[idx,'bare_voie'] 
            try:
                A, mu, sigma, offset = util.estimate_location(s, **self.kwargs)
            except Exception as ex:
                self.snippets.loc[idx, 'goodsnippet'] = False

            self.snippets.loc[idx, 'pixel_A'] = A,
            self.snippets.loc[idx, 'pixel_mu'] = mu,
            self.snippets.loc[idx, 'pixel_sigma'] = sigma,
            self.snippets.loc[idx, 'pixel_offset'] = offset,
            self.snippets.loc[idx, 'gaussfit'] = True
            self.snippets.loc[idx, 'pixel_mean'] = self.snippets.loc[idx, 'a'] + mu
            mean = self.snippets.loc[idx, 'pixel_mean']
            self.snippets.loc[idx, 'left'] =  mean - self.kwargs['SNIPPETS_PIXEL_DELTA']
            self.snippets.loc[idx, 'right'] =  mean + self.kwargs['SNIPPETS_PIXEL_DELTA']
 
        self.extractor.end_logging()
        
        if self.kwargs['n_bootstrap'] > 1:
            
            self.extractor.logging('bootstrapping ')
            index =  self.snippets[
                np.logical_and(
                    np.logical_not(self.snippets['bootstraped']),
                    self.snippets['goodsnippet']
                )
            ].index
            
            for i, idx in enumerate(index):
                util.progress(i, len(index))
                
                intens = self.snippets.loc[idx, 'bare_voie']
                mu, s, sample = util.bootstrap_estimate_location(intens, **self.kwargs)

                self.snippets.loc[idx, 'pixel_std'] = s
                self.snippets.loc[idx, 'bootstraped'] = True

            self.extractor.end_logging()        

        
    def _compute_est_lambda(self):
        """
        computing the estimated wavelength of the snippet bounds using the
        existing wave map
        """
        sn = self.snippets
        # getattr(self, "snippets_voie"+str(voie))
        for o in self.ORDERS:
            p = self.extractor.pix_to_lambda_map_voie[self.voie][o]
            I = sn['true_order_number'] == o
            if (sum(I)) == 0:
                print('no snippets in order', o)
                continue
            sn.loc[I, 'est_lambda'] = p( sn.loc[I, 'pixel_mean'].values )
            sn.loc[I, 'lambda_left'] = p( sn.loc[I, 'left'].values )
            sn.loc[I, 'lambda_right'] = p( sn.loc[I, 'right'].values )
     
    #---------------------
    # FILTER methods for snippets
    #---------------------

    def _filter_snippets_not_excluded(self, o):
        """
        return index of snippets that are not excluded because of
        masking Argon lines
        """

        exclusion = np.loadtxt(self.kwargs['EXCLUSION'])
        I, = np.where(exclusion[:,0] == o)
        exc = exclusion[I]

        res = self.snippets["true_order_number"] == o
        sn = self.snippets

        for ex in exc:
           res = res  &  ( (sn['est_lambda'] < ex[1]) | (sn['est_lambda'] > ex[2]) )
        
        if res.sum() == 0:
            raise Exception(f'no snippets found in method _not_excluded in order {o}')
        return res
        
        
    def _filter_snippets_min_length(self, o):
        alpha = self.kwargs.get("FILTER_MIN_LENGTH", 4)
        # filter says true or false for each snippe
        res = self.snippets["true_order_number"] == o
        res = res & ((self.snippets["b"] - self.snippets["a"]) > alpha)
        if res.sum() == 0:
            raise Exception(f'no snippets found in method min_length in order {o}')
        return res
        
    def _filter_snippets_max_amplitude(self, o):
        alpha = self.kwargs.get("FILTER_AMPLITUDE_QUANTILE", 0.1)

        res = self.snippets["true_order_number"] == o
        q = np.quantile(self.snippets.loc[res, "pixel_max_intens"], alpha)
        res = res & (self.snippets["pixel_max_intens"] > q)
        if res.sum() == 0:
            raise Exception(f'no snippets found in method _max_amplitude in order {o}')
 
        return res
    
       
    def _filter_snippets(self):
        II = False
        for o in self.ORDERS:
            I = True
            I = I & self._filter_snippets_not_excluded(o)
            I = I & self._filter_snippets_max_amplitude(o)
            I = I & self._filter_snippets_min_length(o)
            II = I | II 
        self.snippets.loc[II,'usablesnippet'] = True
        return II
    
    def _match_snippets(self):
        VACUUM_AIR = self.kwargs['WAVEMAP_IN_VACUUM_AIR']
        self.snippets.loc[:, 'goodsnippet'] = False
        ls = self.snippets['est_lambda'].to_numpy()
        lsl = self.snippets['lambda_left'].to_numpy()
        lsr = self.snippets['lambda_right'].to_numpy()

        lc = self.atlas.lambda_ref(VACUUM_AIR)
        lcl = self.atlas.lmin(VACUUM_AIR)
        lcr = self.atlas.lmax(VACUUM_AIR)

        for o in self.ORDERS:
            I = self.snippets['true_order_number'] == o
            I = I & self.snippets['usablesnippet']
            InI = self.snippets.index[I]

            
            J = self.atlas.atlasext(o)
            JnJ, = np.where(J)

            II, JJ = util.match_intervals_centers (
                ls[I],
                lsl[I],
                lsr[I],
                lc[J],
                lcl[J],
                lcr[J]
            )

            self.snippets.loc[InI[II], 'goodsnippet'] = True
            self.snippets.loc[InI[II], 'ref_lambda'] = lc[JnJ[JJ]]

       
