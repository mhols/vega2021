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

class Atlas:
    """
    Default Atlas class. Should be overwritten for specific catalogs (see below)

    :param file_of_lines: path to file containing list of wavelength. defaults to None
    :param **kwargs: list of options. Typically provided by some setting module. defaults to None
    """
    def __init__(self, extract, file_of_lines=None, **kwargs):
        self.kwargs = kwargs
        self.extract = extract
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
        l = (1.-self.kwargs['VRANGE']/self.kwargs['C_LIGHT'])
        return l * self.lambda_ref(VACUUM_AIR)

    def lmax(self, VACUUM_AIR='VACUUM'):
        r = (1.+self.kwargs['VRANGE']/self.kwargs['C_LIGHT'])
        return r * self.lambda_ref(VACUUM_AIR)
    
    def atlasext(self, voie, o):
        """
        the extracted lines of atlas
        """
        minlambda, maxlambda = self.extract.lambda_range_voie(voie, o)
        VACUUM_AIR = self.kwargs['WAVEMAP_IN_VACUUM_AIR']

        I = (minlambda < self.lambda_ref(VACUUM_AIR)) & (self.lambda_ref(VACUUM_AIR) < maxlambda)
        
        l = self.lambda_ref(VACUUM_AIR)

        exclusion = np.loadtxt(self.kwargs['EXCLUSION'])

        J = exclusion[:,0] == o
        exc = exclusion[J]

        #pd.Series(True, index=tmp.index)
        for e in exc:
            I = I & ( (l < e[1]) | (l > e[2]) )
        return I


class Atlas_UVES(Atlas):

    def __init__(self, extractor, **kwargs):
        super(Atlas_UVES, self).__init__(extractor, **kwargs)
        uves = pd.read_csv(kwargs['REF_ATLASLINES_UVES'],  
                 header=None, 
                     sep=r"\s+", 
                     usecols=[0,1,2,3,4,5], 
                     engine='python')

        
        self._lambda_vacuum = 1e8 / uves.loc[:, 0]
        self._lambda_air = uves.loc[:,1]


class Atlas_Redman(Atlas):
    pass

class Atlas_Clicked(Atlas):
    def __init__(self, extractor, **kwargs):
        super(Atlas_Clicked, self).__init__(extractor, **kwargs)
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


class SnippetsMixin:

    def __init__(self):

        for voie in [1,2]:
            self.logging('computing all potential snippets for voie '+str(voie))
            self._prepare_snippets(voie)
            self.end_logging()
        for voie in [1,2]:
            self._compute_est_lambda(voie)
        
    def _prepare_snippets(self, voie):
        
        # generating all possible snippets
        tmp = []
        for o in self.ORDERS:
            tmp.append(self._potential_snippets(voie, o))
        tmp =  pd.concat(tmp, ignore_index=True, axis=0)

        setattr(self, 'snippets_voie'+str(voie), tmp)


    def _potential_snippets(self, voie, o):
        """
        collect all information that describes the snippets
        at order o
        this feature vector is later used for the matching
        with the catalog
        """


        # the signal used to define the snippets
    
        v = self.bare_voie[voie][o]

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
                'left': a, 
                'right': b, 
                'pixel_A': 0.0, #A,
                'pixel_mu': x-a, #mu,
                'pixel_sigma': 1.0, #sigma,
                'pixel_offset': 0.0, #offset,
                'pixel_mean': x, #a + mu,   # TODO change name to pixel_lage
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

    @property
    def snippets_voie(self):
        return {1: self.snippets_voie1, 2: self.snippets_voie2, 3: None}

    def update_snippets(self):
        for voie in [1,2]:
            self._filter_snippets(voie)
            self._match_snippets(voie)
            self._update_snippet_parameters(voie)
            self._compute_est_lambda(voie)


    @lazyproperty
    def _atlas(self):
        atlas_dict = {
            # 'REDMAN': self.atlasline_redman,
            'UVES': Atlas_UVES(self, **self.kwargs),
            'CLICKED': Atlas_Clicked(self, **self.kwargs)
        }
        return atlas_dict

    @property
    def atlas(self):
        return self._atlas[self.kwargs['ATLAS_FOR_SNIPPETS']]


  
    def find_snippet(self, lam):
        """
        :returns: dataframe of snippets containing wavelength lam
        """
        df = self.snippets_voie[voie]
        I = (df['lambda_left'] <= lam) & (lam <= df['lambda_right'])
        return df[I]
    
   
    def _update_snippet_parameters(self, voie):
        """
        For "good" snippets that have not been treated for Gaussian fit
        compute the relevant paramters. If bootstrap is required do it now.
        """


        self.logging('updating snippet parameters')

        index = self.snippets_voie[voie][
            np.logical_and(
                np.logical_not (self.snippets_voie[voie]['gaussfit']),
                self.snippets_voie[voie]['goodsnippet']
            )].index

        for i, idx in enumerate(index):

            util.progress(i, len(index))

            s = self.snippets_voie[voie].loc[idx,'bare_voie'] 
            try:
                A, mu, sigma, offset = util.estimate_location(s, **self.kwargs)
            except Exception as ex:
                continue
            self.snippets_voie[voie].loc[idx, 'pixel_A'] = A,
            self.snippets_voie[voie].loc[idx, 'pixel_mu'] = mu,
            self.snippets_voie[voie].loc[idx, 'pixel_sigma'] = sigma,
            self.snippets_voie[voie].loc[idx, 'pixel_offset'] = offset,
            self.snippets_voie[voie].loc[idx, 'gaussfit'] = True
            self.snippets_voie[voie].loc[idx, 'pixel_mean'] = self.snippets_voie[voie].loc[idx, 'left'] + mu


        self.end_logging()
        
        if self.kwargs['n_bootstrap'] > 1:
            
            self.logging('bootstrapping ')
            index =  self.snippets_voie[voie][
                np.logical_and(
                    np.logical_not(self.snippets_voie[voie]['bootstraped']),
                    self.snippets_voie[voie]['goodsnippet']
                )
            ].index
            
            for i, idx in enumerate(index):
                util.progress(i, len(index))
                
                intens = self.snippets_voie[voie].loc[idx, 'bare_voie']
                mu, s, sample = util.bootstrap_estimate_location(intens, **self.kwargs)

                self.snippets_voie[voie].loc[idx, 'pixel_std'] = s
                self.snippets_voie[voie].loc[idx, 'bootstraped'] = True

            self.end_logging()        

        
    def _compute_est_lambda(self, voie):
        """
        computing the estimated wavelength of the snippet bounds using the
        existing wave map
        """
        sn = self.snippets_voie[voie] 
        # getattr(self, "snippets_voie"+str(voie))
        for o in self.ORDERS:
            p = self.pix_to_lambda_map_voie[voie][o]
            I = sn['true_order_number'] == o 
            sn.loc[I, 'est_lambda'] = \
                p( sn.loc[I, 'pixel_mean'].values )
            sn.loc[I, 'lambda_left'] = \
                p( sn.loc[I, 'left'].values )
            sn.loc[I, 'lambda_right'] = \
                p( sn.loc[I, 'right'].values )
        

    #---------------------
    # FILTER methods for snippets
    #---------------------

    def _filter_snippets_not_excluded(self, voie, o):
        """
        return index of snippets that are not excluded because of
        masking Argon lines
        """

        exclusion = np.loadtxt(self.kwargs['EXCLUSION'])
        I, = np.where(exclusion[:,0] == o)
        exc = exclusion[I]

        res = self.snippets_voie[voie]["true_order_number"] == o
        sn = self.snippets_voie[voie]

        for ex in exc:
           res = res  &  ( (sn['est_lambda'] < ex[1]) | (sn['est_lambda'] > ex[2]) )

        return res
        
        
    def _filter_snippets_min_length(self,voie,o):
        alpha = self.kwargs.get("FILTER_MIN_LENGTH", 3)
        # filter says true or false for each snippe
        res = self.snippets_voie[voie]["true_order_number"] == o
        res = res & ((self.snippets_voie[voie]["right"] - self.snippets_voie[voie]["left"]) > alpha)
        return res
        
    def _filter_snippets_max_amplitude(self,voie, o):
        alpha = self.kwargs.get("FILTER_AMPLITUDE_QUANTILE", 0.1)

        res = self.snippets_voie[voie]["true_order_number"] == o
        q = np.quantile(self.snippets_voie[voie].loc[res, "pixel_max_intens"], alpha)
        res = res & (self.snippets_voie[voie]["pixel_max_intens"] > q)

        return res
    
       
    def _filter_snippets(self, voie):
        II = False
        for o in self.ORDERS:
            I = True
            I = I & self._filter_snippets_not_excluded(voie, o)
            I = I & self._filter_snippets_max_amplitude(voie, o)
            I = I & self._filter_snippets_min_length(voie, o)
            II = I | II 
        self.snippets_voie[voie].loc[II,'usablesnippet'] = True
        return II
    
    def _match_snippets(self, voie):
        VACUUM_AIR = self.kwargs['WAVEMAP_IN_VACUUM_AIR']
        self.snippets_voie[voie].loc[:, 'goodsnippet'] = False
        ls = self.snippets_voie[voie]['est_lambda'].to_numpy(copy=True)
        lsl = self.snippets_voie[voie]['lambda_left'].to_numpy(copy=True)
        lsr = self.snippets_voie[voie]['lambda_right'].to_numpy(copy=True)

        lc = self.atlas.lambda_ref(VACUUM_AIR)
        lcl = self.atlas.lmin(VACUUM_AIR)
        lcr = self.atlas.lmax(VACUUM_AIR)

        for o in self.ORDERS:
            I = self.snippets_voie[voie]['true_order_number'] == o
            I = I & self.snippets_voie[voie]['usablesnippet']
            InI = self.snippets_voie[voie].index[I]

            J = self.atlas.atlasext(voie, o)
            # possible snippets
            #M1 = (lcl <= ls[I]) & (ls[I] <= lcr)
            #M2 = (lsl[I] <= lc)

            II, JJ = util.match_intervals_centers (
                ls[I],
                lsl[I],
                lsr[I],
                lc[J],
                lcl[J],
                lcr[J]
            )

            self.snippets_voie[voie].loc[InI[II], 'goodsnippet'] = True
            self.snippets_voie[voie].loc[InI[II], 'ref_lambda'] = lc[J][JJ]

              
    def update_snippets_voie(self, voie):

        # estimate lamda for snippet
        # self.compute_est_lambda()

        self.logging('matching snippets for voie '+str(self.voie))
        
        # nobody is good (memoryless matching) 
        self.snippets_voie[voie].loc[:,'goodnsippets'] = False


        for i, o in enumerate(self.ORDERS):
            util.progress(i, len(self.ORDERS)-1)
            self._match_snippets_voie[voie](o)

            self.filter_matched_snippets_max_number(o)
        
        self.end_logging()

        # compute the Gaussfit and bootstrap for the used snippets
        self._update_snippet_parameters()

        # recompute the estimated lambda based on the present wavemap
        # estimate lamda for snippet
        self.compute_est_lambda()


        return self.snippets_voie[voie]


    @lazyproperty
    def overlapping_snippets(self):
        tmp = []
        for o, oo in zip(self.ORDERS[:-1], self.ORDERS[1:]):

            Io = self.snippets_voie[voie]['true_order_number'] == o
            Ioo = self.snippets_voie[voie]['true_order_number'] == oo
            pp1 = self.snippets_voie[voie].loc[Io, 'pixel_mean']
            pp2 = self.snippets_voie[voie].loc[Ioo,'pixel_mean']


            p1 = self.pix_to_lambda_map_voie[self.voie][o](pp1)
            p2 = self.pix_to_lambda_map_voie[self.voie][oo](pp2)

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
        return self.snippets_voie[voie]['true_order_number']==o & self.good_snippet

   
    @property
    def good_snippet(self):
        return self.snippets_voie[voie]['goodsnippet'] == True
    
    @property
    def _x(self):
        return self.snippets_voie[voie]['pixel_mean']
    
    @property
    def x(self):
        return self._x[self.good_snippet]
    
    @property
    def _o(self):
        return self.snippets_voie[voie]['true_order_number']
    
    @property
    def o(self):
        return self._o[self.good_snippet]

    @property
    def _l(self):
        return self.snippets_voie[voie]['ref_lambda']
    
    @property
    def l(self):
        return self._l[self.good_snippet]
      
    @property
    def _sigma(self):
        return self.snippets_voie[voie]['pixel_std']
    
    @property
    def sigma(self):
        return self._sigma[self.good_snippet]

    @property
    def _sn(self):
        return self.snippets_voie[voie]

    @property
    def sn(self):
        return self.snippets_voie[voie][self.good_snippet]
    
    @property
    def ci(self):
        return self.snippets_voie[voie][self.good_snippet]['catalog_index']    


    @property
    def ngood(self):
        return sum(self.good_snippet)

