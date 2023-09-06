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

        self.prepare_snippets()

    @property
    def kwargs(self):
        return self.extractor.kwargs


    def prepare_snippets(self):
        if not self._snippets is None:
            return self._snippets
        self.extractor.logging('computing all potential snippets for voie '+str(self.voie))
        
        # generating all possible snippets
        tmp = []
        for i, o in enumerate(self.ORDERS):
            tmp.append(self._potential_snippets(o))
        self._snippets = pd.concat(tmp, ignore_index=True, axis=0)

        self.extractor.end_logging()


        self._snippets.loc[:,'usablesnippet'] = False
        I = self.filter_snippets()
        self._snippets.loc[I, 'usablesnippet'] = True


    # TODO implement specific atlas routines
    @lazyproperty
    def atlasline_uves(self):

        # print('using atlas UVES')
        
        with open(self.kwargs['REF_ATLASLINES_UVES'], 'r') as f:
            alines = f.readlines()

        # extract information...
        if self.kwargs['WAVEMAP_IN_VACUUM_AIR'] == 'VACUUM':
            atlaslines =  np.array([1e8 / float(l.split()[0]) for l in alines])
        else:
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

        f = self.kwargs['REF_ATLASLINES_REEDMAN']

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
    def atlasline_clicked(self):
        tmp = pd.read_csv(self.kwargs['REF_ATLASLINES_CLICKED'])

        uves = pd.read_csv(self.kwargs['REF_ATLASLINES_UVES'],  
                    header=None, 
                    sep=r"\s+", 
                    usecols=[0,1,2,3,4,5], 
                    engine='python')

        if self.kwargs['WAVEMAP_IN_VACUUM_AIR'] == 'VACUUM':
            tmp.loc[:,'ref_lambda'] = 1e8 / uves.loc[:, 0]
        else:
            tmp.loc[:,'ref_lambda'] = uves.loc[:, 1]
        
        I = tmp['click_selected']

        return tmp.loc[I, :]


    @property
    def atlasline(self):
        atlas_dict = {
            # 'REDMAN': self.atlasline_redman,
            'UVES': self.atlasline_uves,
            'CLICKED': self.atlasline_clicked
        }
        key = self.kwargs.get('ATLAS_FOR_SNIPPETS', 'UVES')
        
        # print('using now atlas', key)
        return atlas_dict[key]


    def lambda_range(self, o):
        return  self.extractor.lambda_range_voie(self.voie, o)

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


  
    def find_snippet(self, lam):
        df = self._snippets
        I = (df['lambda_left'] <= lam) & (lam <= df['lambda_right'])
        return df[I]
    
    def _potential_snippets(self, o):
        """
        collect all information that describes the snippets
        at order o
        this feature vector is later used for the matching
        with the catalog
        """


        # the signal used to define the snippets
        # TODO wrong name and use additional voices... singal to noise
    
        v = self.bare_voie(o)

        # hierarchical decomposition of signal into chuncs
        # the snippets are in decreasing order with respect to value of maximum
        lm = util.local_maxima(v)


        # collect all information about the chuncs
        bs = []


        for i,xab in enumerate (lm):
            
            x, a, b = xab  # the position of max and the limits
            
                   
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
                'pixel_std' : 1./ np.sqrt(np.sum(s)),  
                'pixel_sum_intens': np.sum(s),  # TODO A???
                'pixel_max_intens': v[x]-np.min(v),
                'pixel_range' : np.arange(a, b+1),
                'bare_voie': s,      # TODO include true bare voie
                'bootstraped': False,
                'selected': False,
                'goodsnippet': False,
                'usablesnippet': False,
                'gaussfit': False
            })
        bs = pd.DataFrame(bs)
        return bs

    def bare_voie(self, o):
        """
        convinience method. TODO: put into Extractor
        """
        l, v, I = self.extractor.get_lambda_intens(self.voie, o)
        
        res = np.zeros(self.NROWS)
        res[I] = v[I]
        return util.clean_nans(res)
    
    @lazyproperty
    def bare_voies(self):
        return {o: self.bare_voie(o) for o in self.ORDERS}
    
    def _match_snippets_old(self, o):
        """
        central part of algorithm. Here we construct tha match bewteen
        snippets and catalog lines
        """

        # the dataframe of the cleaned catalog at order o
        atlasext = self.atlasext(o)
        
        # select potential snippet lines at order o
        pot_snip = self._snippets.loc[self._snippets['true_order_number'] == o]


        # put your matchin logic here....
        # possibly use a matching table but for the moment we
        # have implemented a one way search starting from the catalog lines
        # and finding the matching snippets

        # take 50 strongest in each order

        n = len(atlasext)

        NCATAL = self.kwargs.get('NCATAL', 500000)

        if NCATAL < n:
            try:
                atlasext = atlasext[atlasext.relative_intensity >= 
                        np.quantile(atlasext.relative_intensity, 1 - NCATAL / n)].copy()
            except:
                pass


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

            self._snippets.loc[idx,'goodsnippet'] = True
            self._snippets.loc[idx,'ref_lambda'] = c
            self._snippets.loc[idx,'catalog_index'] = int(ic)
        

    def update_snippet_parameters(self):
        """
        For "good" snippets that have not been treated for Gaussian fit
        compute the relevant paramters. If bootstrap is required do it now.
        """


        self.extractor.logging('updating snippet parameters')

        index = self._snippets[
            np.logical_and(
                np.logical_not (self._snippets['gaussfit']),
                self._snippets['goodsnippet']
            )].index

        for i, idx in enumerate(index):

            util.progress(i, len(index))

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


        self.extractor.end_logging()
        
        if self.kwargs['n_bootstrap'] > 1:
            
            self.extractor.logging('bootstrapping ')
            index =  self._snippets[
                np.logical_and(
                    np.logical_not(self._snippets['bootstraped']),
                    self._snippets['goodsnippet']
                )
            ].index
            
            for i, idx in enumerate(index):
                util.progress(i, len(index))
                
                intens = self._snippets.loc[idx, 'bare_voie']
                mu, s, sample = util.bootstrap_estimate_location(intens, **self.kwargs)

                self._snippets.loc[idx, 'pixel_std'] = s
                self._snippets.loc[idx, 'bootstraped'] = True

            self.extractor.end_logging()        

        
    def compute_est_lambda(self):
        """
        computing the estimated wavelength of the snippet bounds using the
        existing wave map
        """
        for o in self.ORDERS:
            p = self.extractor.pix_to_lambda_map_voie[self.voie][o]
            I = self._snippets['true_order_number'] == o 
            self._snippets.loc[I, 'est_lambda'] = \
                p( self._snippets.loc[I, 'pixel_mean'].values )
            self._snippets.loc[I, 'lambda_left'] = \
                p( self._snippets.loc[I, 'left'].values )
            self._snippets.loc[I, 'lambda_right'] = \
                p( self._snippets.loc[I, 'right'].values )
        

    #---------------------
    # FILTER methods for snippets
    #---------------------

    def filter_snippets_not_excluded(self, o):
        """
        return index of snippets that are not excluded because of
        @Torsten: please explain
        """

        exclusion = np.loadtxt(self.kwargs['EXCLUSION'])

        I, = np.where(exclusion[:,0] == o)
        
        exc = exclusion[I]

        res = self._snippets["true_order_number"] == o
        sn = self._snippets

        for ex in exc:
           res = res  &  ( (sn['est_lambda'] < ex[1]) | (sn['est_lambda'] > ex[2]) )

        return res
        
        
    def filter_snippets_min_length(self,o):
        alpha = self.kwargs.get("FILTER_MIN_LENGTH", 1.)
        # filter says true or false for each snippe
        res = pd.Series(False,index=self._snippets.index)
        res = self._snippets["true_order_number"] == o
        res = res & (self._snippets["left"] - self._snippets["right"]) > alpha
        return res
        
    def filter_snippets_max_amplitude(self,o):
        alpha = self.kwargs.get("FILTER_AMPLITUDE_QUANTILE", 0.8)

        res = self._snippets["true_order_number"] == o
        q = np.quantile(self._snippets.loc[res, "pixel_max_intens"], alpha)
        res = res & self._snippets["pixel_max_intens"] > q

        return res
    
       
    def filter_snippets(self):
        I = True
        for o in self.ORDERS:
            I = I & self.filter_snippets_not_excluded(o)
            I = I & self.filter_snippets_max_amplitude(o)
            I = I & self.filter_snippets_min_length(o)
        return I
    
    def match_snippets(self, o):
        return self._snippet_match_tor(o)

    def _snippet_match_tor(self,o):
        #selection of uves atlas lines in order o
        #cu=self.atlasline_uves
        
        cu = self.atlasext(o)
        lamlimits=self.extractor.lambda_range_voie(self.voie, o)
       
        I=(cu["ref_lambda"] > lamlimits[0]) & (cu["ref_lambda"] < lamlimits[1])
        cu=cu[I]
        
        #selection of snippets in order o, after filter with intensity
        I = self._snippets['usablesnippet'] #self.filter_snippets_max_amplitude(o) & self.filter_snippets_not_excluded(o)

        # selection=self._snippets.loc[I]
    
        vrange = self.kwargs["VRANGE"]
        
        for i in cu.index:
            c = cu.loc[i,"ref_lambda"]
            l = c * (1 - vrange/C_LIGHT)
            r = c * (1 + vrange/C_LIGHT)
            J = (l <= self._snippets["est_lambda"]) & (self._snippets["est_lambda"] <=r)
            J = I & J
            if sum( J ) != 1:  # check if only one snippet is OK
                continue

            ind=J[J].index[0]
            self._snippets.loc[ind,"goodsnippet"] = True
            self._snippets.loc[ind,'ref_lambda'] = c
            self._snippets.loc[ind,'catalog_index'] = int(i)        

    def filter_matched_snippets_max_number(self, o):
        alpha = self.kwargs.get("FILTER_MAX_NUMBER", 20)

               
    def update_snippets(self):

        # estimate lamda for snippet
        self.compute_est_lambda()

        self.extractor.logging('matching snippets for voie '+str(self.voie))
        
        # nobody is good (memoryless matching) 
        self._snippets.loc[:,'goodnsippets'] = False


        for i, o in enumerate(self.ORDERS):
            util.progress(i, len(self.ORDERS)-1)
            self.match_snippets(o)
            self.filter_matched_snippets_max_number(o)
        
        self.extractor.end_logging()

        # compute the Gaussfit and bootstrap for the used snippets
        self.update_snippet_parameters()

        # recompute the estimated lambda based on the present wavemap
        # estimate lamda for snippet
        self.compute_est_lambda()


        return self._snippets


    @lazyproperty
    def overlapping_snippets(self):
        tmp = []
        for o, oo in zip(self.ORDERS[:-1], self.ORDERS[1:]):

            Io = self._snippets['true_order_number'] == o
            Ioo = self._snippets['true_order_number'] == oo
            pp1 = self._snippets.loc[Io, 'pixel_mean']
            pp2 = self._snippets.loc[Ioo,'pixel_mean']


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
        return self._snippets['true_order_number']==o & self.good_snippet

   
    @property
    def good_snippet(self):
        return self._snippets['goodsnippet'] == True
    
    @property
    def _x(self):
        return self._snippets['pixel_mean']
    
    @property
    def x(self):
        return self._x[self.good_snippet]
    
    @property
    def _o(self):
        return self._snippets['true_order_number']
    
    @property
    def o(self):
        return self._o[self.good_snippet]

    @property
    def _l(self):
        return self._snippets['ref_lambda']
    
    @property
    def l(self):
        return self._l[self.good_snippet]
      
    @property
    def _sigma(self):
        return self._snippets['pixel_std']
    
    @property
    def sigma(self):
        return self._sigma[self.good_snippet]

    @property
    def _sn(self):
        return self._snippets

    @property
    def sn(self):
        return self._snippets[self.good_snippet]
    
    @property
    def ci(self):
        return self._snippets[self.good_snippet]['catalog_index']    


    @property
    def ngood(self):
        return sum(self.good_snippet)

