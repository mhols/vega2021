__doc__="""
Models all observations in a night
"""

import os
import extract
import numpy as np

class Night:
    def __init__(self, dir, **kwargs):
        """
        It is assumed that all fitsfiles in the 
        directory dir belong to a single night

        :dir: (String) name of directory
        """
        self._dir = dir
        self.kwargs = kwargs

        self._sort_fitsies()
    
    def _sort_fitsies(self):
           
        res = []
        for f in os.listdir(self._dir):
            ff = os.path.join(self._dir, f)
            res.append({
                'fitsfile': ff, 
                'time': extract.gettimestamp(ff),
                'is_thorium': extract.is_thorium(ff),
                'is_star': extract.is_star(ff, 'pagasus'),
                'is_bias': extract.is_bias(ff),
                'is_flat': extract.is_flatfield(ff),
            })
        
        self._unclassified = []
        for f in res:
            if not (f['is_thorium'] or f['is_star'] or f['is_bias'] or f['is_flat']):
                self._unclassified.append(f)

        self._fitsies = sorted (res, key=lambda l: l['time'])
                   

    @property
    def thoriums(self):
        return [l for l in self._fitsies if l['is_thorium']]
    
    @property
    def stars(self):
        return [l for l in self._fitsies if l['is_star']]
     
    @property
    def flats(self):
        return [l for l in self._fitsies if l['is_flat']]
     
    @property
    def bias(self):
        return [l for l in self._fitsies if l['is_bias']]
     
    @property
    def time_interval(self):
        return self._fitsies[0]['time'], self._fitsies[-1]['time']
    
    @property
    def time_interval_stars(self):
        return self.stars[0]['time'], self.stars[-1]['time']

    @property
    def four_groups(self):
        MAXTIMELEAP = 0.02
        tt = np.array([l['time'] for l in self.stars])
        dt = tt[1:]-tt[:-1]

        groups = []
        group = [self.stars[0]]
        lt = group[-1]['time']

        for l in self.stars[1:]:
            if l['time'] - lt < MAXTIMELEAP:
                group.append(l)
            else:
                groups.append(group)
                group=[]
        
        groups.append(group)
        
        return groups
    
    def closest_tha(self, star):
        i = np.argmin([np.abs(star['time'] -t ) for t in [ff['time'] for ff in self.thoriums]])
        return self.thoriums[i]
    
    def prepare_extracts(self):
        for f in self.thoriums:
            f['extractor'] = extract.Extractor(f['fitsfile'], **self.kwargs)
            for i in range(3):
                f['extractor'].update()
            f['extractor'].save_to_store()

        for s in self.stars:
            s['extractor'] = extract.get_ext(self.closest_tha(s)['fitsfile'])
            s['extractor'].set_fitsfile(s['fitsfile'])
        


        


