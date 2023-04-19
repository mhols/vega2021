import matplotlib.pyplot as plt
from matplotlib import cm
import util
import numpy as np

class PlotExtractMixin:
    def __init__(self):
        pass

    #---------PLOT methods----------------#
    def color_1(self, voie, o):
        cmap = cm.get_cmap(self.kwargs.get('palette_order', 'rgb'))
        orders = self.ORDERS
        colors = cmap(np.linspace(0,1, max(orders) - min(orders) +1))
        return colors[o-min(orders)]

    def color_2(self, voie, o):
        if voie==1:
            return 'r' if o%2==0 else 'b'
        if voie==2:
            return 'k' if o%2==0 else 'y'

    def plot_voie(self, voie, oo, **kwargs):
        if type(oo) is int:
            oo = list(oo)
        for o in oo:
            l =  self.lambdas_per_order_voie[voie][o]
            I = self.I[o]
            v = self.voie[voie][o]
            plt.plot(
               l[I], v[I], '-', color=self.color_2(voie, o)
        )
            
    def plot_snippets(self, voie, oo):
        if type(oo) is int:
            oo = list(oo)
        
        for o in oo:
            s = self.snippets_voie()[voie].sn
            I = s['true_order_number'] == o
            for i in s[I].index:
                pix = s.loc[i, 'pixel_range']
                lams = self.pix_to_lambda_map_voie[voie][o](pix)
                plt.fill_between(lams, self.voie[voie][o][pix], color='gray') 
 

    def plot_catalog(self, voie, oo):
        if type(oo) is int:
            oo = list(oo)
        I = self.I[33]
        v = util.clean_nans(self.voie1[33][I])
        m = np.max(v)
        dashes = np.linspace(0, 0.5 * m, 50)
        for o in oo:
            s = self.snippets_voie()[voie].sn
            I = s['true_order_number'] == o
            d = dashes if o%2 == 0 else dashes[1:]
            for a, b in zip(d[:-1:2], d[1::2]):
                plt.vlines(
                    s[I]['ref_lambda'], a, b,
                    color=self.color_2(voie, o))
        sn = self.snippets_voie()[voie].snippets

        for o in (ooo for ooo in oo if ooo%2==0):
            s = self.snippets_voie()[voie].atlasext(o)
            l =s['ref_lambda']
            if len(l)>0:
                plt.plot(l, len(l)*[0], 'o', color=self.color_2(voie, o))

        for o in (ooo for ooo in oo if ooo%2==1):
            s = self.snippets_voie()[voie].atlasext(o)
            l =s['ref_lambda']
            if len(l)>0:
                plt.plot(l, len(l)*[0], '.', color=self.color_2(voie, o))

    def plot_catalog_redman(self, oo):
        if type(oo) is int:
            oo = list(oo)
        for o in oo:
            I = self.I[33]
            v =util.clean_nans(self.voie1[33][I])
            m = np.max(v)
            cr=self.snippets_voie1.atlasline_redman
            lamlimits=self.lambda_range_voie1(o)
            I=(cr["ref_lambda"] > lamlimits[0]) & (cr["ref_lambda"] < lamlimits[1])
            cr=cr[I]
            cm = np.max(cr["relative_intensity"])
            plt.vlines(cr["ref_lambda"],0,cr["relative_intensity"]*0.5*m/cm,'k')
        

    def plot_polynomial(self, voie, oo):
        if type(oo) is int:
            oo = list(oo)

        for o in oo:
            I = self.ccd_voie[voie].index_order(o)
            plt.plot(self.ccd_voie[voie].ol[I], self.ccd_voie[voie].x[I], 'o')