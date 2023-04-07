#from settingspegasus import kwargs
from settingsmoon import kwargs
from numpy import *
import matplotlib.pyplot as plt
import extract as ex


myext = ex.get_ext(
        #'/home/hols/vega2021/thar/51Peg_raw/2020_1029/NEO_20201029_172538_th0.fits',
        "/home/hols/vega2021/thar/lune_raw/NEO_20200202_173811_th0.fits",
        finalize=False,
        **kwargs)

#del myext.snippets_voie1
#del myext.snippets_voie2

s1 = myext.snippets_voie1
s2 = myext.snippets_voie2

Ig1 = s1.good_snippet
Ig2 = s2.good_snippet

I = unique(intersect1d(s1.ci.values, s2.ci.values))

def plot_scatter():
    plt.plot(s1.ci, s1.x, 'b.')
    plt.plot(s2.ci, s2.x, 'r.' )

    # common snippets


    tmp = []

    for i in I:
        sn1 = s1.sn[s1.ci==i]
        sn2 = s2.sn[s2.ci==i]

        plt.plot(len(sn1)*[i-0.1], sn1['pixel_mean'], 'go')
        plt.plot(len(sn2)*[i], sn2['pixel_mean'], 'co')

        if len(sn1) == 1 and len(sn2) == 1:
            plt.plot([i,i], [sn1['pixel_mean'].iloc[0], sn2['pixel_mean'].iloc[0]], 'k')
            tmp.append(sn1.iloc[0]['pixel_mean']-sn2.iloc[0]['pixel_mean'])

    return tmp

def plot_snippet(ic):
    sn1 = s1.sn[s1.ci==ic]
    sn2 = s2.sn[s2.ci==ic]

    for s in [sn1, sn2]:
        for i in s.index:
            try:
                a, b = s.loc[i].left, s.loc[i].right
                o = s.loc[i].true_order_number
    
                ab = arange(a,b+1)
                intens = s.loc[i].bare_voie
                plt.plot(ab, intens)
            except:
                pass

def plot_snippet_voie(ss, o):

    s = ss.sn[ss.o==o]
    l = myext.lambdas_per_order_voie1[o]
    plt.vlines(ss.atlasext(o), 0, 100, 'y')
    plt.vlines(s['ref_lambda'], 0, 100, 'g')
    plt.plot(l, ss.bare_voie(o), '-k')

    for i in s.index:
        a,b = s.loc[i].left, s.loc[i].right
        ab = arange(a,b+1)
        intens = s.loc[i].bare_voie
        plt.plot(l[ab], intens, '-r')

                 

myext.save_to_store()