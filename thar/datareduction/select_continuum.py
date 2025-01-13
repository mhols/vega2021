import extract
from extract import *  #  for the pickle to work
import pathlib
import matplotlib.pyplot as plt
import numpy as np
import units
import pandas as pd
import os
import util

picklename = 'catselectpick.pickle'

lockfile = '.lock_cklick'

tharfile =  '/Users/boehm/Desktop/vega2021/thar/06apr23_Moon/NEO_20230406_190457_th0.fits'

print(tharfile)

#tharfile = pathlib.Path(__file__).parents[0] /  'vega_reference/NEO_20220903_191404_th0.fits'

class ClickSelector:

    def __init__(self, o=35, voie=1):

        self.fig = plt.figure('selector')
        self.ax = self.fig.add_subplot(111)
        self.fig.canvas.mpl_connect("button_press_event",self.toggle_click_select)
        self.myext = extract.get_ext(tharfile)

        if not os.path.exists(lockfile):
            self.myext.update()
            self.myext.update()
            self.myext.update()
            self.myext.save_to_store()

        if not os.path.exists(lockfile):
            self.cr = self.myext._snippets_manager_voie1._atlas['UVES']
            self.cr['click_selected'] = False
            # self.cr['true_order_number'] = 0
        else:
            self.cr = pd.read_pickle(picklename)

        if not os.path.exists(lockfile):
            pathlib.Path(lockfile).touch()
        
        self.o = o
        self.voie=voie

        self.plot_catalog_selected()

    @property
    def v(self):
        return extract.util.clean_nans(self.myext.voie[self.voie][self.o][self.myext.I[self.o]])

    @property
    def l(self):
        return self.myext.lambdas_per_order_voie[self.voie][self.o]

    def restart(self):
        os.path.rm(lockfile)

    def set_order(self, o):
        self.ax.clear()
        self.o = o
        self.plot_catalog_selected()

    def index_range_catalog(self):
        lam1, lam2 = self.myext.lambda_range_voie1(self.o)

        J = (self.cr["ref_lambda"] > lam1) & (self.cr["ref_lambda"] < lam2)

        return J
    
    def _plot_selected(self):

        J = self.index_range_catalog()
    
        I = (self.cr['click_selected'] == True ) & J
        
        plt.vlines(self.cr.loc[I, 'ref_lambda'], 0, np.max(self.v), 'g')
        plt.plot(self.cr.loc[I, 'ref_lambda'], sum(I) * [0], 'go')
         
    def _plot_removed(self):
        lam1, lam2 = self.myext.lambda_range_voie1(self.o)

        J = (self.cr["ref_lambda"] > lam1) & (self.cr["ref_lambda"] < lam2)
    
        I = (self.cr['click_selected'] == False) & J
        plt.vlines(self.cr.loc[I, 'ref_lambda'], 0, np.max(self.v), colors='b', linestyles='dashed' )
        plt.plot(self.cr.loc[I, 'ref_lambda'], sum(I) * [0], 'bo')
        

    def plot_catalog_selected(self):

        #self.clear_axis()
        self._plot_selected()
        self._plot_removed()
        self.myext.plot_voie(self.voie,[self.o])

    def plot_selected(self):
        self.clear_axis()
        self._plot_selected()
        self.myext.plot_voie(self.voie, [self.o])    
    
    def plot_removed(self):
        self.clear_axis()
        self._plot_removed()
        self.myext.plot_voie(self.voie, [self.o])    

    def clear_axis(self):
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        self.ax.clear()
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)

    def toggle_click_select(self, event):
        lam = event.xdata

        I = np.abs(self.cr['ref_lambda']/lam - 1)*units.C_LIGHT < 10 * units.KM / units.S
        if sum(I) == 1:
            self.clear_axis()
            i = self.cr.index[I][0]

            selected = self.cr.loc[i,'ref_lambda']
            self.cr.loc[i, 'click_selected'] = not self.cr.loc[i, 'click_selected']
 
            print (selected)

            self.plot_catalog_selected()
            self.cr.to_pickle(picklename)

            self.report(selected)

    def report(self, selected):

        x = []
        lam = []
        w = []

        J = self.index_range_catalog()

        I = self.cr['click_selected'] & J
        s = self.myext.snippets_voie[self.voie].find_snippet(selected)
        self.myext.snippets_voie[self.voie]._snippets.loc[s.index, 'ref_lambda'] = selected
        for l in self.cr.loc[I, 'ref_lambda']:
            s = self.myext.snippets_voie[self.voie].find_snippet(l)
            s = s[s['true_order_number'] == self.o]

            a = s.iloc[0]['left']
            b = s.iloc[0]['right']


            v = self.myext.snippets_voie[self.voie].bare_voies[self.o][a:b+1]

            try:
                A, mu, sigma, offset = util.estimate_location(v, **self.myext.kwargs)
            except Exception as ex:
                print('could not estimate Gaussian fit', ex)
                continue

            x.append(mu + a)
            lam.append(l)
            w.append(A)
        
        x = np.array(x)
        lam = np.array(lam)
        w = np.array(w)
        try:
            p = np.polynomial.Polynomial.fit(x, lam, w=np.sqrt(w), deg=7)
        except:
            print ('no fit possible')
            return
        print ( (1-p(x)/lam) * units.C_LIGHT)
        # print(x, lam)


        
        # for n in range(5, 10):
            

#plot_catalog_selected()
#myext.plot_voie(1, [o])
   
if __name__=="__main__":
    myclick = ClickSelector()




