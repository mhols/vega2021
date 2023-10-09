import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import nextra as nx

reffile =  os.path.join(os.path.dirname(__file__),"../06apr23_Moon/NEO_20230406_190457_th0.fits")
kwargs = nx.settings_reference.get_kwargs()
kwargs['USE_PICKED_CENTRAL_POSITIONS']=False

myref = nx.get_ext(reffile,**kwargs)
myref.save_to_store()

"""
Index(['true_order_number', 'index_pixel_snippet', 'posmax', 'a', 'b', 'left',
       'right', 'pixel_A', 'pixel_mu', 'pixel_sigma', 'pixel_offset',
       'pixel_mean', 'pixel_std', 'pixel_sum_intens', 'pixel_max_intens',
       'total_flux', 'pixel_range', 'bare_voie', 'bootstraped', 'selected',
       'goodsnippet', 'usablesnippet', 'gaussfit', 'ref_lambda', 'lambda_left',
       'lambda_right', 'est_lambda', 'x_2D_ol_o', 'l_2D_x_o', 'dvrad_2D'],
      dtype='object')
"""
def _I(myref,voie,o):
    snippets = myref.snippets_voie[voie]
    good = snippets.goodsnippet
    Io = (snippets.true_order_number == o)
    I = (good & Io)
    return I


def spectral_resolution(myref,voie,o):
    snippets = myref.snippets_voie[voie]
    p = myref.pix_to_lambda_map_voie[voie][o]
    
    #all good snippets in trueorder o
    I = _I(myref,voie,o)
    snippets = snippets[I]
    
    pix_sig = snippets.pixel_sigma
    pixel_mean = snippets.pixel_mean
    dl_dx = p.df(pixel_mean)
    lam_sig = pix_sig * dl_dx
    
    ref_lam = snippets.ref_lambda
    
    #spectral resolution of the ThAr spectrum
    
    resol = ref_lam / (2*np.sqrt(2*np.log(2))*lam_sig)
    return resol
    
def plot_resolution(myref,voie):
    for o in myref.ORDERS:
        res = spectral_resolution(myref,voie,o)
        title = f"median spectral resolution of voie {voie}"
        plt.title(title)
        plt.xlabel("true order number")
        plt.ylabel("R - spectral resolution")
        plt.plot(o,np.median(res),"o")
        
def pix_frac_err(myref,voie,o):
    snippets = myref.snippets_voie[voie]
    I = _I(myref,voie,o)
    snippets = snippets[I]
    delta_vrad = (1 - snippets.ref_lambda/snippets.est_lambda)*nx.units.C_LIGHT
    pix_frac = snippets.pixel_mean-np.floor(snippets.pixel_mean)
    return pix_frac,delta_vrad
    
def plot_frac(myref,voie):
    for o in myref.ORDERS:
        pix_frac,delta_vrad = pix_frac_err(myref,voie,o)
        title = f"radial velocity error (m/s)"
        plt.title(title)
        plt.xlabel("fractionary pixel position of gaussian centroid")
        plt.ylabel("m/s")
        plt.plot(pix_frac,delta_vrad,"o")
        
def position_error(myref,voie,o):
    snippets = myref.snippets_voie[voie]
    I = _I(myref,voie,o)
    snippets = snippets[I]
    delta_vrad = (1 - snippets.ref_lambda/snippets.est_lambda)*nx.units.C_LIGHT
    return delta_vrad
    
def plot_position_error(myref,voie):
    deltavtot = []
    for o in myref.ORDERS:
        delta_vrad = position_error(myref,voie,o)
        deltavtot.extend(delta_vrad)
        
        plt.plot(o,np.quantile(np.abs(delta_vrad),0.68),"o")

    onesigmaquantile_vrad = np.quantile(np.abs(deltavtot),0.68)
    total_lines = np.array(deltavtot).size
    title = f"total of {total_lines} lines provide a one sigma error of {onesigmaquantile_vrad} (m/s)"
    plt.title(title)
    plt.xlabel("true order number")
    plt.ylabel("radial velocity error (m/s)")
    

        
    
    


