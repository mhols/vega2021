import spectrograph
import snippets
import extract
# import thar_write
import os
import re
import pickle
import numpy as np
from astropy.io import fits

settingsmodule = os.environ.get('SETTINGSMODULE', 'settings')

try:
    exec('from %s import *'%(settingsmodule,))
except:
    raise Exception('could not import {settingsmodule}')
print('using settings from %s.py'%(settingsmodule,))


# collection of CCD2d ojects (one for each thar file and for each voie)

if RECOMPUTE_2D_POLYNOMIAL:
    thars = []

    # generate 2-d polynomial for ThAr spectra in DATADIR
    for f_thar in extract.getallthoriumfits(dirname=DATADIR):
        myext = extract.Extractor(**kwargs) # every thar gets its own extractor
        myext.set_fitsfile(f_thar)

        snips = [ snippets.Snippets(voie=i, tharfits=f_thar) for i in [1, 2, ]]  # TODO: voie3
        snippets_voie = [
            s.snippets for s in snips
            #snippets.snippets(myext, i, ORDERS)
            #for i in [1, 2]   # [1, 2, 3]
        ]

        ccd = [
            spectrograph.CCD2d( data=snip, **kwargs)
                for snip in snippets_voie
        ]

        thars.append( {
            'fitsfile':os.path.basename(f_thar), 
            'DATE_JUL': extract.header_info_from_fits(f_thar, 'DATE_JUL'),
            'ccd': ccd,
            'snippets': snippets_voie,
            'extract': myext,
        })

    with open('thars.pickle', 'wb') as f:
        pickle.dump(thars, f)
else:
    with open('thars.pickle', 'rb') as f:
        thars = pickle.load(f)

# computing average poly2d as new basis
if SAFE_LAMS:

    tmp = []
    for ta in thars:
        for c in ta['ccd']:
            tmp.append(c.get_lambda_list())
    np.savetxt(os.path.join(REFFILES, 'hobo.txt'), np.column_stack(tmp).mean(axis=1))

list_of_stars = []
list_of_jd = []

for f in extract.getallstarfits(DATADIR, STARNAME):
    with fits.open(f) as fi:
        header = fi[0].header
        jd = header['DATE_JUL']
        if header.get(PREFIX+"LEVEL", 0)>0:
            continue
    list_of_stars.append(f)
    list_of_jd.append(jd)

list_of_jdc = [c['DATE_JUL'] for c in thars ]


for d, starfits in zip(list_of_jd, list_of_stars):
    i = np.argmin(np.abs(np.array(list_of_jdc) - d))


    myext = thars[i]['extract']
    poly2 = thars[i]['ccd']
    myext.set_fitsfile(starfits)

    filename = os.path.basename(starfits)
    filename = PREFIX + filename[:-8]+'st1.fits'


    order = []
    for o in ORDERS:
        for n in range(NROWS):
            order.append(o)

    # collecting data
    fitstable = fits.BinTableHDU.from_columns(
        [
        fits.Column(
            name="true_order_number",
            format="I",
            array=np.array(order)
        ),
        fits.Column(
            name="wavelength_1",
            format="E",
            array=poly2[0].get_lambda_list()
        ),
        fits.Column(
            name="wavelength_2",
            format="E",
            array=poly2[1].get_lambda_list()
        ),
        fits.Column(
            name="wavelength_3",
            format="E",
            array=np.zeros(len(ORDERS)*NROWS)
        ),
        fits.Column(
            name='flux_1',
            format='E',
            array=myext.voie1_all()
        ),
        fits.Column(
            name='flux_2',
            format='E',
            array=myext.voie2_all()
        ),
        fits.Column(
            name='flux_3',
            format='E',
            array=myext.voie2_all() ### TODO voie3
        ),
        fits.Column(
            name='noise_1',
            format='E',
            array=myext.noise_1
        ),
        fits.Column(
            name='noise_2',
            format='E',
            array=myext.noise_2
        ),
        fits.Column(
            name='noise_3',
            format='E',
            array=myext.noise_3
        )])

        # fits.Column(
        #     name="blaze_1",
        #     format='E',
        #     array=myext.blaze_1
        # ),
        # fits.Column(
        #     name="blaze_2",
        #     format='E',
        #     array=myext.blaze_2
        # ),
        # fits.Column(
        #     name="blaze_3",
        #     format='E',
        #     array=myext.blaze_3
        # ),
        # fits.Column(
        #     name="continuum_1_1",
        #     format='E',
        #     array=myext.continuum_1_1
        # ),
        # fits.Column(
        #     name="continuum_1_2",
        #     format='E',
        #     array=myext.continuum_1_2
        # ),
        # fits.Column(
        #     name="continuum_1_3",
        #     format='E',
        #     array=myext.continuum_1_3
        # )

    page1 = fits.PrimaryHDU()

    with fits.open(starfits) as sfi:

        page1.header = sfi[0].header

        page1.header.append((PREFIX+"LEVEL", 1))
        for key in HEADER_ITEMS:
            page1.header.append((PREFIX+key, globals()[key]))

    newfits = fits.HDUList(
        [ page1, fitstable ]
    )

    newfits.writeto(os.path.join(DATADIR, filename), overwrite=True)

