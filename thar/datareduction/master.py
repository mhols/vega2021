from settings import *
import spectrograph
import snippets
import extract
# import thar_write
import os
import re


myext = extract.Extractor(**kwargs)
# collection of CCD2d ojects (one for each thar file and for each voie)

if RECOMPUTE_2D_POLYNOMIAL:
    ccds = {}

    # generate 2-d polynomial for ThAr spectria in DATADIR
    for f_thar in extract.getallthoriumfits(dirname=DATADIR):
        try:
            myext.set_fitsfile(f_thar)
        except Exception as ex:
            print (ex)
            continue

        snippets_voie = [
            snippets.snippets(myext, i, ORDERS)
            for i in [1, 2]   # [1, 2, 3]
        ]

        ccd = [
            spectrograph.CCD2d(
                data=snip,
                bootstrap_data = True,
                file_lambda_list = os.path.join(
                    RESFILES, os.path.basename(f_thar)+".wave2D_{}.dat".format(i)
                ),
                **kwargs
            )
            for i, snip in enumerate(snippets_voie)
        ]

        ccds[f_thar] = {'ccd': ccd, 'snippets': snippets_voie}
