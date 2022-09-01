import settings
import spectrograph
import thar_write
import os
import re

kwargs = settings.kwargs

fitsfiles = os.listdir(settings.DATADIR) # so werden all files in DATADIR bearbeitet
refname = os.path.join(settings.REFFILES, "thar_spec_MM201006.dat")
atlasname = os.path.join(settings.REFFILES, "thar_UVES_MM090311.dat")


# bootstrapping the data
for f in fitsfiles:
    if not re.match(r'^fo.*\.fits$', f):
        continue
    f = os.path.join(settings.DATADIR, f)
    print ('working on :', f)

    basename = os.path.basename(f)
    basename_stripped = basename.replace('.fits', '')


    snippet_files = [
        os.path.join(
            settings.TMPFILES, 
            basename_stripped+'.snippet.voie_{}.json'.format(i)) for i in [1,] 
    ]

    ### making the snippets for the voices 
    thar_write.thar_write(refname, atlasname, os.path.abspath(f), *snippet_files)

    for voie, snippet in zip ([1, 2, 3], snippet_files):
        print('voie ', voie)
        bootstrapped_file = os.path.join(
            settings.TMPFILES, 
            basename_stripped+'.bstr.voie_{}.json'.format(voie)
        )
        file_lambda_list = os.path.join(
            settings.RESFILES, 
            basename_stripped+'.wave2D.voie_{}.dat'.format(voie)
        )
        kwargs.update({
            'datafile': snippet,
            'bootstrap_data': True,
            'save_bootstraped_data': True,
            'bootstraped_file': bootstrapped_file,
            'file_lambda_list': file_lambda_list
        })

        data = spectrograph.CCD2d(**kwargs)
        data.sigma_clipping()
        data.get_lambda_list()



