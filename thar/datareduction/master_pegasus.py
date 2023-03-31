import settingspegasus as sp
import extract
import os
import pathlib
import numpy as np


## the list of subdirectories to be treated
list_of_dirs = [
    '2020_0912',
    '2020_0915',
    '2020_0904',
    '2020_1017',
    '2020_0917',
    '2020_1029',
    '2020_0905',
    '2020_0911',
    '2022_0529',
    # '2020_0930',   containes only thoriums.......
    '2022_0613',
    '2020_0907',
    # '2020_1008',   contains 100 thoriums  only thoriums ........
    '2022_0617',
    '2020_0913',
    '2020_1016',
    # '2020_0908',
]
badlist = []     # muss haendisch gepflegt werden....

### generate all thorium extracts polynpomes etc...


# use kwargs from pegasus
kwargs = sp.kwargs


PEGDIR = os.path.join(kwargs['BASEDIR'], '51Peg_raw')

# reporting on bad and good
oopsies = {}
goodies = []
oopsiesstar = {}


def step1():
    nfiles = 0
    for d in list_of_dirs:
        kwargs['DATADIR'] = os.path.join(PEGDIR, d)
        for f_thar in extract.getallthoriumfits(kwargs['DATADIR']):
            nfiles += 1
        print (d, nfiles)

    nf = 0
    for d in list_of_dirs:
        kwargs['DATADIR'] = os.path.join(PEGDIR, d)
        for f_thar in extract.getallthoriumfits(kwargs['DATADIR']):
            nf += 1
            print ('\n================================')
            print (nf, 'out of', nfiles, ': ', f_thar)
            print('====================================\n')

            ## retrieve or create Extractor
            try:
                myext = extract.get_ext(f_thar, **kwargs)
                myext.ccd_voie1.get_lambda_list()
                myext.ccd_voie2.get_lambda_list()
            except Exception as ex:
                print('ooops', ex)
                oopsies[f_thar] = ex
                try:
                    del myext
                except Exception as ex:
                    print('could not delete...', ex)
                continue
            goodies.append(f_thar)
            try:
                extract.store.store(f_thar, myext)
            except:
                print ('could not save')
                oopsies[f_thar] = 'could not save'

            try:
                del myext
            except Exception as ex:
                print('problems deleting...beware of memory leak', ex)
                pass

def step2():
    nfiles = 0
    for d in list_of_dirs:
        kwargs['DATADIR'] = os.path.join(PEGDIR, d)
        for f_star in extract.getallstarfits(kwargs['DATADIR']):
            nfiles += 1
        print (d, nfiles)
    nf = 0
    for d in list_of_dirs:
        kwargs['DATADIR'] = os.path.join(PEGDIR, d)

        times = {extract.gettimestamp(f): f for f in extract.getallthoriumfits(kwargs['DATADIR'])}
        ttt = list(times.keys())

        for f_star in extract.getallstarfits(kwargs['DATADIR']):
            nf += 1
            print ('\n================================')
            print (nf, 'out of', nfiles, ': ', f_star)
            print('===================================\n')

            t = extract.gettimestamp(f_star)

            f_thar = times[
                ttt[np.argmin([abs(t - tt) for tt in ttt])]
            ]

            try:
                myext = extract.get_ext(f_thar)
                myext.set_fitsfile(f_star)
                if pathlib.Path(myext.result_path).exists():
                    continue
                myext.save_fits()
            except Exception as ex:
                oopsiesstar[f_star]=[f_thar, ex]
            try:
                del myext
            except:
                pass


if __name__=="__main__":
    step1()
    step2()

    print ( "goodies", goodies)
    print('oopsies', oopsies)
    print('oopsiesstar', oopsiesstar)


