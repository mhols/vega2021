import settingsHD185144
import extract
import os
import pathlib
import numpy as np


## the list of subdirectories to be treated
list_of_dirs = [
"10may22",
"12aug22",
"14apr22",
"14aug22",
"15apr22",
"15aug22",
"15jun22",
"16apr22",
"17jun22",
"18apr22",
"21jul22",
"22aug22",
"22jul22",
"22mar22",
"23aug22",
"23jul22",
"23mar22",
"25apr22",
"25jul22",
"25may22",
"26may22",
"27mar22",
"28may22",
"29may22",
]


badlist = []     # muss haendisch gepflegt werden....

### generate all thorium extracts polynpomes etc...


# use kwargs from pegasus
kwargs = settingsHD185144.kwargs


DIRDIR = os.path.join(kwargs['BASEDIR'], 'HD185144_raw')

# reporting on bad and good
oopsies = {}
goodies = []
oopsiesstar = {}


def step1():
    nfiles = 0
    for d in list_of_dirs:
        kwargs['DATADIR'] = os.path.join(DIRDIR, d)
        for f_thar in extract.getallthoriumfits(kwargs['DATADIR']):
            nfiles += 1
        print (d, nfiles)

    nf = 0
    for d in list_of_dirs:
        kwargs['DATADIR'] = os.path.join(DIRDIR, d)
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
    """
    iterate over stars and look for closest ThAr....
    """
    nfiles = 0
    for d in list_of_dirs:
        kwargs['DATADIR'] = os.path.join(DIRDIR, d)
        for f_star in extract.getallstarfits(kwargs['DATADIR']):
            nfiles += 1
        print (d, nfiles)
    nf = 0
    for d in list_of_dirs:
        kwargs['DATADIR'] = os.path.join(DIRDIR, d)

        times = {extract.gettimestamp(f): f for f in extract.getallthoriumfits(kwargs['DATADIR'])}
        ttt = list(times.keys())

        for f_star in extract.getallstarfits(kwargs['DATADIR']):
            nf += 1
            print ('\n===================================')
            print (nf, 'out of', nfiles, ': ', f_star)
            print(   '===================================\n')

            t = extract.gettimestamp(f_star)

            f_thar = times[
                ttt[np.argmin([abs(t - tt) for tt in ttt])]
            ]

            try:
                myext = extract.get_ext(f_thar, **kwargs)
                # myext.finalize()
                myext.set_fitsfile(f_star)
                if pathlib.Path(myext.result_path).exists():
                    continue
                myext.save_fits()
            except Exception as ex:
                print(ex)
                oopsiesstar[f_star]=[f_thar, ex]
            try:
                del myext
            except:
                pass


if __name__=="__main__":
    #step1()
    step2()

    print ( "goodies", goodies)
    print('oopsies', oopsies)
    print('oopsiesstar', oopsiesstar)


