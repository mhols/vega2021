from extract import *
from settings_reference import kwargs as kwargs_reference
import pickle
import util
from units import *
import matplotlib.pyplot as plt

fitsfile_reference = os.path.join(kwargs_reference['BASEDIR'], 'vega_reference/NEO_20220903_191404_th0.fits')
fitsfile_moon = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230406_190457_th0.fits')

sun_1 = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230407_010642_st0.fits')
sun_2 = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230407_010707_st0.fits')
sun_3 = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230407_010731_st0.fits')
sun_4 = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230407_010755_st0.fits')

# the following lines may be commented out when done..
mymoon = Extractor(fitsfile_moon, **kwargs_reference) # TODO moon need own setting
mymoon.update()
mymoon.update()
mymoon.update()

mymoon.save_to_store()

# mymoon is now ready to reduce any star 
mysun_1 = get_ext(fitsfile_moon)

# reduce the first sun sprectrum
mysun_1.set_fitsfile(sun_1)

# mymoon is now ready to reduce any star 
mysun_2 = get_ext(fitsfile_moon)

# reduce the first sun sprectrum
mysun_2.set_fitsfile(sun_2)

# mymoon is now ready to reduce any star 
mysun_3 = get_ext(fitsfile_moon)

# reduce the first sun sprectrum
mysun_3.set_fitsfile(sun_3)

# mymoon is now ready to reduce any star 
mysun_4 = get_ext(fitsfile_moon)

# reduce the first sun sprectrum
mysun_4.set_fitsfile(sun_4)

#----------- now do what you want with the etracted voices etc... -----#
sl = pickle.load(open(os.path.join(kwargs_reference['BASEDIR'], 
                    'reffiles/selectedlines_fit.pkl3'), 'rb'))['sellines']
sl = np.array(list(set(sl)))

ss = pickle.load(open(os.path.join(kwargs_reference['BASEDIR'], 
                    'reffiles/SolarSpectrum.pkl3'), 'rb'))

refl = ss['Wavelength']
refi = ss['Intensity']

"""
plt.figure()
res = {}
for o in mysun_1.ORDERS:
    lams, v, I = mysun_1.get_lambda_intens(1,o) # logical range
    lm = util.local_maxima(-v)  # absorption lines
    lam1, lam2 = mysun_1.lambda_range_voie1(o)

    lines = sl[ (sl>lam1) & (sl < lam2) ] 
    pix_lines = mysun_1.ccd_voie1._map_2D_x_ol_o[o](o * lines)

    mI, mJ = util.matching(lm[:, 0], pix_lines, lambda a, b: np.abs(a-b) <= 1.0 )

    tmp = []
    for i in mI:
        x, a, b = lm[i]
        try:
            A, mu, sigma, y_offset = util.estimate_location(v[a:b+1], absorption=True, **kwargs_reference)
            pixel_mean = a + mu
        except:
            pixel_mean = x

        tmp.append(pixel_mean)

    try:
        est_lambda = mysun_1.ccd_voie1._map_2D_ol_x_o[o](np.array(tmp)) / o
    except:
        continue

    mysun_1.plot_voie(1, [o])
    plt.vlines(est_lambda, 0, 1e6, 'k')
    plt.vlines(lines[mJ], 0, 1e6, 'r')

    res[o] = C_LIGHT * (1-est_lambda/ lines[mJ])
"""

res = {}
for o in mysun_1.ORDERS:
    lams, v, I = mysun_1.get_lambda_intens(1,o) # logical range
    lams = lams[I]
    v = v[I]
    lm = util.local_maxima(-v)  # absorption lines
    lam1, lam2 = mysun_1.lambda_range_voie1(o)

    extl = refl[ (refl>lam1) & (refl < lam2) ] 
    exti = refi[ (refl>lam1) & (refl < lam2) ] 
    reflm = util.local_maxima(-exti)

    pix_lines = mysun_1.ccd_voie1._map_2D_x_ol_o[o](o * reflm[:, 0])

    mI, mJ = util.matching(lm[:, 0], pix_lines, lambda a, b: np.abs(a-b) <= 2.0 )

    print(o, len(mI))

    tmp = []
    for i, j in zip(mI, mJ):
        x, a, b = lm[i]
        xx, aa, bb = reflm[j]
        try:
            A, mu, sigma, y_offset = util.estimate_location(v[a:b+1], absorption=True, **kwargs_reference)
            pixel_mean = a + mu
            AA, mumu, sigsig, yy_offset = util.estimate_location(exti[aa:bb+1], absorption=True, **kwargs_reference)
            pixel_meanmean = aa + mumu
            tmp.append([A, a + mu, sigma, y_offset, AA, aa + mumu, sigsig, yy_offset])
        except Exception as ex:
            print('cucou', ex)
            continue

    tmp = np.array(tmp)
    try:
        est_lambda = mysun_1.ccd_voie1._map_2D_ol_x_o[o](tmp[:,1]) / o
        est_lambdalam =  mysun_1.ccd_voie1._map_2D_ol_x_o[o](tmp[:,5]) / o
    except Exception as ex:
        print (ex)
        continue


    plt.plot(est_lambdalam, tmp[:,0]*tmp[:,7]/(tmp[:,3]*tmp[:,4]), '.b')


