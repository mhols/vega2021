from extract import *
from settings_reference import kwargs as kwargs_reference
import pickle
import util
from units import *

fitsfile_reference = os.path.join(kwargs_reference['BASEDIR'], 'vega_reference/NEO_20220903_191404_th0.fits')
fitsfile_moon = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230406_190457_th0.fits')

sun_1 = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230407_010642_st0.fits')
sun_2 = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230407_010707_st0.fits')
sun_3 = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230407_010731_st0.fits')
sun_4 = os.path.join(kwargs_reference['BASEDIR'], '06apr23_Moon/NEO_20230407_010755_st0.fits')

# the following lines may be commented out when done..
#mymoon = Extractor(fitsfile_moon, **kwargs_reference) # TODO moon need own setting
#mymoon.update()
#mymoon.update()
#mymoon.update()

#mymoon.save_to_store()

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
res = {}
for o in [33,34]: #mysun_1.ORDERS:
    lams, v, I = mysun_1.get_lambda_intens(1,o) # logical range
    lm = util.local_maxima(-v)  # absorption lines
    lam1, lam2 = mysun_1.lambda_range_voie1(o)

    lines = sl[ (sl>lam1) & (sl < lam2) ] 
    pix_lines = mysun_1.ccd_voie1._map_2D_x_ol_o[o](o * lines)

    mI, mJ = util.matching(lm[:, 0], pix_lines, lambda a, b: np.abs(a-b) <= 2.0 )

    tmp = []
    for i in mI:
        x, a, b = lm[i]
        A, mu, sigma, y_offset = util.estimate_location(-v[a:b+1], **kwargs_reference)
        pixel_mean = a + mu

        tmp.append(pixel_mean)

    est_lambda = mysun_1.ccd_voie1._map_2D_ol_x_o[o](np.array(tmp)) / o
    res[o] = C_LIGHT * (1-est_lambda/ sl[mJ])