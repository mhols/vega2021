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


def doit():
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
    f = 1 + 12 * KM / S / C_LIGHT
    for o in mysun_1.ORDERS:

        lmin, lmax = mysun_1.lambda_range_voie1(o)

        reflams = sl[ (sl >= lmin) & (sl <= lmax) ]

        l, v, I = mysun_1.get_lambda_intens1(o)
        l=l[I]
        v=v[I]
        bg = util.background(l, v, nnodes=10)
        v = v/bg(l)

        sunl, sunv = util.extract_snippets(l, v, reflams/f, reflams*f)
        rl, rv =  util.extract_snippets(refl, refi, reflams/f, reflams*f)

        N = 128

        """
        for l, v in zip(sunl, sunv):
            plt.plot(l, v, 'r-')

        for l, v in zip(rl, rv):
            plt.plot(l, 100000 * v, 'b-')

        """
        for l, v, ll, vv in zip(sunl, sunv, rl, rv):
            try:
                A1, mu1, sigma1, y_offset1 = util.estimate_igauss(l, v)
                A2, mu2, sigma2, y_offset2 = util.estimate_igauss(ll, vv)

                ew1 = A1 / y_offset1
                ew2 = A2 / y_offset2
                plt.plot( [mu1],  [ew1 / ew2], 'o', color=mysun_1.ccd_voie1.color_of_order(o))
            except Exception as ex:
                print (l, v)
                print (ex)
       