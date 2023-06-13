from extract import *
from settings_reference import kwargs as kwargs_reference
import pickle
import util
from units import *
import matplotlib.pyplot as plt

fitsfile_reference = os.path.join(kwargs_reference['BASEDIR'], 'vega_reference/NEO_20220903_191404_th0.fits')
fitsfile_hd108 = os.path.join(kwargs_reference['BASEDIR'], 'hd108/RAW/NEO_20211007_175136_th0.fits')

hd108_1 = os.path.join(kwargs_reference['BASEDIR'], 'hd108/RAW/NEO_20211007_232240_st0.fits')
hd108_2 = os.path.join(kwargs_reference['BASEDIR'], 'hd108/RAW/NEO_20211007_234602_st0.fits')
hd108_3 = os.path.join(kwargs_reference['BASEDIR'], 'hd108/RAW/NEO_20211008_000925_st0.fits')
hd108_4 = os.path.join(kwargs_reference['BASEDIR'], 'hd108/RAW/NEO_20211008_003247_st0.fits')

# the following lines may be commented out when done..

myhd108 = Extractor(fitsfile_hd108, **kwargs_reference) # TODO hd108 need own setting
myhd108.update()
myhd108.update()
myhd108.update()

myhd108.save_to_store()


# myhd108 is now ready to reduce any star
myhd108_1 = get_ext(fitsfile_hd108)
myhd108_1.update_kwargs(VOIE_METHOD ='SUM_DIVIDE')

# reduce the first hd108 sprectrum
myhd108_1.set_fitsfile(hd108_1)

# myhd108 is now ready to reduce any star
myhd108_2 = get_ext(fitsfile_hd108)
myhd108_2.update_kwargs(VOIE_METHOD ='SUM_DIVIDE')

# reduce the first hd108 sprectrum
myhd108_2.set_fitsfile(hd108_2)

# myhd108 is now ready to reduce any star
myhd108_3 = get_ext(fitsfile_hd108)
myhd108_3.update_kwargs(VOIE_METHOD ='SUM_DIVIDE')

# reduce the first hd108 sprectrum
myhd108_3.set_fitsfile(hd108_3)

# myhd108 is now ready to reduce any star
myhd108_4 = get_ext(fitsfile_hd108)
myhd108_4.update_kwargs(VOIE_METHOD ='SUM_DIVIDE')

# reduce the first hd108 sprectrum
myhd108_4.set_fitsfile(hd108_4)

myhd108_1.voie
myhd108_2.voie
myhd108_3.voie
myhd108_4.voie

myhd108_1.save_to_store()
myhd108_2.save_to_store()
myhd108_3.save_to_store()
myhd108_4.save_to_store()





myhd108_1 = get_ext(hd108_1)
myhd108_2 = get_ext(hd108_2)
myhd108_3 = get_ext(hd108_3)
myhd108_4 = get_ext(hd108_4)

myhds = [myhd108_1, myhd108_2, myhd108_3, myhd108_4]

lam_voie1={}
lam_voie2={}
intens_voie1={}
intens_voie2={}
I = {}
intens_voie2_interp = {}
weight_voie1={}
weight_voie2={}


for i, mh in enumerate(myhds):
    lam_voie1[i+1], intens_voie1[i+1], I[i+1] = {}, {}, {}
    lam_voie2[i+1], intens_voie2[i+1], I[i+1] = {}, {}, {}
    intens_voie2_interp[i+1] = {}
    weight_voie1[i+1] = {}
    weight_voie2[i+1] = {}
    for o in mh.ORDERS:
        lam_voie1[i+1][o], intens_voie1[i+1][o], I[o] = mh.get_lambda_intens1(o)
        lam_voie2[i+1][o], intens_voie2[i+1][o], I[o] = mh.get_lambda_intens2(o)
        intens_voie2_interp[i+1][o] = mh.voie2_lam1(o)

        weight_voie1[i+1][o] = np.median(mh.bare_voie1[o][mh.CENTRALROW-100:mh.CENTRALROW+100])
        sumweight = sumweight + weight_voie1[i+1][o]
"""
we now have z.B. lam_voie1[3][23]
"""
bv11 = myhd108_1.bare_voie1
bv12 = myhd108_1.bare_voie2
bv21 = myhd108_2.bare_voie1
bv22 = myhd108_2.bare_voie2
bv31 = myhd108_3.bare_voie1
bv32 = myhd108_3.bare_voie2
bv41 = myhd108_4.bare_voie1
bv42 = myhd108_4.bare_voie2


xlim=[5884,5898]
ylim=[0,0.15]


"""
def plot1():
    plt.figure()
    plt.title("spectrum 1")
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    plt.plot(lam_1_voie1,intens_1_voie1)
    plt.plot(lam_1_voie2,intens_1_voie2)

def plot2():
    plt.figure()
    plt.title("spectrum 2")
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    plt.plot(lam_2_voie1,intens_2_voie1)
    plt.plot(lam_2_voie2,intens_2_voie2)

def plot3():
    plt.figure()
    plt.title("spectrum 3")
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    plt.plot(lam_3_voie1,intens_3_voie1)
    plt.plot(lam_3_voie2,intens_3_voie2)

def plot4():
    plt.figure()
    plt.title("spectrum 4")
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    plt.plot(lam_4_voie1,intens_4_voie1)
    plt.plot(lam_4_voie2,intens_4_voie2)

def plotq():
    plt.figure(figsize=(16,6))
    plt.title("q")
    plt.xlim(*xlim)
    plt.plot(lam_1_voie1,q1)
    plt.plot(lam_2_voie1,q2)
    plt.plot(lam_3_voie1,q3)
    plt.plot(lam_4_voie1,q4)

def allplot():
    plot1()
    plot2()
    plot3()
    plot4()

def plot_meanv():
    plt.figure(figsize=(16,6))
    plt.title("meanV_over_I")
    plt.xlim(*xlim)
    plt.plot(lam_1_voie1,mean_V_over_I)
"""

def outputfile():
    res1 = []
    res2 = []
    res3 = []
    res4 = []
    res5 = []
    mask = []
    NROWS = myhd108_1.NROWS

    for o in myhd108_1.ORDERS[::-1]:
    #lambda
        lam_1_voie1,intens_1_voie1,I = myhd108_1.get_lambda_intens1(o)
        res1.extend(lam_1_voie1/10.)
        for i in range(NROWS):
            if i in I:
                mask.append(1)
            else:
                mask.append(0)

    #normierte total intensity (4 * 2 voies)
        intens_total=(weight_voie1[1][o]*intens_voie1[1][o]+weight_voie2[1][o]*intens_voie2_interp[1][o]\
        +weight_voie1[2][o]*intens_voie1[2][o]+weight_voie2[2][o]*intens_voie2_interp[2][o]\
        +weight_voie1[3][o]*intens_voie1[3][o]+weight_voie2[3][o]*intens_voie2_interp[3][o]\
        +weight_voie1[4][o]*intens_voie1[4][o]+weight_voie2[4][o]*intens_voie2_interp[4][o])\
        (weight_voie1[1][o]+weight_voie2[1]+weight_voie1[2][o]+weight_voie2[2][o]+\
        weight_voie1[3][o]+weight_voie2[3]+weight_voie1[4][o]+weight_voie2[4][o])
        #--------------------------------
        l = np.arange(len(intens_total))
        nnodes=10
        q=0.3
        qq=0.7
        qqq=0.95
        p = util.continuum(l,intens_total,nnodes,q,qq,qqq)
        intens_total_norm = intens_total/p(l)
        #mask of usable order
        #--------------------------------

        intens_total_norm = intens_total_norm * mask
        res2.extend(intens_total_norm)


    #Stokes V/I
        print("jetzt Stokes V")

        q1 = myhd108_1.voie1[o]/myhd108_1.voie2[o]
        q2 = myhd108_2.voie1[o]/myhd108_2.voie2[o]
        q3 = myhd108_3.voie1[o]/myhd108_3.voie2[o]
        q4 = myhd108_4.voie1[o]/myhd108_4.voie2[o]

        R4 = (q1/q2)*(q4/q3)
        Rnull4 = (q1/q4) * (q2/q3)

        R = R4**(1/4.)
        Rnull = Rnull4**(1/4.)

        mean_V_over_I = (R-1.)/(R+1.)
        mean_N_over_I = (Rnull-1.)/(Rnull+1.)

        mean_V_over_I = mean_V_over_I * mask
        mean_N_over_I = mean_N_over_I * mask
        res3.extend(mean_V_over_I)
        res4.extend(mean_N_over_I)


        print("jetzt noise")



    #noise
        inte=(bv11[o]+bv12[o]+bv21[o]+bv22[o]+\
        bv31[o]+bv32[o]+bv41[o]+bv42[o])

        print(o)

        l = np.arange(len(inte))
        nnodes=10
        q=0.3
        qq=0.7
        qqq=0.95
        p = util.continuum(l,inte,nnodes,q,qq,qqq)

        lam_1_voie1,intens_1_voie1,I = myhd108_1.get_lambda_intens1(o)
        mask=[]
        NROWS = kwargs_reference['NROWS']
        for i in range(NROWS):
            if i in I:
                mask.append(1)
            else:
                mask.append(0)

        noise = np.sqrt(inte)/p(l)
        noise=noise*mask
        res5.extend(noise)


    print("im here")




    #return col1,col2,col3,col4,col5,col6
    with open("hd108_lsdtest.s","w") as f:
        f.write("dies ist ein Versuch eines polarisationsspektrums von hd108\n")
        f.write(str(len(col1)) + "   5\n")
        np.savetxt(f,np.column_stack([res1,res2,res3,res4,res5,res6]))





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
    for o in myhd108_1.ORDERS:
        lams, v, I = myhd108_1.get_lambda_intens(1,o) # logical range
        lm = util.local_maxima(-v)  # absorption lines
        lam1, lam2 = myhd108_1.lambda_range_voie1(o)

        lines = sl[ (sl>lam1) & (sl < lam2) ]
        pix_lines = myhd108_1.ccd_voie1._map_2D_x_ol_o[o](o * lines)

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
            est_lambda = myhd108_1.ccd_voie1._map_2D_ol_x_o[o](np.array(tmp)) / o
        except:
            continue

        myhd108_1.plot_voie(1, [o])
        plt.vlines(est_lambda, 0, 1e6, 'k')
        plt.vlines(lines[mJ], 0, 1e6, 'r')

        res[o] = C_LIGHT * (1-est_lambda/ lines[mJ])
    """

    res = {}
    f = 1 + 12 * KM / S / C_LIGHT
    for o in myhd108_1.ORDERS:

        lmin, lmax = myhd108_1.lambda_range_voie1(o)

        reflams = sl[ (sl >= lmin) & (sl <= lmax) ]

        l, v, I = myhd108_1.get_lambda_intens1(o)
        l=l[I]
        v=v[I]
        bg = util.background(l, v, nnodes=10)
        v = v/bg(l)

        hd108l, hd108v = util.extract_snippets(l, v, reflams/f, reflams*f)
        rl, rv =  util.extract_snippets(refl, refi, reflams/f, reflams*f)

        N = 128

        """
        for l, v in zip(hd108l, hd108v):
            plt.plot(l, v, 'r-')

        for l, v in zip(rl, rv):
            plt.plot(l, 100000 * v, 'b-')

        """
        for l, v, ll, vv in zip(hd108l, hd108v, rl, rv):
            try:
                A1, mu1, sigma1, y_offset1 = util.estimate_igauss(l, v)
                A2, mu2, sigma2, y_offset2 = util.estimate_igauss(ll, vv)

                ew1 = A1 / y_offset1
                ew2 = A2 / y_offset2
                plt.plot( [mu1],  [ew1 / ew2], 'o', color=myhd108_1.ccd_voie1.color_of_order(o))
            except Exception as ex:
                print (l, v)
                print (ex)

