from extract import *
from settings_reference import kwargs as kwargs_moon
import pickle
import util
from units import *
import matplotlib.pyplot as plt
import sys

#fitsfile_reference = os.path.join(kwargs_moon['BASEDIR'], 'vega_reference/NEO_20220903_191404_th0.fits')

starname = 'moon'
# fitsfile_tharname = os.path.join(kwargs_moon['BASEDIR'], '06apr23_Moon/NEO_20230406_190457_th0.fits')
starname_1 = os.path.join(kwargs_moon['BASEDIR'], '06apr23_Moon/NEO_20230407_010642_st0.fits')
starname_2 = os.path.join(kwargs_moon['BASEDIR'], '06apr23_Moon/NEO_20230407_010707_st0.fits')
starname_3 = os.path.join(kwargs_moon['BASEDIR'], '06apr23_Moon/NEO_20230407_010731_st0.fits')
starname_4 = os.path.join(kwargs_moon['BASEDIR'], '06apr23_Moon/NEO_20230407_010755_st0.fits')

# the following lines may be commented out when done..


kwargs_moon['SETTING_ID'] = 'MOON'

myhds = []
for starname in [starname_1, starname_2, starname_3, starname_4]:
    star = reduce_star(starname, **kwargs_moon)
    star.voie
    star.save_to_store()
    myhds.append(star)

lam_voie1={}
lam_voie2={}
intens_voie1={}
intens_voie2={}
I = {}
intens_voie2_interp = {}
weight_voie1={}
weight_voie2={}

for i, mh in enumerate(myhds):
    print("spectrum", i)
    lam_voie1[i+1], intens_voie1[i+1], I[i+1] = {}, {}, {}
    lam_voie2[i+1], intens_voie2[i+1], I[i+1] = {}, {}, {}
    intens_voie2_interp[i+1] = {}
    weight_voie1[i+1] = {}
    weight_voie2[i+1] = {}
    for o in mh.ORDERS:
        print("order", o)
        lam_voie1[i+1][o], intens_voie1[i+1][o], I[o] = mh.get_lambda_intens1(o)
        lam_voie2[i+1][o], intens_voie2[i+1][o], I[o] = mh.get_lambda_intens2(o)
        lam_interp,intens_voie2_interp[i+1][o] = mh.voie2_lam1(o)
        weight_voie1[i+1][o] = np.median(mh.bare_voie1[o][mh.CENTRALROW-100:mh.CENTRALROW+100])
        weight_voie2[i+1][o] = np.median(mh.bare_voie2[o][mh.CENTRALROW-100:mh.CENTRALROW+100])

for mh in myhds:
    mh.save_to_store()

# we now have z.B. lam_voie1[3][23]
bv11 = mystarname_1.bare_voie1
bv12 = mystarname_1.bare_voie2
bv21 = mystarname_2.bare_voie1
bv22 = mystarname_2.bare_voie2
bv31 = mystarname_3.bare_voie1
bv32 = mystarname_3.bare_voie2
bv41 = mystarname_4.bare_voie1
bv42 = mystarname_4.bare_voie2


xlim=[5884,5898]
ylim=[0,0.15]


#def outputfile():

res1 = []
res1a = []
res2 = []
res3 = []
res4 = []
res5 = []
res6 = []
mask = []
NROWS = mystarname_1.NROWS

for o in mystarname_1.ORDERS[::-1]:
#for o in mystarname_1.ORDERS:
#lambda
    lam_1_voie1,intens_1_voie1,I = mystarname_1.get_lambda_intens1(o)
    res1.extend(lam_1_voie1/10.)
    res1a.extend(lam_1_voie1)
    mask = []
    for i in range(NROWS):
        if i in I:
            mask.append(1)
        else:
            mask.append(0)


#sum of all flatfielded intensities, weighted by central flux, of all 4 exposures * 2voies
#   I
    intens_total=(weight_voie1[1][o]*intens_voie1[1][o]+weight_voie2[1][o]*intens_voie2_interp[1][o]\
    +weight_voie1[2][o]*intens_voie1[2][o]+weight_voie2[2][o]*intens_voie2_interp[2][o]\
    +weight_voie1[3][o]*intens_voie1[3][o]+weight_voie2[3][o]*intens_voie2_interp[3][o]\
    +weight_voie1[4][o]*intens_voie1[4][o]+weight_voie2[4][o]*intens_voie2_interp[4][o])\
    /(weight_voie1[1][o]+weight_voie2[1][o]+weight_voie1[2][o]+weight_voie2[2][o]+\
    weight_voie1[3][o]+weight_voie2[3][o]+weight_voie1[4][o]+weight_voie2[4][o])

    #Ic    continuum fit of above summed intens_total
    l = np.arange(len(intens_total))
    nnodes=10
    q=0.3
    qq=0.7
    qqq=0.95
    p = util.continuum(l,intens_total,nnodes,q,qq,qqq)

    #I/Ic
    intens_total_norm = intens_total/p(l)
    #--------------------------------

    #selection of the useful part of the order
    intens_total_norm = intens_total_norm * mask
    res2.extend(intens_total_norm)
    """
    if (o % 2) == 0:
        plt.plot(lam_1_voie1,intens_total,'r')
    else:
        plt.plot(lam_1_voie1,intens_total,'b')
    """

#Stokes V/I
    print("jetzt Stokes V")

    q1 = mystarname_1.voie1[o]/mystarname_1.voie2[o]
    q2 = mystarname_2.voie1[o]/mystarname_2.voie2[o]
    q3 = mystarname_3.voie1[o]/mystarname_3.voie2[o]
    q4 = mystarname_4.voie1[o]/mystarname_4.voie2[o]

    R4 = (q1/q2)*(q4/q3)
    Rnull4 = (q1/q4) * (q2/q3)
    Rnull24 = (q1/q2) * (q3/q4)

    R = R4**(1/4.)
    Rnull = Rnull4**(1/4.)
    Rnull2 = Rnull24**(1/4.)

    mean_V_over_I = (R-1.)/(R+1.)
    mean_N_over_I = (Rnull-1.)/(Rnull+1.)
    mean_N2_over_I = (Rnull2-1.)/(Rnull2+1.)

#Stokes V/Ic
    mean_V_over_I = mean_V_over_I * intens_total_norm * mask
#Null1/Ic
    mean_N_over_I = mean_N_over_I * intens_total_norm * mask
    mean_N2_over_I = mean_N2_over_I * intens_total_norm * mask

    res3.extend(mean_V_over_I)
    res4.extend(mean_N_over_I)
    res5.extend(mean_N2_over_I)


    print("jetzt noise")



#noise
    inte=(bv11[o]+bv12[o]+bv21[o]+bv22[o]+\
    bv31[o]+bv32[o]+bv41[o]+bv42[o])

    if (o % 2) == 0:
        plt.plot(lam_1_voie1[I],inte[I],'r')
    else:
        plt.plot(lam_1_voie1[I],inte[I],'b')

    print(o)

    l = np.arange(len(inte))
    nnodes=10
    q=0.3
    qq=0.7
    qqq=0.95
    p = util.continuum(l,inte,nnodes,q,qq,qqq)

    lam_1_voie1,intens_1_voie1,I = mystarname_1.get_lambda_intens1(o)

    NROWS = kwargs_moon['NROWS']

    noise = np.sqrt(inte)/p(l)
    noise=noise*mask
    res6.extend(noise)


print("im here")

plt.show()


#return col1,col2,col3,col4,col5,col6
with open(str(starname)+"_pol3.s","w") as f:
    f.write("Polarisationspektrum of "  + starname +"\n")
    f.write(str(len(res1)) + "   5\n")
    np.savetxt(f,np.column_stack([res1,res2,res3,res4,res5,res6]))

with open(str(starname)+"_int.s","w") as f:
    f.write("Intensityspektrum of " + starname + "\n")
    f.write(str(len(res1)) + "   2\n")
    np.savetxt(f,np.column_stack([res1,res2,res6]))




