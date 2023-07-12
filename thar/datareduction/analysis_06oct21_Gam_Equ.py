from extract import *
from settings_reference import kwargs as kwargs_reference
import pickle
import util
from units import *
import matplotlib.pyplot as plt

fitsfile_reference = os.path.join(kwargs_reference['BASEDIR'], 'vega_reference/NEO_20220903_191404_th0.fits')

starname = 'gammaequ'
fitsfile_starname = os.path.join(kwargs_reference['BASEDIR'], '06oct21_Gam_Equ/NEO_20211006_174327_th0.fits')
starname_1 = os.path.join(kwargs_reference['BASEDIR'], '06oct21_Gam_Equ/NEO_20211006_214054_st0.fits')
starname_2 = os.path.join(kwargs_reference['BASEDIR'], '06oct21_Gam_Equ/NEO_20211006_214255_st0.fits')
starname_3 = os.path.join(kwargs_reference['BASEDIR'], '06oct21_Gam_Equ/NEO_20211006_214454_st0.fits')
starname_4 = os.path.join(kwargs_reference['BASEDIR'], '06oct21_Gam_Equ/NEO_20211006_214656_st0.fits')

# the following lines may be commented out when done..

"""
mystarname = Extractor(fitsfile_starname, **kwargs_reference) # TODO starname need own setting
mystarname.update()
mystarname.update()
mystarname.update()

mystarname.save_to_store()
"""




# mystarname is now ready to reduce any star
mystarname_1 = get_ext(fitsfile_starname)
mystarname_1.update_kwargs(VOIE_METHOD ='SUM_DIVIDE')

# reduce the first starname sprectrum
mystarname_1.set_fitsfile(starname_1)

# mystarname is now ready to reduce any star
mystarname_2 = get_ext(fitsfile_starname)
mystarname_2.update_kwargs(VOIE_METHOD ='SUM_DIVIDE')

# reduce the first starname sprectrum
mystarname_2.set_fitsfile(starname_2)

# mystarname is now ready to reduce any star
mystarname_3 = get_ext(fitsfile_starname)
mystarname_3.update_kwargs(VOIE_METHOD ='SUM_DIVIDE')

# reduce the first starname sprectrum
mystarname_3.set_fitsfile(starname_3)

# mystarname is now ready to reduce any star
mystarname_4 = get_ext(fitsfile_starname)
mystarname_4.update_kwargs(VOIE_METHOD ='SUM_DIVIDE')

# reduce the first starname sprectrum
mystarname_4.set_fitsfile(starname_4)

mystarname_1.voie
mystarname_2.voie
mystarname_3.voie
mystarname_4.voie

mystarname_1.save_to_store()
mystarname_2.save_to_store()
mystarname_3.save_to_store()
mystarname_4.save_to_store()




mystarname_1 = get_ext(starname_1)
mystarname_2 = get_ext(starname_2)
mystarname_3 = get_ext(starname_3)
mystarname_4 = get_ext(starname_4)

myhds = [mystarname_1, mystarname_2, mystarname_3, mystarname_4]


lam_voie1={}
lam_voie2={}
intens_voie1={}
intens_voie2={}
I = {}
intens_voie2_interp = {}
weight_voie1={}
weight_voie2={}
ff_voie1={}
ff_voie2={}

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
        
for o in mystarname_1.ORDERS:
    ff_voie1[o] = mystarname_1.beams[o].beam_sum_voie1(mystarname_1.masterflat)
    ff_voie2[o] = mystarname_1.beams[o].beam_sum_voie2(mystarname_1.masterflat)

"""
for mh in myhds:
    mh.save_to_store()
"""

"""
we now have z.B. lam_voie1[3][23]
"""
bv11 = mystarname_1.bare_voie1
bv12 = mystarname_1.bare_voie2
bv21 = mystarname_2.bare_voie1
bv22 = mystarname_2.bare_voie2
bv31 = mystarname_3.bare_voie1
bv32 = mystarname_3.bare_voie2
bv41 = mystarname_4.bare_voie1
bv42 = mystarname_4.bare_voie2

bf11 = mystarname_1.beams[o].beam_sum_voie1(mystarname_1.masterflat)


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
    
#def outputfile():

res1 = []
res2 = []
res3 = []
res4 = []
res5 = []
res6 = []
int1 = []
ff1 = []
mask = []
masktotal = []
NROWS = mystarname_1.NROWS

for o in mystarname_1.ORDERS[::-1]:
#lambda
    lam_1_voie1,intens_1_voie1,I = mystarname_1.get_lambda_intens1(o)
    res1.extend(lam_1_voie1/10.)
    mask = []
    for i in range(NROWS):
        if i in I:
            mask.append(1)
        else:
            mask.append(0)
    masktotal.extend(mask)

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
    `
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
    inte1 = bv11[o]+bv21[o]+bv31[o]+bv41[o]
    inte2 = bv12[o]+bv22[o]+bv32[o]+bv42[o]
    
    
    inte = inte1 + inte2
    
    int1.extend(inte1)
    ff1.extend(ff_voie1[o])
    """
    if (o % 2) == 0:
        plt.plot(lam_1_voie1[I],inte1[I]/ff_voie1[o][I],'r')
    else:
        plt.plot(lam_1_voie1[I],inte1[I]/ff_voie1[o][I],'b')
    """
    
    if (o % 2) == 0:
        plt.plot(lam_1_voie1[I],inte1[I],'r')
        plt.plot(lam_1_voie1[I],ff_voie1[o][I],'r')
    else:
        plt.plot(lam_1_voie1[I],inte1[I],'b')
        plt.plot(lam_1_voie1[I],ff_voie1[o][I],'b')
    
    
    print(o)
    
    l = np.arange(len(inte))
    nnodes=10
    q=0.3
    qq=0.7
    qqq=0.95
    p = util.continuum(l,inte,nnodes,q,qq,qqq)
    
    lam_1_voie1,intens_1_voie1,I = mystarname_1.get_lambda_intens1(o)

    NROWS = kwargs_reference['NROWS']
    
    noise = np.sqrt(inte)/p(l)
    noise=noise*mask
    res6.extend(noise)


print("im here")

plt.show()


#return col1,col2,col3,col4,col5,col6
with open(str(starname)+".s","w") as f:
    f.write("dies ist ein Versuch eines polarisationsspektrums von starname\n")
    f.write(str(len(res1)) + "   5\n")
    np.savetxt(f,np.column_stack([res1,res2,res3,res4,res5,res6]))
    
    
with open("adrian.s","w") as f:
    f.write("spectre test de beta crB pour la concatenation des ordres ...stage Adrian\n")
    f.write(str(len(res1)) + "   2\n")
    np.savetxt(f,np.column_stack([res1,int1,ff1,masktotal]))




