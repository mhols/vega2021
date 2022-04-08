import json
from util import *
import sys
import matplotlib.pyplot as plt


data1=None
data2=None

class PolySets:

    def __init__(self):
        n=5
        o=7
        crossmax=6

        res=[]


        count = 0
        for i in range (n+1):
            for j in range (o+1):
                count+=1
                if (i+j <= crossmax) | (i==0) | (j==0):
                    #print(count,i,j)
                    res.append((i,j))
        self._monome = res
        #print(res)
    

    def func(self, i, ol, o):
        n,m = self._monome[i]
        return ol**n*o**m
        
# hier nur Ableitung nach lambda (ol), da o bekannt sind
    def dxfunc(self,i, ol, o):
        n,m = self._monome[i]
        return n*ol**(n-1)*o**m


def prepare_jsons():
    global data1, data2
    with open("./ThAr2D_voie_1.dat.json", "r") as f:
        data1 = json.load(f)

    with open("./ThAr2D_voie_2.dat.json", "r") as f:
        data2 = json.load(f)

    res = []
    shit = 0
    for l in data1:
        try:
            p = l['pixels_extract'][0] +\
                estimate_location(
                    np.array(l['flux_values_extract']), 
                    loss_1, # the loss function
                    gauss   # the flux model
                )
        except Exception as ex:
            print(ex)
            shit += 1
            plt.figure()
            plt.plot(l['flux_values_extract'])
            plt.show()
            continue

        l['new_mean_pixel_pos'] = p[1]
        res.append(l)

    with open('ThAr2D_voie_1_new.json', 'w') as f:
        json.dump(res, f)

    print ('n negative ', shit)

def estimate_polynome():
    with open("./ThAr2D_voie_1_new.json", "r") as f:
        data1 = json.load(f)

    



def plot_1():
    with open("./ThAr2D_voie_1_new.json", "r") as f:
        data1 = json.load(f)
    plt.figure()
    old = np.array([l['meas_pixel_position'] for l in data1] )
    new = np.array([l['new_mean_pixel_pos'] for l in data1])
    print (np.where(np.abs(new-old) > 0.5),)
    plt.plot(new-old, 'ro')


def plot_2(n):
    with open("./ThAr2D_voie_1_new.json", "r") as f:
        data1 = json.load(f)
    l = data1[n]
    intens = np.array(l['flux_values_extract'])
    pix = np.array(l['pixels_extract'])
    params_1 = estimate_location(
                    intens, 
                    loss_2, 
                    igauss
                )
    params_2 = estimate_location(
                    intens, 
                    loss_2, 
                    gauss
                )
    params_1[1] += pix[0]
    params_2[1] += pix[0]
                
    plt.figure()
    plt.ylim([0, 1.2*np.max(intens)])
    plt.plot (pix, intens, label='intens')
    plt.plot(pix, igauss(pix, *params_1), label='igauss')
    plt.plot(pix, gauss(pix, *params_2), label='gauss')
    plt.legend()
if __name__=='__main__':
    """    prepare_jsons()
    plot_1()
    #plot_2(int(sys.argv[1]))
    plt.show()
    """
    P = PolySets()
    x=P.func(5,30*5000.,30)
    print(x)
    #P._monome[5]qq
