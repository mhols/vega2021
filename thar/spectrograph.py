import json
from util import *
import sys
import matplotlib.pyplot as plt


data1=None
data2=None

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
                    loss_1, 
                    gauss
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
    #prepare_jsons()
    #plot_1()
    plot_2(int(sys.argv[1]))
    plt.show()
