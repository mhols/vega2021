from extract import *
from settings_reference import kwargs as kwargs_moon
import pickle
import util
from units import *
import matplotlib.pyplot as plt
import sys

deltav = 6. #km/s
clum = 3e5 #km/s

lis=np.loadtxt('../reffiles/sun_mask.5750.45.p00.correct')
cata=10*lis[:,0]


f=open('../reffiles/Selectedlines_fit.pkl3','rb')
c=pickle.load(f)
tamp=[]
tamp=list(set(c['sellines']))
lines=np.sort(np.array(tamp))
aa = lines*(1-deltav/clum)
bb = lines*(1+deltav/clum)

M = (aa[:,None] <= cata[None,:]) & (bb[:,None] >= cata[None,:])
I,J = np.where(M)
#I.shape


sun=np.loadtxt('moon_06apr_int.s',skiprows=2)
sunlam=sun[:,0]
sunint=sun[:,1]


plt.plot(sunlam*10.,sunint)
plt.vlines(util.air_to_vac(cata)*cata,0,1,"y")
plt.vlines(util.air_to_vac(lines)*lines,0,0.5,"r")
plt.vlines(util.air_to_vac(lines[I])*lines[I],0,0.25,"b")
plt.show()


f=open("selectedlines_match_solarkurucz.pkl","wb")
lamselect=lines[I].tolist()
donnees={}
donnees['sellines'] = lamselect
pickle.dump(donnees,f)
f.close()

