#from settingspegasus import kwargs
from settingsmoon import kwargs
from numpy import *
import matplotlib.pyplot as plt
from extract import *


myext = get_ext('/home/hols/vega2021/thar/vega_reference/NEO_20220903_191404_th0.fits')

o = 36
s = myext.snippets_voie[1]
I = s['true_order_number'] == o
x = s.loc[I, 'est_lambda'].to_numpy()
a= s.loc[I, "lambda_left"].to_numpy()
b= s.loc[I, "lambda_right"].to_numpy()

xx= myext._atlas.lambda_ref()
aa= myext._atlas.lmin()
bb= myext._atlas.lmax()

M = (aa[None,:] <= x[:,None]) & (x[:,None] <= bb[None,:])
N = (a[:,None] <= xx[None, :]) & (xx[None,:] <= b[:, None])

print( where(M & N) )

