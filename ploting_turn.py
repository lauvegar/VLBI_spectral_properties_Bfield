import threading, time
import warnings
import sys
import sip
import codecs
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from math import *
from functools import *
import numpy as np
import astropy.io.fits as pf
from pylab import *
import pickle
import iminuit, probfit
from scipy.optimize import curve_fit
from functions_conv import order_by_nu, read_conv_params
from functions_align import find_same_beam,beam_array,check_map_params, cuttingMAP,cross_correlation_shifts_FITS, checking_shift
from functions_turnover import cuttingTURN, synchrotron, synchrotron_v1, guesses_turnover,guesses_turnoverPoint, guesses_PL, powerLaw, powerLawPlot
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, get_ellipse_coords, Annotate
from functions_Bfield import searchNEDnoGUI, B_field, N_UeUb_sigma
import os,glob
import subprocess as sub
from astropy.nddata import Cutout2D
from correlate2d import *
#from fast_ftts import *




filefreq = 'bo/shifted15.fits'
header = take_header(filefreq,False)
cellsize = header[0]
mapdata = read_map(filefreq,False)
realDAT = mapdata[0]
ext = [mapdata[1],mapdata[2],mapdata[3],mapdata[4]]

res = open('turnoverdata.p','rb')
pick = pickle.load(res)
res.close()
nu = pick[0]
s = pick[1]#*1.4213481404770085
alpha0 = pick[2]
alphathick = pick[3]

plt.figure(1)
ax = plt.gca()
cset = plt.contour(realDAT,0.008*np.array([2.,4.,16.,64.,256.,1020.,2050.]),inline=1,colors=['grey'],aspect=1.0,extent=ext)
p1 = ax.imshow(nu,origin='bottom',extent=ext,vmin=6,vmax=12)

plt.axis('scaled')
plt.xlabel('Right Ascension [mas]')
plt.ylabel('Relative Declination [mas]')
plt.xlim(1.5,-3.1)
plt.ylim(-4,1.5)

divider = make_axes_locatable(ax)
cax1 = divider.append_axes("right", size="4%", pad="0.5%")

cb1 = plt.colorbar(p1,cax=cax1,cmap='jet')

cb1.set_label(r'$\nu$ [GHz]')

plt.savefig('turnoverNu.png', bbox_inches='tight')

plt.close()

plt.figure(1)
ax = plt.gca()
cset = plt.contour(realDAT,0.008*np.array([2.,4.,16.,64.,256.,1020.,2050.]),inline=1,colors=['grey'],aspect=1.0,extent=ext)
p1 = ax.imshow(s,origin='bottom',extent=ext)

plt.axis('scaled')
plt.xlabel('Right Ascension [mas]')
plt.ylabel('Relative Declination [mas]')
plt.xlim(1.5,-3.1)
plt.ylim(-4,1.5)

divider = make_axes_locatable(ax)
cax1 = divider.append_axes("right", size="4%", pad="0.5%")

cb1 = plt.colorbar(p1,cax=cax1,cmap='jet')

cb1.set_label(r'$S_y$ [Jy]')

plt.savefig('turnoverS.png', bbox_inches='tight')

plt.close()

plt.figure(1)
ax = plt.gca()
cset = plt.contour(realDAT,0.008*np.array([2.,4.,16.,64.,256.,1020.,2050.]),inline=1,colors=['grey'],aspect=1.0,extent=ext)
p1 = ax.imshow(alpha0,origin='bottom',extent=ext)

plt.axis('scaled')
plt.xlabel('Right Ascension [mas]')
plt.ylabel('Relative Declination [mas]')
plt.xlim(1.5,-3.1)
plt.ylim(-4,1.5)

divider = make_axes_locatable(ax)
cax1 = divider.append_axes("right", size="4%", pad="0.5%")

cb1 = plt.colorbar(p1,cax=cax1,cmap='jet')

cb1.set_label(r'$\alpha_0$')

plt.close()



shapee = np.shape(alpha0)
Ball = np.zeros(shapee)
Kall = np.zeros(shapee)

Ball[:] = np.nan
Kall[:] = np.nan


#cellsize *2 for getting the diameter
for i in xrange (0,shapee[0]):
    for j in xrange(0,shapee[1]):
	if np.isnan(s[i][j]) == False:
	    tempBK = B_field(alpha0[i][j],alphathick[i][j],16945,2.17,8.164,cellsize*2*12,16.5,nu[i][j],s[i][j])
	    Ball[i][j] = tempBK[0]
	    Kall[i][j] = tempBK[1]



Bmasked = np.ma.masked_where(Ball > 0.15, Ball)
Bmasked = np.ma.masked_where(Bmasked < 0.001, Bmasked)	


plt.figure(1)
ax = plt.gca()
cset = plt.contour(realDAT,0.008*np.array([2.,4.,16.,64.,256.,1020.,2050.]),inline=1,colors=['grey'],aspect=1.0,extent=ext)
p1 = plt.imshow(Bmasked,origin='bottom',extent=ext)
plt.axis('scaled')
plt.xlabel('Right Ascension [mas]')
plt.ylabel('Relative Declination [mas]')
plt.xlim(1.5,-3.1)
plt.ylim(-4,1.5)

divider = make_axes_locatable(ax)
cax1 = divider.append_axes("right", size="4%", pad="0.5%")

cb1 = plt.colorbar(p1,cax=cax1,cmap='jet')

cb1.set_label(r'$B$ [G]')

plt.savefig('turnoverB', bbox_inches='tight')


plt.show()







