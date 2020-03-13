from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pylab import *
#import pyspeckit as ps
from scipy import io
from scipy import stats
from scipy.optimize import leastsq
#from lmfit import minimize, Parameters, Parameter, report_fit
#from lmfit.models import GaussianModel
import scipy.special as spc
import scipy.optimize as optimization
import matplotlib.ticker as ticker
import cmath as math
import pickle
import iminuit
import astropy.io.fits as pf
import os
import subprocess as sub
from functions2 import take_header, read_map
from functions2 import convolve_difmap, Annotate
from astropy.nddata import Cutout2D
from correlate2d import *
import warnings

def cuttingTURN(image_concat,final_data,tmp_for_plot,cells):

	#because imshow only show images with 3 or 4 dimensions (each image is counted as one dimension)
	#check dimensions and plot up to a maximum of four to select the zone to cut
	if len(image_concat) == 3 or len(image_concat) == 4: 
		#plotting the maps
		plt.figure(1)
		plt.imshow(np.log10(tmp_for_plot), cmap='afmhot')
		#plt.axis('scaled')
		plt.ylim(0,len(final_data[0]))
		plt.xlabel('Right Ascension [pixels]')
		plt.ylabel('Relative Declination [pixels]')
	if len(image_concat) > 4: 
		temp_image_concat = [image_concat[0],image_concat[1],image_concat[2],image_concat[3]]
     		temp_final_data = np.dstack(temp_image_concat)
		temp2 = ma.masked_less_equal(temp_final_data,0)
		plt.figure(1)
		plt.imshow(np.log10(temp2), cmap='afmhot')
		#plt.axis('scaled')
		plt.ylim(0,len(final_data[0]))
		plt.xlabel('Right Ascension [pixels]')
		plt.ylabel('Relative Declination [pixels]')
		

		
	#select with a box the part to be analyzed
	a = Annotate()
	plt.show()
		
	[limplot_x1,limplot_x2,limplot_y1,limplot_y2] = a()
                                          

	#parameters to crop the map (center position of the rectangle and size)
	centx = ((limplot_x1+limplot_x2)/2.)
	centy = ((limplot_y1+limplot_y2)/2.) 
	position = (centx,centy)
	
	height = round(np.abs(limplot_y2-limplot_y1))+1
	width = round(np.abs(limplot_x2-limplot_x1))+1
	size = (height,width)
	
	cutout_data = []
	image_sigmacut = []

	for i in xrange(0,len(image_concat)):
		#cropping	
		cutout = Cutout2D(image_concat[i],position,size)
		cutout_data.append(cutout.data)

		cutdata = np.asarray(cutout_data)

	"""plt.figure(4)
	plt.imshow(realDAT2, origin='bottom')#,extent=ext)#,extent=ext)
	cutout_v2.plot_on_original(color='white')

	"""

	plt.show()


	return cutdata, position, size, [limplot_x1,limplot_x2,limplot_y1,limplot_y2]


"""def synchrotron(x,S_m, v_m, alpha0):
    alpha1=5./2
    Zm=3./2*(np.sqrt(1-8.*alpha0/(3*alpha1))-1)
    
    #S_v=S_m*(v/v_m)**alpha_1*(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1)*(v/v_m)**(-alpha0-alpha1)/(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1))   
    spectrum = S_m*(x/v_m)**(5./2)*(1-np.exp(-Zm*(x/v_m)**(alpha0-5/2.)))/(1-np.exp(-Zm)) 
    return spectrum"""

def synchrotron(x,S_m, v_m, alpha0,alphathick):
    #alpha1=5./2
    Zm=3./2*(np.sqrt(1-8.*alpha0/(3*alphathick))-1)
    
    #S_v=S_m*(v/v_m)**alpha_1*(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1)*(v/v_m)**(-alpha0-alpha1)/(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1))   
    spectrum = S_m*(x/v_m)**(alphathick)*(1-np.exp(-Zm*(x/v_m)**(alpha0-alphathick)))/(1-np.exp(-Zm)) 
    return spectrum

def synchrotron_v1(x,S_1, v_1, alpha0,alphathick):
    #alpha1=5./2
    #Zm=3./2*(np.sqrt(1-8.*alpha0/(3*alphathick))-1)
    
    #S_v=S_m*(v/v_m)**alpha_1*(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1)*(v/v_m)**(-alpha0-alpha1)/(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1))   
    spectrum = S_1*(x/v_1)**(alphathick)*(1-np.exp(-(x/v_1)**(alpha0-alphathick)))/(1-np.exp(-1)) 
    return spectrum

def powerLaw(x,cte,alpha0):
    #alpha1=5./2
    #Zm=3./2*(np.sqrt(1-8.*alpha0/(3*alphathick))-1)
    
    #S_v=S_m*(v/v_m)**alpha_1*(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1)*(v/v_m)**(-alpha0-alpha1)/(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1))   
    y = cte + alpha0*x
    return y

def powerLawPlot(x,cte,alpha0):
    #alpha1=5./2
    #Zm=3./2*(np.sqrt(1-8.*alpha0/(3*alphathick))-1)
    
    #S_v=S_m*(v/v_m)**alpha_1*(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1)*(v/v_m)**(-alpha0-alpha1)/(1-np.exp(-3./2*(np.sqrt(1+(8*alpha0/3*alpha1))-1))   
    spectrum = cte*x**alpha0
    return spectrum


def guesses_turnover(xdata,ydata):
	'''
	given i = slice number
	num_gauss = number of gaussians used for the fitting
	xdata = array containing the lists for each slice where the X POSITION in the map is stored
	ydata = array containing the lists for each slice where the FLUX in the map is stored
	'''

	plt.ion()

	#plot the profile for the slice	
	plt.figure(1)
	plt.plot(xdata,ydata, '.')
	plt.xlabel('log ' r'$\nu $' ' [GHz]')
	plt.ylabel('log f [Jy]')
        plt.xscale('log')
        plt.yscale('log')
	print("Please, click on the profile.")
	print("The first click will take the guess value of the turnover frequency and flux.")
	print("The first and second will take two values of flux and frequency of the syncrotron part to obtain the slope.")


	#set the number of clics needed in total and store them in an array called param
	param = ginput(2,0) 

	#store the input parameters given in the clics in the correspondent guess parameter type 	
	guess_Sm= np.abs(float(param[0][1]))	
	guess_vm = float(param[0][0])
	guess_alpha0 = (float(param[0][1])-float(param[1][1]))/(float(param[0][0])-float(param[1][0]))

	print("These are your guessing parameters:")
	print("Sm: ", guess_Sm)
	print("vm: ", guess_vm)
	print("alpha0: ", guess_alpha0)



	return guess_Sm,guess_vm,guess_alpha0

def guesses_turnoverPoint(xdata,ydata):
	'''
	given i = slice number
	num_gauss = number of gaussians used for the fitting
	xdata = array containing the lists for each slice where the X POSITION in the map is stored
	ydata = array containing the lists for each slice where the FLUX in the map is stored
	'''

	plt.ion()

	#plot the profile for the slice	
	plt.figure(1)
	plt.plot(xdata,ydata, '.')
	plt.xlabel('log ' r'$\nu $' ' [GHz]')
	plt.ylabel('log f [Jy]')
        plt.xscale('log')
        plt.yscale('log')
	print("Please, click on the profile.")
	print("The first click will take the guess value of the turnover flux.")
	print("The second click will take the guess value of the turnover frequency.")
	print("The third and four will take two values of flux and frequency of the syncrotron part to obtain the slope.")


	#set the number of clics needed in total and store them in an array called param
	param = ginput(1,0) 

	#store the input parameters given in the clics in the correspondent guess parameter type 	
	guess_Sm= np.abs(float(param[0][1]))	
	guess_vm = float(param[0][0])


	print("These are your guessing parameters:")
	print("Sm: ", guess_Sm)
	print("vm: ", guess_vm)



	return guess_Sm,guess_vm

def guesses_PL(xdata,ydata):
	'''
	given i = slice number
	num_gauss = number of gaussians used for the fitting
	xdata = array containing the lists for each slice where the X POSITION in the map is stored
	ydata = array containing the lists for each slice where the FLUX in the map is stored
	'''

	plt.ion()

	#plot the profile for the slice	
	plt.figure(1)
	plt.plot(xdata,ydata, '.')
	plt.xlabel('log ' r'$\nu $' ' [GHz]')
	plt.ylabel('log f [Jy]')
        plt.xscale('log')
        plt.yscale('log')
	print("Please, click on the profile.")
	print("The first click will take the guess value of the turnover flux.")
	print("The second click will take the guess value of the turnover frequency.")
	print("The third and four will take two values of flux and frequency of the syncrotron part to obtain the slope.")


	#set the number of clics needed in total and store them in an array called param
	param = ginput(3,0) 

	#store the input parameters given in the clics in the correspondent guess parameter type 	
	guess_cte= np.abs(float(param[0][0]))	
	#guess_vm = float(param[1][0])
	guess_alpha0 = (float(param[1][1])-float(param[2][1]))/(float(param[1][0])-float(param[2][0]))

	print("These are your guessing parameters:")
	print("cte: ", guess_cte)
	#print("vm: ", guess_vm)
	print("alpha0: ", guess_alpha0)



	return guess_cte,guess_alpha0


def B_field(alpha0,alphathick,D_l,z,scale,R,delta,vm,Sm):

	#Dl 0836 5.5 Gpc  1Jy = 10**(-23) erg

	D_l = D_l*10**6.*3.08*10**(18) #cm
	R = R*scale*3.08*10**(18) #cm
	vm = vm*10**9 #Hz
	Sm = Sm*10**(-23) #erg/(s*Hz*cm**2)

	me = 9.11*10**(-28) #g
 	e = 4.8*10**(-10) #
	c = 3*10**10 #cm/s

	s = -2*alpha0 + 1
	taum = 3./2.*(np.sqrt(1.-8*alpha0/(3*alphathick))-1)


	factor_c = np.sqrt(3)*e**3/(16*np.pi*me*c**2)
	factor_k = np.sqrt(3)*np.pi/72.*e*me**5*c**10
	factor_2 = 3*e/(2*np.pi*me**3*c**5)	
	c_eps = factor_c*factor_2**((s-1)/2)*((s+7./3)/(s+1.))*spc.gamma((3*s-1.)/12.)*spc.gamma((3*s+7.)/12.)
	c_k = factor_k*factor_2**((s+4)/2)*((s+10.)/(3.))*spc.gamma((3*s+2.)/12.)*spc.gamma((3*s+10.)/12.)
	c_eps_b = np.sqrt(np.pi)/2*spc.gamma((s+5.)/4.)*(spc.gamma((s+7.)/4.))**(-1.)
	c_k_b = np.sqrt(np.pi)/2.*spc.gamma((s+6.)/4.)*(spc.gamma((s+8.)/4.))**(-1.)
	
	B = np.pi**2/D_l**4*((c_eps*c_eps_b)/(c_k*c_k_b))**2*(1+z)**7*R**4*delta*vm**5*Sm**(-2)*taum**2

	#normalization factor
	K = (c_eps*c_eps_b)**(-(s+2))*(c_k*c_k_b)**(s+1)*(1.+z)**(-(3*s+5))*delta**(-(s+3))*taum**(-(s+1))*Sm**(s+2)*D_l**(2*s+4)/(np.pi**(s+2))*R**(-(2*s+5))*vm**(-(2*s+3))

	if np.isinf(K) == False:
		print('lll',K)

	return B,K

def N_UeUb_sigma(B,K,gammamin,gammamax,alpha0):

	me = 9.11*10**(-28) #g
 	e = 4.8*10**(-10) #
	c = 3*10**10 #cm/s

	s = -2*alpha0 + 1

	#number of particles
	N = K/(s-1.)*(me*c**2)**(1.-s)*gammamin**(1.-s)*(1.-(gammamax/gammamin)**(1.-s))

	#Ue (total energy distribution of the relativistic particles) done integrating N = KE**-s between the limits E(gammamin) and E(gammamax)
	if s > 1. and s < 2.:
		Ue = K/(2.-s)*(me*c**2)**(2.-s)*gammamax**(2.-s)*(1.-(gammamax/gammamin)**(2.-s))
	elif s > 2.:
		Ue = K/(s-2.)*(me*c**2)**(2.-s)*gammamin**(2.-s)*(1.-(gammamax/gammamin)**(2.-s))
	elif s == 2.:
		Ue = K*math.log(gammamax/gammamin)

	#magnetic energy density
	Ub = B**2/(8*np.pi)

	#ration between magnetic energy density and energy density of the relativistic particles
	sigma = Ub/Ue

	return N,Ue,Ub,sigma

import urllib2
from astropy import units as u

def searchNEDnoGUI(source_name):

	if source_name.find('+') > 0:
		splitedName = source_name.split('+')
		name1 = splitedName[0]
		name2 = splitedName[1]
		response = urllib2.urlopen('https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname='+str(name1)+'%2B'+str(name2)+'&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES')
	elif source_name.find('-') > 0:
		splitedName = source_name.split('-')
		name1 = splitedName[0]
		name2 = splitedName[1]
		response = urllib2.urlopen('https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname='+str(name1)+'-'+str(name2)+'&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES')


	else:
		response = urllib2.urlopen('https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname='+str(source_name)+'&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES')

	html = response.read()

	DLpos = html.find('Luminosity Distance')
	DLpos1st = DLpos+29
	afterDLpos = html[DLpos:]
	DLposEnd = DLpos+afterDLpos.find('Mpc')-2

	#print html[DLpos1st], html[DLposEnd]

	DL= ''
	for i in xrange (DLpos1st,DLposEnd+1):
		DL= DL + html[i]

	try:		
		DL = float(DL)
		print(DL)

		Scalepos = html.find('Scale (Cosmology Corrected)')
		Scalepos1st = Scalepos+31
		afterScalepos = html[Scalepos:]
		ScaleposEnd = Scalepos+afterScalepos.find('pc')-2

		Scale= ''
		print(Scalepos1st,ScaleposEnd+1)
		for i in xrange (Scalepos1st,ScaleposEnd+1):
			Scale= Scale + html[i]

		Scale = float(Scale)/10**3


		z1 = html.find('HREF="#ObjNo1"')+13
		z2 = html[z1:]
		z3 = z1 + z2.find('&nbsp')
		z4 = html[z1:z3].split(' ')
		z5 = []
		for item in z4:
			if item:
				z5.append(item)

		z = float(z5[len(z5)-1])

		print('DL =', DL, 'z=',z)

		ra = z5[len(z5)-5]
		dec = z5[len(z5)-4]

		raSplit1 = ra.split('h')
		raSplit2 = raSplit1[1].split('m')
		raH = float(raSplit1[0])
		raM = float(raSplit2[0])
		raS = float(raSplit2[1].split('s')[0])	
		raHours = raH + raM/60 + raS/(60*60)
		raDeg = raHours*15

		decSplit1 = dec.split('d')
		decSplit2 = decSplit1[1].split('m')
		decD = float(decSplit1[0])
		decM = float(decSplit2[0])
		decS = float(decSplit2[1].split('s')[0])	
		if source_name.find('+') > 0:
			decDeg = decD + decM/(60) + decS/(60*60)
		elif source_name.find('-') > 0:
			decDeg = decD - decM/(60) - decS/(60*60)
		elif decD < 0:
			decDeg = decD - decM/(60) - decS/(60*60)

		#c = SkyCoord(ra=raHours*u.degree,dec=decDays*u.degree)
	
		print('RA =', ra, 'DEC =', dec)
		print('RA =', raDeg, 'degrees, DEC =', decDeg, 'degrees' )

	except ValueError:
		DL = 0.
		z = 0.
		Scale = 0.

	return DL,z, Scale



"""filefreq = 'bo/shifted15.fits'
header = take_header(filefreq)
mapdata = read_map(filefreq)
realDAT = mapdata[0]
ext = [mapdata[1],mapdata[2],mapdata[3],mapdata[4]]

res = open('turnoverdata.p','rb')
pick = pickle.load(res)
res.close()
vmall = pick[0]
small = pick[1]
alpha0all = pick[2]
shapee = np.shape(alpha0all)
ball = np.zeros(shapee)


ball[:] = np.nan


for i in xrange (0,shapee[0]):
    for j in xrange(0,shapee[1]):
        if np.isnan(small[i][j]) == False:
            ball[i][j] = B_field(alpha0all[i][j],2.5,16945,2.17,0.01*12,16.5,vmall[i][j],small[i][j])

plt.figure(1)
cset = plt.contour(realDAT,0.008*np.array([2.,4.,16.,64.,256.,1020.,2050.]),inline=1,colors=['grey'],aspect=1.0,extent=ext)
plt.imshow(ball,origin='bottom',extent=ext,vmax=0.09)
plt.colorbar()
plt.show()"""


