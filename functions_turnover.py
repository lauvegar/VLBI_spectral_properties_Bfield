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

def cuttingTURN(image_concat,final_data,tmp_for_plot,cells,turnoverNu,turnoverNuupper,checkedsync):

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
		
	if checkedsync == True:
		p1 = plt.imshow(turnoverNu, origin='bottom')#,extent=self.ext)#, vmin=-2.5, vmax=1.7)
		p2 = plt.imshow(turnoverNuupper, origin='bottom',cmap='spring')#,extent=self.ext)#, vmin=-2.5, vmax=1.7)

		
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


def B_field(s,D_l,z,R,delta,vm,Sm,taum):

	#Dl 0836 5.5 Gpc

	me = 9.11*10**(-28) #g
 	e = 4.8*10**(-10) #
	c = 3*10**10 #cm/s
	factor_c = np.sqrt(3)*e**3/(16*np.pi*me*c**2)
	factor_k = np.sqrt(3)*np.pi/72.*e*me**5*c**10
	factor_2 = 3*e/(2*np.pi*me**3*c**5)	
	c_eps = factor_c*factor_2**((s-1)/2)*((s+7./3)/(s+1.))*gamma((3*s-1.)/12.)*gamma((3*s+7.)/12.)
	c_k = factor_k*factor_2**((s+4)/2)*((s+10.)/(3.))*gamma((3*s+2.)/12.)*gamma((3*s+10.)/12.)
	c_eps_b = np.sqrt(np.pi)/2*gamma((s+5.)/4.)*(gamma((s+7.)/4.))**(-1.)
	c_k_b = np.sqrt(np.pi)/2.*gamma((s+6.)/4.)*(gamma((s+8.)/4.))**(-1.)
	
	B = np.pi**2/D_l**4*((c_eps*c_eps_b)/(c_k*c_k_b))**2*(1+z)**7*R**4*delta*vm**5*Sm**(-2)*taum**2

	return B


