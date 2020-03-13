import numpy as np
import matplotlib.pyplot as plt
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
import os,glob
#import string,math,sys,fileinput,glob,time
#load modules
#from pylab import *
import subprocess as sub
import re
#from plot_components import get_ellipse_coords, ellipse_axis
import urllib2
from astropy import units as u
#from astropy.coordinates import SkyCoord

#FUNCTION TO READ THE HEADER AND TAKE IMPORTANT PARAMETERS AS
#cell
#BMAJ, BMIN, BPA
#date, freq and epoch

def find_nearest(array,value):
	index = (np.abs(array-value)).argmin()
	return array[index], index

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def get_ellipse_coords(a=0.0, b=0.0, x=0.0, y=0.0, angle=0.0, k=2):
    """ Draws an ellipse using (360*k + 1) discrete points; based on pseudo code
    given at http://en.wikipedia.org/wiki/Ellipse
    k = 1 means 361 points (degree by degree)
    a = major axis distance,
    b = minor axis distance,
    x = offset along the x-axis
    y = offset along the y-axis
    angle = clockwise rotation [in degrees] of the ellipse;
        * angle=0  : the ellipse is aligned with the positive x-axis
        * angle=30 : rotated 30 degrees clockwise from positive x-axis
    """
    pts = np.zeros((360*k+1, 2))

    beta = -angle * np.pi/180.0
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)
    alpha = np.radians(np.r_[0.:360.:1j*(360*k+1)])
 
    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    
    pts[:, 0] = x + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = y + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    return pts

def ellipse_axis(x, y,s):
    x1=x-s
    x2=x+s
   
    if x1<x2:
        xaxis=np.linspace(x1,x2,50)    
    else:
        xaxis=np.linspace(x2,x1,50)           
    
    y1=y-s
    y2=y+s

    if y1<y2:
        yaxis=np.linspace(y1,y2,50)    
    else:
        yaxis=np.linspace(y2,y1,50)  

    return xaxis,yaxis

def ellipse_axis_lines(x,y,size):
	pts_arr=[]
	pt_arr=[]
	x_el_arr=[]
	x_elH_arr=[]
	y_el_arr=[]
	y_elH_arr=[]

	for i in xrange(0,len(x)):
		n = len(x[i])
		pts, pt = [], []
		x_el, y_el = [], []
		x_elH, y_elH = [], []
		for k in xrange(0,n):
			pts.append(get_ellipse_coords(a=size[i][k], b=size[i][k], x=x[i][k],y=y[i][k], angle=0))
			pt.append(get_ellipse_coords(a=0.01, b=0.01, x=x[i][k],y=y[i][k], angle=0))
			#lines axis ellipses      
			x_el.append(ellipse_axis(x=float(x[i][k]),y=float(y[i][k]),s=float(size[i][k]))[0])
			y_el.append(ellipse_axis(x=x[i][k],y=y[i][k],s=size[i][k])[1])
			x_elH.append(np.linspace(x[i][k],x[i][k],50))
			y_elH.append(np.linspace(y[i][k],y[i][k],50))

		pts_arr.append(pts)
		pt_arr.append(pt)
		x_el_arr.append(x_el)
		y_el_arr.append(y_el)
		x_elH_arr.append(x_elH)
		y_elH_arr.append(y_elH)

	return pts_arr,pt_arr,x_el_arr,y_el_arr,x_elH_arr,y_elH_arr

def read_modfile(file1,beam,errors):

	nfiles = len(file1)

	r_arr = []
	errr_arr = [] #np.array([0.]*nfiles)
	psi_arr =  []
	errpsi_arr =  []
	size_arr =  []
	errsize_arr =  []
	flux_arr =  []
	errflux_arr =  []

	ntot=0
	for k in xrange (0,nfiles):	
		with open(file1[k]) as myfile:
		    	count = sum(1 for line in myfile if line.rstrip('\n'))
	
		count = count-4

		#n = len(rms[k])

		n = count

		split_f=[]
		c=[]
		r=np.array([0.]*n)
		errr=np.array([0.]*n)
		psi=np.array([0.]*n)
		errpsi=np.array([0.]*n)
		size=np.array([0.]*n)
		errsize=np.array([0.]*n)
		tb=np.array([0.]*n)
		errtb=np.array([0.]*n)
		flux=np.array([0.]*n)
		fluxpeak = np.array([0.]*n)
		rms = np.array([0.]*n)
		errflux=np.array([0.]*n)
		lim_resol=np.array([0.]*n)
		errlim_resol=np.array([0.]*n)
		
		temp=file1[k]
		temp_file=open(temp,mode='r')
		temp_file.readline()
		temp_file.readline()
		temp_file.readline()
		temp_file.readline()
	        for i in xrange(0,n):
        		split_f = temp_file.readline().split()
        		flux[i] = (float(split_f[0][:-1]))
        		r[i] = (float(split_f[1][:-1]))
        		psi[i] = (float(split_f[2][:-1])*np.pi/180.)
        		size[i] = (float(split_f[3][:-1])/2.)
        		#tb[i] = (float(split_f[7])) 

		if errors == True:
			temp_file2=open('pos_errors.dat',mode='r')
			temp_file2.readline()
			temp_file2.readline()
			for i in xrange(0,ntot):
				temp_file2.readline()
			for i in xrange(0,n):		
				split_f = temp_file2.readline().split()
				fluxpeak[i] = (float(split_f[2][:-1]))	
				rms[i] = (float(split_f[1][:-1]))			

			for i in xrange(0,n):
				errflux[i] = rms[i]
				snr = fluxpeak[i]/rms[i]#[k][i] #change to flux_peak

				dlim = 4/np.pi*np.sqrt(np.pi*np.log(2)*beam[k]*np.log((snr)/(snr-1.))) #np.log((snr+1.)/(snr))) 4/np.pi*beam
				if size[i] > beam[k]:
					ddec=np.sqrt(size[i]**2-beam[k]**2)
				else:
					ddec=0.

				y=[dlim,ddec]
				dg=np.max(y)   
				err_size = rms[i]*dlim/fluxpeak[i]
				err_r = err_size/2.
				if r[i] > 0.:
					err_psi = np.real(math.atan(err_r*180./(np.pi*r[i])))
				else:
					err_psi = 1./5*beam[k]
				
				if err_size < 2./5.*beam[k]:	
					errsize[i] = 2./5.*beam[k]
				else:
					errsize[i] = (err_size)
				if err_r < 1./5*beam:	
					errr[i] =  1./5*beam
					if errr[i] < 1./2.*size[i]:
						errr[i] = 1./2.*size[i]
				else:
					errr[i] = (err_r)
				errpsi[i] = (err_psi)	
		elif errors == 'Done':
			print 'done'
		else:
			for i in xrange(0,n):
				errflux[i] = 0.1*flux[i]
				errr[i] = 1./5.*beam[k]
				errpsi[i] = 0.	
				errsize[i] = 2./5*beam[k]		

		r_arr.append(r)
		errr_arr.append(errr)
		psi_arr.append(psi)
		errpsi_arr.append(errpsi)
		size_arr.append(size)
		errsize_arr.append(errsize)
		flux_arr.append(flux)
		errflux_arr.append(errflux)
	
		ntot = n + ntot + 1

	return r_arr,errr_arr,psi_arr,errpsi_arr,size_arr,errsize_arr,tb,flux_arr,errflux_arr

def x_y(r,errr,psi,errpsi,errors):
	n = len(r)
	x,errx = np.array([0.]*n),np.array([0.]*n)
	y,erry = np.array([0.]*n),np.array([0.]*n)
	x_arr, errx_arr = [], []
	y_arr, erry_arr = [], []
	for i in xrange (0,n):
		x=r[i]*np.sin(psi[i])
	        y=r[i]*np.cos(psi[i])
		if errors == True:
			errx=np.sqrt((errr[i]*np.cos(psi[i]))**2+(r[i]*np.sin(psi[i])*errpsi[i])**2)
			erry=np.sqrt((errr[i]*np.sin(psi[i]))**2+(r[i]*np.cos(psi[i])*errpsi[i])**2)
		else:
			errx = errr[i]	
			erry = errr[i]

		x_arr.append(x)
		errx_arr.append(errx)
		y_arr.append(y)
		erry_arr.append(erry)

	x_arr = np.asarray(x_arr)
	errx_arr = np.asarray(errx_arr)
	y_arr = np.asarray(y_arr)
	erry_arr = np.asarray(erry_arr)

    	return x_arr,errx_arr,y_arr,erry_arr

def r_psi(x,errx,y,erry):
	n = len(r)
	r,errr = np.array([0.]*n),np.array([0.]*n)
	psi,errpsi = np.array([0.]*n),np.array([0.]*n)
	r_arr, errr_arr = [], []
	psi_arr, errpsi_arr = [], []
	for i in xrange (0,n):
		r=np.sqrt(x[i]**2+y[i]**2)
	        psi=np.atan(y[i]/x[i])
		#errr=np.sqrt((1/(2*r)*2*x[i]*errx[i])**2+(1/(2*r)*2*y[i]*erry[i])**2)
		#errpsi=np.sqrt(((y[i]/([x[i]**2+y[i])**2])*errx[i])**2+((x[i]/([x[i]**2+y[i])**2])*erry[i])**2)

		r_arr.append(r)
		#errr_arr.append(errr)
		psi_arr.append(psi)
		#errpsi_arr.append(errpsi)

    	return r_arr,psi_arr

def selectComponent(realDAT,realDAT2, first_contour, pts_arr,x_el_arr,x_elH_arr,y_elH_arr,y_el_arr,ext,freq1,freq2,x,y,numComp,orientation):

	levels = first_contour[0]*np.array([-1., 1., 1.41,2.,2.83,4.,5.66,8.,11.3,16.,
	                                22.6,32.,45.3,64.,90.5,128.,181.,256.,362.,512.,
	                                724.,1020.,1450.,2050.])
	plt.figure(10)
	plt.subplot(121)
	cset = plt.contour(realDAT, levels, inline=1,
		           colors=['grey'],
			   extent=ext, aspect=1.0
		           )
	for j in xrange(0,len(x_el_arr[0])):
		plt.plot(pts_arr[0][j][:,0], pts_arr[0][j][:,1], color='blue',linewidth=4)
		plt.plot(x_el_arr[0][j], y_elH_arr[0][j], color='blue',linewidth=4) 
		plt.plot(x_elH_arr[0][j], y_el_arr[0][j], color='blue',linewidth=4)
	plt.xlim(ext[0],ext[1])
	plt.ylim(ext[2],ext[3])	
	plt.axis('scaled')
	plt.xlabel('Right Ascension [pixels]')
	plt.ylabel('Relative Declination [pixels]')
	plt.title(str('%1.3f' %(freq1))+' GHz')

	levels = first_contour[1]*np.array([-1., 1., 1.41,2.,2.83,4.,5.66,8.,11.3,16.,
	                                22.6,32.,45.3,64.,90.5,128.,181.,256.,362.,512.,
	                                724.,1020.,1450.,2050.])

	#plt.figure(2)
	plt.subplot(122)
	cset = plt.contour(realDAT2, levels, inline=1,
		           colors=['grey'],
			   extent=ext, aspect=1.0
		           )
	for j in xrange(0,len(x_el_arr[1])):
		plt.plot(pts_arr[1][j][:,0], pts_arr[1][j][:,1], color='blue',linewidth=4)
		plt.plot(x_el_arr[1][j], y_elH_arr[1][j], color='blue',linewidth=4) 
		plt.plot(x_elH_arr[1][j], y_el_arr[1][j], color='blue',linewidth=4)
	plt.xlim(ext[0],ext[1])
	plt.ylim(ext[2],ext[3])	
	plt.axis('scaled')
	plt.xlabel('Right Ascension [pixels]')
	plt.title(str('%1.3f' %(freq2))+' GHz')

	param = ginput(4*numComp,0) 

	near_comp1 = []
	near_comp2 = []
	a = 0
	if orientation == 'h':
		for i in xrange(0,numComp):
			x_c = float(param[1+a][0])
			near_comp1.append(int(find_nearest(x[0],x_c)[1]))
			x_c = float(param[3+a][0])
			near_comp2.append(int(find_nearest(x[1],x_c)[1]))
			a = a + 4
	if orientation == 'v':
		for i in xrange(0,numComp):
			y_c = float(param[1+a][1])
			near_comp1.append(int(find_nearest(y[0],y_c)[1]))
			y_c = float(param[3+a][1])
			near_comp2.append(int(find_nearest(y[1],y_c)[1]))
			a = a + 4

	plt.show()

	return near_comp1, near_comp2


def CoreShiftCalculation(indexes,x,y,errx,erry,numComp):
	#indexes[0] low freq, indexes[1] high frequency
	#shift high freq - low freq
	if numComp == 1:
		RaShift = x[1][indexes[1][0]]-x[0][indexes[0][0]]
		DecShift = y[1][indexes[1][0]]-y[0][indexes[0][0]]

		errRaShift = np.sqrt((errx[1][indexes[1][0]])**2+(errx[0][indexes[0][0]])**2)
		errDecShift = np.sqrt((erry[1][indexes[1][0]])**2+(erry[0][indexes[0][0]])**2)

	if numComp > 1:
		#calculate all the Ra and Dec shifts and do an average
		RaShiftArr = np.asarray([0.]*numComp)
		DecShiftArr = np.asarray([0.]*numComp)
		for i in xrange(0,numComp):
			RaShiftArr[i] = x[1][indexes[1][i]]-x[0][indexes[0][i]]			
			DecShiftArr[i] = y[1][indexes[1][i]]-y[0][indexes[0][i]]

		RaShift = np.sum(RaShiftArr)/len(RaShiftArr)
		DecShift = np.sum(DecShiftArr)/len(DecShiftArr)

		if numComp < 4:
			#not enough values to do a proper dispersion, I consider the values' error as more reliable	
			errRaShiftArr = np.asarray([0.]*numComp)
			errDecShiftArr = np.asarray([0.]*numComp)	
			for i in xrange(0,numComp):
				#no square root because I need to square them later in the sum, so i avoid unnecessary calculations
				errRaShiftArr[i] = (errx[1][indexes[1][i]])**2+(errx[0][indexes[0][i]])**2			
				errDecShiftArr[i] = (erry[1][indexes[1][i]])**2+(erry[0][indexes[0][i]])**2

			errRaShift = np.sqrt(np.sum(errRaShiftArr))/numComp
			errDecShift = np.sqrt(np.sum(errDecShiftArr))/numComp
		else:
			#statistical error
			errRaShift = np.sqrt(np.sum((RaShiftArr-RaShift)**2))/(np.sqrt(numComp-1))
			errDecShift = np.sqrt(np.sum((DecShiftArr-DecShift)**2))/(np.sqrt(numComp-1))	

	return RaShift, DecShift, errRaShift, errDecShift








