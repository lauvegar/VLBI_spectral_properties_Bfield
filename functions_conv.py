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
#import iminuit
import astropy.io.fits as pf
import os
import subprocess as sub
from functions2 import take_header, read_map

#create several arrays of map parameters in order (low freq first, high frequency last)

def order_by_nu(files,models,fits,ifhdu):

	#initialize arrays
	cell = np.array([0.]*len(fits))
	bmaj = np.array([0.]*len(fits))
	bmin = np.array([0.]*len(fits))
	bpa = np.array([0.]*len(fits))
	freq = np.array([0.]*len(fits))
	beam = np.array([0.]*len(fits))
	size_map = np.array([0.]*len(fits))
	size_map_y =np.array([0.]*len(fits))
	
	#read some parameters of the map
	for i in xrange(0,len(fits)):
		header = take_header(fits[i],ifhdu)
		cell[i] = header[0]
		bmaj[i] = header[1]
		bmin[i] = header[2]
		bpa[i] = header[3]
		freq[i] = header[5]
		beam[i] = header[7]
		map_details = read_map(fits[i],ifhdu)
		if map_details[7] == map_details[8]:
			size_map[i] = map_details[7]
		else:
			size_map[i] = map_details[7]
			size_map_y[i] = map_details[8]	

	#store them in ascendent order related with the frequency
	n=0
	for z in freq:
		n+=1
	for i in range(n):
		for j in range(1,n-i):
			if freq[j-1] > freq[j]:
				(freq[j-1],freq[j]) = (freq[j],freq[j-1])
				(cell[j-1],cell[j]) = (cell[j],cell[j-1])
				(bmaj[j-1],bmaj[j]) = (bmaj[j],bmaj[j-1])
				(bmin[j-1],bmin[j]) = (bmin[j],bmin[j-1])
				(bpa[j-1],bpa[j]) = (bpa[j],bpa[j-1])
				(beam[j-1],beam[j]) = (beam[j],beam[j-1])
				(size_map[j-1],size_map[j]) = (size_map[j],size_map[j-1])
				(size_map_y[j-1],size_map_y[j]) = (size_map_y[j],size_map_y[j-1])
				(files[j-1],files[j]) = (files[j],files[j-1])
				(models[j-1],models[j]) = (models[j],models[j-1])
				(fits[j-1],fits[j]) = (fits[j],fits[j-1])

	return freq,cell,bmaj,bmin,bpa,beam,size_map,size_map_y,files,models,fits

#read the convolution parameters stored in a pickle file
def read_conv_params(freq,positions):

	freq_selected = []

	for i in xrange(0,len(positions)):
		index = positions[i]
		freq_selected.append(int(round(freq[index])))

	if len(positions) == 2:
		res=open('convolve_param'+str(freq_selected[0])+'and'+str(freq_selected[1])+'.p','rb')
		pick = pickle.load(res)
		res.close() 
	else:
		index = len(freq_selected)-1
		res=open('convolve_param'+str(freq_selected[0])+'to'+str(freq_selected[index])+'.p','rb')
		pick = pickle.load(res)
		res.close() 

	freq2conv = pick[0]
	cell2conv = pick[1]
	bmaj2conv = pick[2]
	bmin2conv = pick[3]
	bpa2conv = pick[4]
	size_map2conv = pick[5]
	size_map_y2conv = pick[6]
	files2conv = pick[7]
	models2conv = pick[8]	
	files_out = pick[9]
	mapsize = pick[10]
	pixelsize = pick[11]

	print 'The files have'
	print 'bmaj = ', bmaj2conv, 'bmin = ', bmin2conv, 'bpa = ', bpa2conv
	print 'Mapsize = ', mapsize, 'cell = ', pixelsize


	return freq2conv,cell2conv,bmaj2conv,bmin2conv,bpa2conv,size_map2conv,size_map_y2conv,files2conv,models2conv,files_out,mapsize,pixelsize 

#check if the cell and mapsize of the images are the same
def check_map_params(image1,image2,ifhdu):
	"""v1 < v2"""

	header1 = take_header(image1,ifhdu)
	header2 = take_header(image2,ifhdu)

	cell1 = header1[0]
	cell2 = header2[0]

	if cell1 == cell2:
		cells = header1[0]
		print 'Same cellsize'
	else:
		print 'The images do not have the same cellsize'
		print 'Convolve them using the same cell'
		raise Exception('Stopping!')

	map_data1 = read_map(image1,ifhdu)		
	realDAT = map_data1[0]
	map_data2 = read_map(image2,ifhdu)		
	realDAT2 = map_data2[0]

	if realDAT.shape == realDAT2.shape:
		print 'Same mapsize' 
	else:
		print 'The images do not have the same mapsize'
		print 'Convolve them using the same mapsize'
		raise Exception('Stopping!')

	beam1 = header1[7]
	beam2 = header2[7]

	if beam1 == beam2:
		beam = header1[7]
		print 'Same beam'
	else:
		print 'The images do not have the same beam'
		print 'Convolve them using the same beam'
		raise Exception('Stopping!')

	#obtaining frequencies
	freq1 = header1[5]
	freq2 = header2[5]

	#obtaining map centers in pixels
	cent_mapx = map_data1[5]
	cent_mapy = map_data1[6]

	#obtaining the four corners of the maps in mas
	x1 = map_data1[1]
	x2 = map_data1[2]
	y1 = map_data1[3]
	y2 = map_data1[4]


	print 'Proceeding with the next step'

	return realDAT,realDAT2, freq1,freq2, cells, beam, cent_mapx, cent_mapy, x1,x2,y1,y2






