import sys
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from math import *
from functools import *
import numpy as np
import math
import astropy.io.fits as pf
import pickle
#from functions_conv import order_by_nu, read_conv_params
#from functions_align import find_same_beam,beam_array,check_map_params, cuttingMAP,cross_correlation_shifts_FITS, checking_shift, search_rms
from functions2 import take_header, read_map, saver
#from functions2 import convolve_difmap, Annotate
import glob
import subprocess as sub
#from astropy.nddata import Cutout2D
#from correlate2d import *
from scipy.fftpack import fftfreq
#from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial
import os, h5py

from collections import Mapping, Container
import pyximport; pyximport.install()
pyximport.install(pyimport = True)
import pixelerror

#from memory_profiler import profile
 
#@profile
 
def deep_getsizeof(o, ids):
    """Find the memory footprint of a Python object
 
    This is a recursive function that drills down a Python object graph
    like a dictionary holding nested dictionaries with lists of lists
    and tuples and sets.
 
    The sys.getsizeof function does a shallow size of only. It counts each
    object inside a container as pointer only regardless of how big it
    really is.
 
    :param o: the object
    :param ids:
    :return:
    """
    d = deep_getsizeof
    if id(o) in ids:
        return 0
 
    r = sys.getsizeof(o)
    ids.add(id(o))
 
    if isinstance(o, str) or isinstance(0, unicode):
        return r
 
    if isinstance(o, Mapping):
        return r + sum(d(k, ids) + d(v, ids) for k, v in o.iteritems())
 
    if isinstance(o, Container):
        return r + sum(d(x, ids) for x in o)
 
    return r 

"""def shift_fft(realDAT,shift):
    shift_rows,shift_cols = shift
    print shift_rows
    print shift_cols
    nr,nc = realDAT.shape
    Nr, Nc = fftfreq(nr), fftfreq(nc)
    Nc,Nr = np.meshgrid(Nc,Nr)
    fft_inputarray = np.fft.fft2(realDAT)
    fourier_shift = np.exp(1j*2*np.pi*((shift_rows*Nr)+(shift_cols*Nc)))
    output_array = np.fft.ifft2(fft_inputarray*fourier_shift)

    return np.real(output_array)

def shiftArr(x,y):
    shiftsArr = [x,y]
    return shiftsArr


class hallo():

    	def __init__(self,*args):

		self.fits1 = 'test22.fits'
		self.fits2 = 'test43.fits'

		header1 = take_header(self.fits1)
		header2 = take_header(self.fits2)
		map_data1 = read_map(self.fits1)		
		self.realDAT = map_data1[0]
		map_data2 = read_map(self.fits2)		
		self.realDAT2 = map_data2[0]

		#obtaining the beam, cell and mapsize from the header
		self.bmaj_files = header1[1]
		self.bmin_files = header1[2]
		self.bpa_files = header1[3]
		self.beam_files = header1[7]
		self.cells_files = header1[0]
		self.mapsize_files = 2*map_data1[7]

		#obtaining frequencies from the header
		self.freq1 = header1[5]
		self.freq2 = header2[5]

		#obtaining map centers in pixels from the header
		self.cent_mapx = map_data1[5]
		self.cent_mapy = map_data1[6]

		#obtaining the four corners of the maps in mas from the header
		x1 = map_data1[1]
		x2 = map_data1[2]
		y1 = map_data1[3]
		y2 = map_data1[4]

		header1 = []
		header2 = []
		map_data1 = []
		map_data2 = []

		self.rms1 = 0.0018
		self.rms2 = 0.0068

		self.shiftsX = np.random.normal(-5.73,1,10**4)

		self.shiftsY = np.random.normal(0.2,1,10**4)

		self.shiftALL = map(shiftArr,self.shiftsX,self.shiftsY)

		a = 0
		self.shiftPartly = []
		for i in xrange(0,10**4+100,100):
			self.shiftPartly.append(self.shiftALL[a:i])
			a = i

		print len(self.shiftPartly)
		self.shiftPartly = self.shiftPartly[0:4]
		#del self.shiftALL

		#shiftedImage = shift_fft(self.realDAT,[0.2,-5.73])
		#image1 = shiftedImage*(shiftedImage > 0.0018)
		#image2 = self.realDAT2*(self.realDAT2 > 0.0068)
		#a = np.log10(image1/image2)/np.log10(self.freq1/self.freq2)

		#plt.imshow(a,origin='bottom')
		#print a
		#plt.show()

		aALL = []

		self.i = 0
		self.i2 = 0

	
		pool = ThreadPool(1)
		self.x = pool.map(self.shift_fft,self.shiftALL)
		#self.x = self.shift_fft(self.shiftALL)
		#x = pool.map(self.hola,self.shiftPartly)
		for i in xrange(0,len(self.x)):
			temp1 = self.x[i].next()
			#if math.isnan(temp1[500][500]) == False:
			print temp1[1024][1024]
		#pool.close()
		#pool.join()



        def mean(self,data):

		n = 0
		mean = 0.0

		for x in data:
			n += 1
			mean += (x-mean)/n

		return mean 
		


	def hola(self,shiftPartly):

		pool2 = ThreadPool(4)

		#print self.realDAT
		#for i in xrange(0,len(shiftPartly)):
		#	 pool = ThreadPool(4)
		print 'i',self.i
		a = pool2.map(self.shift_fft,shiftPartly)
		#pool.close()
		#pool.join()
		#res=open('pickle'+str(self.i)+'.p','wb')
		#pickle.dump(a,res)
		#res.close() 

		#h = h5py.File(str(self.i)+'.h5', 'w')
		#h.create_dataset('data', data=a)
		#h.close()
		#del a


		self.i = self.i +1   


	def shift_fft(self,shift):

	    #i = 0
	    #while i < len(shift): 
	    #print self.realDAT
	    shift_rows,shift_cols = shift
	    #print shift
	    nr,nc = self.realDAT.shape
	    Nr, Nc = fftfreq(nr), fftfreq(nc)
	    Nc,Nr = np.meshgrid(Nc,Nr)
	    fft_inputarray = np.fft.fft2(self.realDAT)
	    fourier_shift = np.exp(1j*2*np.pi*((shift_rows*Nr)+(shift_cols*Nc)))
	    output_array = np.fft.ifft2(fft_inputarray*fourier_shift)

	    shiftedImage = np.real(output_array)

	    image1 = shiftedImage*(shiftedImage > self.rms1)
	    image2 = self.realDAT2*(self.realDAT2 > self.rms2)

	    print self.i2

	    a = np.log10(image1/image2)/np.log10(self.freq1/self.freq2)

	    print 'a', deep_getsizeof(a, set())

	    self.i2 = self.i2 +1   
	    del image1
	    del image2
	    del shiftedImage
	    del output_array

	    yield a

	    del a
	    #i += 1

	def shiftArr(self):
	    shiftsArr = [y,x]
	    return shiftsArr


hallo()"""


import scipy.stats as stats

	
#n, bins,patches = plt.hist(example,100)

#l = plt.plot(bins)

#plt.show()


pixelerror.errorsPixel(758,684)


 
