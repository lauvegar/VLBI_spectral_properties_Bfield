import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import os, h5py
import scipy.stats as stats


def gauss(x,A, x0, sigma ):
	return A*np.exp(-8*np.log(2)*(x-x0)**2./(2.*sigma**2)) #the sigma here is the fwhm

def errorsPixel(int rowsNumber, int columnsNumber):

	cdef int mapsNumber = 10**4
	cdef int n, rowNum, a, columnNum,j
	up = np.zeros((columnsNumber,rowsNumber)) 
	down = np.zeros((columnsNumber,rowsNumber)) 

	for rowNum in xrange(0,rowsNumber):
		filespixel = []
		for filename in sorted(glob.glob('*pixel'+str(rowNum)+'.h5')):
			filespixel.append(filename)

		n = len(filespixel)
		print rowNum
		rows = []
		for j in xrange(0,n):
			h5f = h5py.File(filespixel[j])
			s = h5f['data'][:]
			h5f.close()
			if len(rows) > 1:
				rows=np.concatenate((rows,s),axis=0)
			else:
				rows = s

				
		for colNum in xrange(0,columnsNumber):
		# for a in xrange(0,mapsNumber):
		# # pixelValues contiene todos los pixels de la columna de row
		# pixelValues.append(rows[a][0][colNum])
			pixelValues = []
			for a in xrange(0,mapsNumber):
				pixelValues.append(rows[a][0][colNum]) #si mapsNumber son todos basta con [:][0][colNum]

			print 'lenght',len(pixelValues)

			# Si hay algun nan en el pixel [i,b] de rows NO entra
			if np.isnan(pixelValues).any() == False and np.isinf(pixelValues).any() == False:
				mean = np.mean(pixelValues)
				std = np.std(pixelValues)
				lim1 = mean -5*std
				lim2 = mean +5*std
				nn, bins,patches = plt.hist(pixelValues,bins=101,normed=True)
				xs = np.arange(lim1,lim2,0.01)

				print stats.normaltest(pixelValues)

				plt.plot(xs,gauss(xs,nn.max(),mean,2.355*std),color='black',ls='-',linewidth=2)
				plt.xlim(lim1,lim2)
				median = 0.
				pixelValues = np.sort(pixelValues)
				median = pixelValues[len(pixelValues)/2]
				upp = pixelValues[int(len(pixelValues)*0.84)]
				downn = pixelValues[int(len(pixelValues)*0.16)]

				up[(columnsNumber-1-colNum),rowNum] = abs(median - upp)
				down[(columnsNumber-1-colNum),rowNum] = abs(median - downn)
				plt.axvline(median,c='r',ls='-',linewidth = 2)
				plt.axvline(upp,c='r',ls='--',linewidth = 2)
				plt.axvline(downn,c='r',ls='--',linewidth = 2)
				plt.ylabel('counts')
				plt.xlabel(r'$\alpha$')
				print colNum

				plt.savefig('FIG/histogram'+str(rowNum)+'and'+str(colNum)+'.png')
				plt.close('all')

			else:
				up[(columnsNumber-1-colNum),rowNum] = 0.
				down[(columnsNumber-1-colNum),rowNum] = 0.



	return up,down
