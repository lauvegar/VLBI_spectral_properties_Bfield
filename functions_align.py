#crosscorrelation_shits is part of the code of the image_registration repository

"""Copyright (c) 2012 Adam Ginsburg

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

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
#import iminuit
import astropy.io.fits as pf
import os
import subprocess as sub
from functions2 import take_header, read_map
from functions2 import convolve_difmap, Annotate
from astropy.nddata import Cutout2D
from correlate2d import *
#from fast_ftts import *
import warnings
from collections import Counter

def search_rms():
	search_string = 'Estimated noise='
	rms = []

	with open('difmap.log','r') as infile:
		for line in infile:
			if search_string in line:
				s=line
				s1 = s.split()
				s2 = s1[2].split('=')
				rms_value = float(s2[1])
				if s1[3][0]== 'm':
					rms.append(rms_value/1000.)
					
				elif s1[3][0] == 'J':
					rms.append(rms_value)

	return rms

def checking_shift(shifted_files,position,size,position_feature,size_feature,ifhdu):
	check2 = check_map_params(shifted_files[0],shifted_files[1],ifhdu)
	map_data1 = read_map(shifted_files[0],ifhdu)		
	realDATshift = map_data1[0]
	map_data2 = read_map(shifted_files[1],ifhdu)		
	realDAT2shift = map_data2[0]
	
	cutout_v1shift = Cutout2D(realDATshift,position,size)
	cutout_v2shift = Cutout2D(realDAT2shift,position,size)
	cutout_v1shift_feature = Cutout2D(cutout_v1shift.data,position_feature,size_feature)
	cutout_v2shift_feature = Cutout2D(cutout_v2shift.data,position_feature,size_feature)
	
	image1 = cutout_v1shift_feature.data
	image2 = cutout_v2shift_feature.data
		
	offset_new = cross_correlation_shifts_FITS(image1, image2, sigma_cut=0.004)	

	return offset_new
		

def find_same_beam(beam):
	#finding if some beams are repeated in the array
	beams = [item for item, count in Counter(beam).iteritems() if count > 1] 
	index_b = []	

	if len(beams) == 1:
		#if there is only one number repeated, take the indexes where that number is in the array
		index_b.append(np.where(beam == beams[0])) #output a list containing tuple
		#converting the list into array and the tuples in the list into an array
		index_beam = np.asarray([x for xs in index_b for x in xs]) 

	else:
		#if none or more than one is repeated, take the indexes where that happens in pairs of two
		for i in xrange(0,len(beams)):
			index_b.append(np.where(beam == beams[i]))
		index_beam = np.asarray([x for xs in index_b for x in xs]) 

	return index_beam

#create an array with beams for the two selected frequencies, to later find the same beam in both
def beam_array(checkBOXes,freq,freq_conv,beam_conv):
	freqs = []
	freq1_index_l = []
	beam1_l = []
	freq2_index_l = []
	beam2_l = []
	xshift = 0.
	yshift = 0.

	#obtaining the selected frequencies
	for i in xrange(len(checkBOXes)):
		if checkBOXes[i].isChecked():
			freqs.append(freq[i])

	#getting the index of the selected frequencies in the array with the convolved files 
	#(it probably has more than one frequency matching)	
	freq1_index_l.append(np.where(freq_conv == freqs[0]))
	freq2_index_l.append(np.where(freq_conv == freqs[1]))

	#converting tuples obtained in the previous line into arrays
	freq1_index = np.asarray([x for xs in freq1_index_l for x in xs])
	freq2_index = np.asarray([x for xs in freq2_index_l for x in xs])

	#obtaining the corresponding value of the beam for the indexes obtained for the two frequencies
	for i in xrange(0,len(freq1_index[0])):
		beam1_l.append(beam_conv[freq1_index[0][i]])
	for i in xrange(0,len(freq2_index[0])):
		beam2_l.append(beam_conv[freq2_index[0][i]])
	beam1 = np.asarray(beam1_l)
	beam2 = np.asarray(beam2_l)		

	#concatenate the beam of both frequencies in the same array
	beam12=np.concatenate((beam1,beam2))

	#see if there is a common number in the array, which will be the common beam
	index_beam12=find_same_beam(beam12)

	return freq1_index,freq2_index,index_beam12,beam1,beam2


def check_map_params(image1,image2,ifhdu):
	"""v1 < v2"""

	header1 = take_header(image1,ifhdu)
	header2 = take_header(image2,ifhdu)

	cell1 = header1[0]
	cell2 = header2[0]

	if cell1 == cell2:
		cells = header1[0]
		OK = True
		print('Same cellsize')
	else:
		cells = 0.
		print('The images do not have the same cellsize')
		print('Convolve them using the same cell')
		OK = False

	map_data1 = read_map(image1,ifhdu)		
	realDAT = map_data1[0]
	map_data2 = read_map(image2,ifhdu)		
	realDAT2 = map_data2[0]

	if realDAT.shape == realDAT2.shape:
		if OK == False:		
			OK = False
		else:
			OK = True
		print('Same mapsize')
	else:
		print('The images do not have the same mapsize')
		print('Convolve them using the same mapsize')
		OK = False

	beam1 = header1[7]
	beam2 = header2[7]

	if beam1 == beam2:
		beam = header1[7]
		if OK == False:		
			OK = False
		else:
			OK = True
		print('Same beam')
	else:
		beam = 0.
		print('The images do not have the same beam')
		print('Convolve them using the same beam')
		OK = False
		#raise IOError('Stopping!') compatible with python 3, same that raise Exception

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


	print('Proceeding with the next step')

	return OK, realDAT,realDAT2, freq1,freq2, cells, beam, cent_mapx, cent_mapy, x1,x2,y1,y2


def cuttingMAP(realDAT1,realDAT2,cent_mapx,cent_mapy,cells,freq1,freq2,freq1name,freq2name,freq1unit,freq2unit,iteration):

	#plotting the maps

	if iteration == 0:
		realDAT2 = realDAT2*(realDAT2 > realDAT2.std()*0.01) 
		realDAT2[realDAT2 == 0.0] = np.nan  
		realDAT1 = realDAT1*(realDAT1 > realDAT1.std()*0.01)  
		realDAT1[realDAT1 == 0.0] = np.nan  

	plt.figure(1)
	plt.subplot(121)
	plt.imshow(realDAT2, origin='bottom')
	#plt.axis('scaled')
	plt.xlabel('Right Ascension [pixels]')
	plt.ylabel('Relative Declination [pixels]')
	plt.title(freq2name+freq2unit)

	#plt.figure(2)
	plt.subplot(122)
	plt.imshow(realDAT1, origin='bottom')
	#plt.axis('scaled')
	plt.xlabel('Right Ascension [pixels]')
	plt.title(freq1name+freq1unit)

		
	#select with a box the part of the map to keep
	a = Annotate()
	plt.show()
		
	[limplot_x1,limplot_x2,limplot_y1,limplot_y2] = a()
	[limplot_x1mas,limplot_x2mas,limplot_y1mas,limplot_y2mas] = [(cent_mapx-a()[0])*cells,(cent_mapx-a()[1])*cells,(a()[2]-cent_mapy)*cells,(a()[3]-cent_mapy)*cells]
                                          
	ext_new = [limplot_x1mas,limplot_x2mas,limplot_y2mas,limplot_y1mas]

	#parameters to crop the map (center position of the rectangle and size)
	centx = ((limplot_x1+limplot_x2)/2.)
	centy = ((limplot_y1+limplot_y2)/2.) 
	position = (centx,centy)
	
	height = round(np.abs(limplot_y2-limplot_y1))+1
	width = round(np.abs(limplot_x2-limplot_x1))+1
	size = (height,width)

	#cropping
	cutout_v1 = Cutout2D(realDAT1,position,size)
	cutout_v2 = Cutout2D(realDAT2,position,size)

	"""
	plt.figure(3)
	plt.imshow(realDAT1, origin='bottom')#,extent=ext)#,extent=ext)
	cutout_v1.plot_on_original(color='white')

	plt.figure(4)
	plt.imshow(realDAT2, origin='bottom')#,extent=ext)#,extent=ext)
	cutout_v2.plot_on_original(color='white')

	"""

	plt.show()


	return cutout_v1, cutout_v2, position, size, ext_new



def cross_correlation_shifts(image1, image2, errim1=None, errim2=None,
        maxoff=None, verbose=False, gaussfit=False, return_error=False,
        zeromean=True, **kwargs):
    """ Use cross-correlation and a 2nd order taylor expansion to measure the
    offset between two images
    Given two images, calculate the amount image2 is offset from image1 to
    sub-pixel accuracy using 2nd order taylor expansion.
    Parameters
    ----------
    image1: np.ndarray
        The reference image
    image2: np.ndarray
        The offset image.  Must have the same shape as image1
    errim1: np.ndarray [optional]
        The pixel-by-pixel error on the reference image
    errim2: np.ndarray [optional]
        The pixel-by-pixel error on the offset image.  
    maxoff: int
        Maximum allowed offset (in pixels).  Useful for low s/n images that you
        know are reasonably well-aligned, but might find incorrect offsets due to 
        edge noise
    zeromean : bool
        Subtract the mean from each image before performing cross-correlation?
    verbose: bool
        Print out extra messages?
    gaussfit : bool
        Use a Gaussian fitter to fit the peak of the cross-correlation?
    return_error: bool
        Return an estimate of the error on the shifts.  WARNING: I still don't
        understand how to make these agree with simulations.
        The analytic estimate comes from
        http://adsabs.harvard.edu/abs/2003MNRAS.342.1291Z
        At high signal-to-noise, the analytic version overestimates the error
        by a factor of about 1.8, while the gaussian version overestimates
        error by about 1.15.  At low s/n, they both UNDERestimate the error.
        The transition zone occurs at a *total* S/N ~ 1000 (i.e., the total
        signal in the map divided by the standard deviation of the map - 
        it depends on how many pixels have signal)
    **kwargs are passed to correlate2d, which in turn passes them to convolve.
    The available options include image padding for speed and ignoring NaNs.
    References
    ----------
    From http://solarmuri.ssl.berkeley.edu/~welsch/public/software/cross_cor_taylor.pro
    Examples
    --------
    >>> import numpy as np
    >>> im1 = np.zeros([10,10])
    >>> im2 = np.zeros([10,10])
    >>> im1[4,3] = 1
    >>> im2[5,5] = 1
    >>> import image_registration
    >>> yoff,xoff = image_registration.cross_correlation_shifts(im1,im2)
    >>> im1_aligned_to_im2 = np.roll(np.roll(im1,int(yoff),1),int(xoff),0)
    >>> assert (im1_aligned_to_im2-im2).sum() == 0
    
    """

    if zeromean:
        image1 = image1 - (image1[image1==image1].mean())
        image2 = image2 - (image2[image2==image2].mean())

    image1 = np.nan_to_num(image1)
    image2 = np.nan_to_num(image2)

    quiet = kwargs.pop('quiet') if 'quiet' in kwargs else not verbose
    ccorr = (correlate2d(image1,image2,quiet=quiet,**kwargs) / image1.size)
    # allow for NaNs set by convolve (i.e., ignored pixels)
    ccorr[ccorr!=ccorr] = 0
    if ccorr.shape != image1.shape:
        raise ValueError("Cross-correlation image must have same shape as input images.  This can only be violated if you pass a strange kwarg to correlate2d.")

    ylen,xlen = image1.shape
    xcen = xlen/2-(1-xlen%2) 
    ycen = ylen/2-(1-ylen%2) 

    if ccorr.max() == 0:
        warnings.warn("WARNING: No signal found!  Offset is defaulting to 0,0")
        return 0,0

    if maxoff is not None:
        if verbose: print("Limiting maximum offset to %i" % maxoff)
        subccorr = ccorr[ycen-maxoff:ycen+maxoff+1,xcen-maxoff:xcen+maxoff+1]
        ymax,xmax = np.unravel_index(subccorr.argmax(), subccorr.shape)
        xmax = xmax+xcen-maxoff
        ymax = ymax+ycen-maxoff
    else:
        ymax,xmax = np.unravel_index(ccorr.argmax(), ccorr.shape)
        subccorr = ccorr

    if return_error:
        #if errim1 is None:
        #    errim1 = np.ones(ccorr.shape) * image1[image1==image1].std() 
        #if errim2 is None:
        #    errim2 = np.ones(ccorr.shape) * image2[image2==image2].std() 
        #eccorr =( (correlate2d(errim1**2, image2**2,quiet=quiet,**kwargs)+
                   #correlate2d(errim2**2, image1**2,quiet=quiet,**kwargs))**0.5 
                  # / image1.size)
        eccorr =( (correlate2d((image1*0.15)**2, image2**2,quiet=quiet,**kwargs)+
                   correlate2d((image2*0.15)**2, image1**2,quiet=quiet,**kwargs))**0.5 
                   / image1.size)
        if maxoff is not None:
            subeccorr = eccorr[ycen-maxoff:ycen+maxoff+1,xcen-maxoff:xcen+maxoff+1]
        else:
            subeccorr = eccorr

    if gaussfit:
        try: 
            from agpy import gaussfitter
        except ImportError:
            raise ImportError("Couldn't import agpy.gaussfitter; try using cross_correlation_shifts with gaussfit=False")
        if return_error:
            pars,epars = gaussfitter.gaussfit(subccorr,err=subeccorr,return_all=True)
            exshift = epars[2]
            eyshift = epars[3]
        else:
            pars,epars = gaussfitter.gaussfit(subccorr,return_all=True)
        xshift = maxoff - pars[2] if maxoff is not None else xcen - pars[2]
        yshift = maxoff - pars[3] if maxoff is not None else ycen - pars[3]
        if verbose: 
            print("Gaussian fit pars: ",xshift,yshift,epars[2],epars[3],pars[4],pars[5],epars[4],epars[5])

    else:

        xshift_int = xmax-xcen
        yshift_int = ymax-ycen

        local_values = ccorr[ymax-1:ymax+2,xmax-1:xmax+2]

        d1y,d1x = np.gradient(local_values)
        d2y,d2x,dxy = second_derivative(local_values)

        fx,fy,fxx,fyy,fxy = d1x[1,1],d1y[1,1],d2x[1,1],d2y[1,1],dxy[1,1]

        shiftsubx=(fyy*fx-fy*fxy)/(fxy**2-fxx*fyy)
        shiftsuby=(fxx*fy-fx*fxy)/(fxy**2-fxx*fyy)

        xshift = -(xshift_int+shiftsubx)
        yshift = -(yshift_int+shiftsuby)

        # http://adsabs.harvard.edu/abs/2003MNRAS.342.1291Z
        # Zucker error

        if return_error:
            #acorr1 = (correlate2d(image1,image1,quiet=quiet,**kwargs) / image1.size)
            #acorr2 = (correlate2d(image2,image2,quiet=quiet,**kwargs) / image2.size)
            #ccorrn = ccorr / eccorr**2 / ccorr.size #/ (errim1.mean()*errim2.mean()) #/ eccorr**2
            normalization = 1. / ((image1**2).sum()/image1.size) / ((image2**2).sum()/image2.size) 
            ccorrn = ccorr * normalization
            exshift = (np.abs(-1 * ccorrn.size * fxx*normalization/ccorrn[ymax,xmax] *
                    (ccorrn[ymax,xmax]**2/(1-ccorrn[ymax,xmax]**2)))**-0.5) 
            eyshift = (np.abs(-1 * ccorrn.size * fyy*normalization/ccorrn[ymax,xmax] *
                    (ccorrn[ymax,xmax]**2/(1-ccorrn[ymax,xmax]**2)))**-0.5) 
            if np.isnan(exshift):
                raise ValueError("Error: NAN error!")

    if return_error:
        return xshift,yshift,exshift,eyshift
    else:
        return xshift,yshift



def cross_correlation_shifts_FITS(fitsfile1, fitsfile2,
        return_cropped_images=False, quiet=True, sigma_cut=False,
        register_method=cross_correlation_shifts, **kwargs):
    """
    Determine the shift between two FITS images using the cross-correlation
    technique.  Requires montage or hcongrid.
    Parameters
    ----------
    fitsfile1: str
        Reference fits file name
    fitsfile2: str
        Offset fits file name
    return_cropped_images: bool
        Returns the images used for the analysis in addition to the measured
        offsets
    quiet: bool
        Silence messages?
    sigma_cut: bool or int
        Perform a sigma-cut before cross-correlating the images to minimize
        noise correlation?
    """
    #import montage
    try:
        import astropy.io.fits as pyfits
        #import astropy.wcs as pywcs
    except ImportError:
        import pyfits
        #import pywcs
    #import tempfile

    image2 = fitsfile2
    image1 = fitsfile1
    

    if sigma_cut:
        corr_image1 = image1*(image1 > image1.std()*sigma_cut)
        corr_image2 = image2*(image2 > image2.std()*sigma_cut)
        OK = (corr_image1==corr_image1)*(corr_image2==corr_image2) 
        if (corr_image1[OK]*corr_image2[OK]).sum() == 0:
            print("Could not use sigma_cut of %f because it excluded all valid data" % sigma_cut)
            corr_image1 = image1
            corr_image2 = image2
    else:
        corr_image1 = image1
        corr_image2 = image2

    verbose = kwargs.pop('verbose') if 'verbose' in kwargs else not quiet
    xoff,yoff = register_method(corr_image1, corr_image2, verbose=verbose,**kwargs)
    
    return xoff,yoff


def second_derivative(image):
    """
    Compute the second derivative of an image
    The derivatives are set to zero at the edges
    Parameters
    ----------
    image: np.ndarray
    Returns
    -------
    d/dx^2, d/dy^2, d/dxdy
    All three are np.ndarrays with the same shape as image.
    """
    shift_right = np.roll(image,1,1)
    shift_right[:,0] = 0
    shift_left = np.roll(image,-1,1)
    shift_left[:,-1] = 0
    shift_down = np.roll(image,1,0)
    shift_down[0,:] = 0
    shift_up = np.roll(image,-1,0)
    shift_up[-1,:] = 0

    shift_up_right = np.roll(shift_up,1,1)
    shift_up_right[:,0] = 0
    shift_down_left = np.roll(shift_down,-1,1)
    shift_down_left[:,-1] = 0
    shift_down_right = np.roll(shift_right,1,0)
    shift_down_right[0,:] = 0
    shift_up_left = np.roll(shift_left,-1,0)
    shift_up_left[-1,:] = 0

    dxx = shift_right+shift_left-2*image
    dyy = shift_up   +shift_down-2*image
    dxy=0.25*(shift_up_right+shift_down_left-shift_up_left-shift_down_right)

    return dxx,dyy,dxy




