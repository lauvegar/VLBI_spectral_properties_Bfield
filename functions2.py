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

#########################################################################################################
#########################################################################################################
############################################FUNCTIONS####################################################
#########################################################################################################
#########################################################################################################

#FUNCTION TO SAVE DATA TO A TEXT FILE
def saver(filename, header, data, FORMAT='%1.3e'):
    '''
    Given a header and data array, save to file 'filename'
    '''
    f = open(filename,'w')
    np.savetxt(f, header, delimiter='\t',fmt="%s")
    f.close()
    f_handle = open(filename, 'a')
    np.savetxt(f_handle, data, delimiter='\t',fmt=FORMAT)
    f_handle.close()  

#CREATE THE NECESSARY DIRECTORIES TO RUN THE CODE
def create_dir(path,num_gauss):
	if not os.path.exists(path+'Profiles_plot'):
		os.makedirs(path+'Profiles_plot')

	if not os.path.exists(path+'Pickle'):
		os.makedirs(path+'Pickle')
	for i in xrange(1,num_gauss+1):
		if not os.path.exists(path+'Pickle/'+str(num_gauss)+'gaussian'):
			os.makedirs(path+'Pickle/'+str(num_gauss)+'gaussian')

	if not os.path.exists(path+'Plots_Minuit'):
		os.makedirs(path+'Plots_Minuit')
	for i in xrange(1,num_gauss+1):
		if not os.path.exists(path+'Plots_Minuit/'+str(num_gauss)+'gaussian'):
			os.makedirs(path+'Plots_Minuit/'+str(num_gauss)+'gaussian')
	for i in xrange(1,num_gauss+1):
		if not os.path.exists(path+'Plots_Minuit/'+str(num_gauss)+'gaussian/PNG'):
			os.makedirs(path+'Plots_Minuit/'+str(num_gauss)+'gaussian/PNG')

	if not os.path.exists(path+'Residuals'):
		os.makedirs(path+'Residuals')
	for i in xrange(1,num_gauss+1):
		if not os.path.exists(path+'Residuals/'+str(num_gauss)+'gaussian'):
			os.makedirs(path+'Residuals/'+str(num_gauss)+'gaussian')

	if not os.path.exists(path+'Ridge_line_plots'):
		os.makedirs(path+'Ridge_line_plots')
	for i in xrange(1,num_gauss+1):
		if not os.path.exists(path+'Ridge_line_plots/'+str(num_gauss)+'gaussian'):
			os.makedirs(path+'Ridge_line_plots/'+str(num_gauss)+'gaussian')
	for i in xrange(1,num_gauss+1):
		if not os.path.exists(path+'Ridge_line_plots/'+str(num_gauss)+'gaussian/Amplitude'):
			os.makedirs(path+'Ridge_line_plots/'+str(num_gauss)+'gaussian/Amplitude')
	for i in xrange(1,num_gauss+1):
		if not os.path.exists(path+'Ridge_line_plots/'+str(num_gauss)+'gaussian/Sigma'):
			os.makedirs(path+'Ridge_line_plots/'+str(num_gauss)+'gaussian/Sigma')

#READ THE NECESSARY FILES, CONTAINING THE FLUX, RA AND DEC POSITION IN THE MAP (in profile files)
#AND THE MID POINT OF EACH SLICE (in other file)
#POSITIONS IN BOTH FILES ARE IN PIXELS
def read_profiles(path,directory,mid_point_file,cell,number_read):
	'''
	given path where the directory of the slices file are and where the mid_point file is
        (same path for both)
	directory = directory with the slice profiles
	mid_point_file = name of the file containing the position in pixels of the mid poing x (ra) and y  (dec) of the slice
	cell = cell size of the map
	number_read = total number of slice profiles
	'''

	#read the mid point of the slices

	#initializing lists
	split_f = []
	mid_point_x = []
	mid_point_y = []

	#read the mid_point file and creating two LISTS with x and y positions	
	temp2= path+mid_point_file
	temp_file2=open(temp2,mode='r')
	temp_file2.readline()
	for line in temp_file2:
		split_f= line.split()
		mid_point_x.append(float(split_f[1]))
		mid_point_y.append(float(split_f[2]))

	#read the x,z position and flux for each point in the slices

	#initializing lists
	split=[]
	x_data=[]
	x_all_list=[]
	x_all=[]
	x_data_cent=[]
	x_all_cent_list=[]
	x_all_cent=[]
	x_all_array=[]
	x_all_cent_array=[]

	z_data_cent=[]
	z_all_cent_list=[]

	y_data=[]
	y_all_list=[]
	y_all =[]

	#read the profiles files and creating LISTS with 
	#x position (x_all_list), y position (z_all_list)
        #and flux (y_all_list)
	for i in xrange(0,number_read+1):
		temp = path+directory+'/width_'+str(i)+'.txt'
		temp_file=open(temp,mode='r')
		temp_file.readline()
	
		for line in temp_file:
			split = line.split()
			y = float(split[2])*1000
#			if (y) > 0.5:
			if (y) > -100.0:
				x_data.append(float(split[0]))
				x_data_cent.append(-(float(split[0])-mid_point_x[0])*cell)
				z_data_cent.append((float(split[1])-mid_point_y[0])*cell)   
				y_data.append(float(split[2])*1000)

		x_all_list.append(x_data)
		x_all_cent_list.append(x_data_cent)
		z_all_cent_list.append(z_data_cent)
		y_all_list.append(y_data)

		x_data=[]
		x_data_cent=[]
		z_data_cent=[]	
		y_data=[]

	#creating an ARRAY for x_all_list, z_all_list and y_all_list lists
	x_all=np.array(x_all_list)
	x_all_cent=np.array(x_all_cent_list)
	y_all=np.array(y_all_list)
	z_all_cent=np.array(z_all_cent_list)

	for i in xrange(0,number_read+1):
		x_all_cent[i]=np.array(x_all_cent[i])
	
	y_all=np.array(y_all_list)
	for i in xrange(0,number_read+1):
		y_all[i]=np.array(y_all[i])

	#calculating the x and y position in the map in mas for each slice

	#initializing lists
	x_pos=[]
	y_pos=[]

	#creating two LIST for the estimation of x (x_pos_array) and y (y_pos_array) positions in the map
	for i in xrange(0,number_read+1):
		x_pos.append(-(mid_point_x[i]-mid_point_x[0])*cell)
		y_pos.append((mid_point_y[i]-mid_point_y[0])*cell)

	return x_all_cent,z_all_cent,y_all,x_pos,y_pos,mid_point_x,mid_point_y


#LOAD GAUSSIAN FIT PARAMETERS FROM FILES
def load_gauss_file(gauss_file,gauss_file_rot, number_read, num_gauss):
	'''
	given gauss_file = file with the arrays of gaussian parameters result of the fitting
	gauss_file_rot =  file with the arrays of gaussian parameters result of the fitting, with x and y rotated with respect of the jet P.A.
	number_read = total number of slice profiles
	num_gauss = number of gaussians used for the fitting
	'''

	#initializing the arrays
	slices=np.array([[0.]*(number_read+1)])

	A_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	x0_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	y0_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	sigma_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]

	Aerr_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	x0err_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	sigmaerr_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]

	red_chi_2_arr = np.array([[0.]*(number_read+1)])

	x_data_rot=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	y_data_rot=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	y_data_rot_err=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]

	#loading non rotated values
	print gauss_file
	temp3 = gauss_file
	data1 = np.loadtxt(temp3)

	slices[0] = data1[:,0]
	red_chi_2_arr[0] = data1[:,1]

	counter = 0.
	for i in xrange(0,num_gauss):
		y0_arr[i] = data1[:,(2 + counter)]
		x0err_arr[i] = data1[:,(3 + counter)]		
		x0_arr[i] = data1[:,(4 + counter)]
		A_arr[i] = data1[:,(5 + counter)]
		Aerr_arr[i] = data1[:,(6 + counter)]
		sigma_arr[i] = data1[:,(7 + counter)]
		sigmaerr_arr[i] = data1[:,(8 + counter)]
		counter = counter + 7
	
	#loading rotated values
	print gauss_file_rot
	temp3 = gauss_file_rot 
	data1 = np.loadtxt(temp3)

	counter = 0.
	for i in xrange(0,num_gauss):	
		x_data_rot[i] = data1[:,(2 + counter)]
		y_data_rot[i] = data1[:,(3 + counter)]
		y_data_rot_err[i] = data1[:,(4 + counter)]
		counter = counter + 7

	#calculating the min and max profile used (min profile set to 0)
	min_profile = int(slices[0][0])
	counter = False
	for i in xrange(1,len(slices[0])):
		if counter == False:			
			if slices[0][i] == 0. :
				max_profile_index = i-1
				counter = True
				max_profile = int(slices[0][max_profile_index])
			else:
				max_profile = int(slices[0][len(slices[0])-1])

	return slices,red_chi_2_arr,y0_arr,x0err_arr,x0_arr, A_arr, Aerr_arr, sigma_arr, sigmaerr_arr, x_data_rot, y_data_rot, y_data_rot_err, min_profile, max_profile

#FUNCTION TO READ THE HEADER AND TAKE IMPORTANT PARAMETERS AS
#cell
#BMAJ, BMIN, BPA
#date, freq and epoch
def take_header(x,ifhdu):
	'''
	given x MAP file to be read
	'''

	#open the fits file and read the header
	#FITS = pf.open(x)
	if ifhdu == False:
		HEADER = pf.getheader(x)
	if ifhdu == True:
		hdul = pf.open(x)
		HEADER = hdul[0].header

	#read the parameters that contain the cell, bmaj, bmin and bpa 
	cell = np.abs(HEADER['CDELT1']*3600*1000)  #= -9.7222222253244E-08 /Pixel increment        
	BMAJ =HEADER['BMAJ']*3600*1000
	BMIN = HEADER['BMIN']*3600*1000  # /Reference pixel    
	BPA = HEADER['BPA']

	#read the date as it is in the header
	date = HEADER['DATE-OBS']

	#for dates in the format year-month-day
	if date[2].isdigit():
		d=date.split('-')
		date=str(d[0])+'_'+str(d[1])+'_'+str(d[2])
		epoch=(float(d[0])+((float(d[1])-1)*30.4375+float(d[2]))/365.25)
	#for dates in the format day/month/year (with only the two last digits of the year)
	else:
		d=date.split('/')
		print d[2][0]
		if d[2][0]=='9':
                	d[2]=1900+float(d[2])
                	date=str("%0.0f" % (d[2],))+'_'+str(d[1])+'_'+str(d[0])
                	epoch=(float(d[2])+((float(d[1])-1)*30.4375+float(d[0]))/365.25)        
		else:
                	d[2]=2000+float(d[2])        
                	date=str(d[2])+'_'+str(d[1])+'_'+str(d[1])
                	epoch=(float(d[2])+((float(d[1])-1)*30.4375+float(d[0]))/365.25)
 	freq = float(HEADER['CRVAL3'])/1e9

	#equivalent cicular beam
	circ_beam = np.sqrt(BMAJ*BMIN)

	source_name = HEADER['OBJECT']

	rms = HEADER['NOISE']
	datamax = HEADER['DATAMAX']

	return cell, BMAJ,BMIN,BPA,date,freq,epoch,circ_beam, source_name,rms,datamax

#READ MAP DATA AND DEFINE THE BOX WHERE THE DATA ARE (difmap fits files)
def read_map(x,ifhdu):
	'''
	given x MAP file to be read
	'''

	#open the fits file and read the header. Get the data where the map value are
	#FITS = pf.open(x)

	if ifhdu == False:
		HEADER = pf.getheader(x)
		DATA=pf.getdata(x)
		realDAT=DATA[0][0]
	if ifhdu == True:
		hdul = pf.open(x)
		HEADER = hdul[0].header
		realDAT = hdul[0].data
	
	

	#read from the header for RA and DEC axis: 
	#the lenght of the axis    
	#the reference pixel
	#the reference value
	#the pixel increment (cell)
	raAXIS =HEADER['NAXIS1'] #mapsize
	raPIX = HEADER['CRPIX1']  #=  2048.00000000000  /Reference pixel                       
	raVAL = HEADER['CRVAL1']  # /Reference value                       
	raDEL = HEADER['CDELT1']  # /Pixel increment    (cellsize)       
    
	dcAXIS =HEADER['NAXIS2']  #mapsize
	dcPIX = HEADER['CRPIX2']  #=  2049.00000000000  /Reference pixel                       
	dcVAL = HEADER['CRVAL2']  # /Reference value                       
	dcDEL = HEADER['CDELT2']  # /Pixel increment    (cellsize)       
        
	#calculate the box where the map is
	x1= (raVAL - raPIX*raDEL)*60*60*1000-raVAL*60*60*1000
	x2= (raVAL + raPIX*raDEL)*60*60*1000-raVAL*60*60*1000
	y1= (dcVAL - dcPIX*dcDEL)*60*60*1000-dcVAL*60*60*1000
	y2= (dcVAL + dcPIX*dcDEL)*60*60*1000-dcVAL*60*60*1000

	return realDAT, x1, x2, y1, y2, raPIX, dcPIX, raAXIS, dcAXIS, (raVAL,dcVAL)

#1,2,3 GAUSSIANS
def gauss(x,A, x0, sigma ):
	'''
	given given x = [[set of x values],...[]]
	a = peak flux value
	x0 = peak position value
	sigma = fwhm value
	'''

	func = A*np.exp(-8*np.log(2)*(x-x0)**2./(2.*sigma**2))

	return func

def gauss2(x,A, x0, sigma, A2, x02, sigma2):

	return A*np.exp(8*np.log(2)*(x-x0)**2./(-2.*sigma**2))+A2*np.exp(8*np.log(2)*(x-x02)**2./(-2.*sigma2**2))

def gauss3(x,A, x0, sigma, A2, x02, sigma2,A3,x03,sigma3):
    
	return A*np.exp(8*np.log(2)*(x-x0)**2./(-2.*sigma**2))+A2*np.exp(8*np.log(2)*(x-x02)**2./(-2.*sigma2**2))+A3*np.exp(8*np.log(2)*(x-x03)**2./(-2.*sigma3**2))

#OBTAINING THE GUESS PARAMETERS INTERACTIVELY
def guesses(i,num_gauss,xdata,ydata):
	'''
	given i = slice number
	num_gauss = number of gaussians used for the fitting
	xdata = array containing the lists for each slice where the X POSITION in the map is stored
	ydata = array containing the lists for each slice where the FLUX in the map is stored
	'''

	plt.ion()

	#initialize the arrays
	guess_A = np.array([0.]*num_gauss)
	guess_x0 = np.array([0.]*num_gauss)
	guess_sigma = np.array([0.]*num_gauss)

	#plot the profile for the slice	
	plt.figure(1)
	plt.plot(xdata[i],ydata[i], '.')
	plt.xlabel('d [mas]')
	plt.ylabel('f [mJy]')
	print 'Please, click on the profile.'
	print 'The first click will take the guess value of the amplitude.'
	print 'The second click will take the guess value of the center position.'
	print 'The difference of the third and fourth click will take the sigma of the gaussian.'
	print 'The next four clicks will do the same for the second gaussian, if there is one'
	print 'Same procedure for any other gaussian'

	#set the number of clics needed in total and store them in an array called param
	param = ginput(4*num_gauss,0) 

	#store the input parameters given in the clics in the correspondent guess parameter type 
	counter = 0
	for j in xrange(0,num_gauss*4,4):		
		guess_A[counter] = float(param[j][1])	
		guess_x0[counter] = float(param[j+1][0])
		guess_sigma[counter] = np.abs(float(param[j+2][0])-float(param[j+3][0])) 

		print 'These are your guessing parameters:'
		print 'A'+str(counter+1)+': ', guess_A[counter]
		print 'x0'+str(counter+1)+': ', guess_x0[counter]
		print 'sigma'+str(counter+1)+': ', guess_sigma[counter]

		counter = counter + 1

	return guess_A,guess_x0,guess_sigma

#FUNCTION TO PLOT THE DIFFERENT PROFILES WITH ONE GAUSSIAN FIT AND ITS RESIDUALS
def plot(i,xdata,ydata,args,xpos,ypos,residuals,num_gauss):
	'''
	given i = slice number
	xdata = array containing the lists for each slice where the X POSITION in the map is stored
	ydata = array containing the lists for each slice where the FLUX in the map is stored
	args = parameters obtained by fitting the gaussians: A, ... An, x0 ... x0n, sigma... sigman
	xpos = list with the estimated x(ra) position of the slice in the map
	ypos = list with the estimated y(dec) position of the slice in the map
	residuals = residuals obtained for the fit
	num_gauss = number of gaussians used for the fitting
	'''

	plt.figure(i)
	plt.plot(xdata,ydata, '.')
	if num_gauss == 1:	    
		yfit=gauss(xdata,args[0],args[1],args[2])
		plt.plot(xdata, yfit,'r-')
	if num_gauss == 2:	    
		yfit=gauss2(xdata,args[0],args[1],args[2],args[3],args[4],args[5])
		peak1=gauss(xdata,args[0],args[1],args[2])
		peak2=gauss(xdata,args[3],args[4],args[5])
		plt.plot(xdata, yfit,'r-')
		plt.plot(xdata, peak1,'c--')
		plt.plot(xdata, peak2,'c-.')
	if num_gauss == 3:	    
		yfit=gauss2(xdata,args[0],args[1],args[2],args[3],args[4],
                            args[5],args[6],args[7],args[8])
		peak1=gauss(xdata,args[0],args[1],args[2])
		peak2=gauss(xdata,args[3],args[4],args[5])
		peak3=gauss(xdata,args[6],args[7],args[8])
		plt.plot(xdata, yfit,'r-')
		plt.plot(xdata, peak1,'c--')
		plt.plot(xdata, peak2,'c-.')
		plt.plot(xdata, peak3,'c-')
	plt.plot(xdata, residuals,'g-') 
	plt.xlabel('d [mas]')
	plt.ylabel('f [mJy]')
	plt.title('x='+str("%0.2f" % (xpos,))+' y='+str("%0.2f" % (ypos,)))
	plt.savefig('Plots_Minuit/'+str(num_gauss)+'gaussian/'+str(num_gauss)+'gaussian'+str(i)+'.ps')
	plt.savefig('Plots_Minuit/'+str(num_gauss)+'gaussian/PNG/'+str(num_gauss)+'gaussian'+str(i)+'.png')


#RECOVER THE INITIAL PARAMETERS OF A FIT RESTORED IN A PICKLE (one gaussian)
def pickle_guess(slice_num,num_gauss):
	'''
	given slice_num = slice number
	num_gauss = number of gaussians used for the fitting
	'''

	#read the pickle file. Store the parameters in an array called pick
	res = open('Pickle/'+str(num_gauss)+'gaussian/pickle_minuit'+str(slice_num)+'.p','rb')
	pick = pickle.load(res)
	res.close()

	#initialize the arrays for each parameter
	guess_A = np.array([0.]*num_gauss)
	guess_x0 = np.array([0.]*num_gauss)
	guess_sigma = np.array([0.]*num_gauss)

	#store the pickle parameters in the correspondent guess parameter type 
	for j in xrange(0,num_gauss):
		guess_A[j] = float(pick[0][j])	
		guess_x0[j] = float(pick[1][j])
		guess_sigma[j] = float(pick[2][j]) 

		print 'These are your guessing parameters:'
		print 'A'+str(j+1)+': ', guess_A[j]
		print 'x0'+str(j+1)+': ', guess_x0[j]
		print 'sigma'+str(j+1)+': ', guess_sigma[j]

	return guess_A,guess_x0,guess_sigma

#RECOVER THE LIMITS OF FITS, RESTORED IN A PICKLE (this may the default)

def pickle_lim(slice_num,num_gauss):
	'''
	given slice_num = slice number
	num_gauss = number of gaussians used for the fitting
	'''

	#read the pickle file. Store the parameters in an array called pick
	res=open('Pickle/'+str(num_gauss)+'gaussian/pickle_minuit_limits'+str(slice_num)+'.p','rb')
	pick = pickle.load(res)
	res.close()  

	#initialize the arrays for each limit needed
	min_A = np.array([0.]*num_gauss)
	min_x0 = np.array([0.]*num_gauss)
	max_x0 = np.array([0.]*num_gauss)
	min_sigma = np.array([0.]*num_gauss)
	max_sigma = np.array([0.]*num_gauss)

	#store the pickle parameters in the correspondent limit array 
	for j in xrange(0,num_gauss):
		min_A[j] = float(pick[0][j])	
		min_x0[j] = float(pick[1][j])
		max_x0[j] = float(pick[2][j]) 
		min_sigma[j] = float(pick[3][j])
		max_sigma[j] = float(pick[4][j]) 

		print 'These are your limits:'
		print 'min_A'+str(j+1)+': ', min_A[j]
		print 'min_x0'+str(j+1)+': ', min_x0[j]
		print 'max_x0'+str(j+1)+': ', max_x0[j]
		print 'min_sigma'+str(j+1)+': ', min_sigma[j]
		print 'max_sigma'+str(j+1)+': ', max_sigma[j]
	
	return min_A,min_x0,max_x0,min_sigma,max_sigma

#FUNCTION TO CALCULATE THE LIMITS OF THE VARIABLES WHEN FITTING (CONTINUITY CONDITIONS)
def limits_1gauss(A,x0,sigma,beam,num_gauss):
	'''
	given A = [A_1gauss,...,A_ngauss] array
	x0 = [x0_1gauss,...,x0_ngauss] array
	sigma = [sigma_1gauss,...,sigma_ngauss] array
	beam = beam of the image
	num_gauss = number of gaussians used for the fitting
	'''

	#initialize the arrays for each limit needed
	min_A = np.array([0.]*num_gauss)
	max_A = np.array([0.]*num_gauss)
	min_x0 = np.array([0.]*num_gauss)
	max_x0 = np.array([0.]*num_gauss)
	min_sigma = np.array([0.]*num_gauss)
	max_sigma = np.array([0.]*num_gauss)

	#calculate the limits for each parameter
	#A is within the 30% of the previous value
	#x0 is within 1/4 of the beam of the previous value
	#sigma is within the 20% of the previous value
	for j in xrange(0,num_gauss):
		if A[j]-0.3*A[j] > 0.:
			min_A[j] = A[j]-0.3*A[j]
		else:
			min_A[j] = 0.
  
		max_A[j] = A[j]+0.3*A[j]
		print A[j], max_A[j]

#		if m.values['x0'] > 0:								
		min_x0[j] = x0[j]-1/4.*beam
#		else:								
#			min_x0 = m.values['x0']+1/4.*beam	

#		if m.values['x0'] > 0:								
		max_x0[j] = x0[j]+1/4.*beam
#		else:								
#			max_x0 = m.values['x0']-1/4.*beam

		min_sigma[j] = np.abs(sigma[j])-0.2*np.abs(sigma[j])
		max_sigma[j] = np.abs(sigma[j])+0.2*np.abs(sigma[j])

	return min_A, max_A, min_x0, max_x0, min_sigma, max_sigma

#WHEN A FIT FAILS AND YOU START IT OVER WITHOUT CONTINUITY CONDITIONS, SOMETIMES
#THE RESULTS OF THE VARIABLES ARE NON LOGICAL. THIS ALLOWS TO RESTRICT THOSE VALUES
#ONLY USE IT WHEN IT IS STRICTLY NECESSARY, IF NOT, AVOID!!!!     
def questions_1gauss(xdata,i,num_gauss,sigma1g):				
	'''
	given xdata = array containing the lists for each slice where the X POSITION in the map is stored
	i = slice number
	num_gauss = number of gaussians used for the fitting
	sigma1g = fwhm for the 1gaussian fit, if it exist, if not an array containing a zero
	'''

	#initialize the arrays for each limit needed
	min_A = np.array([0.]*num_gauss)
	max_A = np.array([0.]*num_gauss)
	min_x0 = np.array([0.]*num_gauss)
	max_x0 = np.array([0.]*num_gauss)
	min_sigma = np.array([0.]*num_gauss)
	max_sigma = np.array([0.]*num_gauss)

	#ask if the limit need to be changed. 
	#If yes, it stores the new value in the correspondent limit array
	#If not, it takes the default values
	for j in xrange(0,num_gauss):

		track_change = False
		while track_change == False:
			prompt = 'The lower limit for A'+str(j+1)+' is 0. Do you want to change it? (y/N) \n'
			change = str(raw_input(prompt))
			if 'N' in change or 'y' in change or 'n' in change or len(change)==0:
				track_change = True
		
		if change == 'y':					
			prompt = 'Which lower limit do you want A'+str(j+1)+' to take? \n'
			min_A[j] = float(raw_input(prompt))
		else:
			min_A[j] = float(0)

		track_change = False
		while track_change == False:
			prompt = 'The lower limit for x0'+str(j+1)+' is the minimum value the profile takes in the x-axis. Do you want to change it? (y/N) \n'
			change = str(raw_input(prompt))
			if 'N' in change or 'y' in change or 'n' in change or len(change)==0:
				track_change = True
		
		if change == 'y':		
			prompt = 'Which lower limit do you want x0'+str(j+1)+' to take? \n'
			min_x0[j] = float(raw_input(prompt))
		else:
			min_x0[j] = np.min(xdata)

		track_change = False
		while track_change == False:
			prompt = 'The upper limit for x0'+str(j+1)+' is the maximum value the profile takes in the x-axis. Do you want to change it? (y/N) \n'
			change = str(raw_input(prompt))
			if 'N' in change or 'y' in change or 'n' in change or len(change)==0:
				track_change = True
		
		if change == 'y':									
			prompt = 'Which upper limit do you want x0'+str(j+1)+' to take? \n'
			max_x0[j] = float(raw_input(prompt))
		else:
		    	max_x0[j] = np.max(xdata)

		track_change = False
		while track_change == False:
			prompt = 'The lower limit for sigma'+str(j+1)+' is zero. Do you want to change it? (y/N) \n'
			change = str(raw_input(prompt))
			if 'N' in change or 'y' in change or 'n' in change or len(change)==0:
				track_change = True
		
		if change == 'y':					
			prompt = 'Which lower limit do you want sigma'+str(j+1)+' to take? \n'
			min_sigma[j] = float(raw_input(prompt))
		else:
		   	min_sigma[j] = float(0)

		track_change = False
		while track_change == False:
			prompt = 'The upper limit for sigma'+str(j+1)+' is the maximum value the profile takes in the x-axis. Do you want to change it? (y/N) \n'
			change = str(raw_input(prompt))
			if 'N' in change or 'y' in change or 'n' in change or len(change)==0:
				track_change = True
		
		if change == 'y':									
			prompt = 'Which upper limit do you want sigma'+str(j+1)+' to take? \n'
			max_sigma[j] = float(raw_input(prompt))
		else:
			if num_gauss > 1:
				max_sigma[j] = sigma1g[i]
			else:
				max_sigma[j] = np.abs(2*np.max(xdata[i]))

	return min_A,min_x0,max_x0,min_sigma,max_sigma
  

#FUNCTION TO ROTATE THE RIDGE LINE WITH RESPECT TO A P.A.    
def rotate_ridge_line(x0_arr,y0_arr,x0err_arr,beta,number_read,num_gauss):
	'''
	given x0_arr = position of the central position of the gaussian for the RA axis
	y0_arr = position of the central position of the gaussian for the DEC axis
	x0err_arr = error of the positions
	beta = position angle of the jet
	number_read = total number of slice profiles
	num_gauss = number of gaussians used for the fitting
	'''

	#initialize the arrays 
	x_data_rot=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	y_data_rot=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	y_err=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]

	#calculate the rotation
	for i in xrange(0,num_gauss):
		for j in xrange(0,len(x0_arr[0])):
			y_data_rot[i][j]=y0_arr[i][j]*np.sin(2*math.pi-beta)+x0_arr[i][j]*np.cos(2*math.pi-beta)
			x_data_rot[i][j]=y0_arr[i][j]*np.cos(2*math.pi-beta)-x0_arr[i][j]*np.sin(2*math.pi-beta)
   			y_err[i][j]=np.sqrt((np.sin(2*math.pi-beta)*x0err_arr[i][j])**2+(np.cos(2*math.pi-beta)*x0err_arr[i][j])**2)
        
	return x_data_rot,y_data_rot,y_err

#FUNCTION TO CALCULATE THE ERRORS OF THE RIDGE LINE THROUGHT SNR
def errors(yi,ye,beam,A,fwhm,number_read,num_gauss):
	'''
	given yi = rms close to the first slice
	ye = rms close to the last slice
	beam = beam of the image
	A = [[A_1gauss_1, .., A_1gauss_n],...,[A_ngauss, .., A_ngauss_n ] array and subarrays
	fwhm = [[sigma_1gauss_1, .., sigma_1gauss_n],...,[sigma_ngauss, .., sigma_ngauss_n ] 
               array and subarrays
	number_read = total number of slice profiles
	num_gauss = number of gaussians used for the fitting
	'''

	x1=1.
	x2=number_read+1.
	y1=yi
	y2=ye

	step=(y1-y2)/(x1-x2)
	b=y1-step*x1

	#initialize the arrays needed
	rms=np.array([[0.]*(number_read+1)])
	Aerr_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	x0err_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	sigmaerr_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]
	dlim_arr=[[0.]*(number_read+1),[0.]*(number_read+1),[0.]*(number_read+1)]

	#calculates the rms for each slice
	for i in xrange(0,number_read+1):
		rms[0][i] = step*(i+1) + b

	#calculates the error for each parameter
	for i in xrange(0,num_gauss):
		for j in xrange(0,number_read+1):  
			if A[i][j] != 0:
				snr = A[i][j]/rms[0][j]
    
				errA1 = rms[0][j]    

				dlim = 4/np.pi*beam*np.sqrt(np.pi*np.log(2)*np.log((snr+1.)/(snr)))

				if fwhm[i][j] > beam:
					ddec=np.sqrt(fwhm[i][j]**2-beam**2)
				else:
					ddec=0.

				y=[dlim,ddec]
				dg=np.max(y)        
    
				errx1 = rms[0][j]/A[i][j]*dg/(2.*np.sqrt(2))*np.sqrt(1+dg/dlim) 
				errFwhm1 = rms[0][j]/A[i][j]*dg/2.*(1+dg/dlim)

				Aerr_arr[i][j] = errA1
				sigmaerr_arr[i][j] = errFwhm1      
				x0err_arr[i][j] = errx1
				dlim_arr[i][j] = dlim
    
	return Aerr_arr,x0err_arr,sigmaerr_arr,dlim_arr

#FUNCTION TO OBTAIN THE COORDINATE OF A ELLIPSE, USED TO REPRESENT BEAMS IN THE IMAGES
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


def convolve_difmap(files,models,bmaj,bmin,bpa,xshift,yshift,mapsize,pixelsize,uvw1,uvw2,files2):
	#open difmap and convolve maps with new beam
	from distutils.spawn import find_executable
	difmap_path = find_executable('difmap')
	for i in xrange(0,len(files)): #/usr/local/pgplot.old/uvf_difmap/difmap 
		path = difmap_path #'/usr/local/bin/difmap'
		p=sub.Popen([path],stdin=sub.PIPE,stdout=sub.PIPE,close_fds=True,universal_newlines=True)
		(child_stdout,child_stdin)=(p.stdout,p.stdin)
		child_stdin.write('observe %s\n' %files[i])
		child_stdin.write('select PI\n')
		child_stdin.write('uvweight %d,%d\n'%(uvw1,uvw2)) 
		child_stdin.write('mapsize %s,%s\n' %(mapsize,pixelsize))
		child_stdin.write('rmod %s\n' %models[i])
		child_stdin.write('restore %s,%s,%s\n' %(bmaj,bmin,bpa))
		if float(xshift)!=0 and float(yshift)!=0:
			child_stdin.write('shift %3.3f, %3.3f\n' %(xshift,yshift))
		child_stdin.write('wmap %s\n' %files2[i])
		child_stdin.write('quit\n')
		p.wait()
	
	#remove log file
	#os.system('rm difmap.log*\n')


class Annotate(object):
    def __init__(self):
        self.ax = plt.gca()
        self.rect = Rectangle((0,0), 1, 1, facecolor='None', edgecolor='red')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
	self.x2 = None
	self.y2 = None 
	self.is_pressed = False
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        self.is_pressed = True
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_motion(self, event):
        self.x1, self.y1 = event.xdata, event.ydata
        if (self.is_pressed is True and
                self.x1 is not None and
                self.y1 is not None):
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.rect.set_linestyle('dashed')
            self.ax.figure.canvas.draw()

    def on_release(self, event):
        self.is_pressed = False
        self.x1, self.y1 = event.xdata, event.ydata
        try:
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
        except TypeError:
            if (self.x1 is None or self.y1 is None):
                return
            else:
                raise
        self.rect.set_linestyle('solid')
        self.ax.figure.canvas.draw()
        print self.x0,self.x1,self.y0,self.y1
	self.x2 = self.x1
	self.y2 = self.y1

    def __call__(self):
	return [self.x0,self.x2,self.y0,self.y2]
  

#########################################################################################################
#########################################################################################################
#######################################END_OF_THE_FUNCTIONS##############################################
#########################################################################################################
#########################################################################################################

def swap(x,y):

	temp = x
	x = y
	y = temp
	
	return x,y







### Crops image ###
#def crop(map_data_v1,map_data_v2):
	
    #my $im = $_[0];
    #my $x1 = $_[1];
    #my $y1 = $_[2];
    #my $x2 = $_[3];
    #my $y2 = $_[4];
    #my $xmin = 0;
    #my $ymin = 0;

    # Crop the image
    #if ($x1 > $x2) {
	#$xmin = $x2;
	#if ($y1 > $y2) {
	#    $im = $im($x2:$x1,$y2:$y1);
        #    $ymin = $y2;
	#} else {
	#    $ymin = $y1;
	#    $im = $im($x2:$x1,$y1:$y2);
	#}
    #} else {
	#$xmin = $x1;
	#if ($y1 > $y2) {
	 #   $im = $im($x1:$x2,$y2:$y1);
	 #   $ymin = $y2;
	#} else {
	#    $im = $im($x1:$x2,$y1:$y2);
	#    $ymin = $y1;
	#}
    #}

    # Fix the header
    #my @dim = $im->dims;
    #my $head = $im->gethdr;
   # my $centx = $$head{CRPIX1}; 
    #my $centy = $$head{CRPIX2};       
    #$$head{NAXIS1} = $dim[0]; 
    #$$head{NAXIS2} = $dim[1];
    #$$head{CRPIX1} = $centx-$xmin;
    #$$head{CRPIX2} = $centy-$ymin;
    #$im->sethdr($head);
	
    # Return cropped image
    #return $im;"""







