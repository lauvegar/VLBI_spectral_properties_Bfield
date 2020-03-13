import numpy as np
import scipy.special as scp
import cmath as math
import decimal
import math
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
from random import randint
from pylatex import Document, LongTabu, Tabu, Center, Tabular, Math, NewPage
from pylatex.utils import bold

def round_sigfigs(num, sig_figs=2):
	"""Round to specified number of sigfigs.
	"""
	if num != 0:
		return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
	else:
		return 0  # Can't take the log of 0


#def roundErrorVariable(variable,error):
#	newerror = round_sigfigs(error)
#	d = decimal.Decimal(str(newerror))
#	numberdecimals = d.as_tuple().exponent   
#	newvariable = round(variable,-numberdecimals)

#	return newvariable,newerror

def roundErrorVariable(x,y):
    ''' Function to round up a variable (x) to one significant digit of its error (y). 
        If the first digit of the error is 1, two digits are used. If the last digits 
        of the error are 10 (e.g., 0.10, 0.0010,  0.0000010), the last zero will not be printed and you will have to fix that yourself.'''
    if np.isnan(x):
        print 'Oops! The value you gave as input appears to be NaN! A zero will be returned'
	newX = 0
	newY = 0
    elif np.isnan(y):
        print 'Oops! The error you gave as input appears to be NaN! A zero will be returned'
	newX = 0
	newY = 0
    else:
	    a = int(y)
	    k = 0
	    if a > 20:
		c = int(np.log10(y))
		newY = int(np.round(y, decimals = -c))
		newX = int(np.round(x, decimals = -c))
	    elif a > 1 and a < 20:
		newY = int(np.round(y))
		newX = int(np.round(x))
	    elif a == 1:
		newY = float('%.1f'%(np.round(y, decimals = 1)))
		newX = float('%.1f'%(np.round(x, decimals = 1)))
	    elif a == 0 and y != 0.0:
		j = y
		while j < 1.:
		    j = j*10
		    k += 1
		if int(j) > 1:
		    newY = float('%.*f'%(k,np.round(y, decimals = k)))
		    newX = float('%.*f'%(k,np.round(x, decimals = k)))
		elif int(j) == 1:
		    newY = float('%.*f'%(k+1,np.round(y, decimals = k+1)))
		    newX = float('%.*f'%(k+1,np.round(x, decimals = k+1)))
                    
	    return newX, newY


def Tbformula(freq,size,flux,z):
	h = 6.6260755*10**(-27) #erg/s = g*cm2*s-2/s
	kb = 1.380658*10**(-16) #erg/k = g*cm2*s-2/k
	me = 9.11*10**(-28) #g
	e = 4.8*10**(-10) # g**(1/2)cm**(3/2)/s
	c = 3*10**10 #cm/s

	#pctocm = 3.08*10**(18) #cm 1 pc

	flux =flux*10**(-23)

	freq = freq*10**9
	size = size*np.pi/(6.48*10**8)

	logarithm = np.log(2*h*freq**3/c**2*size**2/flux+1)
	Tb = h*freq/kb*1./logarithm*(1.+z)

	#cte = 1.22*10**(12)

	#Tbapprox = cte*(1.+z)/(freq**2)*flux[i]/(dlim1[i]**2)
	Tbapprox = c**2/(2*kb)*flux/(freq**2*size**2)*(1.+z)


	return Tb, Tbapprox

def Bfield_Tb(Tb,z,gamma,delta,freq):
	#freq in GHz
	#Tb in 10**(12) K
	#only in the case p=2, if not it is kb**(-2)*m**3*c**5/e*k0(p)
	#whre k0(p) is a function depending of the p (p = -2*alpha0 + 1), how i dont know, i did not find it

	B = 7.4*10**(-4)*gamma*delta/(1.+z)*freq*Tb**(-2)

	return B

def constants(alpha0):

	p = -2*alpha0 + 1

	a_alpha =  2**(-alpha0)*np.sqrt(3)/(8.*np.sqrt(np.pi)*(2-2*alpha0))*scp.gamma((6.-2*alpha0)/4.)*scp.gamma((2.-6.*alpha0)/12.)*scp.gamma((22.-6.*alpha0)/12.)*(scp.gamma((8.-2.*alpha0)/4.))**(-1) 

	b_alpha = np.pi**2*a_alpha*3**((p-1.)/2.)*2**((p+9.)/2.)*(p**2+4*p+11.)*((p+3.)**2*(p+1.)*(p+5.))**(-1)	

	c_alpha = 3**(1.-alpha0)/8.*np.sqrt(np.pi)*scp.gamma((7.-2*alpha0)/4.)*scp.gamma((5.-6.*alpha0)/12.)*scp.gamma((25.-6.*alpha0)/12.)*(scp.gamma((9.-2.*alpha0)/4.))**(-1)

	return a_alpha,b_alpha,c_alpha


def Magnetization(alpha0,gammamin,gammamax,B,scale,viewing,opening,z,gamma,delta,kr,freq,Tb,coreshiftmeas):


	me = 9.11*10**(-28) #g
 	e = 4.8*10**(-10) #
	c = 3*10**10 #cm/s

	a_alpha,b_alpha,c_alpha = constants(alpha0)

	#p = -2*alpha0 + 1
	#freq in GHz
	#Tb in 10**(12) K

	p = -2*alpha0 + 1

	if p == 2:

		f_p = np.log(gammamax/gammamin)
	else:
		print 'p>2'
		f_p = 1./(2.-p)*(gammamax**(2.-p)-gammamin**(2.-p))


	#scale = 4.8*DL/(1+z)**2 # DL in Gpc
	#coreshiftmeas in mas GHz
	coreshiftmeas = coreshiftmeas/scale
	rcore = scale*coreshiftmeas/np.sin(viewing)*freq**(-1./kr) 

	#magnetization Sigma
	if p == 2:
		print 'p=2'
		F_Sigma = f_p/(4*np.pi*np.log(gammamax/gammamin))
		Sigma = 1.58*10**(-5)*2*opening*gamma**2*delta**6/(np.sin(viewing)*(1.+z)**7.)*F_Sigma/f_p*rcore*freq*Tb**(-8)

		nrad = 4*10**4*gamma*np.sin(viewing)*(1.+z)**5/(2*opening*delta**4)*f_p*rcore**(-1)*freq*Tb**4
		#Sigma = gamma*B**2*f_p/(4*np.pi*m*c**2*nrad*np.log(gammamax/gammamin))
		

	else:
		F_Sigma = f_p*(2.-p)/(4*np.pi*(gammamax**(2.-p)-gammamin**(2.-p)))
		C_Sigma = F_Sigma/f_p*c_alpha/np.sqrt(5.*(4.+p))*(2*np.pi)**2*(2.8*(1.5)**((p-1.)/2.)*a_alpha/c_alpha)**(p+6.)
		Sigma = 4.1*10**3*(1.7*10**2)**(-p)*C_Sigma*2*opening*gamma**2/(delta*np.sin(viewing))*(delta/(1.+z))**(p+5.)*rcore*freq*Tb**(-(p+6.))

		Cn = f_p*np.sqrt(5*(p+4.))/c_alpha*(2.8*(1.5)**((p-1.)/2.)*a_alpha/c_alpha)**(-(p+2.))
		nrad = 1.1*10**(-3)*(1.7*10**2)**p*Cn*gamma*np.sin(viewing)*(1.+z)**(p+3.)/(2.*opening*delta**(p+2.))*rcore**(-1.)*freq*Tb**(p+2.)  #in cm**(-3)
		#Sigma = gamma*(2.-p)*B**2*f_p/(4*np.pi*m*c**2*nrad*(gammamax**(2.-p)-gammamin**(2.-p)))


	return Sigma,nrad


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

"""for i in xrange(0,len(x_l)):
		pts_l.append(get_ellipse_coords(a=size_l[i], b=size_l[i], x=x_l[i],y=y_l[i], angle=0))
		pt_l.append(get_ellipse_coords(a=0.01, b=0.01, x=x_l[i],y=y_l[i], angle=0))
		#lines axis ellipses      
		x_el_l.append(ellipse_axis(x=x_l[i],y=y_l[i],s=size_l[i])[0])
		y_el_l.append(ellipse_axis(x=x_l[i],y=y_l[i],s=size_l[i])[1])
		x_elH_l.append(np.linspace(x_l[i],x_l[i],50))
"""

def ellipse_axis_lines(x,y,size):

	pts_arr=[]
	pt_arr=[]
	x_el_arr=[]
	x_elH_arr=[]
	y_el_arr=[]
	y_elH_arr=[]

	n = len(x)
	pts, pt = [], []
	x_el, y_el = [], []
	x_elH, y_elH = [], []
	for k in xrange(0,n):
		pts.append(get_ellipse_coords(a=size[k], b=size[k], x=x[k],y=y[k], angle=0))
		pt.append(get_ellipse_coords(a=0.01, b=0.01, x=x[k],y=y[k], angle=0))
		#lines axis ellipses      
		x_el.append(ellipse_axis(x=x[k],y=y[k],s=size[k])[0])
		y_el.append(ellipse_axis(x=x[k],y=y[k],s=size[k])[1])
		x_elH.append(np.linspace(x[k],x[k],50))
		y_elH.append(np.linspace(y[k],y[k],50))

	return pts,pt,x_el,y_el,x_elH,y_elH

def plot_components(pts,x_el,x_elH,y_elH,y_el):
    for i in xrange(0,len(x_el)):
        plt.plot(pts[i][:,0], pts[i][:,1], color='blue',linewidth=4)
        plt.plot(x_el[i], y_elH[i], color='blue',linewidth=4) 
        plt.plot(x_elH[i], y_el[i], color='blue',linewidth=4) 
    #fill(pts3[:,0], pts3[:,1], alpha=0.5, facecolor='blue', edgecolor

def plot_maps(realDAT,ext,first_contour):
	levels = first_contour*np.array([-1., 1., 1.41,2.,2.83,4.,5.66,8.,11.3,16.,
                                22.6,32.,45.3,64.,90.5,128.,181.,256.,362.,512.,
                                724.,1020.,1450.,2050.])
	cset = plt.contour(realDAT, levels, inline=1,
	                  colors=['grey']
	                  ,extent=[ext[0],ext[1],ext[2],ext[3]], aspect=1.0)
	plt.axis('scaled')
	plt.xlabel('Right Ascension [mas]',fontsize=14)
	plt.ylabel('Declination [mas]',fontsize=14)

#file, list with all the modelfit files
def read_modfile(file1,beam,z,freq,errors,errfile):

	nfiles = len(file1)

	r_arr = []
	errr_arr = [] #np.array([0.]*nfiles)
	psi_arr =  []
	errpsi_arr =  []
	size_arr =  []
	errsize_arr =  []
	flux_arr =  []
	errflux_arr =  []
	tb1_arr =  []
	errtb1_arr =  []
	tb2_arr =  []
	errtb2_arr =  []
	tb1_arrPlanck =  []
	errtb1_arrPlanck =  []
	tb2_arrPlanck =  []
	errtb2_arrPlanck =  []
	dlim1_arr = []
	dlim2_arr = []

	cte = 1.22*10**(12)

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
		tbNew1=np.array([0.]*n)
		errtbNew1=np.array([0.]*n)
		tbNew2=np.array([0.]*n)
		errtbNew2=np.array([0.]*n)
		tbNew1Planck=np.array([0.]*n)
		errtbNew1Planck=np.array([0.]*n)
		tbNew2Planck=np.array([0.]*n)
		errtbNew2Planck=np.array([0.]*n)
		dlim1=np.array([0.]*n)
		dlim2=np.array([0.]*n)
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
			fluxValue = (float(split_f[0][:-1]))
        		flux[i] = fluxValue
        		r[i] = (float(split_f[1][:-1]))
        		psi[i] = (float(split_f[2][:-1])*np.pi/180.)
			sizeValue = (float(split_f[3][:-1])/2.)
        		size[i] = sizeValue
        		tb[i] = (float(split_f[7])) 
			#tbNew[i] = cte*(1.+z)/(freq**2)*fluxValue*10**3/(sizeValue**2)
			#print flux[i], size[i]
			#print tbNew[i]/10**(12), '10**(12)'

		if errors == True:
			temp_file2=open(errfile,mode='r')
			temp_file2.readline()
			temp_file2.readline()
			#for i in xrange(0,ntot):
			#	temp_file2.readline()
			for i in xrange(0,n):		
				split_f = temp_file2.readline().split()
				#print split_f
				fluxpeak[i] = (float(split_f[2][:-1]))	
				rms[i] = (float(split_f[1][:-1]))			

			for i in xrange(0,n):
				errflux[i] = rms[i]
				snr = flux[i]/rms[i]#[k][i] #change to flux_peak

				dlim1[i] = 4/np.pi*np.sqrt(np.pi*np.log(2)*beam**2*np.log((snr+1.)/(snr))) #np.log((snr+1.)/(snr))) 4/np.pi*beam
				dlim2[i] = 4/np.pi*np.sqrt(np.pi*np.log(2)*beam**2*np.log((snr)/(snr-1.))) #np.log((snr+1.)/(snr))) 4/np.pi*beam
				if size[i] > beam:
					ddec=np.sqrt(size[i]**2-beam**2)
				else:
					ddec=0.

				y=[dlim1[0],ddec]
				dg=np.max(y)   
				err_size = rms[i]*dlim1[i]/fluxpeak[i]
				err_r = err_size/2.
				if r[i] > 0.:
					err_psi = np.real(math.atan(err_r*180./(np.pi*r[i])))
				else:
					err_psi = 1./5*beam
				
				if err_size < 1./5*beam:	
					errsize[i] = 1./5*beam
				else:
					errsize[i] = (err_size)
				if err_r < 1./5*beam:	
					errr[i] =  1./5*beam
					if errr[i] < 1./2*size[i]:
						errr[i] = 1./2*size[i]
				else:
					errr[i] = (err_r)
				errpsi[i] = (err_psi)	
				if dlim1[i] > size[i]:
					tbvalues = Tbformula(freq,dlim1[i],flux[i],z)
					tbNew1[i] = tbvalues[1]#cte*(1.+z)/(freq**2)*flux[i]/(dlim1[i]**2)
					errtbNew1[i] = tbNew1[i]*np.sqrt((rms[i]/fluxpeak[i])**2)
					tbNew1Planck[i] = tbvalues[0]#cte*(1.+z)/(freq**2)*flux[i]/(dlim1[i]**2)
					errtbNew1Planck[i] = tbNew1Planck[i]*np.sqrt((rms[i]/fluxpeak[i])**2)
				else:
					tbvalues = Tbformula(freq,size[i],flux[i],z)
					tbNew1[i] = tbvalues[1]#cte*(1.+z)/(freq**2)*flux[i]/(dlim1[i]**2)
					errtbNew1[i] = tbNew1[i]*np.sqrt((rms[i]/fluxpeak[i])**2)
					tbNew1Planck[i] = tbvalues[0]#cte*(1.+z)/(freq**2)*flux[i]/(dlim1[i]**2)
					errtbNew1Planck[i] = tbNew1Planck[i]*np.sqrt((rms[i]/fluxpeak[i])**2)

				if dlim2[i] > size[i]:
					tbvalues = Tbformula(freq,dlim2[i],flux[i],z)
					tbNew2[i] = tbvalues[1]#cte*(1.+z)/(freq**2)*flux[i]/(dlim1[i]**2)
					errtbNew2[i] = tbNew2[i]*np.sqrt((rms[i]/fluxpeak[i])**2)
					tbNew2Planck[i] = tbvalues[0]#cte*(1.+z)/(freq**2)*flux[i]/(dlim1[i]**2)
					errtbNew2Planck[i] = tbNew2Planck[i]*np.sqrt((rms[i]/fluxpeak[i])**2)
				else:
					tbvalues = Tbformula(freq,size[i],flux[i],z)
					tbNew2[i] = tbvalues[1]#cte*(1.+z)/(freq**2)*flux[i]/(dlim1[i]**2)
					errtbNew2[i] = tbNew2[i]*np.sqrt((rms[i]/fluxpeak[i])**2)
					tbNew2Planck[i] = tbvalues[0]#cte*(1.+z)/(freq**2)*flux[i]/(dlim1[i]**2)
					errtbNew2Planck[i] = tbNew2Planck[i]*np.sqrt((rms[i]/fluxpeak[i])**2)


				#print tbNew[i]/10**(12), '10**(12)'
		elif errors == 'Done':
			print 'done'
		else:
			for i in xrange(0,n):
				errflux[i] = 0.1*flux[i]
				errr[i] = 1./5*beam
				if errr[i] < 1/2*size[i]:
					errr[i] = 1/2*size[i]
				errpsi[i] = 0.	
				errsize[i] = 2./5*beam			

		r_arr = r
		errr_arr = errr
		psi_arr = psi
		errpsi_arr = errpsi
		size_arr = size
		errsize_arr = errsize
		flux_arr = flux
		errflux_arr = errflux
		if errors == True:
			dlim1_arr.append(dlim1)
			tb1_arr.append(tbNew1)
			errtb1_arr.append(errtbNew1)
			dlim2_arr.append(dlim1)
			tb2_arr.append(tbNew2)
			errtb2_arr.append(errtbNew2)
			tb1_arrPlanck.append(tbNew1Planck)
			errtb1_arrPlanck.append(errtbNew1Planck)
			tb2_arrPlanck.append(tbNew2Planck)
			errtb2_arrPlanck.append(errtbNew2Planck)
		else:
			tb1_arr.append(tb)
			errtb1_arr.append(tb)
			tb2_arr.append(tb)
			errtb2_arr.append(tb)
			tb1_arrPlanck.append(tb)
			errtb1_arrPlanck.append(tb)
			tb2_arrPlanck.append(tb)
			errtb2_arrPlanck.append(tb)
	
		ntot = n + ntot + 1

	genenerate_tabus(r_arr,errr_arr,psi_arr,errpsi_arr,size_arr,errsize_arr,flux_arr,errflux_arr,tb1_arr[0],errtb1_arr[0],dlim1_arr[0],tb2_arr[0],errtb2_arr[0],dlim2_arr[0],tb1_arrPlanck[0],errtb1_arrPlanck[0],tb2_arrPlanck[0],errtb2_arrPlanck[0],freq)

	return r_arr,errr_arr,psi_arr,errpsi_arr,size_arr,errsize_arr,flux_arr,errflux_arr,tb1_arr,errtb1_arr,dlim1_arr,tb2_arr,errtb2_arr,dlim2_arr

def genenerate_tabus(r_arr,errr_arr,psi_arr,errpsi_arr,size_arr,errsize_arr,flux_arr,errflux_arr,tb1_arr,errtb1_arr,dlim1_arr,tb2_arr,errtb2_arr,dlim2_arr,tb1_arrPlanck,errtb1_arrPlanck,tb2_arrPlanck,errtb2_arrPlanck,freq):
	geometry_options = {
		"landscape": False,
		"margin": "0.5in",
		"headheight": "20pt",
		"headsep": "10pt",
		"includeheadfoot": True
	}
	doc = Document(page_numbers=True, geometry_options=geometry_options)

	# Generate data table with 'tight' columns
	table1 = Tabular('c|ccccccccccccccccccc')
	table1.add_row(["Component", "S_y","","", "r","","", "psi","","","x","","","y","","", "theta","","",  "theta_lim"], mapper=[bold])
	table1.add_hline()
	table1.add_hline()

	x_arr,errx_arr,y_arr,erry_arr = x_y(r_arr,errr_arr,psi_arr,errpsi_arr,True)

	for i in range(len(r_arr)):
		flux, errflux = roundErrorVariable(flux_arr[i],errflux_arr[i])
		r, errr = roundErrorVariable(r_arr[i],errr_arr[i])
		psi, errpsi = roundErrorVariable(psi_arr[i]*180./np.pi,errpsi_arr[i]*180./np.pi)
		size, errsize = roundErrorVariable(size_arr[i],errsize_arr[i])
		x, errx = roundErrorVariable(x_arr[i],errx_arr[i])
		y, erry = roundErrorVariable(y_arr[i],erry_arr[i])

		row = [i, flux,"%0.0f" %(0),errflux, r,"%0.0f" %(0),errr,psi,"%0.0f" %(0),errpsi,x,"%0.0f" %(0),errx,y,"%0.0f" %(0),erry,size,"%0.0f" %(0),errsize,"%0.3f" %(dlim1_arr[i])]
		table1.add_row(row)


	table2 = Tabular('c|cccccccccccc')
	table2.add_row(["Component", "Tb_1", "","","Tb_1Planck", "","", "Tb_2", "","","Tb_2Planck", "",""], mapper=[bold])
	table2.add_hline()
	table2.add_hline()


	for i in range(len(r_arr)):
	    row = [i, "%0.3f" %(tb1_arr[i]/10**(11)),"%0.0f" %(0),"%0.3f" %(errtb1_arr[i]/10**(11)),"%0.3f" %(tb1_arrPlanck[i]/10**(11)),"%0.0f" %(0),"%0.3f" %(errtb1_arrPlanck[i]/10**(11)),"%0.3f" %(tb2_arr[i]/10**(11)),"%0.0f" %(0), "%0.3f" %(errtb2_arr[i]/10**(11)),"%0.3f" %(tb2_arrPlanck[i]/10**(11)),"%0.0f" %(0),"%0.3f" %(errtb2_arrPlanck[i]/10**(11))]
	    table2.add_row(row)

	doc.append(table1)
	doc.append(NewPage())
	doc.append(table2)
	doc.append(NewPage())

	doc.generate_pdf("Tb"+str(round(freq)),clean_tex=False)


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

