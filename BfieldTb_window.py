import threading, time
import warnings
import sys
import sip
import codecs
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from math import *
from functools import *
import numpy as np
import astropy.io.fits as pf
from pylab import *
import pickle
import iminuit, probfit
from functions_TbB import Bfield_Tb,constants,Magnetization,get_ellipse_coords,ellipse_axis,ellipse_axis_lines,plot_components,plot_maps,read_modfile,x_y 
from scipy.optimize import curve_fit
from functions_conv import order_by_nu, read_conv_params
from functions_alignComp import natural_keys
from functions_align import find_same_beam,beam_array,check_map_params, cuttingMAP,cross_correlation_shifts_FITS, checking_shift
from functions_turnover import cuttingTURN, synchrotron, synchrotron_v1, guesses_turnover,guesses_turnoverPoint, guesses_PL, powerLaw, powerLawPlot
from functions_Bfield import searchNEDnoGUI
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, get_ellipse_coords, Annotate
import os,glob
import subprocess as sub
from astropy.nddata import Cutout2D
from correlate2d import *
import scipy.special as scp
#from fast_ftts import *

class BfieldTbWindow(QWidget):

    def __init__(self,*args):

        QWidget.__init__(self)

        layout = QGridLayout()


	#(alpha0,gammamin,gammamax,B,scale,viewing,z,gamma,delta,kr,freq,Tb,a_alpha,b_alpha,c_alpha)

	self.ifhdu = False
	self.Tball = []
	self.Tballfreq = []
	self.errTballfreq = []
	self.rallfreq = []
	self.Tbcore = np.zeros(len(needed_param.freq))
	self.fluxcore = np.zeros(len(needed_param.freq))
	self.sizecore = np.zeros(len(needed_param.freq))
	self.Btb = np.zeros(len(needed_param.freq))
	self.nrad = np.zeros(len(needed_param.freq))
	self.Sigma = np.zeros(len(needed_param.freq))
	self.coreshiftmeas = []
	self.a_alpha = 0.
	self.b_alpha = 0.
	self.c_alpha = 0.
	self.kr = 0.
	self.freq = 0.
	self.x = 0.
	self.y = 0.
	self.zvalue = 0.
	self.DLvalue = 0.
	self.scalevalue = 0.
	self.gammaminvalue = 0.
	self.gammamaxvalue = 0.
	self.alpha0value = 0.
	self.gammaavalue = 0.
	self.deltavalue = 0.
	self.beta = 0.
	self.openingangle = 0.
	self.viewingangle = 0.	

	self.fits1 = ''
	self.models1 = ''
	self.models1errors = ''
	self.freq1 = 0.
	self.rms1 = 0.
	self.fits2 = ''
	self.models2 = ''
	self.bmaj_files = 0.
	self.bmin_files = 0.
	self.bpa_files = 0.
	self.beam_files = 0.
	self.cells_files = 0.
	self.cent_mapx = 0.
	self.cent_mapy = 0.
	self.mapsize_files = 0.
	self.extmas=[] 
	self.extmasnew=[] 
	self.image1 = np.asarray([])
	self.image2 = np.asarray([])


	self.labelTEXT = QLabel()
	self.labelTEXT.setText('Select the frequency for which you want to calculate the core Bfield:')
	self.labelTEXT.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

	self.checks = []

        for i in xrange(0,len(needed_param.freq)):
            c = QCheckBox('%s GHz' % ('%1.2f' % (needed_param.freq[i])),self)  #add the MHz thing
	    c.setFixedSize(100,25)
            layout.addWidget(c,2,i+1)
            self.checks.append(c)


	self.labelempty = QLabel()
	self.labelempty.setFixedSize(25,25)
	self.labelempty2 = QLabel()


	self.labelsigmaCut = QLabel()
	self.labelsigmaCut.setText("SigmaCut: ")
	self.labelsigmaCut.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.SigmaCut = QLineEdit()
	self.SigmaCut.setValidator(QDoubleValidator())
	self.SigmaCut.textChanged.connect(self.check_state)
	self.SigmaCut.textChanged.emit(self.SigmaCut.text())
	self.SigmaCut.setFixedSize(100,25)

	self.labelDL = QLabel()
	self.labelDL.setText("DL   : ")
	self.labelDL.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.DL = QLineEdit()
	self.DL.setValidator(QDoubleValidator())
	self.DL.textChanged.connect(self.check_state)
	self.DL.textChanged.emit(self.DL.text())
	self.DL.setFixedSize(100,25)
	self.labelDLunit = QLabel()
	self.labelDLunit.setText("Mpc")
	self.labelDLunit.setAlignment(Qt.AlignLeft | Qt.AlignVCenter) 
	self.labelDLunit.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

	self.labelz = QLabel()
	self.labelz.setText("z     : ")
	self.labelz.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.z = QLineEdit()
	self.z.setValidator(QDoubleValidator())
	self.z.textChanged.connect(self.check_state)
	self.z.textChanged.emit(self.z.text())
	self.z.setFixedSize(100,25)

	self.labelscale = QLabel()
	self.labelscale.setText("Scale : ")
	self.labelscale.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.scale = QLineEdit()
	self.scale.setValidator(QDoubleValidator())
	self.scale.textChanged.connect(self.check_state)
	self.scale.textChanged.emit(self.scale.text())
	self.scale.setFixedSize(100,25)
	self.labelscaleunit = QLabel()
	self.labelscaleunit.setText("pc/mas")
	self.labelscaleunit.setAlignment(Qt.AlignLeft | Qt.AlignVCenter) 
	self.labelscaleunit.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')


	self.labelgammamin = QLabel()
	self.labelgammamin.setText("Gamma min: ")
	self.labelgammamin.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.gammamin = QLineEdit()
	self.gammamin.setValidator(QDoubleValidator())
	self.gammamin.textChanged.connect(self.check_state)
	self.gammamin.textChanged.emit(self.gammamin.text())
	self.gammamin.setFixedSize(100,25)

	self.labelgammamax = QLabel()
	self.labelgammamax.setText("Gamma max: ")
	self.labelgammamax.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.gammamax = QLineEdit()
	self.gammamax.setValidator(QDoubleValidator())
	self.gammamax.textChanged.connect(self.check_state)
	self.gammamax.textChanged.emit(self.gammamax.text())
	self.gammamax.setFixedSize(100,25)

	self.labelalpha0 = QLabel()
	self.labelalpha0.setText("Alpha0 : ")
	self.labelalpha0.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.alpha0 = QLineEdit()
	self.alpha0.setValidator(QDoubleValidator())
	self.alpha0.textChanged.connect(self.check_state)
	self.alpha0.textChanged.emit(self.alpha0.text())
	self.alpha0.setFixedSize(100,25)


	self.labelgamma = QLabel()
	self.labelgamma.setText("Gamma : ")
	self.labelgamma.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.gamma = QLineEdit()
	self.gamma.setValidator(QDoubleValidator())
	self.gamma.textChanged.connect(self.check_state)
	self.gamma.textChanged.emit(self.gamma.text())
	self.gamma.setFixedSize(100,25)

	self.labelviewing = QLabel()
	self.labelviewing.setText("Viewing angle : ")
	self.labelviewing.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.viewing = QLineEdit()
	self.viewing.setValidator(QDoubleValidator())
	self.viewing.textChanged.connect(self.check_state)
	self.viewing.textChanged.emit(self.viewing.text())
	self.viewing.setFixedSize(100,25)

	self.labelopening = QLabel()
	self.labelopening.setText("Opening angle : ")
	self.labelopening.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.opening = QLineEdit()
	self.opening.setValidator(QDoubleValidator())
	self.opening.textChanged.connect(self.check_state)
	self.opening.textChanged.emit(self.opening.text())
	self.opening.setFixedSize(100,25)

        self.Tbbutton = QPushButton("&Brightness \n temperature")
	#self.Bfieldbutton.setFixedSize(100,25)
        self.Tbbutton.clicked.connect(lambda: self.Tbcalculation())
	self.Tbbutton.setAutoDefault(True)

        self.Bfieldbutton = QPushButton("&Magnetic \n field")
	#self.Bfieldbutton.setFixedSize(100,25)
        self.Bfieldbutton.clicked.connect(lambda: self.Bfieldcalculation())
	self.Bfieldbutton.setAutoDefault(True)

        self.Sigmabutton = QPushButton("&Magnetization \n and  \n Particle density")
	#self.Nbutton.setFixedSize(100,25)
        self.Sigmabutton.clicked.connect(lambda: self.SigmaNcalculation())
	self.Sigmabutton.setAutoDefault(True)


	temp = searchNEDnoGUI(needed_param.source_name) 
	#omegarnu coreshift measurement

	self.DL.setText('%1.3f' % (temp[0]))
	self.z.setText('%1.3f' % (temp[1]))
	self.scale.setText('%1.3f' % (temp[2]))

	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 0)	

	layout.addWidget(self.labelTEXT, 1, 1,1,5)

	layout.addWidget(self.labelsigmaCut, 4,1)
	layout.addWidget(self.labelDL, 6,1)
	layout.addWidget(self.labelz, 7,1)
	layout.addWidget(self.labelscale, 8,1)
	layout.addWidget(self.labelgammamin, 10,1)
	layout.addWidget(self.labelgammamax, 11,1)
	layout.addWidget(self.labelalpha0, 6,4)
	layout.addWidget(self.labelgamma, 7,4)
	layout.addWidget(self.labelviewing, 9,4)
	layout.addWidget(self.labelopening, 10,4)
	layout.addWidget(self.SigmaCut, 4,2)
	layout.addWidget(self.DL, 6,2)
	layout.addWidget(self.z, 7,2)
	layout.addWidget(self.scale, 8,2)
	layout.addWidget(self.gammamin, 10,2)
	layout.addWidget(self.gammamax, 11,2)
	layout.addWidget(self.alpha0, 6,5)
	layout.addWidget(self.gamma, 7,5)
	layout.addWidget(self.viewing, 9,5)
	layout.addWidget(self.opening, 10,5)
	layout.addWidget(self.labelDLunit, 6,3)
	layout.addWidget(self.labelscaleunit, 8,3)


	#assigning buddies for elements in the layout that are tied together
	self.labelsigmaCut.setBuddy(self.SigmaCut)
	self.labelDL.setBuddy(self.DL)
	self.labelz.setBuddy(self.z)
	self.labelscale.setBuddy(self.scale)
	self.labelgammamin.setBuddy(self.gammamin)
	self.labelgammamax.setBuddy(self.gammamax)
	self.labelalpha0.setBuddy(self.alpha0)
	self.labelgamma.setBuddy(self.gamma)
	self.labelviewing.setBuddy(self.viewing)
	self.labelopening.setBuddy(self.opening)

	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 6)	

	layout.addWidget(self.Tbbutton, 5,7,2,1)
	layout.addWidget(self.Bfieldbutton, 7,7,2,1)
	layout.addWidget(self.Sigmabutton, 9,7,3,1)

	self.Tbbutton.setEnabled(False)
	self.Bfieldbutton.setEnabled(False)
	self.Sigmabutton.setEnabled(False)

	for i in xrange(len(self.checks)):
		self.checks[i].toggled.connect(lambda checked: self.checksState())


        self.setLayout(layout)

	#put the window in the center of the desktop
	qr = self.frameGeometry()
	cp = QDesktopWidget().availableGeometry().center()
	qr.moveCenter(cp)
	self.move(qr.topLeft())

    #function to check if the arguments given in the text boxes are fine
    #if the value is not the kind of parameter wanted, 
    #for example, if the box requires a double and you give an integrer, it get red
    #if the input value is still the kind of parameter wanted, but outside a range, it gets yellow
    #if the input value is fine, it gets green
    def check_state(self,*args,**kwargs):
		sender = self.sender()
		validator = sender.validator()
		state = validator.validate(sender.text(),0)[0]
		if state == QValidator.Acceptable:
			color = '#c4df9b' #green
		elif state == QValidator.Intermediate:
			color = '#fff79a' #yellow
		else:
			color = '#f6989d' #red
		sender.setStyleSheet('QLineEdit { background-color: %s }' %color)

    #disables buttons
    def checksState(self):
	checked = []
	checked2 = []
	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			checked.append(i)

	if len(checked) == 1:
		self.Tbbutton.setEnabled(True)
		self.Bfieldbutton.setEnabled(True)
	else:
		self.Tbbutton.setEnabled(False)
		self.Bfieldbutton.setEnabled(False)

	if len(checked) == 1 and len(needed_param.coreshiftfile) == 1:
		self.Sigmabutton.setEnabled(True)

	else:
		self.Sigmabutton.setEnabled(False)


    def Tbcalculation(self):

	self.DLvalue = float(self.DL.text())
	self.zvalue = float(self.z.text())
	self.scalevalue = float(self.scale.text())

	a = 0
	for i in xrange(len(self.checks)):
		if self.checks[i].isChecked():
			#if a == 0:
			self.fits1 = needed_param.fits[i]
			self.models1 = needed_param.modelfit[i]
			self.models1errors = needed_param.modelfiterror[i]
			tempi = i
			#	a = 1
			#if a == 1:
			#	self.fits2 = needed_param.fits[i]

	#self.reading_rms_param()
	header1 = take_header(self.fits1,self.ifhdu)
	#header2 = take_header(self.fits2,self.ifhdu)
	map_data1 = read_map(self.fits1,self.ifhdu)		
	self.image1 = map_data1[0]
	#map_data2 = read_map(self.fits2,self.ifhdu)		
	#self.image2 = map_data2[0]

	#obtaining the beam and cell
	self.bmaj_files = header1[1]
	self.bmin_files = header1[2]
	self.bpa_files = header1[3]
	self.beam_files = header1[7]
	self.cells_files = header1[0]
	self.freq1 = header1[5]

	#obtaining map centers in pixels
	self.cent_mapx = map_data1[5]
	self.cent_mapy = map_data1[6]
	self.mapsize_files = 2*map_data1[7]

	#obtaining the four corners of the maps in mas
	x1 = map_data1[1]
	x2 = map_data1[2]
	y1 = map_data1[3]
	y2 = map_data1[4]

	self.extmas=[x1,x2,y1,y2] 
	
		
	self.sigma_cut = float(self.SigmaCut.text())

	if self.freq1 < 0.5:
		self.freq1name = str('%1.0f' %(self.freq1*1000))
		self.freq1unit = 'MHz'
	else:
		self.freq1name = str('%1.2f' %(self.freq1))
		self.freq1unit = 'GHz'

	res=open(needed_param.path+'/rms'+self.freq1name+'.p','rb')
	pick = pickle.load(res)
	res.close() 

	self.rms1 = pick	

	modelfitparameters = read_modfile([self.models1],self.beam_files,self.zvalue,self.freq1,True,self.models1errors)

	#r_arr,errr_arr,
	#2psi_arr,errpsi_arr,
	#4size_arr,errsize_arr,
	#6flux_arr,errflux_arr,
	#8tb1_arr,errtb1_arr,
	#dlim1_arr,
	#11tb2_arr,errtb2_arr,
	#dlim2_arr

	x_and_y = x_y(modelfitparameters[0],modelfitparameters[1],modelfitparameters[2],modelfitparameters[3],True)
	x, errx = np.asarray(x_and_y[0]), np.asarray(x_and_y[1])
	y, erry = np.asarray(x_and_y[2]), np.asarray(x_and_y[3])

	#for plotting the components in the map
	pts_arr=[]
	pt_arr=[]
	x_el_arr=[]
	x_elH_arr=[]
	y_el_arr=[]
	y_elH_arr=[]

	#r_arr 0,errr_arr 1 ,psi_arr 2,errpsi_arr 3,size_arr 4,errsize_arr 5,tb_arr 6,errtb_arr 7,flux_arr 8,errflux_arr 9,tbNew 10

	r = modelfitparameters[0]
	errr = modelfitparameters[1]
	psi = modelfitparameters[2]
	size = modelfitparameters[4]
	errsize = modelfitparameters[5]
	Tb = modelfitparameters[8][0]
	errTb = modelfitparameters[9][0]
	self.Tbcore[tempi] = Tb[0]
	self.sizecore[tempi] = modelfitparameters[4][0]
	self.fluxcore[tempi] = modelfitparameters[6][0]
	print 'Tb', self.Tbcore
	print 'size', self.sizecore
	print 'flux', self.fluxcore

	ellipse_plot = ellipse_axis_lines(x,y,modelfitparameters[4])
	pts_arr,pt_arr = ellipse_plot[0], ellipse_plot[1]
	x_el_arr,y_el_arr = ellipse_plot[2], ellipse_plot[3]
	x_elH_arr,y_elH_arr = ellipse_plot[4], ellipse_plot[5]  

	plt.figure(1)
	plot_components(pts_arr,x_el_arr,x_elH_arr,y_elH_arr,y_el_arr)
	plot_maps(self.image1,self.extmas,self.rms1*self.sigma_cut)
	plt.xlim(x1,x2)
	#plt.ylim(-1.5,1.5)
	#plt.savefig('1642CBAND.png', bbox_inches='tight')

	plt.show()

	"""limits = Annotate()
	plt.show()
		
	ext_new = []
	[self.limplot_x1,self.limplot_x2,self.limplot_y1,self.limplot_y2] = limits()

	self.extmasnew = [self.limplot_x1,self.limplot_x2,self.limplot_y2,self.limplot_y1]

	plt.figure(1)
	plot_components(pts_arr,x_el_arr,x_elH_arr,y_elH_arr,y_el_arr)
	plot_maps(self.image1,self.extmas,self.rms1*self.sigma_cut)
	plt.xlim(self.limplot_x1,self.limplot_x2)
	plt.ylim(self.limplot_y2,self.limplot_y1)
	plt.savefig('modelfit'+str('%1.1f' % (self.freq1))+'.png', bbox_inches='tight')

	plt.show()"""



	plt.close('all')

	plt.figure(2)
	plt.plot(needed_param.freq,self.Tbcore,'r.',markersize=12)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylabel(r'$T_b$ [K]')
	plt.xlabel(r'$\nu$ [GHz]')
	plt.savefig('logTb_ground.png')

	plt.figure(3)
	plt.plot(needed_param.freq,self.sizecore,'r.',markersize=12)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylabel(r'$\theta$ [mas]')
	plt.xlabel(r'$\nu$ [GHz]')
	plt.savefig('logTheta_ground.png')

	plt.figure(4)
	plt.plot(needed_param.freq,self.fluxcore,'r.',markersize=12)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylabel(r'$S_y$ [Jy]')
	plt.xlabel(r'$\nu$ [GHz]')
	plt.savefig('logS_ground.png')
	#plt.xlim(x1,x2)
	#plt.ylim(-1.5,1.5)
	#plt.savefig('1642CBAND.png', bbox_inches='tight')

	plt.close('all')

	plt.figure(2)
	plt.xscale('log')
	plt.yscale('log')
	plt.errorbar(r,Tb,yerr=0.434*errTb,fmt='ro',ecolor='r', capthick=2)
	plt.ylabel(r'$T_b$ [K]')
	plt.xlabel('r [mas]')
	plt.savefig('Tbfreq'+str('%1.1f' % (self.freq1))+'.png')

	plt.close('all')


	res=open(needed_param.path+'/pickle'+str('%1.1f' % (self.freq1))+'.p','wb')
	pickle.dump([r,Tb,0.434*errTb,psi,x,y,size,errsize],res)
	res.close()     


    def Bfieldcalculation(self):

	self.viewingangle = float(self.viewing.text())
	self.viewingangle = np.deg2rad(self.viewingangle)
	self.openingangle = np.deg2rad(self.openingangle)
	self.openingangle = float(self.opening.text())
	self.gammavalue = float(self.gamma.text())

	self.beta = np.sqrt(1. - 1./self.gammavalue**2)
	self.deltavalue = 1./(self.gammavalue*(1-self.beta*np.cos(self.viewingangle)))

	for i in xrange(len(self.checks)):
		if self.checks[i].isChecked():
			#if a == 0:
			self.fits1 = needed_param.fits[i]
			self.models1 = needed_param.modelfit[i]
			self.models1errors = needed_param.modelfiterror[i]
			tempi = i


	header1 = take_header(self.fits1,self.ifhdu)
	self.freq1 = header1[5]

	self.Btb[tempi] = Bfield_Tb(self.Tbcore[tempi]/10**(12),self.zvalue,self.gammavalue,self.deltavalue,self.freq1)

	print self.Btb

    def SigmaNcalculation(self):

	self.gammamaxvalue = float(self.gammamax.text())
	self.gammaminvalue = float(self.gammamin.text())
	self.alpha0value = float(self.alpha0.text())

	for i in xrange(len(self.checks)):
		if self.checks[i].isChecked():
			#if a == 0:
			self.fits1 = needed_param.fits[i]
			self.models1 = needed_param.modelfit[i]
			self.models1errors = needed_param.modelfiterror[i]
			tempi = i

	header1 = take_header(self.fits1,self.ifhdu)
	self.freq1 = header1[5]

	#kr=1
	#coreshiftmeas=48

	self.Sigma[tempi], self.nrad[tempi] = Magnetization(self.alpha0value,self.gammaminvalue,self.gammamaxvalue,self.Btb[tempi],self.scalevalue,self.viewingangle,self.openingangle,self.zvalue,self.gammavalue,self.deltavalue,1.,self.freq1,self.Tbcore[tempi]/10**(12),48.)

	print self.Sigma
	print self.nrad






class needed_param():

	path = os.getcwd()



	#store the uvf, mod and fits files of the maps in a list
	files = []
	for filename in sorted(glob.glob(path+'/UVF/*.uvf*')):   
		files.append(filename)     #for the moment as well, you can also read the file with header = pf.getheader(uvffile) and then freq = header['CRVAL4'], all that would be easier in general as i dont need to depend of having a similar modification date in all of them ----> leads to changing the function order_by_nu
			
	models = []
	for filename in sorted(glob.glob(path+'/MODELS/*.mod*')):   
		models.append(filename)  #for the moment, in the modelfit file it is posible to read the frequency, which will simplify how to get it in general
		
	fits = []
	for filename in sorted(glob.glob(path+'/FITS/*.fits*')):   
		fits.append(filename)  

	modelfit = []
	for filename in sorted(glob.glob(path+'/modelfit/*.mod*')):   
		modelfit.append(filename)  #for the moment, in the modelfit file it is posible to read the frequency, which will simplify how to get it in general

	modelfit.sort(key=natural_keys)

	modelfiterror = []
	for filename in sorted(glob.glob(path+'/modelfit/*.dat*')):   
		modelfiterror.append(filename)  #for the moment, in the modelfit file it is posible to read the frequency, which will simplify how to get it in general

	modelfiterror.sort(key=natural_keys)

	coreshiftfile = []
	for filename in sorted(glob.glob(path+'/coreshiftmeas/*.meas*')):   
		coreshiftfile.append(filename)   #for the moment, in the modelfit file it is posible to read the frequency, which will simplify how to get it in general

					
	#initialize arrays
	cell = np.array([0.]*len(fits))
	bmaj = np.array([0.]*len(fits))
	bmin = np.array([0.]*len(fits))
	bpa = np.array([0.]*len(fits))
	freq = np.array([0.]*len(fits))
	beam = np.array([0.]*len(fits))
	size_map = np.array([0.]*len(fits))
	size_map_y =np.array([0.]*len(fits))

	#order the the list by frequency and 
	#getting the corresponding values of the previous initialized arrays ordered by frequency 
	#(lower to higher)
	ordered_params = order_by_nu(files,models,fits,False)
	
	freq = ordered_params[0]
	files = ordered_params[8]
	models = ordered_params[9]
	fits = ordered_params[10]
	
	#source name
	header = take_header(fits[0],False)
	source_name = header[8]

def main():
	app = QApplication(sys.argv)

	w = BfieldcoreshiftWindow()
	w.show()

	app.exec_()

