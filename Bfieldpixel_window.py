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
from scipy.optimize import curve_fit
from functions_conv import order_by_nu, read_conv_params
from functions_align import find_same_beam,beam_array,check_map_params, cuttingMAP,cross_correlation_shifts_FITS, checking_shift
from functions_turnover import cuttingTURN, synchrotron, synchrotron_v1, guesses_turnover,guesses_turnoverPoint, guesses_PL, powerLaw, powerLawPlot
from functions_Bfield import searchNEDnoGUI, B_field, N_UeUb_sigma
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, get_ellipse_coords, Annotate
import os,glob
import subprocess as sub
from astropy.nddata import Cutout2D
from correlate2d import *
#from fast_ftts import *


def chi2_sync(S_m, v_m, alpha0,alphathick):
    return np.sum((synchrotron(x,S_m, v_m, alpha0,alphathick)-y)**2)

def chi2_syncv1(S_1, v_1, alpha0,alphathick): 
    return np.sum((synchrotron_v1(x,S_1, v_1, alpha0,alphathick)-y)**2)

def chi2_PL(cte,alpha0):
    return np.sum((powerLaw(x,cte,alpha0)-y)**2)


class popup_sourcename(QWidget):
	def __init__(self):
     	 	QWidget.__init__(self,None, Qt.WindowStaysOnTopHint)
	   	#self.layout = QGridLayout()
	   	self.layout = QGridLayout()


		self.sourceName = needed_param.source_name
		#layout.addStretch(1)
		#layout.addLayout(hbox)
		self.labelempty = QLabel()
		self.labelTEXT = QLabel()
		self.labelSOURCE = QLabel()
		self.labelTEXT2 = QLabel()
		self.labelTEXT.setText("Sorry the source was not found in NED by ")
		self.labelTEXT.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.labelSOURCE.setText(self.sourceName)
		self.labelSOURCE.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
		self.labelTEXT2.setText("Please give a more common name or set DL, z and the scale manually in the main window")
		self.labelTEXT2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

		self.layout.addWidget(self.labelempty,1,0)
		self.layout.addWidget(self.labelempty,2,0)
		self.layout.addWidget(self.labelTEXT,1,1,1,3)
		self.layout.addWidget(self.labelSOURCE,1,4)
		self.layout.addWidget(self.labelTEXT2,2,1,1,5)

		self.labelSourceName = QLabel()
		self.labelSourceName.setText("Source name    : ")
		self.labelSourceName.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.SourceName = QLineEdit()
		#self.SigmaCut.setValidator(QDoubleValidator())
		#self.SigmaCut.textChanged.connect(self.check_state)
		#self.SigmaCut.textChanged.emit(self.SigmaCut.text())
		self.SourceName.setFixedSize(100,25)


		self.layout.addWidget(self.labelSourceName,3,1)
		self.layout.addWidget(self.SourceName,3,2)

		self.labelSourceName.setBuddy(self.SourceName)

		self.getSource = QPushButton("&Get")
		self.getSource.setFixedSize(100,25)

		self.layout.addWidget(self.getSource,6,2)

		self.setLayout(self.layout)

		self.adjustSize()

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())



class Bfieldpixel(QWidget):

    def __init__(self,*args):

        QWidget.__init__(self)

        layout = QGridLayout()

	self.DLvalue = 0.
	self.zvalue = 0.
	self.scalevalue = 0.
	self.gammavalue = 0.
	self.deltavalue = 0.
	self.gammaminvalue = 0.
	self.gammamaxvalue = 0.

	self.vmall = np.asarray([])
	self.small = np.asarray([])
	self.alpha0all = np.asarray([])
	self.alphathickall = np.asarray([])
	self.Ball = np.asarray([])
	self.Kall = np.asarray([])
	self.Ueall = np.asarray([])
	self.Uball = np.asarray([])
	self.sigmaall = np.asarray([])

	self.ext = []
	self.cellsize = 0.

	self.labelempty = QLabel()
	self.labelempty.setFixedSize(25,25)
	self.labelempty2 = QLabel()

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


	self.labelgamma = QLabel()
	self.labelgamma.setText("Gamma : ")
	self.labelgamma.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.gamma = QLineEdit()
	self.gamma.setValidator(QDoubleValidator())
	self.gamma.textChanged.connect(self.check_state)
	self.gamma.textChanged.emit(self.gamma.text())
	self.gamma.setFixedSize(100,25)


	self.labeldelta = QLabel()
	self.labeldelta.setText("Delta : ")
	self.labeldelta.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.delta = QLineEdit()
	self.delta.setValidator(QDoubleValidator())
	self.delta.textChanged.connect(self.check_state)
	self.delta.textChanged.emit(self.delta.text())
	self.delta.setFixedSize(100,25)


	temp = searchNEDnoGUI(needed_param.source_name) 


	if temp[0] == 0:
		self.wi = popup_sourcename()
		self.wi.show()
		self.wi.getSource.clicked.connect(lambda: self.getParameters())

	else:
		self.DL.setText('%1.3f' % (temp[0]))
		self.z.setText('%1.3f' % (temp[1]))
		self.scale.setText('%1.3f' % (temp[2]))


        self.Bfieldbutton = QPushButton("&Magnetic \n field")
	#self.Bfieldbutton.setFixedSize(100,25)
        self.Bfieldbutton.clicked.connect(lambda: self.Bfieldcalculation())
	self.Bfieldbutton.setAutoDefault(True)

        self.Nbutton = QPushButton("&Particle \n density")
	#self.Nbutton.setFixedSize(100,25)
        self.Nbutton.clicked.connect(lambda: self.Ncalculation())
	self.Nbutton.setAutoDefault(True)

	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 0)	
	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 4)	

	layout.addWidget(self.labelDL, 4,1)
	layout.addWidget(self.labelz, 5,1)
	layout.addWidget(self.labelscale, 6,1)
	layout.addWidget(self.labelgammamin, 8,1)
	layout.addWidget(self.labelgammamax, 9,1)
	layout.addWidget(self.labelgamma, 4,5)
	layout.addWidget(self.labeldelta, 5,5)
	layout.addWidget(self.DL, 4,2)
	layout.addWidget(self.z, 5,2)
	layout.addWidget(self.scale, 6,2)
	layout.addWidget(self.gammamin, 8,2)
	layout.addWidget(self.gammamax, 9,2)
	layout.addWidget(self.gamma, 4,6)
	layout.addWidget(self.delta, 5,6)
	layout.addWidget(self.labelDLunit, 4,3)
	layout.addWidget(self.labelscaleunit, 6,3)

	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 7)	

	#assigning buddies for elements in the layout that are tied together
	self.labelDL.setBuddy(self.DL)
	self.labelz.setBuddy(self.z)
	self.labelscale.setBuddy(self.scale)
	self.labelgammamin.setBuddy(self.gammamin)
	self.labelgammamax.setBuddy(self.gammamax)
	self.labelgamma.setBuddy(self.gamma)
	self.labeldelta.setBuddy(self.delta)

	layout.addWidget(self.Bfieldbutton, 5,8,2,1)
	layout.addWidget(self.Nbutton, 7,8,2,1)

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


    def Bfieldcalculation(self):

	self.DLvalue = float(self.DL.text())
	self.zvalue = float(self.z.text())
	self.scalevalue = float(self.scale.text())
	self.gammavalue = float(self.gamma.text())
	self.deltavalue = float(self.delta.text())

	filefreq = 'bo/shifted15.fits'
	header = take_header(filefreq,False)
	self.cellsize = header[0]
	mapdata = read_map(filefreq,False)
	realDAT = mapdata[0]
	self.ext = [mapdata[1],mapdata[2],mapdata[3],mapdata[4]]

	res = open('turnoverdata.p','rb')
	pick = pickle.load(res)
	res.close()
	self.vmall = pick[0]
	self.small = pick[1]#*1.4213481404770085
	self.alpha0all = pick[2]
	self.alphathickall = pick[3]
	shapee = np.shape(self.alpha0all)
	self.Ball = np.zeros(shapee)
	self.Kall = np.zeros(shapee)


	self.Ball[:] = np.nan
	self.Kall[:] = np.nan


	#cellsize *2 for getting the diameter
	for i in xrange (0,shapee[0]):
	    for j in xrange(0,shapee[1]):
		if np.isnan(self.small[i][j]) == False:
		    tempBK = B_field(self.alpha0all[i][j],self.alphathickall[i][j],self.DLvalue,self.zvalue,self.scalevalue,self.cellsize*2*self.gammavalue,self.deltavalue,self.vmall[i][j],self.small[i][j])
		    self.Ball[i][j] = tempBK[0]
		    self.Kall[i][j] = tempBK[1]

	plt.figure(1)
	cset = plt.contour(realDAT,0.008*np.array([2.,4.,16.,64.,256.,1020.,2050.]),inline=1,colors=['grey'],aspect=1.0,extent=self.ext)
	plt.imshow(self.Ball,origin='bottom',extent=self.ext,vmin=0,vmax=0.1)
	plt.colorbar()
	plt.show()

    def Ncalculation(self):

	self.gammaminvalue = float(self.gammamin.text())
	self.gammamaxvalue = float(self.gammamax.text())

	shapee = np.shape(self.Ball)

	self.Nall = np.zeros(shapee)
	self.Ueall = np.zeros(shapee)
	self.Uball = np.zeros(shapee)
	self.sigmaall = np.zeros(shapee)

	self.Nall[:] = np.nan
	self.Ueall[:] = np.nan
	self.Uball[:] = np.nan
	self.sigmaall[:] = np.nan


	for i in xrange (0,shapee[0]):
	    for j in xrange(0,shapee[1]):
		if np.isnan(self.Ball[i][j]) == False:
		    tempBUsigma = N_UeUb_sigma(self.Ball[i][j],self.Kall[i][j],self.gammaminvalue,self.gammamaxvalue,self.alpha0all[i][j])
		    self.Nall[i][j] = tempBUsigma[0]
		    self.Ueall[i][j] = tempBUsigma[1]
		    self.Uball[i][j] = tempBUsigma[2]
		    self.sigmaall[i][j] = tempBUsigma[3]
		    print self.Nall[i][j],self.sigmaall[i][j]

	
	#N_UeUb_sigma

    def getParameters(self):
		if self.wi.SourceName.text() != '':
			self.wi.sourceName = str(self.wi.SourceName.text())

			temp = searchNEDnoGUI(self.wi.sourceName)

			self.DL.setText('%1.3f' % (temp[0]))
			self.z.setText('%1.3f' % (temp[1]))
			self.scale.setText('%1.3f' % (temp[2]))

			self.wi.close()





class needed_param():

	path = os.getcwd()

	if not os.path.exists('Plot_fitted_synchrotron'):
		os.makedirs('Plot_fitted_synchrotron')
	if not os.path.exists('Plot_fitted_PL'):
		os.makedirs('Plot_fitted_PL')
	if not os.path.exists('Plot_fitted_synchrotronExtrapolated'):
		os.makedirs('Plot_fitted_synchrotronExtrapolated')


	#store the uvf, mod and fits files of the maps in a list
	files = []
	for filename in sorted(glob.glob(path+'/UVF/*.uvf*')):   
		files.append(filename)   
			
	models = []
	for filename in sorted(glob.glob(path+'/MODELS/*.mod*')):   
		models.append(filename)  
		
	fits = []
	for filename in sorted(glob.glob(path+'/FITS/*.fits*')):   
		fits.append(filename)  
					
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

	w = Bfieldpixel()
	w.show()

	app.exec_()

#main()

#x=self.freqs_conv 
#y=spec

#def chi2_syncv1(S_1, v_1, alpha0,alphathick): 
#    return np.sum((synchrotron_v1(x,S_1, v_1, alpha0,alphathick)-y)**2)


#m=iminuit.Minuit(chi2_syncv1,
  #S_1 = guess_S1, v_1 = guess_v1, alpha0 = guess_alpha0, alphathick = guess_alphathick)

#pfinal, covx = curve_fit(synchrotron_v1, self.freqs_conv, spec, parameters,maxfev=5000)
#fitted_spec = synchrotron(freq, pfinal[0], pfinal[1], pfinal[2])
#guess_alphathick = pfinal[3] #m.values['alphathick']
#guess_alpha0 = pfinal[2]#m.values['alpha0']
#guess_Sm = pfinal[0]#m.values['S_1']
#guess_vm = pfinal[1]#m.values['v_1']

#x = self.freqs_conv
#y = spec


