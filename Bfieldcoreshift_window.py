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
from functions_Bfield import searchNEDnoGUI
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, get_ellipse_coords, Annotate
import os,glob
import subprocess as sub
from astropy.nddata import Cutout2D
from correlate2d import *
import scipy.special as scp
#from fast_ftts import *


def rcoreFunct(freqs,A,kr):
	deltar = A*(freqs**(-1./kr)-43.14**(-1./kr))

	return deltar

def CoreShiftMeasure(coreshift12,v1,v2,kr,z,DL):
	DL = DL*10**6
	dem1 = (v2**(1./kr)-v1**(1./kr))
	csmeasure = 4.85*10**(-9)*coreshift12*DL*v1**(1./kr)*v2**(1./kr)/((1.+z)**2*dem1)

	return csmeasure

def B1nosimplification(alpha0,gammamin,gammamax,theta,psi,csmeasure,kr,delta,z):
	#1pc = 3*10**(18)cm

	betaapp = 17. #for 0836

	r1 = 3.08*10**(18) #cm
	csmeasure = csmeasure*3.08*10**(18)*(10**9)**(1./kr) #cmHz

	me = 9.11*10**(-28) #g
 	e = 4.8*10**(-10) #
	c = 3*10**10 #cm/s
	
	psirad = np.deg2rad(psi)
	thetarad = np.deg2rad(theta)
	
	cte = 2.*np.pi*me**2*c**4/e**3

	####
	cte2 = e**2/(me*c**3)
	var1_thetarad = csmeasure/(r1*np.sin(thetarad))
	paren1 = cte2*var1_thetarad**kr
	####

	var1_betaapp = csmeasure/(r1*(1+betaapp**2)**(1./2))

	####
	"""	
	def C2(alpha0):
		return 2**(-alpha0)*np.sqrt(3)/(8.*np.sqrt(np.pi)*(2-2*alpha0))*scp.gamma((6.-2*alpha0)/4.)*scp.gamma((2.-6.*alpha0)/12.)*scp.gamma((22.-6.*alpha0)/12.)*(scp.gamma((8.-2.*alpha0)/4.))**(-1) #right in christians paper

	"""

	C_alpha = 3**(1.-alpha0)/8.*np.sqrt(np.pi)*scp.gamma((7.-2*alpha0)/4.)*scp.gamma((5.-6.*alpha0)/12.)*scp.gamma((25.-6.*alpha0)/12.)*(scp.gamma((9.-2.*alpha0)/4.))**(-1) #right in christians paper

	numK = (gammamax/gammamin)**(2.*alpha0)-1.
	demK = (gammamax/gammamin)**(2.*alpha0+1.)-1.
	if alpha0 == -0.5:
		K_alphagamma = 0.1
	else:
		K_alphagamma = (2.*alpha0+1.)/(2.*alpha0)*numK/demK


	cte3 = r1*me*c**2/e**2
	var_gammamin = -2*alpha0/(gammamin**(2.*alpha0+1.))

	paren2 = np.pi*C_alpha*cte3*var_gammamin*psirad/np.sin(thetarad)*K_alphagamma*(delta/(1.+z))**(3./2.-alpha0)
	####

	#values of the exponential of the two big parenthesis
	exp1 = (5.-2.*alpha0)/(7.-2.*alpha0) #first parenthesis
	exp2 = -2./(7.-2.*alpha0) #second parenthesis

	B1 = cte*paren1**exp1*paren2**exp2

	#using betaapp
	#B12 = cte*(cte2*var1_betaapp**kr)**exp1*(np.pi*C_alpha*cte3*var_gammamin*K_alphagamma*(1.+betaapp)**(-(1.+2*alpha0)/4.)*(1.+z)**(alpha0-3./2.))**exp2

	return B1

def B1approximation(theta,psi,csmeasure,delta,z):

	#1pc = 3*10**(18)cm

	betaapp = 17. #for 0836

	r1 = 3.08*10**(18) #cm
	#csmeasure pcGHz

	me = 9.11*10**(-28) #g
 	e = 4.8*10**(-10) #
	c = 3*10**10 #cm/s

	psirad = np.deg2rad(psi)
	thetarad = np.deg2rad(theta)
	
	B1 = 0.025*(csmeasure**3*(1+z)**2/(delta**2*psirad*(np.sin(thetarad))**2))**(1./4.)

	return B1

def Bcorenoapprox(csmeasure,theta,freqcore,kr,B1):

	thetarad = np.deg2rad(theta)

	rcore = csmeasure/(np.sin(thetarad))*freqcore**(-1./kr)

	Bcore = B1*rcore**(-1)
	print 'rcore', rcore

	return Bcore, rcore

def Bcoreapprox(csmeasure,theta,freqcore,B1):

	thetarad = np.deg2rad(theta)

	rcore = csmeasure/(np.sin(thetarad))*freqcore**(-1.)

	Bcore = B1*rcore**(-1)

	print 'rcore', rcore

	return Bcore, rcore

def Particledensity(gammamin,gammamax,alpha0,B):

	me = 9.11*10**(-28) #g
 	e = 4.8*10**(-10) #
	c = 3*10**10 #cm/s

	numK = (gammamax/gammamin)**(2.*alpha0)-1.
	demK = (gammamax/gammamin)**(2.*alpha0+1.)-1.
	if alpha0 == -0.5:
		K_alphagamma = 0.1
	else:
		K_alphagamma = (2.*alpha0+1.)/(2.*alpha0)*numK/demK

	cte = (8.*np.pi*me*c**2)

	N = K_alphagamma/cte*gammamin**(-1)*B**2

	return N

class popup_beams(QWidget):
	def __init__(self,parent=None,widget=None):
     	 	QWidget.__init__(self,parent)
	   	layout = QGridLayout(self)
		#layout.addStretch(1)
		#layout.addLayout(hbox)
     	  	self.checks2 = []

           	for i in xrange(0,len(SpixWindow.diff_beams)):
          		c = QRadioButton(str('%1.3f' % (SpixWindow.diff_beams[i]))+' mas')
	  		c.setFixedSize(100,25)
          		layout.addWidget(c,0,i+1)
          		self.checks2.append(c)

	   
        	self.selectButton = QPushButton("&Select")
		self.selectButton.setAutoDefault(True)
		
		layout.addWidget(self.selectButton, 1, trunc(len(SpixWindow.diff_beams)/2.+1.))

		self.setLayout(layout)
		self.adjustSize()

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())


class BfieldcoreshiftWindow(QWidget):

    def __init__(self,*args):

        QWidget.__init__(self)

        layout = QGridLayout()

	self.xshift = []
	self.yshift = []
	self.coreshift = []
	self.freqs = []
	self.kr = 0.
	self.krerr = 0.
	self.freqv1 = 0.
	self.freqv2 = 0.
	self.freqcore = 0.
	self.csMeas = 0.
	self.zvalue = 0.
	self.DLvalue = 0.
	self.scalevalue = 0.
	self.gammaminvalue = 0.
	self.gammamaxvalue = 0.
	self.alpha0value = 0.
	self.deltavalue = 0.
	self.openingangle = 0.
	self.viewingangle = 0.
	self.B1 = 0.
	self.Bcore = 0.
	self.rcore = 0.
	self.N1 = 0.
	self.Ncore = 0.

	mathtext = u'\u03bd\u2081'
	mathtext2 = u'\u03bd\u2082'

	self.labelTEXT = QLabel()
	self.labelTEXT.setText('Select '+mathtext+' and '+mathtext2+':')
	self.labelTEXT.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

	self.checksv1v2 = []

        for i in xrange(0,len(needed_param.freq)):
            c = QCheckBox('%s GHz' % ('%1.2f' % (needed_param.freq[i])),self)
	    c.setFixedSize(100,25)
            layout.addWidget(c,2,i+1)
            self.checksv1v2.append(c)

	self.labelTEXT2 = QLabel()
	self.labelTEXT2.setText('Select the frequency for which you want to calculate the core Bfield:')
	self.labelTEXT2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

	self.checksvcore = []

        for i in xrange(0,len(needed_param.freq)):
            c = QCheckBox('%s GHz' % ('%1.2f' % (needed_param.freq[i])),self)
	    c.setFixedSize(100,25)
            layout.addWidget(c,4,i+1)
            self.checksvcore.append(c)

	self.labelempty = QLabel()
	self.labelempty.setFixedSize(25,25)
	self.labelempty2 = QLabel()

	self.approx = QRadioButton('approximation')
	self.noapprox = QRadioButton('no approximation')

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


	self.labeldelta = QLabel()
	self.labeldelta.setText("Delta : ")
	self.labeldelta.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.delta = QLineEdit()
	self.delta.setValidator(QDoubleValidator())
	self.delta.textChanged.connect(self.check_state)
	self.delta.textChanged.emit(self.delta.text())
	self.delta.setFixedSize(100,25)

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

        self.coreshiftmeasurementbutton = QPushButton("&Coreshift \n Measurement")
	#self.coreshiftmeasurementbutton.setFixedSize(100,25)
        self.coreshiftmeasurementbutton.clicked.connect(lambda: self.coremeasurement())
	self.coreshiftmeasurementbutton.setAutoDefault(True)

        self.Bfieldbutton = QPushButton("&Magnetic \n field")
	#self.Bfieldbutton.setFixedSize(100,25)
        self.Bfieldbutton.clicked.connect(lambda: self.Bfieldcalculation())
	self.Bfieldbutton.setAutoDefault(True)

        self.Nbutton = QPushButton("&Particle \n density")
	#self.Nbutton.setFixedSize(100,25)
        self.Nbutton.clicked.connect(lambda: self.Ncalculation())
	self.Nbutton.setAutoDefault(True)


	self.labelkr = QLabel()
	self.labelkr.setText(u"k\u1d63     : ")
	self.labelkr.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelkr.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
	self.krscreen = QLineEdit()
	self.krscreen.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.krscreen.setFixedSize(50,25)
	self.krerrscreen = QLabel()
	self.krerrscreen.setAlignment(Qt.AlignLeft | Qt.AlignVCenter) 
	self.krerrscreen.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.krerrscreen.setFixedSize(100,25)

	self.labelcsmeas = QLabel()
	self.labelcsmeas.setText(u"\u03a9 \u1d63 \u1d65 : ")
	self.labelcsmeas.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelcsmeas.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
	self.csmeasscreen = QLabel()
	self.csmeasscreen.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.csmeasscreen.setFixedSize(125,25)


	self.labelB1 = QLabel()
	self.labelB1.setText("B1 : ")
	self.labelB1.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelB1.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
	self.B1screen = QLabel()
	self.B1screen.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.B1screen.setFixedSize(125,25)
	self.B1unit = QLabel()
	self.B1unit.setText("G")
	self.B1unit.setAlignment(Qt.AlignLeft | Qt.AlignVCenter) 
	self.B1unit.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.B1unit.setFixedSize(25,25)

	self.labelB = QLabel()
	self.labelB.setText("Bc : ")
	self.labelB.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelB.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
	self.Bscreen = QLabel()
	self.Bscreen.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.Bscreen.setFixedSize(125,25)
	self.Bcunit = QLabel()
	self.Bcunit.setText("G")
	self.Bcunit.setAlignment(Qt.AlignLeft | Qt.AlignVCenter) 
	self.Bcunit.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.Bcunit.setFixedSize(25,25)


	self.labelrc = QLabel()
	self.labelrc.setText("rc : ")
	self.labelrc.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelrc.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
	self.rcscreen = QLabel()
	self.rcscreen.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.rcscreen.setFixedSize(125,25)
	self.rcunit = QLabel()
	self.rcunit.setText("pc")
	self.rcunit.setAlignment(Qt.AlignLeft | Qt.AlignVCenter) 
	self.rcunit.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.rcunit.setFixedSize(25,25)

	self.labelN1 = QLabel()
	self.labelN1.setText("N1 : ")
	self.labelN1.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelN1.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
	self.N1screen = QLabel()
	self.N1screen.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.N1screen.setFixedSize(125,25)

	self.labelN = QLabel()
	self.labelN.setText("Nc : ")
	self.labelN.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelN.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
	self.Nscreen = QLabel()
	self.Nscreen.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.Nscreen.setFixedSize(125,25)
	temp = searchNEDnoGUI(needed_param.source_name) 
	#omegarnu coreshift measurement

	self.DL.setText('%1.3f' % (temp[0]))
	self.z.setText('%1.3f' % (temp[1]))
	self.scale.setText('%1.3f' % (temp[2]))

	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 0)	

	layout.addWidget(self.labelTEXT, 1, 1,1,5)
	layout.addWidget(self.labelTEXT2, 3, 1,1,5)

	layout.addWidget(self.approx, 5, trunc(len(needed_param.freq)/2))
	layout.addWidget(self.noapprox, 5, trunc(len(needed_param.freq)/2)+1)
	layout.addWidget(self.labelDL, 7,1)
	layout.addWidget(self.labelz, 8,1)
	layout.addWidget(self.labelscale, 9,1)
	layout.addWidget(self.labelgammamin, 11,1)
	layout.addWidget(self.labelgammamax, 12,1)
	layout.addWidget(self.labelalpha0, 7,4)
	layout.addWidget(self.labeldelta, 8,4)
	layout.addWidget(self.labelviewing, 10,4)
	layout.addWidget(self.labelopening, 11,4)
	layout.addWidget(self.DL, 7,2)
	layout.addWidget(self.z, 8,2)
	layout.addWidget(self.scale, 9,2)
	layout.addWidget(self.gammamin, 11,2)
	layout.addWidget(self.gammamax, 12,2)
	layout.addWidget(self.alpha0, 7,5)
	layout.addWidget(self.delta, 8,5)
	layout.addWidget(self.viewing, 10,5)
	layout.addWidget(self.opening, 11,5)
	layout.addWidget(self.labelDLunit, 7,3)
	layout.addWidget(self.labelscaleunit, 9,3)

	layout.addWidget(self.labelkr, 14,1)
	layout.addWidget(self.krscreen, 14,2,1,3)
	layout.addWidget(self.krerrscreen, 14,3,1,3)
	layout.addWidget(self.labelcsmeas, 15,1)
	layout.addWidget(self.csmeasscreen, 15,2,1,3)
	layout.addWidget(self.labelB1, 16,1)
	layout.addWidget(self.B1screen, 16,2)
	layout.addWidget(self.B1unit, 16,3)
	layout.addWidget(self.labelB, 16,4)
	layout.addWidget(self.Bscreen, 16,5)
	layout.addWidget(self.Bcunit, 16,6)
	layout.addWidget(self.labelrc, 16,7)
	layout.addWidget(self.rcscreen, 16,8)
	layout.addWidget(self.rcunit, 16,9)
	layout.addWidget(self.labelN1, 17,1)
	layout.addWidget(self.N1screen, 17,2)
	layout.addWidget(self.labelN, 17,4)
	layout.addWidget(self.Nscreen, 17,5)

	#assigning buddies for elements in the layout that are tied together
	self.labelDL.setBuddy(self.DL)
	self.labelz.setBuddy(self.z)
	self.labelscale.setBuddy(self.scale)
	self.labelgammamin.setBuddy(self.gammamin)
	self.labelgammamax.setBuddy(self.gammamax)
	self.labelalpha0.setBuddy(self.alpha0)
	self.labeldelta.setBuddy(self.delta)
	self.labelviewing.setBuddy(self.viewing)
	self.labelopening.setBuddy(self.opening)
	self.labelkr.setBuddy(self.krscreen)
	self.labelcsmeas.setBuddy(self.csmeasscreen)
	self.labelB1.setBuddy(self.B1screen)
	self.labelB.setBuddy(self.Bscreen)
	self.labelrc.setBuddy(self.rcscreen)
	self.labelN1.setBuddy(self.N1screen)
	self.labelN.setBuddy(self.Nscreen)

	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 6)	

	layout.addWidget(self.coreshiftmeasurementbutton, 7,7,2,1)
	layout.addWidget(self.Bfieldbutton, 9,7,2,1)
	layout.addWidget(self.Nbutton, 11,7,2,1)

	self.coreshiftmeasurementbutton.setEnabled(False)
	self.Bfieldbutton.setEnabled(False)
	self.Nbutton.setEnabled(False)

	for i in xrange(len(self.checksv1v2)):
		self.checksv1v2[i].toggled.connect(lambda checked: self.checksState())
	for i in xrange(len(self.checksvcore)):
		self.checksvcore[i].toggled.connect(lambda checked: self.checksState())


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
	for i in xrange(0,len(self.checksv1v2)):
		if self.checksv1v2[i].isChecked():
			checked.append(i)
	for i in xrange(0,len(self.checksvcore)):
		if self.checksvcore[i].isChecked():
			checked2.append(i)

	if len(checked) == 2:
		self.coreshiftmeasurementbutton.setEnabled(True)
	else:
		self.coreshiftmeasurementbutton.setEnabled(False)

	if len(checked) == 2 and len(checked2)==1:
		self.Bfieldbutton.setEnabled(True)
		self.Nbutton.setEnabled(True)

	else:
		self.Bfieldbutton.setEnabled(False)
		self.Nbutton.setEnabled(False)


    def coremeasurement(self):


	# needed_param.path+'/Shift_parameters/shift_param'+self.freq1name+'and'+self.freq2name+'.p','wb')
	xshifttemp = []
	yshifttemp = []
	#change the int round , do  not forget, the name of the files changed after last update %'%1.2f' %(self.freq2)
	for i in xrange(0,len(needed_param.freq)-1):

		if needed_param.freq[i+1] < 0.5:
			self.freq2name = str('%1.0f' %(needed_param.freq[i+1]*1000))
			self.freq2unit = 'MHz'
		else:
			self.freq2name = str('%1.2f' %(needed_param.freq[i+1]))
			self.freq2unit = 'GHz'
		if needed_param.freq[i] < 0.5:
			self.freq1name = str('%1.0f' %(needed_param.freq[i]*1000))
			self.freq1unit = 'MHz'
		else:
			self.freq1name = str('%1.2f' %(needed_param.freq[i]))
			self.freq1unit = 'GHz'

		for filename in sorted(glob.glob('Shift_parameters/shift_param'+str(self.freq1name)+'and'+str(self.freq2name)+'.p')):   
			print filename
			res=open(filename,'rb')
			pick = pickle.load(res)
			res.close() 			
			xshifttemp.append(pick[0])
			yshifttemp.append(pick[1])
			
			#self.noise.append(pick[2])


	self.xshift = []
	self.yshift = []
	self.coreshift = [0.]
	self.coreshift2 = [0.]
	a = -1
	for i in xrange(len(needed_param.freq)-2,-1,-1):
		if i == len(needed_param.freq)-2:
			#tempx2 = 0. + xshifttemp[i]
			tempx = -0.05 + xshifttemp[i]
			tempy = 0 + yshifttemp[i]		
		else:
			tempx = self.xshift[a] + xshifttemp[i]
			tempy = self.yshift[a] + yshifttemp[i]
		self.xshift.append(tempx)
		self.yshift.append(tempy)
		self.coreshift2.append(np.sqrt(tempx**2+tempy**2))
		#self.coreshift2.append(np.sqrt(tempx2**2+tempy**2))
		a = a+1

	self.coreshift2 = np.asarray(self.coreshift2)

	#self.coreshift2[len(self.coreshift2)-1] = 0.7

	self.coreshift = self.coreshift2.copy()
	for i in xrange(len(self.coreshift)-2,0,-1):
		self.coreshift[i] = self.coreshift[i]-0.05

	res=open(needed_param.path+'/coreshiftValuesAll.p','wb')
	pickle.dump(self.coreshift,res)
	res.close()    

	errCoreshift = np.asarray([0.1,0.07,0.05,0.01])[::-1]
	errCoreshift = None


	n = len(needed_param.freq)
	freqs =  needed_param.freq[0:n]#-1]
	self.freqs = freqs[::-1]

	chi2rcore = probfit.Chi2Regression(rcoreFunct, self.freqs,self.coreshift2 , error=errCoreshift, weights=errCoreshift)

	try:

		PLFit = iminuit.Minuit(chi2rcore, A = 1., kr = 1.)

		PLFit.migrad()
		PLFit.hesse()

		print PLFit.values
		print PLFit.errors

		errCoreshift = np.asarray([0.1,0.07,0.05,0.01])[::-1]
		plt.figure(1)

		plt.errorbar(self.freqs,self.coreshift,yerr=errCoreshift,marker='.',color='r',markersize=8,linestyle='',linewidth=2)
		#plt.errorbar(self.freqs,self.coreshift2,yerr=errCoreshift,marker='.',color='g',markersize=8,linestyle='',linewidth=2)
		xaxis=np.linspace(self.freqs[0]+0.2*self.freqs[0],self.freqs[n-1]-0.2*self.freqs[n-1],1000)
		plt.plot(xaxis,rcoreFunct(xaxis,PLFit.values.get('A'), 1.),'b--',label='kr=1')
		#plt.plot(xaxis,rcoreFunct(xaxis,PLFit.values.get('A'), 0.8),'g--',label='kr=0.8')
		plt.plot(xaxis,rcoreFunct(xaxis,PLFit.values.get('A'), PLFit.values.get('kr')),'r-',label='kr= '+str(('%1.2f' % (PLFit.values.get('kr')))))
		ax = plt.gca()
		ax.minorticks_on()
		ax.tick_params('both',length=10,width=2,which='major')
		ax.tick_params('both',length=5,width=1,which='minor')
		ax.set_xticks(ticks=freqs)
		ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
		plt.legend(loc='upper right',numpoints=1,frameon=True,handletextpad=0,markerscale=1,fontsize = 12)
		plt.ylabel('offset [mas]')
		plt.xlabel(r'$\nu$ [GHz]')
		plt.ylim(-0.1,0.8)
		#plt.xlim(0,30)
		plt.show()

		print PLFit.values.get('kr')

	except RuntimeError:
		print 'Covariance is not valid. May be the last Hesse call failed?'   

	self.kr = PLFit.values.get('kr')
	self.krerr = PLFit.errors.get('kr')

	freqsv1v2 = []
	itemp = []
	for i in xrange(0,len(self.checksv1v2)):
		if self.checksv1v2[i].isChecked():
			itemp.append(i)
			freqsv1v2.append(needed_param.freq[i])

	self.freqv1 = freqsv1v2[0]
	self.freqv2 = freqsv1v2[1]

	self.DLvalue = float(self.DL.text())
	self.zvalue = float(self.z.text())
	self.scalevalue = float(self.scale.text())

	if self.approx.isChecked():
		self.csMeas = CoreShiftMeasure(self.coreshift2[::-1][itemp[0]],self.freqv1,self.freqv2,1.,self.zvalue,self.DLvalue)


	if self.noapprox.isChecked():
		self.csMeas = CoreShiftMeasure(self.coreshift2[::-1][itemp[0]],self.freqv1,self.freqv2,self.kr,self.zvalue,self.DLvalue)

	self.krscreen.setText(("%s " % ('%1.3f' % (self.kr))))
	self.krerrscreen.setText((u" \u00b1 %s" % ('%1.3f' % (self.krerr))))
	self.csmeasscreen.setText((u"%s \u00b1 %s" % ('%1.3f' % (self.csMeas), '%1.3f' % (self.krerr))))


    def Bfieldcalculation(self):
	self.kr = float(self.krscreen.text())
	self.viewingangle = float(self.viewing.text())
	self.openingangle = float(self.opening.text())
	self.deltavalue = float(self.delta.text())

	for i in xrange(0,len(self.checksvcore)):
		if self.checksvcore[i].isChecked():
			self.freqcore = needed_param.freq[i]

	if self.approx.isChecked():
		self.B1 = B1approximation(self.viewingangle,self.openingangle,self.csMeas,self.deltavalue,self.zvalue)
		self.B1screen.setText(("%s " % ('%1.5f' % (self.B1))))

		coreValues = Bcoreapprox(self.csMeas,self.viewingangle,self.freqcore,self.B1)
		self.Bcore = coreValues[0]
		self.rcore = coreValues[1]
		self.Bscreen.setText(("%s " % ('%1.5f' % (self.Bcore))))
		self.rcscreen.setText(("%s " % ('%1.2f' % (self.rcore))))

	if self.noapprox.isChecked():
		self.alpha0value = float(self.alpha0.text())
		self.gammaminvalue = float(self.gammamin.text())
		self.gammamaxvalue = float(self.gammamax.text())

		self.B1 = B1nosimplification(self.alpha0value,self.gammaminvalue,self.gammamaxvalue,self.viewingangle,self.openingangle,self.csMeas,self.kr,self.deltavalue,self.zvalue)
		self.B1screen.setText(("%s " % ('%1.5f' % (self.B1))))

		coreValues = Bcorenoapprox(self.csMeas,self.viewingangle,self.freqcore,self.kr,self.B1)
		self.Bcore = coreValues[0]
		self.rcore = coreValues[1]
		self.Bscreen.setText(("%s " % ('%1.5f' % (self.Bcore))))
		self.rcscreen.setText(("%s " % ('%1.2f' % (self.rcore))))

	print self.B1, self.Bcore

    def Ncalculation(self):
	self.alpha0value = float(self.alpha0.text())
	self.gammaminvalue = float(self.gammamin.text())
	self.gammamaxvalue = float(self.gammamax.text())

	self.N1 = Particledensity(self.gammaminvalue,self.gammamaxvalue,self.alpha0value,self.B1)
	self.Ncore = Particledensity(self.gammaminvalue,self.gammamaxvalue,self.alpha0value,self.Bcore)

	self.N1screen.setText(("%s " % ('%1.5f' % (self.N1))))
	self.Nscreen.setText(("%s " % ('%1.5f' % (self.Ncore))))
	




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

	w = BfieldcoreshiftWindow()
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


