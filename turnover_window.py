import threading, time
import warnings
import sys
import sip
import codecs
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
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

class SpectrumWindow(QWidget):

    diff_beams = []
    conv_files = []

    def __init__(self,*args):
        QWidget.__init__(self)

        self.layout = QGridLayout()
	#layout.addStretch(1)
	#layout.addLayout(hbox)

	#initializing parameters
        self.checks = []
        self.fits_freqs = []
	self.fitsFreqsFile = []

	self.files_chosen = []
	self.fits_chosen = []
	self.models_chosen = [] 
	self.shifted_files = []

	self.bmaj_files = 0.
	self.bmin_files = 0.
	self.bpa_files = 0.
	self.beam_files = 0.
	self.beam_area = 0.
	self.cells_files = 0.
	self.mapsize_file = 0.

	self.limplot_x1 = 0.
	self.limplot_x2 = 0.
	self.limplot_y1 = 0.	
	self.limplot_y2 = 0.

	self.freqs_conv = []
	self.dataFits = []
	self.noise = []
	self.datamax = []

	self.dataTurnoverS = []
	self.dataTurnoverNu = []
	self.dataTurnoveralpha0 = []
	self.dataTurnoveralphathick = []

	self.dataTurnoverSupper = []
	self.dataTurnoverNuupper = []
	self.dataTurnoveralpha0upper = []
	self.dataTurnoveralphathickupper = []

	self.ext = []

	self.rms1 = 0.
	self.rms2 = 0.
	self.ext = []

	self.position = []
	self.size = []

        for i in xrange(0,len(needed_param.freq)):
	    if needed_param.units[i] == 'MHz':
           	 c = QCheckBox('%s %s' % ('%1.0f' % (needed_param.freq[i]), needed_param.units[i]),self)
	    elif needed_param.units[i] == 'GHz':
           	 c = QCheckBox('%s %s' % ('%1.2f' % (needed_param.freq[i]), needed_param.units[i]),self)
	    c.setFixedSize(100,25)
            self.layout.addWidget(c,1,i+1)
            self.checks.append(c)

        for i in xrange(0,len(needed_param.freq)):
            c = QLabel()
	    c.setText("FITS %d: " %(1+i))
	    c.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
            c2file = QLabel()
            self.layout.addWidget(c,i+5,1)
            self.layout.addWidget(c2file,i+5,2,1,8)
            self.fits_freqs.append(c)
	    self.fitsFreqsFile.append(c2file)
	    temp = i

        #for i in xrange(0,len(freq)):
        #    c = QPushButton("Beam "+str('%1.2f' % (freq[i])))
        #    layout.addWidget(c,1,i+1)
        #    self.buttons_freqs.append(c)

	#self.lineEditBMAJ = QLineEdit(self.gridLayoutWidget)
        #self.lineEditBMAJ.setObjectName(_fromUtf8("lineEditBMAJ"))

	self.labelempty = QLabel()
	self.labelempty.setFixedSize(25,25)
	self.labelempty2 = QLabel()

	#self.withshift = QRadioButton('shifted')
	#self.withoutshift = QRadioButton('unshifted')

	self.labelOK = QLabel()
	#self.labelOK.setText("hi")

	self.labelSigmaCut = QLabel()
	self.labelSigmaCut.setText("Sigma Cut    : ")
	self.labelSigmaCut.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.SigmaCut = QLineEdit()
	self.SigmaCut.setValidator(QDoubleValidator())
	self.SigmaCut.textChanged.connect(self.check_state)
	self.SigmaCut.textChanged.emit(self.SigmaCut.text())
	self.SigmaCut.setText('%1.1f' % (3.0))
	self.SigmaCut.setFixedSize(100,25)




	for i in xrange(0,17):
		self.layout.addWidget(self.labelempty, i, 0)	

	self.layout.addWidget(self.labelOK, len(needed_param.freq)+5,1,1,5)

	self.powerLawFit = QRadioButton('Power Law')
	self.synchrotronFit = QRadioButton('Synchrotron')

	self.mathText = u'\u03B1\u1D57\u02B0\u2071\u1D9C\u1D4F = ' #\u209C\u2095\u1D62\u2096 thik as subscript

	self.group = QButtonGroup()
        self.group.addButton(self.powerLawFit)
        self.group.addButton(self.synchrotronFit)  

	self.labelMissingPickle = QLabel()
	self.layout.addWidget(self.labelMissingPickle,len(needed_param.freq)+6,1,2,5)

	self.layout.addWidget(self.labelSigmaCut, temp+8,1)
	self.layout.addWidget(self.SigmaCut, temp+8,2)

	self.labelSigmaCut.setBuddy(self.SigmaCut)

	self.synchrotronFit.toggled.connect(lambda checked: self.syncFitState())
	#self.powerLawFit.toggled.connect(lambda checked: self.PLFitState())

	self.layout.addWidget(self.powerLawFit, 2, trunc(len(needed_param.freq)/2))
	self.layout.addWidget(self.synchrotronFit, 2,trunc(len(needed_param.freq)/2)+1)

	"""layout.addWidget(self.withshift, 2, trunc(len(needed_param.freq)/2))
	layout.addWidget(self.withoutshift, 2, trunc(len(needed_param.freq)/2)+1)

	self.labelSigmaCut.setBuddy(self.SigmaCut)"""

	#for disable or enable buttons
	for i in xrange(len(self.checks)):
		self.checks[i].toggled.connect(lambda checked: self.checksState())


	#self.progressBar = QProgressBar(self)
        #self.progressBar.setRange(0,100)
        self.ShiftAllbutton = QPushButton("&Shift all")
	#self.ShiftAllbutton.setFixedSize(100,25)
        self.ShiftAllbutton.clicked.connect(lambda: self.shiftingALL())
	self.ShiftAllbutton.setAutoDefault(True)
        self.PickShiftbutton = QPushButton("&Pick files")
	#self.ShiftAllbutton.setFixedSize(100,25)
        self.PickShiftbutton.clicked.connect(lambda: self.PickShift())
	self.PickShiftbutton.setAutoDefault(True)
        self.Spectrumbutton = QPushButton("&Spectrum")
	#self.ShiftAllbutton.setFixedSize(100,25)
        self.Spectrumbutton.clicked.connect(lambda: self.calculateSpect())
	self.Spectrumbutton.setAutoDefault(True)

	#disable the buttons when opening the GUI
	self.ShiftAllbutton.setEnabled(False)
	self.PickShiftbutton.setEnabled(False)
	self.Spectrumbutton.setEnabled(False)

	if len(needed_param.freq) > 2:		
		self.layout.addWidget(self.ShiftAllbutton, 8, len(needed_param.freq)+2)
		self.layout.addWidget(self.PickShiftbutton, 9, len(needed_param.freq)+2)
		self.layout.addWidget(self.Spectrumbutton, 10, len(needed_param.freq)+2)
		for i in xrange(0,17):
			self.layout.addWidget(self.labelempty2, i, len(needed_param.freq)+2)	
	else:
		self.layout.addWidget(self.ShiftAllbutton, 8, len(needed_param.freq)+4)
		self.layout.addWidget(self.PickShiftbutton, 9, len(needed_param.freq)+4)
		self.layout.addWidget(self.Spectrumbutton, 10, len(needed_param.freq)+4)
		for i in xrange(0,17):
			self.layout.addWidget(self.labelempty2, i, len(needed_param.freq)+4)	

	#layout.addWidget(self.progressBar,17,1,1,len(freq))
	self.alphaThickLabelValue = None 
	self.alphaThickValue = None 

        self.setLayout(self.layout)

	#put the window in the center of the desktop
	qr = self.frameGeometry()
	cp = QDesktopWidget().availableGeometry().center()
	qr.moveCenter(cp)
	self.move(qr.topLeft())

	"""self.setWindowTitle("Shifting")

	"""

    """
    #function to shift an image using fft
    """

    def shift_fft(self,shiftValues,realDAT):

	#print self.realDAT
	shift_rows,shift_cols = shiftValues
	#print shift
	nr,nc = realDAT.shape
	Nr, Nc = fftfreq(nr), fftfreq(nc)
	Nc,Nr = np.meshgrid(Nc,Nr)
	fft_inputarray = np.fft.fft2(realDAT)
	fourier_shift = np.exp(1j*2*np.pi*((shift_rows*Nr)+(shift_cols*Nc)))
	output_array = np.fft.ifft2(fft_inputarray*fourier_shift)

	shiftedImage = np.real(output_array)

	return shiftedImage


    #disables buttons
    def checksState(self):
	checked = []
	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			checked.append(i)
	if len(checked) > 2:
		self.ShiftAllbutton.setEnabled(True)
		self.PickShiftbutton.setEnabled(True)
		self.Spectrumbutton.setEnabled(True)
	else:
		self.ShiftAllbutton.setEnabled(False)
		self.PickShiftbutton.setEnabled(False)
		self.Spectrumbutton.setEnabled(False)

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


    def removeButtons(self):
    	for cnt in reversed(range(self.dvbox.count())):
    	    # takeAt does both the jobs of itemAt and removeWidget
    	    # namely it removes an item and returns it
    	    widget = self.dvbox.takeAt(cnt).widget()

    	    if widget is not None: 
    	        # widget will be None if the item is a layout
    	        widget.deleteLater()

    def syncFitState(self):
	if self.synchrotronFit.isChecked() == True:
		#create the radiobuttons here, because they are deleted if radiobutton synchrotronfit not checked
		self.alphaThickLabel = QLabel()
		self.alphaThickLabel.setFixedSize(100,25)
		self.alphaThickLabel.setText(self.mathText)
		self.alphaThickLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
		#label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
		self.alphaThick = QRadioButton('2.5')
		self.alphaThickFree = QRadioButton('Free')
		self.alphaThickCustom = QRadioButton('Custom')
		self.layout.addWidget(self.alphaThickLabel, 3, trunc(len(needed_param.freq)/2)-1)
		self.layout.addWidget(self.alphaThick, 3, trunc(len(needed_param.freq)/2))
		self.layout.addWidget(self.alphaThickFree, 3,trunc(len(needed_param.freq)/2)+1)		
		self.layout.addWidget(self.alphaThickCustom, 3,trunc(len(needed_param.freq)/2)+2)

		self.alphaThickCustom.toggled.connect(lambda: self.CustomState())

	else:
		self.layout.removeWidget(self.alphaThickLabel)
		self.alphaThickLabel.deleteLater()
		self.alphaThickLabel = None
		self.layout.removeWidget(self.alphaThick)
		self.alphaThick.deleteLater()
		self.alphaThick = None
		self.layout.removeWidget(self.alphaThickFree)
		self.alphaThickFree.deleteLater()
		self.alphaThickFree = None
		self.layout.removeWidget(self.alphaThickCustom)
		self.alphaThickCustom.deleteLater()
		#sip.delete(self.alphaThickCustom)
		self.alphaThickCustom = None

		"""#create alphaThickLabelValue and alphaThickValue because i want to remove the slot 
		#if goes directly to power law after pressing custom
		self.alphaThickLabelValue = QLabel()
		self.alphaThickValue = QLineEdit()
		self.layout.addWidget(self.alphaThickLabelValue, 4, trunc(len(needed_param.freq)/2))
		self.layout.addWidget(self.alphaThickValue, 4, trunc(len(needed_param.freq)/2)+1)

		self.layout.removeWidget(self.alphaThickLabelValue)
		self.alphaThickLabelValue.deleteLater()
		self.alphaThickLabelValue = None
		self.layout.removeWidget(self.alphaThickValue)
		self.alphaThickValue.deleteLater()
		self.alphaThickValue = None"""
		if self.alphaThickLabelValue:
			self.layout.removeWidget(self.alphaThickLabelValue)
			self.alphaThickLabelValue.deleteLater()
			self.alphaThickLabelValue = None
		if self.alphaThickValue:
			self.layout.removeWidget(self.alphaThickValue)
			self.alphaThickValue.deleteLater()
			self.alphaThickValue = None

    def CustomState(self):
	if self.alphaThickCustom.isChecked() == True:
		self.alphaThickLabelValue = QLabel()
		self.alphaThickLabelValue.setFixedSize(100,25)
		self.alphaThickLabelValue.setText(self.mathText)
		self.alphaThickLabelValue.setAlignment(Qt.AlignRight | Qt.AlignVCenter) 
		self.alphaThickValue = QLineEdit()
		self.alphaThickValue.setValidator(QDoubleValidator())
		self.alphaThickValue.textChanged.connect(self.check_state)
		self.alphaThickValue.setText('%1.1f' % (1.0))
		self.alphaThickValue.setFixedSize(100,25)
		self.layout.addWidget(self.alphaThickLabelValue, 4, trunc(len(needed_param.freq)/2))
		self.layout.addWidget(self.alphaThickValue, 4, trunc(len(needed_param.freq)/2)+1)

	if self.alphaThickCustom.isChecked() == False:
		self.layout.removeWidget(self.alphaThickLabelValue)
		self.alphaThickLabelValue.deleteLater()
		self.alphaThickLabelValue = None
		self.layout.removeWidget(self.alphaThickValue)
		self.alphaThickValue.deleteLater()
		self.alphaThickValue = None




    def shiftingALL(self):

	SpectrumWindow.conv_files = []
	for filename in sorted(glob.glob(needed_param.path+'/CONV_ALL/*')):   
			SpectrumWindow.conv_files.append(filename)   

	###ORDER THEM BY FREQUENCY!
	ordered_params_shift = order_by_nu(SpectrumWindow.conv_files,SpectrumWindow.conv_files,SpectrumWindow.conv_files,False)
	SpectrumWindow.conv_files = ordered_params_shift[10]

	TempChecked = []

	for i in xrange(0,len(self.fitsFreqsFile)):
		self.fitsFreqsFile[i].setText(" ")

	self.freqs_conv = []
	self.files_chosen = []
	self.fits_chosen = []
	self.models_chosen = []
	index = []
	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			index.append(i)
			header = take_header(SpectrumWindow.conv_files[i],False)
			self.freqs_conv.append(header[5])
			self.datamax = header[10]
			TempText = SpectrumWindow.conv_files[i]
			TempChecked.append(SpectrumWindow.conv_files[i])
			self.fits_chosen.append(SpectrumWindow.conv_files[i])
			self.files_chosen.append(needed_param.files[i])
			self.models_chosen.append(needed_param.models[i])
			self.fitsFreqsFile[i].setText(TempText[len(needed_param.path):])
			
	self.fitsFreqsFile[len(TempChecked)-1].setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
	allOK = True

	for i in xrange(0,len(TempChecked)-1):
		OK = check_map_params(TempChecked[i],TempChecked[len(TempChecked)-1],False)[0]
		if OK == True:
			self.fitsFreqsFile[index[i]].setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
		else:
			self.fitsFreqsFile[index[i]].setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
			allOK = False

	if allOK == False:
		self.labelOK.setText("Convolve the files in red with the same beam,cell and mapsize than the others")
		self.labelOK.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
	elif allOK == True:
		self.labelOK.setText("Files OK")
		self.labelOK.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')

		header1 = take_header(TempChecked[0],False)

		#obtaining the beam and cell
		self.bmaj_files = header1[1]
		self.bmin_files = header1[2]
		self.bpa_files = header1[3]
		self.beam_files = header1[7]
		self.cells_files = header1[0]

		self.dataFits = []
		for i in xrange(len(self.freqs_conv)):
			map_data1 = read_map(TempChecked[i],False)
			#not in append because it should be the same for all of them
			x1 = map_data1[1]
			x2 = map_data1[2]
			y1 = map_data1[3]
			y2 = map_data1[4]
	
			self.ext=[x1,x2,y1,y2] 

			#self.noise.append(take_header(TempChecked[i],False)[9])
			self.mapsize_files = 2*map_data1[7]		
			self.dataFits.append(map_data1[0])

			if self.freqs_conv[i] < 0.5:
				self.freqname = str('%1.0f' %(self.freqs_conv[i]*1000))
				self.frequnit = 'MHz'
			else:
				self.freqname = str('%1.2f' %(self.freqs_conv[i]))
				self.frequnit = 'GHz'

	
			res=open(needed_param.path+'/rms'+self.freqname+'.p','rb')
			pick = pickle.load(res)
			res.close()  

			self.noise.append(pick) 

		xshift_l = []
		yshift_l = []


		tempPickleFiles = []

		for filename in sorted(glob.glob(needed_param.path+'/Shift_parameters/*.p')):   
			tempPickleFiles.append(filename)   

		if len(tempPickleFiles) == len(self.freqs_conv)-1:

			#self.noise = []
			for i in xrange(0,len(self.freqs_conv)-1):

				if self.freqs_conv[i+1] < 0.5:
					self.freq2name = str('%1.0f' %(self.freqs_conv[i+1]*1000))
					self.freq2unit = 'MHz'
				else:
					self.freq2name = str('%1.2f' %(self.freqs_conv[i+1]))
					self.freq2unit = 'GHz'
				if self.freqs_conv[i] < 0.5:
					self.freq1name = str('%1.0f' %(self.freqs_conv[i]*1000))
					self.freq1unit = 'MHz'
				else:
					self.freq1name = str('%1.2f' %(self.freqs_conv[i]))
					self.freq1unit = 'GHz'

				for filename in sorted(glob.glob('Shift_parameters/shift_param'+str(self.freq1name)+'and'+str(self.freq2name)+'.p')):  
					#for filename in sorted(glob.glob('Shift_parameters/shift_param'+str(int(round(self.freqs_conv[i])))+'and'+str(int(round(self.freqs_conv[i+1])))+'.p')):   
					res=open(filename,'rb')
					pick = pickle.load(res)
					res.close() 			
					xshift_l.append(pick[0])
					yshift_l.append(pick[1])
					#self.noise.append(pick[2])
	
			xshift_l.append(0.) #for taking in account the zero shift of the highest frequency
			yshift_l.append(0.)

			self.xoffset = np.zeros(np.shape(xshift_l))
			self.yoffset = np.zeros(np.shape(yshift_l))

			for i in xrange(0,len(self.xoffset)):
				self.xoffset[i] = -(-xshift_l[i]/self.cells_files)
				self.yoffset[i] = -(yshift_l[i]/self.cells_files)


			self.shifted_files = []
			for i in xrange(0,len(self.freqs_conv)):
				self.shifted_files.append('SHIFT_ALL/'+needed_param.source_name+'_'+str(round(self.freqs_conv[i]))+'convolved_shifted.fits')

			#os.system('rm difmap.log*\n')

			shiftX = 0.
			shiftY = 0.
			for i in xrange(len(self.freqs_conv)-1,-1,-1):
				shiftX = shiftX + self.xoffset[i] #xshift_l[i]
				shiftY = shiftY + self.yoffset[i]#yshift_l[i]
				print shiftX,shiftY

				shiftedfreq = self.shift_fft(np.asarray([shiftY,shiftX]),self.dataFits[i])
				os.system('rm '+self.shifted_files[i]+'\n')

				headerfits = pf.getheader(self.fits_chosen[i])

				hdu = pf.PrimaryHDU(data=shiftedfreq,header=headerfits)
				hdulist = pf.HDUList([hdu])
				hdulist.writeto(self.shifted_files[i])#,overwrite=True) #CHANGE SAVE NAME
	

				#convolve_difmap([self.files_chosen[i]],[self.models_chosen[i]],self.bmaj_files,self.bmin_files,self.bpa_files,shiftX,shiftY,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[i]])

				#os.system('rm difmap.log*\n')
	
		else:
			self.labelMissingPickle.setText("You have one or more shifts remaining. \nPlease shift the remaining frequency pairs")
			self.labelMissingPickle.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')

    def PickShift(self):

	for i in xrange(0,len(self.fitsFreqsFile)):
		self.fitsFreqsFile[i].setText(" ")

	SpectrumWindow.conv_files = []
	for filename in sorted(glob.glob(needed_param.path+'/CONV_ALL/*')):   
			SpectrumWindow.conv_files.append(filename)   

	ordered_params_shift = order_by_nu(SpectrumWindow.conv_files,SpectrumWindow.conv_files,SpectrumWindow.conv_files,False)
	SpectrumWindow.conv_files = ordered_params_shift[10]

	self.freqs_conv = []
	for i in xrange(0,len(self.checks)):
		#if self.checks[i].isChecked():
		header = take_header(SpectrumWindow.conv_files[i],False)
		self.freqs_conv.append(header[5])

	self.shifted_files = []
	for i in xrange(0,len(SpectrumWindow.conv_files)):
		self.shifted_files.append('SHIFT_ALL/'+needed_param.source_name+'_'+str(round(self.freqs_conv[i]))+'convolved_shifted.fits')

	print self.shifted_files

	allOK = True
	TempChecked = []
	self.dataFits = []
	self.freqs_conv = []
	self.noise = []
	self.datamax = []
	index = []
	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			index.append(i)
			map_data1 = read_map(self.shifted_files[i],True)
			TempText = self.shifted_files[i]
			TempChecked.append(self.shifted_files[i])
			self.fitsFreqsFile[i].setText(TempText)
			self.dataFits.append(map_data1[0])
			x1 = map_data1[1]
			x2 = map_data1[2]
			y1 = map_data1[3]
			y2 = map_data1[4]
	
			
			header = take_header(TempText,True)
			self.ext=[x1,x2,y1,y2] 
			#self.noise.append(header[9])
			self.datamax.append(header[10])
			self.mapsize_files = 2*map_data1[7]
			header = take_header(SpectrumWindow.conv_files[i],False)
			self.freqs_conv.append(header[5])

			if header[5] < 0.5:
				self.freqname = str('%1.0f' %(header[5]*1000))
				self.frequnit = 'MHz'
			else:
				self.freqname = str('%1.2f' %(header[5]))
				self.frequnit = 'GHz'

			res=open(needed_param.path+'/rms'+self.freqname+'.p','rb')
			pick = pickle.load(res)
			res.close()  

			self.noise.append(pick) 


	for i in xrange(0,len(TempChecked)):
		OK = check_map_params(TempChecked[i],TempChecked[len(TempChecked)-1],True)[0]
		if OK == True:
			self.fitsFreqsFile[index[i]].setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
		else:
			self.fitsFreqsFile[index[i]].setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
			allOK = False

	if allOK == True:
		header1 = take_header(TempChecked[0],True)

		#obtaining the beam and cell
		self.bmaj_files = header1[1]
		self.bmin_files = header1[2]
		self.bpa_files = header1[3]
		self.beam_files = header1[7]
		self.beam_area = 2*np.pi/(8*np.log(2))*self.bmaj_files*self.bmin_files  #necesary to convert Jy/beam in Jy
		self.cells_files = header1[0]


		self.dataTurnoverS = np.zeros(np.shape(self.dataFits[0]))
		self.dataTurnoverNu = np.zeros(np.shape(self.dataFits[0]))
		self.dataTurnoveralpha0 = np.zeros(np.shape(self.dataFits[0]))
		self.dataTurnoveralphathick = np.zeros(np.shape(self.dataFits[0]))

		self.dataTurnoverSupper = np.zeros(np.shape(self.dataFits[0]))
		self.dataTurnoverNuupper = np.zeros(np.shape(self.dataFits[0]))
		self.dataTurnoveralpha0upper = np.zeros(np.shape(self.dataFits[0]))
		self.dataTurnoveralphathickupper = np.zeros(np.shape(self.dataFits[0]))


		self.dataTurnoverS[:] = np.nan
		self.dataTurnoverNu[:] = np.nan
		self.dataTurnoveralpha0[:] = np.nan
		self.dataTurnoveralphathick[:] = np.nan

		self.dataTurnoverSupper[:] = np.nan
		self.dataTurnoverNuupper[:] = np.nan
		self.dataTurnoveralpha0upper[:] = np.nan
		self.dataTurnoveralphathickupper[:] = np.nan

	
	

    def calculateSpect(self):

	#get the value of alphaThick if synchrotron is selected

	guess_alphathick = 0.
	if self.synchrotronFit.isChecked():
		checked_synch = True
		if self.alphaThick.isChecked() or self.alphaThickFree.isChecked():
			guess_alphathick = 2.5
	
		if self.alphaThickCustom.isChecked():
			if self.alphaThickValue:
				guess_alphathick = float(self.alphaThickValue.text())

	else:
		checked_synch = False
		#print guess_alphathick


	final_data=np.dstack(self.dataFits)
	tmp1=ma.masked_less_equal(final_data,0)

	#print final_data.shape

	tmpres=zeros_like(tmp1)
	tmpres[:,:,:]=float('nan')

	
	#plt.figure(1)
	#plt.imshow(np.log10(tmp1), cmap='afmhot')
	#plt.xlim(860,1336)
	#plt.ylim(0,2000)
	#show()
	cutout = cuttingTURN(self.dataFits,final_data,tmp1,self.cells_files,self.dataTurnoverNu,self.dataTurnoverNuupper,checked_synch)

	cutout_data= []
	if len(cutout[0]) > 4:
		cutout_data = [cutout[0][0],cutout[0][1],cutout[0][2],cutout[0][3]]
	else:
		for i in xrange(len(cutout[0])):
			cutout_data.append(cutout[0][i])


	final_data1=np.dstack(cutout_data)
	tmp2=ma.masked_less_equal(final_data1,0)

	#plt.figure(1)
	#plt.imshow(tmp2, cmap='afmhot')
	#plt.ylim(0,final_data1.shape[0])
	#show() 
	
	sigma_cut = float(self.SigmaCut.text())
	xstart = cutout[3][0]
	xf=cutout[3][1]
	yf=cutout[3][2]
	ypix=cutout[3][3]
	print xstart,xf,ypix,yf
	"""xstart = 1049
	xf=1061
	yf=1000
	ypix=990"""
	fail=0
	j=0

	xstart = int(trunc(xstart))
	xf = int(trunc(xf)+1)
	yf = int(trunc(yf)+1)
	ypix = int(trunc(ypix))

	print xstart,xf,ypix,yf

	ydata = final_data[ypix,xstart,:]

	guess = 0



	dataNU = []
	dataS = []
	dataSpix = []
	dataSpixThick = []


	dataNUupper = []
	dataSupper = []
	dataSpixupper = []
	dataSpixThickupper = []


	xpixL = []
	ypixL = []


	xpixLupper = []
	ypixLupper = []

	if self.synchrotronFit.isChecked():
		while ypix<yf:
			xpix=xstart
			while xpix<xf:            
        			spec=final_data[ypix,xpix,:]        
        			if (spec > np.asarray(self.noise)*sigma_cut).all(): #and spec[3] >  self.noise[3]*sigma_cut:
	

        				spec=np.asarray(spec)*self.beam_area
               				#plt.figure(j)
               				#plt.plot(self.freqs_conv,spec,'ro')
               				#plt.xscale('log')
		       			#plt.yscale('log')
		       			#plt.savefig('Plot/y_'+str('%1.2f' % (ypix))+'x_'+str('%1.2f' % (xpix))+'.png')
		       			#j=j+1
		       			#plt.close('all')
		       			tmpres[ypix,xpix,0]=5
		       			print 'y = ', ypix, 'x = ', xpix
		       			print 'spec = ', spec

					plt.close('all')


					if np.max(spec) == spec[0]:
						guess_alpha0 = (spec[len(spec)-1]-spec[0])/(self.freqs_conv[len(self.freqs_conv)-1]-self.freqs_conv[0])
						guess_cte = spec[len(spec)-1]-guess_alpha0*self.freqs_conv[len(self.freqs_conv)-1]


						x=log10(np.asarray(self.freqs_conv))
						y=log10(spec)


						chi2PL = probfit.Chi2Regression(powerLawPlot, np.asarray(self.freqs_conv),spec , error=None, weights=None)

						try:
							PLValues = []

							PLFit = iminuit.Minuit(chi2PL, cte = guess_cte, alpha0 = guess_alpha0)

							PLFit.migrad()
							PLFit.hesse()
							#print synchrotronFit.values  
							PLValues.append(PLFit.values)

							constant = 5
							nuCutHigh = constant*self.freqs_conv[len(self.freqs_conv)-1]
							ScutHigh = PLValues[0].get('cte')*nuCutHigh**PLValues[0].get('alpha0')
							print 'nuCutHigh =', nuCutHigh, 'ScutHigh = ', ScutHigh

							nuCutLow = self.freqs_conv[0]/constant
							ScutLow = ScutHigh*spec[0]/spec[len(spec)-1]
							print 'nuCutLow =', nuCutLow, 'ScutLow = ', ScutLow

							specNew = np.asarray([0.]*(len(spec)+2))
							freqsNew = np.asarray([0.]*(len(self.freqs_conv)+2))
							specNew[0] = ScutLow
							freqsNew[0] = nuCutLow
					
							for i in xrange(0,len(spec)):	
								specNew[i+1] = spec[i]
								freqsNew[i+1] = self.freqs_conv[i]

							specNew[len(specNew)-1] = ScutHigh
							freqsNew[len(freqsNew)-1] = nuCutHigh

							print specNew, freqsNew

							if guess == 0:
								guess_params = guesses_turnoverPoint(freqsNew,specNew)
								#guess_alphathick = 2.5
								guess_Sm = guess_params[0]
								guess_vm = guess_params[1]
								#taum = 3./2*np.sqrt(1-guess_alpha0/guess_alphathick)-3./2
								guess_v1 = 5.615 #np.exp(np.log(guess_vm)-np.log(taum)/(guess_alpha0-guess_alphathick))
								guess_S1 = 1.543 #synchrotron(guess_v1,guess_Sm, guess_vm, guess_alpha0,guess_alphathick)
		
								parameters = [guess_S1,guess_v1,guess_alpha0,guess_alphathick]

								guess = 1

							plt.ioff()	
							chi2Sync = probfit.Chi2Regression(synchrotron, freqsNew, specNew, error=None, weights=None)

							try:
								synchrotronValues = []

								if self.alphaThickFree.isChecked():
									synchrotronFit = iminuit.Minuit(chi2Sync, S_m = guess_Sm, v_m = guess_vm , alpha0 = PLValues[0].get('alpha0'), 
										alphathick = guess_alphathick,limit_S_m = (0.,3*np.max(final_data)),limit_v_m=(0.,10**4), 
										limit_alphathick = (0.5,3.5),fix_alpha0 = True)

								if self.alphaThick.isChecked() or self.alphaThickCustom.isChecked():
									synchrotronFit = iminuit.Minuit(chi2Sync, S_m = guess_Sm, v_m = guess_vm , 
										alphathick = guess_alphathick, alpha0 = PLValues[0].get('alpha0'), 
										limit_S_m = (0.,3*np.max(final_data)),limit_v_m=(0.,10**4), fix_alphathick = True,fix_alpha0 = True)#,limit_v_m=(0.,10**4))

								#limit_S_m = (0.,3*np.max(self.dataFits))
								synchrotronFit.migrad()
								synchrotronFit.hesse()
								synchrotronValues.append(synchrotronFit.values)

								j=j+1
								print freqsNew,specNew
				       				f = plt.figure(j)
								ax = f.add_subplot(111)
								#plt.ylim(10**(-2),10**1)
								plt.xlim(freqsNew[0]-0.4*freqsNew[0],freqsNew[len(freqsNew)-1]+0.4*freqsNew[len(freqsNew)-1])
					       			plt.xscale('log')
					       			plt.yscale('log')
								plt.errorbar(freqsNew[1:(len(freqsNew)-1)],specNew[1:(len(freqsNew)-1)],yerr=0.1*specNew[1:(len(freqsNew)-1)],fmt='ro',ecolor='r', capthick=2)
								plt.errorbar([nuCutLow,nuCutHigh],[ScutLow,ScutHigh],yerr=[0.3*ScutLow,0.3*ScutHigh],uplims=True, color = 'b', marker ='o', linestyle='none')
								xaxis=np.linspace(freqsNew[0]-0.2*freqsNew[0],freqsNew[len(freqsNew)-1]+0.2*freqsNew[len(freqsNew)-1],1000)
								plt.plot(xaxis,synchrotron(xaxis,synchrotronValues[0].get('S_m'), synchrotronValues[0].get('v_m'), 
									 synchrotronValues[0].get('alpha0'), synchrotronValues[0].get('alphathick')),'g-')
								plt.ylabel(r'$S_y$ [Jy]')
								plt.xlabel(r'$\nu$ [GHz]')
					       			#plt.plot(xaxis,synchrotron(xaxis, Sm, vm, guess_alpha0,guess_alphathick),'g-')
								ax = plt.gca()
								#ax.get_xaxis().get_major_formatter().set_useOffset(False)
								ax.minorticks_on()
								ax.tick_params('both', length=10, width=2, which='major')
								ax.tick_params('both', length=5, width=1, which='minor')
								ax.set_xticks(ticks=freqsNew)
								ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())	
					       			plt.savefig('Plot_fitted_synchrotronExtrapolated/y_'+str('%1.2f' % (ypix))+'x_'+str('%1.2f' % (xpix))+'.png')

								print 'Sm = ',synchrotronValues[0].get('S_m'), 'vm = ', synchrotronValues[0].get('v_m')
								print 'alpha =', synchrotronValues[0].get('alpha0'), 'alpha_thick =', synchrotronValues[0].get('alphathick')

								plt.close('all')

								#guess_alphathick = synchrotronValues[0].get('alphathick')
								#guess_alpha0 = synchrotronValues[0].get('alpha0')
								#guess_Sm = synchrotronValues[0].get('S_m')
								#guess_vm = synchrotronValues[0].get('v_m')

								self.dataTurnoverSupper[ypix,xpix] = synchrotronValues[0].get('S_m')
								self.dataTurnoverNuupper[ypix,xpix] = synchrotronValues[0].get('v_m')
								self.dataTurnoveralpha0upper[ypix,xpix] = synchrotronValues[0].get('alpha0')
								self.dataTurnoveralphathickupper[ypix,xpix] = synchrotronValues[0].get('alphathick')


								xpixLupper.append(xpix)
								ypixLupper.append(ypix)
								dataNUupper.append(synchrotronValues[0].get('v_m'))
								dataSupper.append(synchrotronValues[0].get('S_m'))
								dataSpixupper.append(synchrotronValues[0].get('alpha0'))
								dataSpixThickupper.append(synchrotronValues[0].get('alphathick'))

								file_nameGuess = 'turnoverdataGuess.txt'
								figname = 'turnovernu.png'
								fignameS = 'turnoverS.png'

							except (RuntimeError,OverflowError,ZeroDivisionError) as error:
								print 'Covariance is not valid. May be the last Hesse call failed?'   

						except (RuntimeError,OverflowError,ZeroDivisionError) as error:
							print 'Covariance is not valid. May be the last Hesse call failed?'   

					elif np.max(spec) == spec[len(spec)-1]:
						guess_alpha0 = (spec[len(spec)-1]-spec[0])/(self.freqs_conv[len(self.freqs_conv)-1]-self.freqs_conv[0])
						guess_cte = spec[len(spec)-1]-guess_alpha0*self.freqs_conv[len(self.freqs_conv)-1]

						x=log10(np.asarray(self.freqs_conv))
						y=log10(spec)


						chi2PL = probfit.Chi2Regression(powerLaw, np.asarray(self.freqs_conv),spec , error=None, weights=None)

						try:
							PLValues = []

							PLFit = iminuit.Minuit(chi2PL, cte = guess_cte, alpha0 = guess_alpha0)

							PLFit.migrad()
							PLFit.hesse()
							#print synchrotronFit.values  
							PLValues.append(PLFit.values)

							nuCutmax = 5*self.freqs_conv[0]

						except (RuntimeError,OverflowError,ZeroDivisionError) as error:
							print 'Covariance is not valid. May be the last Hesse call failed?'   

					else:
						chi2Sync = probfit.Chi2Regression(synchrotron, np.asarray(self.freqs_conv), spec, error=None, weights=None)

						if guess == 0:
							guess_params = guesses_turnover(self.freqs_conv,spec)
							guess_alpha0 = guess_params[2]
							#guess_alphathick = 2.5
							guess_Sm = guess_params[0]
							guess_vm = guess_params[1]
							#taum = 3./2*np.sqrt(1-guess_alpha0/guess_alphathick)-3./2
							guess_v1 = 5.615 #np.exp(np.log(guess_vm)-np.log(taum)/(guess_alpha0-guess_alphathick))
							guess_S1 = 1.543 #synchrotron(guess_v1,guess_Sm, guess_vm, guess_alpha0,guess_alphathick)
	
							parameters = [guess_S1,guess_v1,guess_alpha0,guess_alphathick]

							guess = 1

						plt.ioff()	
						try:
							try:
								synchrotronValues = []

								if self.alphaThickFree.isChecked():
									synchrotronFit = iminuit.Minuit(chi2Sync, S_m = guess_Sm, v_m = guess_vm , alpha0 = guess_alpha0, alphathick = guess_alphathick, 
										limit_S_m = (0.,3*np.max(self.dataFits)),limit_v_m=(0.,np.max(np.asarray(self.freqs_conv))), limit_alphathick = (0.,3.5))
	
								if self.alphaThick.isChecked() or self.alphaThickCustom.isChecked():
									synchrotronFit = iminuit.Minuit(chi2Sync, S_m = guess_Sm, v_m = guess_vm , alpha0 = guess_alpha0, alphathick = guess_alphathick, 
										limit_S_m = (0.,3*np.max(self.dataFits)),limit_v_m=(0.,np.max(np.asarray(self.freqs_conv))), fix_alphathick = True)

								synchrotronFit.migrad()
								synchrotronFit.hesse()
								#print synchrotronFit.values  
								synchrotronValues.append(synchrotronFit.values)


								#taum = 3./2*np.sqrt(1-guess_alpha0/guess_alphathick)-3./2
								#vm = np.exp(np.log(taum)/(guess_alpha0-guess_alphathick)+np.log(guess_v1))
								#Sm = synchrotron_v1(vm,guess_S1, guess_v1, guess_alpha0,guess_alphathick)
								#parameters = [guess_Sm,guess_vm,guess_alpha0,guess_alphathick]

								j=j+1
				       				f = plt.figure(j)
								ax = f.add_subplot(111)
					       			plt.xscale('log')
					       			plt.yscale('log')
								#plt.ylim(10**(-2),10**1)
								plt.xlim(self.freqs_conv[0]-0.4*self.freqs_conv[0],self.freqs_conv[len(self.freqs_conv)-1]+0.4*self.freqs_conv[len(self.freqs_conv)-1])
								plt.errorbar(self.freqs_conv,spec,yerr=0.1*np.asarray(self.datamax),fmt='ro',ecolor='r', capthick=2)
								xaxis=np.linspace(self.freqs_conv[0]-0.2*self.freqs_conv[0],self.freqs_conv[len(self.freqs_conv)-1]+0.2*self.freqs_conv[len(self.freqs_conv)-1],1000)
								plt.plot(xaxis,synchrotron(xaxis,synchrotronValues[0].get('S_m'), synchrotronValues[0].get('v_m'), 
									 synchrotronValues[0].get('alpha0'), synchrotronValues[0].get('alphathick')),'g-')
								plt.ylabel(r'$S_y$ [Jy]')
								plt.xlabel(r'$\nu$ [GHz]')
					       			#plt.plot(xaxis,synchrotron(xaxis, Sm, vm, guess_alpha0,guess_alphathick),'g-')
								ax = plt.gca()
								#ax.get_xaxis().get_major_formatter().set_useOffset(False)
								ax.minorticks_on()
								ax.tick_params('both', length=10, width=2, which='major')
								ax.tick_params('both', length=5, width=1, which='minor')
								ax.set_xticks(ticks=self.freqs_conv)
								ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())	
					       			plt.savefig('Plot_fitted_synchrotron/y_'+str('%1.2f' % (ypix))+'x_'+str('%1.2f' % (xpix))+'.png')

								print 'Sm = ',synchrotronValues[0].get('S_m'), 'vm = ', synchrotronValues[0].get('v_m')
								print 'alpha =', synchrotronValues[0].get('alpha0'), 'alpha_thick =', synchrotronValues[0].get('alphathick')

								plt.close('all')

								guess_alphathick = synchrotronValues[0].get('alphathick')
								guess_alpha0 = synchrotronValues[0].get('alpha0')
								guess_Sm = synchrotronValues[0].get('S_m')
								guess_vm = synchrotronValues[0].get('v_m')

								self.dataTurnoverS[ypix,xpix] = synchrotronValues[0].get('S_m')
								self.dataTurnoverNu[ypix,xpix] = synchrotronValues[0].get('v_m')
								self.dataTurnoveralpha0[ypix,xpix] = synchrotronValues[0].get('alpha0')
								self.dataTurnoveralphathick[ypix,xpix] = synchrotronValues[0].get('alphathick')

								xpixL.append(xpix)
								ypixL.append(ypix)
								dataNU.append(synchrotronValues[0].get('v_m'))
								dataS.append(synchrotronValues[0].get('S_m'))
								dataSpix.append(synchrotronValues[0].get('alpha0'))
								dataSpixThick.append(synchrotronValues[0].get('alphathick'))
								file_name = 'turnoverdata.txt'
								figname = 'turnovernu.png'
								fignameS = 'turnoverS.png'

							except (RuntimeError,OverflowError,ZeroDivisionError) as error:
								print 'Covariance is not valid. May be the last Hesse call failed?'   

						except RuntimeWarning:
							print 'Covariance is not valid. May be the last Hesse call failed?'  

		    		#else:
		        	#	tmpres[ypix,xpix,0]=100
		    		#print fail
		    		xpix=xpix+1
			ypix=ypix+1


		final_txt1 = [xpixL, ypixL,dataNU,dataS,dataSpix,dataSpixThick]

		final_txt2 = [xpixLupper, ypixLupper,dataNUupper,dataSupper,dataSpixupper,dataSpixThickupper]

		if len(xpixL) > 0:
			header = np.array([['#x y nu S alpha alphathick']])
			final_txt=np.swapaxes(final_txt1, 0,1)
			saver(file_name, header, final_txt, FORMAT='%1.4f')

		if len(xpixLupper) > 0:
			header = np.array([['#x y nu S alpha alphathick']])
			final_txt=np.swapaxes(final_txt2, 0,1)
			saver(file_nameGuess, header, final_txt, FORMAT='%1.4f')

		first_contour = self.dataFits[len(self.dataFits)-1].std()
		levels = first_contour*np.array([2.,4.,8.,16.,32.,64.,128.,256.,512.,1024.,2048.])

		plt.figure(1)
		ax = plt.gca()
		cset = plt.contour(self.dataFits[len(self.dataFits)-1], levels, inline=1,
			          colors=['grey'],extent=self.ext, aspect=1.0)
		p1 = ax.imshow(self.dataTurnoverS, origin='bottom',extent=self.ext)#, vmin=-2.5, vmax=1.7)

		p2 = ax.imshow(self.dataTurnoverSupper, origin='bottom',extent=self.ext,cmap='spring')#, vmin=-2.5, vmax=1.7)
	
		plt.axis('scaled')
		plt.xlabel('Right Ascension [mas]')
		plt.ylabel('Relative Declination [mas]')

		divider = make_axes_locatable(ax)
		cax1 = divider.append_axes("right", size="4%", pad="0.5%")

		cax2 = divider.append_axes("right", size="4%", pad="15%")

		cb1 = plt.colorbar(p1,cax=cax1,cmap='jet')

		cb2 = plt.colorbar(p2,cax=cax2,cmap='spring')


		cb1.set_label(r'$S_y$ [Jy]')
		cb2.set_label(r'$S_y$ upper [Jy]')

	    	#plt.title('Spectral index between %s GHz and %s GHz \n %s' % ('%1.1f' % (self.freq1), '%1.1f' % (self.freq2),needed_param.source_name ))
		#plt.xlim(self.ext[0], self.ext[1])
		#plt.ylim(self.ext[2],self.ext[3])



		plt.figure(2)
		ax = plt.gca()
		cset = plt.contour(self.dataFits[len(self.dataFits)-1], levels, inline=1,
			          colors=['grey'],extent=self.ext, aspect=1.0)
		p1 = ax.imshow(self.dataTurnoverNu, origin='bottom',extent=self.ext)#, vmin=-2.5, vmax=1.7)
		p2 = ax.imshow(self.dataTurnoverNuupper, origin='bottom',extent=self.ext,cmap='spring')#, vmin=-2.5, vmax=1.7)

		plt.axis('scaled')
		plt.xlabel('Right Ascension [mas]')
		plt.ylabel('Relative Declination [mas]')

		divider = make_axes_locatable(ax)
		cax1 = divider.append_axes("right", size="4%", pad="0.5%")

		cax2 = divider.append_axes("right", size="4%", pad="15%")

		cb1 = plt.colorbar(p1,cax=cax1,cmap='jet')

		cb2 = plt.colorbar(p2,cax=cax2,cmap='spring')

		cb1.set_label(r'$\nu$ [GHz]')
		cb2.set_label(r'$\nu$ upper [GHz]')



		limits = Annotate()
		plt.show()

		[self.limplot_x1,self.limplot_x2,self.limplot_y1,self.limplot_y2] = limits()

		plt.figure(1)
		ax = plt.gca()
		cset = plt.contour(self.dataFits[len(self.dataFits)-1], levels, inline=1,
			          colors=['grey'],extent=self.ext, aspect=1.0)
		p1 = ax.imshow(self.dataTurnoverS, origin='bottom',extent=self.ext)#, vmin=-2.5, vmax=1.7)
		p2 = ax.imshow(self.dataTurnoverSupper, origin='bottom',extent=self.ext,cmap='spring')#, vmin=-2.5, vmax=1.7)

		plt.axis('scaled')
		plt.xlabel('Right Ascension [mas]')
		plt.ylabel('Relative Declination [mas]')
		plt.xlim(self.limplot_x1,self.limplot_x2)
		plt.ylim(self.limplot_y2,self.limplot_y1)

		divider = make_axes_locatable(ax)
		cax1 = divider.append_axes("right", size="4%", pad="0.5%")

		cax2 = divider.append_axes("right", size="4%", pad="15%")

		cb1 = plt.colorbar(p1,cax=cax1,cmap='jet')

		cb2 = plt.colorbar(p2,cax=cax2,cmap='spring')

		cb1.set_label(r'$S_y$ [Jy]')
		cb2.set_label(r'$S_y$ upper [Jy]')


		plt.savefig(fignameS)

		plt.close('all')


		plt.figure(2)
		ax = plt.gca()
		cset = plt.contour(self.dataFits[len(self.dataFits)-1], levels, inline=1,
			          colors=['grey'],extent=self.ext, aspect=1.0)
		p1 = ax.imshow(self.dataTurnoverNu, origin='bottom',extent=self.ext)#, vmin=-2.5, vmax=1.7)
		p2 = ax.imshow(self.dataTurnoverNuupper, origin='bottom',extent=self.ext,cmap='spring')#, vmin=-2.5, vmax=1.7)

		plt.axis('scaled')
		plt.xlabel('Right Ascension [mas]')
		plt.ylabel('Relative Declination [mas]')
		plt.xlim(self.limplot_x1,self.limplot_x2)
		plt.ylim(self.limplot_y2,self.limplot_y1)

		divider = make_axes_locatable(ax)
		cax1 = divider.append_axes("right", size="4%", pad="0.5%")

		cax2 = divider.append_axes("right", size="4%", pad="15%")

		cb1 = plt.colorbar(p1,cax=cax1,cmap='jet')

		cb2 = plt.colorbar(p2,cax=cax2,cmap='spring')

		cb1.set_label(r'$\nu$ [GHz]')
		cb2.set_label(r'$\nu$ upper [GHz]')


		plt.savefig(figname)

		plt.close('all')

		final = [self.dataTurnoverNu,self.dataTurnoverS,self.dataTurnoveralpha0,self.dataTurnoveralphathick]
		final2 = [self.dataTurnoverNuupper,self.dataTurnoverSupper,self.dataTurnoveralpha0upper,self.dataTurnoveralphathickupper]
		res=open(needed_param.path+'/turnoverdata.p','wb')
		pickle.dump(final,res)
		res.close()     
		res=open(needed_param.path+'/UpperLimturnoverdata.p','wb')
		pickle.dump(final2,res)
		res.close()     

	if self.powerLawFit.isChecked():

		while ypix<yf:
			xpix=xstart
			while xpix<xf:            
        			spec=final_data[ypix,xpix,:]        
        			if (spec > np.asarray(self.noise)*sigma_cut).all(): 

					spec = np.asarray(spec)*self.beam_area
	
					guess_alpha0 = (spec[len(spec)-1]-spec[0])/(self.freqs_conv[len(self.freqs_conv)-1]-self.freqs_conv[0])
					guess_cte = spec[len(spec)-1]-guess_alpha0*self.freqs_conv[len(self.freqs_conv)-1]


					x=log10(np.asarray(self.freqs_conv))
					y=log10(spec)


						

					plt.ioff()	

		       			tmpres[ypix,xpix,0]=5
		       			print 'y = ', ypix, 'x = ', xpix
		       			print 'spec = ', spec

					plt.close('all')

					x=log10(np.asarray(self.freqs_conv))
					y=log10(spec)


					chi2PL = probfit.Chi2Regression(powerLawPlot, np.asarray(self.freqs_conv),spec , error=None, weights=None)

					try:
						PLValues = []

						PLFit = iminuit.Minuit(chi2PL, cte = guess_cte, alpha0 = guess_alpha0)

						PLFit.migrad()
						PLFit.hesse()
						#print synchrotronFit.values  
						PLValues.append(PLFit.values)

						#taum = 3./2*np.sqrt(1-guess_alpha0/guess_alphathick)-3./2
						#vm = np.exp(np.log(taum)/(guess_alpha0-guess_alphathick)+np.log(guess_v1))
						#Sm = synchrotron_v1(vm,guess_S1, guess_v1, guess_alpha0,guess_alphathick)

						plt.ioff()
						j=j+1
			       			f = plt.figure(j)
						ax = f.add_subplot(111)
					       	plt.xscale('log')
					  	plt.yscale('log')
						#plt.ylim(10**(-2),10**1)
						#plt.xlim(10**(0),300)
						plt.errorbar(self.freqs_conv,spec,yerr=0.1*np.asarray(self.datamax),fmt='ro',ecolor='r', capthick=2)
			       			#plt.plot(self.freqs_conv,spec,'ro')
						xaxis=np.linspace(self.freqs_conv[0]-0.2*self.freqs_conv[0],self.freqs_conv[len(self.freqs_conv)-1]+0.2*self.freqs_conv[len(self.freqs_conv)-1],1000)
						yaxis = powerLawPlot(xaxis,PLValues[0].get('cte'),PLValues[0].get('alpha0'))
						plt.plot(xaxis,powerLawPlot(xaxis,PLValues[0].get('cte'),PLValues[0].get('alpha0')),'g-')
						plt.ylabel(r'$S_y$ [Jy]')
						plt.xlabel(r'$\nu$ [GHz]')
			       			#plt.plot(xaxis,synchrotron(xaxis, Sm, vm, guess_alpha0,guess_alphathick),'g-')
			       			#plt.xscale('log',nonposx='clip')
			       			#plt.yscale('log',nonposy='clip')
						ax = plt.gca()
						#ax.get_xaxis().get_major_formatter().set_useOffset(False)
						ax.minorticks_on()
						ax.tick_params('both', length=10, width=2, which='major')
						ax.tick_params('both', length=5, width=1, which='minor')
						ax.set_xticks(ticks=self.freqs_conv)
						ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())	
			       			plt.savefig('Plot_fitted_PL/y_'+str('%1.2f' % (ypix))+'x_'+str('%1.2f' % (xpix))+'.png')

						print 'cte = ',PLValues[0].get('cte')
						print 'alpha =', PLValues[0].get('alpha0')

						plt.close('all')

						guess_alpha0 = PLValues[0].get('alpha0')
						guess_cte = PLValues[0].get('cte')

						print PLValues

					except (RuntimeError,OverflowError,ZeroDivisionError) as error:
						print 'Covariance is not valid. May be the last Hesse call failed?'   
					#PLFit = iminuit(chi2_PL,cte,alpha0))
		    		else:
		        		tmpres[ypix,xpix,0]=100
		    		print fail
		    		xpix=xpix+1
			ypix=ypix+1


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
	
	freqOrig = ordered_params[0]
	files = ordered_params[8]
	models = ordered_params[9]
	fits = ordered_params[10]

	freq = freqOrig.copy()

	units = []

	for i in xrange(0,len(freq)):
		if freq[i] < 0.5:
			freq[i] = freq[i]*1000
			units.append('MHz')
		else:
			units.append('GHz')
	
	units = np.asarray(units)

	
	#source name
	header = take_header(fits[0],False)
	source_name = header[8]

def main():
	app = QApplication(sys.argv)

	w = SpectrumWindow()
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


