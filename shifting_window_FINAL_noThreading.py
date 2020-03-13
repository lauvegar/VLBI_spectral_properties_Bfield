import threading, time
import sys
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import matplotlib.pyplot as plt
from math import *
from scipy.fftpack import fftfreq
from functools import *
import numpy as np
import astropy.io.fits as pf
import pickle
from functions_conv import order_by_nu, read_conv_params
from functions_align import find_same_beam,beam_array,check_map_params, cuttingMAP,cross_correlation_shifts_FITS, checking_shift, search_rms
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, Annotate
import os,glob
import subprocess as sub
from astropy.nddata import Cutout2D
from correlate2d import *
#from fast_ftts import *


#lineedit = QLineEdit(self)
#validator = QDoubleValidator()
#lineedit.setValidator(QDoubleValidator())
#lineedit.textChanged.connect(self.check_state)
#lineedit.textChanged.emit(lineedit.text())

class popup_beams(QWidget):
	def __init__(self,parent=None,widget=None):
     	 	QWidget.__init__(self,parent)
	   	layout = QGridLayout(self)
		#layout.addStretch(1)
		#layout.addLayout(hbox)
     	  	self.checks2 = []

           	for i in xrange(0,len(ShiftWindow.diff_beams)):
          		c = QRadioButton(str('%1.3f' % (ShiftWindow.diff_beams[i]))+' mas')
	  		c.setFixedSize(100,25)
          		layout.addWidget(c,0,i+1)
          		self.checks2.append(c)

	   
        	self.selectButton = QPushButton("&Select")
		self.selectButton.setAutoDefault(True)
		
		layout.addWidget(self.selectButton, 1, trunc(len(ShiftWindow.diff_beams)/2.+1.))

		self.setLayout(layout)
		self.adjustSize()

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

class popup_shift(QWidget):
	def __init__(self):
     	 	QWidget.__init__(self)
	   	layout = QGridLayout()
		#layout.addStretch(1)
		#layout.addLayout(hbox)
     	  	
		self.labelempty = QLabel()
		self.labelTEXT = QLabel()
		self.labelTEXT2 = QLabel()
		self.labelNewSHIFT = QLabel()
		self.labelQuestion = QLabel()
          	
		self.labelTEXT.setText("There is still a residual shift in one or both direction")
		self.labelTEXT.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
		self.labelTEXT2.setText("The values of this shift are:")
		self.labelTEXT2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.labelNewSHIFT.setText(("Shift RA: %s mas,  Shift DEC: %s mas" % ('%1.4f' % (ShiftWindow.xshift_new), '%1.4f' % (ShiftWindow.yshift_new))))
		self.labelNewSHIFT.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
		self.labelQuestion.setText("Do you want to shift the image again?")
		self.labelQuestion.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')

		layout.addWidget(self.labelempty,3,0)
		layout.addWidget(self.labelTEXT,0,1,1,5)
		layout.addWidget(self.labelTEXT2,1,1,1,5)
		layout.addWidget(self.labelNewSHIFT,2,1,1,5)
		layout.addWidget(self.labelQuestion,3,1,1,5)

	   
        	self.shiftButton = QPushButton("&Shift")
        	self.noshiftButton = QPushButton("&No Shift")
		self.shiftButton.setAutoDefault(True)
		self.noshiftButton.setAutoDefault(True)
		
		layout.addWidget(self.shiftButton, 4,1,1,2)
		layout.addWidget(self.noshiftButton, 4, 3,1,2)

		self.setLayout(layout)

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

class popup_manualshift(QWidget):
	def __init__(self):
     	 	QWidget.__init__(self)
	   	layout = QGridLayout()
		#layout.addStretch(1)
		#layout.addLayout(hbox)
     	  	
		self.labelempty = QLabel()
		self.labelmSHIFTra = QLabel()
		self.labelmSHIFTra.setFixedSize(100,25)
		self.labelmSHIFTraMAS = QLineEdit()
		self.labelmSHIFTraMAS.setValidator(QDoubleValidator())
		self.labelmSHIFTraMAS.textChanged.connect(self.check_stateQLINE)
		self.labelmSHIFTraMAS.textChanged.emit(self.labelmSHIFTraMAS.text())
		self.labelmSHIFTraMAS.setFixedSize(100,25)
		self.labelmSHIFTra2 = QLabel()
		self.labelmSHIFTra2.setFixedSize(100,25)
		self.labelmSHIFTdec = QLabel()
		self.labelmSHIFTdec.setFixedSize(100,25)
		self.labelmSHIFTdecMAS = QLineEdit()
		self.labelmSHIFTdecMAS.setValidator(QDoubleValidator())
		self.labelmSHIFTdecMAS.textChanged.connect(self.check_stateQLINE)
		self.labelmSHIFTdecMAS.textChanged.emit(self.labelmSHIFTdecMAS.text())
		self.labelmSHIFTdecMAS.setFixedSize(100,25)
		self.labelmSHIFTdec2 = QLabel()
		self.labelmSHIFTdec2.setFixedSize(100,25)
          	
		self.labelmSHIFTra.setText("SHIFT RA : ")
		self.labelmSHIFTra.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.labelmSHIFTdec.setText("SHIFT DEC: ")
		self.labelmSHIFTdec.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.labelmSHIFTra2.setText(" mas ")
		self.labelmSHIFTra2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.labelmSHIFTdec2.setText(" mas ")
		self.labelmSHIFTdec2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

		layout.addWidget(self.labelempty,3,0)
		layout.addWidget(self.labelmSHIFTra,1,1)
		layout.addWidget(self.labelmSHIFTraMAS,1,2)
		layout.addWidget(self.labelmSHIFTra2,1,3)
		layout.addWidget(self.labelmSHIFTdec,2,1)
		layout.addWidget(self.labelmSHIFTdecMAS,2,2)
		layout.addWidget(self.labelmSHIFTdec2,2,3)
	   
        	self.ManualShiftButton = QPushButton("&Shift")
		self.ManualShiftButton.setAutoDefault(True)
		
		layout.addWidget(self.ManualShiftButton, 4,2)

		self.setLayout(layout)

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

	def check_stateQLINE(self,*args,**kwargs):
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

class ShiftWindow(QWidget):

    # declare class parameters that can be call in another classes
    diff_beams = []
    xshift = 0.
    yshift = 0.
    xshift_new = 0.
    yshift_new = 0.
    offset = []
    offset_new = []
    ext = []
    rms = []
    freq1Manual = 0.
    freq2Manual = 0.

    def __init__(self,*args):
        QWidget.__init__(self)

        layout = QGridLayout()
	#layout.addStretch(1)
	#layout.addLayout(hbox)

	#the number of frequencies will appear dynamically after reading the input fits files
	#for that, a list is defined which will contain the elements that will show this frequencies
	#these elements are, a row with checkboxes and a row with buttons (one per frequency)
        self.checks = []
        self.buttons_freqs = []
	self.name_button_freq = []

	#declare list of the convolved  fits files
	self.fits_conv = []

	#declare parameters of the convolved files that might change if more than one beam is present
	self.freq_conv = []
	self.beam_conv = []

	#declare the two fits files, its rms and frequency
	self.fits1 = 'l'
	self.fits2 = 'r'
	self.freq1 = 0.
	self.freq2 = 0.
	self.freq1name = ''
	self.freq2name = ''
	self.freq1unit = ''
	self.freq2unit = ''
	self.rms1 = 0.
	self.rms2 = 0.

	#declare some list to do the ordering later
	self.files_chosen = []
	self.models_chosen = [] 
	self.shifted_files = []

	#declare parameters of the convolved files once the proper one is selected
	self.bmaj_files = 0.
	self.bmin_files = 0.
	self.bpa_files = 0.
	self.beam_files = 0.
	self.cells_files = 0.
	self.mapsize_file = 0.

	#declare list containing the information of the cut of the fits file
	self.position = []
	self.size = []
	self.position_feature = []
	self.size_feature = []

	#given the number of frequencies read in the class needed_files
	#checkboxes are created for each frequency and added to the correspondent list
        for i in xrange(0,len(needed_param.freq)):
	    if needed_param.units[i] == 'MHz':
           	 c = QCheckBox('%s %s' % ('%1.0f' % (needed_param.freq[i]), needed_param.units[i]),self)
	    elif needed_param.units[i] == 'GHz':
           	 c = QCheckBox('%s %s' % ('%1.2f' % (needed_param.freq[i]), needed_param.units[i]),self)
	    c.setFixedSize(100,25)
            layout.addWidget(c,1,i+1)
            self.checks.append(c)

	#empty labels for aesthetic reasons
	self.labelempty = QLabel()
	self.labelempty.setFixedSize(25,25)
	self.labelempty2 = QLabel()

	#declare the label for the FITS1/2
	self.labelFITS1 = QLabel()
	self.labelFITS2 = QLabel()
	#fix a size for the fits1/2 label
	self.labelFITS1.setFixedSize(100,25)
	self.labelFITS2.setFixedSize(100,25)
	#declare the label for the name of the file in FITS1/2
	self.labelFITS1file = QLabel()
	self.labelFITS2file = QLabel()
	#self.labelFITS1file.setFixedSize(100,25)
	#self.labelFITS2file.setFixedSize(100,25)

	#declare label that indicates if the files are fine
	self.labelOK = QLabel()

	#writing the text of the label of FITS1/2
	self.labelFITS1.setText("FITS1 : ")
	self.labelFITS2.setText("FITS2 : ")
	#writing the label of FITS1/2 in bold
	self.labelFITS1.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelFITS2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

	#declare the label for:
	#the line where it says which frequencies you are shifting
	#for the label of the shift in RA, for the value in mas and the value in pixels
	self.freqsSHIFTING = QLabel()
	#self.freqsSHIFTING.setFixedSize(100,25)
	self.labelSHIFTra = QLabel()
	self.labelSHIFTraMAS = QLabel()
	self.labelSHIFTraPIXEL = QLabel()
	#for the label of the shift in DEC, for the value in mas and the value in pixels
	self.labelSHIFTdec = QLabel()
	self.labelSHIFTdecMAS = QLabel()
	self.labelSHIFTdecPIXEL = QLabel()

	#fixing the size of all the RA shift labels
	self.labelSHIFTra.setFixedSize(100,25)
	self.labelSHIFTraMAS.setFixedSize(100,25)
	self.labelSHIFTraPIXEL.setFixedSize(100,25)
	#fixing the size of all the DEC shift labels
	self.labelSHIFTdec.setFixedSize(100,25)
	self.labelSHIFTdecMAS.setFixedSize(100,25)
	self.labelSHIFTdecPIXEL.setFixedSize(100,25)

	#writing the text of the label of RA/DEC (no the shift values)
	self.labelSHIFTra.setText("SHIFT RA : ")
	self.labelSHIFTdec.setText("SHIFT DEC : ")
	#writing the label of RA/DEC (no the shift values) in bold
	self.labelSHIFTra.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelSHIFTdec.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 0)	

	#assigning the previous elements a position in the grid layout
	layout.addWidget(self.labelFITS1, 3,1)
	layout.addWidget(self.labelFITS2, 4,1)
	layout.addWidget(self.labelFITS1file, 3,2, 1, 5)
	layout.addWidget(self.labelFITS2file, 4,2, 1, 5)
	layout.addWidget(self.labelOK, 6,1,1,4)
	layout.addWidget(self.freqsSHIFTING, 8,1,1,5)
	layout.addWidget(self.labelSHIFTra, 9,1)
	layout.addWidget(self.labelSHIFTraMAS, 9,2)
	layout.addWidget(self.labelSHIFTraPIXEL, 9,3)
	layout.addWidget(self.labelSHIFTdec, 10,1)
	layout.addWidget(self.labelSHIFTdecMAS, 10,2)
	layout.addWidget(self.labelSHIFTdecPIXEL, 10,3)

	#assigning buddies for elements in the layout that are tied together
	self.labelFITS1.setBuddy(self.labelFITS1file)
	self.labelFITS2.setBuddy(self.labelFITS2file)
	self.labelSHIFTra.setBuddy(self.labelSHIFTraMAS)
	self.labelSHIFTdec.setBuddy(self.labelSHIFTdecMAS)
	self.labelSHIFTraMAS.setBuddy(self.labelSHIFTraPIXEL)
	self.labelSHIFTraMAS.setBuddy(self.labelSHIFTdecPIXEL)

	#able or disable buttons depending if the checkboxes conditions are fulfilled
	for i in xrange(len(self.checks)):
		self.checks[i].toggled.connect(lambda checked: self.checksState())

	#declaring the buttons, putting them a fix size, activating the button while enter pressed when the user is on it, 
	#and assing the functions for them to work, e.g, findingBEAM,shifting,ManualShifting respectively
        self.findingBEAMbutton = QPushButton("&Check")
	self.findingBEAMbutton.setFixedSize(100,25)
	self.findingBEAMbutton.setAutoDefault(True) #button clicked by pressing enter
        self.findingBEAMbutton.clicked.connect(lambda: self.findingBEAM(self.checks,needed_param.freqOrig))
        self.shiftingButton = QPushButton("&Shifting")
	self.shiftingButton.setFixedSize(100,25)
	self.shiftingButton.setAutoDefault(True)
        self.shiftingButton.clicked.connect(lambda: self.shifting(self.checks,self.fits1,self.fits2,needed_param.files,needed_param.models))
        self.ManualShiftingButton = QPushButton("&Manual \n Shifting")
        self.ManualShiftingButton.clicked.connect(lambda: self.ManualShifting())
	self.ManualShiftingButton.setAutoDefault(True)
	#shiftingButton.setFixedSize(100,25)

	#disable the buttons when opening the GUI
	self.findingBEAMbutton.setEnabled(False)
	self.shiftingButton.setEnabled(False)
	self.ManualShiftingButton.setEnabled(False)

	#assigning the previous buttons a position in the grid layout, depending the number of existing frequencies 
	if len(needed_param.freq) > 2:		
		layout.addWidget(self.findingBEAMbutton, 5, len(needed_param.freq)+1)
		layout.addWidget(self.shiftingButton, 6, len(needed_param.freq)+1)
		layout.addWidget(self.ManualShiftingButton, 7, len(needed_param.freq)+1,2,1)
		for i in xrange(0,17):
			layout.addWidget(self.labelempty2, i, len(needed_param.freq)+2)	
	else:
		layout.addWidget(self.findingBEAMbutton, 5, len(needed_param.freq)+4)
		layout.addWidget(self.shiftingButton, 6, len(needed_param.freq)+4)
		layout.addWidget(self.ManualShiftingButton, 7, len(needed_param.freq)+4,2,1)
		for i in xrange(0,17):
			layout.addWidget(self.labelempty2, i, len(needed_param.freq)+4)	

	#setting the layout
        self.setLayout(layout)

	#put the window in the center of the desktop
	qr = self.frameGeometry()
	cp = QDesktopWidget().availableGeometry().center()
	qr.moveCenter(cp)
	self.move(qr.topLeft())

	self.setWindowTitle("Shifting")

   	#for i in xrange(0,len(freq)):
	#	self.buttons_freqs[i].clicked.connect(partial(self.BEAMparam, bmaj,bmin,bpa,beam,i,self.circ)) 
		#Python only introduces new bindings in namespace through assignment and through parameter lists of functions. i is therefore not actually defined in the namespace of the lambda, but in the namespace of __init__(). The name lookup for i in the lambda consequently ends up in the namespace of __init__(), where i is eventually bound to 9. This is called "closure". You are creating closures. Closures really capture a variable, not the value of a variable. At the end of __init__, i is the last element of range(0, 10), i.e. 9. All the lambdas you created in this scope refer to this i and only when they are invoked, they get the value of i at the time they are at invoked (however, seperate invocations of __init__ create lambdas referring to seperate variables!).

    ###############################################################################################
    ###############################################################################################
    ######################                     FUNCTIONS                     ######################    
    ###############################################################################################
    ###############################################################################################
    """
    #function to shift an image using fft
    """

    def shift_fft(self,shiftValues):

	#print self.realDAT
	shift_rows,shift_cols = shiftValues
	#print shift
	nr,nc = self.realDAT.shape
	Nr, Nc = fftfreq(nr), fftfreq(nc)
	Nc,Nr = np.meshgrid(Nc,Nr)
	fft_inputarray = np.fft.fft2(self.realDAT)
	fourier_shift = np.exp(1j*2*np.pi*((shift_rows*Nr)+(shift_cols*Nc)))
	output_array = np.fft.ifft2(fft_inputarray*fourier_shift)

	shiftedImage = np.real(output_array)

	return shiftedImage


    """
    #function to check if the arguments given in the text boxes are fine
    #if the value is not the kind of parameter wanted, 
    #for example, if the box requires a double and you give an integrer, it get red
    #if the input value is still the kind of parameter wanted, but outside a range, it gets yellow
    #if the input value is fine, it gets green
    """
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

    """
    #function to disable buttons. 
    #In this case, the buttons will get enable if two frequencies are checked 
    #and disabled if less or more than two frequencies are check (or none is checked)
    """
    def checksState(self):
	checked = []
	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			checked.append(i)
	if len(checked) == 2:
		self.ManualShiftingButton.setEnabled(True)
		self.shiftingButton.setEnabled(True)
		self.findingBEAMbutton.setEnabled(True)
	else:
		self.ManualShiftingButton.setEnabled(False)
		self.shiftingButton.setEnabled(False)
		self.findingBEAMbutton.setEnabled(False)

    #select the mapsize (mapsize of the low frequency image) and cellsize (half of the cellsize of the high frequency image)

    """
    #function called by findingBEAM
    #this function calls a pop up that allows the user to select one beam of the valid ones for the frequency pair
    #select the files containing that beam for the frequency pair
    parameters used:
    RadioButtons =  self.wi.checks2 ---> radiobuttons containing the different valid beams for the frequency pair
    diff_beams = ShiftWindow.diff_beams ---> list with the values of the valid beams (if more than one is common for both frequencies)
						created and given from findingBEAM
    beam1 = self.beam1 ---> array containing the beams of all the files of the first frequency
						created and given from findingBEAM
    beam2 = self.beam2 ---> arrays containing the beams of all the files of the second frequency
						created and given from findingBEAM
    freq1_index = self.freq1_index ---> arrays containing the indexes of all the files of the first frequency
						created and given from findingBEAM
    freq2_index = self.freq2_index --->  arrays containing the indexes of all the files of the second frequency
						created and given from findingBEAM
    fits_conv = self.fits_conv ---> list with all the fits files convolved,
						created and given from findingBEAM

    self.selected_beam ---> float, value of the beam selected by the user

    self.fits1,self.fits2 --->  str, final two fits files selected
    self.labelFITS1file, self.labelFITS2file, self.labelOK --->  to update the interface 

    """
    def getBEAM(self,RadioButtons,diff_beams,beam1,beam2,freq1_index,freq2_index,fits_conv):

			#checks with radiobutton is selected
			for i in xrange(len(RadioButtons)):
				if RadioButtons[i].isChecked():
					self.selected_beam = ShiftWindow.diff_beams[i]

			#finds the index of selected beam in the beam arrays for both frequencies
			index_b1 = np.where(beam1 == self.selected_beam)
			index_b2 = np.where(beam2 == self.selected_beam) #output a list containing tuple
			#converting the list into array and the tuples in the list into an array
			index_beam1 = np.asarray([x for xs in index_b1 for x in xs]) 
			index_beam2 = np.asarray([x for xs in index_b2 for x in xs]) 
			#search the index of the fits file for both frequencies that has the selected beam 
			index_fits1 = freq1_index[0][index_beam1]
			index_fits2 = freq2_index[0][index_beam2]

			#gets the name of the fits file
			self.fits1 = fits_conv[index_fits1[0]]
			self.fits2 = fits_conv[index_fits2[0]]
	
			#terminal updated with the name of the files
			self.labelFITS1file.setText(self.fits1[len(needed_param.path):])
			self.labelFITS2file.setText(self.fits2[len(needed_param.path):])

			#parameters of the files checked (beam, cell, mapsize)	
			see_if_OK = check_map_params(str(self.fits1),str(self.fits2),False)
			OK = see_if_OK[0]

			#if the parameters are fine, names appear in green, OK label says OK
			if OK == True:
				self.labelFITS1file.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
				self.labelFITS2file.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
				self.labelOK.setText("Files OK")
				self.labelOK.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
			#if the parameters are not fine, names appear in red, OK label says reconvolve
			else:
				self.labelFITS1file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
				self.labelFITS2file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
				self.labelOK.setText("Convolve the files with the same beam,cell and mapsize")
				self.labelOK.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')

			#close popup window
			self.wi.close()


    """
    #function to find the fits files of the selected frequencies that have a valid beam
    #if there is only one valid beam (same beam for both of them) for the pair of frequencies, they get selected
    #if there is no valid beam, tells the user to convolve
    #if there is more than one common beam, a pop up appears for the user to choose the beam
    parameters used:
    checkBOXes = self.checks ---> list with the check boxes
    freq = needed_param.freq ---> list with all the existing frequencies 
    self.fits_conv ---> list with all the fits files convolved, redefined at the beginning (set to []), to allow the user to use the function again
    
    self.freq_conv ---> list with the frequencies of all fits files convolved, redefined at the beginning (by means of another function)
    self.beam_conv ---> list with the beams of all fits files convolved, redefined at the beginning (by means of another function)
    self.freq1_index, self.freq2_index  --->  arrays containing the indexes of all the files with the two indicated frequencies 
						, redefined at the beginning (by means of another function)
    self.beam1,self.beam2 --->  arrays containing the beams of all the files with the two indicated frequencies 
				, redefined at the beginning (by means of another function)
    self.index_beam12  ---> array containing the beams indexes (beams that have the same value for both of them) for both frequencies
				index corresponding to  self.freq1_index and self.freq2_index respectively
				, redefined at the beginning (by means of another function)
    self.fits1,self.fits2 --->  str, final two fits files selected, gets the new value each time it is equal to something

    ShiftWindow.diff_beams ---> list with the values of the valid beams (if more than one is common for both frequencies)
				redefined in the function
    
    self.labelFITS1file, self.labelFITS2file, self.labelOK --->  to update the interface 
    
    """
    def findingBEAM(self,checkBOXes,freq): #self.checks,freq ###stop_event,checkBOXes,freq 
	#self.myLongTask.start()

	#filling the array with all the existing convolved fits files
	#it is redefined to allow the user to repeat the process and not have stored the values of the first run as well
	self.fits_conv = []

	for filename in sorted(glob.glob(needed_param.path+'/CONV/*.fits*')):   
		self.fits_conv.append(filename)   
	files_temp = self.fits_conv
	models_temp = self.fits_conv


	#orders the convolved fits file by frequency (lower to higher)
	ordered_params_conv = order_by_nu(files_temp,files_temp,self.fits_conv,False)
	self.freq_conv = ordered_params_conv[0]
	self.beam_conv = ordered_params_conv[5]
	self.fits_conv = ordered_params_conv[10]

	#creates an array with beams and indexes of those beams for the two selected frequencies
	#finds the common beams for both frequencies
	#stores the index of the common beams for each subarray in an array, e.g. if two valid beams are found, indexbeam12 = [[0,2],[1,3]]
	#[0,2] ---> indexes in the subarray containing the beams for both frequencies for the first valid beam (in the example)
	#[1,3] ---> indexes in the subarray for the second valid beam (in the example)
	finding_beam = beam_array(checkBOXes,freq,self.freq_conv,self.beam_conv)
	self.freq1_index,self.freq2_index,self.index_beam12 = finding_beam[0],finding_beam[1],finding_beam[2]
	self.beam1,self.beam2 = finding_beam[3],finding_beam[4]

	#len(self.index_beam12)== 1 ----> only one common beam found
	#the fits files corresponding to that beam are stored, the parameters of the file checked
	#and the terminal updated
	if len(self.index_beam12)== 1:
		#the fits files corresponding to that beam are stored
		track_freq1 = self.index_beam12[0][0]
		index_freq1 = self.freq1_index[0][track_freq1]
		self.fits1 = self.fits_conv[index_freq1]

		track_freq2 = self.index_beam12[0][1]-len(self.freq1_index[0])
		index_freq2 = self.freq2_index[0][track_freq2]
		self.fits2 = self.fits_conv[index_freq2]

		#terminal updated with the name of the files
		self.labelFITS1file.setText(self.fits1[len(needed_param.path):])
		self.labelFITS2file.setText(self.fits2[len(needed_param.path):])
		
		#parameters of the files checked (beam, cell, mapsize)	
		see_if_OK = check_map_params(str(self.fits1),str(self.fits2),False)
		OK = see_if_OK[0]

		#if the parameters are fine, names appear in green, OK label says OK
		if OK == True:
			self.labelFITS1file.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
			self.labelFITS2file.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
			self.labelOK.setText("Files OK")
			self.labelOK.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
		#if the parameters are not fine, names appear in red, OK label says reconvolve
		else:
			self.labelFITS1file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
			self.labelFITS2file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
			self.labelOK.setText("Convolve the files with the same beam,cell and mapsize")
			self.labelOK.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')

	#len(self.index_beam12) > 1 ----> more than one common beam found
	#the user selects the desired beam using a popup window
	#it calls the popup class popup_beams() and the function of this class getBEAM controlling the button in the popup
	elif len(self.index_beam12) > 1:
		ShiftWindow.diff_beams = []

		for i in xrange(0,len(self.index_beam12)):
			ShiftWindow.diff_beams.append(self.beam1[self.index_beam12[i][0]])
	
		self.wi = popup_beams()
		self.wi.show()
		self.wi.selectButton.clicked.connect(lambda: self.getBEAM(self.wi.checks2, ShiftWindow.diff_beams,self.beam1,self.beam2,self.freq1_index,self.freq2_index,self.fits_conv))

	#len(self.index_beam12) == 0 ----> no valid files beam related found
	#two X appear in the fits1/2 interface widget and tells the user to convolve the files in OK label
	elif len(self.index_beam12) == 0:
		self.labelOK.setText("No match found. Convolve the frequencies with the same beam")
		self.labelOK.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
		self.labelFITS1file.setText("X")
		self.labelFITS2file.setText("X")
		self.labelFITS1file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
		self.labelFITS2file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')

	self.update()

	#return popup_beams

    """
    #puts the values of the interface in blank, for using when threading if used
    """
    def update(self):
	self.freqsSHIFTING.setText("")
	self.labelSHIFTraMAS.setText("")
	self.labelSHIFTraPIXEL.setText("")
	self.labelSHIFTdecMAS.setText("")
	self.labelSHIFTdecPIXEL.setText("")	

    """
    #
    #
    parameters used:
    checkBOXes = self.checks ---> list with the check boxes
    fits1 = self.fits1 ---> chosen fits file for freq1
    fits2 = self.fits2 ---> chosen fits file for freq2
    files = needed_param.files ---> list with all the uvf files (= total number of frequencies)
    models = needed_param.models ---> list with all the models (= total number of frequencies)

    self.bmaj_files ---> float, common bmaj for both frequencies
    self.bmin_files ---> float, common bmin for both frequencies
    self.bpa_files ---> float, common bpa for both frequencies
    self.beam_files ---> float, common circular beam for both frequencies
    self.cells_files ---> float, common cell for both frequencies
    self.mapsize_files ---> float, common mapsize for both frequencies

    self.freq1 ---> float, frequency of fits1 (lower frequency)
    self.freq2 ---> float, frequency of fits2 (higher frequency)

    ShiftWindow.ext ---> list, BLC and TRC coordinates of the user selected region for future display (first selection)
				values redefined in the function, first as the BLC and TRC of the whole map
				later values redefined by another function as the BLC and TRC of the selected region
    self.position ---> tupple, center position (x and y) of the region selected by the user (region for future display)
				values redefined inside this function using another function
    self.size ---> tupple, size in x and y of the region for future display selected by the user
				values redefined inside this function using another function
    self.position_feature ---> tupple, center position (x and y) of the region selected by the user (region of the optically thin feature)
				values redefined inside this function using another function 
    self.size_feature ---> tupple, size in x and y of the region of the optically thin feature selected by the user
				values redefined inside this function using another function		
    self.shifted_files ---> list, names (created in the function) of the shifted fits files (for the two frequencies selected)
				redefined in the function (set to []), to allow the user to use the function again
    self.models_chosen ---> list, names of the mod files of the two frequencies selected
				redefined in the function (set to []), to allow the user to use the function again
    self.files_chosen ---> list, names of the uvf files of the two frequencies selected
				redefined in the function (set to []), to allow the user to use the function again

    ShiftWindow.offset ---> list, offset value in RA and DEC in pixel
    ShiftWindow.xshift ---> float, offset value in RA in mas
    ShiftWindow.yshift ---> float, offset value in DEC in mas

    self.rms1 ---> float, rms of the first file
    self.rms2 ---> float, rms of the second file

    self.freqsSHIFTING --->  to update the interface 
    self.labelSHIFTraMAS, self.labelSHIFTraPIXEL --->  to update the interface 
    self.labelSHIFTdecMAS, self.labelSHIFTdecPIXEL --->  to update the interface 
    """
    def shifting(self,checkBOXes,fits1,fits2,files,models):

	#reads header
	header1 = take_header(fits1,False)
	header2 = take_header(fits2,False)
	map_data1 = read_map(fits1,False)		
	self.realDAT = map_data1[0]
	map_data2 = read_map(fits2,False)		
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

	if self.freq2 < 0.5:
		self.freq2name = str('%1.0f' %(self.freq2*1000))
		self.freq2unit = 'MHz'
	else:
		self.freq2name = str('%1.2f' %(self.freq2))
		self.freq2unit = 'GHz'
	if self.freq1 < 0.5:
		self.freq1name = str('%1.0f' %(self.freq1*1000))
		self.freq1unit = 'MHz'
	else:
		self.freq1name = str('%1.2f' %(self.freq1))
		self.freq1unit = 'GHz'


	#obtaining map centers in pixels from the header
	cent_mapx = map_data1[5]
	cent_mapy = map_data1[6]

	#obtaining the four corners of the maps in mas from the header
	x1 = map_data1[1]
	x2 = map_data1[2]
	y1 = map_data1[3]
	y2 = map_data1[4]

	ShiftWindow.ext=[x1,x2,y1,y2] 

	#cut maps v1 < v2, map region to study
	cut = cuttingMAP(self.realDAT,self.realDAT2,cent_mapx,cent_mapy,self.cells_files,self.freq1,self.freq2,self.freq1name,self.freq2name,self.freq1unit,self.freq2unit,iteration=0)

	#data of the region cutted
	cutout_v1 = cut[0]
	cutout_v2 = cut[1]
	 
	#position, size and BLC and TRC of the region cutted   
	self.position = cut[2]
	self.size = cut[3]
	ShiftWindow.ext = cut[4]


	#cut maps v1 < v2, feature
	cut = cuttingMAP(cutout_v1.data,cutout_v2.data,cent_mapx,cent_mapy,self.cells_files,self.freq1,self.freq2,self.freq1name,self.freq2name,self.freq1unit,self.freq2unit,iteration=1)

	#data of the optically thin feature cutted
	cutout_v1feature = cut[0]
	cutout_v2feature = cut[1]
    
	#position, size of the optically thin feature cutted 
	self.position_feature = cut[2]
	self.size_feature = cut[3]

	#plot optically thin feature for both frequencies
	plt.figure(1)
	plt.subplot(121)
	plt.imshow(cutout_v2feature.data, origin='bottom')
	#plt.axis('scaled')
	plt.xlabel('Right Ascension [pixels]')
	plt.ylabel('Relative Declination [pixels]')

	plt.title(self.freq2name+self.freq2unit)

	plt.subplot(122)
	plt.imshow(cutout_v1feature.data, origin='bottom')
	#plt.axis('scaled')
	plt.xlabel('Right Ascension [pixels]')
	#plt.ylabel('Relative Declination [pixels]')
	plt.title(self.freq1name+self.freq1unit)

	plt.show()

	#2d cross-correlation
	image1 = cutout_v1feature.data #data of the feature in map1 (freq1)
	image2 = cutout_v2feature.data #data of the feature in map2 (freq2)

	#fast spix calculation without shifting
	a = np.log10(image1 /image2)/np.log10(self.freq1/self.freq2)

	if np.any(a > 0.)==True:
		print 'Some parts of the component selected are not optically thin, you may not want that'


	#calculating the 2d crosscorrelation and setting the offset values
	ShiftWindow.offset = cross_correlation_shifts_FITS(image1, image2, sigma_cut=0.)

	ShiftWindow.xshift=-ShiftWindow.offset[0]*self.cells_files #mas
	ShiftWindow.yshift=ShiftWindow.offset[1]*self.cells_files #mas


	#creating the names of the shifted fits files
	self.shifted_files=[]
	self.shifted_files.append(needed_param.path+'/SHIFT/'+needed_param.source_name+'-'+self.freq1name+self.freq1unit+'convolved_with_beam'+str('%1.2f' % (self.beam_files))+'shifted.fits')
	self.shifted_files.append(needed_param.path+'/SHIFT/'+needed_param.source_name+'-'+self.freq2name+self.freq1unit+'convolved_with_beam'+str('%1.2f' % (self.beam_files))+'shifted.fits')
	
	#getting the mod and uvf files corresponding to the selected frequencies
	self.models_chosen = []
	self.files_chosen = []
	for i in xrange(len(checkBOXes)):
		if checkBOXes[i].isChecked():
			self.models_chosen.append(models[i])
			self.files_chosen.append(files[i])


	#updating the terminal with the shift values
	self.freqsSHIFTING.setText(("Image2 ( %s GHz)  - image1 ( %s GHz)" % ('%1.3f' % (self.freq2), '%1.3f' % (self.freq1))))
	self.freqsSHIFTING.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTraMAS.setText(("%s mas" % ('%1.4f' % (ShiftWindow.xshift))))
	self.labelSHIFTraMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTraPIXEL.setText(("%s pix" % ('%1.2f' % (-ShiftWindow.offset[0]))))
	self.labelSHIFTraPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecMAS.setText(("%s mas" % ('%1.4f' % (ShiftWindow.yshift))))
	self.labelSHIFTdecMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecPIXEL.setText(("%s pix" % ('%1.2f' % (-ShiftWindow.offset[1]))))
	self.labelSHIFTdecPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')

	#shift the files using difmap
	#os.system('rm difmap.log*\n')

	shiftedfreq1 = self.shift_fft(np.asarray([-ShiftWindow.offset[1],-ShiftWindow.offset[0]]))
	os.system('rm '+self.shifted_files[0]+'\n')

	headerfits1 = pf.getheader(fits1)

	hdu = pf.PrimaryHDU(data=shiftedfreq1,header=headerfits1)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto(self.shifted_files[0])#,overwrite=True) #CHANGE SAVE NAME

	res=open(needed_param.path+'/rms'+self.freq1name+'.p','rb')
	pick = pickle.load(res)
	res.close() 

	self.rms1 = pick	

	#convolve_difmap([self.files_chosen[0]],[self.models_chosen[0]],self.bmaj_files,self.bmin_files,self.bpa_files,ShiftWindow.xshift,ShiftWindow.yshift,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[0]])

	#rms1 (lower freq)
	#self.rms1 = search_rms()[0]
	#os.system('rm difmap.log*\n')

	headerfits2 = pf.getheader(fits2)

	#shiftedfreq2 = self.shift_fft(self.realDAT2,np.asarray([0.,0.]))
	os.system('rm '+self.shifted_files[1]+'\n')
	hdu = pf.PrimaryHDU(data=self.realDAT2,header=headerfits2)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto(self.shifted_files[1])#,overwrite=True) #CHANGE SAVE NAME

	res=open(needed_param.path+'/rms'+self.freq2name+'.p','rb')
	pick = pickle.load(res)
	res.close() 

	self.rms2 = pick

	#convolve_difmap([self.files_chosen[1]],[self.models_chosen[1]],self.bmaj_files,self.bmin_files,self.bpa_files,0,0,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[1]])

	#rms2 (higher freq)
	#self.rms2 = search_rms()[0]
	#os.system('rm difmap.log*\n')


	#check if the shift was done properly
	self.checkingShift()

    """
    #
    #
    parameters used:
    ShiftWindow.offset_new ---> list, new offset value in RA and DEC in pixel
				to be added to ShiftWindow.offset if the user decides so (clicking shift button)
				if the user does not decide to shift again, this value will not be added
    ShiftWindow.xshift_new --->  float, new offset value in RA in mas
				to be added to ShiftWindow.xshift if the user decides so (clicking shift button)
    ShiftWindow.yshift_new --->  float, new offset value in DEC in mas
				to be added to ShiftWindow.yshift if the user decides so (clicking shift button)
    self.shifted_files ---> list coming from shifting function
    self.position ---> tupple coming from shifting function
    self.size ---> tupple coming from shifting function
    self.position_feature ---> tupple coming from shifting function
    self.size_feature ---> tupple coming from shifting function
    popup class is popup_shift
    """
    def checkingShift(self):
	#checking the shift using the already shifted files and the values of the region (position,size) containing the optically thin component
	ShiftWindow.offset_new = checking_shift(self.shifted_files,self.position,self.size,self.position_feature,self.size_feature,True)
	
	#new shift values	
	ShiftWindow.xshift_new=-ShiftWindow.offset_new[0]*self.cells_files #mas
	ShiftWindow.yshift_new=ShiftWindow.offset_new[1]*self.cells_files #mas
	
	#if the residual shift is large, ask the user if a new shift is needed
	#in that case, shows the residual shift in a popup and asks the user what he/she wants to do
	if np.abs(ShiftWindow.xshift_new) > 0.005*np.abs(ShiftWindow.xshift) or np.abs(ShiftWindow.yshift_new) > 0.005*np.abs(ShiftWindow.yshift):
		self.wi = popup_shift()
		self.wi.show()
		self.wi.shiftButton.clicked.connect(lambda: self.shiftAgain())
		self.wi.noshiftButton.clicked.connect(lambda: self.noShiftAgain())

	#if the residual shift is too small, the values are directly stored in a file
	else:
		self.write_result_file()

    """
    #
    #
    parameters used:
    ShiftWindow.offset ---> list, offset value in RA and DEC in pixel
    ShiftWindow.offset_new ---> list, new offset value in RA and DEC in pixel
				to be added to ShiftWindow.offset
    ShiftWindow.xshift ---> float, offset value in RA in mas
    ShiftWindow.yshift ---> float, offset value in DEC in mas
    ShiftWindow.xshift_new --->  float, new offset value in RA in mas
				to be added to ShiftWindow.xshift
    ShiftWindow.yshift_new --->  float, new offset value in DEC in mas
				to be added to ShiftWindow.yshift
    """
    def shiftAgain(self):
	#calculate the new total shift
	final_shiftx = ShiftWindow.xshift+ShiftWindow.xshift_new
	final_shifty = ShiftWindow.yshift+ShiftWindow.yshift_new

	#set the default shift as the old shift
	ShiftWindow.xshift = final_shiftx
	ShiftWindow.yshift = final_shifty

	#updating interface
	self.labelSHIFTraMAS.setText(("%s mas" % ('%1.4f' % (ShiftWindow.xshift))))
	self.labelSHIFTraMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTraPIXEL.setText(("%s pix" % ('%1.2f' % (-ShiftWindow.offset[0]-ShiftWindow.offset_new[0]))))
	self.labelSHIFTraPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecMAS.setText(("%s mas" % ('%1.4f' % (ShiftWindow.yshift))))
	self.labelSHIFTdecMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecPIXEL.setText(("%s pix" % ('%1.2f' % (-ShiftWindow.offset[1]-ShiftWindow.offset_new[1]))))
	self.labelSHIFTdecPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')

	#shifting using python

	final_offsetx = -ShiftWindow.offset[0]-ShiftWindow.offset_new[0]
	final_offsety = -ShiftWindow.offset[1]-ShiftWindow.offset_new[1]

	shiftedfreq1 = self.shift_fft(np.asarray([final_offsety,final_offsetx]))
	os.system('rm '+self.shifted_files[0]+'\n')

	headerfits1 = pf.getheader(self.fits1)

	hdu = pf.PrimaryHDU(data=shiftedfreq1,header=headerfits1)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto(self.shifted_files[0])#,overwrite=True) #CHANGE SAVE NAME

	#convolve_difmap([self.files_chosen[0]],[self.models_chosen[0]],self.bmaj_files,self.bmin_files,self.bpa_files,final_shiftx,final_shifty,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[0]])

	#closes popup window
	self.wi.close()
	#checks the shift again
	self.checkingShift()

    """
    #if the user does not want to re-shift
    #this function just saves the important parameters in a text and pickle file
    #using the function write_result_file and closes the popup
    """
    def noShiftAgain(self):

	self.write_result_file()
	self.wi.close()

    """
    #function to write the important parameters of the shifting in a pickle and text file
    parameters used:
    ShiftWindow.xshift ---> float, shift in RA in mas
    ShiftWindow.yshift ---> float, shift in DEC in mas
    self.rms1 ---> float, rms of the first file
    self.rms2 ---> float, rms of the second file
    ShiftWindow.ext ---> list coming from shifting or manualshifting function, 
			BLC and TRC coordinates of the user selected region for future display (first selection)
			if manual shift, BLC and TRC coordinates of the whole map  
    """
    def write_result_file(self):

	#rms1 (lower freq)
	#self.rms1 = search_rms()[0]
	#os.system('rm difmap.log*\n')
	#saves the pickle file
	final = [ShiftWindow.xshift, ShiftWindow.yshift, self.rms1, self.rms2 ,ShiftWindow.ext]
	res=open(needed_param.path+'/Shift_parameters/shift_param'+self.freq1name+'and'+self.freq2name+'.p','wb')
	pickle.dump(final,res)
	res.close()     

	#saves the ascii file
	final_txt = [[ShiftWindow.xshift, ShiftWindow.yshift,self.rms1, self.rms2,ShiftWindow.ext[0],ShiftWindow.ext[1],ShiftWindow.ext[2],ShiftWindow.ext[3]]]
	header = np.array([['#RAshift DECshift rms(low freq image) rms(high freq image) coordinates of the map window -----> lower freq shifted']])
	saver(needed_param.path+'/Shift_parameters/shift_param'+self.freq1name+'and'+self.freq2name+'.txt', header, final_txt, FORMAT='%1.6f')

    """
    #function that shows a popup for the user to select the values in mas of the desired shift
    popup class is popup_manualshift
    """
    def ManualShifting(self):
	self.wi = popup_manualshift()
	self.wi.show()
	self.wi.labelmSHIFTraMAS.editingFinished.connect(self.getManualShiftX)
	self.wi.labelmSHIFTdecMAS.editingFinished.connect(self.getManualShiftY) 

    """
    #function to get the value introduced for the user in the RA shift EDITLINE in the manualshift popup
    it saves the value the user writes in --->  ShiftWindow.xshift
    """
    def getManualShiftX(self):
	ShiftWindow.xshift = float(self.sender().text()) #gives the value to the parameter after the edition of Qlineedit is done

    """
    #function to get the value introduced for the user in the DEC shift EDITLINE in the manualshift popup
    it saves the value the user writes in --->  ShiftWindow.yshift
    it also connects the button of the popup, does not matter if connecting here or in shiftX
    """
    def getManualShiftY(self):
	ShiftWindow.yshift = float(self.sender().text() )
	self.wi.ManualShiftButton.clicked.connect(lambda: self.DoManualShift(self.checks,self.fits1,self.fits2,needed_param.files,needed_param.models))

    """
    #
    #
    parameters used:
    checkBOXes = self.checks ---> list with the check boxes
    fits1 = self.fits1 ---> chosen fits file for freq1
    fits2 = self.fits2 ---> chosen fits file for freq2
    files = needed_param.files ---> list with all the uvf files (= total number of frequencies)
    models = needed_param.models ---> list with all the models (= total number of frequencies)

    self.bmaj_files ---> float, common bmaj for both frequencies
    self.bmin_files ---> float, common bmin for both frequencies
    self.bpa_files ---> float, common bpa for both frequencies
    self.beam_files ---> float, common circular beam for both frequencies
    self.cells_files ---> float, common cell for both frequencies
    self.mapsize_files ---> float, common mapsize for both frequencies

    self.freq1 ---> float, frequency of fits1 (lower frequency)
    self.freq2 ---> float, frequency of fits2 (higher frequency)

    ShiftWindow.ext ---> list, BLC and TRC coordinates of the user selected region for future display (first selection)
				values redefined in the function as the BLC and TRC of the whole map

    self.shifted_files ---> list, names (created in the function) of the shifted fits files (for the two frequencies selected)
				redefined in the function (set to []), to allow the user to use the function again
    self.models_chosen ---> list, names of the mod files of the two frequencies selected
				redefined in the function (set to []), to allow the user to use the function again
    self.files_chosen ---> list, names of the uvf files of the two frequencies selected
				redefined in the function (set to []), to allow the user to use the function again

    self.rms1 ---> float, rms of the first file
    self.rms2 ---> float, rms of the second file
    """
    def DoManualShift(self,checkBOXes,fits1,fits2,files,models):
	#self.wi.labelmSHIFTraMAS.returnPressed.connect(self.getManual)
	#reading the header
	header1 = take_header(fits1,False)
	header2 = take_header(fits2,False)
	map_data1 = read_map(fits1,False)		
	self.realDAT = map_data1[0]
	map_data2 = read_map(fits2,False)		
	self.realDAT2 = map_data2[0]		

	#obtaining the beam and cell from the header
	self.bmaj_files = header1[1]
	self.bmin_files = header1[2]
	self.bpa_files = header1[3]
	self.beam_files = header1[7]
	self.cells_files = header1[0]

	#obtaining frequencies from the header
	self.freq1 = header1[5]
	self.freq2 = header2[5]

	if self.freq2 < 0.5:
		self.freq2name = str('%1.0f' %(self.freq2*1000))
		self.freq2unit = 'MHz'
	else:
		self.freq2name = str('%1.2f' %(self.freq2))
		self.freq2unit = 'GHz'
	if self.freq1 < 0.5:
		self.freq1name = str('%1.0f' %(self.freq1*1000))
		self.freq1unit = 'MHz'
	else:
		self.freq1name = str('%1.2f' %(self.freq1))
		self.freq1unit = 'GHz'

	#obtaining mapsize from the header
	self.mapsize_files = 2*map_data1[7]

	#window limits (BLC, TRC)  from the header
	x1 = map_data1[1]
	x2 = map_data1[2]
	y1 = map_data1[3]
	y2 = map_data1[4]

	ShiftWindow.ext=[x1,x2,y1,y2] 

	#creating the names of the shifted fits files
	self.shifted_files=[]
	self.shifted_files.append(needed_param.path+'/SHIFT/'+needed_param.source_name+'-'+self.freq2name+self.freq2unit+'convolved_with_beam'+str('%1.2f' % (self.beam_files))+'shifted.fits')
	self.shifted_files.append(needed_param.path+'/SHIFT/'+needed_param.source_name+'-'+self.freq1name+self.freq1unit+'convolved_with_beam'+str('%1.2f' % (self.beam_files))+'shifted.fits')
	
	#getting the mod and uvf files corresponding to the selected frequencies
	self.models_chosen = []
	self.files_chosen = []
	for i in xrange(len(checkBOXes)):
		if checkBOXes[i].isChecked():
			self.models_chosen.append(models[i])
			self.files_chosen.append(files[i])


	

	#ShiftWindow.xshift_new=-ShiftWindow.offset_new[0]*self.cells_files #mas
	#ShiftWindow.yshift_new=ShiftWindow.offset_new[1]*self.cells_files #mas

	ShiftWindow.offset = np.asarray([0.,0.])

	ShiftWindow.offset[0] = (-ShiftWindow.xshift/self.cells_files)
	ShiftWindow.offset[1] = (ShiftWindow.yshift/self.cells_files)

	shiftedfreq1 = self.shift_fft(np.asarray([-ShiftWindow.offset[1],-ShiftWindow.offset[0]]))
	os.system('rm '+self.shifted_files[0]+'\n')

	headerfits1 = pf.getheader(fits1)

	hdu = pf.PrimaryHDU(data=shiftedfreq1,header=headerfits1)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto(self.shifted_files[0])#,overwrite=True) #CHANGE SAVE NAME

	res=open(needed_param.path+'/rms'+self.freq1name+'.p','rb')
	pick = pickle.load(res)
	res.close() 

	self.rms1 = pick	


	headerfits2 = pf.getheader(fits2)

	#shiftedfreq2 = self.shift_fft(self.realDAT2,np.asarray([0.,0.]))
	os.system('rm '+self.shifted_files[1]+'\n')
	hdu = pf.PrimaryHDU(data=self.realDAT2,header=headerfits2)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto(self.shifted_files[1])#,overwrite=True) #CHANGE SAVE NAME

	res=open(needed_param.path+'/rms'+self.freq2name+'.p','rb')
	pick = pickle.load(res)
	res.close() 

	self.rms2 = pick


	#shifting using difmap
	#os.system('rm difmap.log*\n')

	#convolve_difmap([self.files_chosen[0]],[self.models_chosen[0]],self.bmaj_files,self.bmin_files,self.bpa_files,ShiftWindow.xshift,ShiftWindow.yshift,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[0]])

	#rms1 (lower freq)
	#self.rms1 = search_rms()[0]
	#os.system('rm difmap.log*\n')

	#convolve_difmap([self.files_chosen[1]],[self.models_chosen[1]],self.bmaj_files,self.bmin_files,self.bpa_files,0,0,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[1]])

	#rms2 (higher freq)
	#self.rms2 = search_rms()[0]
	#os.system('rm difmap.log*\n')

	self.write_result_file()

	#close the popup
	self.wi.close()


class needed_param():

	path = os.getcwd()

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

	units = []

	freq = freqOrig.copy()

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

	w = ShiftWindow()
	w.show()

	app.exec_()


