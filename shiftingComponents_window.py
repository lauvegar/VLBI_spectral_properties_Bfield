import threading, time
import sys
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from math import *
from functools import *
import numpy as np
import astropy.io.fits as pf
import pickle
from functions_conv import order_by_nu, read_conv_params
from functions_align import find_same_beam,beam_array,check_map_params, search_rms
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, Annotate
from functions_alignComp import natural_keys, read_modfile, ellipse_axis_lines, x_y, selectComponent, CoreShiftCalculation
import os,glob
import subprocess as sub
from astropy.nddata import Cutout2D
from correlate2d import *
#from fast_ftts import *

"""def search_rms():
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

	return rms"""

class popup_beams(QWidget):
	def __init__(self,parent=None,widget=None):
     	 	QWidget.__init__(self,parent)
	   	layout = QGridLayout(self)
		#layout.addStretch(1)
		#layout.addLayout(hbox)
     	  	self.checks2 = []

           	for i in xrange(0,len(ShiftComponentsWindow.diff_beams)):
          		c = QRadioButton(str('%1.3f' % (ShiftComponentsWindow.diff_beams[i]))+' mas')
	  		c.setFixedSize(100,25)
          		layout.addWidget(c,0,i+1)
          		self.checks2.append(c)

	   
        	self.selectButton = QPushButton("&Select")
		self.selectButton.setAutoDefault(True)
		
		layout.addWidget(self.selectButton, 1, trunc(len(ShiftComponentsWindow.diff_beams)/2.+1.))

		self.setLayout(layout)
		self.adjustSize()

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

class popupAutomaticSelection(QWidget):
	def __init__(self):
     	 	QWidget.__init__(self)
	   	layout = QGridLayout()
		#layout.addStretch(1)
		#layout.addLayout(hbox)
     	  	
		self.labelempty = QLabel()
		self.labelTEXT = QLabel()
		self.labelTEXT2 = QLabel()
          	
		self.labelTEXT.setText("The components for the rest of the maps have been selected automatically")
		self.labelTEXT.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
		self.labelTEXT2.setText("Are you happy with the automatic selection?")
		self.labelTEXT2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

		layout.addWidget(self.labelempty,3,0)
		layout.addWidget(self.labelTEXT,0,1,1,5)
		layout.addWidget(self.labelTEXT2,1,1,1,5)
	   
        	self.YesButton = QPushButton("&Yes")
        	self.NoButton = QPushButton("&No")
		self.YesButton.setAutoDefault(True)
		self.NoButton.setAutoDefault(True)
		
		layout.addWidget(self.YesButton, 3,2)
		layout.addWidget(self.NoButton, 3,3)

		#for the plot in the window 
		plt.ioff()
		self.plot = plt.figure()		
		self.canvas = FigureCanvas(self.plot)
		self.toolbar = NavigationToolbar(self.canvas,self)
		self.Plot_All(ShiftComponentsWindow.firstCont,ShiftComponentsWindow.pts_arr,
				ShiftComponentsWindow.x_el_arr, ShiftComponentsWindow.y_elH_arr, 
				ShiftComponentsWindow.x_elH_arr, ShiftComponentsWindow.y_el_arr,ShiftComponentsWindow.ext,
				ShiftComponentsWindow.realDATA,ShiftComponentsWindow.freq1,ShiftComponentsWindow.freq2,
				ShiftComponentsWindow.indexes)
		layout.addWidget(self.canvas,4,0,9,9)
		layout.addWidget(self.toolbar,15,0,1,9)

		#self.canvas = FigureCanvas(KinematicsWindow.plot)
       		#layout.addWidget(self.canvas)

		self.setLayout(layout)

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

	def Plot_All(self,first_contour,pts_arr,x_el_arr, y_elH_arr, x_elH_arr, y_el_arr,ext,realDATA,freq1,freq2,indexes):

		print 'indexes', indexes
		levels = first_contour[0]*np.array([-1., 1., 1.41,2.,2.83,4.,5.66,8.,11.3,16.,
			                        22.6,32.,45.3,64.,90.5,128.,181.,256.,362.,512.,
			                        724.,1020.,1450.,2050.])
		#plt.figure(1)
		self.plot.add_subplot(121)
		cset = plt.contour(realDATA[0], levels, inline=1,
				   colors=['grey'],
				   extent=ext, aspect=1.0
				   )
		for j in xrange(0,len(x_el_arr[0])):
			plt.plot(pts_arr[0][j][:,0], pts_arr[0][j][:,1], color='blue',linewidth=4)
			plt.plot(x_el_arr[0][j], y_elH_arr[0][j], color='blue',linewidth=4) 
			plt.plot(x_elH_arr[0][j], y_el_arr[0][j], color='blue',linewidth=4)
		for j in xrange(0,len(indexes[0])):
			plt.plot(pts_arr[0][indexes[0][j]][:,0], pts_arr[0][indexes[0][j]][:,1], color='red',linewidth=4)
			plt.plot(x_el_arr[0][indexes[0][j]], y_elH_arr[0][indexes[0][j]], color='red',linewidth=4) 
			plt.plot(x_elH_arr[0][indexes[0][j]], y_el_arr[0][indexes[0][j]], color='red',linewidth=4)
		plt.xlim(ext[0],ext[1])
		plt.ylim(ext[2],ext[3])	
		plt.axis('scaled')
		plt.xlabel('Right Ascension [pixels]')
		plt.ylabel('Relative Declination [pixels]')
		plt.title(str('%1.3f' %(freq1))+' GHz')

		levels = first_contour[1]*np.array([-1., 1., 1.41,2.,2.83,4.,5.66,8.,11.3,16.,
			                        22.6,32.,45.3,64.,90.5,128.,181.,256.,362.,512.,
			                        724.,1020.,1450.,2050.])

		#plt.figure(2)
		self.plot.add_subplot(122)
		cset = plt.contour(realDATA[1], levels, inline=1,
				   colors=['grey'],
				   extent=ext, aspect=1.0
				   )
		for j in xrange(0,len(x_el_arr[1])):
			plt.plot(pts_arr[1][j][:,0], pts_arr[1][j][:,1], color='blue',linewidth=4)
			plt.plot(x_el_arr[1][j], y_elH_arr[1][j], color='blue',linewidth=4) 
			plt.plot(x_elH_arr[1][j], y_el_arr[1][j], color='blue',linewidth=4)
		for j in xrange(0,len(indexes[0])):
			plt.plot(pts_arr[1][indexes[1][j]][:,0], pts_arr[1][indexes[1][j]][:,1], color='red',linewidth=4)
			plt.plot(x_el_arr[1][indexes[1][j]], y_elH_arr[1][indexes[1][j]], color='red',linewidth=4) 
			plt.plot(x_elH_arr[1][indexes[1][j]], y_el_arr[1][indexes[1][j]], color='red',linewidth=4)


		plt.xlim(ext[0],ext[1])
		plt.ylim(ext[2],ext[3])	
		plt.axis('scaled')
		plt.xlabel('Right Ascension [pixels]')
		plt.title(str('%1.3f' %(freq2))+' GHz')

		self.canvas.draw()

class ShiftComponentsWindow(QWidget):

    # declare class parameters that can be call in another classes
    diff_beams = []
    xshift = 0.
    yshift = 0.
    errxshift = 0.
    erryshift = 0.
    offset = [0.,0.]
    offset_new = []

    freq1 = 0.
    freq2 = 0.
    realDATA = []
    ext = []
    rms = []
    firstCont = []
    numComp = []

    pts_arr,x_el_arr,x_elH_arr,y_elH_arr,y_el_arr = [],[],[],[],[]

    xSelected, errxSelected = [], []
    ySelected, errySelected = [], []

    indexes = []

    def __init__(self,*args):
        QWidget.__init__(self)

        layout = QGridLayout()
	#layout.addStretch(1)
	#layout.addLayout(hbox)

	#declar class parameters to use inside the class
        self.checks = []
        self.buttons_freqs = []
	self.name_button_freq = []

	self.freq_conv = []
	self.beam_conv = []
	self.fits_conv = []

	self.fits1 = 'l'
	self.fits2 = 'r'
	self.freq1 = 0.
	self.freq2 = 0.
	self.rms1 = 0.
	self.rms2 = 0.
	self.files_chosen = []
	self.models_chosen = [] 
	self.shifted_files = []

	self.bmaj_files = 0.
	self.bmin_files = 0.
	self.bpa_files = 0.
	self.beam_files = 0.
	self.cells_files = 0.
	self.mapsize_file = 0.

	self.rms = []
	self.first_contour = []

	self.orientation = 'v'

        for i in xrange(0,len(needed_param.freq)):
	    if needed_param.units[i] == 'MHz':
           	 c = QCheckBox('%s %s' % ('%1.0f' % (needed_param.freq[i]), needed_param.units[i]),self)
	    elif needed_param.units[i] == 'GHz':
           	 c = QCheckBox('%s %s' % ('%1.2f' % (needed_param.freq[i]), needed_param.units[i]),self)
	    c.setFixedSize(100,25)
            layout.addWidget(c,1,i+1)
            self.checks.append(c)

        #for i in xrange(0,len(freq)):
        #    c = QPushButton("Beam "+str('%1.2f' % (freq[i])))
        #    layout.addWidget(c,1,i+1)
        #    self.buttons_freqs.append(c)

	#self.lineEditBMAJ = QLineEdit(self.gridLayoutWidget)
        #self.lineEditBMAJ.setObjectName(_fromUtf8("lineEditBMAJ"))

	self.labelempty = QLabel()
	self.labelempty.setFixedSize(25,25)
	self.labelempty2 = QLabel()

	#label and lineEdit for BMIN
  	self.labelnumComp = QLabel()
   	self.labelnumComp.setText("# Components: ")
	self.labelnumComp.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelnumComp.setFixedSize(100,25)
	self.numComp = QLineEdit()
	self.numComp.setValidator(QIntValidator())#QDoubleValidator(0.#min,3.#max,2.#number decimals))
	self.numComp.textChanged.connect(self.check_state)
	self.numComp.textChanged.emit(self.numComp.text())
	self.numComp.setFixedSize(100,25)
	self.numComp.setText('%d' % 1)

	self.labelnumComp.setBuddy(self.numComp)

	layout.addWidget(self.labelnumComp, 3,1)
	layout.addWidget(self.numComp, 3,2)

   	self.labelSigmaCut = QLabel()
   	self.labelSigmaCut.setText("SigmaCut Map:")
	self.labelSigmaCut.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelSigmaCut.setFixedSize(100,25)
	self.SigmaCut = QLineEdit()
	self.SigmaCut.setValidator(QDoubleValidator())
	self.SigmaCut.textChanged.connect(self.check_state)
	self.SigmaCut.textChanged.emit(self.SigmaCut.text())
	self.SigmaCut.setFixedSize(100,25)
	self.SigmaCut.setText('%1.1f' % (1.0))

	self.labelSigmaCut.setBuddy(self.SigmaCut)

	layout.addWidget(self.labelSigmaCut, 3,4)
	layout.addWidget(self.SigmaCut, 3,5)

	self.labelFITS1 = QLabel()
	self.labelFITS1.setFixedSize(100,25)
	self.labelFITS2 = QLabel()
	self.labelFITS2.setFixedSize(100,25)
	self.labelFITS1file = QLabel()
	#self.labelFITS1file.setFixedSize(100,25)
	self.labelFITS2file = QLabel()
	#self.labelFITS2file.setFixedSize(100,25)

	self.freqsSHIFTING = QLabel()
	#self.freqsSHIFTING.setFixedSize(100,25)
	self.labelSHIFTra = QLabel()
	self.labelSHIFTra.setFixedSize(100,25)
	self.labelSHIFTraMAS = QLabel()
	self.labelSHIFTraMAS.setFixedSize(100,25)
	self.labelSHIFTraPIXEL = QLabel()
	self.labelSHIFTraPIXEL.setFixedSize(100,25)
	self.labelSHIFTdec = QLabel()
	self.labelSHIFTdec.setFixedSize(100,25)
	self.labelSHIFTdecMAS = QLabel()
	self.labelSHIFTdecMAS.setFixedSize(100,25)
	self.labelSHIFTdecPIXEL = QLabel()
	self.labelSHIFTdecPIXEL.setFixedSize(100,25)

	self.labelOK = QLabel()

	self.labelOpticallyThick = QLabel()

	self.labelFITS1.setText("FITS1 : ")
	self.labelFITS1.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelFITS2.setText("FITS2 : ")
	self.labelFITS2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

	self.labelSHIFTra.setText("SHIFT RA : ")
	self.labelSHIFTra.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelSHIFTdec.setText("SHIFT DEC : ")
	self.labelSHIFTdec.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 0)	

	layout.addWidget(self.labelFITS1, 5,1)
	layout.addWidget(self.labelFITS2, 6,1)
	layout.addWidget(self.labelFITS1file, 5,2, 1, 5)
	layout.addWidget(self.labelFITS2file, 6,2, 1, 5)
	layout.addWidget(self.labelOK, 8,1,1,4)
	layout.addWidget(self.freqsSHIFTING, 10,1,1,5)
	layout.addWidget(self.labelSHIFTra, 11,1)
	layout.addWidget(self.labelSHIFTraMAS, 11,2)
	layout.addWidget(self.labelSHIFTraPIXEL, 11,3)
	layout.addWidget(self.labelSHIFTdec, 12,1)
	layout.addWidget(self.labelSHIFTdecMAS, 12,2)
	layout.addWidget(self.labelSHIFTdecPIXEL, 12,3)

	layout.addWidget(self.labelOpticallyThick, 14,1,3,5)

	self.labelFITS1.setBuddy(self.labelFITS1file)
	self.labelFITS2.setBuddy(self.labelFITS2file)
	self.labelSHIFTra.setBuddy(self.labelSHIFTraMAS)
	self.labelSHIFTdec.setBuddy(self.labelSHIFTdecMAS)
	self.labelSHIFTraMAS.setBuddy(self.labelSHIFTraPIXEL)
	self.labelSHIFTraMAS.setBuddy(self.labelSHIFTdecPIXEL)

	for i in xrange(len(self.checks)):
		self.checks[i].toggled.connect(lambda checked: self.checksState())

	#self.progressBar = QProgressBar(self)
        #self.progressBar.setRange(0,100)
        self.findingBEAMbutton = QPushButton("&Check")
	self.findingBEAMbutton.setFixedSize(100,25)
	self.findingBEAMbutton.setAutoDefault(True) #button clicked by pressing enter
        self.findingBEAMbutton.clicked.connect(lambda: self.findingBEAM(self.checks,needed_param.freqOrig))
        #self.selectCoreButton = QPushButton("&Select Core")
	#self.selectCoreButton.setAutoDefault(True)
	#self.selectCoreButton.setFixedSize(100,25)
        #self.selectCoreButton.clicked.connect(lambda: self.selectCore(self.fits1,self.fits2))
        self.selectButton = QPushButton("&Select \n Component")
	self.selectButton.setAutoDefault(True)
        self.selectButton.clicked.connect(lambda: self.selectComp(self.fits1,self.fits2))
        self.ShiftingButton = QPushButton("&Shifting")
        #self.ManualShiftingButton.clicked.connect(lambda: self.ManualShifting())
	self.ShiftingButton.setAutoDefault(True)
	self.ShiftingButton.setFixedSize(100,25)
        self.ShiftingButton.clicked.connect(lambda: self.shifting())

	#disable the buttons when opening the GUI
	self.findingBEAMbutton.setEnabled(False)
	#self.selectCoreButton.setEnabled(False)
	self.selectButton.setEnabled(False)
	self.ShiftingButton.setEnabled(False)


	if len(needed_param.freq) > 2:		
		layout.addWidget(self.findingBEAMbutton, 7, len(needed_param.freq)+2)
		#layout.addWidget(self.selectCoreButton, 8, len(needed_param.freq)+2)
		layout.addWidget(self.selectButton, 8, len(needed_param.freq)+2)
		layout.addWidget(self.ShiftingButton, 9, len(needed_param.freq)+2,2,1)
		for i in xrange(0,17):
			layout.addWidget(self.labelempty2, i, len(needed_param.freq)+3)	
	else:
		layout.addWidget(self.findingBEAMbutton, 7, len(needed_param.freq)+4)
		#layout.addWidget(self.selectCoreButton, 8, len(needed_param.freq)+2)
		layout.addWidget(self.selectButton, 8, len(needed_param.freq)+4)
		layout.addWidget(self.ShiftingButton, 9, len(needed_param.freq)+4,2,1)
		for i in xrange(0,17):
			layout.addWidget(self.labelempty2, i, len(needed_param.freq)+5)	

	#layout.addWidget(self.progressBar,17,1,1,len(freq))

        self.setLayout(layout)

	#put the window in the center of the desktop
	qr = self.frameGeometry()
	cp = QDesktopWidget().availableGeometry().center()
	qr.moveCenter(cp)
	self.move(qr.topLeft())

	self.setWindowTitle("Shifting")

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
	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			checked.append(i)
	if len(checked) == 2:
		if len(needed_param.modelfit) > 0:
			self.findingBEAMbutton.setEnabled(True)
			#self.selectCoreButton.setEnabled(True)
			self.selectButton.setEnabled(True)
			self.ShiftingButton.setEnabled(True)

	else:
		self.findingBEAMbutton.setEnabled(False)
		#self.selectCoreButton.setEnabled(False)
		self.selectButton.setEnabled(False)
		self.ShiftingButton.setEnabled(False)

    def getBEAM(self,RadioButtons,diff_beams,beam1,beam2,freq1_index,freq2_index,fits_conv):

			for i in xrange(len(RadioButtons)):
				if RadioButtons[i].isChecked():
					self.selected_beam = ShiftComponentsWindow.diff_beams[i]

			index_b1 = np.where(beam1 == self.selected_beam)
			index_b2 = np.where(beam2 == self.selected_beam) #output a list containing tuple
			#converting the list into array and the tuples in the list into an array
			index_beam1 = np.asarray([x for xs in index_b1 for x in xs]) 
			index_beam2 = np.asarray([x for xs in index_b2 for x in xs]) 
			index_fits1 = freq1_index[0][index_beam1]
			index_fits2 = freq2_index[0][index_beam2]

			self.fits1 = fits_conv[index_fits1[0]]
			self.fits2 = fits_conv[index_fits2[0]]

			self.labelFITS1file.setText(self.fits1[len(needed_param.path):])
			self.labelFITS2file.setText(self.fits2[len(needed_param.path):])

			see_if_OK = check_map_params(str(self.fits1),str(self.fits2),False)
			OK = see_if_OK[0]

			if OK == True:
				self.labelFITS1file.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
				self.labelFITS2file.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
				self.labelOK.setText("Files OK")
				self.labelOK.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
			else:
				self.labelFITS1file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
				self.labelFITS2file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
				self.labelOK.setText("Convolve the files with the same beam,cell and mapsize")
				self.labelOK.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')

			self.wi.close()


    def findingBEAM(self,checkBOXes,freq): #self.checks,freq ###stop_event,checkBOXes,freq 
	#self.myLongTask.start()

	self.fits_conv = []

	for filename in sorted(glob.glob(needed_param.path+'/CONV/*.fits*')):   
		self.fits_conv.append(filename)   
	files_temp = self.fits_conv
	models_temp = self.fits_conv


	###ORDER THEM BY FREQUENCY
	ordered_params_conv = order_by_nu(files_temp,files_temp,self.fits_conv,False)
	self.freq_conv = ordered_params_conv[0]
	self.beam_conv = ordered_params_conv[5]
	self.fits_conv = ordered_params_conv[10]

	finding_beam = beam_array(checkBOXes,freq,self.freq_conv,self.beam_conv)
	self.freq1_index,self.freq2_index,self.index_beam12 = finding_beam[0],finding_beam[1],finding_beam[2]
	self.beam1,self.beam2 = finding_beam[3],finding_beam[4]

	if len(self.index_beam12)== 1:
		track_freq1 = self.index_beam12[0][0]
		index_freq1 = self.freq1_index[0][track_freq1]
		self.fits1 = self.fits_conv[index_freq1]

		track_freq2 = self.index_beam12[0][1]-len(self.freq1_index[0])
		index_freq2 = self.freq2_index[0][track_freq2]
		self.fits2 = self.fits_conv[index_freq2]

		self.labelFITS1file.setText(self.fits1[len(needed_param.path):])
		self.labelFITS2file.setText(self.fits2[len(needed_param.path):])
			
		see_if_OK = check_map_params(str(self.fits1),str(self.fits2),False)
		OK = see_if_OK[0]

		if OK == True:
			self.labelFITS1file.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
			self.labelFITS2file.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
			self.labelOK.setText("Files OK")
			self.labelOK.setStyleSheet('QLabel {color: green } QLabel {font: Bold }')
		else:
			self.labelFITS1file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
			self.labelFITS2file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
			self.labelOK.setText("Convolve the files with the same beam,cell and mapsize")
			self.labelOK.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')

	elif len(self.index_beam12) > 1:
		ShiftComponentsWindow.diff_beams = []

		for i in xrange(0,len(self.index_beam12)):
			ShiftComponentsWindow.diff_beams.append(self.beam1[self.index_beam12[i][0]])
	
		self.wi = popup_beams()
		self.wi.show()
		self.wi.selectButton.clicked.connect(lambda: self.getBEAM(self.wi.checks2, ShiftComponentsWindow.diff_beams,self.beam1,self.beam2,self.freq1_index,self.freq2_index,self.fits_conv))

	elif len(self.index_beam12) == 0:
		self.labelOK.setText("No match found. Convolve the frequencies with the same beam")
		self.labelOK.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
		self.labelFITS1file.setText("X")
		self.labelFITS2file.setText("X")
		self.labelFITS1file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
		self.labelFITS2file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')

	self.update()

	return popup_beams

    def update(self):
	self.freqsSHIFTING.setText("")
	self.labelSHIFTraMAS.setText("")
	self.labelSHIFTraPIXEL.setText("")
	self.labelSHIFTdecMAS.setText("")
	self.labelSHIFTdecPIXEL.setText("")

    def selectComp(self,fits1,fits2):

	self.rms = []

	header1 = take_header(fits1,False)
	header2 = take_header(fits2,False)
	#self.rms.append(header1[9])
	#self.rms.append(header2[9])
	map_data1 = read_map(fits1,False)		
	realDAT = map_data1[0]
	map_data2 = read_map(fits2,False)		
	realDAT2 = map_data2[0]


	#picks the value of the first contour for plotting the maps
	ShiftComponentsWindow.realDATA = [realDAT,realDAT2]

	ShiftComponentsWindow.numComp = int(self.numComp.text())


	#obtaining the beam and cell
	self.bmaj_files = header1[1]
	self.bmin_files = header1[2]
	self.bpa_files = header1[3]
	self.beam_files = header1[7]
	self.cells_files = header1[0]

	#obtaining frequencies
	self.freq1 = header1[5]
	self.freq2 = header2[5]


	if self.freq1 < 0.5:
		self.freq1name = str('%1.0f' %(self.freq1*1000))
		self.freq1unit = 'MHz'
	else:
		self.freq1name = str('%1.2f' %(self.freq1))
		self.freq1unit = 'GHz'

	if self.freq2 < 0.5:
		self.freq2name = str('%1.0f' %(self.freq2*1000))
		self.freq2unit = 'MHz'
	else:
		self.freq2name = str('%1.2f' %(self.freq2))
		self.freq2unit = 'GHz'

	res=open(needed_param.path+'/rms'+self.freq1name+'.p','rb')
	pick = pickle.load(res)
	res.close()  

	self.rms.append(pick) 

	res=open(needed_param.path+'/rms'+self.freq2name+'.p','rb')
	pick = pickle.load(res)
	res.close()  

	self.rms.append(pick) 

	ShiftComponentsWindow.freq1 = self.freq1
	ShiftComponentsWindow.freq2 = self.freq2

	self.first_contour = []
	sigma_cut = float(self.SigmaCut.text())
	for i in xrange(0,len(self.rms)):
		self.first_contour.append(sigma_cut*self.rms[i])

	ShiftComponentsWindow.firstCont = self.first_contour

	#obtaining map centers in pixels
	cent_mapx = map_data1[5]
	cent_mapy = map_data1[6]
	self.mapsize_files = 2*map_data1[7]

	#obtaining the four corners of the maps in mas
	x1 = map_data1[1]
	x2 = map_data1[2]
	y1 = map_data1[3]
	y2 = map_data1[4]

	ShiftComponentsWindow.ext=[x1,x2,y1,y2] 

	#for plotting the components
	ShiftComponentsWindow.pts_arr,ShiftComponentsWindow.x_el_arr,ShiftComponentsWindow.x_elH_arr = [],[],[]
	ShiftComponentsWindow.y_elH_arr,ShiftComponentsWindow.y_el_arr = [],[]
	ShiftComponentsWindow.xSelected,ShiftComponentsWindow.ySelected = [],[]
	ShiftComponentsWindow.errxSelected,ShiftComponentsWindow.errySelected = [],[]
	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			ShiftComponentsWindow.pts_arr.append(needed_param.pts_arr[i])
			ShiftComponentsWindow.x_el_arr.append(needed_param.x_el_arr[i])
			ShiftComponentsWindow.x_elH_arr.append(needed_param.x_elH_arr[i])
			ShiftComponentsWindow.y_elH_arr.append(needed_param.y_elH_arr[i])
			ShiftComponentsWindow.y_el_arr.append(needed_param.y_el_arr[i]) 
			ShiftComponentsWindow.xSelected.append(needed_param.x[i]) 
			ShiftComponentsWindow.ySelected.append(needed_param.y[i]) 
			ShiftComponentsWindow.errxSelected.append(needed_param.errx[i]) 
			ShiftComponentsWindow.errySelected.append(needed_param.erry[i]) 

	indexComp = selectComponent(realDAT,realDAT2,self.first_contour, ShiftComponentsWindow.pts_arr,ShiftComponentsWindow.x_el_arr,
			ShiftComponentsWindow.x_elH_arr,ShiftComponentsWindow.y_elH_arr,ShiftComponentsWindow.y_el_arr,
			ShiftComponentsWindow.ext,self.freq1,self.freq2,
			ShiftComponentsWindow.xSelected,ShiftComponentsWindow.ySelected,ShiftComponentsWindow.numComp,self.orientation)

	ShiftComponentsWindow.indexes = [indexComp[0],indexComp[1]]
	
	plt.close('all')
	
	self.wi = popupAutomaticSelection()
	self.wi.show()
	self.wi.YesButton.clicked.connect(lambda: self.GoodCompSelection())
        self.wi.NoButton.clicked.connect(lambda: self.BadCompSelection())

    def GoodCompSelection(self):

	plt.close('all')

	ShiftValues = CoreShiftCalculation(ShiftComponentsWindow.indexes,ShiftComponentsWindow.xSelected,ShiftComponentsWindow.ySelected,
					ShiftComponentsWindow.errxSelected,ShiftComponentsWindow.errySelected,ShiftComponentsWindow.numComp)
	ShiftComponentsWindow.xshift = ShiftValues[0]
	ShiftComponentsWindow.yshift = ShiftValues[1]
	ShiftComponentsWindow.errxshift = ShiftValues[2]
	ShiftComponentsWindow.erryshift = ShiftValues[3]
	ShiftComponentsWindow.offset[0]=-ShiftComponentsWindow.xshift/self.cells_files #pix
	ShiftComponentsWindow.offset[1]=ShiftComponentsWindow.yshift/self.cells_files #pix

	self.freqsSHIFTING.setText(("Image2 ( %s GHz)  - image1 ( %s GHz)" % ('%1.2f' % (self.freq2), '%1.2f' % (self.freq1))))
	self.freqsSHIFTING.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTraMAS.setText(("%s mas" % ('%1.4f' % (ShiftComponentsWindow.xshift))))
	self.labelSHIFTraMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTraPIXEL.setText(("%s pix" % ('%1.2f' % (ShiftComponentsWindow.offset[0]))))
	self.labelSHIFTraPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecMAS.setText(("%s mas" % ('%1.4f' % (ShiftComponentsWindow.yshift))))
	self.labelSHIFTdecMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecPIXEL.setText(("%s pix" % ('%1.2f' % (ShiftComponentsWindow.offset[1]))))
	self.labelSHIFTdecPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')

	a = []

	track = 0


	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			if track == 0:
				flux1 = []
				for j in xrange(0,len(ShiftComponentsWindow.indexes[0])):
					flux1.append(needed_param.flux[i][ShiftComponentsWindow.indexes[0][j]])

				track = track + 1

			if track == 1:
				flux2 = []
				for j in xrange(0,len(ShiftComponentsWindow.indexes[1])):
					flux2.append(needed_param.flux[i][ShiftComponentsWindow.indexes[1][j]])

	a = np.zeros(np.shape(flux1))


	for i in xrange(0,len(flux1)):
		a[i] = (np.log10(flux1[i]/flux2[i])/np.log10(self.freq1/self.freq2))


	s = None
	if np.any(a > 0.)==True:
		for i in xrange(0,len(a)):
			if s != None:
				s = s+','+str(i+1)
			if s == None:
				s = str(i+1)

		#self.labelOpticallyThick.setText('The components '+str(s)+', in the order of selection, are optically THICK.')
		#self.labelOpticallyThick.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')

	#self.labelOptically thin or not

	self.wi.close()

    def BadCompSelection(self):
	print 'hi'

    def shifting(self):
	self.shifted_files=[]
	self.shifted_files.append(needed_param.path+'/SHIFT/'+needed_param.source_name+'-'+str(round(self.freq1))+'convolved_with_beam'+str('%1.2f' % (self.beam_files))+'shifted.fits')
	self.shifted_files.append(needed_param.path+'/SHIFT/'+needed_param.source_name+'-'+str(round(self.freq2))+'convolved_with_beam'+str('%1.2f' % (self.beam_files))+'shifted.fits')
	

	self.models_chosen = []
	self.files_chosen = []
	for i in xrange(len(self.checks)):
		if self.checks[i].isChecked():
			self.models_chosen.append(needed_param.models[i])
			self.files_chosen.append(needed_param.files[i])

	os.system('rm difmap.log*\n')

	convolve_difmap([self.files_chosen[0]],[self.models_chosen[0]],self.bmaj_files,self.bmin_files,self.bpa_files,ShiftComponentsWindow.xshift,ShiftComponentsWindow.yshift,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[0]])

	#rms1 (lower freq)
	self.rms1 = search_rms()[0]
	os.system('rm difmap.log*\n')

	convolve_difmap([self.files_chosen[1]],[self.models_chosen[1]],self.bmaj_files,self.bmin_files,self.bpa_files,0,0,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[1]])

	self.rms2 = search_rms()[0]
	os.system('rm difmap.log*\n')

	self.write_result_file()


    def write_result_file(self):
	final = [ShiftComponentsWindow.xshift, ShiftComponentsWindow.yshift, self.rms1, self.rms2 ,ShiftComponentsWindow.ext]
	res=open(needed_param.path+'/Shift_parameters/shift_param'+str(int(round(self.freq1)))+'and'+str(int(round(self.freq2)))+'.p','wb')
	pickle.dump(final,res)
	res.close()     

	final_txt = [[ShiftComponentsWindow.xshift, ShiftComponentsWindow.errxshift, ShiftComponentsWindow.yshift, ShiftComponentsWindow.erryshift, self.rms1, self.rms2,ShiftComponentsWindow.ext[0],ShiftComponentsWindow.ext[1],ShiftComponentsWindow.ext[2],ShiftComponentsWindow.ext[3]]]
	header = np.array([['#RAshift errRAshift DECshift errDECshift rms(low freq image) rms(high freq image) coordinates of the map window -----> lower freq shifted']])
	saver(needed_param.path+'/Shift_parameters/shift_param'+str(round(self.freq1))+'and'+str(round(self.freq2))+'.txt', header, final_txt, FORMAT='%1.6f')

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

	modelfit = []
	for filename in sorted(glob.glob(path+'/modelfit/*.mod*')):   
		modelfit.append(filename)  

	modelfit.sort(key=natural_keys)
					
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
	beam = ordered_params[5]
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

	#reads the modelfit files and obtains the modelfit values and errors

	if len(modelfit) > 0:
		if os.path.isfile('pos_errors.dat'):
			errors = True
		else:
			errors = None

		mod_parameters = read_modfile(modelfit,beam,errors)
		"""
		r, errr = radial distance and error of the component
		psi, errpsi = position angle and error of the component
		size, errsize =  size and error of the component
		"""
		r, errr = mod_parameters[0], mod_parameters[1]     
		psi, errpsi = mod_parameters[2], mod_parameters[3]     
		size, errsize = mod_parameters[4], mod_parameters[5]     
		flux, errflux = mod_parameters[7], mod_parameters[8]    

		#with the radial distance and position angle, the central positions of the components in RA and DEC are calculated
		"""
		x, errx = position in RA and error of the component
		y, erry = position in DEC and error of the component
		"""
		x_and_y = x_y(r,errr,psi,errpsi,errors)
		x, errx = np.asarray(x_and_y[0]), np.asarray(x_and_y[1])
		y, erry = np.asarray(x_and_y[2]), np.asarray(x_and_y[3])

		#for plotting the components in the map
		"""
		pts_arr = points for drawing the external countour of the ellipse, i.e., the ellipse itself

		x_el_arr = points for the x axis in the x direction of the ellipse. They are between (x_cent_component - size component) and (x_cent_component + size component).They are a total of 50
		y_elH_arr = points for the y axis in the x direction of the ellipse. It is a constant, so it is the same value 50 times, for using it with x_el_arr 

		y_el_arr = points for the y axis in the y direction of the ellipse. They are between (y_cent_component - size component) and (y_cent_component + size component). They are a total of 50
		x_elH_arr = points for the x axis in the y direction of the ellipse. It is a constant, so it is the same value 50 times, for using it with y_el_arr 

		"""
		pts_arr=[]
		pt_arr=[]
		x_el_arr=[]
		x_elH_arr=[]
		y_el_arr=[]
		y_elH_arr=[]

		ellipse_plot = ellipse_axis_lines(x,y,size)
		pts_arr,pt_arr = ellipse_plot[0], ellipse_plot[1]
		x_el_arr,y_el_arr = ellipse_plot[2], ellipse_plot[3]
		x_elH_arr,y_elH_arr = ellipse_plot[4], ellipse_plot[5]  	



def main():
	app = QApplication(sys.argv)

	w = ShiftComponentsWindow()
	w.show()

	app.exec_()

#main()
