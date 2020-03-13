import threading, time
import sys
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from math import *
from functools import *
import numpy as np
import astropy.io.fits as pf
from pylab import *
import pickle
from functions_conv import order_by_nu, read_conv_params
from functions_align import find_same_beam,beam_array,check_map_params, cuttingMAP,cross_correlation_shifts_FITS, checking_shift
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, get_ellipse_coords, Annotate
import os,glob
import subprocess as sub
from astropy.nddata import Cutout2D
from correlate2d import *
from scipy.fftpack import fftfreq
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial
import h5py, gc, psutil
import scipy.stats as stats
#import pyximport; pyximport.install()
#pyximport.install(pyimport = True)
#import pixelerror

import matplotlib as mpl
from packaging.version import Version, LegacyVersion

if Version(mpl.__version__) > Version('2.0.0'):
	mpl.style.use('classic')
	mpl.rc('image', cmap='jet')


#from fast_ftts import *

def deep_getsizeof(o, ids):
    """Find the memory footprint of a Python object
 
    This is a recursive function that drills down a Python object graph
    like a dictionary holding nested dictionaries with lists of lists
    and tuples and sets.
 
    The sys.getsizeof function does a shallow size of only. It counts each
    object inside a container as pointer only regardless of how big it
    really is.
 
    :param o: the object
    :param ids:
    :return:
    """
    d = deep_getsizeof
    if id(o) in ids:
        return 0
 
    r = sys.getsizeof(o)
    ids.add(id(o))
 
    if isinstance(o, str) or isinstance(0, unicode):
        return r
 
    if isinstance(o, Mapping):
        return r + sum(d(k, ids) + d(v, ids) for k, v in o.iteritems())
 
    if isinstance(o, Container):
        return r + sum(d(x, ids) for x in o)
 
    return r 


def shiftArr(x,y):
    shiftsArr = [x,y]
    return shiftsArr

#lineedit = QLineEdit(self)
#validator = QDoubleValidator()
#lineedit.setValidator(QDoubleValidator())
#lineedit.textChanged.connect(self.check_state)
#lineedit.textChanged.emit(lineedit.text())


def get_line(x1, y1, x2, y2):
    points = []
    issteep = abs(y2-y1) > abs(x2-x1)
    if issteep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2
    rev = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        rev = True
    deltax = x2 - x1
    deltay = abs(y2-y1)
    error = int(deltax / 2)
    y = y1
    ystep = None
    if y1 < y2:
        ystep = 1
    else:
        ystep = -1
    for x in range(x1, x2 + 1):
        if issteep:
            points.append((y, x))
        else:
            points.append((x, y))
        error -= deltay
        if error < 0:
            y += ystep
            error += deltax
    # Reverse the list if the coordinates were reversed
    if rev:
        points.reverse()
    return points


def div0( a, b ):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide( a, b )
        c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c

class LineDrawer(QWidget):
    lines = []
    x = []
    y = []
    def draw_line(self):
        ax = plt.gca()
	#plt.ginput(1)
        xy = plt.ginput(2)

        x = [p[0] for p in xy]
        y = [p[1] for p in xy]
        line = plt.plot(x,y,'r-')
        #ax.figure.canvas.draw()
	
	LineDrawer.x = x
	LineDrawer.y = y
        self.lines.append(line)

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

class SpixWindow(QWidget):

    proc = psutil.Process(os.getpid())
    gc.collect()

    diff_beams = []

    def __init__(self,*args):
        QWidget.__init__(self)

        layout = QGridLayout()
	#layout.addStretch(1)
	#layout.addLayout(hbox)

	#initializing parameters
        self.checks = []
        self.buttons_freqs = []
	self.name_button_freq = []

	self.fits1 = 'l'
	self.fits2 = 'r'
	self.freq1 = 0.
	self.freq2 = 0.
	self.freq1name = ''
	self.freq2name = ''
	self.freq1unit = ''
	self.freq2unit = ''
	self.files_chosen = []
	self.models_chosen = [] 
	self.shifted_files = []

	self.bmaj_files = 0.
	self.bmin_files = 0.
	self.bpa_files = 0.
	self.beam_files = 0.
	self.cells_files = 0.
	self.mapsize_file = 0.
	self.cent_mapx = 0.
	self.cent_mapy = 0.

	self.ifhdu = False

	self.rms1 = 0.
	self.rms2 = 0.
	self.ext = []
	self.extContours = []

	self.limplot_x1 = 0.
	self.limplot_x2 = 0.
	self.limplot_y1 = 0.
	self.limplot_y2 = 0.

	self.data1 = np.asarray([])
	self.data2 = np.asarray([])
	self.spix = np.asarray([])
	self.spixMasked = np.asarray([])
	self.cutoutSpix = np.asarray([])
	self.up = np.asarray([])
	self.down = np.asarray([])
	self.sigma_cut = 0.
	self.first_contour = 0.
	self.vmin = None
	self.vmax = None

	self.position = []
	self.size = []

	self.i = 0
	self.x = []
	self.ERRcheck = False

	self.aMasked = np.asarray([])

	self.shiftsX = np.asarray([])

	self.shiftsY = np.asarray([])

	self.shiftALL = np.asarray([])

	self.shiftPartly = np.asarray([])

	self.temp = np.asarray([])

	self.len1 = 0
	self.len2 = 0
	self.aaa = np.asarray([])
	self.a1 = np.asarray([])

	self.totalLines = 0

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

	self.withshift = QRadioButton('shifted')
	self.withoutshift = QRadioButton('unshifted')

	self.labelFITS1 = QLabel()
	self.labelFITS1.setFixedSize(100,25)
	self.labelFITS2 = QLabel()
	self.labelFITS2.setFixedSize(100,25)
	self.labelFITS1file = QLabel()
	self.labelFITS2file = QLabel()


  	self.labelnumLines = QLabel()
   	self.labelnumLines.setText("# Lines: ")
	self.labelnumLines.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelnumLines.setFixedSize(100,25)
	self.numLines = QLineEdit()
	self.numLines.setValidator(QIntValidator())#QDoubleValidator(0.#min,3.#max,2.#number decimals))
	self.numLines.textChanged.connect(self.check_state)
	self.numLines.textChanged.emit(self.numLines.text())
	self.numLines.setFixedSize(100,25)
	self.numLines.setText('%d' % 1)

	self.labelnumLines.setBuddy(self.numLines)

	self.labelOK = QLabel()

	self.labelFITS1.setText("FITS1 : ")
	self.labelFITS1.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.labelFITS2.setText("FITS2 : ")
	self.labelFITS2.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')


	self.labelRMS1 = QLabel()
	#self.labelRMS1.setFixedSize(100,25)
	self.labelRMS2 = QLabel()
	#self.labelRMS2.setFixedSize(100,25)

	self.labelSigmaCut = QLabel()
	self.labelSigmaCut.setText("Sigma Cut    : ")
	self.labelSigmaCut.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.SigmaCut = QLineEdit()
	self.SigmaCut.setValidator(QDoubleValidator())
	self.SigmaCut.textChanged.connect(self.check_state)
	self.SigmaCut.textChanged.emit(self.SigmaCut.text())
	self.SigmaCut.setFixedSize(100,25)

	self.labelVMIN = QLabel()
	self.labelVMIN.setText("Min scale   : ")
	self.labelVMIN.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.VMIN = QLineEdit()
	self.VMIN.setValidator(QDoubleValidator())
	self.VMIN.textChanged.connect(self.check_state)
	self.VMIN.textChanged.emit(self.VMIN.text())
	self.VMIN.setFixedSize(100,25)

	self.labelVMAX = QLabel()
	self.labelVMAX.setText("Max scale   : ")
	self.labelVMAX.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.VMAX = QLineEdit()
	self.VMAX.setValidator(QDoubleValidator())
	self.VMAX.textChanged.connect(self.check_state)
	self.VMAX.textChanged.emit(self.VMAX.text())
	self.VMAX.setFixedSize(100,25)


	self.labelErrText = QLabel()
	self.labelErrText.setText("If scale includes more values")
	self.labelErrVMIN = QLabel()
	self.labelErrVMIN.setText("MinErr scale   : ")
	self.labelErrVMIN.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.ErrVMIN = QLineEdit()
	self.ErrVMIN.setValidator(QDoubleValidator())
	self.ErrVMIN.textChanged.connect(self.check_state)
	self.ErrVMIN.textChanged.emit(self.ErrVMIN.text())
	self.ErrVMIN.setFixedSize(100,25)

	self.labelErrVMAX = QLabel()
	self.labelErrVMAX.setText("MaxErr scale   : ")
	self.labelErrVMAX.setStyleSheet('QLabel {color: black } QLabel {font: Bold }')
	self.ErrVMAX = QLineEdit()
	self.ErrVMAX.setValidator(QDoubleValidator())
	self.ErrVMAX.textChanged.connect(self.check_state)
	self.ErrVMAX.textChanged.emit(self.ErrVMAX.text())
	self.ErrVMAX.setFixedSize(100,25)



	for i in xrange(0,17):
		layout.addWidget(self.labelempty, i, 0)	

	layout.addWidget(self.labelnumLines, 3,1)
	layout.addWidget(self.numLines, 3,2)

	layout.addWidget(self.labelFITS1, 4,1)
	layout.addWidget(self.labelFITS2, 5,1)
	layout.addWidget(self.labelFITS1file, 4,2, 1, 8)
	layout.addWidget(self.labelFITS2file, 5,2, 1, 8)
	layout.addWidget(self.labelOK, 7,1,1,4)
	layout.addWidget(self.labelRMS1, 9,1,1,5)
	layout.addWidget(self.labelRMS2, 10,1,1,5)
	layout.addWidget(self.labelSigmaCut, 11,1)
	layout.addWidget(self.SigmaCut, 11,2)
	layout.addWidget(self.labelVMIN, 12,1)
	layout.addWidget(self.VMIN, 12,2)
	layout.addWidget(self.labelVMAX, 13,1)
	layout.addWidget(self.VMAX, 13,2)
	layout.addWidget(self.labelErrVMIN, 12,3)
	layout.addWidget(self.ErrVMIN, 12,4)
	layout.addWidget(self.labelErrVMAX, 13,3)
	layout.addWidget(self.ErrVMAX, 13,4)

	layout.addWidget(self.withshift, 2, trunc(len(needed_param.freq)/2))
	layout.addWidget(self.withoutshift, 2, trunc(len(needed_param.freq)/2)+1)

	self.labelSigmaCut.setBuddy(self.SigmaCut)
	self.labelVMIN.setBuddy(self.VMIN)
	self.labelVMAX.setBuddy(self.VMAX)
	self.labelErrVMIN.setBuddy(self.ErrVMIN)
	self.labelErrVMAX.setBuddy(self.ErrVMAX)

	"""	for i in xrange(len(self.checks)):	
		if self.checks[i].isChecked():
			self.tracker.append(i)
			print len(self.tracker)
	if len(self.tracker) == 0:
		setPIXELandMAP.setEnabled(False)
	else:
		setPIXELandMAP.setEnabled(True)"""
	#self.progressBar = QProgressBar(self)
        #self.progressBar.setRange(0,100)

	#for disable or enable buttons
	for i in xrange(len(self.checks)):
		self.checks[i].toggled.connect(lambda checked: self.checksState())


        self.findingBEAMbutton = QPushButton("&Check")
	self.findingBEAMbutton.setFixedSize(100,25)
        self.findingBEAMbutton.clicked.connect(lambda: self.findingBEAM(self.checks,needed_param.freqOrig))
	self.findingBEAMbutton.setAutoDefault(True)
        self.SPIXbutton = QPushButton("&Spix Map")
	self.SPIXbutton.setFixedSize(100,25)
        self.SPIXbutton.clicked.connect(lambda: self.spixCalculation())
	self.SPIXbutton.setAutoDefault(True)
        self.ErrSPIXbutton = QPushButton("&Error Spix")
	self.ErrSPIXbutton.setFixedSize(100,25)
        self.ErrSPIXbutton.clicked.connect(lambda: self.ErrorSpixMap())
	self.ErrSPIXbutton.setAutoDefault(True)
        self.CutSPIXbutton = QPushButton("&Cut Spix")
	self.CutSPIXbutton.setFixedSize(100,25)
        self.CutSPIXbutton.clicked.connect(lambda: self.cutSpix())
	self.CutSPIXbutton.setAutoDefault(True)

	#disable the buttons when opening the GUI
	self.findingBEAMbutton.setEnabled(False)
	self.SPIXbutton.setEnabled(False)
	self.ErrSPIXbutton.setEnabled(False)
	self.CutSPIXbutton.setEnabled(False)

	if len(needed_param.freq) > 2:		
		layout.addWidget(self.findingBEAMbutton, 6, len(needed_param.freq)+1)
		layout.addWidget(self.SPIXbutton, 7, len(needed_param.freq)+1)
		layout.addWidget(self.ErrSPIXbutton, 8, len(needed_param.freq)+1)
		layout.addWidget(self.CutSPIXbutton, 9, len(needed_param.freq)+1)
		for i in xrange(0,17):
			layout.addWidget(self.labelempty2, i, len(needed_param.freq)+2)	
	else:
		layout.addWidget(self.findingBEAMbutton, 6, len(needed_param.freq)+4)
		layout.addWidget(self.SPIXbutton, 7, len(needed_param.freq)+4)
		layout.addWidget(self.ErrSPIXbutton, 8, len(needed_param.freq)+4)
		layout.addWidget(self.CutSPIXbutton, 9, len(needed_param.freq)+4)
		for i in xrange(0,17):
			layout.addWidget(self.labelempty2, i, len(needed_param.freq)+4)	

	#layout.addWidget(self.progressBar,17,1,1,len(freq))

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

    #disables buttons
    def checksState(self):
	checked = []
	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			checked.append(i)
	if len(checked) == 2:
		self.findingBEAMbutton.setEnabled(True)
		self.SPIXbutton.setEnabled(True)
		self.ErrSPIXbutton.setEnabled(True)
		self.CutSPIXbutton.setEnabled(True)
	else:
		self.findingBEAMbutton.setEnabled(False)
		self.SPIXbutton.setEnabled(False)
		self.ErrSPIXbutton.setEnabled(False)
		self.CutSPIXbutton.setEnabled(False)

	#if self.ERRcheck == True:
	#	self.ErrSPIXbutton.setEnabled(True)
	#	self.CutSPIXbutton.setEnabled(True)
	#else:
	#	self.ErrSPIXbutton.setEnabled(False)
	#	self.CutSPIXbutton.setEnabled(False)



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
    fits_shift = self.fits_shift ---> list with all the fits files shifted,
						created and given from findingBEAM

    self.selected_beam ---> float, value of the beam selected by the user

    self.fits1,self.fits2 --->  str, final two fits files selected
    self.labelFITS1file, self.labelFITS2file, self.labelOK --->  to update the interface 

    """
    def getBEAM(self,RadioButtons,diff_beams,beam1,beam2,freq1_index,freq2_index,fits_shift):

			for i in xrange(len(RadioButtons)):
				if RadioButtons[i].isChecked():
					self.selected_beam = SpixWindow.diff_beams[i]

			index_b1 = np.where(beam1 == self.selected_beam)
			index_b2 = np.where(beam2 == self.selected_beam) #output a list containing tuple
			#converting the list into array and the tuples in the list into an array
			index_beam1 = np.asarray([x for xs in index_b1 for x in xs]) 
			index_beam2 = np.asarray([x for xs in index_b2 for x in xs]) 
			index_fits1 = freq1_index[0][index_beam1]
			index_fits2 = freq2_index[0][index_beam2]

			self.fits1 = fits_shift[index_fits1[0]]
			self.fits2 = fits_shift[index_fits2[0]]

			self.labelFITS1file.setText(self.fits1[len(needed_param.path):])
			self.labelFITS2file.setText(self.fits2[len(needed_param.path):])

			see_if_OK = check_map_params(str(self.fits1),str(self.fits2),self.ifhdu)
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


    """
    #function to find the fits files of the selected frequencies that have a valid beam
    #if there is only one valid beam (same beam for both of them) for the pair of frequencies, they get selected
    #if there is no valid beam, tells the user to convolve
    #if there is more than one common beam, a pop up appears for the user to choose the beam
    parameters used:
    checkBOXes = self.checks ---> list with the check boxes
    freq = needed_param.freq ---> list with all the existing frequencies 
    self.fits_shift ---> list with all the fits files shifted, redefined at the beginning (set to []), to allow the user to use the function again
    
    freq_shift ---> list with the frequencies of all fits files shifted, redefined at the beginning (by means of another function)
    beam_shift ---> list with the beams of all fits files shifted, redefined at the beginning (by means of another function)
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

	fits_shift = []

	if self.withshift.isChecked():
		for filename in sorted(glob.glob(needed_param.path+'/SHIFT/*.fits*')):   
			fits_shift.append(filename)   
		files_shift = fits_shift
		models_shift = fits_shift
		self.ifhdu = True

	if self.withoutshift.isChecked():
		for filename in sorted(glob.glob(needed_param.path+'/CONV/*.fits*')):   
			fits_shift.append(filename)   
		files_shift = fits_shift
		models_shift = fits_shift
		self.ifhdu = False

	#orders the shifted fits file by frequency (lower to higher)
	ordered_params_shift = order_by_nu(files_shift,files_shift,fits_shift,self.ifhdu)
	freq_shift = ordered_params_shift[0]
	beam_shift = ordered_params_shift[5]
	fits_shift = ordered_params_shift[10]

	#creates an array with beams and indexes of those beams for the two selected frequencies
	#finds the common beams for both frequencies
	#stores the index of the common beams for each subarray in an array, e.g. if two valid beams are found, indexbeam12 = [[0,2],[1,3]]
	#[0,2] ---> indexes in the subarray containing the beams for both frequencies for the first valid beam (in the example)
	#[1,3] ---> indexes in the subarray for the second valid beam (in the example)
	finding_beam = beam_array(checkBOXes,freq,freq_shift,beam_shift)
	self.freq1_index,self.freq2_index,self.index_beam12 = finding_beam[0],finding_beam[1],finding_beam[2]
	self.beam1,self.beam2 = finding_beam[3],finding_beam[4]

	#len(self.index_beam12)== 1 ----> only one common beam found
	#the fits files corresponding to that beam are stored, the parameters of the file checked
	#and the terminal updated
	if len(self.index_beam12)== 1:
		track_freq1 = self.index_beam12[0][0]
		index_freq1 = self.freq1_index[0][track_freq1]
		self.fits1 = fits_shift[index_freq1]

		track_freq2 = self.index_beam12[0][1]-len(self.freq1_index[0])
		index_freq2 = self.freq2_index[0][track_freq2]
		self.fits2 = fits_shift[index_freq2]

		self.labelFITS1file.setText(self.fits1[len(needed_param.path):])
		self.labelFITS2file.setText(self.fits2[len(needed_param.path):])
			
		see_if_OK = check_map_params(str(self.fits1),str(self.fits2),self.ifhdu)
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

		self.SigmaCut.setText('%1.1f' % (1.0))
		self.VMIN.setText('')
		self.VMAX.setText('')
		self.ErrVMIN.setText('')
		self.ErrVMAX.setText('')

	#len(self.index_beam12) > 1 ----> more than one common beam found
	#the user selects the desired beam using a popup window
	#it calls the popup class popup_beams() and the function of this class getBEAM controlling the button in the popup
	elif len(self.index_beam12) > 1:
		SpixWindow.diff_beams = []

		for i in xrange(0,len(self.index_beam12)):
			SpixWindow.diff_beams.append(self.beam1[self.index_beam12[i][0]])
	
		self.wi = popup_beams()
		self.wi.show()
		self.wi.selectButton.clicked.connect(lambda: self.getBEAM(self.wi.checks2, SpixWindow.diff_beams,self.beam1,self.beam2,self.freq1_index,self.freq2_index,fits_shift))
		self.SigmaCut.setText('%1.1f' % (1.0))
		self.VMIN.setText('')
		self.VMAX.setText('')
		self.ErrVMIN.setText('')
		self.ErrVMAX.setText('')


	#len(self.index_beam12) == 0 ----> no valid files beam related found
	#two X appear in the fits1/2 interface widget and tells the user to convolve the files in OK label
	elif len(self.index_beam12) == 0:
		self.labelOK.setText("No match found. Convolve the frequencies with the same beam")
		self.labelOK.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
		self.labelFITS1file.setText("X")
		self.labelFITS2file.setText("X")
		self.labelFITS1file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')
		self.labelFITS2file.setStyleSheet('QLabel {color: red } QLabel {font: Bold }')



    def reading_rms_param(self):


	header1 = take_header(self.fits1,self.ifhdu)
	header2 = take_header(self.fits2,self.ifhdu)

	#obtaining frequencies
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

	res=open(needed_param.path+'/Shift_parameters/shift_param'+self.freq1name+'and'+self.freq2name+'.p','rb')
	pick = pickle.load(res)
	res.close() 

	self.rms1 = pick[2]
	self.rms2 = pick[3]
	self.ext = pick[4]

	self.labelRMS1.setText(("rms noise Image1 ( %s GHz) : %s mJy/beam" % ('%1.3f' % (self.freq1), '%1.4f' % (self.rms1*1000))))
	self.labelRMS1.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelRMS2.setText(("rms noise Image2 ( %s GHz) : %s mJy/beam" % ('%1.3f' % (self.freq2), '%1.4f' % (self.rms2*1000))))
	self.labelRMS2.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')

    def spixCalculation(self):

	self.reading_rms_param()
	header1 = take_header(self.fits1,self.ifhdu)
	header2 = take_header(self.fits2,self.ifhdu)
	map_data1 = read_map(self.fits1,self.ifhdu)		
	realDAT = map_data1[0]
	map_data2 = read_map(self.fits2,self.ifhdu)		
	realDAT2 = map_data2[0]

	#obtaining the beam and cell
	self.bmaj_files = header1[1]
	self.bmin_files = header1[2]
	self.bpa_files = header1[3]
	self.beam_files = header1[7]
	self.cells_files = header1[0]

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
	
	self.image1 = realDAT
	self.image2 = realDAT2
		
	self.sigma_cut = float(self.SigmaCut.text())
	self.vmin = self.VMIN.text()
	self.vmax = self.VMAX.text()


	if self.sigma_cut:
		spix_image1 = self.image1*(self.image1 > self.rms1*self.sigma_cut)
		spix_image2 = self.image2*(self.image2 > self.rms2*self.sigma_cut)
		OK = (spix_image1==spix_image1)*(spix_image2==spix_image2) 
		if (spix_image1[OK]*spix_image2[OK]).sum() == 0:
			print("Could not use sigma_cut of %f because it excluded all valid data" % self.sigma_cut)
			spix_image1 = self.image1
			spix_image2 = self.image2
		else:
			print 'Sigma cut was used'

	else:
	        spix_image1 = self.image1
	        spix_image2 = self.image2


	a = np.log10(spix_image1 /spix_image2)/np.log10(self.freq1/self.freq2)

	if self.vmin != '':
		self.vmin = float(self.vmin)
		print self.vmin, self.vmax
		self.aMasked = np.ma.masked_where(a > None, a)
		self.aMasked = np.ma.masked_where(self.aMasked < self.vmin, self.aMasked)	
		print self.aMasked
	else:
		self.vmin = None
		self.aMasked = a
	if self.vmax != '':
		self.vmax = float(self.vmax)
		self.aMasked = np.ma.masked_where(a > self.vmax, a)
		self.aMasked = np.ma.masked_where(self.aMasked < self.vmin, self.aMasked)	
	else:
		self.vmax = None
		self.aMasked = a



	cent_ellips = self.limits_plot(self.ext[1],self.ext[0],self.ext[3],self.ext[2])
	self.x_cent, self.y_cent = cent_ellips[0],cent_ellips[1]
	

	if self.rms1 > self.rms2:
		self.first_contour = self.rms1*1.5*self.sigma_cut
	else:
		self.first_contour = self.rms2*1.5*self.sigma_cut

	self.levels = self.first_contour*np.array([2.,4.,16.,64.,256.,1020.,2050.])

	"""plt.figure(1)
	cset = plt.contour(image1, self.levels, inline=1,
	                  colors=['black'],extent=self.extmas, aspect=1.0)
	plt.figure(1)
	plt.imshow(self.aMasked, origin='bottom',extent=self.extmas, vmin=self.vmin, vmax=self.vmax)
	plt.clim(self.vmin,self.vmax)
	#bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
	#width, height = bbox.width*fig.dpi, bbox.height*fig.dpi

	pts = get_ellipse_coords(a=self.bmaj_files/2, b=self.bmin_files/2, x=self.x_cent ,y=self.y_cent, angle=90+self.bpa_files)
	fill(pts[:,0], pts[:,1], alpha=0.2, facecolor='black', edgecolor='black', linewidth=1, zorder=1)

	plt.axis('scaled')
	plt.xlabel('Right Ascension [mas]')
	plt.ylabel('Relative Declination [mas]')
    	plt.title('Spectral index between %s GHz and %s GHz \n %s' % ('%1.1f' % (self.freq1), '%1.1f' % (self.freq2),needed_param.source_name ))
	plt.xlim(self.ext[0], self.ext[1])
	plt.ylim(self.ext[2],self.ext[3])

	plt.colorbar()

    	#plt.savefig('spectral_index.ps')
    	#plt.savefig('spectral_index.png')

	limits = Annotate()
	plt.show()
		
	ext_new = []
	[self.limplot_x1,self.limplot_x2,self.limplot_y1,self.limplot_y2] = limits()"""

	#[limplot_x1mas,limplot_x2mas,limplot_y1mas,limplot_y2mas] = [(cent_mapx-limplot_x1)*cells,(cent_mapx-limplot_x2)*cells,(limplot_y1-cent_mapy)*cells,(limplot_y2-cent_mapy)*cells]
	
	self.cutMap()

	cent_ellips = self.limits_plot(self.limplot_x2,self.limplot_x1,self.limplot_y1,self.limplot_y2)
	self.x_cent, self.y_cent = cent_ellips[0],cent_ellips[1]


	plt.figure(2)
	ax = plt.gca()
	cset = plt.contour(self.image1, self.levels, inline=1,
	                  colors=['black'],extent=self.extmas, aspect=1.0)
	plt.figure(2)
	pts = get_ellipse_coords(a=self.bmaj_files/2, b=self.bmin_files/2, x=self.x_cent ,y=self.y_cent, angle=90+self.bpa_files)
	fill(pts[:,0], pts[:,1], alpha=0.2, facecolor='black', edgecolor='black', linewidth=1, zorder=1)

	im = plt.imshow(self.aMasked, origin='bottom',extent=self.extmas, vmin=self.vmin, vmax=self.vmax)
	plt.axis('scaled')
	plt.xlabel('Right Ascension [mas]',fontsize=15)
	plt.ylabel('Relative Declination [mas]',fontsize=15)
    	#plt.title('Spectral index between %s %s and %s %s \n %s' % (self.freq1name,self.freq1unit, self.freq2name,self.freq2unit,needed_param.source_name ))
	plt.xlim(self.limplot_x1,self.limplot_x2)
	plt.ylim(self.limplot_y2,self.limplot_y1)

	ax.minorticks_on()
	ax.tick_params('both', length=8, width=2, which='major') 
	plt.tick_params(axis='both',which='both',direction='in', labelsize=13)

	# of ax and the padding between cax and ax will be fixed at 0.05 inch.
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="4.5%", pad="0.5%")

	cb = plt.colorbar(im,cax=cax,cmap='jet')
	cb.ax.tick_params(labelsize=13, width=2)
	cb.set_label(r'$\alpha$',fontsize=15)

    	plt.savefig('spectral_index_between_'+self.freq1name+self.freq1unit+'_'+self.freq2name+self.freq2unit+'_'+needed_param.source_name+'.ps',bbox_inches='tight')
    	plt.savefig('SPIX_png/spectral_index_between_'+self.freq1name+self.freq1unit+'_'+self.freq2name+self.freq2unit+'_'+needed_param.source_name+'.png',bbox_inches='tight')
    	plt.savefig('SPIX_png/spectral_index_between_'+self.freq1name+self.freq1unit+'_'+self.freq2name+self.freq2unit+'_'+needed_param.source_name+'.pdf',bbox_inches='tight')

	plt.close('all')
	
	self.data1 = self.image1.copy()
	self.spix = a.copy()
	self.spixMasked = self.aMasked.copy()
	self.extContours = self.extmas

	os.system('rm SPIX_MAPS/spix'+self.freq1name+self.freq1unit+'-'+self.freq2name+self.freq2unit+'.fits \n')
	hdu = pf.PrimaryHDU(self.spix)
	hdulist = pf.HDUList([hdu])
	hdulist.writeto('SPIX_MAPS/spix'+self.freq1name+self.freq1unit+'-'+self.freq2name+self.freq2unit+'.fits')#,overwrite=True) #CHANGE SAVE NAME

    def limits_plot(self,limplot_x2,limplot_x1,limplot_y1,limplot_y2):
	width = np.abs(limplot_x2 - limplot_x1)
	height = np.abs(limplot_y1 - limplot_y2)

	x_cent = (limplot_x1 - width*0.15)
	y_cent = (limplot_y2 + height*0.15)

	return x_cent,y_cent

    def cutMap(self):
	plt.figure(1)
	cset = plt.contour(self.image1, self.levels, inline=1,
	                  colors=['black'],extent=self.extmas, aspect=1.0)
	plt.figure(1)
	ax = plt.gca()
	im = plt.imshow(self.aMasked, origin='bottom',extent=self.extmas, vmin=self.vmin, vmax=self.vmax)
	plt.clim(self.vmin,self.vmax)
	#bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
	#width, height = bbox.width*fig.dpi, bbox.height*fig.dpi

	pts = get_ellipse_coords(a=self.bmaj_files/2, b=self.bmin_files/2, x=self.x_cent ,y=self.y_cent, angle=90+self.bpa_files)
	fill(pts[:,0], pts[:,1], alpha=0.2, facecolor='black', edgecolor='black', linewidth=1, zorder=1)

	plt.axis('scaled')
	plt.xlabel('Right Ascension [mas]',fontsize=15)
	plt.ylabel('Relative Declination [mas]',fontsize=15)
    	plt.title('Spectral index between %s GHz and %s GHz \n %s' % (self.freq1name, self.freq2name,needed_param.source_name ))
	plt.xlim(self.ext[0], self.ext[1])
	plt.ylim(self.ext[2],self.ext[3])

	ax.minorticks_on()
	ax.tick_params('both', length=8, width=2, which='major') 
	plt.tick_params(axis='both',which='both',direction='in', labelsize=13)

	# of ax and the padding between cax and ax will be fixed at 0.05 inch.
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="4.5%", pad="0.5%")

	cb = plt.colorbar(im,cax=cax,cmap='jet')
	cb.ax.tick_params(labelsize=13, width=2)
	cb.set_label(r'$\alpha$',fontsize=15)

    	#plt.savefig('spectral_index.ps')
    	#plt.savefig('spectral_index.png')

	limits = Annotate()
	plt.show()
		
	ext_new = []
	[self.limplot_x1,self.limplot_x2,self.limplot_y1,self.limplot_y2] = limits()

	newextent = [self.limplot_x1,self.limplot_x2,self.limplot_y1,self.limplot_y2]


	[limplot_x1pix,limplot_x2pix,limplot_y1pix,limplot_y2pix] = [(self.cent_mapx-limits()[0]/self.cells_files),(self.cent_mapx-limits()[1]/self.cells_files),(limits()[2]/self.cells_files+self.cent_mapy),(limits()[3]/self.cells_files+self.cent_mapy)]

	#self.extmas = [limplot_x1pix,limplot_x2pix,limplot_y1pix,limplot_y2pix]
                                     
	#parameters to crop the map (center position of the rectangle and size)
	centx = ((limplot_x1pix+limplot_x2pix)/2.)
	centy = ((limplot_y1pix+limplot_y2pix)/2.) 
	position = (centx,centy)
	
	height = round(np.abs(limplot_y2pix-limplot_y1pix))+1
	width = round(np.abs(limplot_x2pix-limplot_x1pix))+1
	size = (height,width)

	#cropping
	self.cutoutSpix = Cutout2D(self.aMasked,position,size).data
	self.image1cut = Cutout2D(self.image1,position,size).data
	self.image2cut = Cutout2D(self.image2,position,size).data



    def cutSpix(self):

	self.totalLines = int(self.numLines.text())

	#pu levels with self
	self.levels = self.first_contour*np.array([2.,4.,16.,64.,256.,1020.,2050.])

	plt.figure(2)
	ax = plt.gca()
	cset = plt.contour(self.data1, self.levels, inline=1,
	                  colors=['black'], aspect=1.0)
			  #,extent=self.extContours, aspect=1.0)
	im = plt.imshow(self.spixMasked, origin='bottom')#,extent=self.extContours)#, vmin=-2.5, vmax=1.7)
	plt.axis('scaled')
	plt.xlabel('Right Ascension [mas]',fontsize=15)
	plt.ylabel('Relative Declination [mas]',fontsize=15)
    	plt.title('Spectral index between %s %s and %s %s \n %s' % (self.freq1name,self.freq1unit, self.freq2name,self.freq2unit,needed_param.source_name ))
	plt.xlim(-self.limplot_x1/self.cells_files+self.cent_mapx,-self.limplot_x2/self.cells_files+self.cent_mapx)
	plt.ylim(self.limplot_y2/self.cells_files+self.cent_mapy,self.limplot_y1/self.cells_files+self.cent_mapy)

	ax.minorticks_on()
	ax.tick_params('both', length=8, width=2, which='major') 
	plt.tick_params(axis='both',which='both',direction='in', labelsize=13)

	# of ax and the padding between cax and ax will be fixed at 0.05 inch.
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="4.5%", pad="0.5%")

	cb = plt.colorbar(im,cax=cax,cmap='jet')
	cb.ax.tick_params(labelsize=13, width=2)
	cb.set_label(r'$\alpha$',fontsize=15)

	#ld = LineDrawer()
	#ld.draw_line() # here you click on the plot

	#xx = LineDrawer.x
	#yy = LineDrawer.y


	xxMore = []
	yyMore = []

	#if self.totalLines > 1:
	for i in xrange(0,self.totalLines):
		ld = LineDrawer()
		ld.draw_line() # here you click on the plot
		xxMore.append(LineDrawer.x)
		yyMore.append(LineDrawer.y)
	plt.close('all')
	#else:
	#	plt.close('all')


	# of ax and the padding between cax and ax will be fixed at 0.05 inch.

	plt.figure(3)
	ax = plt.gca()
	cset = plt.contour(self.data1, self.levels, inline=1,
	                  colors=['grey'],extent=self.extmas, aspect=1.0)#[-self.extContours[0]/self.cells_files+self.cent_mapx,-self.extContours[1]/self.cells_files+self.cent_mapx,self.extContours[2]/self.cells_files+self.cent_mapy,self.extContours[3]/self.cells_files+self.cent_mapy]
	im = plt.imshow(self.spixMasked, origin='bottom',extent=self.extmas)#[-self.extContours[0]/self.cells_files+self.cent_mapx,-self.extContours[1]/self.cells_files+self.cent_mapx,self.extContours[2]/self.cells_files+self.cent_mapy,self.extContours[3]/self.cells_files+self.cent_mapy])#, vmin=-2.5, vmax=1.7)
	plt.axis('scaled')
	plt.xlabel('Right Ascension [mas]',fontsize=15)
	plt.ylabel('Relative Declination [mas]',fontsize=15)
    	#plt.title('Spectral index between %s %s and %s %s \n %s' % (self.freq1name,self.freq1unit, self.freq2name,self.freq2unit,needed_param.source_name ))
	plt.xlim(self.limplot_x1,self.limplot_x2)
	plt.ylim(self.limplot_y2,self.limplot_y1)
	for i in xrange(0,self.totalLines):
		plt.plot([-(xxMore[i][0]-self.cent_mapx)*self.cells_files,-(xxMore[i][1]-self.cent_mapx)*self.cells_files],[(yyMore[i][0]-self.cent_mapy)*self.cells_files,(yyMore[i][1]-self.cent_mapy)*self.cells_files],color='black',linewidth=3)
	#plt.xlim(-self.limplot_x1/self.cells_files+self.cent_mapx,-self.limplot_x2/self.cells_files+self.cent_mapx)
	#plt.ylim(self.limplot_y2/self.cells_files+self.cent_mapy,self.limplot_y1/self.cells_files+self.cent_mapy)

				#x_data_cent.append(-(float(split[0])-mid_point_x[0])*cell)
				#z_data_cent.append((float(split[1])-mid_point_y[0])*cell)  
	ax.minorticks_on()
	ax.tick_params('both', length=8, width=2, which='major') 
	plt.tick_params(axis='both',which='both',direction='in', labelsize=13)

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="4.5%", pad="0.5%")

	cb = plt.colorbar(im,cax=cax,cmap='jet')
	cb.ax.tick_params(labelsize=13, width=2)
	cb.set_label(r'$\alpha$',fontsize=15)
	#plt.plot([-(xx[0])/self.cells_files+self.cent_mapx,-(xx[1])/self.cells_files+self.cent_mapx],[(yy[0])*self.cells_files+self.cent_mapy,(yy[1])/self.cells_files+self.cent_mapy],color='black')
	#if self.totalLines > 1:

    	plt.savefig('SPIX_png/spectral_index_between_'+self.freq1name+self.freq1unit+'_'+self.freq2name+self.freq2unit+'_'+needed_param.source_name+'LineCut.pdf',bbox_inches='tight')
    	plt.savefig('SPIX_png/spectral_index_between_'+self.freq1name+self.freq1unit+'_'+self.freq2name+self.freq2unit+'_'+needed_param.source_name+'LineCut.png',bbox_inches='tight')


	spixArr = [] #np.array([0.]*self.totalLines)
	rArr = [] #np.array([0.]*self.totalLines)

	for j in xrange(0,self.totalLines):

		xi = int(trunc(xxMore[j][0]))
		xf = int(round(xxMore[j][1]))
		yi = int(trunc(yyMore[j][0]))
		yf = int(round(yyMore[j][1])) #corresponds to xf
		pointsLine = get_line(xi, yi, xf, yf)


		spixArr.append(np.array([0.]*len(pointsLine)))
		rArr.append(np.array([0.]*len(pointsLine)))
		#r = []

		for i in xrange(0,len(pointsLine)):
			x = -(pointsLine[i][0]-self.cent_mapx)*self.cells_files
			y = (pointsLine[i][1]-self.cent_mapy)*self.cells_files

			if np.abs(xxMore[j][1]-xxMore[j][0]) > np.abs(yyMore[j][1]-yyMore[j][0]):
				check = x
			else:
				check = y
		

			if check > 0:
				rArr[j][i] = (-np.sqrt(x**2+y**2))
			else:
				rArr[j][i] = (np.sqrt(x**2+y**2))
		
			spixArr[j][i] = self.spix[pointsLine[i][1],pointsLine[i][0]]	

	spixArr = np.asarray(spixArr) #np.array([0.]*self.totalLines)
	rArr = np.asarray(rArr) #np.array([0.]*self.totalLines)

	plt.figure(4)	
	for j in xrange(0,self.totalLines):
		plt.plot(rArr[j],spixArr[j],'r-')

	plt.xlabel('r [mas]', fontsize=15)
	plt.ylabel(r'$\alpha$', fontsize=17)
	plt.tick_params(axis='both',which='both',direction='in', labelsize=13)
	ax = plt.gca()
	#ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on()
	ax.tick_params('both', length=10, width=2, which='major')
	ax.tick_params('both', length=5, width=1, which='minor')
    	plt.savefig('SPIX_png/spectral_index_between_'+self.freq1name+self.freq1unit+'_'+self.freq2name+self.freq2unit+'_'+needed_param.source_name+'CUT.ps',bbox_inches='tight')
    	plt.savefig('SPIX_png/spectral_index_between_'+self.freq1name+self.freq1unit+'_'+self.freq2name+self.freq2unit+'_'+needed_param.source_name+'CUT.png',bbox_inches='tight')
    	plt.savefig('SPIX_png/spectral_index_between_'+self.freq1name+self.freq1unit+'_'+self.freq2name+self.freq2unit+'_'+needed_param.source_name+'CUT.pdf',bbox_inches='tight')

		
	plt.show()

	res=open(needed_param.path+'/cutspixdata'+self.freq1name+self.freq1unit+'_'+self.freq2name+self.freq2unit+'.p','wb')
	pickle.dump([spixArr,rArr],res)
	res.close()  



    def ErrorSpixMap(self):

		#meany meanx CHANGE FOR REAL VALUES SOON, save pixel value in txt file as well
		mean = [0.2,-5.73]
		#errors y x and z if there is z
		cov = [[1,0],[0,1]]

		self.shiftALL = np.random.multivariate_normal(mean,cov,10**4)

		a = 0
		self.shiftPartly = []
		for i in xrange(0,10**4+250,250):
			self.shiftPartly.append(self.shiftALL[a:i])
			a = i


		len(self.shiftPartly)
	
		self.i2 = 0

		if self.ERRcheck == False:
			map(self.shiftingMCmaps,self.shiftPartly)

			self.up, self.down = pixelerror.errorsPixel(self.len1,self.len2)

		self.errvmin = self.ErrVMIN.text()
		self.errvmax = self.ErrVMAX.text()

		#self.up = np.ma.masked_where(self.up <= 0, self.up)
		#self.down = np.ma.masked_where(self.down <= 0, self.down)

		if self.errvmin != '':
			self.errvmin = float(self.errvmin)
		else:
			self.errvmin = None
		if self.errvmax != '':
			self.errvmax = float(self.errvmax)
		else:
			self.errvmax = None

		
		


		plt.figure(14)
		plt.imshow(self.up, origin='bottom',extent=self.extmas, vmin=self.errvmin, vmax=self.errvmax)
		plt.axis('scaled')
		plt.xlabel('Right Ascension [mas]')
		plt.ylabel('Relative Declination [mas]')
	    	#plt.title('Spectral index between %s GHz and %s GHz \n %s' % ('%1.1f' % (self.freq1), '%1.1f' % (self.freq2),needed_param.source_name ))
		plt.xlim(self.limplot_x1,self.limplot_x2)
		plt.ylim(self.limplot_y2,self.limplot_y1)
		plt.colorbar()
		plt.figure(15)
		plt.imshow(self.down, origin='bottom',extent=self.extmas, vmin=self.errvmin, vmax=self.errvmax)
		plt.axis('scaled')
		plt.xlabel('Right Ascension [mas]')
		plt.ylabel('Relative Declination [mas]')
	    	#plt.title('Spectral index between %s GHz and %s GHz \n %s' % ('%1.1f' % (self.freq1), '%1.1f' % (self.freq2),needed_param.source_name ))
		plt.xlim(self.limplot_x1,self.limplot_x2)
		plt.ylim(self.limplot_y2,self.limplot_y1)
		plt.colorbar()

		self.ERRcheck = True

		plt.show()


    def shiftingMCmaps(self,shiftPartly):

		#pool2 = ThreadPool(4)

		#print self.realDAT
		#for i in xrange(0,len(shiftPartly)):
		#	 pool = ThreadPool(4)
		#print 'i',self.i
		self.aaa= map(self.shift_fft,shiftPartly)


		#pool2.join()
		aaa2 = np.asarray(self.aaa)

		for i in xrange(0,self.len1):
			self.a1 = aaa2[:,i]
			h = h5py.File(str(self.i)+'pixel'+str(i)+'.h5', 'w')
			h.create_dataset('data', data=self.a1)
			h.close()
			del h



		#pool2.close()
		#pool2.join()


		del self.aaa[:]
		del aaa2
 


    def shift_fft(self,shift):

	    #print self.realDAT
	    shift_rows,shift_cols = shift
	    #print shift
	    nr,nc = self.image1cut.shape
	    Nr, Nc = fftfreq(nr), fftfreq(nc)
	    Nc,Nr = np.meshgrid(Nc,Nr)
	    fft_inputarray = np.fft.fft2(self.image1cut)
	    fourier_shift = np.exp(1j*2*np.pi*((shift_rows*Nr)+(shift_cols*Nc)))
	    output_array = np.fft.ifft2(fft_inputarray*fourier_shift)

	    shiftedImage = np.real(output_array)

	    image1 = shiftedImage*(shiftedImage > self.rms1*self.sigma_cut)
	    image2 = self.image2cut*(self.image2cut > self.rms2*self.sigma_cut)

	    #print self.i2

	    a = np.log10(image1/image2)/np.log10(self.freq1/self.freq2)


            #plt.show()

	    print self.i
	    #print 'a', deep_getsizeof(a, set())

	    #self.i2 = self.i2 +1   
	    del image1
            del image2
            del shiftedImage
            del output_array

	    self.i = self.i + 1	

	    self.len1,self.len2 = np.shape(a)
	    a2 = np.reshape(a,(self.len1,1,self.len2))

	    #self.len1 = numero de filas para cada una de las cuales tengo un fichero
				

	    return a2

    def shiftArr(self):
	    shiftsArr = [y,x]
	    return shiftsArr


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

	w = SpixWindow()
	w.show()

	app.exec_()


