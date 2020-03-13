import threading, time
import sys, pickle
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from math import *
from functools import *
import numpy as np
import astropy.io.fits as pf
from functions_align import search_rms
from functions_conv import order_by_nu, read_conv_params, check_map_params
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, Annotate
import os,glob, shutil
import subprocess as sub
from shutil import copyfile


#lineedit = QLineEdit(self)
#validator = QDoubleValidator()
#lineedit.setValidator(QDoubleValidator())
#lineedit.textChanged.connect(self.check_state)
#lineedit.textChanged.emit(lineedit.text())

class ConvWindow(QWidget):
    ''' 
	creates a window used for the user to convolve the maps at different frequencies with the same beam
	To perform the convolution, the code calls difmap.
	The convolution can be done by frequency pairs or with any number of frequency larger than 2.
    '''

    def __init__(self,*args):
        QWidget.__init__(self)

	#sets up the layout of the window. Grid in this case
        layout = QGridLayout()
	#layout.addStretch(1)
	#layout.addLayout(hbox)

	#the number of frequencies will appear dynamically after reading the input fits files
	#for that, a list is defined which will contain the elements that will show this frequencies
	#these elements are, a row with checkboxes and a row with buttons (one per frequency)
        self.checks = []
        self.buttons_freqs = []
	#self.name_button_freq = []

	#given the number of frequencies read in the class needed_files
	#checkboxes are created for each frequency and added to the correspondent list
        for i in xrange(0,len(needed_files.freq)):
	    if needed_files.units[i] == 'MHz':
           	 c = QCheckBox('%s %s' % ('%1.0f' % (needed_files.freq[i]), needed_files.units[i]),self)
	    elif needed_files.units[i] == 'GHz':
           	 c = QCheckBox('%s %s' % ('%1.2f' % (needed_files.freq[i]), needed_files.units[i]),self)
	    c.setFixedSize(100,25)
	    #c.toogle()
            layout.addWidget(c,1,i+1)
            self.checks.append(c)

	#given the number of frequencies read in the class needed_files
	#buttons are created for each frequency and added to the correspondent list
        for i in xrange(0,len(needed_files.freq)):
            c = QPushButton('Beam %s' % ('%1.2f' % (needed_files.freq[i])),self)
            layout.addWidget(c,2,i+1)
            self.buttons_freqs.append(c)

	#self.lineEditBMAJ = QLineEdit(self.gridLayoutWidget)
        #self.lineEditBMAJ.setObjectName(_fromUtf8("lineEditBMAJ"))

	#empty label, created for stetic reasons
	self.labelempty = QLabel()
	self.labelempty2 = QLabel()
	self.labelempty.setFixedSize(25,25)
	
	#label for the BMAJ
	self.labelBMAJ = QLabel()
	#line Edit for BMAJ, it allows the user to put a value in the box
	self.BMAJ = QLineEdit()
	#validates that the given values in the box are correct, (float in this case)
	self.BMAJ.setValidator(QDoubleValidator())
	self.BMAJ.textChanged.connect(self.check_state)
	self.BMAJ.textChanged.emit(self.BMAJ.text())
	self.BMAJ.setFixedSize(100,25)
	#label and lineEdit for BMIN
  	self.labelBMIN = QLabel()
	self.BMIN = QLineEdit()
	self.BMIN.setValidator(QDoubleValidator())#QDoubleValidator(0.#min,3.#max,2.#number decimals))
	self.BMIN.textChanged.connect(self.check_state)
	self.BMIN.textChanged.emit(self.BMIN.text())
	self.BMIN.setFixedSize(100,25)
	#label and lineEdit for BPA
   	self.labelBPA = QLabel()
	self.BPA = QLineEdit()
	self.BPA.setValidator(QDoubleValidator())
	self.BPA.textChanged.connect(self.check_state)
	self.BPA.textChanged.emit(self.BPA.text())
	self.BPA.setFixedSize(100,25)
	#label and lineEdit for the pixel size
	self.labelPIXsize = QLabel()
	self.PIXsize = QLineEdit()
	self.PIXsize.setValidator(QDoubleValidator())
	self.PIXsize.textChanged.connect(self.check_state)
	self.PIXsize.textChanged.emit(self.PIXsize.text())
	self.PIXsize.setFixedSize(100,25)
	#label and lineEdit for the mapsize
	self.labelMAPsize = QLabel()
	self.MAPsize = QLineEdit()
	self.MAPsize.setValidator(QDoubleValidator())
	self.MAPsize.textChanged.connect(self.check_state)
	self.MAPsize.textChanged.emit(self.MAPsize.text())
	self.MAPsize.setFixedSize(100,25)

	#label and lineEdit for the two parameters that difmap uses for the uv weighting
	self.labelUVW = QLabel()
	self.uvw1 = QLineEdit()
	self.uvw2 = QLineEdit()
	self.uvw1.setValidator(QDoubleValidator())
	self.uvw2.setValidator(QDoubleValidator())
	self.uvw1.textChanged.connect(self.check_state)
	self.uvw1.textChanged.emit(self.uvw1.text())
	self.uvw2.textChanged.connect(self.check_state)
	self.uvw2.textChanged.emit(self.uvw2.text())
	self.uvw1.setFixedSize(100,25)
	self.uvw2.setFixedSize(100,25)

	#set the text of the empty label blank
   	self.labelempty.setText(" ")	
	self.labelempty.setFixedSize(10,25)
	#set the text of the labels for BMAJ,BMIN,BPA,cellsize,mapsize and uvw in the interface
   	self.labelBMAJ.setText("BMAJ    : ")
   	self.labelBMIN.setText("BMIN    : ")
   	self.labelBPA.setText("BPA       : ")
   	self.labelPIXsize.setText("Cell    :  ")
   	self.labelMAPsize.setText("Mapsize    : ")
   	self.labelUVW.setText("UVW     : ")

	#define the labels that will contain the information of the beam for each frequency 
	#once the button for the beam of the correspondent frequency is clicked
	self.labelBEAM = QLabel()
	self.labelBEAM2 = QLabel()
	self.labelBEAM3 = QLabel()
	self.labelBEAM4 = QLabel()

	#assigning the previous elements a position in the grid layout
	for i in xrange(0,18):
		layout.addWidget(self.labelempty, i, 0)	
	layout.addWidget(self.labelempty, 4, 1)	
	layout.addWidget(self.labelempty, 5, 1)
	layout.addWidget(self.labelBMAJ, 6, 1)
	layout.addWidget(self.BMAJ, 6, 2, 1, 1)
	layout.addWidget(self.labelBMIN, 7, 1)
	layout.addWidget(self.BMIN, 7, 2, 1, 1)
	layout.addWidget(self.labelBPA, 8, 1)
	layout.addWidget(self.BPA, 8, 2, 1, 1)
	layout.addWidget(self.labelPIXsize, 9, 1)
	layout.addWidget(self.PIXsize, 9, 2, 1, 1)
	layout.addWidget(self.labelMAPsize, 10, 1)
	layout.addWidget(self.MAPsize, 10, 2, 1, 1)
	layout.addWidget(self.labelUVW, 11, 1)
	layout.addWidget(self.uvw1, 11, 2, 1, 1)
	layout.addWidget(self.uvw2, 11, 3, 1, 1)
	layout.addWidget(self.labelBEAM, 13,1, 1, 4)
	layout.addWidget(self.labelBEAM2, 14,1, 1, 4)
	layout.addWidget(self.labelBEAM3, 15,1, 1, 4)
	layout.addWidget(self.labelBEAM4, 16,1, 1, 4)

	#set buddys for the elements that are related among themselves
	#for instance, the label saying BMAJ and the box where the user sets the value
	self.labelBMAJ.setBuddy(self.BMAJ)
	self.labelBMIN.setBuddy(self.BMIN)
	self.labelBPA.setBuddy(self.BPA)     
	self.labelPIXsize.setBuddy(self.PIXsize)
	self.labelMAPsize.setBuddy(self.MAPsize)

	#checkbox for setting circular or elliptical beam
	self.circ = QCheckBox('circular')
	self.ellip = QCheckBox('elliptical')

	#setting the latter checkboxes in the layout
	layout.addWidget(self.circ, 3, trunc(len(needed_files.freq)/2))
	layout.addWidget(self.ellip, 3, trunc(len(needed_files.freq)/2)+1)

	#for the elliptical and circular checkboxes
	#if the user activates one of them, uncheck the other one if activated
	self.ellip.setEnabled(True)
    	#self.circ.toggled.connect(self.ellip.setDisabled)
    	self.circ.toggled.connect(lambda checked: checked and self.ellip.setChecked(False))

	self.circ.setEnabled(True)
	#self.ellip.toggled.connect(self.circ.setDisabled)
	self.ellip.toggled.connect(lambda checked: checked and self.circ.setChecked(False))

	#set the action that the buttons for each frequency do when clicked. It calls the function BEAMparam
   	for i in xrange(0,len(needed_files.freq)):
		self.buttons_freqs[i].clicked.connect(partial(self.BEAMparam, needed_files.bmaj,needed_files.bmin,needed_files.bpa,needed_files.beam,i,self.circ)) 
		#Python only introduces new bindings in namespace through assignment and through parameter lists of functions. i is therefore not actually defined in the namespace of the lambda, but in the namespace of __init__(). The name lookup for i in the lambda consequently ends up in the namespace of __init__(), where i is eventually bound to 9. This is called "closure". You are creating closures. Closures really capture a variable, not the value of a variable. At the end of __init__, i is the last element of range(0, 10), i.e. 9. All the lambdas you created in this scope refer to this i and only when they are invoked, they get the value of i at the time they are at invoked (however, seperate invocations of __init__ create lambdas referring to seperate variables!).

	#create the main two buttons:
	#one for setting the pixel and mapsize
        self.setPIXELandMAP = QPushButton("&Set Pixelsize \n and Mapsize")
	self.setPIXELandMAP.setEnabled(False)
	self.setPIXELandMAP.setAutoDefault(True)
	#self.setPIXELandMAP.setFixedSize(100,25)
	#self.tracker = []
	"""for i in xrange(len(self.checks)):	
	#	self.checks[i].toogle()
		self.checks[i].stateChanged.connect(lambda checked: self.setPIXELandMAP.setEnabled(True))
	#	self.tracker.append(i)"""

	for i in xrange(len(self.checks)):
		self.checks[i].toggled.connect(lambda checked: self.checksState())
	self.setPIXELandMAP.clicked.connect(lambda: self.setPIXandMAPsize(self.checks,needed_files.freqOrig,needed_files.cell,needed_files.size_map))


	#self.progressBar = QProgressBar(self)
        #self.progressBar.setRange(0,100)

	#other button to perform the convolution
        self.convolvebutton = QPushButton("&Convolve")
	self.convolvebutton.setFixedSize(120,25)
        self.convolvebutton.clicked.connect(lambda: self.CONVOLVE(self.checks,needed_files.freq,needed_files.files,needed_files.models,needed_files.source_name,self.convolvebutton))
	self.convolvebutton.setEnabled(False)
	self.convolvebutton.setAutoDefault(True)

	#self.myLongTask = TaskThread()
        #self.myLongTask.notifyProgress.connect(self.onProgress)

	#setting the buttons in the layout
	if len(needed_files.freq) > 2:		
		layout.addWidget(self.convolvebutton, 8, len(needed_files.freq)+1)
		layout.addWidget(self.setPIXELandMAP, 6, len(needed_files.freq)+1,2,1)
		for i in xrange(0,18):
			layout.addWidget(self.labelempty2, i, len(needed_files.freq)+2)	
	else:
		layout.addWidget(self.convolvebutton, 8, len(needed_files.freq)+4)
		layout.addWidget(self.setPIXELandMAP, 6, len(needed_files.freq)+4,2,1)
		#for i in xrange(0,18):
		#	layout.addWidget(self.labelempty2, i, len(needed_files.freq)+4)	

	#layout.addWidget(self.progressBar,17,1,1,len(freq))

	#set the layout of the window
        self.setLayout(layout)

	#for the window to be in the center of the desktop when the program is opened
	qr = self.frameGeometry()
	cp = QDesktopWidget().availableGeometry().center()
	qr.moveCenter(cp)
	self.move(qr.topLeft())

	self.setWindowTitle("Convolution")


    ###############################################################################################
    ###############################################################################################
    ######################                     FUNCTIONS                     ######################    
    ###############################################################################################
    ###############################################################################################

    """
    #function to disable buttons. 
    #In this case, the buttons will get enable if one or more frequencies are checked
    #and disabled if no frequency is checked
    """
    def checksState(self):
	checked = []
	for i in xrange(0,len(self.checks)):
		if self.checks[i].isChecked():
			checked.append(i)
	if len(checked) < 1:
		self.setPIXELandMAP.setEnabled(False)
		self.convolvebutton.setEnabled(False)
	else:
		self.setPIXELandMAP.setEnabled(True)
		self.convolvebutton.setEnabled(True)

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
    #function to write the beam parameters in the window when the button with beam #GHz is pressed
    #also, if the checkbox with circular is checked, writes a circular beam in the text boxes 
    #which corresponds to the one of the button pressed
    #if the checkbox is not pressed, then, the correspondent elliptical beam is written
    """
    def BEAMparam(self,bmaj,bmin,bpa,beam,i,checkboxCIRC):
        self.labelBEAM.setText("BMAJ = %s" % (
            '%1.3f' % (bmaj[i])))
	self.labelBEAM.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
        self.labelBEAM2.setText("BMIN = %s" % (
            '%1.3f' % (bmin[i])))
	self.labelBEAM2.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
        self.labelBEAM3.setText("BPA = %s" % (
            '%1.3f' % (bpa[i])))
	self.labelBEAM3.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
        self.labelBEAM4.setText("beam_circ = %s" % (
            '%1.3f' % (beam[i])))
	self.labelBEAM4.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	if checkboxCIRC.isChecked():
		self.BMAJ.setText('%1.3f' % (beam[i]))
		self.BMIN.setText('%1.3f' % (beam[i]))
		self.BPA.setText('%1.3f' % (0.))
	else:
		self.BMAJ.setText('%1.3f' % (bmaj[i]))
		self.BMIN.setText('%1.3f' % (bmin[i]))
		self.BPA.setText('%1.3f' % (bpa[i]))

    """
    #select the mapsize (mapsize of the low frequency image) and cellsize (half of the cellsize of the high frequency image)
    """
    def setPIXandMAPsize(self,checkBOXes,freq,cell,size_map): #self.checks,freq ###stop_event,checkBOXes,freq
	self.freq_index = []
	for i in xrange(len(checkBOXes)):
		if checkBOXes[i].isChecked():
			self.freq_index.append(i)

	min_freq_index = np.min(self.freq_index)
	max_freq_index = np.max(self.freq_index)
	cell = cell[max_freq_index]
	mapsize = 4*(size_map[min_freq_index])
	self.PIXsize.setText('%1.3f' % cell)
	self.MAPsize.setText('%1.3f' % mapsize)

    """
    #function to convolve
    """
    def CONVOLVE(self,checkBOXes,freq,files,models,source_name,convolvebutton): #self.checks,freq ###stop_event,checkBOXes,freq 
	#self.myLongTask.start()
	self.freq_index = []
	self.files2conv = []
	self.models2conv= []
	self.files_out = []
	self.xshift = 0.
	self.yshift = 0.
	self.beam = np.sqrt(float(self.BMAJ.text())*float(self.BMIN.text()))
	tempFreqUsed = []

	self.uvw1val = int(self.uvw1.text())
	self.uvw2val = int(self.uvw2.text())

	for i in xrange(len(checkBOXes)):
		if checkBOXes[i].isChecked():
			tempFreqUsed.append(i)
			self.freq_index.append(i)
			self.files2conv.append(files[i])
			self.models2conv.append(models[i])
			self.files_out.append(needed_files.path+'/CONV/'+source_name+'_'+str('%1.2f' % (freq[i]))+'GHz_convolved_with_beam'+str('%1.2f' % (self.beam))+'.fits')

	if len(tempFreqUsed) == len(checkBOXes):
		self.files_out = []
		for i in xrange(0,len(tempFreqUsed)):
			self.files_out.append(needed_files.path+'/CONV_ALL/'+source_name+'_'+str('%1.2f' % (freq[i]))+'GHz_convolved.fits')
		if len(checkBOXes) == 2:
			self.files_out = []
			for i in xrange(0,len(tempFreqUsed)):
				self.files_out.append(needed_files.path+'/CONV/'+source_name+'_'+str('%1.2f' % (freq[i]))+'GHz_convolved.fits')
	
	#print len(tempFreqUsed),self.files_out
	#if len(self.freq_index)==0:
	#	convolvebutton.setEnabled(False)
	#else:
	for i in xrange(0,len(self.files2conv)):
		convolve_difmap([self.files2conv[i]],[self.models2conv[i]], float(self.BMAJ.text()), float(self.BMIN.text()), float(self.BPA.text()),self.xshift,self.yshift, float(self.MAPsize.text()),float(self.PIXsize.text()),self.uvw1val,self.uvw2val,[self.files_out[i]])

		rms = search_rms()[0]
		os.system('rm difmap.log*\n')

		if needed_files.freq[self.freq_index[i]] < 0.5:
			self.freqname = str('%1.0f' %(needed_files.freq[self.freq_index[i]]*1000))
			self.frequnit = 'MHz'
		else:
			self.freqname = str('%1.2f' %(needed_files.freq[self.freq_index[i]]))
			self.frequnit = 'GHz'

	
		res=open(needed_files.path+'/rms'+self.freqname+'.p','wb')
		pickle.dump(rms,res)
		res.close()   
	#params2conv = conv_params(freq,cell,bmaj,bmin,bpa,beam,size_map,size_map_y,files,models,fits)


    #def CONVOLVE(self):
    #    	self.stop_event=threading.Event()
    #    	self.c_thread=threading.Thread(target=self.convolveEvent, args=(self.stop_event,))
    #    	self.c_thread.start()   

class TaskThread(QThread):
    notifyProgress = pyqtSignal(int)
    def run(self):
        for i in range(101):
            self.notifyProgress.emit(i)
            time.sleep(0.1)

class needed_files():

	path = os.getcwd()
	sourcepath = os.getcwd()

	if not os.path.exists('UVF'):
		os.makedirs('UVF')
	if not os.path.exists('MODELS'):
		os.makedirs('MODELS')
	if not os.path.exists('FITS'):
		os.makedirs('FITS')
	if not os.path.exists('CONV'):
		os.makedirs('CONV')
	if not os.path.exists('CONV_ALL'):
		os.makedirs('CONV_ALL')
	if not os.path.exists('SHIFT'):
		os.makedirs('SHIFT')
	if not os.path.exists('Shift_parameters'):
		os.makedirs('Shift_parameters')
	if not os.path.exists('SHIFT_ALL'):
		os.makedirs('SHIFT_ALL')
	if not os.path.exists('Plot_fitted'):
		os.makedirs('Plot_fitted')
	if not os.path.exists('SPIX_MAPS'):
		os.makedirs('SPIX_MAPS')
	if not os.path.exists('coreshiftmeas'):
		os.makedirs('coreshiftmeas')

	source = os.listdir(sourcepath)
	destinationpath_uvf = path+'/UVF/'
	destinationpath_models = path+'/MODELS/'
	destinationpath_fits = path+'/FITS/'

	for files in source:
		if files.endswith('.uvf'):
       			 shutil.move(os.path.join(sourcepath,files), os.path.join(destinationpath_uvf,files))
		if files.endswith('.mod'):
       			 shutil.move(os.path.join(sourcepath,files), os.path.join(destinationpath_models,files))
		if files.endswith('.fits'):
       			 shutil.move(os.path.join(sourcepath,files), os.path.join(destinationpath_fits,files))

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
	cell = ordered_params[1]
	bmaj = ordered_params[2]
	bmin = ordered_params[3]
	bpa = ordered_params[4]
	beam = ordered_params[5]
	size_map = ordered_params[6]
	size_map_y = ordered_params[7]
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
	
	xshift = 0
	yshift = 0

def main():
	app = QApplication(sys.argv)

	w = ConvWindow()
	w.show()

	app.exec_()


"""	
   vbox = QVBoxLayout()
   vbox.addWidget(l1)
   vbox.addStretch()
   vbox.addWidget(l2)
   vbox.addStretch()
   vbox.addWidget(l3)
   vbox.addStretch()
   vbox.addWidget(l4)
	
   l1.setOpenExternalLinks(True)
   l4.linkActivated.connect(clicked)
   l2.linkHovered.connect(hovered)
   l1.setTextInteractionFlags(Qt.TextSelectableByMouse)
   win.setLayout(vbox)
	
   win.setWindowTitle("QLabel Demo")
   win.show()
   sys.exit(app.exec_())
	
def hovered():
   print "hovering"
def clicked():
   print "clicked"
"""
