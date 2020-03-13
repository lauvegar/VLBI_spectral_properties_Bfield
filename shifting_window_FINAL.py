import threading, time
import sys
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import matplotlib.pyplot as plt
from math import *
from functools import *
import numpy as np
import astropy.io.fits as pf
import pickle
from functions_conv import order_by_nu, conv_params, read_conv_params
from functions_align import find_same_beam,beam_array,check_map_params, cuttingMAP,cross_correlation_shifts_FITS, checking_shift
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, Annotate
import os,glob
import subprocess as sub
from astropy.nddata import Cutout2D
from correlate2d import *
from fast_ftts import *


#lineedit = QLineEdit(self)
#validator = QDoubleValidator()
#lineedit.setValidator(QDoubleValidator())
#lineedit.textChanged.connect(self.check_state)
#lineedit.textChanged.emit(lineedit.text())

class updateGUIshiftText(QThread):
	def __init__(self,*args):
	        QThread.__init__(self)

	def __del__(self):
		self.wait()

	"""def _get_top_post(self, subreddit):
		url = "https://www.reddit.com/r/{}.json?limit=1".format(subreddit)
		headers = {'User-Agent': 'nikolak@outlook.com tutorial code'}
		request = urllib2.Request(url, headers=headers)
		response = urllib2.urlopen(request)
		data = json.load(response)
		top_post = data['data']['children'][0]['data']
        return "'{title}' by {author} in {subreddit}".format(**top_post)"""

	def run(self):
		ShiftWindow.labelSHIFTraMAS.__init__.setText("")
		ShiftWindow.labelSHIFTraPIXEL.__init__.setText("")
		ShiftWindow.labelSHIFTdecMAS.__init__.setText("")
		ShiftWindow.labelSHIFTdecPIXEL.__init__.setText("")	

def search_rms():
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

	return rms

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
        	self.noshfitButton = QPushButton("&No Shift")
		
		layout.addWidget(self.shiftButton, 4,1,1,2)
		layout.addWidget(self.noshfitButton, 4, 3,1,2)

		self.setLayout(layout)

		#put the window in the center of the desktop
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

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

	self.position = []
	self.size = []
	self.position_feature = []
	self.size_feature = []

        for i in xrange(0,len(needed_param.freq)):
            c = QCheckBox('%s GHz' % ('%1.2f' % (needed_param.freq[i])),self)
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

	self.labelFITS1.setBuddy(self.labelFITS1file)
	self.labelFITS2.setBuddy(self.labelFITS2file)
	self.labelSHIFTra.setBuddy(self.labelSHIFTraMAS)
	self.labelSHIFTdec.setBuddy(self.labelSHIFTdecMAS)
	self.labelSHIFTraMAS.setBuddy(self.labelSHIFTraPIXEL)
	self.labelSHIFTraMAS.setBuddy(self.labelSHIFTdecPIXEL)


	#self.progressBar = QProgressBar(self)
        #self.progressBar.setRange(0,100)
        findingBEAMbutton = QPushButton("&Check")
	findingBEAMbutton.setFixedSize(100,25)
        findingBEAMbutton.clicked.connect(lambda: self.findingBEAM(self.checks,needed_param.freq))
        shiftingButton = QPushButton("&Shifting")
	shiftingButton.setFixedSize(100,25)
        shiftingButton.clicked.connect(lambda: self.shifting(self.checks,self.fits1,self.fits2,needed_param.files,needed_param.models))


	if len(needed_param.freq) > 2:		
		layout.addWidget(findingBEAMbutton, 6, len(needed_param.freq)+1)
		layout.addWidget(shiftingButton, 7, len(needed_param.freq)+1)
		for i in xrange(0,17):
			layout.addWidget(self.labelempty2, i, len(needed_param.freq)+2)	
	else:
		layout.addWidget(findingBEAMbutton, 6, len(needed_param.freq)+4)
		layout.addWidget(shiftingButton, 7, len(needed_param.freq)+4)
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

    #check checkboxes
    def checkCHECKBOX(self,checkBOXes,setPIXELandMAP): #self.checks,freq ###stop_event,checkBOXes,freq
	self.freq_index = []
	for i in xrange(len(checkBOXes)):
		if checkBOXes[i].isChecked():
			self.freq_index.append(i)
	if len(self.freq_index) == 0:
		setPIXELandMAP.setEnabled(False)
	else:
		setPIXELandMAP.setEnabled(True)

    #select the mapsize (mapsize of the low frequency image) and cellsize (half of the cellsize of the high frequency image)

    def getBEAM(self,RadioButtons,diff_beams,beam1,beam2,freq1_index,freq2_index,fits_conv):

			for i in xrange(len(RadioButtons)):
				if RadioButtons[i].isChecked():
					self.selected_beam = ShiftWindow.diff_beams[i]

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

			see_if_OK = check_map_params(str(self.fits1),str(self.fits2))
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
	ordered_params_conv = order_by_nu(files_temp,files_temp,self.fits_conv)
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
			
		see_if_OK = check_map_params(str(self.fits1),str(self.fits2))
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

	else:
		ShiftWindow.diff_beams = []

		for i in xrange(0,len(self.index_beam12)):
			ShiftWindow.diff_beams.append(self.beam1[self.index_beam12[i][0]])
	
		self.wi = popup_beams()
		self.wi.show()
		self.wi.selectButton.clicked.connect(lambda: self.getBEAM(self.wi.checks2, ShiftWindow.diff_beams,self.beam1,self.beam2,self.freq1_index,self.freq2_index,self.fits_conv))

	return popup_beams


    def shifting(self,checkBOXes,fits1,fits2,files,models):

	self.get_thread = updateGUIshiftText()
	self.get_thread.start()	

	header1 = take_header(fits1)
	header2 = take_header(fits2)
	map_data1 = read_map(fits1)		
	realDAT = map_data1[0]
	map_data2 = read_map(fits2)		
	realDAT2 = map_data2[0]

	#obtaining the beam and cell
	self.bmaj_files = header1[1]
	self.bmin_files = header1[2]
	self.bpa_files = header1[3]
	self.beam_files = header1[7]
	self.cells_files = header1[0]

	#obtaining frequencies
	self.freq1 = header1[5]
	self.freq2 = header2[5]

	#obtaining map centers in pixels
	cent_mapx = map_data1[5]
	cent_mapy = map_data1[6]
	self.mapsize_files = 2*map_data1[7]

	#obtaining the four corners of the maps in mas
	x1 = map_data1[1]
	x2 = map_data1[2]
	y1 = map_data1[3]
	y2 = map_data1[4]

	ShiftWindow.ext=[x1,x2,y1,y2] 

	#cut maps v1 < v2, map region to study
	cut = cuttingMAP(realDAT,realDAT2,cent_mapx,cent_mapy,self.cells_files,self.freq1,self.freq2,iteration=0)

	cutout_v1 = cut[0]
	cutout_v2 = cut[1]
	    
	self.position = cut[2]
	self.size = cut[3]
	ShiftWindow.ext = cut[4]


	#cut maps v1 < v2, feature
	cut = cuttingMAP(cutout_v1.data,cutout_v2.data,cent_mapx,cent_mapy,self.cells_files,self.freq1,self.freq2,iteration=1)

	cutout_v1feature = cut[0]
	cutout_v2feature = cut[1]
    
	self.position_feature = cut[2]
	self.size_feature = cut[3]

	plt.figure(1)
	plt.subplot(121)
	plt.imshow(cutout_v2feature.data, origin='bottom')
	#plt.axis('scaled')
	plt.xlabel('Right Ascension [pixels]')
	plt.ylabel('Relative Declination [pixels]')
	plt.title(str('%1.3f' %(self.freq2))+' GHz')

	plt.subplot(122)
	plt.imshow(cutout_v1feature.data, origin='bottom')
	#plt.axis('scaled')
	plt.xlabel('Right Ascension [pixels]')
	#plt.ylabel('Relative Declination [pixels]')
	plt.title(str('%1.3f' %(self.freq1))+' GHz')

	plt.show()

	#2d cross-correlation

	image1 = cutout_v1feature.data
	image2 = cutout_v2feature.data

	ShiftWindow.offset = cross_correlation_shifts_FITS(image1, image2, sigma_cut=0.004)

	ShiftWindow.xshift=-ShiftWindow.offset[0]*self.cells_files #mas
	ShiftWindow.yshift=ShiftWindow.offset[1]*self.cells_files #mas

	self.shifted_files=[]
	self.shifted_files.append(needed_param.path+'/SHIFT/'+needed_param.source_name+'-'+str(round(self.freq1))+'convolved_with_beam'+str('%1.2f' % (self.beam_files))+'shifted.fits')
	self.shifted_files.append(needed_param.path+'/SHIFT/'+needed_param.source_name+'-'+str(round(self.freq2))+'convolved_with_beam'+str('%1.2f' % (self.beam_files))+'shifted.fits')
	

	self.models_chosen = []
	self.files_chosen = []
	for i in xrange(len(checkBOXes)):
		if checkBOXes[i].isChecked():
			self.models_chosen.append(models[i])
			self.files_chosen.append(files[i])


	self.freqsSHIFTING.setText(("Image2 ( %s GHz)  - image1 ( %s GHz)" % ('%1.2f' % (self.freq2), '%1.2f' % (self.freq1))))
	self.freqsSHIFTING.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTraMAS.setText(("%s mas" % ('%1.4f' % (ShiftWindow.xshift))))
	self.labelSHIFTraMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTraPIXEL.setText(("%s pix" % ('%1.2f' % (ShiftWindow.offset[0]))))
	self.labelSHIFTraPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecMAS.setText(("%s mas" % ('%1.4f' % (ShiftWindow.yshift))))
	self.labelSHIFTdecMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecPIXEL.setText(("%s pix" % ('%1.2f' % (ShiftWindow.offset[1]))))
	self.labelSHIFTdecPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')


	#print 'Image2 (', freq2, ') - image1 (' , freq1, ')' 
	#print 'Shift RA: ', offset[0], 'pixels'
	#print 'Shift DEC: ', offset[1], 'pixels'

	#print 'Shift RA: ', xshift, 'mas'
	#print 'Shift DEC: ', yshift, 'mas'


	os.system('rm difmap.log*\n')

	convolve_difmap([self.files_chosen[0]],[self.models_chosen[0]],self.bmaj_files,self.bmin_files,self.bpa_files,ShiftWindow.xshift,ShiftWindow.yshift,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[0]])

	#rms1 (lower freq)
	self.rms1 = search_rms()[0]
	os.system('rm difmap.log*\n')

	convolve_difmap([self.files_chosen[1]],[self.models_chosen[1]],self.bmaj_files,self.bmin_files,self.bpa_files,0,0,self.mapsize_files,self.cells_files,2,-1,[self.shifted_files[1]])

	self.rms2 = search_rms()[0]
	os.system('rm difmap.log*\n')


	#check if the shift was done properly
	self.checkingShift()

    def checkingShift(self):
	ShiftWindow.offset_new = checking_shift(self.shifted_files,self.position,self.size,self.position_feature,self.size_feature)
		
	ShiftWindow.xshift_new=-ShiftWindow.offset_new[0]*self.cells_files #mas
	ShiftWindow.yshift_new=ShiftWindow.offset_new[1]*self.cells_files #mas
	
	#if the residual shift is large, ask the user if a new shift is needed
	if np.abs(ShiftWindow.xshift_new) > 0.01*np.abs(ShiftWindow.xshift) or np.abs(ShiftWindow.yshift_new) > 0.01*np.abs(ShiftWindow.yshift):
		self.wi = popup_shift()
		self.wi.show()
		self.wi.shiftButton.clicked.connect(lambda: self.shiftAgain())
		self.wi.noshfitButton.clicked.connect(lambda: self.noShiftAgain())

	else:
		self.write_result_file()

    def shiftAgain(self):
	#calculate the new total shift
	final_shiftx = ShiftWindow.xshift+ShiftWindow.xshift_new
	final_shifty = ShiftWindow.yshift+ShiftWindow.yshift_new

	#set the default shift as the old shift
	ShiftWindow.xshift = final_shiftx
	ShiftWindow.yshift = final_shifty

	self.labelSHIFTraMAS.setText(("%s mas" % ('%1.4f' % (ShiftWindow.xshift))))
	self.labelSHIFTraMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTraPIXEL.setText(("%s pix" % ('%1.2f' % (ShiftWindow.offset[0]+ShiftWindow.offset_new[0]))))
	self.labelSHIFTraPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecMAS.setText(("%s mas" % ('%1.4f' % (ShiftWindow.yshift))))
	self.labelSHIFTdecMAS.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')
	self.labelSHIFTdecPIXEL.setText(("%s pix" % ('%1.2f' % (ShiftWindow.offset[1]+ShiftWindow.offset_new[1]))))
	self.labelSHIFTdecPIXEL.setStyleSheet('QLabel {color: blue } QLabel {font: Bold }')

	#shifting
	convolve_difmap([self.files_chosen[0]],[self.models_chosen[0]],self.bmaj_files,self.bmin_files,self.bpa_files,final_shiftx,final_shifty,self.mapsize_files,self.cells_files,2,-1[self.shifted_files[0]])

	self.wi.close()
	self.checkingShift()

    def noShiftAgain(self):

	self.write_result_file()
	self.wi.close()

    def write_result_file(self):
	final = [ShiftWindow.xshift, ShiftWindow.yshift, self.rms1, self.rms2 ,ShiftWindow.ext]
	res=open(needed_param.path+'/Shift_parameters/shift_param'+str(int(round(self.freq1)))+'and'+str(int(round(self.freq2)))+'.p','wb')
	pickle.dump(final,res)
	res.close()     

	final_txt = [[ShiftWindow.xshift, ShiftWindow.yshift,self.rms1, self.rms2,ShiftWindow.ext[0],ShiftWindow.ext[1],ShiftWindow.ext[2],ShiftWindow.ext[3]]]
	header = np.array([['#RAshift DECshift rms(low freq image) rms(high freq image) coordinates of the map window -----> lower freq shifted']])
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
	ordered_params = order_by_nu(files,models,fits)
	
	freq = ordered_params[0]
	files = ordered_params[8]
	models = ordered_params[9]
	fits = ordered_params[10]
	
	#source name
	header = take_header(fits[0])
	source_name = header[8]


def main():
	app = QApplication(sys.argv)

	w = ShiftWindow()
	w.show()

	app.exec_()

main()


"""import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *

def window():
   app = QApplication(sys.argv)
   win = QWidget() 
	
   l1 = QLabel()
   l2 = QLabel()
   l3 = QLabel()
   l4 = QLabel()
	
   l1.setText("Hello World")
   l4.setText("TutorialsPoint")
   l2.setText("welcome to Python GUI Programming")
	
   l1.setAlignment(Qt.AlignCenter)
   l3.setAlignment(Qt.AlignCenter)
   l4.setAlignment(Qt.AlignRight)
   l3.setPixmap(QPixmap("python.jpg"))
	
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
	
if __name__ == '__main__':
   window()"""
