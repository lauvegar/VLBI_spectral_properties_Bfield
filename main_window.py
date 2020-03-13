import threading, time
import warnings
import sys
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import *
from functools import *
import numpy as np
import astropy.io.fits as pf
from functions_conv import order_by_nu, read_conv_params
from functions_align import find_same_beam,beam_array,check_map_params, cuttingMAP,cross_correlation_shifts_FITS, checking_shift
from functions2 import take_header, read_map, saver
from functions2 import convolve_difmap, Annotate
import os,glob, shutil
import subprocess as sub
from astropy.nddata import Cutout2D
from correlate2d import *
#from fast_ftts import *
#from shifting_window import main
from convolution_window_FINAL import ConvWindow
from shifting_window_FINAL_noThreading import ShiftWindow
from shiftingComponents_window import ShiftComponentsWindow
from turnover_window import SpectrumWindow
from spix_window import SpixWindow
from Bfieldpixel_window import Bfieldpixel
from Bfieldcoreshift_window import BfieldcoreshiftWindow
from BfieldTb_window import BfieldTbWindow
from packaging.version import Version, LegacyVersion

if Version(mpl.__version__) > Version('2.0.0'):
	mpl.style.use('classic')
	mpl.rc('image', cmap='jet')


class tooldemo(QMainWindow):
   def __init__(self, parent = None):
      super(tooldemo, self).__init__(parent)
      layout = QVBoxLayout()
      tb = self.addToolBar("File")

      icon = QIcon("images/Beam.png")		
      #icon.iconSize(QSize(24, 24))
      convolution = QAction(icon,"Convolve different images with the same beam",self)
      #tb.setIconSize(QSize,(25,25))
      tb.addAction(convolution)
      convolution.triggered.connect(self.Conv)

      shifting = QAction(QIcon("images/CrossCorrelation.png"),"Shift two different images at two different frequencies using the 2D cross correlation method",self)
      tb.addAction(shifting)
      shifting.triggered.connect(self.Shift)

      shiftingComponent = QAction(QIcon("images/components.png"),"Shift two different images at two different frequencies using an optically thin component",self)
      tb.addAction(shiftingComponent)
      shiftingComponent.triggered.connect(self.ShiftComponent)

      spix = QAction(QIcon("images/Spix.png"),"Creates the spectral index map of two images at two different frequencies ",self)
      tb.addAction(spix)
      spix.triggered.connect(self.Spix)

      tb.addSeparator()

      turnover = QAction(QIcon("images/Turnover.png"),"Calculates de Turnover Frequency and Flux using multifrequency data",self)
      tb.addAction(turnover)
      turnover.triggered.connect(self.Spectrum)

      bfieldcoreshift = QAction(QIcon("images/Bfieldcoreshift.png"),"Calculates the magnetic field using the coreshift",self)
      tb.addAction(bfieldcoreshift)
      bfieldcoreshift.triggered.connect(self.BfieldCoreshift)

      bfieldpixel = QAction(QIcon("images/Bfieldsync.png"),"Calculates the magnetic field using the synchrotron spectrum",self)
      tb.addAction(bfieldpixel)
      bfieldpixel.triggered.connect(self.BfieldPixel)

      bfieldTb = QAction(QIcon("images/BfieldTb.png"),"Calculates the magnetic field using the brightness temperature",self)
      tb.addAction(bfieldTb)
      bfieldTb.triggered.connect(self.BfieldTb)

      tb.setIconSize(QSize(45,45))
		
     # open = QAction(QIcon("open.bmp"),"open",self)
      #tb.addAction(open)
      #save = QAction(QIcon("save.bmp"),"save",self)
      #tb.addAction(save)
      #tb.actionTriggered[QAction].connect(self.toolbtnpressed)
      self.setLayout(layout)
      self.setGeometry(50,50,750,600)
      self.setWindowTitle("Spix+Turnover")

      #tb.iconSize(QSize(24, 24))

      qr = self.frameGeometry()
      cp = QDesktopWidget().availableGeometry().center()
      qr.moveCenter(cp)
      self.move(qr.topLeft())
		
   def toolbtnpressed(self,a):
      print "pressed tool button is",a.text()

   def Conv(self):
	self.w = ConvWindow(self)
	self.setCentralWidget(self.w)
	self.w.show()

   def Shift(self):
	self.w = ShiftWindow(self)
	self.setCentralWidget(self.w)
	self.w.show()

   def ShiftComponent(self):
	self.w = ShiftComponentsWindow(self)
	self.setCentralWidget(self.w)
	self.w.show()

   def Spix(self):
	self.w = SpixWindow(self)
	self.setCentralWidget(self.w)
	self.w.show()

   def Spectrum(self):
	self.w = SpectrumWindow(self)
	self.setCentralWidget(self.w)
	self.w.show()

   def BfieldCoreshift(self):
	self.w = BfieldcoreshiftWindow(self)
	self.setCentralWidget(self.w)
	self.w.show()

   def BfieldPixel(self):
	self.w = Bfieldpixel(self)
	self.setCentralWidget(self.w)
	self.w.show()

   def BfieldTb(self):
	self.w = BfieldTbWindow(self)
	self.setCentralWidget(self.w)
	self.w.show()
		
def main():

	warnings.filterwarnings("ignore")
	"""sourcepath = os.getcwd()

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
	if not os.path.exists('SHIFT_ALL'):
		os.makedirs('SHIFT_ALL')
	if not os.path.exists('Plot_fitted'):
		os.makedirs('Plot_fitted')
	if not os.path.exists('SPIX_MAPS'):
		os.makedirs('SPIX_MAPS')
	if not os.path.exists('SPIX_png'):
		os.makedirs('SPIX_png')

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
       			 shutil.move(os.path.join(sourcepath,files), os.path.join(destinationpath_fits,files))"""


	app = QApplication(sys.argv)
	ex = tooldemo()
	ex.show()	
	sys.exit(app.exec_())
	
if __name__ == '__main__':
   main()
