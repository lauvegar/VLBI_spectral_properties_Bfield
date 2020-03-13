	if self.synchrotronFit.isChecked():
		while ypix<yf:
			xpix=xstart
			while xpix<xf:            
        			spec=final_data[ypix,xpix,:]        
        			if spec[0] > self.noise[0]*sigma_cut and spec[1] >  self.noise[1]*sigma_cut and spec[2] >  self.noise[2]*sigma_cut: #and spec[3] >  self.noise[3]*sigma_cut:
	
					if guess == 0:
						guess_params = guesses_turnover(self.freqs_conv,spec)
						guess_alpha0 = guess_params[2]
						guess_alphathick = 2.5
						guess_Sm = guess_params[0]
						guess_vm = guess_params[1]
						#taum = 3./2*np.sqrt(1-guess_alpha0/guess_alphathick)-3./2
						guess_v1 = 5.615 #np.exp(np.log(guess_vm)-np.log(taum)/(guess_alpha0-guess_alphathick))
						guess_S1 = 1.543 #synchrotron(guess_v1,guess_Sm, guess_vm, guess_alpha0,guess_alphathick)
		
						parameters = [guess_S1,guess_v1,guess_alpha0,guess_alphathick]

						guess = 1

					plt.ioff()	
        				spec=np.asarray(spec)
               				#plt.figure(j)
               				#plt.plot(self.freqs_conv,spec,'ro')
               				#plt.xscale('log')
		       			#plt.yscale('log')
		       			#plt.savefig('Plot/y_'+str('%1.2f' % (ypix))+'x_'+str('%1.2f' % (xpix))+'.png')
		       			#j=j+1
		       			#plt.close('all')
		       			#tmpres[ypix,xpix,0]=5
		       			print 'y = ', ypix, 'x = ', xpix
		       			print 'spec = ', spec

					plt.close('all')



					if np.max(spec) == spec[0]:
						guess_alpha0 = (spec[len(spec)-1]-spec[0])/(self.freqs_conv[len(self.freqs_conv)-1]-self.freqs_conv[0])
						guess_cte = spec[len(spec)-1]-guess_alpha0*self.freqs_conv[len(self.freqs_conv)-1]
spec=np.asarray(spec)

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

							constant = 5
							nuCutHigh = constant*self.freqs_conv[len(self.freqs_conv)-1]
							ScutHigh = PLValues[0].get('alpha0')*nuCutHigh + PLValues[0].get('cte')
							print nuCutHigh, ScutHigh

							nuCutLow = self.freqs_conv[len(self.freqs_conv)-1]/constant
							ScutLow = ScutHigh*spec[0]/spec[len(spec)-1]
							print nuCutLow, ScutLow

							"""plt.ioff()
							j=j+1
				       			f = plt.figure(j)
							ax = f.add_subplot(111)
							plt.errorbar(self.freqs_conv,spec,yerr=0.1*np.asarray(self.datamax),fmt='ro',ecolor='r', capthick=2)
				       			#plt.plot(self.freqs_conv,spec,'ro')
							xaxis=np.linspace(0,self.freqs_conv[len(self.freqs_conv)-1],1000)
							yaxis = powerLaw(xaxis,PLValues[0].get('cte'),PLValues[0].get('alpha0'))
							plt.plot(xaxis,powerLaw(xaxis,PLValues[0].get('cte'),PLValues[0].get('alpha0')),'g-')
							plt.ylabel(r'$S_y$ [Jy]')
							plt.xlabel(r'$\nu$ [GHz]')
				       			#plt.plot(xaxis,synchrotron(xaxis, Sm, vm, guess_alpha0,guess_alphathick),'g-')
				       			plt.xscale('log',nonposx='clip')
				       			plt.yscale('log',nonposy='clip')
							plt.ylim(10**(-2),10**1)
							plt.xlim(10**(0),10**2)
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

							print PLValues"""

						except RuntimeError:
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

					else:


						chi2Sync = probfit.Chi2Regression(synchrotron, np.asarray(self.freqs_conv), spec, error=None, weights=None)

						try:
							try:
								synchrotronValues = []

								if self.alphaThickFree.isChecked():
									synchrotronFit = iminuit.Minuit(chi2Sync, S_m = guess_Sm, v_m = guess_vm , alpha0 = guess_alpha0, alphathick = guess_alphathick, 
										limit_S_m = (0.,3*np.max(self.dataFits)),limit_v_m=(0.,10**4), limit_alphathick = (0.,3.5))
	
								if self.alphaThick.isChecked() or self.alphaThickCustom.isChecked():
									synchrotronFit = iminuit.Minuit(chi2Sync, S_m = guess_Sm, v_m = guess_vm , alpha0 = guess_alpha0, alphathick = guess_alphathick, 
										limit_S_m = (0.,3*np.max(self.dataFits)),limit_v_m=(0.,10**4), fix_alphathick = True)

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
								plt.errorbar(self.freqs_conv,spec,yerr=0.1*np.asarray(self.datamax),fmt='ro',ecolor='r', capthick=2)
								xaxis=np.linspace(0,self.freqs_conv[len(self.freqs_conv)-1],1000)
								plt.plot(xaxis,synchrotron(xaxis,synchrotronValues[0].get('S_m'), synchrotronValues[0].get('v_m'), 
									 synchrotronValues[0].get('alpha0'), synchrotronValues[0].get('alphathick')),'g-')
								plt.ylabel(r'$S_y$ [Jy]')
								plt.xlabel(r'$\nu$ [GHz]')
					       			#plt.plot(xaxis,synchrotron(xaxis, Sm, vm, guess_alpha0,guess_alphathick),'g-')
					       			plt.xscale('log')
					       			plt.yscale('log')
								plt.ylim(10**(-2),10**1)
								plt.xlim(10**(0),10**2)
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

							except RuntimeError:
								print 'Covariance is not valid. May be the last Hesse call failed?'   

						except RuntimeWarning:
							print 'Covariance is not valid. May be the last Hesse call failed?'  

		    		#else:
		        	#	tmpres[ypix,xpix,0]=100
		    		#print fail
		    		xpix=xpix+1
			ypix=ypix+1

		first_contour = self.dataFits[len(self.dataFits)-1].std()
		levels = first_contour*np.array([2.,4.,16.,64.,128.,256.,512.,1024.,2050.])

		plt.figure(1)
		cset = plt.contour(self.dataFits[len(self.dataFits)-1], levels, inline=1,
			          colors=['grey'],extent=self.ext, aspect=1.0)
		plt.figure(1)
		plt.imshow(self.dataTurnoverS, origin='bottom',extent=self.ext)#, vmin=-2.5, vmax=1.7)
		#bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		#width, height = bbox.width*fig.dpi, bbox.height*fig.dpi

		#pts = get_ellipse_coords(a=self.bmaj_files, b=self.bmin_files, x=x_cent ,y=y_cent, angle=90+self.bpa_files)
		#fill(pts[:,0], pts[:,1], alpha=0.2, facecolor='black', edgecolor='black', linewidth=1, zorder=1)

		plt.axis('scaled')
		plt.xlabel('Right Ascension [mas]')
		plt.ylabel('Relative Declination [mas]')
	    	#plt.title('Spectral index between %s GHz and %s GHz \n %s' % ('%1.1f' % (self.freq1), '%1.1f' % (self.freq2),needed_param.source_name ))
		#plt.xlim(self.ext[0], self.ext[1])
		#plt.ylim(self.ext[2],self.ext[3])

		plt.colorbar()

		plt.figure(2)
		cset = plt.contour(self.dataFits[len(self.dataFits)-1], levels, inline=1,
			          colors=['grey'],extent=self.ext, aspect=1.0)
		plt.figure(2)
		plt.imshow(self.dataTurnoverNu, origin='bottom',extent=self.ext)#, vmin=-2.5, vmax=1.7)
		plt.axis('scaled')
		plt.xlabel('Right Ascension [mas]')
		plt.ylabel('Relative Declination [mas]')
		plt.colorbar()

		plt.show()

	if self.powerLawFit.isChecked():

		while ypix<yf:
			xpix=xstart
			while xpix<xf:            
        			spec=final_data[ypix,xpix,:]        
        			if spec[0] > self.noise[0]*sigma_cut and spec[1] >  self.noise[1]*sigma_cut and spec[2] >  self.noise[2]*sigma_cut and spec[3] >  self.noise[3]*sigma_cut:
	
					if guess == 0:
						#guess_params = guesses_PL(self.freqs_conv,spec)
						guess_alpha0 = (spec[len(spec)-1]-spec[0])/(self.freqs_conv[len(self.freqs_conv)-1]-self.freqs_conv[0])
						guess_cte = spec[len(spec)-1]-guess_alpha0*self.freqs_conv[len(self.freqs_conv)-1]

						guess = 1

					plt.ioff()	
        				spec=np.asarray(spec)
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

						#taum = 3./2*np.sqrt(1-guess_alpha0/guess_alphathick)-3./2
						#vm = np.exp(np.log(taum)/(guess_alpha0-guess_alphathick)+np.log(guess_v1))
						#Sm = synchrotron_v1(vm,guess_S1, guess_v1, guess_alpha0,guess_alphathick)

						plt.ioff()
						j=j+1
			       			f = plt.figure(j)
						ax = f.add_subplot(111)
						plt.errorbar(self.freqs_conv,spec,yerr=0.1*np.asarray(self.datamax),fmt='ro',ecolor='r', capthick=2)
			       			#plt.plot(self.freqs_conv,spec,'ro')
						xaxis=np.linspace(0,self.freqs_conv[len(self.freqs_conv)-1],1000)
						yaxis = powerLaw(xaxis,PLValues[0].get('cte'),PLValues[0].get('alpha0'))
						plt.plot(xaxis,powerLaw(xaxis,PLValues[0].get('cte'),PLValues[0].get('alpha0')),'g-')
						plt.ylabel(r'$S_y$ [Jy]')
						plt.xlabel(r'$\nu$ [GHz]')
			       			#plt.plot(xaxis,synchrotron(xaxis, Sm, vm, guess_alpha0,guess_alphathick),'g-')
			       			plt.xscale('log',nonposx='clip')
			       			plt.yscale('log',nonposy='clip')
						plt.ylim(10**(-2),10**1)
						plt.xlim(10**(0),10**2)
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

					except RuntimeError:
						print 'Covariance is not valid. May be the last Hesse call failed?'   
					#PLFit = iminuit(chi2_PL,cte,alpha0))
		    		else:
		        		tmpres[ypix,xpix,0]=100
		    		print fail
		    		xpix=xpix+1
			ypix=ypix+1
