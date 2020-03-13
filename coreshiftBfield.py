import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import iminuit,probfit
import scipy.special as scp

def rcoreFunct(freqs,A,kr):
	deltar = A*(freqs**(-1./kr)-43.14**(-1./kr))

	return deltar

def chi2_rcoreFunct(freqs,A,kr):
	return np.sum((rcoreFunct(freqs,A,kr)-y)**2)	

def CoreShiftMeasure(coreshift12,v1,v2,kr,z,DL):
	dem1 = (v2**(1./kr)-v1**(1./kr))
	csmeasure = 4.85*10**(-9)*coreshift12*DL*v1**(1./kr)*v2**(1./kr)/((1.+z)**2*dem1)

	return csmeasure

def B1nosimplification(alpha0,gammamin,gammamax,theta,psi,csmeasure,kr,delta,z):

	#1pc = 3*10**(18)cm

	scale = 8.2 #for 0836
	betaapp = 17. #for 0836

	r1 = 3.08*10**(18) #cm
	csmeasure = csmeasure*3.08*10**(18)*10**9 #cmHz

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
	C_alpha = 3**(1.-alpha0)/8.*np.sqrt(np.pi)*scp.gamma((7.-2*alpha0)/4.)*scp.gamma((5.-6.*alpha0)/12.)*scp.gamma((25.-6.*alpha0)/12.)*(scp.gamma((9.-2.*alpha0)/4.))**(-1) #right in christians paper

	numK = (gammamax/gammamin)**(2.*alpha0)-1.
	demK = (gammamax/gammamin)**(2.*alpha0+1.)-1.
	if alpha0 == -0.5:
		K_alphagamma = 0.1
	else:
		K_alphagamma = (2.*alpha0+1.)/(2.*alpha0)*numK/demK


	cte3 = r1*me*c**2/e**2
	var_gammamin = -2*alpha0/(gammamin**(2.*alpha0+1.))

	paren2 = np.pi*C_alpha*psirad/np.sin(thetarad)*K_alphagamma*cte3*var_gammamin*(delta/(1.+z))**(3./2.-alpha0)
	####

	#values of the exponential of the two big parenthesis
	exp1 = (5.-2.*alpha0)/(7.-2.*alpha0) #first parenthesis
	exp2 = -2./(7.-2.*alpha0) #second parenthesis

	B1 = cte*paren1**exp1*paren2**exp2

	print 'B1_proff,',B1
	B12 = cte*(cte2*var1_betaapp**kr)**exp1*(np.pi*C_alpha*cte3*var_gammamin*K_alphagamma*(1.+betaapp)**(-(1.+2*alpha0)/4.)*(1.+z)**(alpha0-3./2.))**exp2


	print 'B1,',B12
	return B1




freqRef = 43.14
freqs = np.asarray([1.67,4.84,22.24])
coreshift = np.asarray([1.76,0.26,0.03])
errCoreshift = np.asarray([0.15,0.07,0.01])

plt.errorbar(freqs,coreshift,yerr=errCoreshift,marker='.',color='r',markersize=8,linestyle='',linewidth=2)
ax = plt.gca()
ax.minorticks_on()
ax.tick_params('both',length=10,width=2,which='major')
ax.tick_params('both',length=5,width=1,which='minor')
ax.set_xticks(ticks=freqs)
ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.ylabel('offset [mas]')
plt.xlabel(r'$\nu$ [GHz]')
plt.ylim(-0.15,2)
plt.xlim(0,30)
plt.savefig('coreshift.png')
#plt.show()

chi2rcore = probfit.Chi2Regression(rcoreFunct, freqs,coreshift , error=errCoreshift, weights=None)

try:

	PLFit = iminuit.Minuit(chi2rcore, A = 3., kr = 1.)

	PLFit.migrad()
	PLFit.hesse()

	print PLFit.values
	print PLFit.errors

	plt.errorbar(freqs,coreshift,yerr=errCoreshift,marker='.',color='r',markersize=8,linestyle='',linewidth=2)
	xaxis=np.linspace(freqs[0],freqs[len(freqs)-1],1000)
	plt.plot(xaxis,rcoreFunct(xaxis,PLFit.values.get('A'), PLFit.values.get('kr')),'r-')
	ax = plt.gca()
	ax.minorticks_on()
	ax.tick_params('both',length=10,width=2,which='major')
	ax.tick_params('both',length=5,width=1,which='minor')
	ax.set_xticks(ticks=freqs)
	ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
	plt.ylabel('offset [mas]')
	plt.xlabel(r'$\nu$ [GHz]')
	plt.ylim(-0.15,2)
	plt.xlim(0,30)
	#plt.show()

except RuntimeError:
	print 'Covariance is not valid. May be the last Hesse call failed?'   

kr = PLFit.values.get('kr')

csMeas = CoreShiftMeasure(coreshift[2],freqs[2],freqRef,PLFit.values.get('kr'),2.17,16945*10**6)
print csMeas

B1 = 0.025*(csMeas**3*(1+2.17)**2/(17.**2*np.deg2rad(0.3)*(np.sin(np.deg2rad(3.2)))**2))**(1./4.)
print B1


#B1nosimplification(alpha0,gammamin,gammamax,theta,psi,csmeasure,kr,delta,z)
B1_2 = B1nosimplification(-0.4,100.,10**5,3.2,2.,csMeas,kr,12.,2.17)
print 'b1, printed', B1_2
print B1_2*(csMeas*(1.+17**2)**(1./2.)*43.**(-1./PLFit.values.get('kr')))


rcore43 = csMeas/(np.sin(np.deg2rad(3)))*43.**(-1./PLFit.values.get('kr'))

print rcore43, 'pc', rcore43/8.2, 'mas'

Bcore = B1_2*rcore43**(-1)

print Bcore





