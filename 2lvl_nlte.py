import matplotlib.pyplot as plt 
import numpy as np 
import astropy.io.fits as fits
from numba import jit
import sys

font = {'family' : 'normal',                                   
    	'weight' : 'normal',             
    	'size'   : 18} 

import matplotlib
matplotlib.rc('font', **font)  

#-------------------------------------------
# Short characteristics formal solver. Solves for intensity along the ray using 
# second order polynomial based interpolation between the points. 
# INPUT:
# ------------------------------------------------------------------------------------  
# tau - a 1D array of optical depths 
# S   - a 1D array of source functions
# mu  - scalar cos of the angle w.r.t to atmospheric normal. positive is outgoing rays.
#       negative for incoming ones. 
# I_boundary - scalar, boundary condition (for the given mu)
# ------------------------------------------------------------------------------------
# 
# OUTPUT:
# I - intensity at each point in the atmosphere.
# L - approximate lambda operator at each point in the atmosphere. 
# The function returns these two stacked. 

@jit('float64[:,:](float64[:],float64[:],float64,float64)',nopython=True)
def sc_2nd_order(tau, S, mu, I_boundary):

	ND = S.shape[0]
	
	#first we determine direction, to see how we sweep the grid.
	begin = ND-1
	end = -1
	step = -1
	if (mu<0):
		begin = 0
		end = S.shape[0]
		step = 1

	# allocate space for the output:
	I = np.zeros(ND)
	L = np.zeros(ND)
	I[begin] = I_boundary
	
	# go one point at the time and solve RTE in the integral form on little segments:
	for d in range(begin+step,end-step,step):
		
		delta_u = (tau[d-step] - tau[d])/mu
		delta_d = (tau[d] - tau[d+step])/mu

		expd = np.exp(-delta_u)
		#now some numerical expressions for the contributions:
		if delta_u < 0.01:
			w0=delta_u*(1.-delta_u/2.+delta_u**2/6.-delta_u**3/24.+delta_u**4/120.-delta_u**5/720.+delta_u**6/5040.-delta_u**7/40320.+delta_u**8/362880.)
			w1=delta_u**2*(0.5-delta_u/3.+delta_u**2/8.-delta_u**3/30.+delta_u**4/144.-delta_u**5/840.+delta_u**6/5760.-delta_u**7/45360.+delta_u**8/403200.)
			w2=delta_u**3*(1./3.-delta_u/4.+delta_u**2/10.-delta_u**3/36.+delta_u**4/168.-delta_u**5/960.+delta_u**6/6480.-delta_u**7/50400.+delta_u**8/443520.)
		else:
			w0 = 1.0 - expd
			w1 = w0 - delta_u * expd
			w2 = 2.0 * w1 - delta_u * delta_u * expd

		# coefficients that tells us how local, previous and the next source function
		# contribute to the intensity
		psi0 = w0 + (w1 * (delta_u/delta_d - delta_d / delta_u) - w2 * (1.0 / delta_d + 1.0 / delta_u)) / (delta_u + delta_d)
		psiu = (w2 / delta_u + w1*delta_d/delta_u)/(delta_u+delta_d)
		psid = (w2 / delta_d - w1 * delta_u/delta_d)/(delta_u+delta_d)

		I[d] = I[d-step]*expd + psiu*S[d-step] + psi0*S[d] + psid*S[d+step]
		L[d] = psi0
	
	#last point is linear:
	d = end-step
	delta_u = (tau[d-step] - tau[d])/mu
	expd = np.exp(-delta_u)	
	psi0 = 1.0 - 1.0/delta_u * (1.0 - expd)
	psiu = -expd + 1.0/delta_u * (1.0 - expd)

	if (delta_u < 0.01):
		expd = 1.0 - delta_u + delta_u**2.0 / 2.0 - delta_u**3.0/6.0
		psi0 = delta_u/2. - delta_u*delta_u/6. + delta_u**3.0 / 24.
		psiu = delta_u/2. - delta_u**2.0/3. + delta_u**3.0 / 8.

	I[d] = I[d-step]*expd + psiu*S[d-step] + psi0*S[d]
	L[d] = psi0

	return np.stack((I,L))

@jit('float64[:,:](float64[:],float64[:],float64[:],float64[:],float64[:])',nopython=True)
def calc_lambda_full(tau, mu, wmu, profile, wx):

	ND = tau.shape[0]
	Lambda_full = np.zeros((ND,ND))

	for d in range(0,ND):
		S_mock = np.zeros(ND)
		S_mock[d] = 1.0
		for m in range(0,mu.shape[0]):
			for l in range(0,profile.shape[0]):

				#outward:
				Lambda_monoc = sc_2nd_order(tau*profile[l],S_mock,mu[m],0.)

				Lambda_full[:,d]+=Lambda_monoc[0]*profile[l]*wx[l]*wmu[m]*0.5
				
				#inward
				Lambda_monoc = sc_2nd_order(tau*profile[l],S_mock,-mu[m],0.)

				Lambda_full[:,d]+=Lambda_monoc[0]*profile[l]*wx[l]*wmu[m]*0.5

	return Lambda_full

@jit('float64[:,:](float64[:],float64)',nopython=True)
def calc_lambda_monoc(tau, mu):
	ND = tau.shape[0]
	lambda_monoc = np.zeros((ND,ND))

	for d in range(0,ND):
		S_mock = np.zeros(ND)
		S_mock[d] = 1.0
			
		lambda_monoc[:,d] = sc_2nd_order(tau,S_mock,mu,0.)[0]

	return lambda_monoc

@jit('float64[:](float64[:],float64[:],float64[:],float64[:],float64[:],float64[:])')
def calc_anisotropy(tau,mu,wmu,profile,wx,S):

	ND = tau.shape[0]
	J_02 = np.zeros(ND)

	for m in range(0,mu.shape[0]):
		for l in range(0,profile.shape[0]):

			#outward:
			I = sc_2nd_order(tau*profile[l],S,mu[m],S[-1])
			J_02+=I[0]*profile[l]*wx[l]*wmu[m]*0.5*(1.0-3.*mu[m]*mu[m]) /2.4248
			
			#inward
			I = sc_2nd_order(tau*profile[l],S,-mu[m],0.)
			J_02+=I[0]*profile[l]*wx[l]*wmu[m]*0.5*(1.0-3.*mu[m]*mu[m]) /2.4248
	
	return J_02



#@jit('float64[:,:](float64[:],float64[:],float64[:],float64[:],float64[:],float64[:,:])')
#def calc_anisotropy_response(tau,mu,wmu,profile,wx,dS_dB):




@jit('float64[:](float64[:],float64[:],float64,float64[:],float64)',nopython=True)
def one_full_fs(tau,S,mu,profile,boundary):

	NL = profile.shape[0]
	I = np.zeros(NL)

	for l in range(0,NL):
		I[l] = sc_2nd_order(tau*profile[l],S,mu,boundary)[0,0]

	return I

#-------------------------------------------

# We will define a discrete grid for our spacial coordinate, that is logtau 

ND = 91
logtau = np.linspace(-7,2,ND)
tau = 10.0**logtau

B = np.zeros(ND)
B[:] = 1.0

eps = np.zeros(ND)
eps[:] = 1E-4

line_ratio = 1E3

NL = 21
x = np.linspace(-4,4,NL)
profile = 1./np.sqrt(np.pi) * np.exp(-(x**2.0))
prof_norm = np.sum(profile)
wx = np.zeros(NL)
wx[0] = (x[1] - x[0]) * 0.5
wx[-1] = (x[-1] - x[-2]) * 0.5
wx[1:-1] = (x[2:NL] - x[0:-2]) * 0.5
norm =  (np.sum(profile*wx))
wx/= norm

S = np.copy(B)

mu=([1./np.sqrt(3.0)])
wmu=[1.0]

mu=np.cos([0.4793425352,1.0471975512,1.4578547042])
wmu=[.2777777778,0.4444444444,0.2777777778]
NM = mu.shape[0]

mu = np.asarray(mu)
wmu = np.asarray(wmu)

I_boundary_lower = B[-1];
I_boundary_upper = 0.0;


for iter in range(0,200):
	
	J = np.zeros(ND)
	L = np.zeros(ND)

	for m in range(0,NM):
		for l in range(0,NL):

			#outward
			ILambda = sc_2nd_order(tau*profile[l]*line_ratio,S,mu[m],B[-1])

			J+=ILambda[0]*profile[l]*wx[l]*wmu[m]*0.5
			L+=ILambda[1]*profile[l]*wx[l]*wmu[m]*0.5

			#inward
			ILambda = sc_2nd_order(tau*profile[l]*line_ratio,S,-mu[m],0)

			J+=ILambda[0]*profile[l]*wx[l]*wmu[m]*0.5
			L+=ILambda[1]*profile[l]*wx[l]*wmu[m]*0.5
			
	dS = (eps * B + (1.-eps) * J - S) / (1.-(1.-eps)*L)

	max_change  = np.max(np.abs(dS / S))
	print (max_change)

	S += dS
	if (max_change<1E-4):
		break;

LL = calc_lambda_full(tau*line_ratio,mu,wmu,profile,wx)



# Response function of the source function to the temperature:
dS_dB = np.linalg.inv(np.eye(ND)-(1.-eps)*LL)@ (eps*np.eye(ND))

# Then the emergent intensity is:
x_detailed = np.linspace(-6,6,121)
detailed_profile = 1./np.sqrt(np.pi) * np.exp(-x_detailed**2.0)

spectra = one_full_fs(tau*line_ratio,S,1.0,detailed_profile,B[-1])

rf = np.zeros((ND,detailed_profile.shape[0]))

for d in range(0,ND):
	rf[d] = one_full_fs(tau*line_ratio,S+dS_dB[:,d]*1E-3,1.0,detailed_profile,B[-1])-spectra

CF = np.zeros((ND,detailed_profile.shape[0]))

for l in range(0,detailed_profile.shape[0]):
	CF[:,l] = S * np.exp(-tau*line_ratio*detailed_profile[l] * tau*line_ratio*detailed_profile[l])

#plot the results.
plt.figure(figsize=[14,5])
plt.clf()
plt.cla()
plt.subplot(121)
plt.pcolormesh(x_detailed,np.log10(tau),rf/spectra)
plt.xlabel("Reduced wavelength")
plt.ylabel("$\\mathcal{R}(\\tau,\lambda)/I_\lambda$")
plt.title("Response function")
plt.subplot(122)
plt.pcolormesh(x_detailed,np.log10(tau),CF/spectra)
plt.xlabel("Reduced wavelength")
plt.ylabel("$\\mathcal{C}(\\tau,\lambda)/I_\lambda$")
plt.title("Contribution function")
plt.tight_layout()
plt.savefig("2LVL_NLTE_RF_vs_CF.png",fmt='png',bbox_inches='tight')

#plot more results 
plt.clf()
plt.cla()
plt.subplot(121)
plt.plot(np.log10(tau),np.log10(S),label='Source function')
plt.plot(np.log10(tau),np.log10(B),label='Planck function')
plt.legend()
plt.xlabel("$\\log\\tau_c$")
plt.subplot(122)
plt.plot(x_detailed,spectra)
plt.xlabel("Reduced wavelength")
plt.ylabel("Intensity")
plt.tight_layout()
plt.savefig("2LVL_NLTE.png",fmt='png',bbox_inches='tight')

# ==========================================================================================
# Anisotropy

J_02 = calc_anisotropy(tau*line_ratio, mu, wmu, profile, wx, S)








