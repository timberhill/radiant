import numpy as np
import warnings
from modules.contracts import OrbitFitter
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
from lmfit import Model, Parameter, minimize, Parameters, report_fit
from modules.settings import Settings

class RVOrbitFitter(OrbitFitter):

	def _validatePositiveRange(self, r, name, alt=[0.0, float('inf')]):
		if len(r) == 2:
			if r[0] < r[1] or r[0] >= 0:
				return r
			warnings.warn('RVOrbitFitter: invalid <' + name + '> provided. Value should be positive and left boundary strictly less than right. Ignoring...')
				
		return alt


	"""
	m/s, days, M_earth
	starmass in solar masses

	"""
	def __init__(self, tdata, rvdata, errs=None, starmass=1.0, massrange=[], periodrange=[]):
		if len(tdata) != len(rvdata):
			raise ValueError('RVOrbitFitter: <tdata> and <rvdata> should have same dimensions.')

		if errs != None and len(errs) != len(tdata):
			raise ValueError('RVOrbitFitter: <tdata> and <rvdata> should have same dimensions as <errs>.')

		self._G = 8.88677e-10 # AU^3 Me^-1 day^-2
		self._SunInEarthMasses = 332978.9
		self._ts = tdata
		self._rvs = rvdata
		self._errs = errs
		self._starmass = starmass

		self._massrange = self._validatePositiveRange(massrange, 'massrange', alt=[0, 0.1*starmass*self._SunInEarthMasses])
		self._periodrange = self._validatePositiveRange(periodrange, 'periodrange')
		self.settings = Settings()


	def _curve(self, t, Mp, Ms, P, e, w, v0=0): # assume sin(i) = 1 or Mp is actually Mp*sin(i)
		AUmeters = 1.496e11 # m
		DAYseconds = 86400 # sec

		Ms *= self._SunInEarthMasses
		Mtot = Mp + Ms

		a = np.power( self._G*Mtot*P*P / (4*np.pi*np.pi) , 1.0/3.0)
		# P = np.power(4.0*np.power(np.pi, 2)*np.power(a, 3) / (self._G*Mtot), 0.5)

		t1 = np.power(2*np.pi*self._G / P, 1.0/3.0)
		t2 = Mp / np.power(Mtot, 2.0/3.0)
		t3 = 1.0 / np.power(1 - np.power(e, 2), 0.5)

		K = t1*t2*t3

		keplerEq = lambda E : E - e*np.sin(E) - 2.0*np.pi*t/P

		# compute eccentric anomaly from Kepler equation
		E_initial_guess = 2.0*np.pi*t/P # as for circular orbit
		E = E_initial_guess
		
		if not e == 0.0:
			E = fsolve(keplerEq, E_initial_guess)

		# compute true anomaly
		nu = 2.0 * np.arctan(np.power((1 + e) / (1 - e), 0.5) * np.tan(E/2.0))

		return K*(np.sin(w + nu) + e*np.sin(w)) * (AUmeters/DAYseconds) + v0



	def getParameters(self, initial_guess=None):
		# 0: Mplanet [Earth masses]
		# 1: period [days]
		# 2: e (eccentricity)
		# 3: w (argument of periastron) [rad]
		# 4: v0, barycentric velocity [m/s]

		if initial_guess == None:
			initial_guess = [100, (self._periodrange[0] + self._periodrange[1])/2, 0.0, 0.0, 0.0]

		fun = lambda t, _Mp, _p, _e, _w, _v0: self._curve(t, _Mp, self._starmass, _p, _e, _w, v0=_v0)
		bounds = (	[self._massrange[0], self._periodrange[0], 0.0, 0.0, -1e8], \
					[self._massrange[1], self._periodrange[1], self.settings.max_eccentricity, 2.0*np.pi, 1e8])

		popt, pcov = curve_fit(fun, self._ts, self._rvs, p0=initial_guess, bounds=bounds, sigma=self._errs)
	
		errors = [] 
		for i in range(len(popt)):
			try:
				errors.append(np.absolute(pcov[i][i])**0.5)
			except:
				errors.append( 0.00 )

		return  (popt[0], errors[0]), \
				(popt[1], errors[1]), \
				(popt[2], errors[2]), \
				(popt[3], errors[3]), \
				(popt[4], errors[4])