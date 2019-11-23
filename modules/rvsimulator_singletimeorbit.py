import numpy as np
from modules.contracts import RVSimulator
from scipy.optimize import fsolve

class SingleTimeOrbitRVSimulator(RVSimulator):

	def __init__(self, timepoints, Mstar, Mplanet, a, e=0, w=0, v0=0, periTime=0):
		self._G = 8.88677e-10 # AU^3 Me^-1 day^-2
		self._SunInEarthMasses = 332978.9

		if periTime != 0:
			for i in range(0,len(timepoints)):
				timepoints[i] -= periTime

		self._timepoints = np.asarray(timepoints)
		self._periTime = periTime # time of perihelion
		self._Mstar = Mstar * self._SunInEarthMasses
		self._Mplanet = Mplanet
		self._a = a
		self._e = e
		self._w = w
		self._v0 = v0
		self._period = np.power(4.0*np.power(np.pi, 2)*np.power(self._a, 3) / (self._G*(self._Mstar + self._Mplanet)), 0.5)
		self._K = self._semi_amplitude()

		self._AUmeters = 1.496e11 # m
		self._DAYseconds = 86400 # sec

	"""
	USAGE SAMPLE:
	sim = SingleTimeOrbitRVSimulator(times, Ms, Mps[planet], a, e)
	ps, rvs = sim.getRVCurve(addNoise=True, sigma=3)


	ARGUMENTS:
	timepoints:	array-like list of datapoints [days]
	Mstar:	mass of star in solar masses (NOTE: the value stores in Earth masses instead)
	a:		semi-major axis in AU
	e:		eccentricity
	w:		argument of periastron
	v0:		systemic velocity
	periTime: time of crossing perihelion
	"""

	def _semi_amplitude(self): # assume sin(i) = 1 or Mp is actually Mp*sin(i)
		Mtot = self._Mplanet + self._Mstar

		t1 = np.power(2*np.pi*self._G / self._period, 1.0/3.0)
		t2 = self._Mplanet / np.power(Mtot, 2.0/3.0)
		t3 = 1.0 / np.power(1 - np.power(self._e, 2), 0.5)

		return t1*t2*t3

	def _timeToPhase(self, t):
		keplerEq = lambda E : E - self._e*np.sin(E) - 2.0*np.pi*t/self._period

		# compute eccentric anomaly from Kepler equation
		E_initial_guess = 2.0*np.pi*t/self._period # as for circular orbit
		E = E_initial_guess
		
		if not self._e == 0.0:
			E = fsolve(keplerEq, E_initial_guess)

		# compute true anomaly
		nu = 2.0 * np.arctan(np.power((1 + self._e) / (1 - self._e), 0.5) * np.tan(E/2.0))

		return nu / (2.0*np.pi)

	def _rv(self, time): # returns m/s
		phase = self._timeToPhase(time)
		return self._K*( np.sin(self._w + 2.0*np.pi*phase) + self._e*np.sin(self._w)) * (self._AUmeters/self._DAYseconds)

	def getRVCurve(self, addNoise=False, sigma=0.5, mean=0.0):
		rvs = self._rv(self._timepoints)

		if addNoise:
			# add Gaussian noise
			noise = np.random.normal(mean, sigma, (len(rvs)))
			return np.asarray(self._timepoints) + self._periTime, np.asarray([x+y for x, y in zip(rvs, noise)]) + self._v0

		return np.asarray(self._timepoints) + self._periTime, np.asarray(rvs) + self._v0