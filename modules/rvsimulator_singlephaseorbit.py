import numpy as np
from modules.contracts import RVSimulator

class SinglePhaseOrbitRVSimulator(RVSimulator):

	def __init__(self):
		self._G = 8.88677e-10 # AU^3 Me^-1 day^-2
		self._SunInEarthMasses = 332978.9
		self._valuesSet = False

	"""
	USAGE SAMPLE:
	sim = SinglePhaseOrbitRVSimulator()
	sim.setInput(phases, Ms, Mps[planet], a, e)
	ps, rvs = sim.getRVCurve(addNoise=True, sigma=3)


	ARGUMENTS:
	phases:	array-like list of datapoints - phases of orbit (span of 1.0 is a full orbit)
	Mstar:	mass of star in solar masses (NOTE: the value stores in Earth masses instead)
	a:		semi-major axis in AU
	e:		eccentricity
	w:		argument of periastron
	v0:		systemic velocity
	"""
	def setInput(self, phases, Mstar, Mplanet, a, e=0, w=0, v0=0): # orbits and such
		self._phases = np.asarray(phases)
		self._Mstar = Mstar * self._SunInEarthMasses
		self._Mplanet = Mplanet
		self._a = a
		self._e = e
		self._w = w
		self._v0 = v0
		self._K = self._semi_amplitude()

		self._AUmeters = 1.496e11 # m
		self._DAYseconds = 86400 # sec

		self._valuesSet = True

	def _semi_amplitude(self): # assume sin(i) = 1 or Mp is actually Mp*sin(i)
		Mtot = self._Mplanet + self._Mstar

		P2 = 4.0*np.power(np.pi, 2)*np.power(self._a, 3) / (self._G*Mtot)
		P = np.power(P2, 0.5)

		t1 = np.power(2*np.pi*self._G / P, 1.0/3.0)
		t2 = self._Mplanet / np.power(Mtot, 2.0/3.0)
		t3 = 1.0 / np.power(1 - np.power(self._e, 2), 0.5)

		return t1*t2*t3

	def _rv(self, phase): # returns m/s
		return self._K*( np.sin(self._w + 2.0*np.pi*phase) + self._e*np.sin(self._w)) * (self._AUmeters/self._DAYseconds)

	def getRVCurve(self, addNoise=False, sigma=0.5, mean=0.0):
		if not self._valuesSet:
			raise ValueError('SinglePhaseOrbitRVSimulator.getRVCurve: provide data first via SinglePhaseOrbitRVSimulator.setInput() method')

		rvs = self._rv(self._phases)

		if addNoise:
			# add Gaussian noise
			noise = np.random.normal(mean, sigma, (len(rvs)))
			return self._phases, np.asarray([x+y for x, y in zip(rvs, noise)]) + self._v0
		else:
			return self._phases, np.asarray(rvs) + self._v0