import numpy as np
import scipy as sp
from modules.settings import Settings
from modules.s1reader import S1Reader
from modules.helpers import getClosest
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from modules.helpers import out
import matplotlib.pyplot as plt

class S1SpectraFetcher():
	def __init__(self):
		settings = Settings()

		self.phase_data = False
		self.fit_type = 'spline'

		self.input_files = settings.input_files
		self.timepoints = []	# time of the observations
		self.filenames = []		# names of the S1 files
		self.readers = []		# readers/streams of the S1 files
		self.mperiod = settings.magnetic_period
		self.rperiod = settings.rotation_period

		for i in range(0, len(settings.input_files)):
			x = settings.input_files[i]
			self.timepoints.append(x[0])
			self.filenames.append(x[1])
			self.readers.append(S1Reader(x[1]))

		mphases = np.asarray(self.timepoints)/self.mperiod
		for j in range(0,len(mphases)):
			while mphases[j] >= 1:
				mphases[j] -= 1

		self.mphases, self.phased_readers = (list(t) for t in zip(*sorted(zip(mphases, self.readers))))

		#self.__initialize()

	def __initialize(self):
		# we'll put here 2D array: x[rotPhases_index][tfit_index]
		self.afit = []
		self.mufit = []
		self.sfit = []
		self.dyfit = []

		# how many rotational phases?
		rot_phases = 100
		self.rotPhases = np.linspace(0, 1, rot_phases)

		# TODO : get resolution of the fit in days ?
		dt = 1.0
		self.tfit = np.arange(min(self.timepoints), max(self.timepoints), dt)

		# fit magnetic cycle for each rotational phase of the star
		for irph, rph in enumerate(self.rotPhases):
			out( str(irph) + ' of ' + str(len(self.rotPhases)) + ' . . .' )
			amps = []
			mus = []
			sigmas = []
			dys = []

			# get profiles of the observations and fit them
			for i in range(0,len(self.timepoints)):
				wls, fls = self.readers[i].getLine(rph)
				if i == 0:
					self.wls = wls
					self.R = self.readers[i].R

				p0 = [1-min(fls), sum(wls)/len(wls), 0.2, 1.]
				coeff, var_matrix = curve_fit(self.__gauss, wls, fls, p0=p0)

				amps.append(coeff[0])
				mus.append(coeff[1])
				sigmas.append(coeff[2])
				dys.append(coeff[3])

			# fit profiles' parameters
			tfit, afit   = self.__fit_spline(self.timepoints, amps,   new_x=self.tfit)
			tfit, mufit  = self.__fit_spline(self.timepoints, mus,    new_x=self.tfit)
			tfit, sfit   = self.__fit_spline(self.timepoints, sigmas, new_x=self.tfit)
			tfit, dyfit  = self.__fit_spline(self.timepoints, dys,    new_x=self.tfit)

			# add to self object
			self.afit.append(afit)
			self.mufit.append(mufit)
			self.sfit.append(sfit)
			self.dyfit.append(dyfit)

		out('DONE')

	def __gauss(self, x, *p):
		A, mu, sigma, dy = p
		return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + dy

	def __fit_sine(self, xs, ys, period, new_x=[]):
		def sine(x, *p):
			A, phase, dy = p
			return A * np.sin(x*2*np.pi/period + phase) + dy

		if len(new_x) == 0:
			mn = min(xs)
			mx = max(xs)
			new_x = np.arange(mn, mx, (mx-mn)/100)

		# initial guess
		a0  = (max(ys)-min(ys))/2
		dy0 = (max(ys)+min(ys))/2
		p0  = [a0, 0, dy0]

		# fit
		coeff, var_matrix = curve_fit(sine, xs, ys, p0=p0)

		yfit = []
		for x in new_x:
			yfit.append(sine(x, *coeff))

		return new_x, yfit

	def __fit_spline(self, xs, ys, new_x=[]):
		if len(new_x) == 0:
			mn = min(xs)
			mx = max(xs)
			new_x = np.arange(mn, mx, (mx-mn)/100)

		new_y = interp1d(xs, ys, kind='cubic')(new_x)
		return new_x, new_y

	def __dopplerShift(self, wls, fluxes, rv, continuum=1.0):
		func = interp1d(wls, fluxes)
		fluxes_rv = []
		wls_max = max(wls)
		wls_min = min(wls)
		for wl in wls:
			dwl = wl * rv / 299792458.0
			new_wl = wl - dwl

			if new_wl < wls_min or new_wl > wls_max:
				newFlux = continuum
			else:
				newFlux = func(new_wl)

			fluxes_rv.append(newFlux)

		return wls, fluxes_rv

	def getProfile(self, time, rphase, rv=0):
		if len(self.input_files) == 1:
			self.R = self.readers[0].R
			return self.readers[0].getLine(rphase, rv=rv)

		# TODO : get time from rphase and gaussian parameters (self.*fit)

		# find index of closest point in tfit
		it, t = getClosest(self.tfit, time)
		irph, rph = getClosest(self.rotPhases, rphase)

		# get the new profile
		fluxes_fit = self.__gauss(self.wls, self.afit[irph][it], self.mufit[irph][it], self.sfit[irph][it], self.dyfit[irph][it])

		# apply RV
		self.wls, fluxes_fit_rv = self.__dopplerShift(self.wls, fluxes_fit, rv, continuum=self.dyfit[irph][it])

		# TODO : convert wls to kmps somehow. No idea what central wavelength to use. 

		return self.wls, fluxes_fit_rv


	def getProfile_old(self, rphase=0, mphase=0, rv=0, timepoint=None):
		if len(self.input_files) == 1:
			self.R = self.readers[0].R
			return self.readers[0].getLine(rphase, rv=rv)

		def gauss(x, *p):
			A, mu, sigma, dy = p
			return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + dy

		amps = []
		mus = []
		sigmas = []
		dys = []

		wls = []

		for i in range(0,len(self.timepoints)):
			wls, fls = [], []
			Rsum = 0
			if self.phase_data:
				wls, fls = self.phased_readers[i].getLine(rphase)
				self.R = self.phased_readers[i].R
			else:
				wls, fls = self.readers[i].getLine(rphase)
				self.R = self.readers[i].R

			p0 = [1-min(fls), sum(wls)/len(wls), 0.2, 1.]

			coeff, var_matrix = curve_fit(gauss, wls, fls, p0=p0)
			amps.append(coeff[0])
			mus.append(coeff[1])
			sigmas.append(coeff[2])
			dys.append(coeff[3])


		tfit = [mphase*self.mperiod] # for non-phased data
		xaxis = self.timepoints
		period = self.mperiod
		if self.phase_data:
			tfit = [mphase] # for phased data
			xaxis = self.mphases
			period = 1
			while mphase >= 1:
				mphase -= 1

		if timepoint != None and not self.phase_data:
			tfit = [timepoint]

		afit, mufit, sfit, dyfit = [], [], [], []

		if tfit[0] > max(xaxis):
			out('Got a point above an interpolation range: ' + str(tfit))
			return self.readers[len(self.readers)-1].getLine(rphase)
		elif tfit[0] < min(xaxis):
			out('Got a point above an interpolation range: ' + str(tfit))
			return self.readers[0].getLine(rphase)

		if self.fit_type == 'sine':
			tfit, afit   = self.__fit_sine(xaxis, amps,   period, new_x=tfit)
			tfit, mufit  = self.__fit_sine(xaxis, mus,    period, new_x=tfit)
			tfit, sfit   = self.__fit_sine(xaxis, sigmas, period, new_x=tfit)
			tfit, dyfit  = self.__fit_sine(xaxis, dys,    period, new_x=tfit)
		elif self.fit_type == 'spline':
			tfit, afit   = self.__fit_spline(xaxis, amps,   new_x=tfit)
			tfit, mufit  = self.__fit_spline(xaxis, mus,    new_x=tfit)
			tfit, sfit   = self.__fit_spline(xaxis, sigmas, new_x=tfit)
			tfit, dyfit  = self.__fit_spline(xaxis, dys,    new_x=tfit)
		else:
			raise Exception('Dunno how to fit this way ("' + self.fit_type + '"). I can fit either "sine" or "spline".')
			
		# get the new profile
		fluxes_fit = gauss(wls, afit[0], mufit[0], sfit[0], dyfit[0])

		# plt.step(wls, fluxes_fit, where='mid')
		# plt.show()

		# apply RV
		wls, fluxes_fit_rv = self.__dopplerShift(wls, fluxes_fit, rv, continuum=dyfit[0])

		# TODO : convert wls to kmps somehow. No idea what central wavelength to use. 

		return wls, fluxes_fit_rv
