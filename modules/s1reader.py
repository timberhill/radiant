import numpy as np
import os.path
from scipy.interpolate import splev, splrep, interp1d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class S1Reader():
	def __init__(self, filename):
		if os.path.isfile(filename):
			self.filename = filename
			self.sfactor = 100 # increase sampling by this value with splines
			self._readFile()
		else:
			raise Exception('Could not find file <' + str(filename) + '>.')

	def _readFile(self):
		self.wls = []
		self.fluxes = []
		self.pbounds = []
		self.centralwls = []

		with open(self.filename) as f:
			nums = f.read().split()

			i = 0
			while i < len(nums):
				if 'Mean' in nums[i]:
					# read block of single line profile
					i += 1 # go to central wavelength
					wl0 = float(nums[i]) * 10.0 # turn to angstroems

					i += 6 # go to number of points
					n = int(nums[i])

					pbound = []
					i += 1 # go to left bound of the phase
					phase_left = float(nums[i])
					i += 1 # go to right bound of the phase
					phase_right = float(nums[i])					
					pbound.append(phase_left)
					pbound.append(phase_right)

					wlbound = []
					i += 1 # go to left wavelength bound
					wlbound.append(float(nums[i]) * 10.0)
					i += 1 # go to right wavelength bound
					wlbound.append(float(nums[i]) * 10.0)

					i += 1 # go to SNR of the LSD profile
					SNR = float(nums[i])

					i += 1 # go to continuum level
					continuum = float(nums[i])

					wlsi = []
					fluxesi = []
					i += 2 # go to the line values

					R_i = n*wl0/(wlbound[1]-wlbound[0])

					if len(self.wls) > 0 and abs(1-R_i/self.R) > 0.01:
						print(' ... WARNING: RESOLUTION MAY BE DIFERENT WITHIN FILE <' + self.filename + '>')

					self.R = R_i

					wl_temp = wlbound[0]
					for j in range(0, n):
						i += 1
						wl_temp += wl_temp / self.R
						wlsi.append(wl_temp)
						fluxesi.append(float(nums[i]))

					self.wls.append(wlsi)
					self.fluxes.append(fluxesi)
					self.pbounds.append(pbound)
					self.centralwls.append(wl0)
					self.continuum = continuum

				i += 1

	def _dopplerShift(self, wls, fluxes, rv, continuum=1.0):
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

	def _correctPhase(self , phase):		
		# correct phase to be [0..1]
		while phase < 0:
			phase += 1
		while phase >= 1:
			phase -= 1

		return phase

	def _getPhaseIndex(self, phase):
		phase = self._correctPhase(phase)
		wl0 = 0
		index = 0
		temp = 1
		for i, p in enumerate(self.pbounds):
			mid = (p[0] + p[1])/2
			if abs(phase - mid) < temp:
				temp = abs(phase - mid)
				index = i

		return index

	def wlToRv(self, phase, wls=None, centralwl=0):
		if wls == None or len(wls) == 0:
			phase_index = self._getPhaseIndex(phase)
			wls = self.wls[phase_index]
		if centralwl <= 0:
			phase_index = self._getPhaseIndex(phase)
			centralwl = self.centralwls[phase_index]

		lsd_rvs = []
		for i in range(0,len(wls)):
			val = 299792.458 * (wls[i] - centralwl)/centralwl
			lsd_rvs.append(val)

		return lsd_rvs

	def getCentralwl(self, phase):
		index = self._getPhaseIndex(phase)
		return self.centralwls[index]

	def shiftLine(self, wls, fluxes, rv=0):
		# interpolate
		data_tck = splrep(wls, fluxes)
		wl_fit = []
		temp = min(wls)
		while temp <= max(wls):
			wl_fit.append(temp)
			temp += (temp / self.R) / self.sfactor
		flux_fit = splev(wl_fit, data_tck)

		# apply RV
		wl, flux = self._dopplerShift(wl_fit, flux_fit, rv, self.continuum)

		# reduce resolution back to original
		f = interp1d(wl, flux, bounds_error=False, fill_value=self.continuum)
		flux_new = f(wls)

		return wls, flux_new


	def getLine(self, phase, rv=0):
		phase_index = self._getPhaseIndex(phase)

		# find profile of the phase
		wl0 = self.wls[phase_index]
		flux0 = self.fluxes[phase_index]

		return self.shiftLine(wl0, flux0, rv)














