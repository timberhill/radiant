import numpy as np
from scipy.special import wofz
#from PyAstronomy import pyasl
from modules.contracts import SpectrumGenerator

class VoigtLinesGenerator(SpectrumGenerator):

	def __init__(self, centers=[0.0], margin=40, rv=0, g_fwhm=[1.0], l_fwhm=[1e-5], R=1e4, fluxes=[1.0], cont=0):
		self._c = 299792458 # speed of light, m/s

		self._rv = rv
		self._centers = centers
		self._margin = margin
		self._dvelocity = self._c/R
		self._g_fwhm = g_fwhm
		self._l_fwhm = l_fwhm
		self._R = R
		self._fluxes = fluxes
		self._cont = cont


	def _normalize(self, ys):
		maximum = max(ys)
		minimum = min(ys)

		ymin = self._cont
		ymax = max(self._fluxes)

		norm = []
		for i, y in enumerate(ys):
			#print('>>  maximum: ' + str(maximum))
			#print('>>  (y - minimum)/(maximum - minimum): ' + str( (y - minimum)/(maximum - minimum) ))
			norm.append(ymax*(y - minimum)/(maximum - minimum) + ymin)

		return norm

	def _Voigt(self, x, g_fwhm, l_fwhm, semi_amp=1, center=0.0):
	    """
	    Return the Voigt line shape at x with Lorentzian component FWHM l_fwhm
	    and Gaussian component FWHM g_fwhm.

	    """
	    sigma = g_fwhm / np.sqrt(2.0 * np.log(2.0))

	    return semi_amp * np.real(wofz(((x-center) + 1j*l_fwhm)/(sigma*np.sqrt(2.0)))) / (sigma * np.sqrt(2.0*np.pi))

	def _Voigts(self, x, g_fwhms, l_fwhms, semi_amps, centers):
		val = 0.0

		for i in range(0,len(g_fwhms)):
			val += self._Voigt(x, g_fwhms[i], l_fwhms[i], semi_amps[i], centers[i])

		return val

	def getSpectrum(self, SNR=0.0, mean=0.0, absorption=False):
		vals = []
		wls = []
		wl = float(min(self._centers) - self._margin)
				
		while wl <= max(self._centers) + self._margin:
			wls.append(wl)
			shift = float(self._rv)*wl/self._c
			val = self._cont + self._Voigts(wl-shift, self._g_fwhm, self._l_fwhm, self._fluxes, self._centers)

			vals.append(val)
			wl += wl / self._R

		vals_norm = np.asarray(self._normalize(vals))

		result = vals_norm

		if absorption:
			# invert profile
			result = vals_norm - 2.0*(vals_norm - self._cont)

		if SNR > 0.0:
			# add Gaussian noise with sigma = 1.0/SNR
			noise = np.random.normal(mean, 1.0/SNR, (len(result)))
			result = np.asarray([x+y for x, y in zip(result, noise)])

		return np.asarray(wls), np.asarray(result)
		