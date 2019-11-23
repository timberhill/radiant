import numpy as np
from scipy.special import wofz
from modules.contracts import SpectrumGenerator

class VoigtLineGenerator(SpectrumGenerator):

	def __init__(self, center=6562.8, margin=40, rv=0, g_fwhm=1.0, l_fwhm=1e-5, R=1e4, flux=[0, 1]):
		self._c = 299792458 # speed of light, m/s
		self._rv = rv
		self._center = center
		self._span = span
		self._dlambda = center/R
		self._dvelocity = self._c/R
		self._g_fwhm = g_fwhm
		self._l_fwhm = l_fwhm
		self._R = R
		self._flux = flux

	def _normalize(self, ys, ymin=0, ymax=1):
		maximum = max(ys)
		minimum = min(ys)

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

	def getSpectrum(self, addNoise=False, sigma=1, mean=0.0):
		wls0 = np.arange(-self._margin, self._margin, self._dlambda)
		#wls = np.arange(self._center - self._span, self._center + self._span, self._dlambda)

		vals = []
		wls = []
		for i, wl0 in enumerate(wls0):
			wls.append( self._center + wl0)

			shift = self._rv*(self._center + wl0)/self._c
			value = self._Voigt( wl0 - shift , self._g_fwhm, self._l_fwhm)
			vals.append(value)

		vals_norm = self._normalize(vals, ymin=self._flux[0], ymax=self._flux[1])

		if addNoise:
			# add Gaussian noise
			noise = np.random.normal(mean, sigma, (len(vals_norm)))
			return wls, np.asarray([x+y for x, y in zip(vals_norm, noise)])

		return wls, vals_norm


		