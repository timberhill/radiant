import numpy as np
from modules.contracts import LombScargle
from PyAstronomy.pyTiming.pyPeriod import Gls, TimeSeries

"""
documentation:
http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyTimingDoc/pyPeriodDoc/periodograms.html#the-generalized-lomb-periodogram
"""

#class LombScarglePyAs(LombScargle):
class LombScarglePyAs():
	"""
	# BASIC SAMPLE:

	lsg = LombScargleGatspy()
	lsg.setInput(xs, ys, dys, periods)
	periods, powers = lsg.getPeriodogram()

	print('period orig = ' + str(period_base))
	print('period norm = ' + str(lsg.bestPeriod))

	lsg._fast = True
	periods, powers = lsg.getPeriodogram()
	print('period fast = ' + str(lsg.bestPeriod))

	"""

	def __init__(self, xs, ys, errs, freqs=None, ofac=10, hifac=1.0, fast=False, norm='HorneBaliunas', stats=False):
		if len(xs) != len(ys) or len(xs) != len(errs) or len(errs) != len(ys):
			raise ValueError('LombScarglePyAs: xs, ys and errs should have same dimensions.')

		self._timeseries = TimeSeries(xs, ys, errs)

		self._xs = xs
		self._ys = ys
		self._errs = errs
		self._freqs = freqs
		self._fast = fast
		self._ofac = ofac
		self._hifac = hifac
		self._norm = norm
		self._stats = stats

	def _getYMax (self, xs, ys):
		i = np.nanargmax(ys)
		return xs[i], ys[i]


	def _getGlsInstance(self, data):
		if len(self._freqs) == 0:
			return Gls(data, ofac=self._ofac, hifac=self._hifac, norm=self._norm, stats=self._stats, plot=False)
		else:
			return Gls(data, freq=self._freqs, norm=self._norm, stats=self._stats, plot=False)

	def getPeriodogram(self, Pn=0.0):
		data = TimeSeries(self._xs, self._ys, self._errs)
		ls = self._getGlsInstance(data)
		fap = ls.FAP(Pn)

		maxf, maxp = self._getYMax(ls.freq, ls.power)
		self.bestPeriod = 1.0 / maxf

		return ls.freq, ls.power, fap 


	def getFAP(self):
		if not self._valuesSet:
			raise ValueError('LombScarglePyAs.getPeriodogram: provide data first via LombScargleGatspy.setInput() method')

		raise NotImplementedError()


