import numpy as np
from gatspy.periodic import LombScargle as GLS
from modules.contracts import LombScargle

class LombScargleGatspy(LombScargle):
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

	def __init__(self):
		self._valuesSet = False

	def _getYMax (self, xs, ys):
		i = np.nanargmax(ys)
		return xs[i], ys[i]

	def setInput(self, xs, ys, errs, periods, fast=False):
		if len(xs) != len(ys) or len(xs) != len(errs) or len(errs) != len(ys):
			raise ValueError('LombScargleGatspy.setInput: xs, ys and errs should have same dimensions.')

		self._xs = xs
		self._ys = ys
		self._errs = errs
		self._periods = periods
		self._fast = fast

		self._valuesSet = True


	def getPeriodogram(self):
		if not self._valuesSet:
			raise ValueError('LombScargleGatspy.getPeriodogram: provide data first via LombScargleGatspy.setInput() method')

		if not self._fast:
			model = GLSFast(fit_offset=True, fit_period=True)
		else:
			model = GLS(fit_offset=True, fit_period=True)

		model.optimizer.period_range = (min(self._periods), max(self._periods))
		model.fit(self._xs, self._ys, self._errs)
		powers = model.score(self._periods)

		x_max, y_max = self._getYMax(self._periods, powers)
		self.bestPeriod = x_max #model.best_period

		return self._periods, powers

	def getFAP(self):
		if not self._valuesSet:
			raise ValueError('LombScargleGatspy.getPeriodogram: provide data first via LombScargleGatspy.setInput() method')

		raise NotImplementedError('Gatspy package has no FAP estimation')


