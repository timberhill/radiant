from abc import ABCMeta, abstractmethod

class LombScargle(object):
	__metaclass__ = ABCMeta

	@abstractmethod
	def setInput(self): # xs, ys, errs, periods/freqs
		pass

	@abstractmethod
	def getPeriodogram(self):
		pass

	@abstractmethod
	def getFAP(self):
		pass

class RVSimulator(object):
	__metaclass__ = ABCMeta
	
	@abstractmethod
	def getRVCurve(self):
		pass

class SpectrumGenerator(object):
	__metaclass__ = ABCMeta

	@abstractmethod
	def getSpectrum(self):
		pass

class OrbitFitter(object):
	__metaclass__ = ABCMeta

	@abstractmethod
	def getParameters(self):
		pass