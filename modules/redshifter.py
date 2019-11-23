from scipy.optimize import curve_fit
#from PyAstronomy import pyasl
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class RedshiftMeasurer():
	def __init__(self):
		pass

	def _polynomial(self, x, p0, p1, p2):
		return p0 + p1*x + p2*x*x


	def _getParabolaVertex(self, x, y):
		denom = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
		A     = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / denom
		B     = (x[2]*x[2] * (y[0] - y[1]) + x[1]*x[1] * (y[2] - y[0]) + x[0]*x[0] * (y[1] - y[2])) / denom
		C     = (x[1] * x[2] * (x[1] - x[2]) * y[0] + x[2] * x[0] * (x[2] - x[0]) * y[1] + x[0] * x[1] * (x[0] - x[1]) * y[2]) / denom

		x_top = -B / (2*A)
		y_top = C - B*B / (4*A)

		return x_top, y_top

	def ccf(self, x, y, maxdelay):
		n = len(y)
		mx = 0
		my = 0
		for i in range(0, n):
			mx += x[i]
			my += y[i]

		mx /= n
		my /= n

		sx = 0
		sy = 0
		for i in range(0, n):
			sx += (x[i] - mx) * (x[i] - mx)
			sy += (y[i] - my) * (y[i] - my)
		
		denom = np.power(sx*sy, 0.5)

		ccf = []
		delay = -maxdelay
		iss = np.arange(0, len(y), 1)
		while delay <= maxdelay:
			sxy = 0
			for i in range(0, n):
				j = i + delay
				if j < maxdelay or j >= n-maxdelay:
				   continue
				else:
				   sxy += (x[i] - mx) * (y[j] - my)
				#/* Or should it be (?)
				#if j < 0 or j >= n:
				#   sxy += (x[i] - mx) * (-my)
				#else:
				#   sxy += (x[i] - mx) * (y[j] - my)
				#*/
			ccf.append(sxy / denom)
			#print('point ' + str(delay) + ' ...')
			#tempiss = np.arange(delay, len(x) + delay, 1)
			#plt.title('cc = ' + str(sxy / denom))
			#plt.step(iss, y, 'r-', where='mid')
			#plt.step(tempiss, x, 'k-', where='mid')
			#plt.show()
			delay += 1

		return np.asarray(ccf)

	def getRV (self, wl, flux, wl_model, flux_model, maxshift=1e3, R=1e5, fit_type='linespline', sfactor=1000):
		if fit_type == 'parabola':
			return self.getRV_parabola(wl, flux, wl_model, flux_model, maxshift, R)

		if fit_type == 'linespline':
			return self.getRV_linespline(wl, flux, wl_model, flux_model, maxshift, R, sfactor)

		raise Exception('Wrong fit type specified.')

	def getRV_linespline (self, wl, flux, wl_model, flux_model, maxshift=1e3, R=1e5, sfactor=1000):
		c = 299792458.0 # speed of light

		wl_fit = []
		flux_fit = []
		wl_model_fit = []
		flux_model_fit = []
		if sfactor > 1:
			data_tck = sp.interpolate.splrep(wl, flux)
			temp = min(wl)
			while temp <= max(wl) :
				wl_fit.append(temp)
				temp += (temp / R) / sfactor
			flux_fit = sp.interpolate.splev(wl_fit, data_tck)

			model_tck = sp.interpolate.splrep(wl_model, flux_model)
			wl_model_fit = wl_fit
			flux_model_fit = sp.interpolate.splev(wl_model_fit, model_tck)
		else:
			wl_fit = wl
			flux_fit = flux
			wl_model_fit = wl_model
			flux_model_fit = flux_model

		maxshift = int(maxshift * R * sfactor / c)
		mcc = self.ccf(flux_model_fit, flux_fit, maxshift)
		xmcc = np.arange(-maxshift, len(mcc)-maxshift) * c / (R * sfactor)

		# this fits parabola to the top 3 points of the CCF and returns peak value
		i = np.argmax(mcc)
		x_top = xmcc[i]
		if i != 0 and i != len(mcc)-1:
			x_top, y_top = self._getParabolaVertex(xmcc[i-1:i+2], mcc[i-1:i+2])


		#from s1reader import S1Reader
		#filename = 'oct12_spots_evenphase.s1'
		#reader = S1Reader(filename)
		#phase = 0
		#shifts = reader.wlToRv(phase, wls=wl)
		#shifts_fit = reader.wlToRv(phase, wls=wl_fit)
		#shifts_model = reader.wlToRv(phase, wls=wl_model)
		#shifts_model_fit = reader.wlToRv(phase, wls=wl_model_fit)

		#plt.figure(figsize=(15,9))
		#gs = gridspec.GridSpec(1, 2)
		#gs.update(hspace=0.05)
		#ax0 = plt.subplot(gs[0])
		#ax1 = plt.subplot(gs[1])
		#ax0.set_xlim([-10, 10])
		#ax0.step(shifts_model,		flux_model, 'k-', where='mid',		label=r'LSD profile template')
		#ax0.step(shifts_model_fit,	flux_model_fit, 'k-', where='mid',)
		#ax0.step(shifts,			flux, 'r-', where='mid', 			label=r'observed LSD profile, $RV=' + str(int(x_top)) + '$ m/s')
		#ax0.step(shifts_fit,		flux_fit, 'r-', where='mid', )
		#ax0.legend(loc=4)
		#ax0.set_xlabel(r'Velocity, $km/s$')
		#ax0.set_ylabel(r'$I/I_c$')
		#ax1.step(xmcc, mcc, 'b-', where='mid', label='CCF from interpolation')
		#ax1.set_xlim([-100, 100])
		#ax1.set_ylim([0.9887, 0.98973])
		#ax1.set_xlabel(r'$RV$, $m/s$')
		#ax1.set_ylabel(r'$CCF$')
		#ax0.legend(loc=4)
		#ax1.legend(loc=4)
		#plt.savefig('sample_ccf.pdf

		# this will just return highest point of the CCF
		return x_top

	def getRV_parabola (self, wl, flux, wl_model, flux_model, maxshift=1e3, R=1e5):
		c = 299792458.0 # spped of light
				
		shift = 200
		mcc = self.ccf(flux_model, flux, shift)
		xmcc = np.arange(-shift, len(mcc)-shift) * c / R
		i = np.argmax(mcc)

		x_top, y_top = self._getParabolaVertex(xmcc[i-1:i+2], mcc[i-1:i+2])

		return x_top