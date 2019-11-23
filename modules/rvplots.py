import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from modules.helpers import ensurePathExists, truncate, insertSubfolder
from modules.settings import Settings
plt.style.use('ggplot')
plt.rc('font', family='serif')
settings = Settings()

def savePeriodogramPlot(freqs, periods, powers, P, fit_p, subfolder=''):
	plt.figure(figsize=(15,5))
	xvals = freqs
	plt.xlabel(r'frequency, day$^{-1}$')

	if settings.periodogram_units == 'f_log':
		plt.xscale('log')
		plt.xlabel(r'log frequency, day$^{-1}$')

	if settings.periodogram_units == 'p':
		xvals = periods
		plt.xlabel(r'period, days')

	if settings.periodogram_units == 'p_log':
		xvals = periods
		plt.xscale('log')
		plt.xlabel(r'log period, days')


	ls_axis = [min(xvals), max(xvals), min(powers), max(powers)]
	ls_axis[3] += (ls_axis[3] - ls_axis[2]) / 10.0
	plt.axis(ls_axis)
	plt.plot(xvals, powers, 'b-')
	plt.plot((P, P), (0, max(powers)), 'b--', label='original period, P = ' + str(truncate(P)))
	if fit_p > 0:
		plt.plot((fit_p, fit_p), (0, max(powers)), 'r-', label='fitted period, P = ' + str(truncate(fit_p)))
	
	plt.ylabel(r'power')
	plt.legend(loc='best')

	if len(settings.periodogram_plot) > 0:
		ensurePathExists(insertSubfolder(settings.periodogram_plot, subfolder))
		plt.savefig(insertSubfolder(settings.periodogram_plot, subfolder))
	if settings.periodogram_show == True:
		plt.show()
	else:
		plt.clf()

def saveRVPlot(ts, ts_full, rvs, rv_full_fit, rv_full_orig, rv_fit, rv_orig, Ms, SNR, params, params_fit, params_err, subfolder=''):

	# params = Mp, p, a, e, w, v0
	Mplanet = [params[0], params_fit[0], params_err[0]]
	P = [params[1], params_fit[1], params_err[1]]
	a = [params[2], params_fit[2], params_err[2]]
	e = [params[3], params_fit[3], params_err[3]]
	w = [params[4], params_fit[4], params_err[4]]
	v0 = [params[5], params_fit[5], params_err[5]]

	plt.figure(figsize=(15,12))
	plt.suptitle(r'$M_\bigstar = ' + str(truncate(Ms)) + ' M_\odot$, $SNR=' + str(SNR) + '$\n' + \
		'Original orbit: $M_p = ' + str(truncate(Mplanet[0])) + ' M_\oplus$, $P = ' + str(truncate(P[0])) + '$ days, $a = ' + str(truncate(a[0])) + \
		'$ AU, $e = ' + str(truncate(e[0])) + '$, $w = '+ str(truncate(w[0])) + '$ rad, $v_0 = '+ str(truncate(v0[0])) + '$ m/s \n' + \
		'Best fit result: $M_p = ' + str(truncate(Mplanet[1])) + '\pm' + str(truncate(Mplanet[2])) + ' M_\oplus$, $P = ' \
			+ str(truncate(P[1])) + '\pm' + str(truncate(P[2])) + '$ days, $a = ' + str(truncate(a[1])) + '\pm'+ str(truncate(a[2])) + \
		'$ AU, $e = ' + str(truncate(e[1])) + '\pm' + str(truncate(e[2])) + '$, $w = '+ str(truncate(w[1])) + '\pm' + str(truncate(w[2])) + \
			'$ rad, $v_0 = '+ str(truncate(v0[1])) + '\pm'+ str(truncate(v0[2])) + '$ m/s',\
		fontsize=16)
	gs = gridspec.GridSpec(3, 1, height_ratios=[4, 1, 1])
	gs.update(hspace=0.05)
	ax0 = plt.subplot(gs[0])
	ax1 = plt.subplot(gs[1], sharex=ax0)
	ax2 = plt.subplot(gs[2], sharex=ax0)

	ax0.plot(ts_full, rv_full_orig, 'b--', label='original curve')
	ax0.plot(ts, rvs, 'bo', label='datapoints', alpha=0.5)
	ax0.plot(ts_full, rv_full_fit, 'r-', label='fit result')
	ax0.legend(loc=1)

	tmin = min(ts)
	tmax = max(ts)

	# residuals
	ax1.plot((tmin, tmax), (0, 0), 'b--', label='original curve')
	ax1.plot(ts_full, rv_full_fit - rv_full_orig, 'r-', label='fit')
	ax1.plot(ts, rvs - rv_orig, 'bo', label='data points', alpha=0.5)

	ax2.plot((tmin, tmax), (0, 0), 'r-', label='fit')
	ax2.plot(ts, rvs - rv_fit, 'bo', label='data points', alpha=0.5)

	ax0.set_ylabel(r'RV, m/s')
	ax2.set_xlabel(r'time, days')

	ax0.locator_params(axis='y',nbins=10)
	ax1.locator_params(axis='y',nbins=4)
	ax2.locator_params(axis='y',nbins=4)

	for tick in ax0.get_xticklabels():
		tick.set_fontsize(0.0)
	for tick in ax1.get_xticklabels():
		tick.set_fontsize(0.0)
	
	if len(settings.rv_plot) > 0:
		ensurePathExists(insertSubfolder(settings.rv_plot, subfolder))
		plt.savefig(insertSubfolder(settings.rv_plot, subfolder))
	if settings.rv_show:
		plt.show()
	else:
		plt.clf()