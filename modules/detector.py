import numpy as np
from scipy.optimize import curve_fit
from modules.lombscargle_pyastronomy import LombScarglePyAs
from modules.rvorbitfitter import RVOrbitFitter
from modules.settings import Settings
from modules.helpers import out, ensurePathExists, savedatafile, truncate, insertSubfolder
from modules.rvplots import savePeriodogramPlot, saveRVPlot
settings = Settings()


def runDetector(ts, rvs, ts0, rvs0, errs, Ms, Mps, a, P, e, w, v0, SNR, subfolder='', fit_rotation=False):
	G = 8.88677e-10 # AU^3 Mearth^-1 day^-2

	if fit_rotation:
		out( 'subtracting stellar rotation ...' )
		def sine(t, period, amp, phase):
			return amp*np.sin(2*np.pi*t/period + phase)

		for x in range(1,4): # fit P, P/2, P/3		
			rotation_period = settings.rotation_period / x
			fit_curve = lambda _t, _amp, _phase: sine(_t, settings.rotation_period, _amp, _phase)
			par, covm = curve_fit(fit_curve, ts, rvs, p0=[max(rvs), 0])
			rv_rot_fit = sine(ts, settings.rotation_period, par[0], par[1])
			rvs = rvs - rv_rot_fit


	out( 'calculating periodogram ...' )
	lsg = LombScarglePyAs(ts, rvs, errs, freqs=np.arange(settings.periodogram_low, settings.periodogram_high, settings.periodogram_step))
	freqs, powers, fap = lsg.getPeriodogram()
	periods = 1.0 / freqs

	# save LSP data
	if len(insertSubfolder( settings.periodogram_data , subfolder)) > 0:
		savedatafile(insertSubfolder( settings.periodogram_data , subfolder), \
			(freqs, periods, powers), \
			('freq', 'period', 'power'))


	out( 'fitting orbit ...' )
	fit_success = True
	try:
		# get nyquist estimate
		mindiff = 1e10 # minimum time difference between points
		for i in range(0, len(ts)-1):
			val = ts[i+1] - ts[i]
			if val < mindiff:
				mindiff = val

		fitter = RVOrbitFitter(ts, rvs, starmass=Ms, periodrange=[2*mindiff, 7000])#periodrange=[0.5*lsg.bestPeriod, 2.0*lsg.bestPeriod])
		fit_Mp, fit_p, fit_e, fit_w, fit_v0 = fitter.getParameters(initial_guess=[Mps, P, e, w, v0])

		def afromp(p, m):
			return np.power( G*(Ms*332978.9 + m)*p*p / (4*np.pi*np.pi) , 1.0/3.0)

		fit_a = (afromp(fit_p[0], fit_Mp[0]), afromp(fit_p[1], fit_Mp[1]))
	except Exception as exc:
		fit_success = False
		print(exc)

	tmin = min(ts)
	tmax = max(ts)
	tmargin = 0 # np.abs(tmax - tmin) / 100
	tmin -= tmargin
	tmax += tmargin
	ts_full = np.arange(tmin, tmax, (tmax-tmin)/settings.rv_curves_points)

	if not fit_success:
		out( 'WARNING : sorry, Master, I could not fit the orbit T_T' )

		# save rv data
		if len(insertSubfolder(settings.rv_data, subfolder)) > 0:
			savedatafile(insertSubfolder(settings.rv_data, subfolder), \
				(ts, rvs, rvs-fitter._curve(ts, Mps, Ms, P, e, w, v0)), \
				('time[days]', 'RV[m/s]', 'original_resid'))

		if len(insertSubfolder( settings.rv_curves , subfolder)) > 0:
			savedatafile(insertSubfolder( settings.rv_curves , subfolder), (ts_full, fitter._curve(ts_full, Mps, Ms, P, e, w, v0)), \
				('time[days]', 'original_RV[m/s]'))

		rv_full_orig = fitter._curve(ts_full, Mps, Ms, P, e, w, v0)
		rv_orig = fitter._curve(ts, Mps, Ms, P, e, w, v0)
		rv_full_fit = np.linspace(0, 0, len(ts_full))
		rv_fit = np.linspace(0, 0, len(ts))

		# params = Mp, p, a, e, w, v0
		out( 'producing plots ...' )
		saveRVPlot(ts, ts_full, rvs, rv_full_fit, rv_full_orig, rv_fit, rv_orig, Ms, SNR, \
			(Mps, P, a, e, w, v0), \
			(0,0,0,0,0,0), \
			(0,0,0,0,0,0), \
			subfolder=subfolder)

		# LS periodigram
		if len(insertSubfolder( settings.periodogram_plot , subfolder)) > 0:
			savePeriodogramPlot(freqs, periods, powers, P, -1, subfolder=subfolder)

		return


	out('[OUTPUT] Mp = '   +   str(truncate(fit_Mp[0])) + u'\t+-\t' + str(truncate(fit_Mp[1])) + '\t(' + str(truncate(100.0*fit_Mp[1]/fit_Mp[0])) + ' %)  [earth masses]')
	out('[OUTPUT] a =  '   +   str(truncate(fit_a[0]) ) + u'\t+-\t' + str(truncate(fit_a[1]) ) + '\t(' + str(truncate(100.0*fit_a[1]/fit_a[0])  ) + ' %)  [AU]')
	out('[OUTPUT] e =  '   +   str(truncate(fit_e[0]) ) + u'\t+-\t' + str(truncate(fit_e[1]) ) + '\t(' + str(truncate(100.0*fit_e[1]/fit_e[0])  ) + ' %)')
	out('[OUTPUT] w =  '   +   str(truncate(fit_w[0]) ) + u'\t+-\t' + str(truncate(fit_w[1]) ) + '\t(' + str(truncate(100.0*fit_w[1]/fit_w[0])  ) + ' %)  [radians]')
	out('[OUTPUT] v0 = '   +   str(truncate(fit_v0[0])) + u'\t+-\t' + str(truncate(fit_v0[1])) + '\t(' + str(truncate(100.0*fit_v0[1]/fit_v0[0])) + ' %)  [m/s]')
	out('[OUTPUT] P =  '   +   str(truncate(fit_p[0]) ) + u'\t+-\t' + str(truncate(fit_p[1]) ) + '\t(' + str(truncate(100.0*fit_p[1]/fit_p[0])  ) + ' %)  [days]')
	out('[OUTPUT] LS P = ' + str(truncate(lsg.bestPeriod)) + ' [days]\n')

	# sigma is assumed to be 1
	err_sum = np.sum((fitter._curve(ts, fit_Mp[0], Ms, fit_p[0], fit_e[0], fit_w[0], fit_v0[0]) - rvs)**2)
	chi2reduced = err_sum / (len(rvs) - 5)

	out( 'chi2 = ' + str(chi2reduced))
	residuals = rvs - fitter._curve(ts, Mps, Ms, P, e, w, v0)
	out( 'RMS  = ' + str(np.sqrt(np.mean(residuals**2))))


	# save parameters
	if len(insertSubfolder(settings.outputfile , subfolder)) > 0:
		ensurePathExists(insertSubfolder(settings.outputfile , subfolder))
		with open(insertSubfolder(settings.outputfile , subfolder), 'w+') as f:
			f.write('Parameter\tunits\t\toriginal\tdetermined\terrorbar')
			f.write('\n\nMplanet\t\t'  + '[Mearth]\t'  + str(Mps)         + '\t\t'   + str(truncate(fit_Mp[0])) + '\t\t' + str(truncate(fit_Mp[1])))
			f.write('\na\t\t\t'        + '[AU]\t\t'    + str(a)           + '\t\t\t' + str(truncate(fit_a[0]))  + '\t\t' + str(truncate(fit_a[1])))
			f.write('\ne\t\t\t'        + '[-]\t\t\t'   + str(e)           + '\t\t'   + str(truncate(fit_e[0]))  + '\t\t' + str(truncate(fit_e[1])))
			f.write('\nw\t\t\t'        + '[radians]\t' + str(w)           + '\t\t'   + str(truncate(fit_w[0]))  + '\t\t' + str(truncate(fit_w[1])))
			f.write('\nv0\t\t\t'       + '[m/s]\t\t'   + str(v0)          + '\t\t\t' + str(truncate(fit_v0[0])) + '\t\t' + str(truncate(fit_v0[1])))
			f.write('\nP\t\t\t'        + '[days]\t\t'  + str(truncate(P)) + '\t\t'   + str(truncate(fit_p[0]))  + '\t\t' + str(truncate(fit_p[1])))

			f.write('\n\nchi2reduced = ' + str(chi2reduced))
			f.write('\n\nMstar = ' + str(Ms) + ' [solar masses]')
			f.write('\nSNR = ' + str(SNR))


	# save rv data
	if len(insertSubfolder( settings.rv_data , subfolder)) > 0:
		savedatafile(insertSubfolder( settings.rv_data , subfolder), \
			(ts, rvs, \
				rvs-fitter._curve(ts, Mps, Ms, P, e, w, v0), \
				rvs-fitter._curve(ts, fit_Mp[0], Ms, fit_p[0], fit_e[0], fit_w[0], fit_v0[0])), \
			('time[days]', 'RV[m/s]', \
				'original_resid', \
				'fit_resid'))

	if len(insertSubfolder( settings.rv_curves , subfolder)) > 0:
		savedatafile(insertSubfolder( settings.rv_curves , subfolder), \
			(ts_full, fitter._curve(ts_full, Mps, Ms, P, e, w, v0), \
				fitter._curve(ts_full, fit_Mp[0], Ms, fit_p[0], fit_e[0], fit_w[0], fit_v0[0]), \
				fitter._curve(ts_full, fit_Mp[0], Ms, fit_p[0], fit_e[0], fit_w[0], fit_v0[0])-fitter._curve(ts_full, Mps, Ms, P, e, w, v0)), \
			('time[days]', 'original_RV[m/s]', \
				'fit_RV[m/s]', \
				'difference'))


	out( 'producing plots ...' )
		
	# LS periodigram
	if len(settings.periodogram_plot) > 0 or settings.periodogram_show:
		savePeriodogramPlot(freqs, periods, powers, P, fit_p[0], subfolder=subfolder)

	# orbit
	if len(settings.rv_plot) > 0 or settings.rv_show:
		rv_full_fit = fitter._curve(ts_full, fit_Mp[0], Ms, fit_p[0], fit_e[0], fit_w[0], fit_v0[0])
		rv_full_orig = fitter._curve(ts_full, Mps, Ms, P, e, w, v0)
		rv_fit = fitter._curve(ts, fit_Mp[0], Ms, fit_p[0], fit_e[0], fit_w[0], fit_v0[0])
		rv_orig = fitter._curve(ts, Mps, Ms, P, e, w, v0)

		# params = Mp, p, a, e, w, v0
		saveRVPlot(ts, ts_full, rvs, rv_full_fit, rv_full_orig, rv_fit, rv_orig, Ms, SNR, \
			(Mps, P, a, e, w, v0), \
			(fit_Mp[0], fit_p[0], fit_a[0], fit_e[0], fit_w[0], fit_v0[0]), \
			(fit_Mp[1], fit_p[1], fit_a[1], fit_e[1], fit_w[1], fit_v0[1]), \
			subfolder=subfolder)