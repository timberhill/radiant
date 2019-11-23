import numpy as np
from modules.containers import VoigtLine
from modules.redshifter import RedshiftMeasurer
from modules.spectrumgen_voigts import VoigtLinesGenerator
from modules.rvsimulator_singletimeorbit import SingleTimeOrbitRVSimulator
from modules.s1reader import S1Reader
from modules.settings import Settings
from modules.helpers import ensurePathExists, savedatafile, truncate
from modules.s1spectrafetcher import S1SpectraFetcher

class RVSpectraSimulator():
	def __init__(self, Mstar=1.0, Mplanet=1.0, a=1.0, e=0.0, w=0.0, v0=0.0, SNR=0.0, R=1e5, times=np.arange(0.0, 365.25, 0.1)):
		self.G = 8.88677e-10 # AU^3 Mearth^-1 day^-2
		self.SunInEarthMasses = 332978.9 # 1 solar mass in Earth masses

		self.Mstar = Mstar # solar masses
		self.Mplanet = Mplanet # Earth masses
		self.a = a # AU
		self.e = e # 0..1
		self.w = w # argument of periastron, radians
		self.v0 = v0 # barycenter velocity (constant component)

		self.SNR = SNR # set to 0 to eliminate all noise
		self.R = R # spectral resolution Lambda/dLambda
		self.times = times

		self.simulator = SingleTimeOrbitRVSimulator(times, Mstar, Mplanet, a, e, w, v0)
		self.redshiftMeasurer = RedshiftMeasurer()

		self.lines = np.asarray([])
		self.continuum = 1.0
		self.absorption = True
		self.margin = 5 # angstroems

		self.maxshift = 1e3 # m/s
		self.sfactor = 300 # for spline interpolation of the line

		self.settings = Settings()


	def _getLinesProperties(self):
		if len(self.lines) == 0:
			raise Exception('Provide at least one item of type <VoigtLine> for the spectrum in <RVSpectraSimulator.lines[]>')

		centers = [VoigtLine.center for VoigtLine in self.lines]
		fluxes = [VoigtLine.flux for VoigtLine in self.lines]
		g_fwhms = [VoigtLine.g_fwhm for VoigtLine in self.lines]
		l_fwhms = [VoigtLine.l_fwhm for VoigtLine in self.lines]

		return centers, fluxes, g_fwhms, l_fwhms

	def _getModelSpectrum(self):
		centers, fluxes, g_fwhms, l_fwhms = self._getLinesProperties()
		model = VoigtLinesGenerator(centers=centers, margin=self.margin, rv=0, g_fwhm=g_fwhms, l_fwhm=l_fwhms, R=self.R, fluxes=fluxes, cont=self.continuum)

		return model.getSpectrum(absorption=self.absorption)


	def addLine(self, voigt):
		self.lines = np.append(self.lines, voigt)

	def removeAllLines(self):
		self.lines = np.asarray([])

	def period(self):
		return np.power(4.0*np.power(np.pi, 2)*np.power(self.a, 3) / (self.G*(self.Mstar * self.SunInEarthMasses + self.Mplanet)), 0.5)

	def getRVCurve(self, times=[], dt=None):
		if dt != None:
			times = np.arange(min(self.times), max(self.times), dt)

		if len(times) > 0:
			simulator = SingleTimeOrbitRVSimulator(times, self.Mstar, self.Mplanet, self.a, self.e, self.w, self.v0)
			return simulator.getRVCurve()

		return self.simulator.getRVCurve()

	def getRVData(self, progressbar=False, template='average', outfunction=None, times_full=None, cycle=True):
		ts0, rvs0 = self.getRVCurve()
		fetcher = S1SpectraFetcher()
		time_from = self.settings.time_from
		
		wl_model, flux_model = [], []
		lines_wls = []
		lines_fluxes = []

			
		if progressbar and outfunction != None:
			outfunction( 'producing profiles ...' )
		# get lines first

		points = []
		if progressbar:
			from progressbar import ProgressBar,SimpleProgress,ETA
			pbar = ProgressBar(widgets=['                        ', SimpleProgress(), ', ', ETA()])
			points = pbar(range(0,len(rvs0)))
		else:
			points = range(0,len(rvs0))
			
		for rvi in points:
			# stellar rotation phase
			start_rphase = self.settings.start_rphase
			rotation_period = self.settings.rotation_period
			rphase = start_rphase + (ts0[rvi] - time_from) / rotation_period

			# magnetic cycle phase
			start_mphase = self.settings.start_mphase
			magnetic_period = self.settings.magnetic_period
			mphase = start_mphase + (ts0[rvi] - time_from) / magnetic_period

			# get profile
			wls, fluxes = fetcher.getProfile_old(rphase, mphase, rv=rvs0[rvi], timepoint=ts0[rvi])
			#wls, fluxes = fetcher.getProfile(time=ts0[rvi], rphase=rphase, rv=rvs0[rvi])

			if self.SNR > 0.0:
				# add Gaussian noise with sigma = 1.0/SNR
				noise = np.random.normal(0, 1.0/self.SNR, (len(fluxes)))
				fluxes = np.asarray([x+y for x, y in zip(fluxes, noise)])

			# save profile to file
			if len(self.settings.profiles_path) > 0:
				path = self.settings.profiles_path.replace(r'{N}', str(rvi)).replace(r'{RV}', str(truncate(rvs0[rvi]))).replace(r'{time}', str(truncate(ts0[rvi])))
				ensurePathExists(path)
				lsd_rvs = wls#reader.wlToRv(phase=rphase, wls=wls)
				savedatafile(path, \
					(wls, lsd_rvs, fluxes), ('wl[A]', 'rv_shift[m/s]', 'flux'), \
					comment='RV = ' + str(rvs0[rvi]) + ' [m/s],\tt = ' + str(ts0[rvi]) + ' [days],\tstellar rotation phase: ' + str(rphase))


			lines_wls.append(wls)
			lines_fluxes.append(fluxes)

		self.R = fetcher.R

		if progressbar and outfunction != None:
			outfunction( 'computing template ...' )

		if template == 'average':
			# average lines and use it as a model
			for rvi in range(0,len(rvs0)):
				if rvi == 0:
					for wli in range(0,len(lines_wls[rvi])):
						wl_model.append(lines_wls[rvi][wli])
						flux_model.append(lines_fluxes[rvi][wli])
				else:
					for wli in range(0,len(lines_wls[rvi])):
						flux_model[wli] += lines_fluxes[rvi][wli]

				lines_wls[rvi] = np.asarray(lines_wls[rvi])
				lines_fluxes[rvi] = np.asarray(lines_fluxes[rvi])
			n = len(rvs0)
			for i in range(0, len(flux_model)):
				flux_model[i] /= n

			wl_model = np.asarray(wl_model)
			flux_model = np.asarray(flux_model)

		elif template == 'model':
			wl_model, flux_model = self._getModelSpectrum()
		else:
			raise Exception('Unknown template selected in RVSpectraSimulator.getRVData(). Available options: "model" and "average".')


		if progressbar and outfunction != None:
			outfunction( 'measuring RVs ...' )

		points = []
		if progressbar:
			from progressbar import ProgressBar,SimpleProgress,ETA
			pbar = ProgressBar(widgets=['                        ', SimpleProgress(), ', ', ETA()])
			points = pbar(rvs0)
		else:
			points = rvs0

		i = 0
		rvdata = []
		for rv in points:
			# measure redshift
			val = self.redshiftMeasurer.getRV(lines_wls[i], lines_fluxes[i], wl_model, flux_model, maxshift=self.maxshift, R=self.R, fit_type='linespline', sfactor=self.sfactor)
			rvdata.append(val)
			i += 1
			
		return ts0, np.asarray(rvdata), ts0, rvs0