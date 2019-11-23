
class Settings():
	def __init__(self, filename='input.config'):
		self.filename = filename

		# default values, Sun-Earth (no reason)
		self.starmass = 1.0
		self.rotation_period = 14.713
		self.start_phase = 0.0
		self.a = 1.0
		self.e = 0.0167
		self.planetmass = 1.0
		self.w = 0.0
		self.v0 = 0.0
		self.SNR = 1e4
		self.digits = 2
		self.input_files = []

		self.__initialize()

	def __initialize(self):
		with open(self.filename) as f:
		    for line in f:
		        self.__processline(line)

	def __processline(self, line):
		if len(line.strip()) == 0:
			return

		splitted = line.split('#')
		data = splitted[0].strip()

		if len(data.strip()) == 0:
			return

		name, value = data.split('=')
		name = name.strip().lower()
		value = value.strip()

		if name == 'starmass':
			self.starmass = float(value)
		elif name == 'rotation_period':
			self.rotation_period = float(value)
		elif name == 'start_rphase':
			self.start_rphase = float(value)
		elif name == 'start_mphase':
			self.start_mphase = float(value)
		elif name == 'a':
			self.a = float(value)
		elif name == 'e':
			self.e = float(value)
		elif name == 'planetmass':
			self.planetmass = float(value)
		elif name == 'w':
			self.w = float(value)
		elif name == 'v0':
			self.v0 = float(value)
		elif name == 'snr':
			self.SNR = float(value)
		elif name == 'digits':
			self.digits = int(value)

		if name == 'time_source':
			self.time_source = value.strip()
		elif name == 'time_from':
			self.time_from = float(value)
		elif name == 'time_to':
			self.time_to = float(value)
		elif name == 'time_step':
			self.time_step = float(value)
		elif name == 'time_randomize':
			self.time_randomize = self.__getBoolean(value)
		elif name == 'time_data':
			self.time_data = value.strip()

		if name == 'outputfile':
			self.outputfile = value.strip() #   TODO: remove restricted characters
		elif name == 'logfile':
			self.logfile = value.strip()
		elif name == 'rv_curves_points':
			self.rv_curves_points = int(value)
		elif name == 'rv_data':
			self.rv_data = value.strip()
		elif name == 'rv_curves':
			self.rv_curves = value.strip()
		elif name == 'rv_plot':
			self.rv_plot = value.strip()
		elif name == 'rv_show':
			self.rv_show = self.__getBoolean(value)
		elif name == 'progressbar':
			self.progressbar = self.__getBoolean(value)
		elif name == 'profiles_path':
			self.profiles_path = value.strip()

		if name == 'periodogram_low':
			self.periodogram_low = float(value)
		elif name == 'periodogram_high':
			self.periodogram_high = float(value)
		elif name == 'periodogram_step':
			self.periodogram_step = float(value)
		elif name == 'periodogram_data':
			self.periodogram_data = value.strip()
		elif name == 'periodogram_plot':
			self.periodogram_plot = value.strip()
		elif name == 'periodogram_show':
			self.periodogram_show = self.__getBoolean(value)
		elif name == 'periodogram_units':
			self.periodogram_units = value.strip()
		elif name == 'max_eccentricity':
			self.max_eccentricity = float(value)
		elif name == 'fit_stellar_rotation':
			self.fit_stellar_rotation = self.__getBoolean(value)

		if name == 'add_input_file':
			spl = value.split(',')
			vals = [float(spl[0].strip()), spl[1].strip()]
			self.input_files.append(vals)
		elif name == 'magnetic_period':
			self.magnetic_period = float(value)


	def __getBoolean(self, value):
		return value.strip().lower() == 'true'