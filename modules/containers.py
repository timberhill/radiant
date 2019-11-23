
class VoigtLine():
	def __init__(self, center, flux, g_fwhm, l_fwhm):
		self.center = center
		self.flux = flux
		self.g_fwhm = g_fwhm
		self.l_fwhm = l_fwhm