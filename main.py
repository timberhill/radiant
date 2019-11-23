import sys
import numpy as np
from modules.datagen import RVSpectraSimulator
from modules.settings import Settings
from modules.helpers import out, ensurePathExists, savedatafile, getObsPoints, truncate
from detector import runDetector
settings = Settings()


# star
Ms = settings.starmass

# orbit
Mps = settings.planetmass
a = settings.a
e = settings.e
w = settings.w
v0 = settings.v0

SNR = settings.SNR

# read or generate observations time points
times = []
if len(settings.time_source) > 0:
	times = np.loadtxt(settings.time_source, usecols=(0,), unpack=True)
	times = np.sort(times)
else:
	times = getObsPoints(settings.time_from, settings.time_to, settings.time_step, settings.time_randomize)

# save time points
ensurePathExists(settings.time_data)
with  open(settings.time_data,'w') as f:
	for t in times:	
		f.write(str(t) + '\n')


tmin = min(times)
tmax = max(times)
ts_full = np.arange(tmin, tmax, (tmax-tmin)/settings.rv_curves_points)


gen = RVSpectraSimulator(Mstar=Ms, Mplanet=Mps, a=a, e=e, w=w, v0=v0, SNR=SNR, times=times)

G = gen.G
Mtot = gen.Mstar*gen.SunInEarthMasses + gen.Mplanet
P = gen.period()
R = gen.R

# estimate semi amplitude
t1 = np.power(2*np.pi*G / P, 1.0/3.0)
t2 = gen.Mplanet / np.power(Mtot, 2.0/3.0)
t3 = 1.0 / np.power(1 - e**2, 0.5)
K = t1*t2*t3*(1.496e11/86400) # (AUmeters/DAYseconds)
gen.maxshift = K * 1.1 + 20 #1500 # m/s
gen.sfactor = 100

out( '[INPUT] star mass = ' + str(Ms) + ' [Msun]')
out( '[INPUT] planet mass = ' + str(Mps) + ' [Mearth]')
out( '[INPUT] a =  ' + str(a) + ' [AU]')
out( '[INPUT] e =  ' + str(e))
out( '[INPUT] w =  ' + str(w/np.pi) + ' * pi')
out( '[INPUT] v0 = ' + str(v0))
out( '[INPUT] period = ' + str(truncate(P)) + ' [days]\n')

out( 'generating orbit ...' )
ts, rvs, ts0, rvs0 = gen.getRVData(progressbar=settings.progressbar, template='average', outfunction=out, times_full=ts_full)
errs = np.asarray([1.0]*len(ts))

# save rv data
if len(settings.rv_data) > 0:
	savedatafile(settings.rv_data, \
		(ts, rvs), \
		('time[days]', 'RV[m/s]'))

runDetector(ts, rvs, ts0, rvs0, errs, Ms, Mps, a, P, e, w, v0, SNR)

if settings.fit_stellar_rotation:
	runDetector(ts, rvs, ts0, rvs0, errs, Ms, Mps, a, P, e, w, v0, SNR, fit_rotation=True, subfolder='subtracted_rotation')


out( 'DONE\n\n' )