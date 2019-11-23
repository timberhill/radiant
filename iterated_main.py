from subprocess import call, PIPE
from numpy import power, pi, sort
from sklearn.datasets.samples_generator import make_blobs
from modules.helpers import truncate, out, ensurePathExists
import sys

# python iterated_main.py <caseX> <cluster> <ppp> (points per period)
# settings
periodsNumber = 3
pointsPerPeriod = int(sys.argv[3]) #30
samples = 10
clustered = sys.argv[2] == "cluster"
template = open('input_template_{}.config'.format(sys.argv[1]), 'r').read()
foldername = "output"


def pfroma(a, Ms, Mp):
	G = 8.88677e-10 # AU^3 Mearth^-1 day^-2
	SunInEarthMasses = 332978.9 # 1 solar mass in Earth masses
	return power( 4.0*pi**2*a**3 / (G*(Mp + Ms*SunInEarthMasses)), 0.5)


def generate_observations(left, right, n, cluster_n=1, cluster_std=0.3, seed=None):
	x, y = make_blobs(n_samples=n, centers=n//cluster_n, n_features=1, cluster_std=cluster_std, center_box=(left, right), random_state=seed)
	return x.T[0]


starmass = 0.856 # solar masses
planet_masses = [1, 2, 5, 10, 20, 50, 100, 159, 318, 636, 1589] # Earth masses
planet_orbits = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0] # AU
tbounds = [130, 2570]
time_from = tbounds[0]

out("Go get a coffee, it's gonna take a while. Trust me\n")

N = len(planet_masses) * len(planet_orbits)
count = 1
for mpi, mp in enumerate(planet_masses):
	for ai, a in enumerate(planet_orbits):
		# get parameters
		p = pfroma(a, starmass, mp)
		time_from = tbounds[0]
		length = periodsNumber * p
		time_to = tbounds[0] + length
		if time_to > tbounds[1]: # fix for observations limit
			time_to = tbounds[1]

		# THIS IS FOR MULTIPLE CALCULATIONS of same orbit
		unused = (tbounds[1] - tbounds[0]) - (time_to - time_from)
		sample_step = unused / samples

		for i in range(0, samples):
			out('\n\nStarting orbit ' + str(count) + '.' + str(i) + '/' + str(N) + ', Mp = ' + str(truncate(mp)) + ' Me, a = ' + str(truncate(a)) + ' AU')
			
			time_from	= tbounds[0] + i * sample_step
			time_to		= tbounds[0] + length + i * sample_step
			if time_to > tbounds[1]:
				time_to = tbounds[1]

			new_config = template
			output_folder = foldername + '/mp_' + str(truncate(mp)) + '_a_' + str(truncate(a)) + '_' + str(i)

			time_step = 2 * (time_to - time_from) / (periodsNumber*pointsPerPeriod)
			n = pointsPerPeriod * periodsNumber

			if clustered:
				n_clusters = (periodsNumber*5) # 5 clusters per period
				cluster_n = n//n_clusters
				cluster_std = 0.4 # ~1 night, OR:
				if (time_to - time_from)/n_clusters <= time_step*5:
					cluster_std = time_step/30
				times = generate_observations(left=time_from, right=time_to, n=n, cluster_n=cluster_n, cluster_std=cluster_std)
			else:
				times = generate_observations(left=time_from, right=time_to, n=n, cluster_n=1)

			# save the time input file
			ensurePathExists(output_folder + '/times.dat')
			with open(output_folder + '/times.dat','w') as f:
				for t in times:	
					f.write(str(t) + '\n')

			if clustered:
				new_config = new_config.replace(r'{time_source}', output_folder + '/times.dat')
			else:
				new_config = new_config.replace(r'{time_source}', '#times.dat')

			# put parameters into the template
			new_config = new_config.replace(r'{a}', str(a))
			new_config = new_config.replace(r'{planetmass}', str(mp))
			new_config = new_config.replace(r'{time_from}', str(time_from))
			new_config = new_config.replace(r'{time_to}',   str(time_to))
			new_config = new_config.replace(r'{time_step}', str(time_step))
			new_config = new_config.replace(r'{output_folder}', str(output_folder))

			# save input file to "root" to be used and to the output folder
			with open("input.config", "w") as f1:
				f1.write(new_config)

			ensurePathExists(output_folder + "/input.config")
			with open(output_folder + "/input.config", "w") as f2:
				f2.write(new_config)

			# run code
			call(['python','main.py'], stdout=PIPE)

		# repeat
		count += 1

out('\nWow, all done. Have a cookie.')
