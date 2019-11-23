
import datetime, time, os, random
from numpy import sort, asarray
from bisect import bisect_left
from modules.settings import Settings

def ensurePathExists(filename):
	if not os.path.exists(os.path.dirname(filename)):
		try:
			os.makedirs(os.path.dirname(filename))
		except OSError as exc: # Guard against race condition
			pass

def out(message):
	st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
	message = '[' + st + ']\t' + message
	print(message)
	settings = Settings()
	if len(settings.logfile) > 0:
		ensurePathExists(settings.logfile)
		with open(settings.logfile, 'a') as f:
			f.write(message + '\n')

def savedatafile(filename, lists, captions, separator='\t', comment=''):
	n = len(lists[0])
	ensurePathExists(filename)
	with open(filename, 'w+') as file:
		if len(comment) > 0:
			file.write('# ' + comment + '\n')

		line = '# '
		for x in captions:
			line += (str(x) + separator)
		file.write(line.strip() + '\n\n')

		for i in range(0, n):
			line = ''
			for x in lists:
				line += (str(x[i]) + separator)
			file.write(line.strip() + '\n')

def getObsPoints(left, right, step, randomize=True):
	xs = []
	phase = left

	def mult(randomize):
		if randomize:
			return random.random()
		return 1.0

	#phase += mult(randomize) * step
	while phase <= right:
		xs.append(phase)
		phase += mult(randomize) * step

	return sort(xs)

def truncate(f, n=2):
	settings = Settings()
	n = settings.digits
	s = '%.12f' % round(f, n)
	i, p, d = s.partition('.')
	return '.'.join([i, (d+'0'*n)[:n]])

def getClosest(myList, myNumber, isSorted=True):
	if isSorted:
		pos = bisect_left(myList, myNumber)

		if pos == 0:
			return pos, myList[0]
		if pos == len(myList):
			return pos-1, myList[-1]

		before = myList[pos - 1]
		after = myList[pos]
		if after - myNumber < myNumber - before:
			return pos, after
		else:
			return pos, before
	else: # bruteforce for unsorted lists
		index = -1
		mindiff = max(myList) - min(myList)
		for i, x in enumerate(myList):
			diff = abs(x - myNumber)
			if diff < mindiff:
				index = i
				mindiff = diff

		return index, myList[index]

def insertSubfolder(filename, subfolder):
	if subfolder == '' or filename == '':
		return filename

	splitted = filename.split('/')
	new = []
	for i in range(0,len(splitted)):
		if i == len(splitted)-1:
			new.append(subfolder)
		new.append(splitted[i])

	return '/'.join(new)
