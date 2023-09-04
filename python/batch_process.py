import sys
import os
import subprocess
import time
import tempfile
import numpy as np


# number of runs
lf = 8

# number of cores
num_cores = 8
nc=min([num_cores, lf])

outputFileName='/tmp/myoutput.npz'


if __name__ == "__main__":
	"""
		submit dummy jobs using subprocess
	"""
	p=[None]*nc
	log=[None]*nc
	submission=[None]*nc
	for i in range(nc):
		p[i] = subprocess.Popen('', shell=True)
	time.sleep(1)

	"""
	   main loop 
	"""
	strOutputFiles=[]
	strOutputFiles_h=[]
	i = 0
	iter1=0
	while True:
		runsProcessing = False
		# go through all cores
		for j in range(nc):
			# submit a new run
			if p[j].poll() == 0 and i < lf:
				# run command
				strOutputFile = tempfile.NamedTemporaryFile()
				strOutputFiles_h.append(strOutputFile)
				strOutputFiles.append(strOutputFile.name)
				runs=sys.executable + ' ' + "main.py '" + strOutputFile.name + "'"
				
				submission[j]=runs
				print(runs)
				# submit the job
				p[j]=subprocess.Popen(runs, shell=True,stdout=subprocess.DEVNULL, \
					stderr=subprocess.STDOUT)
				i += 1

	
			# resubmit if it didn't work
			if (p[j].poll() != 0) and (p[j].poll() != None):
				print('resubmitting')
				p[j].terminate()
				p[j]=subprocess.Popen(submission[j], shell=True,stdout=subprocess.DEVNULL, \
					stderr=subprocess.STDOUT)

			# check if runs are still going
			runsProcessing = runsProcessing or (p[j].poll() != 0)

		# break out if the runs are not processing and we are past the last file
		if not(runsProcessing) and (i>=lf):
			break
		iter1 += 1
		time.sleep(10)
		
	print(strOutputFiles)
	
	
	"""
		OK when run, read in the files and write to a new file
	"""
	temp=np.load(strOutputFiles[0] + '.npz')
	t = temp['t']
	Patms 	= np.zeros((len(t),len(strOutputFiles)))
	dNs 	= np.copy(Patms)
	mole1s 	= np.copy(Patms)
	for i in range(len(strOutputFiles)):
		temp = np.load(strOutputFiles[i] + '.npz')
		Patms[:,i]	=temp['Patm']
		dNs[:,i]	=temp['dN']
		mole1s[:,i]	=temp['mole1']
	"""
		now save it
	"""
	np.savez(outputFileName,t=t,Patm=Patms,dN=dNs,mole1=mole1s)
	
	
	
	
	"""
		close files - this will delete them
	"""
	for i in range(len(strOutputFiles)):
		strOutputFiles_h[i].close()
		os.remove(strOutputFiles[i] + '.npz')
		