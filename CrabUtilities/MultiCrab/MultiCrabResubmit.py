from os import listdir, system
from os.path import isfile, join

workname='2017_data_20190425'

mypath="crab_"+workname+"/"
dirs = [f for f in listdir(mypath) if not isfile(join(mypath, f))]

njobs=len(dirs)

for ind in range(njobs):	
	print "\n==========================\nResubmiting job "+str(ind+1)+" of "+str(njobs)+"\n==========================\n"
	system("crab resubmit -d "+mypath+dirs[ind])
