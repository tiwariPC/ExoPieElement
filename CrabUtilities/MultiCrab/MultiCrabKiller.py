from os import listdir, system
from os.path import isfile, join

workname='monoH_2016DATA_20190411'

mypath="crab_"+workname+"/"
dirs = [f for f in listdir(mypath) if not isfile(join(mypath, f))]

njobs=len(dirs)

for ind in range(njobs):	
	print "\n==========================\nKilling job "+str(ind+1)+" of "+str(njobs)+"\n==========================\n"
	system("crab kill -d "+mypath+dirs[ind])
