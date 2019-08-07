import os
import sys


def submit(datasetdetail):
    print "submitting"
    
    print datasetdetail
    #os.system('crab submit -c crabConfig_2017_MC.py General.requestName='+datasetdetail[0]+' Data.inputDataset='+datasetdetail[1]+' Data.unitsPerJob='+datasetdetail[2])
    print ('crab submit -c crabConfig_2017_MC.py General.requestName='+datasetdetail[0]+' Data.inputDataset='+datasetdetail[1]+' Data.unitsPerJob='+datasetdetail[2])
    


'''
def status(crabdirname):
    import os
    os.system ("./Statusall.sh "+crabdirname)
    

## Add a help or usage function here 
def help() :
    print "this is under progress"

    



## safety check
## improve this
if len(sys.argv) < 2 :
    print "insufficient options provided see help function "
    help()
    exit (1)


## submit jobs 
if len(sys.argv) == 2 :
    if sys.argv[1] == "submit" :
        submit()


## check status of jobs 
## send the crab directory 
if len(sys.argv) == 3 : 
    if sys.argv[1] == "status" :
        crabdir = sys.argv[2]
        status(crabdir)


'''



f = open('backgroundList_2017_miniaodV2.txt','r')
for line in f:
    print line
    a,b,c = line.split()
    datasetdetail=[a,b,c]
    
    submit(datasetdetail)

    
