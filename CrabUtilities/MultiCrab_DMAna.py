import os

### Preparing the dataset, crab and pset configuration files
def prepare(dataset):
    os.system('rm -rf ana_inputdataset.txt')
    f = open(dataset,'r')
    for line in f:
        print line
        a=[]
        b=[]
        a,b = line.split()
        datasetdetail=[a,b]
        print datasetdetail
        name='crab status -d crab_projects_MonoHStep4/crab_step4-'+datasetdetail[1]+' | grep -a "Output dataset:" | awk -v my_var='+datasetdetail[1]+' \'{print my_var" treeMaker_Summer16_cfg.py "$3" 1"}\'>> ana_inputdataset.txt'
        print name
        os.system(name)


def submit(official):
    print "submitting"
    f = open('ana_inputdataset.txt','r')
    for line in f:
        print line
        a=[]
        b=[]
        c=[]
        d=[]
        a,b,c,d = line.split()
        datasetdetail=[a,b,c,d]
        print datasetdetail
        if not official:
            name='crab submit -c crabConfig.py General.requestName='+datasetdetail[0]+' JobType.psetName='+datasetdetail[1]+' Data.inputDataset='+datasetdetail[2]+' Data.unitsPerJob='+datasetdetail[3]+' Data.inputDBS=phys03'
        else:
            name='crab submit -c crabConfig.py General.requestName='+datasetdetail[0]+' JobType.psetName='+datasetdetail[1]+' Data.inputDataset='+datasetdetail[2]+' Data.unitsPerJob='+datasetdetail[3]
        print name
        os.system(name)


def status(crabdirname):
    import os
    os.system ("./Statusall.sh "+crabdirname)
    

## Add a help or usage function here 
def help() :
    print "this is under progress"

    


####################################################################################################################################################
####################################################################################################################################################
## this will control the functions   ##
## convert this to python main.      ##
####################################################################################################################################################
####################################################################################################################################################
import os
import sys
print sys.argv

## safety check
## improve this
if len(sys.argv) < 2 :
    print "insufficient options provided see help function "
    help()
    exit (1)




## prepare jobs
if sys.argv[1] == "prepare" :
    if len(sys.argv) != 3 :
        print "Usuage: python MultiCrab_DMAna.py prepare step4_inputdataset.txt"
        exit (1)
    else:
        print "preparing for the submission"
        datasetTxt = sys.argv[2]
        prepare(datasetTxt)


## submit jobs 
if len(sys.argv) == 2 :
    if sys.argv[1] == "submit" :
        submit(False)

if len(sys.argv) == 3 :
    if sys.argv[1] == "submit" and sys.argv[2] == "official":
        submit(True)


## check status of jobs 
## send the crab directory 
if len(sys.argv) == 3 : 
    if sys.argv[1] == "status" :
        crabdir = sys.argv[2]
        status(crabdir)




