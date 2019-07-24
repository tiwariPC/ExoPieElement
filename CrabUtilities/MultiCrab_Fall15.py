import os

## check if file exist then ask if you want to delete old one and create new one. 
## if answer is yes then deelte the old and create this one new. 
## each line will contain four parameters. 
## taskname   cfg.py  datasetname  numberofdiles
## cfg.py is configurable because data and MC will have different configurations.
## And number of files canbe used as number of lumis in that case. 
postfix=""

##var='/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'

outfile = open('datasetdetails_Fall15.txt','w')

datasetfile = open('datsetsFall15.txt','r')
for dataset in datasetfile:
    a=dataset.split('/')
    b=a[1]
    b=b.rstrip()
    dataset=dataset.rstrip()
    line=b+postfix+" treeMaker_Fall15_Nocleaning_cfg.py "+dataset+" 1 \n"
    print line
    outfile.write(line)


outfile.close()


def submit():
    print "submitting"
    f = open('datasetdetails_Fall15.txt','r')
    for line in f:
        print line
        a=[]
        b=[]
        c=[]
        d=[]
        a,b,c,d = line.split()
        datasetdetail=[a,b,c,d]
        print datasetdetail
        os.system('crab submit  General.requestName='+datasetdetail[0]+' JobType.psetName='+datasetdetail[1]+' Data.inputDataset='+datasetdetail[2]+' Data.unitsPerJob='+datasetdetail[3])
    #name =  'crab submit General.requestName='+datasetdetail[0]+' JobType.psetName='+datasetdetail[1]+' Data.inputDataset='+datasetdetail[2]+' Data.unitsPerJob='+datasetdetail[3]
    #print name 
        
        



STATUS='''
source /cvmfs/cms.cern.ch/crab3/crab.sh

if [ -f crabtaskslist.txt ] ; then
    rm crabtaskslist.txt
else echo "file crabtaskslist.txt does not exist will list the directory and create one now"
fi

if [ $# -lt 1 ] ; then
    echo $crabProjectDir " does not exist"
    echo "give correct crabProjectDir name"
    exit 0;
fi

crabProjectDir=$1
if [ -d $crabProjectDir ] ; then
    ls -1 $crabProjectDir >& crabtaskslist.txt
    for which in `less crabtaskslist.txt` ; do
        echo " "
        echo " "
        echo "---------------- CHECKING STATUS OF THE JOBS ---------------------"
        echo $which
        echo " "
        echo " "
        crab status $crabProjectDir/$which
    done
else
    echo $crabProjectDir " does not exist"
    echo "give correct crabProjectDir name"
fi

'''

statussh = open('Statusall.sh','w')
statussh.write(STATUS)
statussh.close()
os.system('chmod 777 Statusall.sh')

RESUBMIT='''
source /cvmfs/cms.cern.ch/crab3/crab.sh

if [ -f crabtaskslist.txt ] ; then
    rm crabtaskslist.txt
else echo "file crabtaskslist.txt does not exist will list the directory and create one now"
fi

if [ $# -lt 1 ] ; then
    echo $crabProjectDir " does not exist"
    echo "give correct crabProjectDir name"
    exit 0;
fi

crabProjectDir=$1
if [ -d $crabProjectDir ] ; then
    ls -1 $crabProjectDir >& crabtaskslist.txt
    for which in `less crabtaskslist.txt` ; do
        echo " "
        echo " "
        echo "---------------- Resubmiting jobs ---------------------"
        echo $which
        echo " "
        echo " "
        crab resubmit $crabProjectDir/$which
    done
else
    echo $crabProjectDir " does not exist"
    echo "give correct crabProjectDir name"
fi
'''

resubmitsh = open('resubmit.sh','w')
resubmitsh.write(RESUBMIT)
resubmitsh.close()
os.system('chmod 777 resubmit.sh')


def status(crabdirname):
    import os
    os.system ("./Statusall.sh "+crabdirname)
    


def resubmit(crabdirname):
    import os
    os.system ("./resubmit.sh "+crabdirname)
    

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


if len(sys.argv) == 3 : 
    if sys.argv[1] == "resubmit" :
        crabdir = sys.argv[2]
        resubmit(crabdir)





#  LocalWords:  MonoHToBBarMZp
