import os

## check if file exist then ask if you want to delete old one and create new one. 
## if answer is yes then deelte the old and create this one new. 
## each line will contain four parameters. 
## taskname   cfg.py  datasetname  numberofdiles
## cfg.py is configurable because data and MC will have different configurations.
## And number of files canbe used as number of lumis in that case. 
postfix="_MC25ns_LegacyMC_20170328"

##var='/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'

outfile = open('allmonoHsamples.txt','w')
datasetfile = open('allMCMonoHbb2016Legacy.txt','r')
#datasetfile = open('datsetsMonoHbb.txt','r')
for dataset in datasetfile:
    a=dataset.split('/')
    b=a[1]
    b=b.rstrip()
    dataset=dataset.rstrip()
    line=b+postfix+" treeMaker_Summer16_cfg.py "+dataset+" 1 \n"
    print line
    outfile.write(line)


outfile.close()


def submit():
    print "submitting"
    f = open('allmonoHsamples.txt','r')
    for line in f:
        print line
        a=[]
        b=[]
        c=[]
        d=[]
        a,b,c,d = line.split()
        datasetdetail=[a,b,c,d]
        print datasetdetail
        os.system('cp -p crabConfig_MC.py crabConfig.py')
        os.system('crab submit General.requestName='+datasetdetail[0]+' JobType.psetName='+datasetdetail[1]+' Data.inputDataset='+datasetdetail[2]+' Data.unitsPerJob='+datasetdetail[3])
        
        



def status(crabdirname):
    import os
    os.system ("./Statusall.sh "+crabdirname)


def resubmit(crabdirname):
    import os
    os.system ("./Resubmit.sh "+crabdirname)

    
    

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
