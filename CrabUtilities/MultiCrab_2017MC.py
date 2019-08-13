import os
import sys, optparse,argparse

usage = "use only one option at one time, if you use c, r, k then do provide crabdir  "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-s", "--submit", action="store_true",  dest="submit")
parser.add_argument("-c", "--status", action="store_true",  dest="status")
parser.add_argument("-r", "--resubmit", action="store_true",  dest="resubmit")
parser.add_argument("-k", "--kill", action="store_true",  dest="kill")
parser.add_argument("-d", "--crabdir", dest="crabdir" )

args = parser.parse_args()

print args
print len(sys.argv)
def submit(datasetdetail):
    print "submitting"
    
    print datasetdetail
    os.system('crab submit -c crabConfig_2017_MC.py General.requestName='+datasetdetail[0]+' Data.inputDataset='+datasetdetail[1]+' Data.unitsPerJob='+datasetdetail[2])
    #print ('crab submit -c crabConfig_2017_MC.py General.requestName='+datasetdetail[0]+' Data.inputDataset='+datasetdetail[1]+' Data.unitsPerJob='+datasetdetail[2])




def status(crabdirname):
    print "cehcking the status of all jobs "
    os.system("ls -1 "+crabdirname + " >& crabtasklist.txt")
    for itask in open('crabtasklist.txt'):
        print "checking status of ",itask.rstrip()
    
        os.system('crab status -d '+crabdirname+'/'+itask)



def resubmit(crabdirname):
    print "resubmitting all the failed jobs"
    os.system("ls -1 "+crabdirname + " >& crabtasklist.txt")
    for itask in open('crabtasklist.txt'):
        print "resubmiting failed jobs for the tast of ",itask.rstrip()
    
        os.system('crab resubmit -d '+crabdirname+'/'+itask)
        


def kill(crabdirname):
    print "killing all the jobs in the dir you have provided  all the failed jobs: use it only when you are fully sure, this will be used only when there was some error"
    os.system("ls -1 "+crabdirname + " >& crabtasklist.txt")
    for itask in open('crabtasklist.txt'):
        print "killing all the jobsfor the tast of ",itask.rstrip()
    
        os.system('crab kill -d '+crabdirname+'/'+itask)
        


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





textfilename='raman.txt'


if args.submit:
    f = open(textfilename,'r')
    for line in f:
        print line
        a,b,c = line.split()
        datasetdetail=[a,b,c]
        submit(datasetdetail)

if (args.status or args.resubmit or args.kill ) and len(sys.argv)<3 :
    print "insufficienct input provided, please look at the script MultiCrab_2017MC.py for usage. or ask Raman"
    
if (args.status or args.resubmit or args.kill ) and len(sys.argv)==3 :
    print "two argument provided"
    if args.status: status(args.crabdir)

    if args.resubmit: resubmit(args.crabdir)
    
    if args.kill: kill(args.crabdir) 
