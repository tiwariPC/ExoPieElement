import os
import sys, optparse,argparse

usage = "use only one option at one time, if you use c, r, k then do provide crabdir  "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-s", "--submit", action="store_true",  dest="submit")
parser.add_argument("-c", "--status", action="store_true",  dest="status")
parser.add_argument("-r", "--resubmit", action="store_true",  dest="resubmit")
parser.add_argument("-k", "--kill", action="store_true",  dest="kill")
parser.add_argument("-d", "--crabdir", dest="crabdir" )
parser.add_argument("-l", "--ss", action="store_true",  dest="ss") ## status summary

args = parser.parse_args()

print args
print len(sys.argv)
def submit(datasetdetail):
    print "submitting"
    print datasetdetail
    os.system('crab submit -c crabConfig_2018_Signal.py General.requestName='+datasetdetail[0]+' Data.inputDataset='+datasetdetail[1])


def cleannumber(N):
    return str(str(N).replace("(","").replace(")",""))

def statussummary(logfilename):
    ## get the task name
    os.system ("cat "+logfilename+" "+ " | grep \"CRAB project directory\" > tmp")
    ftmp=open("tmp")
    firstline = (ftmp.readline()).rstrip().split("/")[-1]
    ftmp.close()

    totaljobs=0
    finishedjobs=0
    transferringjobs=0
    failedjobs=0
    idlejobs=0
    submittedjobs=0



    os.system("cat statuslog | grep \"finished\" > tmp")
    if not os.stat('tmp').st_size==0:
        cols=(open("tmp")).readline().rstrip().split(" ")[-1]

        ## get the total jobs
        totaljobs=cols.split("/")[-1]


        ## get the finished jobs
        finishedjobs=cols.split("/")[0]


    ## get the transferring jobs
    os.system("cat statuslog | grep \"transferring\" > tmp")
    if not os.stat('tmp').st_size==0:
        cols=(open("tmp")).readline().rstrip().split(" ")[-1]
        transferringjobs=cols.split("/")[0]

    ## get the failed jobs
    os.system("cat statuslog | grep \"failed\" > tmp")
    if not os.stat('tmp').st_size==0:
        cols=(open("tmp")).readline().rstrip().split(" ")[-1]
        failedjobs=cols.split("/")[0]


    ## get the running jobs
    os.system("cat statuslog | grep \"running\" > tmp")
    if not os.stat('tmp').st_size==0:
        cols=(open("tmp")).readline().rstrip().split(" ")[-1]
        runningjobs=cols.split("/")[0]


    ## get the idle jobs
    os.system("cat statuslog | grep \"idle\" > tmp")
    if not os.stat('tmp').st_size==0:
        cols=(open("tmp")).readline().rstrip().split(" ")[-1]
        idlejobs=cols.split("/")[0]

    ## submitted jobs
    os.system("cat statuslog | grep \"submittedjobs\" > tmp")
    if not os.stat('tmp').st_size==0:
        cols=(open("tmp")).readline().rstrip().split(" ")[-1]
        submittedjobsjobs=cols.split("/")[0]



    #print "total     finished     running     failed     idle    submitted     transferring     crabdir "
    #print cleannumber(totaljobs),"       ", cleannumber(finishedjobs), "         ", cleannumber(runningjobs), "      ",cleannumber(failedjobs), "        ", cleannumber(idlejobs),\
    #    "         ", cleannumber(submittedjobs), "            ", cleannumber(transferringjobs), "      ", firstline

    towrite= (cleannumber(totaljobs)+"       "+ cleannumber(finishedjobs)+ "     "+ cleannumber(runningjobs) + "      "+ cleannumber(failedjobs)+ "        "+ cleannumber(idlejobs)+\
        "         "+ cleannumber(submittedjobs)+ "            "+ cleannumber(transferringjobs)+ "      "+ firstline+"\n")
    #print towrite

    fout=open("crabtasksummary.rkl","a")
    fout.write(towrite)
    fout.close()


def status(crabdirname):
    print "cehcking the status of all jobs "
    os.system("ls -1 "+crabdirname + " >& crabtasklist.txt")
    os.system("rm crabtasksummary.rkl")

    for itask in open('crabtasklist.txt'):
        print "checking status of ",itask.rstrip()
        if not args.ss:
            os.system('crab status -d '+crabdirname+'/'+itask.rstrip())
        if args.ss:
            os.system ('crab status -d '+crabdirname+'/'+itask.rstrip() + ' | tee statuslog')
            statussummary("statuslog")




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



textfilename='signal_2018_ZpBaryonic.txt'


if args.submit:
    f = open(textfilename,'r')
    for line in f:
        print line
        a,b,c = line.split()
        datasetdetail=[a,b,c]
        submit(datasetdetail)

if (args.status or args.resubmit or args.kill ) and len(sys.argv)<3 :
    print "insufficienct input provided, please look at the script MultiCrab_2018MC.py for usage. or ask Raman"

if (args.status or args.resubmit or args.kill or args.ss ) and len(sys.argv)>=3 :
    print "two argument provided"
    if args.status:
        status(args.crabdir)
        print "------------- ------------- ------------ ----------- ---------------- --------------"
        print "here is the summary of all the jobs (this summary feature is still in testing mode"
        print "total     finished     running     failed     idle    submitted     transferring     crabdir "
        os.system('cat crabtasksummary.rkl')
        print "------------- ------------- ------------ ----------- ---------------- --------------"
    if args.resubmit: resubmit(args.crabdir)

    if args.kill: kill(args.crabdir)

    #if args.ss: statussummary(args.crabdir)
