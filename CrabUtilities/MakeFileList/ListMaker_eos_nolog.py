import os,sys

errmsg="Usage:\n$ python ListMaker.py crab T2path outfileprefix\nor\n$ python ListMaker.py st T2path outfileprefix"

if len(sys.argv)==4:
    if sys.argv[1]=="crab":
        isCrab=True
    elif sys.argv[1]=="st":
        isCrab=False
    else:
        print errmsg
        sys.exit()
    T2path=sys.argv[2]
    filepref=sys.argv[3]
else:
    print errmsg
    sys.exit()

inpfilename='T2List_'+filepref+'.txt'
os.system('ls -R '+"/hdfs/"+T2path+' | cat &> ' +inpfilename)

f=open(inpfilename,"r")

os.system("mkdir -p Filelist"+"_"+filepref)
pref="root://cmsxrootd.hep.wisc.edu/"
filecount=1
lineminus1=""
lineminus2=""
fileopen=False
failedfile=False
#log=open("log_"+filepref+".txt","w")

for line in f:
    if not line=="\n":
        fname=line.split()[-1]
    else:
        fname=""

    if fname.endswith(".root") and not fileopen:
        folder=pref+lineminus2[:-2]+"/"
        
        folder=folder.replace("/hdfs/","/")
        if 'failed' in lineminus2 or lineminus2.split("/")[-1].strip()=="failed:" or lineminus2.split("/")[-1].strip()=="log:": failedfile=True

        if not failedfile:
            if isCrab:
                #print ('lineminus2: ',lineminus2)
                realname=lineminus2.split("/")[-3]+"_"+lineminus2.split("/")[-1][:-2]
            else:
                #print (lineminus2)
                realname=lineminus2.split("/")[-1][:-2]
            out=open("Filelist"+"_"+filepref+"/"+realname+".txt","w")
            out.write(folder+fname+"\n")
            filecount+=1
        fileopen=True
    elif fname.endswith(".root"):
        if not failedfile: out.write(folder+fname+"\n")
    elif fileopen:
        if not failedfile: out.close()
        fileopen=False
        failedfile=False

    #print ('end_lineminus2',lineminus2)
    if lineminus1=="\n":
        lineminus2=line
    else:
        lineminus2 = lineminus1
    #print ('end_linminus1',lineminus1)
    lineminus1=line
    #print ('end_line', line)

#log.close()
f.close()

print ("Created Filelist_%s directory and log_%s.txt." %(filepref,filepref))
