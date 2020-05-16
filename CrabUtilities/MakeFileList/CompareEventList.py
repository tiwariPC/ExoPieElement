import os 
import sys


store_path1="/hdfs/store/user/khurana/ExoPieElement/setup_2017_2016_v03"
store_path2="/hdfs/store/user/khurana/ExoPieElement/setup_2017_2016_v03_ite2Automatic/"

a=[]
file1=sys.argv[1]
file2=sys.argv[2]

for line in open(file1):
    a.append(line)

b=[]
for line2 in open(file2):
    #bb = line.split()
    b.append(line2)


aminusb=set(a)-set(b)
bminusa=set(b)-set(a)
sameevent=set(a)&set(b)

ab=list(aminusb)
ba=list(bminusa)
same=list(sameevent)

print "-----------------------------------------------------------------------------------------------------------------------------------------------"
to_print_same="list of dataset which are common in {} and {}"
print to_print_same.format(file1,file2)
for iele in same: 
    datasetname =  (iele.rstrip().split(".txt")[0]).split("crab_")[1].split("_00")[0]
    print "ls -ltr "+datasetname+" ; du -sh "+datasetname
    
    

    #print "---------------------------------------------------size in each crab dir--------------------------------------------------------------------------------------------"
    #os.system("du -sh "+store_path1+"/"+datasetname)
    #os.system("du -sh "+store_path2+"/"+datasetname)
    #to_Print_size= "size for {} is {} and size for {} is {}"
    #print to_Print_size.format(store_path1,os.system("du -sh "+store_path1+"/"+datasetname), store_path2, os.system("du -sh "+store_path2+"/"+datasetname) )

print "-----------------------------------------------------------------------------------------------------------------------------------------------"
to_print_diff="list of dataset which are present in {} but not in {}"
print to_print_diff.format(file1,file2)
for jevent in ab:
    print jevent.rstrip()

print "-----------------------------------------------------------------------------------------------------------------------------------------------"
print to_print_diff.format(file2,file1)
for ievent in ba:
    print ievent.rstrip()


print "-----------------------------------------------------------------------------------------------------------------------------------------------"
print "lenght of ab =", len(ab)
print "lenght of ba =", len(ba)
print "length of sameevent =", len(same)
