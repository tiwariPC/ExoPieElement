import os, datetime
datestr = datetime.date.today().strftime("%Y%m%d")

f = open('data_2017.txt','r')
for line in f:

	a=line.split()[0]
	b=line.split()[1]
	c=line.split()[2]
	d=line.split()[3]

	workname='2017_data_'+datestr

	if len(a)>100:
		reqname=a[0:100]
	else:
		reqname=a
	dataset=c


        if '2016B' in c:
            period='BCD'
        elif '2016C' in c:
                period='BCD'
        elif '2016D' in c:
            period='BCD'
        elif '2016E' in c:
            period='EF'
        elif '2016F' in c:
            period='EF'
        elif '2016G' in c:
            period='G'
        else:
            period='H'



	tempfile = open('crabConfigTemp.py','w')
	cfile = open('crabConfig_data.py','r')
	for cline in cfile:
		if cline.startswith("workname="):
			tempfile.write("workname=\'"+workname+"\'\n")
		elif cline.startswith("reqname="):
			tempfile.write("reqname=\'"+reqname+"\'\n")
		elif cline.startswith("dataset="):
			tempfile.write("dataset=\'"+dataset+"\'\n")
		elif cline.startswith("PERIOD="):
			tempfile.write("PERIOD=\'"+period+"\'\n")
		else:
			tempfile.write(cline)
	cfile.close()
	tempfile.close()

	print "\n==========================\nSubmitting "+a+"\n==========================\n"
	os.system("crab submit -c crabConfigTemp.py")
        #print a, "  ", b, "  ", c , "  ", d

f.close()

