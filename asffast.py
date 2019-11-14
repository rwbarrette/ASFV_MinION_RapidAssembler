
import time
import os
import shutil
import sys, getopt
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from time import mktime
from datetime import datetime
from stat import S_ISREG, ST_CTIME, ST_MODE


def main(argv):
   #try: (SourceLOC, BCtarget)
   opts, args = getopt.getopt(argv,"hf:b:",["inFASTQ=","inBC="])

   for opt, arg in opts:
      if opt == '-h':
         #print 'test.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-f", "--FQ"):
         FASTQinFOLDER = arg
      elif opt in ("-b", "--BC"):
         BCtarget = arg

   return FASTQinFOLDER, BCtarget

def RunPORECHOP(FastqInFile, DemultiplexedOUT_Folder):
    import os
    import sys
    print("Running Porechop Demultiplexer...")

    PORECHOPline = "porechop -i %s -b %s" % (FastqInFile, DemultiplexedOUT_Folder)
    status = os.system(PORECHOPline)


def RunQCAT(FastqInFileFolder, DemultiplexedOUT_Folder):
    import os
    import sys
    print("Running QCAT Demultiplexer...")
    

    PORECHOPline = "cat %s/*.fastq | qcat -b %s" % (FastqInFileFolder, DemultiplexedOUT_Folder)
    status = os.system(PORECHOPline)


def RunBLAST_IMAGER_PL(OutPNG):
	BLAST_IMAGERline = "perl /blast-imager.pl /REPORT/VIRCurrentTorrent_DIR.TBL > %s"%OutPNG
	status = os.system(BLAST_IMAGERline)


def Parse_STAT_File(OutSubFolderLOC, TimeSTAMP, SEQ_PROFILE_ALL):
	
	BC_count=0
	listingSTAT = os.listdir(OutSubFolderLOC)
	for infileSTAT in listingSTAT:
		isStat = infileSTAT.find("_STAT")
		if isStat>=0:
				if infileSTAT==".DS_Store":
					print("Ignoring")
				else:
					BC_count=BC_count+1
					QUEpct_ID=0
					print(infileSTAT)

					SeqProfileRECORD=[]

					BCparse = infileSTAT[0:12].rfind(".")
					BC = infileSTAT[0:BCparse].replace("barcode","") #Sample ID
					print("BARCODE_ID:::", BC)
					if BC=="none":
						BC=0
					else:
						if BC=="":
							BC=0
						else:	
							BC=int(BC)
					
					print("Barcode", BC)
					print(isStat,BCparse)
					
					infileSTATloc = OutSubFolderLOC+infileSTAT
					StatFile = open(infileSTATloc, "r")
					for xS in StatFile:
						#print xS
						isTOTAL = xS.find("in total (")
						isMAPPED = xS.find("mapped (")
						isTOTALNT = xS.find("### Total NT:")
						isALN_NT = xS.find("### ALN NT:")
						isPCT_ID = xS.find("### Percent Coverage(%):")

						#QUEtotal=""
						#QUEmapped=""
						#QUEtotalnt=""
						#QUEaln_nt=""
						#QUEpct_ID=""

						if isTOTAL>=0:
							TOTALend = xS.find("+")
							QUEtotal =xS[0:TOTALend-1]
							QUEtotal = QUEtotal.replace("\n","")
							print(QUEtotal)
						if isMAPPED>=0:
							MAPPEDend = xS.find("+")
							QUEmapped = xS[0:MAPPEDend-1]
							QUEmapped = QUEmapped.replace("\n","")
							print(QUEmapped)
						if isTOTALNT>=0:
							QUEtotalnt = xS[14:]
							print(QUEtotalnt)
						if isALN_NT>=0:
							QUEaln_nt = xS[12:]
							print(QUEaln_nt)
						if isPCT_ID>=0:
							QUEpct_ID = xS[25:]
							print(QUEpct_ID)
							QUEpct_ID = QUEpct_ID.replace("\n","")
					
				SeqProfileRECORD.append(BC)
				SeqProfileRECORD.append(TimeSTAMP)
				SeqProfileRECORD.append(QUEpct_ID)
				SeqProfileRECORD.append(QUEtotal)
				SeqProfileRECORD.append(QUEmapped)

				#SeqProfileRECORD.update({"TS":TimeSTAMP})
				##SeqProfileRECORD.update({"BC":BC})
				#SeqProfileRECORD.update({"PCTID":QUEpct_ID})
				#SeqProfileRECORD.update({"TOTREAD":QUEtotal})
				#SeqProfileRECORD.update({"SPECREAD":QUEmapped})
				#print(SeqProfileRECORD)

				SEQ_PROFILE_ALL.append(SeqProfileRECORD)
					
					

	return SEQ_PROFILE_ALL, BC_count


if __name__ == "__main__":
    SourceLOC, BCtarget = main(sys.argv[1:])   
    INIT_Watcher(SourceLOC, BCtarget)

def INIT_Watcher(SourceLOC, BCtarget):    
	FirstPass = 1
	BC_count = 0
	SEQ_PROFILE_ALL = []
	BC_list = []
	BC_used = []
	BC_TARGET = BCtarget
	RootDataFolder = ""
	ReferenceGenome = "/INPUT/ASF_Georgia.fasta"
	SourceFOLDER = SourceLOC #"/INPUT/FASTQ_IN/04232019Killean/"
	QueFOLDER = "%sfastqQUE/"%(SourceFOLDER)
	DemultiplexedFOLDER = "%sfastqDMTPLX/"% (SourceFOLDER)
	DoneFOLDER = "%sfastqDONE/"%(SourceFOLDER)
	ResultsFOLDER = "%sfastqRESULTS/"%(SourceFOLDER)

	shutil.copyfile("/ReportINIT.txt", "/REPORT/FastqQUE_TEST.txt")

	print("Removing Old Directories.")
	try:
		shutil.rmtree(QueFOLDER)
	except:
		print("no QueFolder.")
	try:	
		shutil.rmtree(DemultiplexedFOLDER)
	except:
		print("no DMLTPLX_Folder.")
	try:
		shutil.rmtree(DoneFOLDER)
	except:
		print("no DoneFolder.")
	try:
		shutil.rmtree(ResultsFOLDER)
	except:
		print("no ResultsFolder.")
		
	print("Creating New Directories.")

	os.mkdir(QueFOLDER)
	os.mkdir(DemultiplexedFOLDER)
	os.mkdir(DoneFOLDER)
	os.mkdir(ResultsFOLDER)


	CSVout = "%s/SEQReadData_Profile.csv"%(ResultsFOLDER)

	f = open(CSVout, "w")
	writer = csv.writer(f)
	writer.writerow( ('Date','Time', 'Percent_Coverage', 'ASF_Reads', 'Non_ASF_Reads', 'Tm') )

	print(time.time())
	starttime = time.time()

	t_end = time.time() + 31
	CheckTime = time.time() + 15
	print(starttime, t_end)


	firstRound = True
	firstRange = True

	while time.time() < t_end:


		# Get files in SOURCE directory

		if time.time() > CheckTime:
			print("Check Sequence FOLDER NOW!!")
			listing = os.listdir(SourceFOLDER)

			# Check for new file
				#No file in directory

			if len(listing)==0:
				print("No new sequencing file yet.")

				#New file in directory
			else:

				# path to the directory (relative or absolute)

				dirpath = SourceFOLDER

				# get all entries in the directory w/ stats
				entries = (os.path.join(dirpath, fn) for fn in os.listdir(dirpath))
				entries = ((os.stat(path), path) for path in entries)

				entries = ((stat[ST_CTIME], path)
						for stat, path in entries if S_ISREG(stat[ST_MODE]))


				for cdate, path in sorted(entries):
					print(time.ctime(cdate))
				#for infile in listing:
					infile = os.path.basename(path)

					if infile==".DS_Store":
						print("Ignoring")
					else:
						try:
							shutil.move(SourceFOLDER+infile, QueFOLDER)
						except:
							print("Folder exists...")

						
						print(infile)
						FileCheckTime = time.ctime(cdate) #time.gmtime()
						OutSubFolderNAME = infile.replace('.fastq','/')
						OutSubFolderLOC = DemultiplexedFOLDER+OutSubFolderNAME
						
						try:
							shutil.rmtree(OutSubFolderLOC)
							print("Removing Old Directory.")
						except:
							print("Creating New Directory.")
						print(FileCheckTime)
						DateFound = datetime.strptime(FileCheckTime, "%a %b %d %H:%M:%S %Y")
						#TimeFound = datetime.strptime(FileCheckTime, "%H_%M:%S")
						TimeSTAMP = DateFound
						#timestamp = DateFound
						print(TimeSTAMP)

						#DateFound = time.strftime("%Y-%m-%d", FileCheckTime)
						#TimeSTAMP = time.strftime("%m%d%y_%H_%M_%S")
						#timestamp = time.strftime("%y-%m-%d %H:%M:%S")
		
						dt = datetime.fromtimestamp(mktime(DateFound.timetuple()))

						print(dt.hour, dt.minute, dt.day)
						dDays = int(dt.day)


						if firstRound == True:
							startDay = dDays
							startTime = (int(dt.hour)*60)+int(dt.minute)
						print(startTime)
						print(startDay)
						print(dDays)
						dElapsedTime = (((startDay-dDays)*24)+(int(dt.hour)*60)+int(dt.minute))-int(startTime)


						print(dElapsedTime, "<<<<<<<<<<<<<<<<(((((((((<<<<<<<<<<<<<< Timestamp")  # THIS IS THE STARTTIME IN MINUTES
						
	
						RunQCAT(QueFOLDER, OutSubFolderLOC)
						print(QueFOLDER, OutSubFolderLOC)
						

						SystemLine = ("python ./Multi_Reference_Guided_Assembly_DOCKERb.py -f %s -r %s -t %s")% (OutSubFolderLOC, ReferenceGenome, dElapsedTime)
							
						print(SystemLine)
						os.system(SystemLine)
						
						#RunBLAST_IMAGER_PL(PNGoutRESULT)

						SEQ_PROFILE_ALL, BC_count = Parse_STAT_File(OutSubFolderLOC, dElapsedTime, SEQ_PROFILE_ALL)
						print(SEQ_PROFILE_ALL)

						df = pd.DataFrame(SEQ_PROFILE_ALL, columns = ['barcode','timestamp','pct_coverage','total_reads','specific_reads'])
						print(df)
						df['pct_coverage'] =df['pct_coverage'].astype(float)
						df['barcode'] =df['barcode'].astype(np.int64)
						#df['timestamp']=df['timestamp'].astype('daetime64[ns]')
						#df['timestamp'] = [time.date() for time in df['timestamp']]
					if 1==1:#quit
					###try:
						for bcX in range(0,97):
							print(bcX)
							if bcX==96:
								print("END OF RANGE 1--------------------------")
								firstRange==False
								fig.savefig(OutSubFolderLOC_PNG)
							print("*************")
							if bcX == BC_TARGET:
								Sampling = (df.loc[df['barcode']==bcX])
								if Sampling.empty!=True:

									OutSubFolderLOC_PNG = DemultiplexedFOLDER+"_"+str(dElapsedTime)+"_BC"+str(bcX)+".png"
									OutSubFolderLOC_XLS = DemultiplexedFOLDER+"_"+str(dElapsedTime)+"_BC"+str(bcX)+".xlsx"
									print("*************++")
									print(OutSubFolderLOC_XLS)

									#Sampling['timestamp'] = pd.to_datetime(Sampling['timestamp'])
									Sampling['timestamp'] =Sampling['timestamp'].astype(np.int64)
									print(Sampling)

									Sampling.to_excel(OutSubFolderLOC_XLS)
									#Sampling.plot(kind='scatter',x='timestamp',y='pct_coverage', color='red')
												
									
									fig=plt.figure(1)
									ax=fig.add_subplot(111)

									BarCODE = "BC%s"%str(bcX)
									if BarCODE in BC_used:
										BarCODE = ""
									else:
										BC_used.append(BarCODE)

									#ax.plot(Sampling['timestamp'],Sampling['pct_coverage'], label=BarCODE)

									if firstRound==True:
										print("Current Barcode", str(BarCODE))
										box = ax.get_position()
										ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
										firstRound=False
										print(Sampling['timestamp'])
										print(Sampling['pct_coverage'])
									#if firstRange==True:	
									ax.plot(Sampling['timestamp'],Sampling['pct_coverage'], label=BarCODE)
									ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
										
									print(Sampling['timestamp'])
									print(Sampling['pct_coverage'])
									
									#ax.annotate(s=('blah'),xy=(Sampling['timestamp'],Sampling['pct_coverage']))

									plt.xlabel("Time of Day")
									plt.ylabel("% Coverage")
									plt.title("Genome Sequence Coverage of ASF virus")
									#plt.legend()
									
									#plt.set_ylim(bottom=0, top=100)
									plt.ylim((0,100))
									plt.xticks(rotation='vertical')


						#Need to write a stats parser!! #writer.writerow( ('Date','Time', 'Percent_Coverage', 'ASF_Reads', 'Non_ASF_Reads', 'Tm') )

		          		#BLAST TAB OUT::: "/users/rwbarrettemac/bioinformatics/pythonfolders/TorrentFiles/VIRCurrentTorrent_DIR.TBL"
		          		#perl /Users/rwbarrettemac/Desktop/2015_Py_Scripts/blast-imager.pl /Users/rwbarrettemac/Bioinformatics/PythonFolders/TorrentFiles/VIRCurrentTorrent_DIR.TBL > /Users/rwbarrettemac/Desktop/ASF_WGS/TBLtestOUT.png"

				# Set next file check time
			CheckTime = CheckTime + 15

	