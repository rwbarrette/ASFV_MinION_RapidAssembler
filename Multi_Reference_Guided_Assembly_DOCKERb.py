
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import os
import sys
import string
import sys, getopt

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
from shutil import copy
import shutil

def main(argv):
   #try:
   opts, args = getopt.getopt(argv,"hf:r:t:",["inFASTQ=","inREF=","timeSTAMP="])

   for opt, arg in opts:
      if opt == '-h':
         #print 'test.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-f", "--inFASTQ"):
         FASTQinFOLDER = arg
      elif opt in ("-r", "--inREF"):
         FASTAinREF = arg
      elif opt in ("-t", "--tStamp"):
         timeSTAMP = arg

   return FASTQinFOLDER, FASTAinREF, timeSTAMP

def RunBWA_SE(REFInfile, FastqInFileFWD, FastqInFileREV):
    print("Running BWA Indexing and SAI file build...")
    
    
    BWA_index = "/BWA/bwa index %s" % (REFInfile)
    status = os.system(BWA_index)
    print(BWA_index)
    
    BWA_saiF = "/BWA/bwa mem %s %s >  /BWAfiles/BWAsam.sam" % (REFInfile, FastqInFileFWD)
    print(BWA_saiF)
    
    status = os.system(BWA_saiF)

def RunSAM_SortedBam(REFInfile):
    print("Running SAMTOOLS Index Build (faidx)...")
    
    RunSamIDX = "/samtools/samtools faidx %s" % (REFInfile)
    #print RunSamIDX, "1"
    status = os.system(RunSamIDX)
  
    RunSam_Import = "/samtools/samtools import %s.fai /BWAfiles/BWAsam.sam /BWAfiles/BWAbam.bam" %(REFInfile)
    #print RunSam_Import, "2"
    status = os.system(RunSam_Import)
    


    RunSamSORT = "/samtools/samtools sort /BWAfiles/BWAbam.bam -o /BWAfiles/BWAbam.sorted.bam"
    #print RunSamSORT, "3"
    status = os.system(RunSamSORT)
    
    print("SAM_DEPTH----start")
    RunSam_Depth = """/samtools/samtools depth  /BWAfiles/BWAbam.sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'"""
    RunSam_DepthB = """/samtools/samtools depth  /BWAfiles/BWAbam.sorted.bam > /BWAfiles/FastqQUE_DEPTH.txt""" 
    RunSam_DepthC = """/samtools/samtools depth  /BWAfiles/BWAbam.sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > /BWAfiles/FastqQUE_DEPTHavg.txt"""
     #print RunSam_Import, "2"
    status = os.system(RunSam_Depth)
    status = os.system(RunSam_DepthB)
    status = os.system(RunSam_DepthC)
    print("SAM_DEPTH----end")

    df_COVERAGE = pd.read_csv("/BWAfiles/FastqQUE_DEPTH.txt", sep='\t')


    RunSamIndexBam = "/samtools/samtools index /BWAfiles/BWAbam.sorted.bam"
    #print RunSamIndexBam, "4"
    status = os.system(RunSamIndexBam)

def RunBAM_STATS(OutfileSTAT):
    print("Running SAMTOOLS Flagstat for read statistics...")
    
    RunSTATS_CMDLINE = "/samtools/samtools flagstat /BWAfiles/BWAbam.bam > %s" %(OutfileSTAT)
    
    status = os.system(RunSTATS_CMDLINE)

    # Duplicate Stats file to central LOC
    print(OutfileSTAT.find("none"))
    print("*****&&&&&&&&*******&&&&&&&*****")
    if OutfileSTAT.find("none")>0:
      print("unbarcoded stats...N/A")
    else:
      shutil.copy(OutfileSTAT, "/BWAfiles/CurrentSTAT.txt")


def RunSAM_Consensus(REFInfile, Outfile):

    print("Running SAMTOOLS Consensus Build (mpileup)...")
    
    RunSam_Cons_CMDLINE = "/samtools/samtools mpileup -uf %s /BWAfiles/BWAbam.sorted.bam | /bcftools/bcftools call -c | /bcftools/misc/vcfutils.pl vcf2fq > %s" % (REFInfile, Outfile)
    
    print(RunSam_Cons_CMDLINE)
    
    status = os.system(RunSam_Cons_CMDLINE)

def Convert_CONSENS_to_FASTA(Outfile, NewOutfile, OutfileSTAT, timeSTAMP, RefLEN):
    StatOUTfile = open(OutfileSTAT, "a")
    ModOUTFILE = open(NewOutfile, "w")
    ModINFILE = open(Outfile, "r")
    BLASTout = open("/BWAfiles/CurrentSTAT.txt","a")

    mLineTotal = 0
    mLineNonN = 0
    inlineCNT = 0
    prcntid = 0
    for mLine in ModINFILE:
        inlineCNT=inlineCNT+1
        if "+" in mLine:
            "DONE"
            break
        else:
            if "@" in mLine:
                StartString = ">"+mLine
                ModOUTFILE.write(StartString)
            
            else:
                #print(mLine)
                mLineLEN = len(mLine)-1
                mLineNcnt = mLine.count("n")
                #print(mLineLEN, mLineNcnt, "**********************###############")
                ModOUTFILE.write(mLine)
                mLineTotal = mLineTotal+(mLineLEN)
                mLineNonN = mLineNonN+(mLineLEN-mLineNcnt)
                #print(mLineTotal, mLineNonN)
                prcntid = (float(mLineNonN/float(RefLEN)))*100
    
    StatOUTfile.write("\n")
    StatOUTfile.write("### Total NT: " + str(RefLEN))
    StatOUTfile.write("\n")
    StatOUTfile.write("### Region NT: " + str(mLineTotal))
    StatOUTfile.write("\n")
    StatOUTfile.write("### ALN NT: " + str(mLineNonN))
    StatOUTfile.write("\n")
    StatOUTfile.write("### Percent Coverage(%): "+ str(prcntid))
    StatOUTfile.write("\n")

    BLASTout.write("\n")
    BLASTout.write("### Total NT: " + str(RefLEN))
    BLASTout.write("\n")
    BLASTout.write("### Region NT: " + str(mLineTotal))
    BLASTout.write("\n")
    BLASTout.write("### ALN NT: " + str(mLineNonN))
    BLASTout.write("\n")
    BLASTout.write("### Percent Coverage(%): "+ str(prcntid))
    BLASTout.write("\n")

    ModOUTFILE.close()
    ModINFILE.close()
    StatOUTfile.close()
    BLASTout.close()

    return mLineTotal, mLineNonN     

def ReadFastQinFile(infile_FASTq, outfile_FASTa, New): #Sets up FASTA file from FASTQ
                           
    ID = 0
    seqcount= 0
    seqcountMAX = 100000000
    handleFWD = open(infile_FASTq)
                           
    fq_dictFWD = SeqIO.parse(handleFWD,"fastq")
    if New == "TRUE":

        fastq_file = open(outfile_FASTa, "w")
    else:
        fastq_file = open(outfile_FASTa, "a")

    print("Reading Fastq...")
    
    for seqFWD in fq_dictFWD:
                        if seqcount <=seqcountMAX:   
                               seqcount = seqcount + 1
                               
                               QscoreMin = min(seqFWD.letter_annotations["phred_quality"])
                               #print QscoreMin
                               ID = ID + 1
                               
                               fastq_sequence = seqFWD.seq
                               SeqTitle = seqFWD.description
                               
                               if len(fastq_sequence) >= 0:
                                   
                                   fastq_IDstring = ">"+ SeqTitle + "\n"
                                   fastq_file.write(fastq_IDstring)
                                   
                                   for i in range(0, len(fastq_sequence), 72):
                                       fastq_file.write(fastq_sequence[i:i+72])
                                       fastq_file.write("\n")
    fastq_file.close()
    handleFWD.close()


def RunGsnap(REFInfile, FastqInFileFWD, FastqInFileREV):

    print("Running GSnap Indexing and SAI file build...")
    

    GSnap_index = "gmap_build -d ViralRefGMp %s" % (REFInfile)
    print(GSnap_index)
    status = os.system(GSnap_index)

    GSnap_saiF = "gsnap -d ViralRefGMp %s %s -A sam > /BWAfiles/BWAsam.sam" % (FastqInFileFWD, FastqInFileREV)

    status = os.system(GSnap_saiF)

def BLASTN_v29(outfile, infile, dbfile):

    print("BLASTn (Sequence Read Match)")
    BLASTN_CMDLINE = "/NCBI/bin/blastn -out %s -query %s -db %s -outfmt 5 -num_alignments 1000000"# -num_alignments -task blastn-short 100000" #took out -task blastn-short as it was too slow and didn't increase alignments
    
    status = os.system(BLASTN_CMDLINE % (outfile, infile, dbfile))

def BLASTN_v29_TabularOUT(outfile, infile, dbfile):

    print("BLASTn (Sequence Read Match)")
    BLASTN_CMDLINE = "/NCBI/bin/blastn -out %s -query %s -db %s -outfmt 6 -num_alignments 1000000"# -num_alignments -task blastn-short 100000" #took out -task blastn-short as it was too slow and didn't increase alignments
    
    status = os.system(BLASTN_CMDLINE % (outfile, infile, dbfile))

def BlastSequencesV29_TAB(TemplateInFile, BlastType):

    from Bio.Blast.Applications import NcbiblastnCommandline
     
    BlastNThandle = open(TemplateInFile)
    BlastNTread = SeqIO.parse(BlastNThandle,"fasta")
    
    IDcnt = 0
    INDEX = 0
    preIndexCNT = 0
    DBblast = "/BWAfiles/RawSequence_DIRECT"
    #SaveTitle = "c:\\TorrentFiles\\FastaQue.fasta"
    print("Running BLAST...")
    
    BLASTN_v29_TabularOUT("/BWAfiles/VIRCurrentTorrent_DIR.TBL", TemplateInFile, DBblast)
    
def RunBLAST_IMAGER_PL(OutPNG):
    BLAST_IMAGERline = "perl /blast-imager.pl /BWAfiles/VIRCurrentTorrent_DIR.TBL > %s"%OutPNG
    status = os.system(BLAST_IMAGERline)


def BlastSequencesV29(TemplateInFile, BlastType, OutFileStat):

    from Bio.Blast.Applications import NcbiblastnCommandline
    
    BlastNThandle = open(TemplateInFile)
    BlastNTread = SeqIO.parse(BlastNThandle,"fasta")
    for reads in BlastNTread:
        TestSeq = reads.seq
    
    QueLen = len(TestSeq)
    print("Template Length: ", QueLen)
    
    IDcnt = 0
    INDEX = 0
    preIndexCNT = 0
    DBblast = "/BWAfiles/RawSequence_DIRECT"

    print("Running BLAST...")
    
    BLASTN_v29_TabularOUT("/BWAfiles/VIRCurrentTorrent_DIR.TBL", TemplateInFile, DBblast)
    BLASTN_v29("/BWAfiles/VIRCurrentTorrent_DIR.xml", TemplateInFile, DBblast)
    
    fastq_file = open(OutFileStat, "a")

    print("BLAST complete.")
    
                           
    from Bio.Blast import NCBIXML
    result_handleNEW = open("/BWAfiles/VIRCurrentTorrent_DIR.xml")
    blast_records = NCBIXML.parse(result_handleNEW)
    for blast_record in blast_records:
            
            
                               for alignment in blast_record.alignments:
                                   preIndexCNT = preIndexCNT+1
                #if IDcnt <1:
                                   for hsp in alignment.hsps:
                                       if hsp.identities >=0:
                                        
                                        IDpep = float(hsp.identities)
                                        Lenpep = float(len(hsp.query))
                                        if Lenpep>=0:
                                           prcntid = IDpep/Lenpep
                                              
                                           PreGI = alignment.title
                                           print(PreGI, ".......",Lenpep, "nt")
                                           QueSEQ = hsp.query  # FASTA in file
                                           MatchSEQ = hsp.match
                                           SubjSEQ = hsp.sbjct  # Illumina record
                                           MatchCNT = MatchSEQ.count("|")
                                           #print SubjSEQ

                                           
                                           StartPos = hsp.query_start
                                           EndPos = hsp.query_end
                        
                                           GIstart = 3
                                           PreGIend = PreGI[GIstart:15]
                                           GIend = PreGIend.find("|")+ 3
                                           GIid = PreGI[GIstart:GIend]
                                           INDEX = INDEX + 1

                                           fastq_file.write("\n")
                                           fastq_file.write(PreGI)
                                           fastq_file.write("### Percent ID(%): "+ str(prcntid))
                                           fastq_file.write("\n")
        
    fastq_file.close()
    print("preIndex Count",preIndexCNT)
    print("index",INDEX)

def Parse_STAT_FileBLAST(OutSubFolderLOC, TimeSTAMP, SEQ_PROFILE_ALL):
  
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

          BCparse = infileSTAT[0:6].rfind(".")
          BC = infileSTAT[0:BCparse].replace("BC","") #Sample ID
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
              QUEpct_ID = xS[19:]
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
          
          
  #StatIN = open(OutSubFolderLOC, "r")

  #for sIN in StatIN:
#   print sIN
  return SEQ_PROFILE_ALL, BC_count

def BlastSequencesV29_Xref(TemplateInFile, BlastType, OutFileBLAST, StatsFile):

    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML
    
    
    BlastNThandle = open(TemplateInFile)
    BlastNTread = SeqIO.parse(BlastNThandle,"fasta")
    for reads in BlastNTread:
        TestSeq = reads.seq
    
    QueLen = len(TestSeq)
    print("Template Length: ", QueLen)
    
    IDcnt = 0
    INDEX = 0
    preIndexCNT = 0
    HSPcount = 0
    HSPstart=0
    HSPpreviousstart=1

    DBblast = "/INPUT/ASFV_p72all"

    print("Running BLAST xREF...")
    
    BLASTN_v29("/BWAfiles/VIRCurrentTorrent_DIRBLAST.xml", TemplateInFile, DBblast)
    
    fastq_file = open(OutFileBLAST, "w")

    print("BLAST complete.")
    
                           
    fastq_file.write("\n")
    fastq_file.write(str("-----------------------------------------------\n"))
    fastq_file.write(str("ASF-FAST Real-Time Sequencing Report Results\n"))
    fastq_file.write(str("-----------------------------------------------\n"))
    fastq_file.write("\n")
    fastq_file.write(str("[ALIGNMENT RESULTS (TOP 5) OF ASF GENOME CROSS-REFERENCE] \n"))


    result_handleNEW = open("/BWAfiles/VIRCurrentTorrent_DIRBLAST.xml")
    blast_records = NCBIXML.parse(result_handleNEW)
    for blast_record in blast_records:
            
                               
                               for alignment in blast_record.alignments:
                                  if IDcnt <1:
                                       for hsp in alignment.hsps:
                                          HSPstart = alignment.title[0:10]
                                          if HSPstart!=HSPpreviousstart:
                                             


                                             if hsp.identities >=0:
                                              
                                              IDpep = float(hsp.identities)
                                              Lenpep = float(len(hsp.query))
                                              if Lenpep>=0:
                                                 HSPpreviousstart=HSPstart
                                                 prcntid = (IDpep/Lenpep)*100
                                                    
                                                 PreGI = alignment.title
                                                 print(PreGI, ".......",Lenpep, "nt")
                                                 QueSEQ = hsp.query  # FASTA in file
                                                 MatchSEQ = hsp.match
                                                 SubjSEQ = hsp.sbjct  # Illumina record
                                                 MatchCNT = MatchSEQ.count("|")
                                                 #print SubjSEQ

                                                 
                                                 StartPos = hsp.query_start
                                                 EndPos = hsp.query_end
                              
                                                 GIstart = 3
                                                 PreGIend = PreGI[GIstart:15]
                                                 GIend = PreGIend.find("|")+ 3
                                                 GIid = PreGI[GIstart:GIend]
                                                 INDEX = INDEX + 1

                                                 fastq_file.write("\n")
                                                 fastq_file.write(str(PreGI))
                                                 fastq_file.write("\n")
                                                 fastq_file.write("### Percent ID(%): "+ str(prcntid))
                                                 fastq_file.write("\n\n")
                                                 fastq_file.write(str("%s nucleotides matched")%IDpep)
                                                 fastq_file.write("\n")
                                                 fastq_file.write(str("%s total nucleotides resolved")%Lenpep)
                                                 fastq_file.write("\n")
                                                 HSPcount=HSPcount+1
                                             else:
                                              continue
                                  IDcnt = IDcnt+1

    fastq_file.write("\n")
    fastq_file.write(str("-----------------------------------------------\n"))
    fastq_file.write(str("ASF-FASTA SEQUENCE ASSEMBLY STATISTICS\n"))
    fastq_file.write(str("-----------------------------------------------\n"))
    fastq_file.write("\n")


    STAT_inFILE = open("/BWAfiles/CurrentSTAT.txt", "r")
    for xS in STAT_inFILE:
 
      fastq_file.write(str(xS))

    STAT_inFILE.close()

    fastq_file.close()
    print ("preIndex Count",preIndexCNT)
    print ("index",INDEX)


       

def RunGsnapSE(REFInfile, FastqInFileFWD):
    print("Running GSnap Indexing and SAI file build...")
    

    GSnap_index = "gmap_build -d ViralRefGMp %s" % (REFInfile)
    print(GSnap_index)
    status = os.system(GSnap_index)

    GSnap_saiF = "gsnap -d ViralRefGMp %s -A sam > /BWAfiles/BWAsam.sam" % (FastqInFileFWD)
    status = os.system(GSnap_saiF)

def BWA_SAM(REFInfile,FastqInFileFWD,FastqInFileREV,OutfileFASTa, OutfileSTAT, timeSTAMP, RefLEN, OutPNG, OutBLAST, OutCOV, OutCOVL, OutCOVpng):

    FastqInFileFWDassem = FastqInFileFWD #"/users/rwbarrettemac/bioinformatics/velvet/data/assem/FWDreads/contigs.fa"
    FastqInFileREVassem = FastqInFileREV #"/users/rwbarrettemac/bioinformatics/velvet/data/assem/REVreads/contigs.fa"
    OutfileFASTq = "/BWAfiles/FastqQUEOUT.fasta"
    #OutfileSTAT = OutfileFASTa+"_STATS_OUT.txt"

    RunBWA_SE(REFInfile, FastqInFileFWDassem, FastqInFileREVassem)

    RunSAM_SortedBam(REFInfile)
    #RunBAM_STATS(OutfileSTAT) #Added 2/15/2019 RWB

    RunSAM_Consensus(REFInfile, OutfileFASTa)
    RunBAM_STATS(OutfileSTAT)
  

    mLineTotal, mLineNonN = Convert_CONSENS_to_FASTA(OutfileFASTa, "/BWAfiles/FastqQUE_TEST.fasta", OutfileSTAT, timeSTAMP, RefLEN)


    MakeDB(REFInfile,"nucl","/BWAfiles/RawSequence_DIRECT")
    ###FIX_  BlastSequencesV29("/users/rwbarrettemac/bioinformatics/BWAfiles/FastqQUE_TEST.fasta", "blastn", OutfileSTAT)
    BlastSequencesV29_TAB("/BWAfiles/FastqQUE_TEST.fasta", "blastn")
    BlastSequencesV29_Xref("/BWAfiles/FastqQUE_TEST.fasta", "blastn", "/BWAfiles/FastqQUE_TEST.txt", OutfileSTAT)
    #RunBLAST_IMAGER_PL(OutPNG)

    df_COVERAGE = pd.read_csv("/BWAfiles/FastqQUE_DEPTH.txt", sep='\t', names=['ACC','Loc','Depth'])
    print(df_COVERAGE)
    #plt.figure();
    df_COVlim = df_COVERAGE.groupby(np.arange(len(df_COVERAGE))//100).mean()

    fig, ax = plt.subplots()

    print(df_COVlim)

    #Make copy of /users/rwbarrettemac/bioinformatics/BWAfiles/FastqQUE_DEPTH.txt but BC associated

    shutil.copy("/BWAfiles/FastqQUE_DEPTH.txt", OutCOV)
    shutil.copy("/BWAfiles/FastqQUE_DEPTHavg.txt", OutCOVL)

    df_COVlim.plot.line(x='Loc', y='Depth', style='b', ax=ax)
    plt.savefig(OutCOVpng)
    print("done")


def MakeDB(Inputfile, dbType, DBtitle):
    
    print("Building BLAST database...")
    MAKEDB_CMDLINE = "/NCBI/bin/makeblastdb -in %s -out %s -parse_seqids -dbtype %s"
    status = os.system(MAKEDB_CMDLINE % (Inputfile, DBtitle, dbType))


def File_Iterate(FilePath, RefGenome, timeSTAMP):
    handle = open(RefGenome)
    fa_REFSeq = SeqIO.parse(handle,"fasta")
    for NTseq in fa_REFSeq:
                fasta_sequence = NTseq.seq
                fasta_LEN = len(fasta_sequence)


    path = FilePath
    listing = os.listdir(path)

    for infile in listing:
        #try:

          if ".fastq.gz" in infile:
            continue #print("negatory...")
          else:
            if "R2_001" in infile:
              continue#print("negatory")
            else:
              if "._" in infile:
                continue#print("negatory")
              else:
                if infile!="none.fastq":
                  print("current file is: " + path + infile)
                  #print("POSITORY...")
                  RunFileR2 = infile.replace('_R1_001.', '_R2_001.')
                  RunFileR2x = path+RunFileR2
                  RunFileR1x = path+infile
                  FileRootID = infile.replace('L001_R1_001.fastq','')
                  OutStat = path+FileRootID+"OUT_STAT.txt"
                  OutFASTA = path+FileRootID+"OUT_CONS.fasta"
                  OutPNG = path+FileRootID+timeSTAMP+".png"
                  OutCOV = "/BWAfiles/FastqQUE_DEPTH/"+FileRootID+"COV.txt"
                  OutCOVL = "/BWAfiles/FastqQUE_DEPTH/"+FileRootID+"COVall.txt"
                  OutCOVpng = "/BWAfiles/FastqQUE_DEPTH/"+FileRootID+"COVall.png"

                  OutBLAST = "/OUTPUT/BLASTout.txt"
                  print(OutBLAST, "blasthere<<<<<")

                  print(RunFileR1x, RunFileR2x, OutStat, OutFASTA, timeSTAMP)
                  
                  BWA_SAM(RefGenome,RunFileR1x,RunFileR2x,OutFASTA, OutStat, timeSTAMP, fasta_LEN, OutPNG, OutBLAST, OutCOV, OutCOVL, OutCOVpng)


if __name__ == "__main__":
    FASTQinFOLDER, FASTAinREF, timeSTAMP = main(sys.argv[1:])   

    File_Iterate(FASTQinFOLDER, FASTAinREF, timeSTAMP)



