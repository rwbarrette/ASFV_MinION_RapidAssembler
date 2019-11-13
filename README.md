# ASFV_MinION_RapidAssembler

https://www.ncbi.nlm.nih.gov/pubmed/?term=31694969

O'Donnell VK, Grau FR, Mayr GA, Sturgill Samayoa TL1, Dodd KA, Barrette RW.,
RAPID SEQUENCE-BASED CHARACTERIZATION OF AFRICAN SWINE FEVER VIRUS USING THE OXFORD NANOPORE MINION SEQUENCE SENSING DEVICE AND A COMPANION ANALYSIS SOFTWARE TOOL.
J Clin Microbiol. 2019 Nov 6. pii: JCM.01104-19. doi: 10.1128/JCM.01104-19.

Foreign Animal Disease Diagnostic Laboratory, National Veterinary Services Laboratories, Animal and Plant Health Inspection Service, United States Department of Agriculture, Plum Island Animal Disease Center, New York.

---------------------------------------------------------------------------------------------------------------------------

ASF_FAST Companion Software for rapid sequence assembly in real-time.  Intended for use with the Oxford Nanopore MinION Sequence Sensing device.  Software can be used with or without SMS reporting, and will generate output data within the raw fastq folders.  

Software is provided as both source code and as a Docker Image.  It is recommended that the software be built using the docker image for platform compatability and to ensure that all dependencies are available.   Source code is provided for reference and has not been adapted for easy deployment and would require extensive modification for use.

Please note!  Demultiplexing by Porechop as originally described and implemented for generation of data in the manuscript has been replaced with the QCAT demultiplexer for improved performance (https://github.com/nanoporetech/qcat).


**INSTALLATION INSTRUCTIONS**

1.  Download and install Docker desktop. 

  https://www.docker.com/products/docker-desktop

2.  Download ASF_FAST docker image.

3.  Load Docker Image.  Enter into command line with Docker Desktop running.

  docker load --input ASFfast.tar
  

**OPTIONAL INSTALLATION OF TWILIO FOR SMS RESULTS REPORTING**
 (note: currently generates report for only single barcode of interest for SMS reporting)

5.  Install ngrok

  https://github.com/inconshreveable/ngrok

**Setting up ngrok/Twilio for optional SMS reporting**

Optional for SMS based reporting:

  Setup ngrok service
  
   1.  Start ngrok
    
      ngrok http 5000
      
   2.  Copy from ngrok session window Forwarding http:
   
 ![](/images/ngrok_Screenshot1.png)
 
   3.  Setup Twilio account and environment
   
    https://www.twilio.com/docs/sms/quickstart/python
   
   Log into Twilio and paste Forwarding http: for SMS.  

**RUNNING ASF-FAST**

  1.  Run Docker Image ASFfast and linking docker folder with input/output folders
  
    docker run -ti -v /{Input_Output_Folder_Location}:/OUTPUT python-barcode /bin/bash
    
  2.  Set up Oxford Nanopore MinION to export basecalled FASTQ raw data in  files to a new folder in {Input_Output_Folder_Location}
  
  3. start ASF-FAST software within Docker pointing at the folder where the FASTQ raw data files are being copied to
  
    
  
  
