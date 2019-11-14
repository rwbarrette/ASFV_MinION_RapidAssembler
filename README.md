# ASFV_MinION_RapidAssembler

https://www.ncbi.nlm.nih.gov/pubmed/?term=31694969

O'Donnell VK, Grau FR, Mayr GA, Sturgill Samayoa, Dodd KA, Barrette RW.,
RAPID SEQUENCE-BASED CHARACTERIZATION OF AFRICAN SWINE FEVER VIRUS USING THE OXFORD NANOPORE MINION SEQUENCE SENSING DEVICE AND A COMPANION ANALYSIS SOFTWARE TOOL.
J Clin Microbiol. 2019 Nov 6. pii: JCM.01104-19. doi: 10.1128/JCM.01104-19.

Foreign Animal Disease Diagnostic Laboratory, National Veterinary Services Laboratories, Animal and Plant Health Inspection Service, United States Department of Agriculture, Plum Island Animal Disease Center, New York.

---------------------------------------------------------------------------------------------------------------------------

**!!! PLEASE NOTE: SOURCE CODE FOR ASF_FAST IS CURRENTLY AVAILABLE HERE.  THE SELF-CONTAINED DOCKER IMAGE TO SIMPLIFY TOOL DEPLOYMENT ACROSS PLATFORMS WILL BE MADE AVAILABLE BY DECEMBER 1**

ASF_FAST Companion Software for rapid sequence assembly in real-time.  This software is intended for use with the Oxford Nanopore MinION Sequence Sensing device.  Software can be used with or without SMS reporting, and will generate output data within the raw fastq folders. 

User must indicate location of basecalled FASTQ files generated by the MinION device at runtime, and designate the same folder within the MinKNOW software.

Software is provided currently as source code but will be released as a Docker Image soon to improvel usabilty and provide cross-platform flexibilty.  It is recommended that the software be built using the docker image for platform compatability and to ensure that all dependencies are available.   Source code is provided primarily for reference and has not been adapted for easy deployment and would require extensive modification for use.

Please note!  Demultiplexing by Porechop as originally described and implemented for generation of data in the manuscript has been replaced with the QCAT demultiplexer for improved performance (https://github.com/nanoporetech/qcat).



  
  
