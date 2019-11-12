# ASFV_MinION_RapidAssembler

https://www.ncbi.nlm.nih.gov/pubmed/?term=31694969

O'Donnell VK, Grau FR, Mayr GA, Sturgill Samayoa TL1, Dodd KA, Barrette RW.,
RAPID SEQUENCE-BASED CHARACTERIZATION OF AFRICAN SWINE FEVER VIRUS USING THE OXFORD NANOPORE MINION SEQUENCE SENSING DEVICE AND A COMPANION ANALYSIS SOFTWARE TOOL.
J Clin Microbiol. 2019 Nov 6. pii: JCM.01104-19. doi: 10.1128/JCM.01104-19.

Foreign Animal Disease Diagnostic Laboratory, National Veterinary Services Laboratories, Animal and Plant Health Inspection Service, United States Department of Agriculture, Plum Island Animal Disease Center, New York.

---------------------------------------------------------------------------------------------------------------------------

ASF_FAST Companion Software for rapid sequence assembly in real-time.  Intended for use with the Oxford Nanopore MinION Sequence Sensing device.

Software is provided as both source code and as a Docker Image.  It is recommended that the software be built using the docker image for platform compatability and to ensure that all dependencies are available.   Source code is provided for reference and has not been adapted for easy deployment and would require extensive modification for use.

Please note!  Demultiplexing by Porechop as originally described and implemented for generation of data in the manuscript has been replaced with the QCAT demultiplexer for improved performance (https://github.com/nanoporetech/qcat).
