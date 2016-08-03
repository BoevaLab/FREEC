# Control-FREEC
Copy number and genotype annotation from whole genome and whole exome sequencing data.
  
---------------------------------------------------------------------------------------------------------------------------

Control-FREEC is a tool for detection of copy-number changes and allelic imbalances (including LOH) using deep-sequencing data originally developed by the Bioinformatics Laboratory of Institut Curie (Paris). Since 2016, the project has moved to Insitut Cochin, INSERM U1016 (Paris).

Control-FREEC automatically computes, normalizes, segments copy number and beta allele frequency (BAF) profiles, then calls copy number alterations and LOH. The control (matched normal) sample is optional for whole genome sequencing data but mandatory for whole exome or targeted sequencing data. For whole genome sequencing data analysis, the program can also use mappability data (files created by GEM).

Starting from version 8.0, we provide a possibility to detect subclonal gains and losses and evaluate the likeliest average ploidy of the sample. Also, the evaluation procedure for the level of contamination by normal cells has been improved.

---------------------------------------------------------------------------------------------------------------------------

**Input for CNA detection:** aligned single-end, paired-end or mate-pair data in SAM, BAM, SAMtools pileup.

Control-FREEC accepts .GZ files. Support of Eland, BED, SOAP, arachne, psl (BLAT) and Bowtie formats has been discontinued starting from version 8.0.

**Input for CNA+LOH detection:** There are two options: (a) provide aligned reads in SAMtools pileup format. Files can be GZipped; (b) provide BAM files together with options "makePileup" and "fastaFile" (see "How to create a config file?" in the manual.)

**Output:** Regions of gains, losses and LOH, copy number and BAF profiles.

---------------------------------------------------------------------------------------------------------------------------

**Installation:** To install FREEC, type "make" in the command line. If you are using Linux 32bit, please remove 64bit-tags from the Makefile file before building the program.

---------------------------------------------------------------------------------------------------------------------------

To cite please use:

Boeva V, Zinovyev A, Bleakley K, Vert JP, Janoueix-Lerosey I, Delattre O, Barillot E. (2011) Control-free calling of copy number alterations in deep-sequencing data using GC-content normalization. Bioinformatics 2011; 27(2):268-9. PMID: 21081509.

Boeva V, Popova T, Bleakley K, Chiche P, Cappo J, Schleiermacher G, Janoueix-Lerosey I, Delattre O, Barillot E. (2011) Control-FREEC: a tool for assessing copy number and allelic content using next generation sequencing data. Bioinformatics. 2011 Dec 6. [Epub ahead of print] PubMed PMID: 22155870.

---------------------------------------------------------------------------------------------------------------------------

People who contributed to Control-FREEC:

   - Valentina Boeva
   - Andrei Zinovyev
   - Tatiana Popova
   - Carino Gurjao
   - Kevin Bleakley
   - Pierre Chiche
   - Jean-Philippe Vert
   - Isabelle Janoueix
   - Joern Toedling
   - Emmanuel Barillot
   - Olivier Delattre 

---------------------------------------------------------------------------------------------------------------------------

**Contact:** We will be pleased to address any question or concern you may have with the Control-FREEC software: valentina.boeva %at% inserm.fr

---------------------------------------------------------------------------------------------------------------------------

This work was supported by grants from the Institut National de la Sante et de la Recherche Medicale, the Institut Curie, the Ligue Nationale contre le Cancer (Equipe labellisee and CIT program). 
