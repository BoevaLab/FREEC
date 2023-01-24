# Control-FREEC
Copy number and genotype annotation from whole genome and whole exome sequencing data.

Complete tutorial on how to run FREEC: http://boevalab.inf.ethz.ch/FREEC/tutorial.html
  
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

**Installation:** 

To install FREEC, type "make" in the command line. If you are using Linux 32bit, please remove 64bit-tags from the Makefile file before building the program.

Alternatively, you may experiment with a FREEC executable on any OS, w/o installing make and C++ compiler, by following the steps in the **Using Control-FREEC from a container** section.

***Using Control-FREEC from a container***

You may experiment with FREEC using a [Docker image](https://www.docker.com/resources/what-container/), prepared by the project's development team. By default, the image uses a FREEC executable, compiled from the latest code on the **master** branch. Please see [the list of available versions](https://hub.docker.com/repository/docker/knotnote/control-freec/tags?page=1&ordering=last_updated) in case you want to use another version of FREEC.

1. [Install Docker](https://docs.docker.com/get-docker/).
2. Get the `knotnote/control-freec:latest` docker image:
   - if you are using Docker Desktop, follow [these instructions](https://docs.docker.com/desktop/dashboard/#pull-the-latest-image-from-docker-hub)
   - or from the console: `docker pull knotnote/control-freec:latest`
3. Get the config file that you plan to use. Modify all file paths to be absolute.
   
   Example:

         [general]
         chrLenFile = hs18_chr.len

   becomes:
         
         [general]
         chrLenFile = /home/user/data/hs18_chr.len

4. Run FREEC. Use the following command: 
   
   `docker run --rm -t --mount type=bind,source=abs_path_to_data,target=abs_path_to_data knotnote/control-freec:latest abs_path_to_config_file`

   where:

      - abs_path_to_data - **absolute** path to the directory with input data **and** a configuration file. 

      - abs_path_to_config_file - **absolute** path to the config file.

   The following options are being used:

      - --rm - remove the container after it finishes execution to spare disk space
      
      - --mount - bind-mount a local directory into the container. The data is **not exposed** to the Internet and stays on the private machine.

   **NOTE.** If your files are located in different folders, you need to add the \'--mount type=bind,source=<path>,target=<path>\' for each folder.


   If set up correctly, FREEC outputs will be in the `outputDir` folder, specified in the configuration file.

(optional) 5. Running an interactive experiment.

   If you want to get into the container and modify its contents (for e.g., update the FREEC code and recompile the executable), you may use the following command:

   `docker run -it --name=experiment1 --entrypoint="" --mount type=bind,source=abs_path_to_data,target=abs_path_to_data knotnote/control-freec:latest abs_path_to_config_file bash`

   Here:

   - --name=experiment1 - gives your container a unique name (in case you want to build your own image based on the experiments)
 
   - '--entrypoint=""' and 'bash' - allows to get into the console inside of the container.

   This could also **be useful** for debugging purposes (validate that data is properly mounted, the code version, etc.)

***Creating a new FREEC container version (for developers)***

From the folder with the Dockerfile, please run the command:

   `docker build -t knotnote/control-freec:<FREEC version> .`

To make the image publicly available:

1. Ensure that you have authenticated to the Dockerhub registry:

   `docker login`

2. Push the image:

   `docker push knotnote/control-freec:<FREEC version>`

Note. To push to the *knotnote/control-freec* Dockerhub repository, please ask collaborators to grant you permissions.

---------------------------------------------------------------------------------------------------------------------------

**Running the compiled version:** 

'PathToFREEC/freec -conf myConfig.txt -sample sample.bam -control control.bam'

or

'PathToFREEC/freec -conf myConfig.txt' if BAM files are provided in the config file.

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
