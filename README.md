# Classification-app
The purpose of this application is to create a database and run a machine learning classification algorithm for any gene. Therefore, the application would need information about the gene in order to build such a database. So the information can be entered manually or the application can download it automatically. However, the main idea of this programme was to be as automatic as possible, and then it is recommended to enter the name and protein sequence length of the particular gene, and the application will download all the information.

The information that is needed will be extracted from databases and websites:
+ ClinVar
+ AlphaFold
+ Ensembl
+ Lovd3
+ GenomAD
+ UniProt
+ CCDS

For example, the ClinVar information can be taken from the website REST API or the database, although I strongly recommend using the database as it is much faster. On the other hand, the CCDS information is also taken from the database and the rest from websites using REST APIs or other extraction methods.

In this way, the main idea of this application is to specify the gene name and sequence length, and the program will return a database that will be used by the machine learning algorithm to perform the classification training.
## Requisites
This application uses some python packages that must be installed.

* pandas
* biopandas
* biopython
* numpy
* mdtraj
* requests

I have been using the anaconda environment that can be obtained from https://www.anaconda.com/download. Therefore, all the packages mentioned above can be installed with both pip and conda.

On the other hand, the application uses the hhblits programme from HHSuite to construct the MSA for the residue conservation calculation. This program requires a lot of resources, since you have to enter a FASTA sequence of the gene as input, and it will extract all the sequences related to this one, comparing them with the huge databases it needs, since they contain a large number of sequences, and in this way it is able to relate the input sequence to those in the database. Note that this search may take more or less time depending on the device you have. Thus, I recommend to read the HHSuite user guide (https://github.com/soedinglab/hh-suite/wiki) to avoid any problem you can have.
## Content
The Classification app is a folder which contains different folders and files:

+ Data:

  This folder contains the ClinVar and CCDS databases, and the rest of the information is also downloaded here. However, this information is temporary and will be deleted after the database has been built. **This folder must be inside the Classification app folder**.
+ msas:

  This folder contains more folders that are called as the indicated gene, and contains the msas for the gene. However, there will be five different files:

    + .hhr
    + .a3m
    + .fil.a3m
    + .fil.fas
    + .aligned.gap.filtered.fas

  The first two files correspond to those returned by the hhblits programme. Then we applied a filter with hhfilter to the .a3m file to extract the aligned sequences without insertions, and then we changed the format from .a3m to .fas (FASTA) with reformat.pl, obtaining .fil.fas. Finally, we applied a second filter to avoid sequences with more than 50% gaps. In this way, we obtain the file .aligned.gap.filtered.fas, which we use to calculate the resiude conservation. **hhblits, hhfilter and reformat.pl are scripts from HHSuite.**
+ Packages:

  This folder contains some of the packages that are used in the application.
+ RESULTS:

  This folder is returned by the application. It will contain different folders with the name of the genes, and these ones contain the results for each gene:

   + .xlsx
   + _encoded.xlsx
   + Descriptor_folder

  The first two files correspond to the databases. However, the first one contains all the information related to the gene and the second one is prepared for the machine learning classification algorithm, as it has been performed a hot encoding on it and the information not required has been deleted. On the other hand, the Descriptor_folder (it is named as the descriptors used for the classification) contains the files with the performance statistics, and the GRAPHS folder contains the confusion matrices.

+ Scripts:

  This folder contains the scripts for the database creation and the machine learning classification training:

   + Data_Set_Creation:

     This folder contains all the scripts needed to create the database, from the scripts that download the information to the ones that merge it. It is worth noting that there are three different versions, but the Data_Set_Creation_3 is the latest.

   + Machine_Learning:
 
     This folder contains the required scripts to perform the machine learning classification algorithm.
+ tmp:

  This folder contains the .err and .out files for each gene. These files contain the possible errors during execution and the output information respectively.
+ config.txt:

  **This file must be inside the Classification app folder**. It contains the absolute paths for the Data folder and the databases needed for hhblits:

  + >DATA PATH=
  + >UNIREF PATH=
  + >BFD PATH=
+ gene_app.py:

  This file corresponds to the execution script of the application. It will analyse whether a database already exists for the specified gene and, if so, it will automatically run the machine learning algorithm. Otherwise, it will create the database and then perform the classification.
+ gene_app.slurm
  
