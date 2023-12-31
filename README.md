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
## Databases:
  As mentioned, the ClinVar and CCDS databases are contained in the Data folder, and can be downloaded from:
  + ClinVar: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
  + CCDS:
    1. https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS_protein.current.faa.gz
    2. https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS_nucleotide.current.fna.gz
    3. https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS2Sequence.current.txt   
  
  It is worth noting that the ClinVar database 'variant_summary.txt.gz' is updated on the first Thursday of each month. Therefore, the bash shell file "variant_summary_download.sh" is included in the Data folder and automatically downloads the new ClinVar database.

  On the other hand, the mentioned database must be downloaded and included in the Data folder.
## Content
The Classification app is a folder which contains different folders and files:

+ Data:

  This folder contains the ClinVar and CCDS databases, and the rest of the information is also downloaded here. However, this information is temporary and will be deleted after the database has been built. In addition, it contains the "Amino_Acid_Features.xlsx" and "CodigoGenetico.xlsx" Excel files which are required to calculate some of the descriptors.
  **This folder must be inside the Classification app folder**.
+ msas:

  This folder contains more folders that are called as the indicated gene, and contains the msas for the gene. However, there will be five different files for each folder:

    + .hhr
    + .a3m
    + .fil.a3m
    + .fil.fas
    + .aligned.gap.filtered.fas

  The first two files correspond to those returned by the hhblits programme. Then we applied a filter with hhfilter to the .a3m file to extract the aligned sequences without insertions, and then we changed the format from .a3m to .fas (FASTA) with reformat.pl, obtaining .fil.fas. Finally, we applied a second filter to avoid sequences with more than 50% gaps. In this way, we obtain the file .aligned.gap.filtered.fas, which we use to calculate the resiude conservation. **Keep in mind that hhblits, hhfilter and reformat.pl are scripts from HHSuite.**
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

  + \>DATA PATH=
  + \>UNIREF PATH=
  + \>BFD PATH=

  If these paths do not exit the application will report an error.
+ gene_app.py:

  This file corresponds to the execution script of the application. It will analyse whether a database already exists for the specified gene and, if so, it will automatically run the machine learning algorithm. Otherwise, it will create the database and then perform the classification.
+ gene_app.slurm:

  This file is used for the server queuing system.
## Usage:
  The application is initialized by argument parsing. Therefore, there are required and optional arguments.
  
  usage: gene_app.py [-h] [-i FILENAME FILENAME FILENAME FILENAME FILENAME FILENAME FILENAME FILENAME] [-g GEN_NAM] [-l LENGTH] -o OUTFILE [-s STATUS] [-w WEB_DECISION] -cpu CPU_NUMBER [-n_splits NUMBER_SPLITS] [-n_splits_grid NUMBER_SPLITS_GRID] [-column_names COLUMN_NAMES] -config_path CONFIG_PATH

  In this way, the arguments inside the [] are optional and the others are required. We can divide the arguments in two groups:

  + Database group: Parser arguments for database creation
    + Information:
      + **-i**: Indicating that the information related to the gene will be manually entered.
      + **-g** and **-l**: Otherwise, introducing the gene name (-g) and sequence length (-l). **If the length is not a number it will raise an error.**
      **You only can choose one option, -i or -g and -l. The application will raise an error in case of mixing them.**
    + **-s**: Status for secondary structure analysis. Regular or Extend. Regular will contain three different types of structures while Extend eight. (default = Regular) \<Regular\>
    + **-w**:  Status for the files downloading: -w <Internet> or -w <Internet/DataBase>. Default <Internet/DataBase>. **The Internet/DataBase option is to take the ClinVar information from the database instead of using the REST API. As mentioned before, I strongly recommend this option.**
    + **-cpu**:  Number of CPUs for hhblits. Recommendation: 1-4
    + **-o**: Name for the output file
 + Machine Learning group: Parser arguments for Machine Learning analysis
   + **-n_splits**:  Number of splits for cross validation. <default = 5>
   + **-n_splits_grid**: Number of splits for grid search. <default = 5>
   + **-column_names**: Column names of the database for the Machine Learning. Default: <My_Label> <residue_conserv> <d_size>
+ **-config_path**: Path for config.txt **It must be inside the Classification app folder**
+ **-h**: show this help message and exit


**Note that absolute paths must be specified correctly in the config.txt file.**
## Examples:
  These could be possible examples:
  1. `python gene_app.py -g $SLURM_JOB_NAME -l 872 -o $SLURM_JOB_NAME -cpu $SLURM_CPUS_PER_TASK -config_path "/home/einstein/martin/git_enviroment/Classification_app/config.txt" -column_names "My_Label,pLDDT,MTI_31" -n_splits_grid 14`
  2. `python gene_app.py -g $SLURM_JOB_NAME -l 872 -o $SLURM_JOB_NAME -cpu $SLURM_CPUS_PER_TASK -config_path "/home/einstein/martin/git_enviroment/Classification_app/config.txt" -column_names "My_Label,pLDDT,MTI_31,d_size,d_vol,d_pol_e,d_ip_e,d_hf_e,d_msa,residue_conserv" -n_splits_grid 14`

  where $SLURM_JOB_NAME corresponds to the gene name. In these two examples the column names are different. This argument parsing must be indicated in the gene_app.slurm file.

  On the other hand, the msas and RESULTS folders are created automatically when the analysis is finished. In addition, **the tmp folder must be created before starting the execution**.
  ## Data Download:
  + ClinVar:
    As mentioned above, this information is extracted from the ClinVar REST API or the variant_summary.txt.gz database. For using the REST API, I recommend reading https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/.
  + CCDS:
    The CCDS information is extracted from the databases. However, these databases are ordered by their own id, i.e. each gene has a CCDS id. Therefore, a query is made to Ensembl's REST API to extract both the CCDS and transcript ids for each gene. Once the CCDS id is known, it is used to search for the associated information in the CCDS databases.

    **The Ensembl REST API use: https://rest.ensembl.org/**
  + GenomAD:
    The GenomAD information is extracted using a Python package called Selenium, which simulates the user interacting with the website to click on the download button. This method is used because there is no REST API available and because the site is dynamic, so web scraping is not possible. Therefore, the transcript id extracted from the Ensembl REST API is used to download the information.
  + AlphaFold:
    The AlphaFold file is download after extracting the id that UniRef has used for classifying the particular gene. In this way, the UniRef id is obtained using the REST API: **https://rest.uniprot.org/uniprotkb/stream?format=json&query=%28%28gene%3ASCN7A%29%29+AND+%28model_organism%3A9606%29&sort=length+desc/**. This A url is for the SCN7A gene, but we can change the name to search any other gene.

    Therefore, once the UniRef id is obtained, it must be used to download the pdb file from AlphaFold as the url is constructed as **https://alphafold.ebi.ac.uk/files/AF-{UniRef_ID}-F1-model_v4.pdb**.
  + Lovd3:
    The Lovd3 file is obtained using webscraping, in other words, reading the html code. This method is used as no REST API is available or there is no a database.
## Errors:
  The application may raise an error if the gene information is entered incorrectly or if there is no information in the databases. This is checked for each file download, so it will raise an error if any of the websites or databases fail. In addition, some websites may not have information about the gene, so this will also be reported.

  Once the databases are built, there may be an error in training the machine learning models, as the clinical condition may only be of one class and therefore no classification training can be performed.


  
