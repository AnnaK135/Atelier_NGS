############## SCRIPTS DESCRIPION #############################################

This folder contains the following three scripts.
script1: preparationRNAseq.sh
script2: pipelineRNAseq.sh
script3: quality_controlRNAseq.sh --optional

These scripts perform analysis of the RNA-seq data and produce the table with raw counts
(to be further analyzed for differential expression, eg with DESeq2).

########## RUNNING THE SCRIPTS ################################################

To work on the Biosphere VMs, the scripts should be placed in the directory 
"/home/ubuntu/data/mydatalocal/" and run with the command ./scriptname.sh.

The scripts should be launched strictly in the following order:
0) conda activate
1) preparationRNAseq.sh 
2) pipelineRNAseq.sh 

For an optional (but highly recommended!) quality control step, the script 
quality_controlRNAseq.sh can be run at any timepoint following the step (1).

Path variables introduced at the beginning of each script can be changed
to run the pipeline in the environment other than the Biosphere VM.

###############################################################################
Anna KARPUKHINA
11/12/2020






