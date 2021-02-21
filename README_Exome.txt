############## SCRIPTS DESCRIPION #############################################

This folder contains the following three scripts.
script1: preparationExome.sh
script2: pipelineExome.sh
script3: quality_controlExome.sh --optional

These scripts perform genetic variants analysis using exome sequencing data.

########## RUNNING THE SCRIPTS ################################################

To work on the Biosphere VMs, the scripts should be placed in the directory 
"/home/ubuntu/data/mydatalocal/" and run with the command ./scriptname.sh.

The scripts should be launched strictly in the following order:
0) conda activate
1) preparationExome.sh 
2) pipelineExome.sh 

For an optional (but highly recommended!) quality control step, the script 
quality_controlExome.sh can be run at any timepoint following the step (1).

Path variables introduced at the beginning of each script can be changed
to run the pipeline in the environment other than the Biosphere VM.

#Comment: following the installations in preparationExome.sh, running the 
"varscan somatic" command produced VCF file and not MAF 
(contrary to what have been written in the email from 30 November 2020)

###############################################################################
Anna KARPUKHINA
11/12/2020






