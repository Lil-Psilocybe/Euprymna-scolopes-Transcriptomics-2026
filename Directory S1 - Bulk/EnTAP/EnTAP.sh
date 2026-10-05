#!/bin/bash
#SBATCH --job-name=EnTAP.sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=200G
#SBATCH --exclude=mantis-023
#SBATCH -o /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/EnTAP.sh_%j.out
#SBATCH -e /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/EnTAP.sh_%j.err
echo `hostname`


####################
###Module Loaders###
####################
module load entap/2.3.0-jmykac3

EsV3EnTAP=/labs/Nyholm/Es_V3_Genome/GCF_053919585.1/kb_entap_res/

#################################
#EnTAP for functional annotation#
#################################
#(written by Wegryzyn Lab)
#Combines RSEM, TransDecoder (HMMER), DIAMOND, EggNOG, &, InterProScan woah!

cd ${EsV3EnTAP}
EnTAP --run \
--run-ini entap_run.params \
--entap-ini real_entap_config.ini \
-t 10
#This works! - 1/10-26