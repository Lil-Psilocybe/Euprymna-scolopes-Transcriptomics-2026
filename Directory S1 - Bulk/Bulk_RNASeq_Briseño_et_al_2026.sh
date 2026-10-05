#!/bin/bash
#SBATCH --job-name=Bulk_RNASeq_Briseño_et_al_2026.sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mem=300G
#SBATCH --mail-user=
#SBATCH -o /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/Bulk_RNASeq_Briseño_et_al_2026_%j.out
#SBATCH -e /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/Bulk_RNASeq_Briseño_et_al_2026_%j.err



######################
##Variables & Paths##
#####################

###########
#TRIM & QC#
###########
#Raw reads trim and QC of demultiplexed data recieved in 2 batches from the Simakov lab
#where first batch (sequenced DATE) is 
FIRST_BATCH_RAW_READS=/labs/Nyholm/Oleg_project/HHH72DSX2_1_R11994_20210815/demultiplexed/
#and second batch (sequenced DATE) is 
SECOND_BATCH_RAW_READS=/labs/Nyholm/Oleg_project/HCVFYDSX3_1_R13227_20220410/demultiplexed/
#Raw data #Illumina read variables
R1=L001_R1_001.fastq.gz
R2=L001_R2_001.fastq.gz
#Other QC paths
PRE_TRIM=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/QC/PRE_TRIM/
POST_TRIM=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/QC/POST_TRIM/
MULTIQC=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/QC/MULTIQC/
LOGS=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/QC/LOGS/
#Output cleaned reads
CLEANED_READS=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/READS/CLEANED_READS/
Unpaired_Reads=/labs/Nyholm/JB_Nyholm_Lab/Ceph_Omes/Transcriptomes/E_scolopes/READS/Unpaired_Reads/
FastP=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/FastP/



###############
###ALIGNMENT###
###############
#Adding V3 run on 9/18/2026
Es_V3_dir=/labs/Nyholm/Es_V3_Genome/GCF_053919585.1/
Es_V3_genome=${Es_V3_dir}GCF_053919585.1_ASM5391958v1_genomic.fna
Es_V3_GTF=${Es_V3_dir}genomic.gtf
STAR_EsV3_dir=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/STAR_GenomeV3_BRAKER_BULK/
TAG=Parent
FEATURE=exon
CPUs=20
BAM=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/STAR_GenomeV3_BRAKER_BULK/BAM/
FEATURECOUNTS=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/FEATURECOUNTS/


####################
###Quantification###
####################
#Quantifying reads to genes (exons) with featureCounts 
FEATURECOUNTS=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/FEATURECOUNTS/


#Organ Sample Prefixes
MNTLs="
167717_S2
190853_S1
190910_S9
190922_S21
190931_S28
"
SKINs="
167720_S5
190904_S3
190912_S11
190924_S23
190933_S30
"
HECT0_ARMs="
167719_S4
190854_S2
190911_S10
190923_S22
190932_S29
"
CBs="
190909_S8
190914_S13
190921_S20
167744_S29
190939_S36
190954_S51
"
OLLs="
167722_S7
167728_S13
190906_S5
167736_S21
167745_S30
190935_S32
"
WBLs="
167730_S15
190913_S12
190916_S15
190918_S17
190907_S6
190925_S24
167737_S22
167747_S32
190951_S48
190936_S33
"
OVRYs="
167749_S34
167741_S26
190934_S31
190965_S53
"
TSTs="
167721_S6
167734_S19
190905_S4
190915_S14
"
LO_CCs="
190908_S7
190917_S16
190920_S19
167743_S28
167751_S36
190938_S35
190953_S50
"
#Combining mutiple variables into one
#https://stackoverflow.com/questions/8446146/combining-multiple-variables-into-another-variable-in-unix
ORGANS="${HECT0_ARMs} ${MNTLs} ${SKINs} ${WBLs} ${CBs} ${OLLs} ${OVRYs} ${TSTs} ${LO_CCs}"


####################
###END USER INPUT###
####################




###############
###TRIM & QC###
###############
#For EsV3 - 9.18.2026
#trimming with fastp instead of trimmomatic
#using multiqc on fastp json files collectively 

#loadin' modules
module load MultiQC

###conda###
source ~/.bashrc
conda activate fastp
#done!

cd ${FastP}
for organ in ${ORGANS} #This can be subbed for specific organs -> to use just use the organ! like "for wb in ${WB}"
	do
		#fastp command
		fastp \
		--in1 ${FIRST_BATCH_RAW_READS}"$organ"_${R1} \
		--in2 ${FIRST_BATCH_RAW_READS}"$organ"_${R2} \
		--out1 ${FastP}"$organ"_R1.fq.gz \
		--out2 ${FastP}"$organ"_R2.fq.gz \
		-j ${FastP}"$organ".json
		
		fastp \
		--in1 ${SECOND_BATCH_RAW_READS}"$organ"_${R1} \
		--in2 ${SECOND_BATCH_RAW_READS}"$organ"_${R2} \
		--out1 ${FastP}"$organ"_R1.fq.gz \
		--out2 ${FastP}"$organ"_R2.fq.gz \
		-j ${FastP}"$organ".fastp.json
	done

#multiqc check on fastp json files
multiqc . --force



###############
###ALIGNMENT###
###############
#Indexing genome for Bulk STAR alignment

#loadin' modules
module load star/2.7.11b #aligmment

STAR \
--runThreadN ${CPUs} \
--runMode genomeGenerate \
--genomeDir ${STAR_EsV3_dir} \
--genomeFastaFiles ${Es_V3_genome} \
--sjdbGTFfile ${Es_V3_GTF} \
--sjdbGTFfeatureExon ${FEATURE} \
--sjdbOverhang 149 #150 bp read length so N - 1 = 149 for Overhang
#9.18.26 - completed genome generation

#Aligning cleaned reads from TRIM & QC section to this newly indexed genome
cd ${BAM}
for organ in ${ORGANS} #This can be subbed for specific organs -> to use just use the organ! like "for wb in ${WB}"
	do
		#Mappin to V2 of the genome using STAR 2pass alignment 
		echo "Basic two-pass mapping mode to ID splice junctions"
		STAR \
		--runThreadN ${CPUs} \
		--runMode alignReads \
		--runDirPerm All_RWX \
		--sjdbOverhang 149 \
		--twopassMode Basic \
		--sjdbGTFfeatureExon ${FEATURE} \
		\
		--genomeDir ${STAR_EsV3_dir} \
		--sjdbGTFfile ${Es_V3_GTF} \
		\
		--readFilesIn ${FastP}"$organ"_R1.fq.gz ${FastP}"$organ"_R2.fq.gz \
		--readFilesCommand gunzip -c \
		\
		--outFileNamePrefix "$organ" \
		--outTmpKeep None \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterMultimapNmax 1 
	done



####################
###Quantification###
####################
##Quantification with FeatureCounts##

#loadin' modules
module load subread #quantification via freaturecounts (v2.1.1)

cd $FEATURECOUNTS
featureCounts \
-T 10 \
-a $Es_V3_GTF \
-t exon -g gene_id \
-p --countReadPairs \
-o featureCounts.matrix.txt \
--donotsort \
$BAM/*.bam