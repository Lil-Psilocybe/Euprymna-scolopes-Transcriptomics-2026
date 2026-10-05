#!/bin/bash
#SBATCH --job-name=Single_Nuclei_RNASeq_Briseño_et_al_2026.sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=300G
#SBATCH --mail-user=
#SBATCH -o /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/Single_Nuclei_RNASeq_Briseño_et_al_2026.sh_%j.out
#SBATCH -e /home/FCAM/jbriseno/JB_Nyholm_Lab/Scripts/Error.Output/Single_Nuclei_RNASeq_Briseño_et_al_2026.sh_%j.err



#what is this script:
#produce spliced/unspliced count matrices of single-nucleus RNASeq (snRNASeq) reads
#from adult male E. scolopes white body (MWB, sample was A7) tissue to be analyzed for RNA velocity

#experiment details:
#Conducted in Spring 2024
#ParseBioSciences V2 WT MiniKit was captured ~5000 MWB nuclei
#Sequenced on an Illumina Nova Seq at (insert sequencing depth)

#General workflow here:
#Pre step - extracting male white body reads from shared Parse sublibrary
#1 - clean reads of noisy barcodes with kallisto-bustools (kb)
#2 - align cleaned reads to the genome using Parse's Split-Pipe pipeline to produce a BAM file
#3 - recreate geomic annotation GTF with BAM file to capture 3' UTRs of reads that are canonically missed in snRNASeq with GeneExt
#4 - pseudo-align reads for quantification with updated GTF using kb





###Loading modules###
#module load fastqc

###activating conda envs###
#source ~/.bashrc
#conda activate spipe
#conda activate geneext
#conda activate kb_env



########################################################
###Generating and Setting Variables/Directories/Paths###
########################################################

#Adding V3 run on 9/18/2026

#latest assemby dir/files
Es_V3_dir=/labs/Nyholm/Es_V3_Genome/GCF_053919585.1/
Es_V3_genome=${Es_V3_dir}GCF_053919585.1_ASM5391958v1_genomic.fna
Es_V3_GTF=${Es_V3_dir}genomic.gtf

#working directory
V3_ReAnalyses_Dir=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/ParseV2_GenomeV3/

#Parse directories
V3_parse_genome=/labs/Nyholm/scsnRNASeq/Parse/

#sn reads
SL2_read1=SL2_S3_L002_R1_001.fastq.gz
SL2_read2=SL2_S3_L002_R2_001.fastq.gz
MA7_cleaned_reads=/labs/Nyholm/scsnRNASeq/Parse/MA7.WB_FB4.ANG_ParseMiniV2_2024/kb_cleaned_reads/ 

#GeneExt and kb
geneext_dir=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/ParseV2_GenomeV3/genomes/geneext/
kb_dir=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/kb_dir/
cleaned_dir=/labs/Nyholm/scsnRNASeq/Parse/MA7.WB_FB4.ANG_ParseMiniV2_2024/kb_cleaned_reads/

sorted_bam=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/ParseV2_GenomeV3/analysis/process/sorted_barcode_headAligned_anno.bam
sorted_cleaned_bam=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/ParseV2_GenomeV3/analysis/process/sorted_cleaned_barcode_headAligned_anno.bam
agat_gtf=/labs/Nyholm/Es_V3_Genome/GCF_053919585.1/agat_genomic.gtf
geneext_dir=/home/FCAM/jbriseno/JB_Nyholm_Lab/Sept2026_Ch2_BulkDataV3_ReAnalyses/geneext/

kb_dir=/labs/Nyholm/scsnRNASeq/sn_kb_cleaned_reads/
output_dir=/home/FCAM/jbriseno/JB_Nyholm_Lab/cleaned_fastqs_tmp/raised_atlas_cleaned_fastqs/
bustools_path=/home/FCAM/jbriseno/miniconda3/envs/kb_env/lib/python3.14/site-packages/kb_python/bins/linux/bustools/bustools



####################
###END USER INPUT###
####################



###Pre step###
#GOAL: partitioning raw reads by sample to pull out white body from ParseV2 mini run 
#10K nuclie sequenced - ~5K from Male A7 white body and ~5k from Female B4 ANG
cd /labs/Nyholm/scsnRNASeq/Parse/MA7_WB.FB4_ANG/expdata/SL2/
python /home/FCAM/jbriseno/JB_programs/fastq_sep_groups_v0.5.py \
--chemistry v2 \
--fq1 SL2_S3_L002_R1_001.fastq.gz \
--fq2 SL2_S3_L002_R2_001.fastq.gz \
--opath . \
--group MA7_WB A1-A6 \
--group FB4_ANG A7-A12



###1 - Kallisto-Bustools###
#GOAL: clean reads of unwanted barcodes
#Script from CO at CalTech

#loading conda env
source ~/.bashrc
conda activate kb_env

cd /labs/Nyholm/scsnRNASeq/sn_kb_cleaned_reads/

#creates reference of all possible k-mers (this part only has to be run once)
#echo 'create reference of all possible k-mers (this part only has to be run once)'
#echo {A,T,C,G}{A,T,C,G}{A,T,C,G}| tr ' ' '\n'|awk '{print ">0\n" $0}' > 3mers.txt
#echo "0 0" > 3mers.t2g
#kb ref --workflow=custom --distinguish --overwrite -k 3 -i 3mers.idx 3mers.txt

#actual filtering
f1=/labs/Nyholm/scsnRNASeq/Parse/MA7_WB.FB4_ANG/expdata/SL2/SL2_S3_L002_group_MA7_WB_R1.fastq.gz
f2=/labs/Nyholm/scsnRNASeq/Parse/MA7_WB.FB4_ANG/expdata/SL2/SL2_S3_L002_group_MA7_WB_R2.fastq.gz

echo 'Finding Reads'
#barcodes found at: https://github.com/Yenaled/barcodes
kb count \
-x SPLIT-SEQ \
-w splitseqv2_barcodes.txt \
-r splitseqv2_replace.txt \
-i 3mers.idx \
-g 3mers.t2g \
-o . \
--num $f1 $f2 \
--overwrite
#where do barcodes come from? my udi plate was v2 for sure
#checkout -x technology for parse here -> https://kallisto.readthedocs.io/en/latest/sc/webpages/splitseq_preprocessing.html
#this works
#manually changed names on all the .json files (inspect.json, run_info.json, kb_info.json) for each SL

echo 'Converting to Text'
$bustools_path sort --flags-bc -o output.sorted.bus output.bus
$bustools_path sort --flags -o output_modified.unfiltered.sorted.bus output_modified.unfiltered.bus
$bustools_path text -o output_modified.unfiltered.txt output_modified.unfiltered.sorted.bus
#this works

echo 'Extracting Good Barcodes'
awk '{print $1}' output_modified.unfiltered.txt | uniq | sort | uniq > good_barcodes.txt
#this works

echo 'Filtering Bus by Barcodes'
$bustools_path capture -o good_barcodes.bus -c good_barcodes.txt --barcode output.sorted.bus
#this works

echo 'Extract Good fastqs'
$bustools_path extract \
good_barcodes.bus \
-N 2 \
-o $output_dir \
--include -f $f1,$f2
#this works



###2 - ParseBioSciences SplitPipe###
#GOAL: produce a BAM file of aligned MWB reads

cd $V3_ReAnalyses_Dir

###genome indexing###
source ~/.bashrc
conda activate spipe
module load fastqc

split-pipe \
--mode mkref \
--genome_name Es_V3_Csomes \
--fasta $Es_V3_genome \
--genes $Es_V3_GTF \
--output_dir /labs/Nyholm/scsnRNASeq/Parse/ParseGenomeIndex_V3

###quanitfying cleaned reads with Parse Index###
split-pipe \
--mode all \
--chemistry v2 \
--kit WT_mini \
--genome_dir $V3_parse_genome/ParseGenomeIndex_V3/ \
--fq1 /labs/Nyholm/scsnRNASeq/Parse/MA7.WB_FB4.ANG_ParseMiniV2_2024/kb_cleaned_reads/Cleaned_MA7_WB_R1.fastq.gz \
--fq2 /labs/Nyholm/scsnRNASeq/Parse/MA7.WB_FB4.ANG_ParseMiniV2_2024/kb_cleaned_reads/Cleaned_MA7_WB_R2.fastq.gz \
--output_dir analysis/ \
--sample MA7_WB A1-A6



###3 - GeneExt###
#GOAL: produce a UTR-aware GTF with the snRNASeq alignment itself to inform the annotation


##BAM prep##
#need to resort parse-generated 
module load samtools

cd $V3_ReAnalyses_Dir/analysis/process
samtools sort barcode_headAligned_anno.bam -o sorted_barcode_headAligned_anno.bam
samtools view -H sorted_barcode_headAligned_anno.bam | grep @HD


cd /home/FCAM/jbriseno/JB_programs/GeneExt-main/
#copied GeneExtMain from V3 genome dir


#ah yes, miust strip Es-V3 from bam file csomes!
cd $V3_ReAnalyses_Dir/analysis/process
python /home/FCAM/jbriseno/JB_programs/split_pipe/bam_noprefix_add_sublib.py \
--in_bam $sorted_bam \
--out_bam sorted_cleaned_barcode_headAligned_anno.bam \
--strip_prefix "Es-V3_" \
--sublib_number "2"

#fix "genes with missign exons"
source ~/.bashrc
conda activate agat

agat_convert_sp_gxf2gxf.pl -g $Es_V3_GTF -o $Es_V3_dir/agat_genomic.gtf


#running geneext
source ~/.bashrc
conda activate geneext

python geneext.py \
-g $geneext_dir/agat_genomic.gtf \
-b $sorted_cleaned_bam \
-o $geneext_dir/MWB_geneext.gtf \
-v 1 \
--force \
--peak_perc 0
#this worked!!!
#now rerun split pipe



###4 - KB###
#GOAL: pseudo-align cleaned reads with GeneExt GTF for (un)splice matrices for single-cell velocity


#9.23.26 - kbtools

####building a velocity index###
#https://www.kallistobus.tools/tutorials/kb_velocity_index/python/kb_velocity_index/
cd $kb_dir
kb ref \
--tmp TMP --verbose \
-t 8 --aa --workflow lamanno --overwrite \
-i index.idx \
-g t2g.txt \
-f1 cdna.fa -f2 intron.fa \
-c1 cdna_t2c.txt -c2 intron_t2c.txt \
-t 10 \
$Es_V3_genome $geneext_dir/MWB_geneext.gtf

###pseudo-alignment and quantification###
kb count \
-x SPLIT-SEQ \
-w onlist.txt -r replace.txt \
-i index.idx -g t2g.txt \
-c1 cdna_t2c.txt -c2 intron_t2c.txt \
--h5ad --workflow nucleus \
--num -t 10 -m 300G --report \
-o MWB_kb --overwrite \
$cleaned_dir/Cleaned_MA7_WB_R1.fastq.gz $cleaned_dir/Cleaned_MA7_WB_R2.fastq.gz
#completed 9.23!
#now move to scvelo


#9.19.26 - created scvelo conda env and will follow tutorial for this
#9.22.26 - scvelo works well
