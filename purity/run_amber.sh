#! /bin/sh
#$ -m beas
#$ -M jespejov@umcutrecht.nl
#$ -l h_rt=48:0:0
#$ -l h_vmem=80G
#$ -cwd
#$ -pe threaded 8

REF=$1

TUM=$2

SAMPLE=$3

HET_BED=/hpc/cog_bioinf/kloosterman/tools/hmftools_pipeline_v4/dbSNP.20180423.vcf.2millPositions.bed

sambamba mpileup -t 8 -L $HET_BED $REF --samtools -q 1 -f /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta > REFERENCE.mpileup

sambamba mpileup -t 8 -L $HET_BED $TUM --samtools -q 1 -f /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta > TUMOR.mpileup

module load Java/1.8.0_60

java -jar ~/kloosterman/tools/hmftools_pipeline_v4/amber_v1-5.jar -sample $SAMPLE -output_dir . -reference REFERENCE.mpileup -tumor TUMOR.mpileup
