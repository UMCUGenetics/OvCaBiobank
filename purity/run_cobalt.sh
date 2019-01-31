
#! /bin/sh
#$ -m beas
#$ -M jespejov@umcutrecht.nl
#$ -l h_rt=6:0:0
#$ -l h_vmem=20G
#$ -cwd
#$ -pe threaded 8

REF=$1
REF_BAM=$2

TUM=$3
TUM_BAM=$4

OUT_DIR=.

module load Java/1.8.0_60

java -jar ~/kloosterman/tools/hmftools_pipeline_v4/cobalt_v1-4.jar -reference $REF -reference_bam $REF_BAM -tumor $TUM -tumor_bam $TUM_BAM -output_dir . -threads 8 -gc_profile ~/kloosterman/tools/hmftools_pipeline_v4/GC_profile.1000bp.cnp
