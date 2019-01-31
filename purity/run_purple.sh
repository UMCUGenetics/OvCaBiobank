
#! /bin/sh
#$ -m beas
#$ -M jespejov@umcutrecht.nl
#$ -l h_rt=4:0:0
#$ -l h_vmem=20G
#$ -cwd
#$ -pe threaded 8

REF=$1
TUM=$2

module load Java/1.8.0_60

java -jar ~/kloosterman/tools/hmftools_pipeline_v4/purple_v2-14.jar \
   -run_dir . \
   -ref_sample $REF \
   -tumor_sample $TUM \
   -threads 8 \
   -gc_profile ~/kloosterman/tools/hmftools_pipeline_v4/GC_profile.1000bp.cnp \
 
