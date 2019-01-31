#! /bin/bash

#Example without filenames for anonimity
module load Java/1.8.0_60
java -Xmx16G -jar /hpc/local/CentOS7/cog_bioinf/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T CombineVariants -R /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta --variant:${ORGANOID} ${BLOOD}_${ORGANOID}_merged_somatics_snpEff_dbSNP_Cosmicv76_melted.vcf --variant:${TUMOR} ${BLOOD}_${TUMOR}_merged_somatics_snpEff_dbSNP_Cosmicv76_melted.vcf -o ${SAMPLE}.somatic.merged.vcf
