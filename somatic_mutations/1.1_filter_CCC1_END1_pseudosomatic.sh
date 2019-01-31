#! /bin/sh

#THESE TWO SAMPLES DO NOT HAVE A REFERENCE
#IN ORDER TO ENRICH FOR SOMATIC CALLS, REMOVE ALL CALLS IN 1000GENOMES AND GENOME OF THE NETHERLANDS
#ALSO, SELECT FOR MODERATE OR HIGH IMPACT FROM SNPEFF

module load Java/1.8.0_60

PATIENT=END-1 #HERE GOES END-1 or CCC-1
VCF=${PATIENT}.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf
OUT=${VCF/.vcf/_1kg.vcf}
KG=/hpc/cog_bioinf/kloosterman/resources/1kg/phase3-b37/raw/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
#It had been annotated by the IAP with GoNL already

java -Xmx16g -jar /hpc/local/CentOS7/cog_bioinf/snpEff_v4.1h/SnpSift.jar  filter "(ANN[*].IMPACT = 'MODERATE') | (ANN[*].IMPACT = 'HIGH')" $VCF > ${VCF/.vcf/_moderatehigh.vcf}
java -Xmx16g -jar /hpc/local/CentOS7/cog_bioinf/snpEff_v4.1h/SnpSift.jar annotate -tabix -name 1KG_ -info AF,AN,AC $KG ${VCF/.vcf/_moderatehigh.vcf} > ${VCF/.vcf/_1kg_moderatehigh.vcf}

grep '^#' ${VCF/.vcf/_1kg_moderatehigh.vcf} > ${VCF/.vcf/_1kg_moderatehigh_filtered.vcf}
grep -v '^#' ${VCF/.vcf/_1kg_moderatehigh.vcf} | grep -v '1KG' | grep -v 'GoNL' >> ${PATIENT}.pseudosomatic.merged.filtered.vcf

