#! /bin/bash
module load Java/1.8.0_60

#Filter if it has caller support CSP >= 2, no MT or Y chromosome variants and read count at least 10 in the tumor or 5 in the organoid

for VCF in *.somatic.merged.vcf
do
java -jar /hpc/local/CentOS7/cog_bioinf/snpEff_v4.2/SnpSift.jar filter "(!(exists CSP) | (CSP >= 2)) & (CHROM != 'MT') & (CHROM != 'Y') & ((GEN[0].AD[1] >= 10) | (GEN[1].AD[1] >= 5))" $VCF > ${VCF/.vcf/.filtered.vcf}
done

