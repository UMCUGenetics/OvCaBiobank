#! /bin/sh
###RUN ON LOCAL MACHINE!!
#Directory with all the bam files (or symlinks of course)
BAMDIR=. 

for VCF in *.merged.filtered_genePanel.vcf
do
	python rescue_genePanel.py --vcf $VCF --bamDir $BAMDIR > ${VCF/.vcf/_rescued.vcf}
done

