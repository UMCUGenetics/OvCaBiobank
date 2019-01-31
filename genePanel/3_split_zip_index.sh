#! /bin/sh

module load Java/1.8.0_60

for VCF in *.somatic.merged.filtered_genePanel_rescued.vcf
do
	vcf-sort -c $VCF > ${VCF/.vcf/_sorted.vcf}	

	for SAMPLE in `grep ^#CHROM $VCF | cut -f10-`
do
	if [[ "$SAMPLE" == *"-S" ]]; then
		TYPE="organoidLate"
	elif [[ "$SAMPLE" == "O"* ]]; then
		TYPE="organoid"
	else
		TYPE="tumor"
	fi
	echo $SAMPLE $TYPE
	java -Xmx16G -jar /hpc/local/CentOS7/cog_bioinf/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T SelectVariants \
	-R /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta -V ${VCF/.vcf/_sorted.vcf}  \
	-se $SAMPLE | awk '$10 != "./."' > ${VCF/.vcf/_$TYPE.vcf}
	bgzip ${VCF/.vcf/_$TYPE.vcf}
	tabix -p vcf ${VCF/.vcf/_$TYPE.vcf.gz}
done
done


#CCC-1 and END-1 go on their own due to different regexp needed

for VCF in CCC-1.pseudosomatic.merged.filtered_genePanel_rescued.vcf
do
	vcf-sort -c $VCF > ${VCF/.vcf/_sorted.vcf}	

	for SAMPLE in `grep ^#CHROM $VCF | cut -f10-`
do
	if [[ "$SAMPLE" == *"-S" ]]; then
		TYPE="organoidLate"
	elif [[ "$SAMPLE" == "O"* ]]; then
		TYPE="organoid"
	else
		TYPE="tumor"
	fi
	echo $SAMPLE $TYPE
	java -Xmx16G -jar /hpc/local/CentOS7/cog_bioinf/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T SelectVariants \
	-R /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta -V ${VCF/.vcf/_sorted.vcf}  \
	-se $SAMPLE | awk '$10 != "./."' > ${VCF/.vcf/_$TYPE.vcf}
	bgzip ${VCF/.vcf/_$TYPE.vcf}
	tabix -p vcf ${VCF/.vcf/_$TYPE.vcf.gz}
done
done


for VCF in END-1.pseudosomatic.merged.filtered_genePanel_rescued.vcf
do
	vcf-sort -c $VCF > ${VCF/.vcf/_sorted.vcf}	

	for SAMPLE in O76-A O76-B T76
do
	if [[ "$SAMPLE" == *"-A" ]]; then
		TYPE="organoid.A"
	elif [[ "$SAMPLE" == "O"* ]]; then
		TYPE="organoid.B"
	else
		TYPE="tumor"
	fi
	echo $SAMPLE $TYPE
	java -Xmx16G -jar /hpc/local/CentOS7/cog_bioinf/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T SelectVariants \
	-R /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta -V ${VCF/.vcf/_sorted.vcf}  \
	-se $SAMPLE | awk '$10 != "./."' > ${VCF/.vcf/_$TYPE.vcf}
	bgzip ${VCF/.vcf/_$TYPE.vcf}
	tabix -p vcf ${VCF/.vcf/_$TYPE.vcf.gz}
done
done
