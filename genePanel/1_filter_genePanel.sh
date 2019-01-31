#! /bin/sh



VCFDIR=.
#USES SAME INPUT AS SHARED MUTATIONS

for VCF in $VCFDIR/*.merged.filtered.vcf
do
OUT=`echo $VCF | rev | cut -f1 -d'/' | rev | sed 's:.vcf:_genePanel.vcf:g'`
echo $OUT

#Create an HPC job
cat << EOF > ${OUT}_job.sh
#! /bin/bash
#$ -cwd
#$ -m ea
#$ -M jespejov@umcutrecht.nl
#$ -l h_rt=4:0:0
#$ -N ${OUT/_genePanel.vcf/}

module load Java/1.8.0_60
java -jar /hpc/local/CentOS7/cog_bioinf/snpEff_v4.1h/SnpSift.jar filter " (!(ANN[0].EFFECT = 'sequence_feature')) &&  ((ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')) && \
((ANN[*].GENE = 'USH2A' ) |  \
(ANN[*].GENE = 'RECQL4' ) |  \
(ANN[*].GENE = 'CSMD3' ) |  \
(ANN[*].GENE = 'ATR' ) |  \
(ANN[*].GENE = 'BCORL1' ) |  \
(ANN[*].GENE = 'GABRA6' ) |  \
(ANN[*].GENE = 'MYH2' ) |  \
(ANN[*].GENE = 'FLG2' ) |  \
(ANN[*].GENE = 'MYH1' ) |  \
(ANN[*].GENE = 'PPP1R3A' ) |  \
(ANN[*].GENE = 'MECOM' ) |  \
(ANN[*].GENE = 'PIK3CA' ) |  \
(ANN[*].GENE = 'MAP3K1' ) |  \
(ANN[*].GENE = 'LRP1B' ) |  \
(ANN[*].GENE = 'MAP3K4' ) |  \
(ANN[*].GENE = 'PARK2' ) |  \
(ANN[*].GENE = 'MYC' ) |  \
(ANN[*].GENE = 'FGFR1' ) |  \
(ANN[*].GENE = 'BRCA1' ) |  \
(ANN[*].GENE = 'BRCA2' ) |  \
(ANN[*].GENE = 'EXT1' ) |  \
(ANN[*].GENE = 'SYNE1' ) |  \
(ANN[*].GENE = 'FGFR2' ) |  \
(ANN[*].GENE = 'KRAS' ) |  \
(ANN[*].GENE = 'TRIOBP' ) |  \
(ANN[*].GENE = 'ARID1A' ) |  \
(ANN[*].GENE = 'GTSE1' ) |  \
(ANN[*].GENE = 'FAT3' ) |  \
(ANN[*].GENE = 'PTEN' ) |  \
(ANN[*].GENE = 'FAT4' ) |  \
(ANN[*].GENE = 'SLC19A1' ) |  \
(ANN[*].GENE = 'CIC' ) |  \
(ANN[*].GENE = 'PPP2R1A' ) |  \
(ANN[*].GENE = 'MMP16' ) |  \
(ANN[*].GENE = 'PTK2' ) |  \
(ANN[*].GENE = 'CDKN2A' ) |  \
(ANN[*].GENE = 'PLEC' ) |  \
(ANN[*].GENE = 'DAB2' ) |  \
(ANN[*].GENE = 'NDRG1' ) |  \
(ANN[*].GENE = 'RPTOR' ) |  \
(ANN[*].GENE = 'APOB' ) |  \
(ANN[*].GENE = 'CDK12' ) |  \
(ANN[*].GENE = 'CCNE1' ) |  \
(ANN[*].GENE = 'NRAS' ) |  \
(ANN[*].GENE = 'RB1' ) |  \
(ANN[*].GENE = 'EPPK1' ) |  \
(ANN[*].GENE = 'CXCR1' ) |  \
(ANN[*].GENE = 'WWOX' ) |  \
(ANN[*].GENE = 'AURKA' ) |  \
(ANN[*].GENE = 'NF1' ) |  \
(ANN[*].GENE = 'PRKCI' ) |  \
(ANN[*].GENE = 'MDC1' ) |  \
(ANN[*].GENE = 'PARP1' ) |  \
(ANN[*].GENE = 'TP53BP1' ) |  \
(ANN[*].GENE = 'CDKN2B' ) |  \
(ANN[*].GENE = 'MMP26' ) |  \
(ANN[*].GENE = 'RASA1' ) |  \
(ANN[*].GENE = 'TP53' ) |  \
(ANN[*].GENE = 'E2F1' ) |  \
(ANN[*].GENE = 'BRAF' ) |  \
(ANN[*].GENE = 'AHNAK2' ) |  \
(ANN[*].GENE = 'CREBBP' ) |  \
(ANN[*].GENE = 'AGO2' )\
)" $VCF > $OUT
EOF

qsub ${OUT}_job.sh
done
