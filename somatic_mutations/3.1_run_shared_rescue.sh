
BAMDIR=/home/cog/jvalleinclan/uproov_wd/bams #Directory with BAM files
VCFDIR=.
OUTDIR=.
SCRIPT=3_patient_variants_withrescue.py

for AF in 0.05
do
	for PATIENT in ${ARRAY WITH PATIENT_NAMES}
	do
		python $SCRIPT --noref -f $VCFDIR/$PATIENT.somatic.merged.filtered.vcf -b $BAMDIR -a $AF > $OUTDIR/$PATIENT.shares.af.$AF.txt
	done
done


