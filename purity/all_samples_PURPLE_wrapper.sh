#! /bin/bash

BAM_SRC=/hpc/cog_bioinf/kloosterman/users/jvalleinclan/ovarian/WGS/data/bam/
RES_SRC=/hpc/cog_bioinf/kloosterman/users/jvalleinclan/ovarian/WGS/purple/results/

COBALT=/hpc/cog_bioinf/kloosterman/users/jvalleinclan/ovarian/WGS/purple/run_cobalt.sh
AMBER=/hpc/cog_bioinf/kloosterman/users/jvalleinclan/ovarian/WGS/purple/run_amber.sh
PURPLE=/hpc/cog_bioinf/kloosterman/users/jvalleinclan/ovarian/WGS/purple/run_purple.sh


for PATIENT_DIR in $BAM_SRC/*/
do
	REF_BAM=$PATIENT_DIR/B*.bam
	REF=`echo $REF_BAM | rev | cut -f1 -d'/' | rev | cut -f1 -d'.'`
	for TUM_BAM in $PATIENT_DIR/T*.bam $PATIENT_DIR/O*.bam
	do
		
		TUM=`echo $TUM_BAM | rev | cut -f1 -d'/' | rev | cut -f1 -d'.'`
		cd $RES_SRC
		mkdir $TUM
		cd $TUM

		mkdir amber
		cd amber
		qsub -N AMBER_${TUM} $AMBER $REF_BAM $TUM_BAM $TUM
		cd ..

		mkdir cobalt
		cd cobalt
		qsub -hold_jid AMBER_${TUM} -N COBALT_${TUM} $COBALT $REF $REF_BAM $TUM $TUM_BAM
		cd ..

		qsub -hold_jid COBALT_${TUM} -N PURPLE_${TUM} $PURPLE $REF $TUM
		
	done
done
