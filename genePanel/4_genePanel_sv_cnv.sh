#! /bin/sh


SCRIPT=4.1_genePanel_SV_CNV_revision.py
SNVDIR=/home/cog/jvalleinclan/hpcMount/ovarian/WGS/data/somatic/genePanel_filtered #Dir with genePanel_filtered VCFs
SVDIR=/home/cog/jvalleinclan/hpcMount/ovarian/WGS/data/sv #Dir with SV VCFs
CNVDIR=/home/cog/jvalleinclan/hpcMount/ovarian/WGS/data/cnv #Dir with CNV files
OUTPUT=genePanel_sv_cnv.txt
touch $OUTPUT

 for SAMPLE in HGS-1 HGS-13.3 HGS-13.4 HGS-19 HGS-1-R1 HGS-1-R2 HGS-1-R3 HGS-2 HGS-3.1 HGS-3.2 HGS-4.1 HGS-4.3 HGS-4.5 HGS-5.1 HGS-6 LGS-1.1 LGS-1.2 LGS-1.3 LGS-1.4 LGS-2.2 LGS-3.1 LGS-3.2 LGS-5.2 LGS-5.4 MBT-1 MBT-2.1 MBT-2.2 MC-1.1 MC-1.2 MC-2.1 MC-2.2 SBT-3.1 SBT-3.2 SBT-4.1 SBT-4.2
 do
     for TYPE in tumor organoid
     do
         python $SCRIPT --snv $SNVDIR/${SAMPLE}.somatic.merged.filtered_genePanel_rescued_${TYPE}.vcf.gz --sv $SVDIR/${SAMPLE}_${TYPE}_somaticSV.vcf --cnv $CNVDIR/${SAMPLE}_${TYPE}_CNVs --sample $TYPE --patient $SAMPLE >> $OUTPUT
     done
 done

for SAMPLE in CCC-1
 do
     for TYPE in tumor organoid
     do
         python $SCRIPT --snv $SNVDIR/${SAMPLE}.pseudosomatic.merged.filtered_genePanel_rescued_${TYPE}.vcf.gz --sv $SVDIR/${SAMPLE}_${TYPE}_somaticSV.vcf --cnv $CNVDIR/${SAMPLE}_${TYPE}_CNVs --sample $TYPE --patient $SAMPLE >> $OUTPUT
     done
 done

for SAMPLE in END-1
do
     TYPE="organoid-A"
     python $SCRIPT --snv $SNVDIR/END-1.pseudosomatic.merged.filtered_genePanel_rescued_organoid.A.vcf.gz --sv $SVDIR/END-1_organoid_A_somaticSV.vcf --cnv $CNVDIR/END-1-A_organoid_CNVs --sample $TYPE --patient $SAMPLE >> $OUTPUT
     TYPE="organoid-B"
     python $SCRIPT --snv $SNVDIR/END-1.pseudosomatic.merged.filtered_genePanel_rescued_organoid.B.vcf.gz --sv $SVDIR/END-1_organoid_B_somaticSV.vcf --cnv $CNVDIR/END-1-B_organoid_CNVs --sample $TYPE --patient $SAMPLE >> $OUTPUT
     TYPE="tumor"
     python $SCRIPT --snv $SNVDIR/END-1.pseudosomatic.merged.filtered_genePanel_rescued_tumor.vcf.gz --sv $SVDIR/END-1_tumor_A_somaticSV.vcf --cnv $CNVDIR/END-1_tumor_CNVs --sample $TYPE --patient $SAMPLE >> $OUTPUT
     
 done
 
 
for SAMPLE in HGS-1 HGS-1-R2 HGS-2 HGS-3.1 HGS-3.2 HGS-6
 do
     TYPE="organoidLate"
     python $SCRIPT --snv $SNVDIR/${SAMPLE}.somatic.merged.filtered_genePanel_rescued_${TYPE}.vcf.gz --sv $SVDIR/${SAMPLE}_${TYPE}_somaticSV.vcf --cnv $CNVDIR/${SAMPLE}_${TYPE}_CNVs --sample $TYPE --patient $SAMPLE >> $OUTPUT
 done

#SOME MODS TO MAKE LIFE EASIER LATER
sed -i "s:tumor:T:g" $OUTPUT
sed -i "s:organoidLate:O-S:g" $OUTPUT
sed -i "s:organoid:O:g" $OUTPUT
sed -i "s:disruptive_inframe_insertion:inframe_insertion:g" $OUTPUT
sed -i "s:inframe_insertion:inframe_insertion:g" $OUTPUT
sed -i "s:disruptive_inframe_deletion:inframe_deletion:g" $OUTPUT
sed -i "s:inframe_deletion:inframe_deletion:g" $OUTPUT


Rscript 4.2_genePanel_sv_cnv.R $OUTPUT ${OUTPUT/.txt/.pdf}
