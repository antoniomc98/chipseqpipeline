## Authors: Antonio Melgar de la Cruz, Salvador Priego Poyato
##
## Contact: salvapriego@gmail.com, melgar.cruz.antonio@gmail.com
##
## ChIP-Seq Data Analysis
##
## Chip processing script


#! /bin/bash

# Reading input parameters

WD=$1
NS=$2
NCHIP=$3
CHIP_ID=$4
SCRIPTS=$5
PROMOTER=$6
OUTPUT=$7


# Adress sample folder

cd $WD/samples/chip_sample_${CHIP_ID}


# Sample quality control and read mapping to reference genome

fastqc chip_${CHIP_ID}.fastq
bowtie2 -x $WD/genome/index -U chip_${CHIP_ID}.fastq -S chip_${CHIP_ID}.sam


# Generting sorted bam file

samtools view -@ 2 -S -b chip_${CHIP_ID}.sam > chip_${CHIP_ID}.bam
rm chip_${CHIP_ID}.sam
samtools sort chip_${CHIP_ID}.bam -o chip_${CHIP_ID}_sorted.bam
samtools index chip_${CHIP_ID}_sorted.bam

# Synchronization

echo "chip_${CHIP_ID} DONE" >> $WD/logs/blackboard.txt

DONE_CHIP=$(wc -l $WD/logs/blackboard.txt | awk '{ print $1 }' )

if [ $DONE_CHIP -eq $NS ]
then
  qsub -N callpeak -o $WD/logs/callpeak $SCRIPTS/callpeak.sh $WD $NCHIP $SCRIPTS $PROMOTER $OUTPUT
fi

