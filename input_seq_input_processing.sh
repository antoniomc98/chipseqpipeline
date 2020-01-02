## Authors: Antonio Melgar de la Cruz, Salvador Priego Poyato
##
## Contact: salvapriego@gmail.com, melgar.cruz.antonio@gmail.com
##
## ChIP-Seq Data Analysis
##
## Input processing script


#! /bin/bash

# Input parameters

WD=$1
NS=$2
NINPUT=$3
INPUT_ID=$4
SCRIPTS=$5
PROMOTER=$6
OUTPUT=$7


# Access to sample folder

cd $WD/samples/input_sample_${INPUT_ID}


# Performing quality control and mapping to reference genome

fastqc input_${INPUT_ID}.fastq
bowtie2 -x $WD/genome/index -U input_${INPUT_ID}.fastq -S input_${INPUT_ID}.sam


# Generating bam file

samtools view -@ 2 -S -b input_${INPUT_ID}.sam > input_${INPUT_ID}.bam
rm input_${INPUT_ID}.sam
samtools sort input_${INPUT_ID}.bam -o input_${INPUT_ID}_sorted.bam
samtools index input_${INPUT_ID}_sorted.bam


# Integration

echo "input_${INPUT_ID} DONE" >> $WD/logs/blackboard.txt

DONE_INPUT=$(wc -l $WD/logs/blackboard.txt | awk '{ print $1 }' )

if [ $DONE_INPUT -eq $NS ]
then
   qsub -N callpeak -o $WD/logs/callpeak $SCRIPTS/callpeak.sh $WD $NINPUT $SCRIPTS $PROMOTER $OUTPUT
fi

