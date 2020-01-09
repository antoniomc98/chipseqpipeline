## Authors: Antonio Melgar de la Cruz, Salvador Priego Poyato
##
## Contact: salvapriego@gmail.com, melgar.cruz.antonio@gmail.com
##
## ChIP-Seq Data Analysis
##
## Callpeak script


#! /bin/bash


# Reading input parameters

WD=$1
NCHIP=$2
SCRIPTS=$3
PROMOTER=$4
OUTPUT=$5


# Callpeak function

cd $WD/results

I=1

while [ $I -le $NCHIP ]
do
   macs2 callpeak -t $WD/samples/chip_sample_$I/chip_${I}_sorted.bam -c $WD/samples/input_sample_$I/input_${I}_sorted.bam -n peaks_${I} --outdir . -f BAM
   ((I++))
done

# Run HOMER analysis

findMotifsGenome.pl $WD/results/peaks_1_peaks.narrowPeak tair10 $WD/results/homer_results -size 10


# Run R script analysis

mkdir $WD/results/R_analysis
cd R_analysis
echo "Promoter length is $PROMOTER "
echo "Output file name is $OUTPUT "

Rscript $SCRIPTS/Rchip.R  $PROMOTER $OUTPUT $WD/results/peaks_1_peaks.narrowPeak

mv $HOME/chip* ../../logs/
mv $HOME/callpeak* ../../logs/
mv $HOME/input* ../../logs/

echo "All analyses are finished!" >> $WD/results/finish.txt
