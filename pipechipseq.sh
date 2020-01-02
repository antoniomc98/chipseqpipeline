## Authors: Antonio Melgar de la Cruz, Salvador Priego Poyato
##
## Contact: salvapriego@gmail.com, melgar.cruz.antonio@gmail.com
##
## ChIP-Seq Data Analysis
##
## Pipeline script


#! /bin/bash

if [ $# -eq 0 ]
  then
   echo "This pipeline analyzes ChIP-seq Data"
   echo "Usage: pipechipseq <param_files>"
   echo ""
   echo "param_file: file with the parameters specifications. Please, check params.txt for an example"
   echo "Lay back and let me do the work"

   exit 0
fi

# Introducing parameters

PARAMS=$1

WD=$(grep working_directory: $PARAMS | awk '{ print $2 }' )
INT=$(grep int: $PARAMS | awk '{ print $2 }' )
GENOME=$(grep genome: $PARAMS | awk '{ print $2 }' )
ANNOTATION=$(grep annotation: $PARAMS | awk '{ print $2 }' )
NS=$(grep number_of_samples: $PARAMS | awk '{ print $2 }' )
NCHIP=$(grep number_chip_samples: $PARAMS | awk '{ print $2 }' )
NINPUT=$(grep number_input_samples: $PARAMS | awk '{ print $2 }' )
SCRIPTS=$(grep scripts: $PARAMS | awk '{ print $2 }' )
PROMOTER=$(grep promoter: $PARAMS | awk '{ print $2 }' )
OUTPUT=$(grep output_file: $PARAMS | awk '{ print $2 }' )


echo "Working directory:" $WD
echo "Number of samples:" $NS
echo "Genome downloaded/copied from:" $GENOME
echo "Annotation downloaded/copied from:" $ANNOTATION
echo "Number of chip samples:" $NCHIP
echo "Number of input samples:" $NINPUT




SAMPLES=()


I=0
J=0
K=0

while [ $I -lt $NS ]
do
   if [ $J -lt $NCHIP ]
   then
      SAMPLES[$I]=$(grep chip_sample_$(($J + 1)): $PARAMS | awk '{ print $2 }')
      ((I++))
      ((J++))
   elif [ $K -lt $NINPUT ]
   then
      SAMPLES[$I]=$(grep input_sample_$(($K + 1)): $PARAMS | awk '{ print $2 }')
      ((I++))
      ((K++))
   fi
done


I=0
J=0
K=0

while [ $I -lt $NS ]
do
   if [ $J -lt $NCHIP ]
   then
      echo chip_sample_$(($J+1)) = ${SAMPLES[$I]}
      ((J++))
      ((I++))
   elif [ $K -lt $NINPUT ]
   then
      echo input_sample_$(($K+1)) = ${SAMPLES[$I]}
      ((K++))
      ((I++))
   fi
done


# Creating working directory

mkdir $WD
cd $WD
mkdir genome annotation samples results logs
cd samples


I=1
J=1
K=1

while [ $I -le $NS ]
do
   if [ $J -le $NCHIP ]
   then
      mkdir chip_sample_$J
      ((I++))
      ((J++))
   elif [ $K -le $NINPUT ]
   then
      mkdir input_sample_$K
      ((I++))
      ((K++))
   fi
done

# Getting reference genome

cd $WD/genome

if [ $INT == "YES" ]
then
   wget -O genome.fa.gz $GENOME
   gunzip genome.fa.gz
else
   cp $GENOME genome.fa
fi

# Downloading annotation

cd $WD/annotation

if [ $INT == "YES" ]
then
   wget -O annotation.gtf.gz $ANNOTATION
   gunzip annotation.gtf.gz
   else
   cp $ANNOTATION annotation.gtf
fi


# Creating reference genome index

cd $WD/genome
bowtie2-build genome.fa index

# Downloading/Copying samples

cd ../samples/
if [ $INT == "NO" ]
then
   I=0
   J=0
   K=0
   while [ $I -lt $NS ]
   do
      if [ $J -lt $NCHIP ]
      then
         cp ${SAMPLES[$I]} chip_sample_$(($J + 1))/chip_$(($J + 1)).fastq
         ((I++))
         ((J++))
      elif [ $K -lt $NINPUT ]
      then
         cp ${SAMPLES[$I]} input_sample_$(($K + 1))/input_$(($K + 1)).fastq
         ((I++))
         ((K++))
      fi
   done
else
   I=0
   J=0
   K=0
   while [ $I -lt $NS ]
   do
      if [ $J -lt $NCHIP ]
      then
         fastq-dump --split-files ${SAMPLES[$I]} -O ./chip_sample_$(($J+1))
         mv ./chip_sample_$(($J+1))/${SAMPLES[$I]}* ./chip_sample_$(($J+1))/chip_$(($J+1)).fastq
         ((I++))
         ((J++))
      elif [ $K -lt $NINPUT ]
      then
         fastq-dump --split-files ${SAMPLES[$I]} -O ./input_sample_$(($K+1))
         mv ./input_sample_$(($K+1))/${SAMPLES[$I]}* ./input_sample_$(($K+1))/input_$(($K+1)).fastq
         ((I++))
         ((K++))
      fi
   done
fi


# Executing sample processing scripts

I=1
J=1
K=1
while [ $I -le $NS ]
do
   if [ $J -le $NCHIP ]
   then
      qsub -N chip_$J -o $WD/logs/chip_$J $SCRIPTS/chip_seq_chip_processing.sh $WD $NS $NCHIP $J $SCRIPTS $PROMOTER $OUTPUT
      ((I++))
      ((J++))
   elif [ $K -le $NINPUT ]
   then
      qsub -N input_$K -o $WD/logs/input_$K $SCRIPTS/input_seq_input_processing.sh $WD $NS $NINPUT $K $SCRIPTS $PROMOTER $OUTPUT
      ((I++))
      ((K++))
   fi
done
