
This pipeline analyses CHIP-seq raw data (SRA files) creating peak files that are later analysed using an R script and KEGG to get further information about the studied transcription factor, such as its target genes, functions of those genes and the metabolic pathways in which those genes are involved.

Apart from that, this pipeline gives the transcription factor sequence that binds DNA (the motif) using the findMotifsGenome tool from HOMER.

 

In order to do that, the pipeline has 5 different scripts and an additional parameters file with all input data that is needed for the mentioned analysis. Please, open this file and complete it before going forward with the analysis.

 

Before using this pipeline, please make sure you follow the next requirements:

-       You have installed R in your device and the following packages: ChIPseeker, TxDb.Athaliana.BioMart.plantsmart28, clusterProfiler, AnnotationDbi, org.At.tair.db, topGO, Rgraphviz y pathview.

-       You have installed HOMER and the genome of the organism of study. To install HOMER, it is imperative that you copy the PATH in the file ‘.profile’ in order for this pipeline to work.

 

To run this pipeline, execute the files ‘pipechipseq.sh’ and ‘params.txt’.

If an error is found, please check again you follow all requirements and all data in the parameters file are as expected.

 

Enjoy!

Bioinformática y Análisis Genómico. Universidad de Sevilla.
