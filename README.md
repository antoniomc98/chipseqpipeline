
This pipeline analyses CHIP-seq raw data (SRA files) creating peak files that are later analysed using an R script and KEGG to get further information about the studied transcription factor, such as its target genes, functions of those genes and the metabolic pathways in which those genes are involved.

Apart from that, this pipeline gives the transcription factor sequence that binds DNA (the motif) using the findMotifsGenome tool from HOMER.


In order to do that, the pipeline has 5 different scripts and an additional parameters file with all input data that is needed for the mentioned analysis. Please, open this file and complete it before going forward with the analysis.


Before using this pipeline, please make sure you follow the next requirements:

- You have installed R and the following packages: ChIPseeker, TxDb.Athaliana.BioMart.plantsmart28, clusterProfiler, AnnotationDbi, org.At.tair.db, topGO, Rgraphviz y pathview.

- You have installed HOMER and the genome of the organism of study. To install HOMER, it is imperative that you copy the PATH in the file ‘.profile’ in order for this pipeline to work.


To run this pipeline, execute the files ‘pipechipseq.sh’ with ‘params.txt’.

If an error is found, please check again you follow all requirements and all data in the parameters file are as expected.


Enjoy!

Bioinformática y Análisis Genómico. Universidad de Sevilla.

References

R: R Core Team (2019). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
- CHIPseeker package: Guangchuang Yu, Li-Gen Wang, Qing-Yu He. ChIPseeker: an R/Bioconductor package for ChIP peak     annotation, comparison and visualization. Bioinformatics 2015, 31(14):2382-2383
- Bioconductor: Orchestrating high-throughput genomic analysis with Bioconductor. W.Huber, V.J. Carey, R. Gentleman, ..., M. Morgan Nature Methods, 2015:12,115.
- clusterProfiler package: Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
- AnnotationDbi: Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2019). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.46.1.
- org.At.tair.db: Marc Carlson (2019). org.At.tair.db: Genome wide annotation for Arabidopsis. R package version 3.8.2
- topGO: Adrian Alexa and Jorg Rahnenfuhrer (2019). topGO: Enrichment Analysis for Gene Ontology. R package version 2.36.0.
- Rgraphviz: Kasper Daniel Hansen, Jeff Gentry, Li Long, Robert Gentleman, Seth Falcon, Florian Hahne and Deepayan Sarkar (2019). Rgraphviz: Provides plotting capabilities for R graph objects. R package version 2.28.0.
- pathview: Luo, W. and Brouwer C., Pathview: an R/Bioconductor package for pathway-based data integration and visualization. Bioinformatics, 2013,29(14): 1830-1831, doi: 10.1093/bioinformatics/btt285

KEGG pathways: Kanehisa, M.; Toward pathway engineering: a new database of genetic and molecular pathways. 
Science & Technology Japan, No. 59, pp. 34-38 (1996).
Herramienta HOMER: Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription 
Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432
