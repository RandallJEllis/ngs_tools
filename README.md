# Tools for NGS analysis 🧬

### alignment
1. Generating indices for alignment with STAR from species reference genomes
2. Aligning zipped FASTQ files (.gz) to the reference genome index and quantifying counts on the gene or transcript level
3. Collating counts and running QC on alignment BAM files

### deg
DESeq on gene- or transcript-level counts

### pathways
Run pathway analysis on DEG results using the Enrichr API

### coexprs_network
1. Run MEGENA on expression data to identify gene co-expression modules
2. Calculate eigengenes of specific modules for each sample and correlate them with other sample-level measurements
3. Calculate module differential connectivity between experimental groups
4. Pathway/ontology enrichment analyses on modules
5. Network analysis measures (e.g., betweenness centrality)

### single_cell_deconv
1. Dampened weighted least squares for single-cell deconvolution using a reference single-cell RNA-seq dataset
2. Cell-type aware RNA-seq to run differential expression analysis on deconvolved cell fractions