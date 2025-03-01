![image](https://github.com/user-attachments/assets/4b5c0a80-718e-4b75-b224-8e231e63db06)# Code for analyses described in Doxey et al.

* Data:
  - raw FASTQ files available from the NCBI SRA: accession # PRJNA1212169

* Code:
  -  `/docs/*.sh`: scripts for running Salmon, de-hosting, etc.
  -  `/docs/Rmarkdown.html`: initial version of analysis
  -  `/docs/Rmarkdown-V2.html`: updated analysis including negative controls
  -  `/docs/bacterialTranscriptomicsAnalysis.R`: code for analysis of MCAT/SPN/HFLU gene expression
  -  `/docs/hostResponseClassifier.R`: code for constructing host-response random-forest classifier
* Data required for analysis:
  - `d1-d5_fastp_dehost_krakenmerged_reads_normalized.txt`:  Taxonomic predictions for all samples
  - `gencode.v39.metadata.HGNC`: reference human transcriptome (https://www.gencodegenes.org/)
  - `/quants/`: gene expression data for all samples - generate from fastq files using Salmon
