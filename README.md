# Code for analyses described in Doxey et al.

* Code:
  -  `/docs/Rmarkdown.html`: initial version of analysis
  -  `/docs/Rmarkdown-V2.html`: updated analysis including negative controls
  -  `/docs/bacterialTranscriptomicsAnalysis.R`: code for analysis of MCAT/SPN/HFLU gene expression
  -  `/docs/hostResponseClassifier.R`: code for constructing host-response random-forest classifier
* Data required for analysis:
  - `d1-d5_fastp_dehost_krakenmerged_reads_normalized.txt`:  Taxonomic predictions for all samples
  - `gencode.v39.metadata.HGNC`: reference human transcriptome (https://www.gencodegenes.org/)
  - `/quants/`: gene expression data for all samples - please contact the authors
