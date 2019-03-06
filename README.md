# Cross Population Clustering (CPC)
Population-based Detection and Genotyping of Single Nucleotide Polymorphisms in Repeated Regions of the Genome
Jayon Lihm<sup>1</sup>, Vladimir Makarov<sup>2</sup>, and Seungtai (Chris) Yoon<sup>1</sup>  

1 Cold Spring Harbor Laboratory, Cold Spring Harbor, NY 11724  
2 Memorial Sloan Kettering Cancer Center, New York, NY 10065

---

CPC is a SNP Genotyping algorithm in large-scale population, specialized in repeated regions. The method can detect SNPs in small repeated regions that are hard to be genotyped with conventional methods.
The main algorithm of CPC works on the proportion of alternative allele per position, called "Alternative Allele Proportion (AAP)".
Assuming bi-allelic SNP models, we first identify the reference allele and alternative allele per position.
From pileup files generated from Samtools, the number of reads containing A, G, C, T, and DEL are counted per position.
We then identify reference allele and its corresponding counts, alternative allele and its count. When there are multiple alternative alleles mapped to this position, we consider the alternative allele with the largest count as the alternative allele for the position. If two or more alternative alleles have tie in counts, the script randomly selects the alternative allele out of them.

We filter out low coverage regions are filtered; at least 9 reads should have mapping quality 20 or greater at a position, in order to assume the continuity of AAP values.

