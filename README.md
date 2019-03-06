# Cross Population Clustering (CPC)
Population-based Detection and Genotyping of Single Nucleotide Polymorphisms in Repeated Regions of the Genome
Jayon Lihm<sup>1</sup>, Vladimir Makarov<sup>2</sup>, and Seungtai (Chris) Yoon<sup>1</sup>  

1 Cold Spring Harbor Laboratory, Cold Spring Harbor, NY 11724  
2 Memorial Sloan Kettering Cancer Center, New York, NY 10065

(Manuscript in prep for submission)

---

CPC is a SNP Genotyping algorithm in large-scale population, specialized in repeated regions. The method can detect SNPs in small repeated regions that are hard to be genotyped with conventional methods. The main algorithm of CPC works on the proportion of alternative allele per position, called "Alternative Allele Proportion (AAP)". Assuming bi-allelic SNP models, we first identify the reference allele and alternative allele per position. From pileup files generated from Samtools, the number of reads containing A, G, C, T, and DEL are counted per position. We then identify reference allele and its corresponding counts, alternative allele and its count. When there are multiple alternative alleles mapped to this position, we consider the alternative allele with the largest count as the alternative allele for the position. If two or more alternative alleles have tie in counts, the script randomly selects the alternative allele out of them. We filter out low coverage regions for quality control; Positions with at least 9 reads of mapping quality 20 or greater are selected.  

For regular bi-allelic SNPs in two copy regions, the distribution of AAP from unrelated individuals is assumed to have three peaks at 0 for ref/ref (g0), 0.5 for ref/alt (g1), and 1 for alt/alt (g2) genotypes. For repeted regions, however, we observe that peaks are located at unexpected values, for example 0.25 and 0.5, while maintaining the tri-modal distribution. Our algorithm uses PAM clustering to determine two boundaries between g0 and g1, g1 and g2 and assigns genotypes based on the two thresholds to accomodate complicated AAP structure under repeated regions.

PAM clustering: R library "cluster"

**1. Generate_AAP_forQuads.sh**  
This script processes bam file to pileup files to alternative allele count files that will be the groundwork for CPC. BAM files we used were by family (mostly quads consisting of mother, father, sibling, and proband). Thus the script starts from seprating a single family BAM file to four individual BAM files. Individual BAM files are then passed to python scripts to tabulate pileup files and determine reference and alternative alleles per position and count the corresponding number of reads containing reference and alternative alleles.  


**2. Make_AAP_withD30.sh**  
We deterine boundaries based on high quality positions, that are positions with 30 or more read depth. This script generates a list of AAP values with depth 30 or greater.  

**3. Filter_positions_with_small_samples.sh**  
In order for PAM to work on AAP values, we require at least 50 samples to be in g1 and at least 5 samples in g2 genotype. Thus we initially scan the number of samples with AAP>0 per position and filter out the positions with less than or equal to 55 samples with AAP>0. PAM would not work well on these positions due to small number of samples.

"make_allChr_windows.sh" needs to be run before step #3.  

**4. Clustering_AAP.R**  
This R script is the main body of our method. It collects AAP values from all samples and generate a single AAP matrix (of a given block of chromsome to reduce the size). Then PAM clustering runs to determine the lower and upper thresholds that separates g0, g1, and g2; thus the number of clusters is forced to be 3. PAM is applied on four combinations of two conditions: Include/exclude 0 AAPs, trim 5% on both tails or do not trim. We compare the silhouette values from the four clustering results, and consider the one with the largest silhouette value as the **best** clustering. LT and UT are determined from the best clustering result.  

