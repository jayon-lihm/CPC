# Cross Population Clustering (CPC)
Population-based Detection and Genotyping of Single Nucleotide Polymorphisms in Repeated Regions of the Genome
Jayon Lihm<sup>1</sup>, Vladimir Makarov<sup>2</sup>, and Seungtai (Chris) Yoon<sup>1</sup>  

1 Cold Spring Harbor Laboratory, Cold Spring Harbor, NY 11724  
2 Memorial Sloan Kettering Cancer Center, New York, NY 10065

(Manuscript in prep for submission)

Contact Info: jlihm@cshl.edu  

---

CPC is a SNP genotyping algorithm in large-scale population. The method can detect SNPs in small repeated regions that are difficult to be genotyped with conventional methods. The method requires "cluster" package in R.  

**1. Generate_AAP_forQuads.sh**  
This script processes bam file to pileup files to alternative allele count files that will be the groundwork for CPC. BAM files we used were by family (mostly quads consisting of mother, father, sibling, and proband). Thus the script starts from seprating a single family BAM file to four individual BAM files. Individual BAM files are then passed to python scripts to tabulate pileup files and determine reference and alternative alleles per position.  

**2. Make_AAP_withD30.sh**  
We deterine parameters based on high quality positions, that are positions with 30 or more read depth.  

**3. Filter_positions_with_small_samples.sh**  
In order for PAM to work on our data, we require at least 55 samples to have genotypes other than ref/ref. Thus we initially scan the number of samples per position and filter out the positions with less than 55 samples with data.

"make_allChr_windows.sh" needs to be run before step #3.  

**4. Clustering_AAP.R**  
This R script is the main body of our method. It collects data from all samples and generate a single matrix (of a given block of chromsome to reduce the size). Then PAM clustering is applied.  



