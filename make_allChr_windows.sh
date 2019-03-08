## Generate 100KB windows per chromosome, 0-based

for i in {1..22}
do
    echo $i
    bedtools makewindows -g /seq/SNP_pipeline/hg19_forGATK/hg19.fasta.fai -w 100000 | awk '$1=="chr'${i}'"' | sed 's/chr//g' > ./hg19_chr${i}.100Kbp.windows.bed
done

