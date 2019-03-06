#$ -v PATH
#$ -o [stdout output]
#$ -e [stderr output]
#$ -l m_mem_free=4G
#$ -l tmp_free=2G

#### INPUT: Alternative allele count AAP and REF files
#### OUTPUT: bed file with aap for d30 positions (ref/aap)

## qsub -t 1-22 ./run_collect_allChr_allAAP.sh

famfile="./familyinfo.parents.txt" ## first column: famID, second column: sample ID

chrom=${SGE_TASK_ID}

while read line
do
    famID=$(echo $line | awk '{print $1}')
    sampleID=$(echo $line | awk '{print $2}')
    aap_path="./${famID}/${sampleID}/" ## Specify
    out_path="./" ## Specify

    aapfile="${aap_path}chr${chrom}_q20_d9_mapped.pileup.alt_count.aap.gz"
    reffile="${aap_path}chr${chrom}_q20_d9_mapped.pileup.alt_count.ref.gz"
    
    d30aap="${out_path}chr${chrom}_q20_d30.aap.bed.gz"
    d30ref="${out_path}chr${chrom}_q20_d30.ref.bed.gz"

    zcat ${aapfile} | sed '1d' | awk 'BEGIN{OFS="\t"}{if($12=="SNP" && $13>=30) print $1, $2-1, $2, $14/$13}' | gzip -c > ${d30aap}
    zcat ${reffile} | sed '1d' | awk 'BEGIN{OFS="\t"}{if($13>=30) print $1, $2-1, $2}' | gzip -c > ${d30ref}
done < $famfile


