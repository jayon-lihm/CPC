#!/bin/sh
#$ -v PATH
#$ -o [stdout out]
#$ -e [stderr out]
#$ -l m_mem_free=4G
#$ -l tmp_free=2G

## Count number of samples who have AAP>0 (d30) at the positions
## Record if the total count is greater than 55

## chrom="chr2"
## cutoff=55
## window_bed="/seq/jlihm/AAP/greyAAP2/hg19_${chrom}.100Kbp.windows.bed"
## maxJob=$(cat ${window_bed} | wc -l)
## qsub -t 1-$maxJob -tc 200 ./Filter_positions_with_small_samples.sh ${chrom} ${window_bed} ${cutoff}

chrom=$1
window_bed=$2
cutoff=$3

famfile="./familyinfo.parents.txt"

chrbed="./hg19_${chrom}.100Kbp.windows.sti${SGE_TASK_ID}.bed"
outdir="./${chrom}/"
if [ ! -d $outdir ]; then mkdir -p $outdir; fi

resfile="${outdir}${chrom}_d30AAP_sampleCounts.sti${SGE_TASK_ID}.bedg.gz"

sed -n "${SGE_TASK_ID}p" ${window_bed} > ${chrbed}

count=0
echo "chrom: ${chrom}"
echo "*START: $(date +"%D"), $(date +"%T")"
while read line
do
    famID=$(echo $line | awk '{print $1}')
    sampleID=$(echo $line | awk '{print $2}')

    d30aap="./${famID}/${sampleID}/${chrom}_q20_d30.aap.bed.gz"
    d30aapBG="./${famID}/${sampleID}/tmp_sti${SGE_TASK_ID}.${famID}_${sampleID}.${chrom}_q20_d30.aap.bedg.gz"
    
    bedtools intersect -a ${chrbed} -b ${d30aap} -wa -wb | awk 'BEGIN{OFS="\t"}{print $4, $5, $6, 1}' | gzip -c > ${d30aapBG}
    numAAP=$(zcat ${d30aapBG} | wc -l)
    echo $count $famID $sampleID $numAAP

    if [ $numAAP -gt 0 ]; then
	(( count=$count+1 ))
    fi

    if [ $count -eq 1 ]; then
	afile=${d30aapBG}
    elif [ $count -gt 1 ]; then
	if [ $numAAP -gt 0 ]; then
	    bfile=${d30aapBG}
	    new_afile="./tmp_${chrom}_sti${SGE_TASK_ID}.d30AAP_count${count}.bedg.gz"
	    bedtools unionbedg -i ${afile} ${bfile} | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4+$5}' | gzip -c > ${new_afile}
	    rm $afile
	    rm $bfile
	    afile=${new_afile}
	fi
    fi
done < $famfile

if [ $count -le $cutoff ]; then
    echo -n "" | gzip -c  > $resfile
elif [ $count -gt $cutoff ]; then
    zcat $afile | awk '$4>"'$cutoff'"' | gzip -c > $resfile
    rm ${new_afile}
    rm ${chrbed}
fi

echo "END: $(date +"%D"), $(date +"%T")"
