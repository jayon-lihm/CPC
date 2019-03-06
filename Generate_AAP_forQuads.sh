#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -o [stdout output path]
#$ -e [stderr output path]
#$ -l vf=40G

# This script processes one quad family (four members: mother, father, proband, and sibling)
# It is assumed that the reads for all four members are stored in a single BAM file
# Customize family ID and sample ID to suit your needs.

## Scripts Needed: Count_Alternative_Alleles.py, run_AAcount_forSSC.py

## submit as: qsub ./Generate_AAP_forQuads.sh

###########################################
## ${out_path}: root directory to store output files by family.
## ${bam_path}: directory to BAM file
## ${outFormat}: output as text file 
## ${refFasta}: path to reference FASTA file
## ${MapQual}: minimum mapping quality of reads required
## ${depth}: minimum number of reads required to calculate AAP
## ${chrstart}, ${chrend}: specify chromosome range



SAMTOOLS_PATH="/path/to/samtools/"

MapQual=20
depth=9 
chrstart=1
chrend=22

outFormat="txt"
refFasta="./reference/human_g1k_v37.fasta"

famID="A"
s1="mother"
s2="father"
s3="proband"
s4="sibling"

## Output paths
out_path="./"
fam_aap_path="${out_path}/${famID}/"
tmp_path="${fam_aap_path}/tmp/"

if [ ! -d ${fam_aap_path} ]; then mkdir -p ${fam_aap_path}; fi
if [ ! -d ${tmp_path} ]; then mkdir -p ${tmp_path}; fi

## Files
fam_bam_file=./quad_realigned_recalibrated_cs.bam
header_file=${tmp_path}${famID}_header.txt

echo "########JOB INFO###########"
echo "Running at: ${HOSTNAME}"
echo "JOB ID: ${JOB_ID}"
echo "TASK ID: ${SGE_TASK_ID}"
echo "START: $(date +"%D"), $(date +"%T")"
echo "FAMILY: ${famID}"

fam_start=$(date +%s)
fam_mem=($s1 $s2 $s3 $s4)
numMem=${#fam_mem[@]}

echo "####################################################"

# Export Header of Quad Bam
${SAMTOOLS_PATH}samtools view -H ${bam_path}quad_realigned_recalibrated_cs.bam > ${header_file}

for ((j=0; j<${numMem}; j++))
do
    sample=${fam_mem[$j]}
    samp_start=$(date +%s)

    aap_path="${fam_aap_path}${sample}/"
    if [ ! -d ${aap_path} ]; then mkdir -p ${aap_path}; fi

    rg_file=${tmp_path}${sample}_rglist.txt    
    sample_bam_file=${tmp_path}${sample}_realigned_recalibrated_cs.bam
    
    ## Find Read Group (RG) for a given sample
    grep "SM:${sample}" ${header_file} | awk '{print $2}' | cut -d':' -f2 > ${rg_file}
    ${SAMTOOLS_PATH}samtools view -bhR ${rg_file} ${fam_bam_file} > ${sample_bam_file}
    ${SAMTOOLS_PATH}samtools index ${sample_bam_file}
    
    echo "Sample: $sample"
    echo "Bam File Name: $sample_bam_file"

    if [ -e $sample_bam_file ]; then
	NumReadQ20=$(${SAMTOOLS_PATH}samtools view -q 20 -c $sample_bam_file)
    
	echo "Number of Reads"
	echo -e "ID\tNumReadsQ20"
	echo -e "${sample}\t${NumReadQ20}"
    
	if [ ${NumReadQ20} -ne 0 ]; then
	    for (( i=${chrstart}; i<=${chrend}; i++ ))
	    do
		ichar=$i
		if [ $i == 23 ]; then
		    ichar="X"
		elif [ $i == 24 ]; then
		    ichar="Y"
		fi
	  
		qpileupName=chr${ichar}_q${MapQual}_mapped.pileup
		qd_pileupName=chr${ichar}_q${MapQual}_d${depth}_mapped.pileup
		qd_altcntName=${qd_pileupName}.alt_count
		aapFile=${qd_altcntName}.aap
		refFile=${qd_altcntName}.ref
	  
 	  #******Q20 PILEUP GENERATION******#
		echo "..CHR${ichar} Q${MapQual} pileup"
		${SAMTOOLS_PATH}samtools mpileup -r ${ichar} -f ${refFasta} -q ${MapQual} ${sample_bam_file} > ${tmp_path}${qpileupName}
		gzip -f ${tmp_path}${qpileupName}
	  
      	  #for analysis, take positions with depth >= 9
		zcat ${tmp_path}${qpileupName}.gz | awk '$4>='$depth'' > ${tmp_path}${qd_pileupName}
		rm -f ${tmp_path}${qpileupName}.gz
		gzip -f ${tmp_path}${qd_pileupName}
	  
	  #******AAP GENERATION on Q20 PILEUP******#
		if [ -e ${tmp_path}${qd_pileupName}.gz ]; then
		    echo "..CHR${ichar} Q${MapQual} D${depth} alt count file"
		    python  ./Count_All_Alternative_Alleles.py ${tmp_path}${qd_pileupName}
		    rm -f ${tmp_path}${qd_pileupName}.gz
		    echo "...CHR${ichar} AAP/REF file"
		    python  ./run_AAcount_forSSC.py ${tmp_path}${qd_altcntName}
		    rm -f ${tmp_path}${qd_altcntName}.gz
		    cp -f ${tmp_path}${aapFile}.gz ${aap_path}. # copy aap and ref file to AAP count path
		    cp -f ${tmp_path}${refFile}.gz ${aap_path}.
		    rm -f ${tmp_path}${aapFile}.gz
		    rm -f ${tmp_path}${refFile}.gz
		else
		    echo "q${MapQual}_d${depth}_pileup doesn't exist."                                                                              
		fi # if ${qd_pileupName}.gz exists 
	    done # for chr 1 to 24
	else # if there is a read
	    echo "There is no read in the bam"
	fi # if there is a read
	rm -f ${sample_bam_file}*
    else # if bam exist
	echo "${sample_bam_file} does not exist"
    fi # if bam exists
    
    echo "END: $(date +"%D"), $(date +"%T")"
    samp_end=$(date +%s)
    echo "IT TOOK $(( ${samp_end} - ${samp_start} )) SECONDS FOR THIS SAMPLE."
    echo "####################################################"

done #sample 1-4	    

fam_end=$(date +%s)
echo "***IT TOOK $(( ${fam_end} - ${fam_start} )) SECONDS FOR THIS FAMILY.***"

