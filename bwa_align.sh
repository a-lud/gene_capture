#!/bin/bash

# Arguments
# $1 = reference
# $2 = R1
# $3 = R2 
# $4 = outBAM = all aligned
# $5 = outBAM_flagstat = statistics
# $6 = outBAMFILT = filtered output
# $7 = threads

## Align read using bwa

## Check if index file exists
if [[ ! -f $(dirname ${1})/$(basename ${1} ]).bwt ]]; then
    bwa index ${1}
fi

## Align reads
bwa mem -B 2 -M -t ${7} ${1} ${2} ${3} | samtools view -b -@ ${7} -o ${4} > ${4}

## Statistics on bam
samtools flagstat -@ ${7} ${4} > ${5}
mosdepth -t ${7} $(dirname ${4})/$(basename ${5} .flagstat) ${4}

## Filtering raw alignments
samtools view -b -@ ${7} -F 4 ${4} |
samtools sort -O BAM -@ ${7} -o ${6}

## Index the bam
samtools index -@ ${7} ${6}