#!/usr/bin/env bash

module load seqkit/0.8.1

READS="${FASTDIR}/path/to/supertranscripts-fastas-dir"
MAP="${FASTDIR}/path/to/trinotate-key-value-annotation-dir"
OUT="${FASTDIR}/path/to/output/dir"

mkdir -p ${OUT}

FILES=$(find ${MAP} -type f -name '*.txt' -exec basename {} .annotationMap.txt \;)
TYPE="pep cds"

for f in ${FILES}; do

    for t in ${TYPE}; do
        seqkit replace -p " len:.*" -r "" ${READS}/${f}/${f}.Trinity.SuperTrans.fasta.transdecoder.${t} | \
        seqkit replace -p "\sGENE\\.TRINITY.*ORF\s" -r ' ' | \
        seqkit replace -p "(^.*)\s" -r "{kv}::" -K -k ${MAP}/${f}.annotationMap.txt > ${OUT}/${f}.SuperTrans.anno.${t}
    done

done