#!/bin/bash

nextflow main.nf    \
    --downloadFastq true \
    --downloadGenome true \
    --downloadAnnotation true \
    --createGenome true \
    --doQuality false \
    --doTrimmomatic false \
    --getTrimmomatic false \
    --mapping true \
    --indexBam true \
    --countingReads true \
    --differentialAnalysis true \

    --SRAID ["SRR628589", "SRR628588", "SRR628587", "SRR628586", "SRR628585", "SRR628584", "SRR628583", "SRR628582"] \
    --CHR ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"] \

    --files 'data/seqs/SRR*_{1,2}.fastq' \
    --referenceGenome 'data/genome/referenceGenome.fa' \
    --annotatedGenome 'data/annotation/annotedGenome.gtf' \
    --GenomeDir 'data/index/GenomeDir/' \
    --trimmoFiles 'data/trimmomatic/*{1,2}P.fastq' \
    --mappingFiles 'data/map/*bam' \
    --baiFiles 'data/map/*bai'