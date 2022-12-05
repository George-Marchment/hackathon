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
    --differentialAnalysis true