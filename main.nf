process downloadGenome {

label 'downloadGenome'

output:
file "*.fa"

script:
"""
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
gunzip -c *.fa.gz > referenceGenome.fa
rm *gz
"""
}


process genomeAnnotations {

label 'genomeAnnotations'

output:
file "*.gtf"

script:
"""
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
gunzip -c *.gtf.gz > annotedGenome.gtf
rm *gz
"""
}


process downloadFastqFiles {

label 'downloadFastqFiles'

input:
val SRAID

output:
file "*.sra"

script:
"""
fasterq-dump ${SRAID}
"""
}


process genomeIndex {
    label 'STAR'

    input:
    file fasta
    file genomeAnnot

    script:
    """
    STAR --runMode genomeGenerate --runThreadN 4 \
    --genomeSAindexNbases 12 \
    --genomeFastaFiles ${fasta} \
    --sjdbGTFfile ${genomeAnnot}
    """
}

workflow {
    genomeIndex(downloadGenome(), genomeAnnotations())
}

