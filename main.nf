process downloadGenome {

label 'downloadGenome'

output:
file "*.fa"

//Faut faire en sorte avec les tous les chromosones
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

publishDir '/data/seqs'
label 'downloadFastqFiles'
//Le publishdir ne fonctionne pas -> je ne sais pas pourquoie mais les fastq ce trouve dans le work en attendent

input:
val SRAID

output:
file "*"

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
	
	//Il faudra mettre des options pour l'utilisateur pour télécharger les données (qu'il faudra placé à des endroits précis)
	//Sinon juste pour l'analyse -> une option analyse qui ne télécharge pas les données (pour éviter de le faire tout le temps)

    //genomeIndex(downloadGenome(), genomeAnnotations())
   downloadFastqFiles(Channel.from(params.SRAID)).view()
	//Channel.from(params.SRAID).view()
}

