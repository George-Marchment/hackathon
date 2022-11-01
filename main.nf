process downloadFastqFiles {
    /*
    Download the fastq files with the tool fasterq-dump
    @param : a list of SRAID
    @return : all the fastq download
    */

    publishDir 'data/seqs', mode: 'copy'
    label 'SRA'

    input:
        val SRAID

    output:
        tuple val(SRAID), path ("*1.fastq"), path ("*2.fastq")  

    script:
        """
        fasterq-dump ${SRAID}
        """
}

process downloadGenome {
    /*
    Download the genome of a list of chromosomes
    @return : the file of the chromosome download
    */

    label 'downloadGenome'

    publishDir 'data/genome', mode: 'copy'

    input : 
        val CHR

    output:
        path "*.fa"

    script:
        """
        wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${CHR}.fa.gz
        gunzip -c *.fa.gz > referenceGenomeChr${CHR}.fa
        rm *gz
        """
}

process concatenateGenome {

    label 'concatenateGenome'

    publishDir 'data/genome', mode: 'copy'

    input :
        path(genome)
    
    output : 
        path "*.fa"

    script :
    """
    cat ${genome} > referenceGenome.fa
    """
}

process genomeAnnotations {
    /*
    Download the genome Annotation
    @return : the gtf file with the annotation
    */

    label 'genomeAnnotations'

    output:
        path "*.gtf"

    script:
        """
        wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
        gunzip -c *.gtf.gz > annotedGenome.gtf
        rm *gz
        """
}

process genomeIndex {
    /*
    Create the repository with the genome files created 
    @param : the fasta file (chromosome) and the genome annotation (gft)
    @return : the genome repository
    */
    label 'STAR'

    input:
        path fasta
        path genomeAnnot
    
    output:
        path "GenomeDir"

    script:
        """
        STAR --runMode genomeGenerate --runThreadN 4 \
        --genomeSAindexNbases 12 \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${genomeAnnot}
        """
}

process mappingFastQ {
    /*
    Maaping the read with STAR
    @param : the reads R1 and R2 and the gemome dir
    @return : file bam aligned
    */
    label 'STAR'

    input:
        tuple val(sample), path(fastq1), path(fastq2)
        path GenomeDir

    output: 
        file "*.bam"

    script:
        """
        STAR --outFilterMismatchNmax 4 \
            --outFilterMultimapNmax 10 \
            --genomeDir ${GenomeDir} \
            --readFilesIn ${fastq1} ${fastq2} \
            --runThreadN 4 \
            --outSAMtype BAM SortedByCoordinate \
            --outStd BAM_SortedByCoordinate \
            > ${sample}.bam
        """
}

log.info """\

 H A C K A T H O N  P I P E L I N E
===================================

    downloadFastq      : ${params.downloadFastq}
    downloadGenome     : ${params.downloadGenome}
    downloadAnnotation : ${params.downloadAnnotation}
    createGenome       : ${params.createGenome}
    mapping            : ${params.mapping}

 """

workflow {
    //Il faudra mettre des options pour l'utilisateur pour télécharger les données (qu'il faudra placer à des endroits précis)
    //Sinon juste pour l'analyse -> une option analyse qui ne télécharge pas les données (pour éviter de le faire tout le temps)

    //il faudrait mettre en param : le nb de coeurs et autres parametres generaux -> (George) je ne suis pas sûr qu'il ya besoin comme on fait tourné en local sur le VM

    
    //Download Fastq files
    if (params.downloadFastq == true){
        fastq = downloadFastqFiles(Channel.from(params.SRAID))
    }else{
        fastq = Channel.fromFilePairs('data/seqs/SRR*_{1,2}.fastq', checkIfExists : true, flat: true, followLinks: false)
    }
    
    //Download genome and annotation
    if (params.downloadGenome == true){
        pathG = concatenateGenome(downloadGenome(Channel.from(params.CHR)).toList())
    }else{
        pathG = Channel.fromPath('data/genome/referenceGenome.fa', checkIfExists : true, followLinks: false)  
    }
    
    if (params.downloadAnnotation == true){
        genomeAnnotations()
        pathA = genomeAnnotations.out
    }else{
        pathA = Channel.fromPath('work/**/annotedGenome.gtf', checkIfExists : true, followLinks: false)
    }
    
    //Create gemome dir
    if (params.createGenome == true){
        genomeIndex(pathG, pathA)
        pathGenomeDir = genomeIndex.out
    }else{
        pathGenomeDir = Channel.fromPath('work/**/GenomeDir/', checkIfExists : true, type: 'dir', followLinks: false)
    }
    

    //Mapping 
    if (params.mapping == true){
        bam = mappingFastQ(fastq, pathGenomeDir)
    //}else{
        //channel pour trouver les fichiers si deja telecharges
    }
}