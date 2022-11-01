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
    publishDir 'data/annotation', mode: 'copy'

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
    publishDir 'data/index', mode: 'copy'

    input:
        path fasta
        path genomeAnnot
    
    output:
        path "GenomeDir"

    script:
        """
        STAR --runMode genomeGenerate --runThreadN 16 \
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
    publishDir 'data/map', mode: 'copy'

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
            --runThreadN 16 \
            --outSAMtype BAM SortedByCoordinate \
            --outStd BAM_SortedByCoordinate \
            > ${sample}.bam
        """
}


process indexBam {

	/*
	Indexation of the bam files
	@param : bam file
	@return : bai file
	*/

	label 'samtools'
    publishDir 'data/map', mode: 'copy'
    
    input:
    	path BAM

    output:
    	path "*.bai"

    script:  
		"""
		samtools index ${BAM}
		"""
}

process counting {
	/*
	Counting Reads 
	@param : all the bam files and the genome annotated
	@return : a txt file with all the information about the counting
	*/

	label 'featureCounts'
	publishDir 'data/counting', mode: 'copy'

	input:
		val alignedGene
		path annotedGenome

	output:
		path "countingReads.txt"

	script:
	"""
		featureCounts -T 4 -t gene -g gene_id -s 0 -a ${annotedGenome} -o countingReads.txt ${alignedGene}
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
    countingReads      : ${params.countingReads}

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
        pathA = genomeAnnotations()
    }else{
        pathA = Channel.fromPath('data/annotation/annotedGenome.gtf', checkIfExists : true, followLinks: false)
    }
    
    //Create gemome dir
    if (params.createGenome == true){
        pathGenomeDir = genomeIndex(pathG, pathA).collect()
    }else{
        pathGenomeDir = Channel.fromPath('data/index/GenomeDir/', checkIfExists : true, type: 'dir', followLinks: false).collect()
    }
    

    //Mapping 
    if (params.mapping == true){
        bam = mappingFastQ(fastq, pathGenomeDir)
    }else{
        bam =  Channel.fromPath('data/map/*bam', checkIfExists : true, followLinks: false)
    }
    
    //Index Bam
    if (params.indexBam == true){
        bai = indexBam(bam)
    }else{
        bai =  Channel.fromPath('data/map/*bai', checkIfExists : true, followLinks: false)
    }
    
    //Counting Reads
    if (params.countingReads == true){
    	bam.tolist().join(' ').view()
        //count = counting(bam.tolist().join(' '),pathA)
    //}else{
    	//count = 
    } 
}
