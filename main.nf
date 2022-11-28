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
        tuple val(SRAID), path ("*1.fastq.gz"), path ("*2.fastq.gz")  

    script:
        """
        fasterq-dump ${SRAID}
	gzip ${SRAID}*.fastq
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
    /*
    Take in input all the genome for each chromosome and concatenate it into one file
    @param : all the chr genome 
    @return : the 'global' genome
    */

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

process qualityControl {
    /*
    Verify the reads quality
    @param : the reads (name of the sample - R1 - R2)
    @return : the html file of the quality control
    */

    label 'fastqc'
    publishDir "data/qualityGraph", mode: 'copy'

    input:
        tuple val(SRAID), path(R1), path(R2)

    output:
        path "*.html"

    script:
        """
        fastqc ${R1} ${R2}
        """
}

process trimming {
    /*
    Remove the adaptators
    @param : the reads (name of the sample - R1 - R2)
    @return : the reads trimmed 
    */

    label 'trimmomatic'
    publishDir 'data/trimmomatic/', mode: 'copy'

    input:
        tuple val(SRAID), path(R1), path(R2)

    output:
        tuple val(SRAID), path("*1P.fastq"), path("*2P.fastq")

    script:
        """
        trimmomatic PE ${R1} ${R2} -baseout \
        ${SRAID}.fastq  LEADING:20 TRAILING:20 MINLEN:50
        """
}

process mappingFastQ {
    /*
    Mapping the read with STAR
    @param : the reads R1 and R2 and the genome dir
    @return : file bam aligned
    */
    label 'STAR'
    publishDir 'data/map', mode: 'copy'

    input:
        tuple val(sample), path(R1), path(R2)
        path GenomeDir

    output: 
        path "*.bam"

    script:
        """
        STAR --outFilterMismatchNmax 4 \
            --outFilterMultimapNmax 10 \
            --genomeDir ${GenomeDir} \
            --readFilesIn ${R1} ${R2} \
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
		path alignedGene
		path annotedGenome

	output:
		path "countingReads.txt"

	script:
        """
        featureCounts -p -T 4 -t gene -g gene_id -s 0 -a ${annotedGenome} -o countingReads.txt ${alignedGene}
        """
}

process analyseDifferentielle {
    /*
	Differential analysis
	@param : the feature count file and the metadata file
	@return : analysis repository (a txt file with differentially expressed genes, histogram of p-values)
	*/

	label 'deseq2'
	publishDir 'data/differentialAnalysis', mode: 'copy'

	input:
		path countingReads
		path metadata

	output:
		path "DEgenes.txt"
        path "*.png"

	script:
        """
        # Chargement des données
        countData <- read.table(header=TRUE, row.names = 1, ${countingReads})
        countData <- countData[,6:13]
        metadata <- metadata <- read.table(header=TRUE, sep="\t", ${metadata})

        cond <- c()
        for (ind in colnames(countData)){
            cond = c(cond, metadata[which(metadata8[,1]==ind),2] == "yes")
        }
        cond <- factor(cond)

        # Filtrage des gènes
        countData <- countData[-which(rowSums(countData8) == 0),]

        mutant = metadata[which(metadata[,2]==1),1]
        nonMutant = metadata[-which(metadata[,2]==1),1]

        for (indMutant in mutant){
            for (indNonMutant in nonMutant) {
                countData <- countData[countData[,indMutant]!=0 | countData[,indNonMutant]!=0,]
            }
        }

        # Analyse différentielle
        dds <- DESeqDataSetFromMatrix(countData=countData, colData=DataFrame(cond), design=~cond)
        dds <- DESeq(dds)
        res <- results(dds)

        png("histogram.png")
        hist(res$pvalue, main='Histogram of pvalues')
        dev.off()

        DEgenes <- res[which(res$padj<0.05),]

        rld <- rlog(dds)

        png("PCA.png")
        plotPCA(rld, intgroup='mutant') + geom_text(aes(label=name),vjust=2)
        dev.off()

        write.table(DEgenes, "DEgenes.txt", sep="\t")
        """

}


log.info """\

 H A C K A T H O N  P I P E L I N E
===================================

    downloadFastq      : ${params.downloadFastq}
    downloadGenome     : ${params.downloadGenome}
    downloadAnnotation : ${params.downloadAnnotation}
    createGenome       : ${params.createGenome}
    doQuality          : ${params.doQuality}
    doTrimmomatic      : ${params.doTrimmomatic}
    getTrimmomatic     : ${params.getTrimmomatic}
    mapping            : ${params.mapping}
    indexBam           : ${params.indexBam}
    countingReads      : ${params.countingReads}

 """

workflow {
    //Sinon juste pour l'analyse -> une option analyse qui ne télécharge pas les données (pour éviter de le faire tout le temps)

    //il faudrait mettre en param : le nb de coeurs et autres parametres generaux -> (George) je ne suis pas sûr qu'il ya besoin comme on fait tourné en local sur le VM

    
    //Download Fastq files
    if (params.downloadFastq == true){
        fastq = downloadFastqFiles(Channel.from(params.SRAID))
    }else{
        fastq = Channel.fromFilePairs(params.files, checkIfExists : true, flat: true, followLinks: false)
    }
    
    //Download genome and annotation
    if (params.downloadGenome == true){
        pathG = concatenateGenome(downloadGenome(Channel.from(params.CHR)).toList())
    }else{
        pathG = Channel.fromPath(params.referenceGenome, checkIfExists : true, followLinks: false)  
    }
    
    if (params.downloadAnnotation == true){
        pathA = genomeAnnotations()
    }else{
        pathA = Channel.fromPath(params.annotatedGenome, checkIfExists : true, followLinks: false)
    }
    
    //Create genome dir
    if (params.createGenome == true){
        pathGenomeDir = genomeIndex(pathG, pathA).collect()
    }else{
        pathGenomeDir = Channel.fromPath(params.GenomeDir, checkIfExists : true, type: 'dir', followLinks: false).collect()
    }

    //Quality Control
    if (params.doQuality == true){
        qualityControl(fastq)
    }

    //Trimmomatic
    if (params.doTrimmomatic == true){
        if (params.getTrimmomatic == true){
            new_fastq = trimming(fastq)
        }else{
            new_fastq = Channel.fromFilePairs(params.trimmoFiles, checkIfExists : true, flat: true, followLinks: false)
        }   
    }else{
        new_fastq = fastq
    }

    //Mapping 
    if (params.mapping == true){
        bam = mappingFastQ(new_fastq, pathGenomeDir)
    }else{
        bam =  Channel.fromPath(params.mappingFiles, checkIfExists : true, followLinks: false)
    }
    
    //Index Bam
    if (params.indexBam == true){
        bai = indexBam(bam)
    }else{
        bai =  Channel.fromPath(params.baiFiles, checkIfExists : true, followLinks: false)
    }
    
    //Counting Reads
    if (params.countingReads == true){
        count = counting(bam.toList(),pathA)
    }else{
    	count = Channel.fromPath('data/counting/countingReads.txt', checkIfExists : true, followLinks: false)
    } 

    //Differential analysis
    
}
