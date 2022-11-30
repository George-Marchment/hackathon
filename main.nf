/**
Nextflow workflow performing an RNA-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

/**
Notes :
    - We have chosen to leave the full workflow in one file, for readibility reasons
    - Structure of workflow (code) : define the list of processes then define the main of the workflow 
    - 
**/

process downloadFastqFiles {
    /*
    Download the fastq files with the tool fasterq-dump
    Inputs : a list of SRAID
    Outputs : the fastqs 
    */

    publishDir 'data/seqs', mode: 'copy'
    label 'SRA'

    input:
        val SRAID

    output:
        tuple val(SRAID), path("*1.fastq.gz"), path ("*2.fastq.gz")  

    script:
        """
        fasterq-dump ${SRAID}
	    gzip ${SRAID}*.fastq
        """
}

process downloadGenome {
    /*
    Download the genome of a list of chromosomes
    Inputs : ID chromosome
    Outputs : the file of the chromosome download
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
    Take in input the 'genome' for each chromosome and concatenate it into one file
    Inputs : The chromosones
    Outputs : the 'global' genome
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
    Outputs : the gtf file with the annotation
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
    Inputs : the fasta file (genome) and the genome annotation (gft)
    Outputs : the genome repository
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
        STAR --runMode genomeGenerate --runThreadN ${task.cpus} \
        --genomeSAindexNbases 12 \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${genomeAnnot}
        """
}

process qualityControl {
    /*
    Verify the reads quality
    Inputs : the reads (name of the sample - R1 - R2)
    Outputs : the html file of the quality control
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
    Inputs : the reads (name of the sample - R1 - R2)
    Outputs : the reads trimmed 
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
    Mapping the read to genome with STAR
    Inputs : the reads R1 and R2 and the genome dir
    Outputs : file bam aligned
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
        STAR --outSAMstrandField intronMotif \
            --outFilterMismatchNmax 4 \
            --outFilterMultimapNmax 10 \
            --genomeDir ${GenomeDir} \
            --readFilesIn ${R1} ${R2} \
            --runThreadN ${task.cpus} \
            --outSAMtype BAM SortedByCoordinate \
            --outStd BAM_SortedByCoordinate \
            > ${sample}.bam
        """
}


process indexBam {
	/*
	Indexation of the bam files
	Inputs : bam file
	Outputs : bai file
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
	Inputs : all the bam files and the genome annotated
	Outputs : a txt file with all the information about the counting
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

process test_george {

    publishDir 'George/', mode: 'copy'
    
    output:
        stdout

    script:
        """
        echo 'george' > test.txt
        """
}


process diffAnalysis {
    /*
	Differential analysis
	Inputs : the feature count file and the metadata file
	Outputs : analysis repository (a txt file with differentially expressed genes, histogram of p-values)
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
        Rscript analyse.R
        """
}


log.info """\

===================================
 H A C K A T H O N  W O R K F L O W
===================================

Written by :

    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

Options selected :

    - downloadFastq        : ${params.downloadFastq}
    - downloadGenome       : ${params.downloadGenome}
    - downloadAnnotation   : ${params.downloadAnnotation}
    - createGenome         : ${params.createGenome}
    - doQuality            : ${params.doQuality}
    - doTrimmomatic        : ${params.doTrimmomatic}
    - getTrimmomatic       : ${params.getTrimmomatic}
    - mapping              : ${params.mapping}
    - indexBam             : ${params.indexBam}
    - countingReads        : ${params.countingReads}
    - differentialAnalysis : ${params.differentialAnalysis}

 """

workflow {

    res = test_george()
    res.view()

    counting(bam.toList(),pathA)
    
}
