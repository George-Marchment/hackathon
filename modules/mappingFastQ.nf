/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

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
            --runThreadN 16 \
            --outSAMtype BAM SortedByCoordinate \
            --outStd BAM_SortedByCoordinate \
            > ${sample}.bam
        """
}