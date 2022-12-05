/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

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