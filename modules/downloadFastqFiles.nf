/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
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
        tuple val(SRAID), path("*1.fastq"), path ("*2.fastq")  

    script:
        """
        fasterq-dump ${SRAID}
        """
}