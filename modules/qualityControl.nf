/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

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