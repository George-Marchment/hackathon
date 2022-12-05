/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

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