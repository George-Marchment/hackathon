/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

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