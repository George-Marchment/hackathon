/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

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
        STAR --runMode genomeGenerate --runThreadN ${params.nb_threads_star} \
        --genomeSAindexNbases 12 \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${genomeAnnot}
        """
}