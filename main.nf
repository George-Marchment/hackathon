/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

//Import the different workflows
include { RNA_SEQ } from './workflows/rna-seq'

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

    RNA_SEQ()

}
