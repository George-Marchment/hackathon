/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

//Import the different modules
include { downloadFastqFiles } from '../modules/downloadFastqFiles'
include { downloadGenome } from '../modules/downloadGenome'
include { concatenateGenome } from '../modules/concatenateGenome'

//For more info on the individuals processes, see modules folder


//Definition of the subworkflow qualityControl
workflow downloadFiles{
    main:
        //Download Fastq files
        if (params.downloadFastq == true){
            fastq = downloadFastqFiles(Channel.from(params.SRAID))
        }else{
            //If the user doesn't want to download them, check if they exist
            fastq = Channel.fromFilePairs(params.files, checkIfExists : true, flat: true, followLinks: false)
        }
        
        //Download and concatenate genomes
        if (params.downloadGenome == true){
            pathG = concatenateGenome(downloadGenome(Channel.from(params.CHR)).toList())
        }else{
            //If the user doesn't want to download them, check if they exist
            pathG = Channel.fromPath(params.referenceGenome, followLinks: false)  
        }

    emit:
        //-------
        //Outputs
        //-------
        fastq
        pathG
}