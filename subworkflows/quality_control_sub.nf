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
include { qualityControl } from '../modules/qualityControl'
include { trimming } from '../modules/trimming'

//For more info on the individuals processes, see modules folder


//Definition of the subworkflow qualityControlSub
workflow qualityControlSub{
    main:
        //Quality Control
        if (params.doQuality == true){
            qualityControl(fastq)
        }

        //Trimmomatic
        if (params.doTrimmomatic == true){
            //In the case the user wants to trim the reads (in the case they are poor quality)
            if (params.getTrimmomatic == true){
                //If the user wants to run the process Trimmomatic
                new_fastq = trimming(fastq)
            }else{
                //If not, checks in the files exist
                new_fastq = Channel.fromFilePairs(params.trimmoFiles, checkIfExists : true, flat: true, followLinks: false)
            }   
        }else{
            //In the case the user doesn't want to trim the reads (in the case they are good quality)
            new_fastq = fastq
        }

    emit:
        //-------
        //Outputs
        //-------
        new_fastq

}