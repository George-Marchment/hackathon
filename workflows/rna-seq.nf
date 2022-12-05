/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

//Import the different subworkflows
include { downloadFiles } from '../subworkflows/download_files'
include { qualityControlSub } from '../subworkflows/quality_control_sub'

//Import the different modules
include { genomeAnnotations } from '../modules/genomeAnnotations'
include { genomeIndex } from '../modules/genomeIndex'
include { mappingFastQ } from '../modules/mappingFastQ'
include { indexBam } from '../modules/indexBam'
include { counting } from '../modules/counting'
include { diffAnalysis } from '../modules/diffAnalysis'

//For more info on the individuals processes and subworkflows, see modules and subworflows folder


//Definition of the workflow
workflow RNA_SEQ {
    
    //***************************
    //SUBWORKFLOW DOWNLOAD FILES
    //***************************
    //Call the subworkflow "downloadFiles" which can be found in the folder subworkflows
    downloadFiles()
    //Retrieve the outputs of downloadFiles
    fastq = downloadFiles.out.fastq
    pathG = downloadFiles.out.pathG
    
    
    //Download genome annotation
    if (params.downloadAnnotation == true){
        pathA = genomeAnnotations()
    }else{
        //If the user doesn't want to download it, check if it exists
        pathA = Channel.fromPath(params.annotatedGenome, checkIfExists : true, followLinks: false)
    }
    
    //Create genome dir
    if (params.createGenome == true){
        pathGenomeDir = genomeIndex(pathG, pathA).collect()
    }else{
        //If the user doesn't want to create it, check if it exists
        pathGenomeDir = Channel.fromPath(params.GenomeDir, checkIfExists : true, type: 'dir', followLinks: false).collect()
    }

    
    //***************************
    //SUBWORKFLOW QUALITY CONTROL
    //***************************
    //Call the subworkflow "qualityControlSub" which can be found in the folder subworkflows
    qualityControlSub()
    //Retrieve the outputs of downloadFiles
    new_fastq = qualityControlSub.out.new_fastq

    

    //Mapping 
    if (params.mapping == true){
        bam = mappingFastQ(new_fastq, pathGenomeDir)
    }else{
        //In the case the user doesn't want to run the mapping process -> check that the files exist 
        bam =  Channel.fromPath(params.mappingFiles, checkIfExists : true, followLinks: false)
    }
    
    //Index Bam
    if (params.indexBam == true){
        bai = indexBam(bam)
    }else{
        //In the case the user doesn't want to run the indexBam process -> check that the files exist 
        bai =  Channel.fromPath(params.baiFiles, checkIfExists : true, followLinks: false)
    }
    
    //Counting Reads
    if (params.countingReads == true){
        count = counting(bam.toList(),pathA)
    }else{
        //In the case the user doesn't want to run the countingReads process -> check that the file exist 
    	count = Channel.fromPath('data/counting/countingReads.txt', checkIfExists : true, followLinks: false)
    } 

    //Differential analysis
    if (params.differentialAnalysis == true){
        diffAnalysis(Channel.fromPath('Rscript/analyse.R', checkIfExists : true, followLinks: false), count, Channel.fromPath('Rscript/metadata.txt', checkIfExists : true, followLinks: false))
    }

}
    