[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![Version 1.0.1](https://img.shields.io/badge/version-v1.0.0-blue)]()

# 2022 Hackathon Project 

**Master M2 AMI2B - Universite Paris Saclay**

Contributors :

* [Ambre Baumann](https://github.com/ambrebaumann)
* [Lindsay Goulet](https://github.com/Lindsay-Goulet)
* [George Marchment](https://github.com/George-Marchment)
* [Clémence Sebe](https://github.com/ClemenceS)

This repository contains :
    
* The Nextflow workflow which performs the ARN-seq analysis
* The scripts from which are created the Docker Images
* The final report
* Documentation on how to run the workflow on a VM

Supervisors : 

* Frédéric Lemoine
* Thomas Cokelaer

___

## Dependencies to install

Run these commands to install **Git**, **Singularity** and **Nextflow** onto the VM. 

> Note : The VM used during the project was Biopipes (it is feasible that some of these dependencies are already installed, in the case of an other VM)

```
sudo apt install -y git-all
conda activate
mamba install -y  singularity=3.6.3
conda install -c bioconda -y nextflow
conda update -y nextflow
```

___

## Run 

To run the workflow, run the following commands : 

```
conda activate
nextflow main.nf
```

When running the workflow, the user can specify which processes they want to be executed (by default they are all True in the *nextflow.config*).

Imagine the case, they don't want to run the *downloadFastqFiles* process to download the reads, but want to use the ones have been already downloaed (on a previous run), they would simply have run the following command :
 
```
nextflow main.nf --downloadFastqFiles False
```

> Note :  The reads which are used in this case are the ones found in data/seqs/ and following the pattern SRR*_{1,2}.fastq. If the user wants to use different fastq, they would have to change the *files* parameter in the same way.

___

<img src="pictures/paris-saclay.png">

