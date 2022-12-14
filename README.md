[![MIT licensed](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE) [![Version 1.0.1](https://img.shields.io/badge/version-v1.0-yellow)]() [![Reproducibility](https://img.shields.io/badge/Crucial-Reproducibility-orange)]()


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
./run.sh
```

When running the workflow, the user can specify which processes they want to be executed (by default all the primary processes are all set to True, see *run.sh*).

Imagine the case, the user does not want to run the *downloadFastqFiles* process (to download the reads), but wants to use the ones which have already been downloaded (on a previous run), they would simply have to change the run.sh to :
 
```
#!/bin/bash

nextflow main.nf    \
    --downloadFastq false \
    --downloadGenome true \
    --downloadAnnotation true \
    --createGenome true \
    --doQuality false \
    --doTrimmomatic false \
    --getTrimmomatic false \
    --mapping true \
    --indexBam true \
    --countingReads true \
    --differentialAnalysis true
```

> Note :  The reads which are used in this case are the ones found in *data/seqs/* and following the pattern 'SRR*_{1,2}.fastq'. If the user wants to use different fastqs, they would have to change the *files* parameter int the *nextflow.config*.

> Note : In the same way, the user would have to modify the other parameters in the *nextflow.config* and *run.sh* for more flexibility, in the case they want to do a different analysis.

Here is a representation of the directed acyclic graph corresponding to the workflow :

<img src="pictures/dag.svg">

___

<img align="right" src="pictures/paris-saclay.png">


