/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

process counting {
	/*
	Counting Reads 
	Inputs : all the bam files and the genome annotated
	Outputs : a txt file with all the information about the counting
	*/

	label 'featureCounts'
	publishDir 'data/counting', mode: 'copy'

	input:
		path alignedGene
		path annotedGenome

	output:
		path "countingReads.txt"

	script:
        """
        featureCounts -p -T 4 -t gene -g gene_id -s 0 -a ${annotedGenome} -o countingReads.txt ${alignedGene}
        """
}