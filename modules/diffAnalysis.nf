/**
Nextflow workflow performing an ARN-seq Analysis created by :
    - Ambre Baumann (https://github.com/ambrebaumann)
    - Lindsay Goulet (https://github.com/Lindsay-Goulet)
    - George Marchment (https://github.com/George-Marchment)
    - Clémence Sebe (https://github.com/ClemenceS)

This workflow was created in the context of the 2022 Hackathon project (master M2 AMI2B) 
This project was surpervised by Frédéric Lemoine & Thomas Cokelaer
**/

process diffAnalysis {
    /*
	Differential analysis
	Inputs : the feature count file and the metadata file
	Outputs : analysis repository (a txt file with differentially expressed genes, histogram of p-values)
	*/

	label 'deseq2'
	publishDir 'data/differentialAnalysis', mode: 'copy'

	input:
        path script
        path countingReads
        path metadata

	output:
        path "*.txt"
        path "*.pdf"


	script:
        """
        Rscript ${script} ${countingReads} ${metadata}
        """
    
}