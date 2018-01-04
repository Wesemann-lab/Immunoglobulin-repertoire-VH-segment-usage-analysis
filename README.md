# Immunoglobulin repertoire VH segment usage analysis

This pipeline merges, filters, annotates and calculates the VH gene segment usage in Illumina MiSeq generated paired end B cell immunoglobulin heavy chain sequences. It includes removal of PCR repeats using Unique Molecular Identifiers (UMI).

This is the pipeline used for analysis of data in the publication "Microbial symbionts regulate the primary Ig repertoire".

It is a very simple workflow which works mostly in R. Pear and IgBlast were used in intermediate steps. The functions used and the steps followed are described below.

Requirements:

1) R and R studio
2) Pear (Paired-End reAd mergeR) (https://sco.h-its.org/exelixis/web/software/pear/doc.html)
4) IgBlast (standalone) (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/)


Make sure that all the bio-packages listed at the beginning of each R script are installed in R. 

Work flow:

 1) Merge the paired end reads by running Pear with the following parameters
		i)  minimum overlap for merge (-v) = 20
		ii) quality score threshold for trimming the low quality part of a read (Phred score) (-q) = 20
		iii)minimum possible length of the assembled sequences = 80

 2) The following function in R replaces intermittent low quality reads(Phred score <= 20) by “N” and converts Fastq file to a Fasta file. It also trims the sequences if it finds four consecutive low quality reads. 
	Seq_File_filter.R - The input for this file is the output file from pear (merged reads) and a file containing the Phred scores along with their ascii characters(available with the files) to decode the characters into scores. This function also calls another function “Lower_score_N.R”. Both the functions should be sources into the workspace before running the first function.

 3) In-order to exclude very poor quality sequences, a filter is applied to allow sequences with a maximum of 2 percent ’N’s in them to be included in the next step. All sequences with greater than 2 % “N”s were excluded. 
	Percent_N_apply_file - The input for this function is the output from Seq_File_filter.R. This function calls another function Percent_N_apply_Seq.R present in same file. Source the file to load both function into workspace.

 4) The sequences are run through IgBlast to be aligned with the germline sequences and annotated. This requires the standalone igBlast to be installed and the databases from IMGT to be downloaded and converted to blast data format.

 5) The output results from igBlast are parsed and converted into table format using the R program ParseIgblast.R. The input for this file is the output file from the standalone IgBlast program.

 6) The next step involves extracting the Unique Molecular Identifier(UMI) and the CDR3 region. In this case the CDR3 is described as the 26 nucleotide region just upstream of the J Heavy region. The R function Extract_UMI_CDR3.R extracts the two sequences and adds them to the result table created by from ParseIgBlast.R. 
    The inputs for this function are:
	i) file_seq - name of the file with sequences
	ii) dir_seq - location of the file_seq (path)
	iii) file_blast - name of the file with table
	iv) dir_blast - location of the file_blast (path)
	
   The output of this function is a table with Igblast results along with “UMI” and “CDR3” columns added to it.

 7) The next function (Remove_PCR_repeats.R) in the pipeline identifies all the PCR repeats and removes them. PCR repeats are defined as sequences with same VH , same JH, same UMI and same CDR3. Input for this function is the result table from Extract_UMI_CDR3.R. Output from this function is a table with all PCR repeats removes. The number of PCR repeats in mentioned in the column “copies” of the table.

 8) The function Calculate_clonal_expansion.R calculates the number of clones for each sequence. Two sequences are clones of each other if they have the same VH, same JH and same CDR3 but a different UMI. Input for this function is the result table from Remove_PCR_repeats.R. Output from this function is a table with the count of number of clones for each sequence mentioned in the “copies” column of the table. 

 These counts are then used to calculate the VH segment usage frequency.

