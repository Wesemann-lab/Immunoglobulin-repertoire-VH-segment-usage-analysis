#Two functions to filter the Sequences such that only those are slected which have less than 2 pecent "N"s in them.

#Input file is the fasta output file of Seq_File_filter.R function -  make sure to set the path such that 
#the input file is available to the console or provide the complete location of the file.
library(Biostrings)
library(seqinr)
Percent_N_apply_file <-function(file) {
  # create a new folder to save the output if it does not already exist
  if(!file.exists("2_percent"))
  { dir.create("2_percent")}
  #read the Sequence file into R
  Sequences <- readDNAStringSet(file)
  if(length(Sequences) != 0){
  # create an output file
  file_out <- paste("2_percent/",file, sep = "")
  # for each sequence call the function Percent_N_apply_seq
  x = length(Sequences)
  sapply(1:x,function(x){Percent_N_apply_seq(Sequences[x],file_out)})
  rm(Sequences,file_out,file,x)
  }
}


#input file for Percent_N_apply_seq -  this function is called internally by Percent_N_apply_file.
#The sequences are passed to this function one at a time along with the ooutput file name.

Percent_N_apply_seq <- function(Sequence,file_out) {
  #split each sequnce into a vector
  seq <- strsplit(as.character(Sequence), "")[[1]]

  # count the number of "N"s in the sequence and check whether the number is greater than 2 percent of the sequence length : if true than write the sequence to a file
  if (length(which(seq == "N"))< 0.02*length(seq)) {
    write.fasta(sequences = Sequence, names = names(Sequence), nbchar = 80, file.out = file_out, open = "a")
  }
  rm(Sequence,file_out,seq)
}