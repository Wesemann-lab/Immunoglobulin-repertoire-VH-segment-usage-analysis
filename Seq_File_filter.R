#Two functions Seq_File_filter and Lower_score_N work together to convert a Fastq file to a Fasta file 
#and replace intermittend low score (phred score <= 20) nucleotides by "N". If four consequtive "N"s are 
#found, the sequnce is truncated at the point.

# input files for Seq_File_Filter -
#       i) file - the merged output file from pear
#       ii) phred_file - A file containing the Phred score values with their respective ascii characters

# initailize the libraries required for this function
library(ShortRead)
library(seqinr)
Seq_File_filter <-function(file,phred_file){
   
    Phred_score <- as.matrix(read.table(phred_file, sep="\t",header=TRUE))
    Phred_score[3,2] = '#'
    Phred_score[7,2] = '\''
    
    # create a new folder to save the output if it does not already exist
    if(!file.exists("Fasta_N"))
    { dir.create("Fasta_N")}

    #Read the Sequence file into R
    Seq <- readFastq(file)
    #create the name and path of the output fasta file
    file_out <- gsub("tq", "ta", file)
    file_out <-paste("Fasta_N/",file_out,sep="")
    
    #For each sequence, call function Lower_score_N
    x = length(Seq)
    if(x>0){
        sapply(1:x,function(x){Lower_score_N(Seq[x],file_out,Phred_score)})
        rm(Seq, file_out,file,phred_file,Phred_score,x)
    }
}



#input file for Lower_score_N -  this function is called internally by Seq_File_filter. 
#The sequences are passed to this function one at a time along with the ooutput file name.

Lower_score_N <- function(Sequence,file_out,Phred_score){
  #split the seq into a vector
  Nucl_seq <- strsplit(as.character(sread(Sequence)[[1]]), "")[[1]]
  #split the score into a vector
  Score_seq <- strsplit(as.character(quality(Sequence)[[1]]), "")[[1]]
  #scan each nucleotide for score and replace by "N" if score is <= 20; 
  #break if there are more than 4 "N"s
  low = 0
  for (j in 1: length(Score_seq)) {
    h = match(Score_seq[j], Phred_score[ ,2])
    if(low >= 4){
      j = j - 5
      break
    } else 
      if (h < 21) {
        Nucl_seq[j] = "N"
        low = low + 1
      } else if (h >=  21) {low = 0 }
  }
  #print(j)
  #convert the sequence vector back to DNA string format
  Seq_filtered <- DNAStringSet(paste(Nucl_seq[1:j],sep = "", collapse =""))  
  #write the sequence to the output fasta file
  write.fasta(sequences = (paste(Nucl_seq[1:j],sep = "", collapse ="")) , names = Sequence@id[[1]], nbchar = 80, file.out = file_out, open = "a",as.string = TRUE)
  rm(Nucl_seq,Score_seq,low,j,h,Sequence,file_out)
}

