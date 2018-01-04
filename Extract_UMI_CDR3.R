# The function extracts the UMI and the CDR3 sequence from the sequences and adds these to the result table
# The definition of CDR3 here is the 26 nucletides upstream of the JH region (J1 to J4)

# Input 
# i) file_seq - name of the file with sequences
#ii) dir_seq - location of the file_seq (path)
#iii) file_blast - name of the file with table
#iv) dir_blast - location of the file_blast (path)

#Output
#file with same name as "file_blast" with two additional columns -"UMI" and "CDR3"

library(Biostrings)
Extract_UMI_CDR3  <- function(file_seq, dir_seq, dir_blast,file_blast) {
    
    setwd(dir_seq)  #set directory to dir_seq
    Sequences <- readDNAStringSet(file_seq) #read the sequences from the file_seq
    Seq_new = Sequences[which(attr(attr(Sequences,"ranges"),"width") > 31)]  #exclude sequences less than 30 nts
    NAMES <- names(Seq_new) # extract names of the sequences
    Seq_barcode <- as.character(subseq(Seq_new, start = 23, end = 30)) # extract UMI(23-30) which is after the primer sequence (1-22)
    Barcode <- cbind(NAMES, Seq_barcode) # combine the name of sequence with UMI
    rownames(Barcode) <- (1:dim(Barcode)[1])
    colnames(Barcode) <- c("Query","UMI") # names rows and columns
    Barcode <- as.data.frame(Barcode)   #convert to dataframe
    
    setwd(dir_blast) #set directory to dir_blast
    Data <- read.delim( file_blast,header=T,stringsAsFactors=F)  #read the Igblast result table
    merge_UMI <- merge(Data,Barcode,by="Query")  # merge the Barcode table to the result table
  
    #extract CDR3 as 26 nucleotide region just upstream of JH region
    Add_cdr3 <- function(location,Sequence)
    {
        if(end(location) > 49)
        {
            cdr3_1 <- as.character(subseq(Sequence, start = max(start(location))-26, end = max(start(location))-1)[[1]])
        } else { cdr3_1 <- "NA"}
    
        return(cdr3_1)
    }
    Names_cdr3 <- data.frame(Query=character(), CDR3=character())
    
    #IGHJ1
    n <- vcountPattern("tggggcgcagggaccacggtcacc", Sequences, max.mismatch = 3, fixed = FALSE)
    if(length(which(n==1)) !=0)
    {
        Sequences_J1 <- Sequences[which(n==1)]
        Sequences <- Sequences[which(n==0)]
        Names_j <- names(Sequences_J1)
        m <- vmatchPattern("tggggcgcagggaccacggtcacc", Sequences_J1, max.mismatch = 3, fixed = FALSE)
        cdr3_j <- sapply(1:length(Sequences_J1),function(x){Add_cdr3(m[[x]],Sequences_J1[x])})
        Names_j <- cbind(Names_j,cdr3_j)
        colnames(Names_j) <- c("Query","CDR3")
        Names_j <- as.data.frame(Names_j)
        Names_cdr3 <- rbind(Names_cdr3,Names_j)
        rm(Sequences_J1,Names_j,cdr3_j,m,n)
    }
    #IGHJ2
    n <- vcountPattern("tggggccaaggcaccactctcaca", Sequences, max.mismatch = 3, fixed = FALSE)
    if(length(which(n==1)) !=0)
    {
        Sequences_J1 <- Sequences[which(n==1)]
        Sequences <- Sequences[which(n==0)]
        Names_j <- names(Sequences_J1)
        m <- vmatchPattern("tggggccaaggcaccactctcaca", Sequences_J1, max.mismatch = 3, fixed = FALSE)
        cdr3_j <- sapply(1:length(Sequences_J1),function(x){Add_cdr3(m[[x]],Sequences_J1[x])})
        Names_j <- cbind(Names_j,cdr3_j)
        colnames(Names_j) <- c("Query","CDR3")
        Names_j <- as.data.frame(Names_j)
        Names_cdr3 <- rbind(Names_cdr3,Names_j)
        rm(Sequences_J1,Names_j,cdr3_j,m,n)
    }
    #IGHJ3
    n <- vcountPattern("tggggccaagggactctggtcact", Sequences, max.mismatch = 3, fixed = FALSE)
    if(length(which(n==1)) !=0)
    {
        Sequences_J1 <- Sequences[which(n==1)]
        Sequences <- Sequences[which(n==0)]
        Names_j <- names(Sequences_J1)
        m <- vmatchPattern("tggggccaagggactctggtcact", Sequences_J1, max.mismatch = 3, fixed = FALSE)
        cdr3_j <- sapply(1:length(Sequences_J1),function(x){Add_cdr3(m[[x]],Sequences_J1[x])})
        Names_j <- cbind(Names_j,cdr3_j)
        colnames(Names_j) <- c("Query","CDR3")
        Names_j <- as.data.frame(Names_j)
        Names_cdr3 <- rbind(Names_cdr3,Names_j)
        rm(Sequences_J1,Names_j,cdr3_j,m,n)
    }
    #IGHJ4
    n <- vcountPattern("tggggtcaaggaacctcagtcacc", Sequences, max.mismatch = 3, fixed = FALSE)
    if(length(which(n==1)) !=0)
    {
        Sequences_J1 <- Sequences[which(n==1)]
        Sequences <- Sequences[which(n==0)]
        Names_j <- names(Sequences_J1)
        m <- vmatchPattern("tggggtcaaggaacctcagtcacc", Sequences_J1, max.mismatch = 3, fixed = FALSE)
        cdr3_j <- sapply(1:length(Sequences_J1),function(x){Add_cdr3(m[[x]],Sequences_J1[x])})
        Names_j <- cbind(Names_j,cdr3_j)
        colnames(Names_j) <- c("Query","CDR3")
        Names_j <- as.data.frame(Names_j)
        Names_cdr3 <- rbind(Names_cdr3,Names_j)
        rm(Sequences_J1,Names_j,cdr3_j,m,n)
    }
    #find all possible sequences
    n <- vcountPattern("tggggcgcagggaccacggtcacc", Sequences, max.mismatch = 7, fixed = FALSE)
    if(length(which(n==1)) !=0)
    {
        Sequences_J1 <- Sequences[which(n==1)]
        Names_j <- names(Sequences_J1)
        m <- vmatchPattern("tggggcgcagggaccacggtcacc", Sequences_J1, max.mismatch = 7, fixed = FALSE)
        cdr3_j <- sapply(1:length(Sequences_J1),function(x){Add_cdr3(m[[x]],Sequences_J1[x])})
        Names_j <- cbind(Names_j,cdr3_j)
        colnames(Names_j) <- c("Query","CDR3")
        Names_j <- as.data.frame(Names_j)
        Names_cdr3 <- rbind(Names_cdr3,Names_j)
        rm(Sequences_J1,Names_j,cdr3_j,m,n)
    }
    #merge the merge_UMI table with the CDR3 table
    merge_CDR3 <- merge(merge_UMI,Names_cdr3,by="Query")
    # Keep only the first match of VH, DH and JH
    merge_CDR3$VH_top_match <- sub('\\*.*','',sub( '\\,.*','',merge_CDR3$VH_top_match))
    merge_CDR3$DH_top_match <- sub('\\*.*','',sub( '\\,.*','',merge_CDR3$DH_top_match))
    merge_CDR3$JH_top_match <- sub('\\*.*','',sub( '\\,.*','',merge_CDR3$JH_top_match))
    #write the table
    write.table(merge_CDR3,file_blast, sep = "\t")
    rm(Data, merge_UMI,merge_CDR3,Barcode,Names_cdr3)
}




