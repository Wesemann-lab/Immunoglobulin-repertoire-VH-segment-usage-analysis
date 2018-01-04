# This function calculated the clonal expansion for each CDR3
# Two sequences are clones if they have the same VH, same JH and same CDR3 but a different UMI
# Input for this function is output file from Remove_PCR_repeast.R function
# Output of this function is another tables with the "copies" column giving the number of clones for each sequence.
library(Biostrings)
library(msa)
Calclulate_clonal_expansion <- function(file) {
    # create a folder if it does not already exist
    if(!file.exists("CE"))
    { dir.create("CE")}
    # read the table
    Data <- read.delim(file,header=T,stringsAsFactors=F)
    Data$copies <- 1 #This will give the number of clones for each sequence
    Filter <- Data[0,] #Final table of results
    uniq_VJ <- unique(Data[,c("VH_top_match","JH_top_match")]) find unique VH and JH combinations in the data
    uniq_VJ <- uniq_VJ[!is.na(uniq_VJ$VH_top_match),]
    # for each unique pair
    for(i in 1:length(uniq_VJ[,1]))
    {
        genes <- uniq_VJ[i,]
        reads = Data
        same_VJ <- reads[ reads$VH_top_match == genes$VH_top_match &  reads$JH_top_match == genes$JH_top_match & !is.na(reads$VH_top_match), ]
        uniq_CDR3 <- unique(same_VJ[,c("CDR3")])  #check number of unique CDR3s
        #if the CDR3s are exactly the same, count the number of duplicates and remove duplicates
        cdr3 <- same_VJ[0,]
        for (l in 1: length(uniq_CDR3))
        {
            same_CDR3 <- same_VJ[same_VJ$CDR3 == uniq_CDR3[l],]
            same_CDR3[1,"copies"] = sum(same_CDR3$copies)
            cdr3 <- rbind(cdr3, same_CDR3[1,])
            rm(same_CDR3 )
        }
        same_VJ<- cdr3
        rm(cdr3,uniq_CDR3,l)
        #Convert the CDR3s into string structure for comparison
        AsSeq_CDR3 <- DNAStringSet(as.character(same_VJ[,"CDR3"]))
        while(length(AsSeq_CDR3) != 0)
        {
            if  (length(AsSeq_CDR3) == 1)  #if there is just one sequence - save it
            {
                Filter <- rbind(Filter,same_VJ[1,])
                m = 1
                same_VJ <- same_VJ[!(m),]
                AsSeq_CDR3<-AsSeq_CDR3[!(m)]
                rm(m)
            }
            else
            {   #if there are more than one sequences; check if how many have a miss match of one or less with the first sequence
                m <- vcountPattern(AsSeq_CDR3[[1]], AsSeq_CDR3 ,max.mismatch = 1,fixed = FALSE)
                if(length(which(m != 0)) == 1)
                { #if no other sequence matches save the first sequence
                    Filter <- rbind(Filter,same_VJ[which(m != 0)[1],])
                    same_VJ <- same_VJ[(which(m == 0)),]
                    AsSeq_CDR3<-AsSeq_CDR3[(which(m == 0))]
                    rm(m)
                }
                else
                { #if there are other macthing sequences a consensus is created
                    myaln <- msa(AsSeq_CDR3[which(m != 0)])
                    conMat_seq <- consensusMatrix(myaln)
                    conMat_seq <- conMat_seq[rowSums(conMat_seq)>0,]
                    consensus = character()
                    for (k in 1:dim(conMat_seq)[2])
                    {
                        if(colSums(conMat_seq)[k] > 1)
                        {
                            x <- conMat_seq[,k]
                            rownames(x) <- rownames(conMat_seq[,k])
                            x <- x[order(-x)]
                            if(ROWNAMES(x)[1] != "N")
                            {
                                consensus = rbind(consensus, ROWNAMES(x)[1])
                            }
                            else if(x[2] !=0 && ROWNAMES(x)[2] != "N" )
                            {
                                consensus = rbind(consensus,ROWNAMES(x)[2])
                            } else if( !(is.na(x[3] !=0)) && x[3] !=0 && ROWNAMES(x)[3] != "N" )
                            {
                                consensus = rbind(consensus,ROWNAMES(x)[3])
                            } else
                            {
                                consensus = rbind(consensus,"N")
                            }
                        }
                    }
                    
                    Seq <- paste(consensus,collapse="")
                    rm(myaln,conMat_seq,k,x)
                    #the consensus is added to the table with the number of matching sequences recorded to keep the count
                    noN <-same_VJ[which(m != 0)[1],]
                    noN$CDR3 <- Seq
                    noN$copies = sum(same_VJ$copies[which(m != 0)])
                    Filter <- rbind(Filter,noN)
                    same_VJ <- same_VJ[(which(m == 0)),]
                    AsSeq_CDR3<-AsSeq_CDR3[(which(m == 0))]
                    rm(m,noN,consensus, Seq)
                    
                }
            }
        } 
        rm(same_VJ,AsSeq_CDR3,same_VJ)
    }
    file_out <- paste("CE/",  file,sep = "")
    # write to output file
    write.table(Filter,file_out , sep="\t")        
    rm(Data_matrix,Filter,reads,genes,uniq_VJ,file_out,i,UMI_all_seq,UMI_unique_seq) 
  
}

