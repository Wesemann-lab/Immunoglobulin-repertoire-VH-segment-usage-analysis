# The function removes the pcr repeats
# Definition of PCR repeats is - same VH , same JH, same UMI and same CDR3
# Input for this function is the result table from Extract_UMI_CDR3
# Output from this function is a table with all PCR repeats removes. The number of PCR repeats in menstioned in the column "copies" of the table.
library(Biostrings)
library(msa)
Remove_PCR_repeats <- function(file) {
    # create a folder if it does not already exist
    if(!file.exists("no_pcr"))
    { dir.create("no_pcr")}
    # read the table
    Data <- read.delim(file,header=T,stringsAsFactors=F)
    #check and remove incomplete sequences
    Data <- Data[Data$VH_top_match != "N/A" & Data$JH_top_match != "N/A"  & Data$J.region_start != "N/A" &  !is.na(Data$VH_top_match), ]
    # add a counter for number of pcr repeats for each sequence
    Data$copies <- 1
    # matrix to same the output
    Filter <- Data[0,]
    # find all unique combintaions of VH and JH
    uniq_VJ <- unique(Data[,c("VH_top_match","JH_top_match")])
    uniq_VJ <- uniq_VJ[!is.na(uniq_VJ$VH_top_match),]
    #print(length(uniq_VJ[,1]))
    #for each unique VJ extract all sequences; if they have same UMI and CDR3 consider them as PCR repeats
    for(i in 1:length(uniq_VJ[,1]))
    {
        genes <- uniq_VJ[i,]
        reads = Data
        filter_1 <- Data[0,]
        # all sequences with same V and J
        same_VJ <- reads[ reads$VH_top_match == genes$VH_top_match &  reads$JH_top_match == genes$JH_top_match & !is.na(reads$VH_top_match), ]
        #among the sequences with same VJ, find unique barcodes
        uniq_Barcode <- unique(same_VJ[,c("UMI")])
        # foreach barcode
        for(j in 1: length(uniq_Barcode))
        {
            # extract all seq with saem barcode
            same_barcode <- same_VJ[ same_VJ$UMI == uniq_Barcode[j],]
            # same CDR3 amonth same barcode
            uniq_CDR3 <- unique(same_barcode[,c("CDR3")])
            # remove same umi same cdr3 - exact same
            cdr3 <- same_barcode[0,]
            for (l in 1: length(uniq_CDR3))
            {
                same_CDR3 <- same_barcode[same_barcode$CDR3 == uniq_CDR3[l],]
                same_CDR3[1,"copies"] = length(same_CDR3[,1])
                cdr3 <- rbind(cdr3, same_CDR3[1,])
                rm(same_CDR3 )
            }
            same_barcode <- cdr3
            rm(cdr3,uniq_CDR3,l)
            #match shorter to longer, replace N and remove shorter
            if (length(same_barcode[,1]) < 2)
            {
                filter_1 <- rbind(filter_1,same_barcode)
            }  else
            {
                #if the CDR3s are different; if they have a mismatch of two, they are still considered pcr repeats and a consensus is created and used as CDR3
                # chances are that "N" may be removed when creating consensus
                AsSeq_CDR3_barcode <- DNAStringSet(as.character(same_barcode[,"CDR3"]))  #convert CDR3 into DNA strings
                while(length(AsSeq_CDR3_barcode) != 0)
                 {
                  if  (length(AsSeq_CDR3_barcode) == 1) #if there is just one sequence - it is included
                  {
                    filter_1 <- rbind(filter_1,same_barcode[1,])
                    m = 1
                    same_barcode <- same_barcode[!(m),]
                    AsSeq_CDR3_barcode<-AsSeq_CDR3_barcode[!(m)]
                    rm(m)
                  } else 
                  {
                    #find the sequences which have a mismatch of two or less with the first sequence
                    m <- vcountPattern(AsSeq_CDR3_barcode[[1]], AsSeq_CDR3_barcode ,max.mismatch = 2, fixed = FALSE)
                    #if there is no other matching sequences, this sequence in included
                    if(length(which(m != 0)) == 1)
                    {
                        filter_1 <- rbind(filter_1,same_barcode[which(m != 0)[1],])
                        same_barcode <- same_barcode[(which(m == 0)),]
                        AsSeq_CDR3_barcode<-AsSeq_CDR3_barcode[(which(m == 0))]
                        rm(m)
                    } else 
                    {  #if there are other matching sequences a consensus is made using these sequences
                        myaln <- msa(AsSeq_CDR3_barcode[which(m != 0)])
                        conMat_seq <- consensusMatrix(myaln)
                        conMat_seq <- conMat_seq[rowSums(conMat_seq)>0,]
                        consensus = character()
                        for (k in 1:dim(conMat_seq)[2])
                        {
                            x <- conMat_seq[,k]
                            rownames(x) <- rownames(conMat_seq[,k])
                            x <- x[order(-x)]
                            if(ROWNAMES(x)[1] != "-" && ROWNAMES(x)[1] != "N")
                            {
                                consensus = rbind(consensus, ROWNAMES(x)[1])
                            } else if(x[2] !=0 && ROWNAMES(x)[2] != "-" && ROWNAMES(x)[2] != "N" )
                            {
                                consensus = rbind(consensus,ROWNAMES(x)[2])
                            } else if( !(is.na(x[3] !=0)) && x[3] !=0 && ROWNAMES(x)[3] != "-" && ROWNAMES(x)[3] != "N" )
                            {
                                consensus = rbind(consensus,ROWNAMES(x)[3])
                            } else 
                            { 
                                    consensus = rbind(consensus,"N")
                            }
                        }
                        
                        Seq <- paste(consensus,collapse="")  #consensus sequence is added
                        rm(myaln,conMat_seq, x, k)
                        noN <-same_barcode[which(m != 0)[1],]
                        noN$CDR3 <- Seq
                        noN$copies <-sum(same_barcode$copies[which(m != 0)])  # no of pcr repeats is the total number of matched sequences
                        rm(Seq)
                        filter_1 <- rbind(filter_1,noN)
                        same_barcode <- same_barcode[(which(m == 0)),]
                        AsSeq_CDR3_barcode<-AsSeq_CDR3_barcode[(which(m == 0))]
                        rm(m,noN,consensus)
                    }
                  } 
                 }
            }
            rm(same_barcode)
        }
        Filter <- rbind(Filter, filter_1)
        rm(filter_1,AsSeq_CDR3_barcode,uniq_Barcode,same_VJ)
    }
  #write to file in a different folder
  file_out <- paste("no_pcr/", file ,sep = "")
  write.table(Filter,file_out , sep="\t")
}
  
  
    