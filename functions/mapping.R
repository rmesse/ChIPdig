fastq_list_import <- function(fastq_list_path, second_choice) {
  fastq_list <- read.delim(fastq_list_path, sep = "\t", header = F)
  
  if (second_choice==1) {
    if (ncol(fastq_list)!=1) {
      return("ncol_not_1")
    } 
    fastq_list <- within(fastq_list, SampleName <- seq(1,nrow(fastq_list)))
    fastq_list <- within(fastq_list, SampleName <- paste("Sample_",SampleName,sep = ""))
    colnames(fastq_list) <- c("FileName", "SampleName")
    fastq_list$FileName <- as.character(fastq_list$FileName)
  } else {
    
    fastq_list <- within(fastq_list, SampleName <- seq(1,nrow(fastq_list)))
    fastq_list <- within(fastq_list, SampleName <- paste("Sample_",SampleName,sep = ""))
    colnames(fastq_list) <- c("FileName1","FileName2" ,"SampleName")
    fastq_list$FileName1 <- as.character(fastq_list$FileName1)
    fastq_list$FileName2 <- as.character(fastq_list$FileName2)
  }
  return(fastq_list)
}






preprocessing <- function(fastq_list,second_choice,truncateStartBases,truncateEndBases,Lpattern,Rpattern,minLength,nBases) {
  if (second_choice==1) {
    td <- tempdir()
    infiles <- fastq_list$FileName
    outfiles <- paste("filtered_",infiles,sep = "")
    res <- preprocessReads(filename = infiles,
                           outputFilename = outfiles,
                           truncateStartBases = truncateStartBases,
                           truncateEndBases = truncateEndBases,
                           Lpattern = Lpattern,
                           Rpattern = Rpattern,
                           minLength = minLength,
                           nBases = nBases)
  }
  if (second_choice==2) {
    infiles1 <- fastq_list$FileName1
    infiles2 <- fastq_list$FileName2
    outfiles1 <- paste("filtered_",infiles1,sep = "")
    outfiles2 <- paste("filtered_",infiles2,sep = "")
    res <- preprocessReads(filename = infiles1,
                           filenameMate=infiles2,
                           outputFilename = outfiles1,
                           outputFilenameMate=outfiles2,
                           truncateStartBases = truncateStartBases,
                           truncateEndBases = truncateEndBases,
                           Lpattern = Lpattern,
                           Rpattern = Rpattern,
                           minLength = minLength,
                           nBases = nBases)
  }
  return(res)
}

