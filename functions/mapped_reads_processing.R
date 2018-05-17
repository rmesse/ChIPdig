# 1 - sample loading
sample_list_import <- function(sample_list_path) {
  sample_list <- read.delim(sample_list_path, sep = "\t", header = F)
  if (ncol(sample_list)!=5) {
    return("ncol_not_5")
  } 
  if (sum(sample_list=="")>0) {
    return("missing_values")
  }
  sample_list[,2] <- as.character(sample_list[,2])
  sample_list[,3] <- as.character(sample_list[,3])
  bam_files <- c(sample_list[,2],sample_list[,3])  
  for (i in 1:length(bam_files)) {
    temp <- length(list.files(pattern = bam_files[i] ))
    if (i==1) {
      proceed <- temp
    } else {
      proceed <- temp * proceed
    }
  }
  
  if (proceed==0) {
    return("bam_missing")
  }
  colnames(sample_list) <- c("Sample ID", "Mapped reads file, treatment","Mapped reads file, control","Condition","Color for exported peak and coverage files")
  return(sample_list)
}


# 2a - calculation of number of reads and fragment lengths  
fragment_size_estimation <- function(BAM_file_path) {
  dedup.on <- readParam()
  dedup.on <- reform(dedup.on, dedup=TRUE)
  x <- correlateReads(BAM_file_path, param=dedup.on)
  x <- maximizeCcf(x)
  return(x)
}
bam.files_list <- function(sample_list,list,i,j) {
  if (i==1) {
    list <- rep(NA,nrow(sample_list))
    list[1] <- as.character(sample_list[1,2])
    list[2] <- as.character(sample_list[1,3])
    return(bam.files_list(sample_list,list,2,1))
  }
  if (i>nrow(sample_list)) {
    return(list)
  } else {
    if (j==1) {
      list[i*2-1] <- as.character(sample_list[i,2])
      return(bam.files_list(sample_list,list,i,2))
    } else {
      list[i*2] <- as.character(sample_list[i,3])
      return(bam.files_list(sample_list,list,i+1,1))
    }
  }
}
fragment_size_estimation_batch <- function(list) {
  list2 <- rep(0,length(list))
  for (i in 1:length(list)) {
    list2[i] <- fragment_size_estimation(as.character(list[i]))
  }
  list2 <- as.numeric(as.character(list2))
  return(list2)
}
add_vector_elements_to_last_2_columns_of_frame <- function(frame,vector) {
  for (i in 1:length(vector)) {
    if (i %% 2 != 0) {
      frame[(i+1)/2,(ncol(frame)-1)] <- vector[i]
    } else {
      frame[(i+0)/2,ncol(frame)] <- vector[i]
    }
  }
  return(frame)
 }

# 2b - guessing of chromosome sizes from bam files
bam_import_for_sizes_estimation <- function(bam_file_path,single_or_paired,strandMode) {
  BamFile <- BamFile(bam_file_path)
  param <- ScanBamParam(mapqFilter=20, flag=scanBamFlag(isDuplicate=FALSE))
  if(single_or_paired=="single") {ga <- readGAlignments(BamFile, param=param)} else {ga <- readGAlignmentPairs(BamFile, param=param,strandMode=strandMode)}
  gr <- granges(ga)
  return(gr)
}
GR_to_DF_to_spans <- function(gr) {
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr),
                   ends=end(gr),
                   names=c(rep(".", length(gr))),
                   scores=c(rep(".", length(gr))),
                   strands=strand(gr))
  df <- within(df, starts <- 1)
  chroms <- levels(df$seqnames)
  for (i in 1:length(chroms)) {
    temp <- subset(df, seqnames==chroms[i])
    temp <- subset(temp, ends == max(temp$ends,na.rm = T))
    temp <- temp[!duplicated(temp), ]
    assign(paste("subset_pulp_fiction_",i, sep=""),temp)
  }
  df_summary <- do.call(rbind, mget(ls(pattern="subset_pulp_fiction_")))
  return(df_summary)
}
chrom_sizes_fun <- function(sample_list,single_or_paired,strandMode) {
  for (i in 1:nrow(sample_list)) {
    gc()
    temp1 <- bam_import_for_sizes_estimation(as.character(sample_list[i,2]),single_or_paired,strandMode)
    temp1 <- GR_to_DF_to_spans(temp1)
    temp2 <- bam_import_for_sizes_estimation(as.character(sample_list[i,3]),single_or_paired,strandMode)
    temp2 <- GR_to_DF_to_spans(temp2)
    assign(paste("subset_pulp_fiction_batch_treat_",i, sep=""),temp1)
    assign(paste("subset_pulp_fiction_batch_control_",i, sep=""),temp2)
  }
  df_summary_1 <- do.call(rbind, mget(ls(pattern="subset_pulp_fiction_batch_treat_")))
  df_summary_2 <- do.call(rbind, mget(ls(pattern="subset_pulp_fiction_batch_control_")))
  df_summary <- rbind(df_summary_1,df_summary_2)
  chroms <- levels(df_summary$seqnames)
  for (i in 1:length(chroms)) {
    temp <- subset(df_summary, seqnames==chroms[i])
    temp <- subset(temp, ends == max(temp$ends,na.rm = T))
    temp <- temp[!duplicated(temp), ]
    assign(paste("subset_pulp_fiction_chroms_",i, sep=""),temp)
  }
  df_summary <- do.call(rbind, mget(ls(pattern="subset_pulp_fiction_chroms_")))
  df_summary <- within(df_summary, strands <- "+")
  df_summary <- df_summary[!duplicated(df_summary),]
  row.names(df_summary) <- NULL
  return(df_summary)
}

# 2c - preparation for differential enrichment analysis and coverage export
replicates_present_assessment <- function(sample_list) {
  conditions <- as.character(sample_list$Condition)
  duplicates_out <- conditions[!duplicated(conditions)]
  output <- ifelse(length(duplicates_out)!=length(conditions),TRUE,FALSE)
  return(output)
}
counts_to_edger <- function(counts,sample_list) {
  group <- rep(NA,nrow(sample_list)*2)
  for (i in 1:nrow(sample_list)) {
    group[i*2-1] <- paste(sample_list[i,grep("^Condition$", colnames(sample_list))],"_treatment",sep = "")
    group[i*2-0] <- paste(sample_list[i,grep("^Condition$", colnames(sample_list))],"_control",sep = "")
  }
  samples <- rep(NA,nrow(sample_list)*2)
  for (i in 1:nrow(sample_list)) {
    samples[i*2-1] <- paste(sample_list[i,grep("^Sample ID$", colnames(sample_list))],"_treatment",sep = "")
    samples[i*2-0] <- paste(sample_list[i,grep("^Sample ID$", colnames(sample_list))],"_control",sep = "")
  }
  colnames(counts) <- samples
  edger <- DGEList(assay(counts), group=group)
  return(edger)
}
edger_to_TMM <- function(edger) {
  TMM <- calcNormFactors(edger, method="TMM")
}
coverages_function <- function(sample_list,TMM) {
  coverage_matrix <- matrix(0,nrow=nrow(TMM$counts),ncol=nrow(sample_list)*3)
  for (i in 1:nrow(sample_list)) {
    coverage_matrix[,i*3-2] <- cpm(TMM[,i*2-1], normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)
    coverage_matrix[,i*3-1] <- cpm(TMM[,i*2-0], normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)
    coverage_matrix[,i*3-0] <- coverage_matrix[,i*3-2] - coverage_matrix[,i*3-1]
  }
  colnames <- rep(NA,ncol(coverage_matrix))
  for (i in 1:nrow(sample_list)) {
    colnames[i*3-2] <- paste(sample_list[i,1],"_treatment",sep = "")
    colnames[i*3-1] <- paste(sample_list[i,1],"_control",sep = "")
    colnames[i*3-0] <- paste(sample_list[i,1],"_normalized",sep = "")
  }
  colnames(coverage_matrix) <- colnames
  return(coverage_matrix)
}

# 3 - coverage export
export_bedgraphs <- function(counts,edger,sample_list) {
  ifelse(!dir.exists(file.path("./", "./coverage_files/")), dir.create(file.path("./", "./coverage_files/")), FALSE)
  setwd("./coverage_files/")
  df <- data.frame(seqnames=seqnames(rowRanges(counts)),
                   starts=start(rowRanges(counts)),
                   ends=end(rowRanges(counts)))
  df <- as.data.frame(df)
  sample_list_original <- sample_list
  sample_list <- within(sample_list, background_color <- "gray85")
  sample_list <- within(sample_list, norm_color <- sample_list[,ncol(sample_list)-1])
  sample_list <- sample_list[,(ncol(sample_list)-2):ncol(sample_list)] 
  sample_list <- c(t(sample_list))
  color_bck <- as.vector(col2rgb("gray85" )[,1])
  color_bck <- paste(as.character(color_bck), collapse=",")
  for (j in 1:nrow(sample_list_original)) {
    subset_edger <- edger[,(j*2-1):(j*2)]
    subset_TMM <- edger_to_TMM(subset_edger)
    sample_list_j <- sample_list_original[j,]
    subset_coverage <- coverages_function(sample_list_j,subset_TMM)
    subset_coverage <- as.data.frame(subset_coverage)
    for (i in 1:ncol(subset_coverage)) {
      temp <- cbind(df,subset_coverage[,i])
      track_definition <- temp[1,]
      temp[,ncol(temp)] <- round(temp[,ncol(temp)],3)
      temp <- rbind(track_definition,temp)
      temp[1,] <- rep(NA,ncol(temp))
      header <- "track type=bedGraph name=name_replace description=name_replace visibility=full color=color_main altColor=color_bck windowingFunction=mean smoothingWindow=10"
      temp[,1] <- sapply(temp[,1], as.character)
      temp[1,1] <- header
      temp <- as.data.frame(sapply(temp,gsub,pattern="name_replace",replacement=shQuote(colnames(subset_coverage[i]))))
      if (i==2) {
        color_main <- color_bck
      } else {
        color_main <- as.vector(col2rgb(sample_list_original[j,5])[,1])
        color_main <- paste(as.character(color_main), collapse=",")
      }
      temp <- as.data.frame(sapply(temp,gsub,pattern="color_main",replacement= color_main))
      temp <- as.data.frame(sapply(temp,gsub,pattern="color_bck",replacement= color_bck))
      write.table(temp, file=paste(colnames(subset_coverage[i]),".bedgraph",sep = "") , quote=F, sep="\t", na="", row.names=F, col.names=F)
    }
  }
  setwd("..")
}

# 4 - peak calling
prepareChIPseq_defined <- function(reads,fr_leng){
  reads.extended <- trim(resize(reads, width = as.integer(fr_leng)))
  return(reads.extended)
}
bam_import <- function(bam_file_path,single_or_paired,dupl_reads_remove,strandMode,read_extension_option,fr_leng) {
  BamFile <- BamFile(bam_file_path)
  if(dupl_reads_remove==TRUE) {
    param <- ScanBamParam(mapqFilter=20, flag=scanBamFlag(isDuplicate=FALSE))
  } else {
    if(dupl_reads_remove==FALSE) {
      param <- ScanBamParam(mapqFilter=20, flag=scanBamFlag(isDuplicate=NA))
    }
  }
  if(single_or_paired=="single") {ga <- readGAlignments(BamFile, param=param)} else {ga <- readGAlignmentPairs(BamFile, param=param,strandMode=strandMode)}
  gr <- granges(ga)
  if (read_extension_option=="YES") {gr <- prepareChIPseq_defined(gr,fr_leng)}
  return(gr)
}
peaks_batch <- function(sample_list,sample_list_with_fr_sizes,single_or_paired,dupl_reads_remove,strandMode,read_extension_option,PP,bin_size) {
  ifelse(!dir.exists(file.path("./", "./peaks/")), dir.create(file.path("./", "./peaks/")), FALSE)  
  for (i in 1:nrow(sample_list)) {
    gc()
    if (i==1) {
      output_table <- sample_list
      output_table <- within(output_table, Peak_count <- 0)
      }
    temp1 <- bam_import(as.character(sample_list_with_fr_sizes[i,2]),single_or_paired,dupl_reads_remove,strandMode,read_extension_option,as.numeric(as.character(sample_list_with_fr_sizes[i,(ncol(sample_list_with_fr_sizes)-3)])))
    temp2 <- bam_import(as.character(sample_list_with_fr_sizes[i,3]),single_or_paired,dupl_reads_remove,strandMode,read_extension_option,as.numeric(as.character(sample_list_with_fr_sizes[i,(ncol(sample_list_with_fr_sizes)-2)])))
    df1 <- data.frame(seqnames=seqnames(temp1),
                     starts=start(temp1),
                     ends=end(temp1),
                     names=c(rep(".", length(temp1))),
                     scores=c(rep(".", length(temp1))),
                     strands=strand(temp1))
    df1[2:3] <- lapply(df1[2:3], as.integer)
    df2 <- data.frame(seqnames=seqnames(temp2),
                      starts=start(temp2),
                      ends=end(temp2),
                      names=c(rep(".", length(temp2))),
                      scores=c(rep(".", length(temp2))),
                      strands=strand(temp2))
    df2[2:3] <- lapply(df2[2:3], as.integer)
    setwd("./peaks/")
	write.table(df1, file=paste(as.character(sample_list[i,2]),".bed",sep = ""), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(df2, file=paste(as.character(sample_list[i,3]),".bed",sep = ""), quote=F, sep="\t", row.names=F, col.names=F)
    result <- bayespeak(paste(as.character(sample_list[i,2]),".bed",sep = ""),paste(as.character(sample_list[i,3]),".bed",sep = ""),bin.size = bin_size)
    file.remove(paste(as.character(sample_list[i,2]),".bed",sep = ""))
    file.remove(paste(as.character(sample_list[i,3]),".bed",sep = ""))
    result <- summarize.peaks(result, threshold = PP)
    output_table[i,ncol(output_table)] <- nrow(result)
    write.table(result, file=paste(as.character(sample_list[i,1]),"_",as.character(sample_list[i,4]),"_peaks.bed",sep = ""), quote=F, sep="\t", row.names=F, col.names=F)
	setwd("..")
  }  
  return(output_table)
}

# 5 - differential enrichment analysis
counts_input_filter <- function(coverages,FC,counts) {
  filter_logical <- matrix(TRUE,nrow=nrow(coverages),ncol=ncol(coverages)/3)
  for (i in 1:ncol(coverages)) {
    if (i%%3 == 0) {
      filter_logical[,i/3] <- coverages[,i] > log2(FC)
    } else {}
  }
  filter_logical <- rowSums(filter_logical)
  keep <- filter_logical != 0
  counts_backg_removed <- counts[keep,]
  return(counts_backg_removed)
}
edger_subset_keeping_params_all_samples <- function(TMM,counts_subset,sample_list) {
  libsizes <- TMM$samples$lib.size
  normfacs <- TMM$samples$norm.factors
  group <- TMM$samples$group
  samples <- rep(NA,nrow(sample_list)*2)
  for (i in 1:nrow(sample_list)) {
    samples[i*2-1] <- paste(sample_list[i,grep("^Sample ID$", colnames(sample_list))],"_treatment",sep = "")
    samples[i*2-0] <- paste(sample_list[i,grep("^Sample ID$", colnames(sample_list))],"_control",sep = "")
  }
  colnames(counts_subset) <- samples
  edger_subset <- DGEList(assay(counts_subset), group=group, norm.factors=normfacs, lib.size = libsizes)
}
peaks_import <- function(peaks_path) {
  peaks <- read.delim(peaks_path, header = F)
  peaks <- peaks[,1:3]
  colnames(peaks) <- c("chrom","start","end")
  peaks <- peaks[!grepl("track", peaks[,1]),]
  peaks$chrom <- factor(peaks$chrom)  
  peaks <- makeGRangesFromDataFrame(peaks)
  return(peaks)
}
binranges_fun <- function(coverages,sample_list,contrast_test,contrast_ref,counts_subset) {
  indexes_to_keep <- seq(1,ncol(coverages), by = 3)
  coverages <- coverages[,indexes_to_keep]
  colnames(coverages) <- sample_list[,4]
  coverages_ref <- coverages[,which(colnames(coverages)==contrast_ref)]
  coverages_ref <- rowMeans(as.data.frame(coverages_ref), na.rm = T)
  coverages_test <- coverages[,which(colnames(coverages)==contrast_test)]
  coverages_test <- rowMeans(as.data.frame(coverages_test), na.rm = T)
  FC <- cbind(coverages_test,coverages_ref)
  binranges <- data.frame(chrom=seqnames(rowRanges(counts_subset)),
                          start=start(rowRanges(counts_subset))-1,
                          end=end(rowRanges(counts_subset)))
  result <- cbind(binranges,FC)
  return(result)
}
diff_enrich_tab_fun <- function(peaks,bins_with_FC_values,contrast_test,contrast_ref) {
  peaks_DF <- data.frame(chrom=seqnames(peaks),
                         start=start(peaks),
                         end=end(peaks),
                         PValue=elementMetadata(peaks)$PValue,
                         FDR=elementMetadata(peaks)$FDR)
  total_peak_count <- nrow(peaks_DF)
  peaks_DF <- within(peaks_DF,peak_ID <- seq(1,total_peak_count))
  peaks_DF <- bed_intersect(peaks_DF,bins_with_FC_values)
  peaks_DF <- as.data.frame(peaks_DF)
  colnames(peaks_DF) <- c("chrom","start","end","PValue","FDR","peakID","start","end","CPMt","CPMr","overlap")
  peaks_DF <- peaks_DF[,c(1:6,9:11)]
  for (i in 1:total_peak_count) {
    peak_summary <- subset(peaks_DF,peakID==as.integer(i))
    if (nrow(peak_summary)>1) {
      peak_summary_init <- peak_summary
      peak_summary <- within(peak_summary,prod_cov_test <- peak_summary[,7]*peak_summary[,9])
      peak_summary <- within(peak_summary,prod_cov_ref <- peak_summary[,8]*peak_summary[,9])
      peak_summary_final <- within(peak_summary_init,CPMt <- sum(peak_summary$prod_cov_test)/sum(peak_summary$overlap) )
      peak_summary_final <- within(peak_summary_final,CPMr <- sum(peak_summary$prod_cov_ref)/sum(peak_summary$overlap) )
      peak_summary_final <- within(peak_summary_final,overlap<- sum(peak_summary$overlap) )
      peak_summary <- peak_summary_final[1,]
    }
    if (i==1) {
      result <- peak_summary
    } else {
      result <- rbind(result,peak_summary)
    }
  }
  result <- result[,c(1:5,7:8)]
  colnames(result) <- c("Chromosome","Start","End","P value","FDR",paste("Average Log2 CPM per bin,",contrast_test),paste("Average Log2 CPM per bin,",contrast_ref))
  result <- within(result,FC<-result[,6]-result[,7])
  colnames(result)[8] <- "Fold enrichment"
  return(result)
}

