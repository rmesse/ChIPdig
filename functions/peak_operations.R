# for user-supplied peaks
ext_peaks_import <- function(ext_peaks_path) {
  peak_file_list <- read.delim(ext_peaks_path, sep = "\t", header = F)
  colnames(peak_file_list) <- c("Sample ID", "Peak file","Condition")
  return(peak_file_list)
}
bed_import <- function(bed_path) {
  bed <- read_bed(bed_path, n_fields = 3)
  bed <- bed[!grepl("track", bed$chrom),]
  return(bed)
}
replicated_peaks_batch <- function(peak_table) { 
  ifelse(!dir.exists(file.path("./", "./replicated_peaks/")), dir.create(file.path("./", "./replicated_peaks/")), FALSE)
  conditions <- unique(as.character(peak_table[,3]))
  for (j in 1:length(conditions)) {
    peak_table_subset <- subset(peak_table, peak_table[,3]==as.character(conditions[j]))
    if (nrow(peak_table_subset)==1) { 
      repl_peaks <- bed_import(as.character(peak_table_subset[1,2]))
      setwd("./replicated_peaks/")
      write.table(repl_peaks, file=paste(conditions[j],"_replicated_peaks.bed",sep = ""), quote=F, sep="\t", row.names=F, col.names=F)
      setwd("..")
    } else {
      for (i in 1:nrow(peak_table_subset)) {
        if (i==1) {
          bed1 <- bed_import(as.character(peak_table_subset[i,2]))
          bed2 <- bed_import(as.character(peak_table_subset[i+1,2]))
          repl_peaks <- bed_intersect(bed1,bed2)
          repl_peaks <- dplyr::mutate(repl_peaks, start = pmax(start.x, start.y),end = pmin(end.x, end.y))
          repl_peaks <- repl_peaks[,c(1,which( colnames(repl_peaks)=="start" ),which( colnames(repl_peaks)=="end"))]
        } else {
          if (i==nrow(peak_table_subset)) {
            setwd("./replicated_peaks/")
            write.table(repl_peaks, file=paste(conditions[j],"_replicated_peaks.bed",sep = ""), quote=F, sep="\t", row.names=F, col.names=F)
            setwd("..")
          } else {
            bed2 <- bed_import(as.character(peak_table_subset[i+1,2]))
            repl_peaks <- bed_intersect(repl_peaks,bed2)
            repl_peaks <- dplyr::mutate(repl_peaks, start = pmax(start.x, start.y),end = pmin(end.x, end.y))
            repl_peaks <- repl_peaks[,c(1,which( colnames(repl_peaks)=="start" ),which( colnames(repl_peaks)=="end"))]
          }
        }
      }
    }
  }
}
consensus_peak_set_function <- function(file_list, bin_size) {
  if(length(file_list)==1) {
    consensus_peaks <- bed_import(file_list[1])
    setwd("..")
    write.table(consensus_peaks,"consensus_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)
    setwd("./replicated_peaks/")
  } else {
    for (i in 1:length(file_list)) {
      if (i==1) {
        bed1 <- bed_import(file_list[i])
        bed2 <- bed_import(file_list[i+1])
        consensus_peaks_1 <- bed_intersect(bed1,bed2)
        consensus_peaks_1 <- dplyr::mutate(consensus_peaks_1, start = pmax(start.x, start.y),end = pmin(end.x, end.y))
        consensus_peaks_1 <- consensus_peaks_1[,c(1,which( colnames(consensus_peaks_1)=="start" ),which( colnames(consensus_peaks_1)=="end"))]
        consensus_peaks_2 <- bed_subtract(bed1,bed2)
        consensus_peaks_3 <- bed_subtract(bed2,bed1)
        consensus_peaks <- rbind(consensus_peaks_1,consensus_peaks_2,consensus_peaks_3)
      } else {
        if (i==length(file_list)) {
          setwd("..")
          consensus_peaks <- within(consensus_peaks, width <- consensus_peaks[,3]-consensus_peaks[,2])
          consensus_peaks <- consensus_peaks[consensus_peaks[,4]>bin_size,]
          consensus_peaks <- consensus_peaks[,1:3]
          write.table(consensus_peaks,"consensus_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)
          setwd("./replicated_peaks/")
        } else {
          bed2 <- bed_import(file_list[i+1])
          consensus_peaks_1 <- bed_intersect(consensus_peaks,bed2)
          consensus_peaks_1 <- dplyr::mutate(consensus_peaks_1, start = pmax(start.x, start.y),end = pmin(end.x, end.y))
          consensus_peaks_1 <- consensus_peaks_1[,c(1,which( colnames(consensus_peaks_1)=="start" ),which( colnames(consensus_peaks_1)=="end"))]
          consensus_peaks_2 <- bed_subtract(consensus_peaks,bed2)
          consensus_peaks_3 <- bed_subtract(bed2,consensus_peaks)
          consensus_peaks <- rbind(consensus_peaks_1,consensus_peaks_2,consensus_peaks_3)
          consensus_peaks <- within(consensus_peaks, width <- consensus_peaks[,3]-consensus_peaks[,2])
          consensus_peaks <- consensus_peaks[consensus_peaks[,4]>bin_size,]
          consensus_peaks <- consensus_peaks[,1:3]
        }
      }
    } 
  }
}

# for app-generated peaks
add_track_using_sample_list_colors <- function(sample_list) {
  setwd("./peaks/")
  for (i in 1:nrow(sample_list)) {
    bed <- read_bed(paste(as.character(sample_list[i,1]),"_",as.character(sample_list[i,4]),"_peaks.bed",sep = ""), n_fields = 5)
    bed <- bed[!grepl("track", bed$chrom),]
    bed[,4] <- NA
    track_definition <- bed[1,]
    bed <- rbind(track_definition,bed)
    bed[1,] <- rep(NA,ncol(bed))
    header <- "track name=name_replace description=name_replace color=color_main"
    bed[,1] <- sapply(bed[,1], as.character)
    bed[1,1] <- header
    bed <- as.data.frame(sapply(bed,gsub,pattern="name_replace",replacement=shQuote(as.character(sample_list[i,1]))))
    color_main <- as.vector(col2rgb(sample_list[i,5])[,1])
    color_main <- paste(as.character(color_main), collapse=",")
    bed <- as.data.frame(sapply(bed,gsub,pattern="color_main",replacement= color_main))
    write.table(bed, file=paste(as.character(sample_list[i,1]),"_",as.character(sample_list[i,4]),"_peaks.bed",sep = ""), quote=F, sep="\t",  na = "",row.names=F, col.names=F)
  }
  setwd("..")
}
replicated_peaks_batch_app_peaks <- function(sample_list) { 
  setwd("./peaks/")
  ifelse(!dir.exists(file.path("./", "./replicated_peaks/")), dir.create(file.path("./", "./replicated_peaks/")), FALSE)
  conditions <- unique(as.character(sample_list[,4]))
  for (j in 1:length(conditions)) {
    sample_list_subset <- subset(sample_list, sample_list[,4]==as.character(conditions[j]))
    if (nrow(sample_list_subset)==1) { 
      bed <- read_bed(paste(as.character(sample_list_subset[1,1]),"_",as.character(sample_list_subset[1,4]),"_peaks.bed",sep = ""), n_fields = 5)
      setwd("./replicated_peaks/")
      write.table(repl_peaks, file=paste(conditions[j],"_replicated_peaks.bed",sep = ""), na = "", quote=F, sep="\t", row.names=F, col.names=F)
      setwd("..")
    } else {
      for (i in 1:nrow(sample_list_subset)) {
        if (i==1) {
          bed1 <- read_bed(paste(as.character(sample_list_subset[i,1]),"_",as.character(sample_list_subset[i,4]),"_peaks.bed",sep = ""), n_fields = 5)
          bed2 <- read_bed(paste(as.character(sample_list_subset[i+1,1]),"_",as.character(sample_list_subset[i+1,4]),"_peaks.bed",sep = ""), n_fields = 5)
          bed1 <- bed1[!grepl("track", bed1$chrom),]
          bed2 <- bed2[!grepl("track", bed2$chrom),]
          repl_peaks <- bed_intersect(bed1,bed2)
          temp <- repl_peaks
          repl_peaks <- dplyr::mutate(repl_peaks, start = pmax(start.x, start.y),end = pmin(end.x, end.y))
          repl_peaks <- repl_peaks[,c(1,which( colnames(repl_peaks)=="start" ),which( colnames(repl_peaks)=="end"))]
        } else {
          if (i==nrow(sample_list_subset)) {
            setwd("./replicated_peaks/")
            track_definition <- repl_peaks[1,]
            repl_peaks <- rbind(track_definition,repl_peaks)
            repl_peaks[1,] <- rep(NA,ncol(repl_peaks))
            header <- "track name=name_replace description=name_replace color=color_main"
            repl_peaks[,1] <- sapply(repl_peaks[,1], as.character)
            repl_peaks[1,1] <- header
            repl_peaks <- as.data.frame(sapply(repl_peaks,gsub,pattern="name_replace",replacement=shQuote(as.character(sample_list_subset[i,4]))))
            color_main <- as.vector(col2rgb(sample_list_subset[i,5])[,1])
            color_main <- paste(as.character(color_main), collapse=",")
            repl_peaks <- as.data.frame(sapply(repl_peaks,gsub,pattern="color_main",replacement= color_main))
            write.table(repl_peaks, file=paste(conditions[j],"_replicated_peaks.bed",sep = ""), na = "",quote=F, sep="\t", row.names=F, col.names=F)
            setwd("..")
          } else {
            bed2 <- read_bed(paste(as.character(sample_list_subset[i+1,1]),"_",as.character(sample_list_subset[i+1,4]),"_peaks.bed",sep = ""), n_fields = 5)
            repl_peaks <- bed_intersect(repl_peaks,bed2)
            repl_peaks <- dplyr::mutate(repl_peaks, start = pmax(start.x, start.y),end = pmin(end.x, end.y))
            repl_peaks <- repl_peaks[,c(1,which( colnames(repl_peaks)=="start" ),which( colnames(repl_peaks)=="end"))]
          }
        }
      }
    }
  }
  setwd("..")
}


