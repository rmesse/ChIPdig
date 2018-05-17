# coverages import and chromosome sizes
coverage_import <- function(coverage_path) {
  coverage <- read_bedgraph(coverage_path)
  coverage <- coverage[!grepl("track", coverage$chrom),]
  return(coverage)
}
coverage_import_batch <- function(sample_list) {
  length_list <- nrow(sample_list)
  cov_list <- list()
  for( i in 1:length_list ){
    cov_list[[ i ]] <- coverage_import(as.character(sample_list[i,2]))
  }
  return(cov_list)
}
chrom_sizes_from_coverage <- function(coverage) {
  chroms <- unique(coverage$chrom)
  size <- rep(1,length(chroms))
  for (i in 1:length(size)) {
    temp <- subset(coverage, chrom==chroms[i])
    size[i] <- max(temp$end,na.rm = T)
  }
  genome <- as.data.frame(cbind(chroms,size)); colnames(genome) <- c("chrom","size")
  genome$chrom <- as.character(genome$chrom)
  genome$size <- as.integer(as.character(genome$size))
  return(genome)
}
chrom_sizes_from_coverage_batch <- function(coverages) {
  for (i in 1:length(coverages)) {
    if (i==1) {
      genome <- chrom_sizes_from_coverage(coverages[[1]])
    } else {
      temp <- chrom_sizes_from_coverage(coverages[[i]])
      genome <- rbind(temp,genome)
    }
  }
  chromosomes <- unique(genome[,1])
  genome_unique <- rep(NA,length(chromosomes))
  for (i in 1:length(genome_unique)) {
    temp <- subset(genome, genome[,1]==chromosomes[i])
    genome_unique[i] <- max(temp[,2],na.rm = T)
  }
  genome <- as.data.frame(cbind(chromosomes,genome_unique));colnames(genome) <- c("chrom","size")
  genome$chrom <- as.character(genome$chrom)
  genome$size <- as.integer(as.character(genome$size))
  return(genome)
}

# regions
regions_fun <- function(regions_string) {
  regions <- read_bed(regions_string, n_fields = 6)
  regions <- regions[!grepl("track", regions$chrom),]
  regions <- within(regions, strand <- ifelse(is.na(strand),"+",strand))
  regions <- within(regions, score <- ifelse(is.na(score),0,score))
  regions <- within(regions, name <- ifelse(is.na(name),paste(chrom,start,end, sep = "_"),name))
  return(regions)
}
final_regions_TSS_or_TTS <- function(mode,length_before,length_after,regions) {
  if (mode=="TSS") {
    result <- within(regions, new_start <- ifelse(strand=="+", start-length_before,end-length_after))
    result <- within(result, new_end <- ifelse(strand=="+", start+length_after,end+length_before))
    result <- within(result, center <- ifelse(strand=="+",start,end))
    result <- result[,c(1,7,8,4,5,6,9)]; colnames(result) <- c(colnames(regions),"center")
  } else {
    result <- within(regions, new_start <- ifelse(strand=="+", end-length_before,start-length_after))
    result <- within(result, new_end <- ifelse(strand=="+", end+length_after,start+length_before))
    result <- within(result, center <- ifelse(strand=="+",end,start))
    result <- result[,c(1,7,8,4,5,6,9)]; colnames(result) <- c(colnames(regions),"center")
  }
  return(result)
}

# signal around coordinate
bin_vector_around_coordinate <- function(length_before,length_after,binsize) {
  return(seq(-length_before,length_after,by=binsize))
}
signal <- function(coverage, chrom_sizes, binsize) {
  binned_coverage <- bed_makewindows(coverage, chrom_sizes, win_size = binsize)
  binned_coverage <- binned_coverage[, c(1:ncol(binned_coverage)-1)]
  return(binned_coverage)
}
signal_around_coordinate <- function(final_regions,signal,bin_vector) {
  result <- bed_intersect(final_regions,signal)
  result <- within(result, start.y <- ifelse(start.y<start.x,start.x,start.y))
  result <- within(result, end.y <- ifelse(end.y>end.x,end.x,end.y))
  result <- within(result, start.y <- ifelse(strand.x=="+",start.y-center.x,center.x-start.y))
  result <- within(result, end.y <- ifelse(strand.x=="+",end.y-center.x,center.x-end.y))
  result <- within(result, coordinate <- (start.y+end.y)/2)
  result <- result[,c(4,10,12)]
  result$bin <- .bincode(result$coordinate, bin_vector, include.lowest = TRUE)
  result <- result[,c(1,2,4)]
  colnames(result) <- c("name", "bin_coverage","bin_number")
  return(result)
}

# matrices
matrix_if_TSS_or_TTS <- function(sample_list,mode,chrom_sizes,regions_string,length_before,length_after,binsize) {
  bin_vector_around_coordinate <- bin_vector_around_coordinate(length_before,length_after,binsize)
  regions <- regions_fun(regions_string)
  final_regions_TSS_or_TTS <- final_regions_TSS_or_TTS(mode,length_before,length_after,regions)
  result <- list()
  for (i in 1:nrow(sample_list)) {
    coverage <- coverage_import(sample_list[i,2])
    signal <- signal(coverage, chrom_sizes, binsize)
    temp <- signal_around_coordinate(final_regions_TSS_or_TTS,signal,bin_vector_around_coordinate)
    temp <- dcast(temp, name ~ bin_number, value.var = "bin_coverage", fun.aggregate = mean, na.rm = T)
    temp <- as.matrix(temp[,2:ncol(temp)])
    result[[i]] <- temp
  }
  return(result)
}
matrix_if_body <- function(sample_list,regions_string,chrom_sizes,L,binsize) {
  regions <- regions_fun(regions_string)
  bin_vector_body <- seq(-L/2,L/2,by=binsize)
  result <- list()  
  for (i in 1:nrow(sample_list)) {
    coverage <- coverage_import(sample_list[i,2])
    signal <- signal(coverage, chrom_sizes, binsize)
    signal_body <- bed_intersect(regions,signal)
    signal_body <- within(signal_body, start.y <- ifelse(start.y<start.x,start.x,start.y))
    signal_body <- within(signal_body, end.y <- ifelse(end.y>end.x,end.x,end.y))
    signal_body <- within(signal_body, center <- (start.x+end.x)/2)
    signal_body <- within(signal_body, operation <- ifelse(end.x-start.x<=L,"expand","compress"))
    signal_body_compress <- subset(signal_body,operation=="compress")
    signal_body_expand <- subset(signal_body,operation=="expand")
    signal_body_expand <- within(signal_body_expand, start.y <- ifelse(start.y<1/2*(start.x+end.x),
                                                                       start.y*1/start.x*(start.x-1/2*(L-(end.x-start.x))),
                                                                       start.y*1/end.x*(end.x+1/2*(L-(end.x-start.x)))))
    signal_body_expand <- within(signal_body_expand, end.y <- ifelse(end.y<1/2*(start.x+end.x),
                                                                     end.y*1/start.x*(start.x-1/2*(L-(end.x-start.x))),
                                                                     end.y*1/end.x*(end.x+1/2*(L-(end.x-start.x)))))
    signal_body_expand <- within(signal_body_expand, start.y <- round(start.y))
    signal_body_expand <- within(signal_body_expand, end.y <- round(end.y))
    signal_body_expand <- signal_body_expand[,c(1,7,8,4,9,6,11)]
    colnames(signal_body_expand) <- c("chrom", "start",  "end",    "name",   "score" , "strand", "center")
    signal_body_expand <- bed_makewindows(signal_body_expand, chrom_sizes, win_size = binsize)
    signal_body_expand <- within(signal_body_expand,coordinate <- ifelse(strand=="+", (start+end)/2-center,-(start+end)/2+center) )
    signal_body_expand <- signal_body_expand[,c(4,5,9)]
    signal_body_expand <- within(signal_body_expand, coordinate <- ifelse(coordinate<=-L/2,coordinate+1,ifelse(coordinate>=L/2,coordinate-1,coordinate) ))
    signal_body_expand$bin <- .bincode(signal_body_expand$coordinate, bin_vector_body, include.lowest = TRUE)
    signal_body_expand <- signal_body_expand[,c(1,2,4)]
    colnames(signal_body_expand) <- c("name", "bin_coverage","bin_number")
    signal_body_compress <- within(signal_body_compress, factor <- L/(end.x-start.x))
    signal_body_compress <- signal_body_compress[,c(1,7,8,4,9,6,11,13)]
    colnames(signal_body_compress) <- c("chrom", "start",  "end",    "name",   "score" , "strand", "center","factor")
    signal_body_compress <- within(signal_body_compress,coordinate <- ifelse(strand=="+", (start+end)/2-center,-(start+end)/2+center) )
    signal_body_compress <- within(signal_body_compress, coordinate <- factor * coordinate)
    signal_body_compress<-signal_body_compress[,c(4,5,9)]
    signal_body_compress <- within(signal_body_compress, coordinate <- ifelse(coordinate<=-L/2,coordinate+1,ifelse(coordinate>=L/2,coordinate-1,coordinate) ))
    signal_body_compress$bin <- .bincode(signal_body_compress$coordinate, bin_vector_body, include.lowest = TRUE)
    signal_body_compress <- signal_body_compress[,c(1,2,4)]
    colnames(signal_body_compress) <- c("name", "bin_coverage","bin_number")
    signal_body <- rbind(signal_body_compress,signal_body_expand)
    temp <- dcast(signal_body, name ~ bin_number, value.var = "bin_coverage", fun.aggregate = mean, na.rm = T)
    temp <- as.matrix(temp[,2:ncol(temp)])
    result[[i]] <- temp  
  }
  return(result)
}
matrix_if_all <- function(regions_string,chrom_sizes,sample_list,L,binsize,length_all_before,length_all_after) {
  sizes <- chrom_sizes
  regions <- regions_fun(regions_string)
  bin_vector_body <- seq(-L/2,L/2,by=binsize)
  bin_vector_before <- bin_vector_around_coordinate(length_all_before,0,binsize)
  final_regions_before <- final_regions_TSS_or_TTS("TSS",length_all_before,0,regions)
  bin_vector_after <- bin_vector_around_coordinate(0,length_all_after,binsize)
  final_regions_after <- final_regions_TSS_or_TTS("TTS",0,length_all_after,regions)
  result <- list() 
  for (i in 1:nrow(sample_list)) {
    coverage <- coverage_import(sample_list[i,2])
    signal <- signal(coverage, chrom_sizes, binsize)
    # body
    signal_body <- bed_intersect(regions,signal)
    signal_body <- within(signal_body, start.y <- ifelse(start.y<start.x,start.x,start.y))
    signal_body <- within(signal_body, end.y <- ifelse(end.y>end.x,end.x,end.y))
    signal_body <- within(signal_body, center <- (start.x+end.x)/2)
    signal_body <- within(signal_body, operation <- ifelse(end.x-start.x<=L,"expand","compress"))
    signal_body_compress <- subset(signal_body,operation=="compress")
    signal_body_expand <- subset(signal_body,operation=="expand")
    signal_body_expand <- within(signal_body_expand, start.y <- ifelse(start.y<1/2*(start.x+end.x),
                                                                       start.y*1/start.x*(start.x-1/2*(L-(end.x-start.x))),
                                                                       start.y*1/end.x*(end.x+1/2*(L-(end.x-start.x)))))
    signal_body_expand <- within(signal_body_expand, end.y <- ifelse(end.y<1/2*(start.x+end.x),
                                                                     end.y*1/start.x*(start.x-1/2*(L-(end.x-start.x))),
                                                                     end.y*1/end.x*(end.x+1/2*(L-(end.x-start.x)))))
    signal_body_expand <- within(signal_body_expand, start.y <- round(start.y))
    signal_body_expand <- within(signal_body_expand, end.y <- round(end.y))
    signal_body_expand <- signal_body_expand[,c(1,7,8,4,9,6,11)]
    colnames(signal_body_expand) <- c("chrom", "start",  "end",    "name",   "score" , "strand", "center")
    signal_body_expand <- bed_makewindows(signal_body_expand, sizes, win_size = binsize)
    signal_body_expand <- within(signal_body_expand,coordinate <- ifelse(strand=="+", (start+end)/2-center,-(start+end)/2+center) )
    signal_body_expand <- signal_body_expand[,c(4,5,9)]
    signal_body_expand <- within(signal_body_expand, coordinate <- ifelse(coordinate<=-L/2,coordinate+1,ifelse(coordinate>=L/2,coordinate-1,coordinate) ))
    signal_body_expand$bin <- .bincode(signal_body_expand$coordinate, bin_vector_body, include.lowest = TRUE)
    signal_body_expand <- signal_body_expand[,c(1,2,4)]
    colnames(signal_body_expand) <- c("name", "bin_coverage","bin_number")
    signal_body_compress <- within(signal_body_compress, factor <- L/(end.x-start.x))
    signal_body_compress <- signal_body_compress[,c(1,7,8,4,9,6,11,13)]
    colnames(signal_body_compress) <- c("chrom", "start",  "end",    "name",   "score" , "strand", "center","factor")
    signal_body_compress <- within(signal_body_compress,coordinate <- ifelse(strand=="+", (start+end)/2-center,-(start+end)/2+center) )
    signal_body_compress <- within(signal_body_compress, coordinate <- factor * coordinate)
    signal_body_compress<-signal_body_compress[,c(4,5,9)]
    signal_body_compress <- within(signal_body_compress, coordinate <- ifelse(coordinate<=-L/2,coordinate+1,ifelse(coordinate>=L/2,coordinate-1,coordinate) ))
    signal_body_compress$bin <- .bincode(signal_body_compress$coordinate, bin_vector_body, include.lowest = TRUE)
    signal_body_compress <- signal_body_compress[,c(1,2,4)]
    colnames(signal_body_compress) <- c("name", "bin_coverage","bin_number")
    signal_body <- rbind(signal_body_compress,signal_body_expand)
    #before
    before <- signal_around_coordinate(final_regions_before,signal,bin_vector_before)
    #after
    after <- signal_around_coordinate(final_regions_after,signal,bin_vector_after)
    #sum
    signal_body <- within(signal_body, bin_number <- max(before$bin_number, na.rm = TRUE) + bin_number)
    after <- within(after, bin_number <- max(signal_body$bin_number, na.rm = TRUE) + bin_number)
    temp <- rbind(before,signal_body,after)
    temp <- dcast(temp, name ~ bin_number, value.var = "bin_coverage", fun.aggregate = mean, na.rm = T)
    temp <- as.matrix(temp[,2:ncol(temp)])
    result[[i]] <- temp  
  }
  return(result)
}

# reordering and replacing NANs and NAs with 0
reorder <- function(matrices) {
  ref <- as.data.frame(matrices[[1]])
  ref <- within(ref, max <- rowMeans(ref, na.rm = T))
  ref <- within(ref, ID <- seq(1,nrow(ref)))
  ref <- ref[order(ref$max),]
  ref <- within(ref, order_column <- seq(1,nrow(ref)))
  ref <- ref[,c(ncol(ref)-1,ncol(ref))]
  result <- list()
  for (i in 1:length(matrices)) {
    temp <- as.data.frame(matrices[[i]])
    temp <- within(temp, ID <- seq(1,nrow(temp)))
    temp <- merge(ref,temp,by="ID")
    temp <- temp[order(temp$order_column),]
    temp <- temp[,3:ncol(temp)]
    temp <- as.matrix(temp)
    temp[is.nan(temp)] = 0
    temp[is.na(temp)] = 0
    result[[ i ]] <- temp
  }
  return(result)
}
matrix_fun <- function(sample_list,regions_string, mode,binsize, length_before,length_after, L,  length_all_before,length_all_after,chrom_sizes) {
  if (mode=="TSS") {
    matrices <- matrix_if_TSS_or_TTS(sample_list,"TSS",chrom_sizes,regions_string ,length_before,length_after,binsize)
  }
  if (mode=="TTS") {
    matrices <- matrix_if_TSS_or_TTS(sample_list,"TTS",chrom_sizes,regions_string ,length_before,length_after,binsize)
  }
  if (mode=="body") {
    matrices <- matrix_if_body(sample_list,regions_string,chrom_sizes,L,binsize)
  }
  if (mode=="body_plus") {
    matrices <- matrix_if_all(regions_string,chrom_sizes,sample_list,L,binsize,length_all_before,length_all_after)
  }
  matrices <- reorder(matrices)
  return(matrices)
}

# plotting
multiplot_enriched_heatmaps <- function(sample_list, matrices, mode, binsize, length_before,length_after, L,  length_all_before,length_all_after) {
  if (mode=="TSS") {axis_name = c(paste("-",length_before,"bp",sep = " "), "start",paste(length_after,"bp",sep = " "))}
  if (mode=="TTS") {axis_name = c(paste("-",length_before,"bp",sep = " "), "end",paste(length_after,"bp",sep = " "))} 
  if (mode=="body") {axis_name = c("start","end")} 
  if (mode=="body_plus") {axis_name = c(paste("-",length_all_before,"bp",sep = " "),"start","end",paste(length_all_after,"bp",sep = " "))} 
  mat1 <- matrices[[1]]
  if (mode=="TSS"|mode=="TTS") {
    attr(mat1, "upstream_index") = 1:(length_before/binsize+0.5)
    attr(mat1, "target_index") = integer(0)
    attr(mat1, "downstream_index") =  ((length_before/binsize)+0.5):( (length_before+length_after)/binsize)
  } 
  if (mode=="body") {
    attr(mat1, "upstream_index") = integer(0)
    attr(mat1, "target_index") = 1:(L/binsize)
    attr(mat1, "downstream_index") =  integer(0)
  } 
  if (mode=="body_plus") {
    attr(mat1, "upstream_index") = 1:(length_all_before/binsize+0.5)
    attr(mat1, "target_index") = (length_all_before/binsize+0.5):((length_all_before+L)/binsize+0.5)
    attr(mat1, "downstream_index") =  ((length_all_before+L)/binsize+0.5):((length_all_before+L+length_all_after)/binsize+0.5)
  }
  class(mat1) = c("normalizedMatrix", "matrix")
  if (min(mat1,na.rm = T)<0) {
    col_fun = colorRamp2(c(quantile(mat1, 0.01),0,quantile(mat1, 0.99)), c("lightblue1", "white", "red"))
  } else {
    col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
  }
  ht_global_opt(heatmap_legend_labels_gp = gpar(fontsize = 15), heatmap_legend_title_gp = gpar(fontsize = 18), heatmap_legend_grid_height = unit(0.5, "cm"))
  mat_list <- list()
  for (i in 1:length(matrices)) {
    mat_temp <- matrices[[i]]
    if (mode=="TSS"|mode=="TTS") {
      attr(mat_temp, "upstream_index") = 1:(length_before/binsize+0.5)
      attr(mat_temp, "target_index") = integer(0)
      attr(mat_temp, "downstream_index") =  ((length_before/binsize)+0.5):( (length_before+length_after)/binsize)
    } 
    if (mode=="body") {
      attr(mat_temp, "upstream_index") = integer(0)
      attr(mat_temp, "target_index") = 1:(L/binsize)
      attr(mat_temp, "downstream_index") =  integer(0)
    } 
    if (mode=="body_plus") {
      attr(mat_temp, "upstream_index") = 1:(length_all_before/binsize+0.5)
      attr(mat_temp, "target_index") = (length_all_before/binsize+0.5):((length_all_before+L)/binsize+0.5)
      attr(mat_temp, "downstream_index") =  ((length_all_before+L)/binsize+0.5):((length_all_before+L+length_all_after)/binsize+0.5)
    } 
    class(mat_temp) = c("normalizedMatrix", "matrix")
    mat_list[[i]] <- mat_temp
    if (i==1) {
      ht_list <- EnrichedHeatmap(mat_list[[i]], col = col_fun, column_title= as.character(sample_list[i,1]), name =  " ", axis_name = axis_name, axis_name_gp = gpar(fontsize = 18), top_annotation = HeatmapAnnotation(lines = anno_enriched()), top_annotation_height = unit(2, "cm"), axis_name_rot = 90) 
    } else {
      ht_list <- ht_list + EnrichedHeatmap(mat_list[[i]], col = col_fun, column_title=as.character(sample_list[i,1]),name = as.character(sample_list[i,1]), axis_name = axis_name, axis_name_gp = gpar(fontsize = 18), top_annotation = HeatmapAnnotation(lines = anno_enriched()), top_annotation_height = unit(2, "cm"),show_heatmap_legend = FALSE, axis_name_rot = 90) 
    }
  }
  print(draw(ht_list, gap = unit(rep(12,nrow(sample_list)-1), "mm")))
}
ci95perc_min <- function(x) {mean(x) - qnorm(0.975)*sd(x)/sqrt(length(x))}
ci95perc_max <- function(x) {mean(x) + qnorm(0.975)*sd(x)/sqrt(length(x))}
metaplot <- function(sample_list, matrices, mode, binsize, length_before,length_after, L,  length_all_before,length_all_after) {
  for (i in 1:length(matrices)) {
    if (i==1) {
      new_mat <- matrices[[i]]
      new_mat <- data.frame(matrix(0, ncol = ncol(matrices[[i]]), nrow = 4))
      colnames(new_mat) <- colnames(matrices[[i]])
      rownames(new_mat) <- c("position","average","conf_minus","conf_plus")
      new_mat[1,] <- colnames(new_mat)
      new_mat[2,] <- apply(matrices[[i]], 2, mean)
      new_mat[3,] <- apply(matrices[[i]], 2, ci95perc_min)
      new_mat[4,] <- apply(matrices[[i]], 2, ci95perc_max)
      new_mat <- as.data.frame(t(new_mat))
      new_mat <- within(new_mat, SampleID <- as.character(sample_list[i,1]))
    } else {
      temp <- matrices[[i]]
      temp <- data.frame(matrix(0, ncol = ncol(matrices[[i]]), nrow = 4))
      colnames(temp) <- colnames(matrices[[i]])
      rownames(temp) <- c("position","average","conf_minus","conf_plus")
      temp[1,] <- colnames(temp)
      temp[2,] <- apply(matrices[[i]], 2, mean)
      temp[3,] <- apply(matrices[[i]], 2, ci95perc_min)
      temp[4,] <- apply(matrices[[i]], 2, ci95perc_max)
      temp <- as.data.frame(t(temp))
      temp <- within(temp, SampleID <- as.character(sample_list[i,1]))
      new_mat <- rbind(new_mat,temp)
    }
  }
  new_mat$position <- factor(new_mat$position, levels = colnames(matrices[[1]]))
  new_mat$SampleID <- factor(new_mat$SampleID, levels = unique(as.character(sample_list[,1])))
  new_mat$average <- as.numeric(as.character(new_mat$average))
  new_mat$conf_minus <- as.numeric(as.character(new_mat$conf_minus))
  new_mat$conf_plus <- as.numeric(as.character(new_mat$conf_plus))
  colvector <- as.character(sample_list[,3])
  xintercept1 <- 0
  xintercept2 <- 0
  xintercept3 <- 0
  segm_offset <- (max(new_mat$conf_plus)-min(new_mat$conf_minus))/50
  if (mode=="body") {alpha1<- 0; alpha2<- 0; xintercept3 <- as.integer(L/binsize/2)}
  if (mode=="body_plus") {alpha1<- 1; alpha2<- 1; xintercept3 <- as.integer(0.05*(length_all_before+length_all_after+L)/binsize)}
  if (mode=="TSS"|mode=="TTS") {alpha1<- 1; alpha2<- 0; xintercept1 <- length_before/binsize+0.5; xintercept3 <- as.integer(0.05*(length_before+length_after)/binsize)}
  if (mode=="body_plus") {xintercept1 <- length_all_before/binsize + 0.5; xintercept2 <- (length_all_before+L)/binsize + 0.5}
  p <- ggplot(new_mat,aes(x=position, y=average,group=SampleID)) + geom_line(aes(color=SampleID)) + geom_ribbon(aes(ymin=conf_minus, ymax=conf_plus, fill = SampleID), alpha=0.2)+ scale_color_manual(values=colvector) + scale_fill_manual(values=colvector) + labs(x="",y="Coverage average per bin") + theme_bw() +
    theme(panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(), legend.key=element_blank(),
          panel.grid.minor = element_blank(), legend.title =  element_blank(),legend.text = element_text(size=16),
          axis.text.x = element_blank(),axis.text.y = element_text(size=16), axis.ticks.x=element_blank(),
          axis.title.x =  element_blank(),axis.title.y = element_text(size=16,face = "bold")) + geom_vline(xintercept=xintercept1, alpha=alpha1,linetype="dotted") + geom_vline(xintercept=xintercept2, alpha=alpha2,linetype="dotted") +
    annotate("text", x = c(xintercept1,xintercept2,xintercept3), fontface=c(2,2,1) , y = c(min(new_mat$conf_minus),min(new_mat$conf_minus),min(new_mat$conf_minus)),  alpha=c(alpha1,alpha2,1), label = c("Start","End",paste(binsize*2,"bp")),  angle = c(90,90,0)) +
    geom_segment(x = xintercept3-1, y = min(new_mat$conf_minus)+segm_offset, xend = xintercept3+1, yend = min(new_mat$conf_minus)+segm_offset, colour = "black")
  print(p)
}