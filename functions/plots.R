MDS_plot <- function(TMM,sample_list) {
  p <- plotMDS(TMM)
  data <- within(as.data.frame(p$cmdscale.out), condition <- as.character(TMM$samples$group))
  data$condition <- factor(data$condition, levels = unique(as.character(data$condition)))
  sample_list_colors <- sample_list[,(ncol(sample_list)-1):ncol(sample_list)]
  sample_list_colors <- sample_list_colors[!duplicated(sample_list_colors),]
  sample_list_colors <- within(sample_list_colors, control_color <- sample_list_colors[,ncol(sample_list_colors)])
  sample_list_colors <- sample_list_colors[,(ncol(sample_list_colors)-1):ncol(sample_list_colors)]
  sample_list_colors <- c(t(sample_list_colors))
  shape <- rep(c(16,1),length(sample_list_colors)/2)
  keysize <- rep(5,length(shape))
  plot <- ggplot(data, aes(x=V1, y=V2,group=condition)) +  geom_point(aes(shape=condition, color=condition, size=condition)) +
    scale_shape_manual(values=shape) + scale_color_manual(values=sample_list_colors) + labs(x="Dimension 1",y="Dimension 2") + scale_size_manual(values=keysize) +
    theme_bw() + theme(panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), legend.title =  element_blank(),legend.text = element_text(size=16),
                       axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),
                       axis.title.x = element_text(size=16,face = "bold"),axis.title.y = element_text(size=16,face = "bold"))
  print(plot)
}



diff_result_plot <- function(results,FDR_threshold, string_test,string_ref) {
  diff_result_table <- subset(results, FDR!="Below treatment/control threshold")
  for(i in 4:8) {
    diff_result_table[,i] <- as.numeric(as.character(diff_result_table[,i]))
  }
  diff_result_table <-  within(diff_result_table, sign <- ifelse(FDR<FDR_threshold,1,0))
  diff_result_table <-  within(diff_result_table, class <- ifelse(sign==0,"not changing",ifelse(diff_result_table[,8]>0,"enriched",ifelse(diff_result_table[,8]<0,"depleted","not changing"))))
  diff_result_table$class <- factor(diff_result_table$class, levels = c("enriched", "depleted", "not changing"))
  diff_result_table$class[is.na(diff_result_table$class)] <- "not changing"
  colvector <- c("red","blue","black") 
  if (nrow(subset(diff_result_table, class=="enriched"))==0) {
    if (nrow(subset(diff_result_table, class=="depleted"))==0) {
      colvector <-"black"  
    } else {
      if (nrow(subset(diff_result_table, class=="not changing"))==0) {
        colvector <- "blue"
      } else {
        colvector <- c("blue","black")
      }
    }
  }
  if (nrow(subset(diff_result_table, class=="depleted"))==0) {
    if (nrow(subset(diff_result_table, class=="enriched"))==0) {
      colvector <-"black"  
    } else {
      if (nrow(subset(diff_result_table, class=="not changing"))==0) {
        colvector <- "red"
      } else {
        colvector <- c("red","black")
      }
    }
  }
  if (nrow(subset(diff_result_table, class=="not changing"))==0) {
    if (nrow(subset(diff_result_table, class=="enriched"))==0) {
      colvector <-"blue"  
    } else {
      if (nrow(subset(diff_result_table, class=="depleted"))==0) {
        colvector <- "red"
      } else {
        colvector <- c("red","blue")
      }
    }
  }
  colnames(diff_result_table)[8] <- "FC"
  diff_result_table <- within(diff_result_table, Average_CPM <- (diff_result_table[,6] + diff_result_table[,7])/2)
  p <- ggplot(diff_result_table,aes(x=Average_CPM, y=FC,group=class)) + geom_point(aes(color=class)) + scale_color_manual(values=colvector)  + labs(x=paste("Average Log2 CPM",string_test, "&",string_ref, sep = " "),y=paste("Log2 Fold Enrichment",string_test, "vs",string_ref, sep = " ")) + theme_bw() +
    theme(panel.border = element_rect(colour = "black"), panel.grid.major = element_line(colour = "gray88"),
          panel.grid.minor = element_blank(), legend.title =  element_blank(),legend.text = element_text(size=16),
          axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),
          axis.title.x = element_text(size=16,face = "bold"),axis.title.y = element_text(size=16,face = "bold"))
  print(p)
} 
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = T)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
boxplot_Q95 <- function(results,string_test,string_ref) {
  diff_result_table <- subset(results, FDR!="Below treatment/control threshold")
  diff_result_table_test <- cbind(diff_result_table[,6],rep(string_test,nrow(diff_result_table)))
  diff_result_table_ref <- cbind(diff_result_table[,7],rep(string_ref,nrow(diff_result_table)))
  diff_result_table <- as.data.frame(rbind(diff_result_table_test,diff_result_table_ref))
  colnames(diff_result_table) <- c("CPM", "condition" )
  diff_result_table$condition <- factor(diff_result_table$condition, levels = c(string_ref, string_test))
  diff_result_table$CPM <- as.numeric(as.character(diff_result_table$CPM))
  p <- ggplot(diff_result_table,aes(x=condition, y=CPM)) + stat_summary(fun.data = quantiles_95, geom="boxplot") +  labs(x="",y="Median Log2 CPM") + scale_y_continuous(expand=c(0.1,0)) +
    theme_bw() + theme(panel.border = element_rect(colour = "black"), panel.grid.major = element_line(colour = "gray88"),panel.grid.minor = element_blank(),
                       axis.text.x = element_text(size=16),axis.text.y = element_text(size=16), legend.title =  element_blank(),legend.text = element_text(size=16),
                       axis.title.x = element_text(size=16,face = "bold"),axis.title.y = element_text(size=16,face = "bold")) +
    geom_signif(test= 'wilcox.test',comparisons = list(c(string_test,string_ref)), map_signif_level=TRUE, textsize = 10, vjust=0.25, y_position = (quantiles_95(diff_result_table$CPM)[5]*1.2), tip_length = 0)
  print(p)
  return(diff_result_table)
}

