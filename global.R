list_of_packages = c("shiny","shinyFiles","BSgenome","QuasR","csaw","edgeR","GenomicAlignments","BayesPeak","ChIPseeker",
                     "GenomicFeatures","ggplot2","ggsignif","valr","reshape2","circlize","EnrichedHeatmap")

lapply(list_of_packages, 
       function(x) if(!require(x,character.only = TRUE)) install.packages(x, dependencies = TRUE))
lapply(list_of_packages, 
       function(x) if(!require(x,character.only = TRUE)) {
         source("https://bioconductor.org/biocLite.R")
         biocLite(x, suppressUpdates=TRUE, ask=FALSE)
       } 
       )

library(shiny)
library(shinyFiles)
library(BSgenome) 
library(QuasR) 
library(csaw)
library(edgeR)
library(GenomicAlignments)
library(BayesPeak)
library(ChIPseeker)
library(GenomicFeatures)
library(ggplot2)
library(ggsignif)
library(valr)
library(reshape2) 
library(circlize)
library(EnrichedHeatmap)

source("./functions/mapped_reads_processing.R")
source("./functions/plots.R")
source("./functions/peak_operations.R")
source("./functions/mapping.R")
source("./functions/heatmaps_and_metaplot.R")