# ChIPdig

ChIPdig is a tool designed for the bulk analysis of ChIP-seq data comprising multiple samples. Its capabilities are organized into four analysis modules:
•	Read alignment 
•	Normalization and comparison of ChIP-seq data sets
•	Annotation of genomic regions
•	Visualization of coverage (heatmap and metaplot generation)
To be able to use the tool, you will first need to download the folder to your computer, and then either manually install the necessary packages or let the application do it for you when executed for the first time. Bioconductor packages and their dependencies can be installed by typing the following commands in the R studio console:  
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome")

In the example above, the first command establishes communication with the Bioconductor repository and the second command installs package “edgeR”. Repeat the second command for the remaining Bioconductor packages used by ChIPdig:
biocLite("QuasR")
biocLite("csaw")
biocLite("edgeR")
biocLite("GenomicAlignments")
biocLite("BayesPeak")
biocLite("ChIPseeker")
biocLite("GenomicFeatures")
biocLite("EnrichedHeatmap")

Other necessary packages may be installed with these commands:
install.packages("shiny")
install.packages("shinyFiles")
install.packages("ggplot2")
install.packages("ggsignif")
install.packages("reshape2")
install.packages("valr")
install.packages("circlize")

To load ChIPdig, simply navigate to the folder downloaded from the repository and open any of the following files using RStudio: “server.R”, “global.R”, or “ui.R”. Then, in the RStudio console, click on the “Run App” button. If any of the packages necessary to run ChIPdig is not installed, this will be done automatically.  
A tutorial (ChIPdig_tutorial.pdf) explaining how to use the tool is included, along with a folder containing toy data enabling the user to easily test ChIPdig before analyzing his own samples.

Please contact me for questions and suggestions: rmesse@bu.edu. 
