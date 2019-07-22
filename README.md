 10.5281/zenodo.3345788 

# ChIPdig

ChIPdig is a tool designed for the bulk analysis of ChIP-seq data comprising multiple samples. Its capabilities are organized into four analysis modules: <br>

•	Read alignment <br>
•	Normalization and comparison of ChIP-seq data sets <br>
•	Annotation of genomic regions <br>
•	Visualization of coverage (heatmap and metaplot generation) <br>

To be able to use the tool, you will first need to download the folder to your computer, and then either manually install the necessary packages or let the application do it for you when executed for the first time. Bioconductor packages and their dependencies can be installed by typing the following commands in the R studio console:  

•	source("https://bioconductor.org/biocLite.R") <br>
•	biocLite("BSgenome") <br>

In the example above, the first command establishes communication with the Bioconductor repository and the second command installs package “edgeR”. Repeat the second command for the remaining Bioconductor packages used by ChIPdig:

•	biocLite("QuasR") <br>
•	biocLite("csaw") <br>
•	biocLite("edgeR") <br>
•	biocLite("GenomicAlignments") <br>
•	biocLite("BayesPeak") <br>
•	biocLite("ChIPseeker") <br>
•	biocLite("GenomicFeatures") <br>
•	biocLite("EnrichedHeatmap") <br>

Other necessary packages may be installed with these commands:

•	install.packages("shiny") <br>
•	install.packages("shinyFiles") <br>
•	install.packages("ggplot2") <br>
•	install.packages("ggsignif") <br>
•	install.packages("reshape2") <br>
•	install.packages("valr") <br>
•	install.packages("circlize") <br>

To load ChIPdig, simply navigate to the folder downloaded from the repository and open any of the following files using RStudio: “server.R”, “global.R”, or “ui.R”. Then, in the RStudio console, click on the “Run App” button. If any of the packages necessary to run ChIPdig is not installed, this will be done automatically.  

A tutorial (ChIPdig_tutorial.pdf) explaining how to use the tool is included, along with a folder containing toy data enabling the user to easily test ChIPdig before analyzing his own samples.

Please contact me for questions and suggestions: rmesse@bu.edu. 
