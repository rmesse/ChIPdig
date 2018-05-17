shinyServer(function(input,output,session){
  
  options(warn=-1)
  
  # input directory
  initialdir <- getwd()
  volumes = getVolumes()
  shinyDirChoose(input, 'dir', roots=volumes, session=session, restrictions=system.file(package='base'))
  output$dir <- renderPrint(initialdir)
  observeEvent(input$dir, {
    output$dir <- renderPrint(parseDirPath(volumes, input$dir))
  })
  directory <- reactive(input$dir)
  
  # initial instruction
  output$text1 <- renderUI(h3("ChIPdig provides a user-friendly interface for analyzing multi-sample ChIP-seq data using using packages implemented in R. First, select the analysis module in the left pane and the input directory containing the files to be analyzed."))
  
  # loading of user input options
  observeEvent(input$first_choice, {
    if(input$first_choice=="1" | input$first_choice=="2") {
      output$second_choice <- renderUI(radioButtons("second_choice", label = "Is data single-ended or paired-ended?", choices = list("Single-ended" = 1, "Paired-ended" = 2), selected = 1))
    } else {
      output$second_choice <- renderUI(NULL)
    }
    if(input$first_choice=="2") {
      updateTabsetPanel(session, "inTabset",selected = "initial_proc")
      output$text8 <- renderUI( h5("Processing of mapped reads requires a tab-delimited text file listing mapped read files in BAM format. Each BAM file must be accompanied by a matching indexing file (ending in '.bai'). The text file should not have column names and should be in the format specified in the following example:"))
      data_frame_processing_example <- data.frame(a=c("H3K4me3_wildtype_replicate_1","H3K4me3_wildtype_replicate_2","H3K4me3_mutant_replicate_1","H3K4me3_mutant_replicate_2") , b=c("H3K4me3_wildtype_replicate_1_treatment.bam","H3K4me3_wildtype_replicate_2_treatment.bam","H3K4me3_mutant_replicate_1_treatment.bam","H3K4me3_mutant_replicate_2_treatment.bam"), c=c("H3K4me3_wildtype_replicate_1_control.bam","H3K4me3_wildtype_replicate_2_control.bam","H3K4me3_mutant_replicate_1_control.bam","H3K4me3_mutant_replicate_2_control.bam"),d=c("H3K4me3_wildtype","H3K4me3_wildtype","H3K4me3_mutant","H3K4me3_mutant"),e=c("blue","blue","red","red"))
      rv$data_frame_processing_example <- data_frame_processing_example
      output$data_frame_processing_example <- renderTable(data_frame_processing_example, include.colnames=FALSE)      
      output$BAM_sample_processing_title <- renderUI(h3("Processing of mapped reads"))
      output$file_input <- renderUI(fileInput('file_input', 'Select tab-delimited text file listing aligned read files (see instructions):',multiple = FALSE, accept='text/comma-separated-values,text/plain'))
      output$genome <- renderUI(NULL)
      output$anno_button <- renderUI(NULL)
      output$TSS_def <- renderUI(NULL)
      output$anno_title <- renderUI(NULL)
      output$align_title  <- renderUI(NULL)
      output$heatmaps_title <- renderUI(NULL)
      output$heatmaps_regions <- renderUI(NULL)
      output$heatmaps_mode <- renderUI(NULL)
      output$heatmaps_binsize <- renderUI(NULL)
      output$length_before <- renderUI(NULL)
      output$length_after <- renderUI(NULL)
      output$heatmaps_button <- renderUI(NULL)
    } else {
      if(input$first_choice=="3") {
        updateTabsetPanel(session, "inTabset",selected = "anno")
        output$text7 <- renderUI( h5("Load the peak file in BED format. If the file has a track definition line and column names, these have to be deleted first."))
        output$text1 <- renderUI(NULL)
        output$file_input <- renderUI(fileInput('file_input', 'Select  peak file in BED format:',multiple = FALSE, accept='text/comma-separated-values,text/plain'))
        output$genome <- renderUI(selectInput("genome", "Choose reference genome:", ucscGenomes()[,"db"], multiple=F))
        output$anno_button <- renderUI(actionButton("anno_button", "Annotate peaks"))
        output$TSS_def <- renderUI(sliderInput("TSS_def", label = "Distance upstream (choose negative value) and downstream (choose positive value) of annotated transcription start site to consider:", min = -5000, max = 5000, value = c(-1500, 500)))
        output$anno_title <- renderUI(h3("Annotation of peaks"))
        output$BAM_sample_processing_title <- renderUI(NULL)
        output$align_title  <- renderUI(NULL)
        output$heatmaps_title <- renderUI(NULL)
        output$heatmaps_regions <- renderUI(NULL)
        output$heatmaps_mode <- renderUI(NULL)
        output$heatmaps_binsize <- renderUI(NULL)
        output$length_before <- renderUI(NULL)
        output$length_after <- renderUI(NULL)
        output$heatmaps_button <- renderUI(NULL)
      } else {
        if(input$first_choice=="1") {
          output$align_title <- renderUI(h3("Alignment of reads to reference genome"))
          output$file_input <- renderUI(fileInput('file_input', 'Select text file listing sequence files (see instructions):',multiple = FALSE, accept='text/comma-separated-values,text/plain'))
          updateTabsetPanel(session, "inTabset",selected = "map")
          output$text1 <- renderUI(NULL)
          output$text2 <- renderUI( h5("The alignment module requires a tab-delimited text file listing read files. Files ending with '.fq', '.fastq', '.gz', '.bz2' and 'xz' are supported. For a single-read experiment, the text file is merely a list of the names of files to be analyzed. Do not assign column names. Here's an example:"))
          singleread_exp_input_example <- data.frame(a=c("H3K4me3_IP_replicate_1.fq.gz2","H3K4me3_input_replicate_1.fq.gz2","H3K4me3_IP_replicate_2.fq.gz2","H3K4me3_input_replicate_2.fq.gz2") , c=5)
          singleread_exp_input_example <- singleread_exp_input_example[,1, drop=F]
          output$singleread_exp_input_example <- renderTable(singleread_exp_input_example, include.colnames=FALSE)
          output$text3 <- renderUI(h5("For a paired-end experiment:"))
          pairedend_exp_input_example <- data.frame(a=c("H3K4me3_IP_replicate_1_1.fq.gz2","H3K4me3_input_replicate_1_1.fq.gz2","H3K4me3_IP_replicate_2_1.fq.gz2","H3K4me3_input_replicate_2_1.fq.gz2") , a=c("H3K4me3_IP_replicate_1_2.fq.gz2","H3K4me3_input_replicate_1_2.fq.gz2","H3K4me3_IP_replicate_2_2.fq.gz2","H3K4me3_input_replicate_2_2.fq.gz2") )
          output$pairedend_exp_input_example <- renderTable(pairedend_exp_input_example, include.colnames=FALSE)
          output$text4 <- renderUI(h5("Select \"single-ended\" or \"paired-ended\" and load the tab-delimited text file."))
          output$anno_button <- renderUI(NULL)
          output$TSS_def <- renderUI(NULL)
          output$anno_title <- renderUI(NULL)
          output$BAM_sample_processing_title <- renderUI(NULL)
          output$genome <- renderUI(NULL)
          output$heatmaps_title <- renderUI(NULL)
          output$heatmaps_regions <- renderUI(NULL)
          output$heatmaps_mode <- renderUI(NULL)
          output$heatmaps_binsize <- renderUI(NULL)
          output$length_before <- renderUI(NULL)
          output$length_after <- renderUI(NULL)
          output$heatmaps_button <- renderUI(NULL)
        } else {
          updateTabsetPanel(session, "inTabset",selected = "heatmaps")
          output$text1 <- renderUI(NULL)
          output$text6 <- renderUI( h5("Visualization of coverage in specific regions requires a corresponding file in BED format. At least 3 columns (chromosome, start and end) are needed. If regions have strand specificity, such information should be provided in the 6-th column. Do not assign column names. Load the regions file using the corresponding button in the left pane. In addition, a tab-delimited text file listing coverage files in bedGraph format is required. Do not assign column names. Use the format specified in the following example:"))
          bedgraph_tab_example <- data.frame(a=c("condition_1_replicate_1","condition_1_replicate_2","condition_2_replicate_1","condition_2_replicate_2"),b=c("condition_1_replicate_1.bedgraph","condition_1_replicate_2.bedgraph","condition_2_replicate_1.bedgraph","condition_2_replicate_2.bedgraph"),c=c("blue","blue","red","red"))
          output$bedgraph_tab_example <- renderTable(bedgraph_tab_example, include.colnames=FALSE)
          output$file_input <- renderUI(fileInput('file_input', 'Select  tab-delimited text file listing coverage files (see instructions):',multiple = FALSE, accept='text/comma-separated-values,text/plain'))
          output$heatmaps_title <- renderUI(h3("Visualization of coverage in specific regions"))
          output$heatmaps_regions <- renderUI(fileInput('heatmaps_regions', 'Select  peak file in BED format (do not include header/track definition line):',multiple = FALSE, accept='text/comma-separated-values,text/plain'))
          output$heatmaps_mode <- renderUI(selectInput("heatmaps_mode", "Choose reference coordinate(s):", c("", "start","end","both"), multiple=F))
          output$heatmaps_binsize <- renderUI(numericInput("heatmaps_binsize", "Bin size (bp):", 50, min = 10))
          output$length_before <- renderUI(numericInput("length_before", "Length upstream to show (bp):", 500, min = 0))
          output$length_after <- renderUI(numericInput("length_after", "Length downstream to show (bp):", 500, min = 0))
          output$genome <- renderUI(NULL)
          output$anno_button <- renderUI(NULL)
          output$TSS_def <- renderUI(NULL)
          output$anno_title <- renderUI(NULL)
          output$BAM_sample_processing_title <- renderUI(NULL)
          output$align_title  <- renderUI(NULL)
        }
      }
    }
  })
  
  # reactive values
  rv <- reactiveValues()
  
  # sample loading
  observeEvent(input$file_input, {
    if (input$first_choice == "2") {
      updateTabsetPanel(session, "inTabset",selected = "initial_proc")
      output$text8 <- renderUI(NULL)
      output$data_frame_processing_example <- renderUI(NULL)
      if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
      inFile <- input$file_input
      sample_list_path <- inFile$datapath
      sample_list <- sample_list_import(sample_list_path)
      if (sample_list=="ncol_not_5") {
        output$warning1 <- renderText("Error loading sample list. Tab-delimited text file must have 5 columns.")
        output$warning2 <- renderUI( h5("Processing of mapped reads requires a tab-delimited text file listing mapped read files in BAM format. Each BAM file must be accompanied by a matching indexing file (ending in '.bai'). The text file should not have column names and should be in the format specified in the following example:"))
        output$warning3 <- renderTable(rv$data_frame_processing_example, include.colnames=FALSE) 
      } 
      if (sample_list=="missing_values") {
        output$warning1 <- renderText("Error loading sample list. Missing values. All fields in the tab-delimited text file describing samples must be non-empty.")
        output$warning2 <- renderUI( h5("Processing of mapped reads requires a tab-delimited text file listing mapped read files in BAM format. Each BAM file must be accompanied by a matching indexing file (ending in '.bai'). The text file should not have column names and should be in the format specified in the following example:"))
        output$warning3 <- renderTable(rv$data_frame_processing_example, include.colnames=FALSE) 
      } 
      if (sample_list=="bam_missing") {
        output$warning1 <- renderText("Error loading sample list. At least one file not found in the input directory. Make sure to provide the correct file names.")
        output$warning2 <- renderUI( h5("Processing of mapped reads requires a tab-delimited text file listing mapped read files in BAM format. Each BAM file must be accompanied by a matching indexing file (ending in '.bai'). The text file should not have column names and should be in the format specified in the following example:"))
        output$warning3 <- renderTable(rv$data_frame_processing_example, include.colnames=FALSE) 
      } 
      if (sample_list!="ncol_not_5" & sample_list!="missing_values" & sample_list!="bam_missing") {
        output$sample_list <- renderTable(sample_list)
        rv$sample_list <- sample_list
        output$dupl_reads_remove <- renderUI(checkboxInput("dupl_reads_remove", "Remove duplicate reads" , value = FALSE))
        output$process_reads <- renderUI(actionButton("process_reads", "Load mapped read files"))
        output$process_reads_note <- renderUI("NOTE: Loading may take a while. Be patient.")
        output$bin_size <- renderUI(numericInput("bin_size", "Bin size (bp):", 50, min = 25)) 
        setwd(initialdir)
        if (input$second_choice == "1") {
          output$second_choice <- renderUI(NULL)
          output$sample_list_title <- renderUI(h3("Sample list (single-end reads):"))
          rv$single_or_paired <- "single"
          output$read_extension_single <- renderUI(radioButtons("read_extension_single", label = NULL, choices = list("Extend reads to computationally estimated median fragment length" = 1, "Extend reads to experimentally observed median fragment length" = 2), selected = FALSE))
          output$max_distance_pairs <- renderUI(NULL)
        } else {
          output$second_choice <- renderUI(NULL)
          output$sample_list_title <- renderUI(h3("Sample list (paired-end reads):"))
          rv$single_or_paired <- "paired"
          output$strandMode <- renderUI(radioButtons("strandMode", label = NULL, choices = list("Strand of the pair is always +." = 0, "Strand of the pair is strand of its first alignment. This mode should be used when the paired-end data was generated using one of the following stranded protocols: Directional Illumina (Ligation), Standard SOLiD." = 1, "Strand of the pair is strand of its last alignment. This mode should be used when the paired-end data was generated using one of the following stranded protocols: dUTP, NSR, NNSR, Illumina stranded TruSeq PE protocol." = 2), selected = 1))      
          output$read_extension_single <- renderUI(NULL)
          output$max_distance_pairs <- renderUI(numericInput("max_distance_pairs", "For differential peak calling, only consider paired reads separated by distance less than (bp):", 400, min = 0)) 
        }
        output$file_input <- renderUI(NULL)
      }
    }
    if (input$first_choice == "1") {
      updateTabsetPanel(session, "inTabset",selected = "map")
      if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
      
      inFile <- input$file_input
      path <- inFile$datapath
      fastq_list <- fastq_list_import(path,input$second_choice)
      
      
      if (fastq_list=="ncol_not_1") {
        output$warning_mapping <- renderText("Error loading sample list. Text file must have only 1 columns.")
        
      }
      
      if (fastq_list!="ncol_not_1" ) {
        output$text2 <- renderUI(NULL); output$text3 <- renderUI(NULL); output$text4 <- renderUI(NULL)
        output$singleread_exp_input_example <- renderUI(NULL);output$pairedend_exp_input_example <- renderUI(NULL)
        output$text5 <- renderUI(h5("Read preprocessing can be used to prepare the input sequence files prior to alignment. Removal of sequence segnments corresponding to adapters and full low quality reads is supported. Preprocessing is not required."))
        output$preprocessing_title <- renderUI(h4("Pre-processing"))
        output$preprocessing_button <- renderUI(actionButton("preprocessing_button", "Pre-process reads"))
        output$preprocessing_no_button <- renderUI(actionButton("preprocessing_no_button", "Skip read pre-processing"))
        output$truncateEndBases <- renderUI(numericInput("truncateEndBases", "Number of bases to be removed from the end of each sequence:", 0, min = 0))
        output$truncateStartBases <- renderUI(numericInput("truncateStartBases", "Number of bases to be removed from the beginning of each sequence:", 0, min = 0))
        if (input$second_choice == "1") {output$Lpattern <- renderUI(textInput("Lpattern", "Left (5'-end) adapter sequence (leave empty if not to trim):",""))} 
        if (input$second_choice == "1") {output$Rpattern <- renderUI(textInput("Rpattern", "Right (3'-end) adapter sequence (leave empty if not to trim):",""))}
        output$minLength <- renderUI(numericInput("minLength", "Minimal allowed sequence length:",15))
        output$nBases <- renderUI(numericInput("nBases", "Maximal number of Ns allowed per sequence:",5))
        if (input$second_choice == "1") {output$fastq_list_title <- renderUI(h3("Sample list (single-end reads):"))} else {output$fastq_list_title <- renderUI(h3("Sample list (paired-end reads):"))}
        output$fastq_list <- renderTable(fastq_list)
        rv$fastq_list <- fastq_list
        output$file_input <- renderUI(NULL)
      }
      
       
      
    }
    if (input$first_choice == "4") {
      updateTabsetPanel(session, "inTabset",selected = "heatmaps")
      if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
      inFile <- input$file_input
      sample_list_path <- inFile$datapath
      sample_list <- read.delim(sample_list_path, header = F)
      if (ncol(sample_list)!=3) {
        output$warning1_h <- renderText("Error loading sample list. Tab-delimited text file must have 3 columns.")
      } else {
        if (sum(sample_list=="")>0) {
          output$warning1_h <- renderText("Error loading sample list. Missing values. All fields in the tab-delimited text file describing samples must be non-empty.")
        } else {
          cov_files <- as.character(sample_list[,2])  
          for (i in 1:length(cov_files)) {
            temp <- length(list.files(pattern = cov_files[i] ))
            if (i==1) {
              proceed <- temp
            } else {
              proceed <- temp * proceed
            }     
          }
          if (proceed==0) {
            output$warning1_h <- renderText("Error loading sample list. At least one file not found in the input directory. Make sure to provide the correct file names.")
          } else {
            output$text6 <- renderUI(NULL); output$bedgraph_tab_example <- renderUI(NULL)
            colnames(sample_list) <- c("SampleID","File","Color")
            sample_list$File <- as.character(sample_list$File)
            sample_list$SampleID <- as.character(sample_list$SampleID)
            output$sample_list_heatmaps <- renderTable(sample_list)
            rv$sample_list <- sample_list
            output$heatmaps_button <- renderUI(actionButton("heatmaps_button", "Get heatmaps and metaplot"))
          }
        }
      }
    }
  })
 
  # mapping
  observeEvent(input$preprocessing_button, {
    if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
    output$text5 <- renderUI(NULL)
    fastq_list <- rv$fastq_list
    truncateStartBases <- input$truncateStartBases 
    truncateEndBases <- input$truncateEndBases 
    minLength <- input$minLength
    nBases <- input$nBases
    Lpattern <- "" 
    Rpattern <- ""
    if (input$second_choice==1) {
      filter_results <- preprocessing(fastq_list,1,truncateStartBases,truncateEndBases,input$Lpattern,input$Rpattern ,minLength,nBases)
      fastq_list <- within(fastq_list, totalSequences <- filter_results[1,])
      fastq_list <- within(fastq_list, totalPassed <- filter_results[nrow(filter_results),])
      fastq_list$totalSequences <- as.integer(as.character(fastq_list$totalSequences)); fastq_list$totalPassed <- as.integer(as.character(fastq_list$totalPassed))
      fastq_list <- within(fastq_list, FileName <- paste("filtered_",FileName,sep = ""))
    }
    if (input$second_choice==2) {
      filter_results <- preprocessing(fastq_list,2,truncateStartBases,truncateEndBases,Lpattern,Rpattern,minLength,nBases)
      fastq_list <- within(fastq_list, totalSequences <- filter_results[1,])
      fastq_list <- within(fastq_list, totalPassed <- filter_results[nrow(filter_results),])
      fastq_list$totalSequences <- as.integer(as.character(fastq_list$totalSequences)); fastq_list$totalPassed <- as.integer(as.character(fastq_list$totalPassed))
      fastq_list <- within(fastq_list, FileName1 <- paste("filtered_",FileName1,sep = ""))
      fastq_list <- within(fastq_list, FileName2 <- paste("filtered_",FileName2,sep = ""))
      output$qalign_pair <- renderUI(radioButtons("qalign_pair", label = "Type of paired-end library:", choices = list("Paired mates are in forward/reverse position" = 1, "Paired mates are in forward/forward position" = 2, "Paired mates are in reverse/forward position" = 3), selected = 1))
    }
    fastq_list_2 <- fastq_list
    output$fastq_list <- renderTable(fastq_list_2)
    if (input$second_choice==1) {
      fastq_list <- fastq_list[,1:2]
      write.table(fastq_list, file = "ChIPdig_sample_list_single_end_reads.txt", sep = "\t", row.names = F, quote = F, na="", col.names = T)
    }
    if (input$second_choice==2) {
      fastq_list <- fastq_list[,1:3]
      write.table(fastq_list, file = "ChIPdig_sample_list_paired_end_reads.txt", sep = "\t", row.names = F, quote = F, na="", col.names = T)
    } 
    rv$fastq_list <- fastq_list
    output$preprocessing_title <- renderUI(NULL)
    output$preprocessing_no_button <- renderUI(NULL)
    output$preprocessing_button <- renderUI(NULL)
    output$truncateStartBases <- renderUI(NULL)
    output$truncateEndBases <- renderUI(NULL)
    output$minLength <- renderUI(NULL)
    output$nBases <- renderUI(NULL)
    output$Lpattern <- renderUI(NULL)
    output$Rpattern <- renderUI(NULL)
    setwd(initialdir)
    output$mapping_title <- renderUI(h4("Alignment"))
    output$mapping_button <- renderUI(actionButton("mapping_button", "Align reads to reference genome"))
    output$genome <- renderUI(selectInput("genome", "Choose reference genome:", available.genomes(), multiple=F))
  })
  observeEvent(input$preprocessing_no_button, {
    output$text5 <- renderUI(NULL)
    output$preprocessing_title <- renderUI(NULL)
    output$preprocessing_no_button <- renderUI(NULL)
    output$preprocessing_button <- renderUI(NULL)
    output$truncateStartBases <- renderUI(NULL)
    output$truncateEndBases <- renderUI(NULL)
    output$minLength <- renderUI(NULL)
    output$nBases <- renderUI(NULL)
    output$Lpattern <- renderUI(NULL)
    output$Rpattern <- renderUI(NULL)
    if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
    if (input$second_choice==1) {
      fastq_list <- rv$fastq_list
      write.table(fastq_list, file = "ChIPdig_sample_list_single_end_reads.txt", sep = "\t", row.names = F, quote = F, na="", col.names = T)
    }
    if (input$second_choice==2) {
      fastq_list <- rv$fastq_list
      write.table(fastq_list, file = "ChIPdig_sample_list_paired_end_reads.txt", sep = "\t", row.names = F, quote = F, na="", col.names = T)
      output$qalign_pair <- renderUI(radioButtons("qalign_pair", label = "Type of paired-end library:", choices = list("Paired mates are in forward/reverse position" = 1, "Paired mates are in forward/forward position" = 2, "Paired mates are in reverse/forward position" = 3), selected = 1))
    } 
    setwd(initialdir)
    output$mapping_title <- renderUI(h4("Alignment"))
    output$mapping_button <- renderUI(actionButton("mapping_button", "Align reads to reference genome"))
    output$genome <- renderUI(selectInput("genome", "Choose reference genome:", available.genomes(), multiple=F))
  })
  observeEvent(input$mapping_button, {
    if (input$second_choice==1) {
      if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
      if (is.null(input$dir)) {
        path <- file.path(paste(initialdir,"/ChIPdig_sample_list_single_end_reads.txt", sep = ""))
      } else {
        path <- file.path(paste(parseDirPath(volumes, directory()),"/ChIPdig_sample_list_single_end_reads.txt", sep = ""))
      }
      mapping_project <- qAlign(path, genome=input$genome) 
    }
    if (input$second_choice==2) {
      if (input$qalign_pair==1) {
        paired="fr"
      } else {
        if (input$qalign_pair==2) {
          paired="ff"
        } else {
          if (input$qalign_pair==3) {
            paired="rf"
          }
        }
      } 
      if (is.null(input$dir)) {
        path <- file.path(paste(initialdir,"/ChIPdig_sample_list_paired_end_reads.txt", sep = ""))
      } else {
        path <- file.path(paste(parseDirPath(volumes, directory()),"/ChIPdig_sample_list_paired_end_reads.txt", sep = ""))
      }
      mapping_project <- qAlign(path, genome=input$genome, paired=paired) 
    }
    qQCReport(mapping_project, pdfFilename="qc_report.pdf")
    output$mapping_success <- renderUI(h3("Genomic alignments have been created successfully! A quality control report (PDF file) in the form of a series of diagnostic plots with details on sequences and alignments was created in the input directory."))
    setwd(initialdir)
  })

  # annotation
  observeEvent(input$anno_button, {
    updateTabsetPanel(session, "inTabset",selected = "anno")
    inFile <- input$file_input
    path <- inFile$datapath
    name <- inFile$name
    txdb <- makeTxDbFromUCSC(genome=input$genome, tablename="ensGene")
    peakAnno <- annotatePeak(path, tssRegion=input$TSS_def,TxDb=txdb)
    annotated_table <- as.data.frame(peakAnno)
    output$annotated_table <- renderTable(annotated_table)
    output$annotated_table_download <- downloadHandler(
      filename <- function() {
        paste(name, "_annotation_table.csv", sep = "") 
      },
      content <- function(file) {
        write.csv(annotated_table, file, row.names = F)
      }
    ) 
    output$annotated_pie <- renderPlot(plotAnnoPie(peakAnno))
    output$annotated_pie_download <- downloadHandler(
      filename <- function() {
        paste(name, "_annotation_pie.tiff", sep = "") 
      },
      content <- function(file2) {
        tiff(file2, width = 20, height = 14, units = "cm", res = 300)
        print(plotAnnoPie(peakAnno))
        dev.off()
      }
    )    
  })
  
  # processing of mapped reads - initial processing
  observeEvent(input$read_extension_single, {
    if (input$read_extension_single == 2) {
      output$read_extension_numeric <- renderUI(numericInput("read_extension_numeric", "Extend reads to this length:", 250, min = 0)) 
    } else {
      output$read_extension_numeric <- renderUI(NULL)
    }
  })
  observeEvent(input$process_reads, {
    updateTabsetPanel(session, "inTabset",selected = "initial_proc")
    if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
    sample_list <- rv$sample_list 
    single_or_paired <- rv$single_or_paired
    bam.files <- bam.files_list(sample_list,sample_list,1,1)
    ext<-0
    rv$ext <- "NO"
    sample_list_with_parameters <- sample_list
    sample_list_with_parameters <- within(sample_list_with_parameters, ext_length_treat <- 0)
    sample_list_with_parameters <- within(sample_list_with_parameters, ext_length_control <- 0)
    if (as.character(single_or_paired) == "single") {
      if(is.null(input$read_extension_single)) {} else {
        if (input$read_extension_single==1) {
          sizes <- fragment_size_estimation_batch(bam.files)
          ext <- list(sizes,sizes)
          sample_list_with_parameters <- add_vector_elements_to_last_2_columns_of_frame(sample_list_with_parameters,sizes)
          output$read_extension_decision <- renderUI(NULL)
          rv$ext <- "YES"
        } else {
          ext <- input$read_extension_numeric
          sample_list_with_parameters <- within(sample_list_with_parameters, ext_length_treat <- ext)
          sample_list_with_parameters <- within(sample_list_with_parameters, ext_length_control <- ext)
          output$read_extension_decision <- renderUI(NULL)
          rv$ext <- "YES"
        }
      }
    }
    output$read_extension_single <- renderUI(NULL)
    output$read_extension_numeric <- renderUI(NULL)
    param <- readParam()
    param <- reform(param, minq=20)
    if (input$dupl_reads_remove==TRUE) {
      param <- reform(param, dedup=TRUE)
      } 
    rv$dupl_reads_remove <- input$dupl_reads_remove
    output$dupl_reads_remove <- renderUI(NULL)
    if (input$second_choice=="2") {param <- reform(param, max.frag=input$max_distance_pairs, pe="both")}
    output$max_distance_pairs <- renderUI(NULL)
    width <- input$bin_size
    rv$bin_size <- input$bin_size
    output$bin_size <- renderUI(NULL)
    counts <- windowCounts(bam.files, width=width, spacing=width, ext=ext, shift=0, filter=0, bin=TRUE, param=param)
    edger <- counts_to_edger(counts,sample_list)
    TMM <- edger_to_TMM(edger)
    coverages <- coverages_function(sample_list,TMM)
    rv$counts <- counts
    rv$edger <- edger
    rv$TMM <- TMM
    rv$coverages <- coverages
    output$process_reads <- renderUI(NULL)
    output$process_reads_note <- renderUI(NULL)
    libsizes <- counts$totals
    sample_list_with_parameters <- within(sample_list_with_parameters, libsize_treat <- 0)
    sample_list_with_parameters <- within(sample_list_with_parameters, libsize_control <- 0)
    sample_list_with_parameters <- add_vector_elements_to_last_2_columns_of_frame(sample_list_with_parameters,libsizes)
    rv$sample_list_with_parameters <- sample_list_with_parameters
    sample_list_with_parameters <- sample_list_with_parameters[,c(1,ncol(sample_list_with_parameters)-3,ncol(sample_list_with_parameters)-2,ncol(sample_list_with_parameters)-1,ncol(sample_list_with_parameters)-0)]
    sample_list_with_parameters <- within(sample_list_with_parameters, ext_length_treat <- as.integer(ext_length_treat))
    sample_list_with_parameters <- within(sample_list_with_parameters, ext_length_control <- as.integer(ext_length_control))
    sample_list_with_parameters <- within(sample_list_with_parameters, libsize_treat <- as.integer(libsize_treat))
    sample_list_with_parameters <- within(sample_list_with_parameters, libsize_control <- as.integer(libsize_control))
    colnames(sample_list_with_parameters) <- c("Sample ID","Treatment reads extended to:","Control reads extended to:","Library size, treatment:","Library size, control:")
    output$sample_list_with_parameters_title <- renderUI(h3("Library sizes and read extension sizes:"))
    output$sample_list_with_parameters <- renderTable(sample_list_with_parameters)
    output$peak_calling_title <- renderUI(h4("Peak calling"))
    output$peak_calling <- renderUI(actionButton("peak_calling", "Call peaks"))
    output$peak_calling_note <- renderUI("NOTE: Peak files will be created in a subfolder created in the input directory. Peak calling may take a while. Be patient.")
    output$user_suppl_peaks_button <- renderUI(actionButton("user_suppl_peaks_button", "Load peaks called externally"))
    output$PP <- renderUI(numericInput("PP", "Posterior probability threshold:", 0.5, min = 0)) 
    if(length(unique(sample_list$Condition)) > 1 ) {
      output$diff_title <- renderUI(h4("Differential enrichment analysis"))
      if(replicates_present_assessment(sample_list)==FALSE) {output$diff_note <- renderUI("NOTE: No replicates were detected in the sample list provided. Negative binomial distribution coefficient set to 0.05 (a typical value based on past experience).")}
      if(replicates_present_assessment(sample_list)==FALSE) {rv$repl_present <- FALSE} else {rv$repl_present <- TRUE}
      output$diff_action <- renderUI(actionButton("diff_action", "Perform differential enrichment analysis on chosen contrast."))
      output$FDR <- renderUI(numericInput("FDR", "False discovery rate threshold (for MA plot):", 0.1, min = 0))
      output$contrast_test <- renderUI(selectInput("contrast_test", "Choose condition to test:", unique(sample_list$Condition), multiple=F))
      output$contrast_ref <- renderUI(selectInput("contrast_ref", "Choose reference condition:", unique(sample_list$Condition), multiple=F))
      output$bed_input_for_diff <- renderUI(fileInput('bed_input_for_diff', 'Select regions to limit differential enrichment analysis to (BED format) (please do not include column names):',multiple = FALSE, accept='text/comma-separated-values,text/plain'))
      output$FC_cutoff_all_samples <- renderUI(numericInput("FC_cutoff_all_samples", "Discard bins in which fold enrichment (treatment/control) in all samples is less than (bp):", 1, min = 0))
    }
    strandMode <- 1
    if (as.character(single_or_paired) == "paired") {strandMode <- input$strandMode}
    rv$strandMode <- strandMode
    read_extension_option <- rv$ext
    dupl_reads_remove <- rv$dupl_reads_remove
    chrom_sizes <- chrom_sizes_fun(sample_list,single_or_paired,strandMode)
    output$chrom_sizes_title <- renderUI(h3("Chromosome sizes:"))
    output$chrom_sizes <- renderTable(chrom_sizes)
    chrom_sizes <- chrom_sizes[,c(1,3)]
    colnames(chrom_sizes) <- c("Chromosome","Size")
    rv$chrom_sizes <- chrom_sizes
    setwd(initialdir)
    output$generate_coverage_title <- renderUI(h4("Coverage files export"))
    output$generate_coverage <- renderUI(actionButton("generate_coverage", "Export coverage files"))
    output$generate_coverage_note <- renderUI("NOTE: Coverage files (treatment, control and normalized) will be created in a subfolder created in the input directory. Be patient.")
    output$MDS_all_title <- renderUI(h3("Multidimensional scaling plot (genome-wide):"))
    output$MDS_all <- renderPlot(MDS_plot(TMM,sample_list))
    output$MDS_all_download <- downloadHandler(
      filename <- function() {
        "MDS_plot.tiff" 
      },
      content <- function(file) {
        tiff(file, width = 20, height = 14, units = "cm", res = 300)
        print(MDS_plot(TMM,sample_list))
        dev.off()
      }
    )
  })
  
  # peak calling and peak operations
  observeEvent(input$peak_calling, {
    updateTabsetPanel(session, "inTabset",selected = "peak_calling")
    if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
    chrom_sizes <- rv$chrom_sizes
    sample_list <- rv$sample_list
    sample_list_with_fr_sizes <- rv$sample_list_with_parameters
    single_or_paired <- rv$single_or_paired
    dupl_reads_remove <- rv$dupl_reads_remove
    strandMode <- rv$strandMode
    read_extension_option <- rv$ext
    bin_size <- rv$bin_size
    PP <- input$PP
    table_with_peaks <- peaks_batch(sample_list,sample_list_with_fr_sizes,single_or_paired,dupl_reads_remove,strandMode,read_extension_option,PP,bin_size)
    table_with_peaks <- within(table_with_peaks, Peak_count <- as.integer(as.character(Peak_count)))
    colnames(table_with_peaks)[ncol(table_with_peaks)] <- "Peak count"
    table_with_peaks <- table_with_peaks[,c(1,ncol(table_with_peaks))]
    output$table_with_peaks_title <- renderUI(h3("Number of peaks called for each treatment-control pair:"))
    output$table_with_peaks <- renderTable(table_with_peaks)
    #setwd(initialdir)
    output$add_track_line <- renderUI(actionButton("add_track_line", "Add UCSC track to peak files"))
    output$add_track_line_note <- renderUI("NOTE: The color defined in the sample list will be used.")
    output$repl_peaks <- renderUI(actionButton("repl_peaks", "Generate file with replicated peaks in each condition"))
    output$repl_peaks_note <- renderUI("NOTE: A file with replicated peaks for each condition will be created in a subfolder created in peak files folder.")
  })
  observeEvent(input$add_track_line, {
    sample_list <- rv$sample_list
    add_track_using_sample_list_colors(sample_list)
  })
  observeEvent(input$repl_peaks, {
    sample_list <- rv$sample_list
    replicated_peaks_batch_app_peaks(sample_list)
    output$consensus_peaks <- renderUI(actionButton("consensus_peaks", "Generate consensus peaks"))
    output$consensus_peaks_note <- renderUI("NOTE: A consensus peak set (peaks present in all files) will be created in the peaks directory. Peaks smaller than the input bin size will be discarded.")
  })
  observeEvent(input$consensus_peaks,{
    setwd("./peaks/")
    setwd("./replicated_peaks/")
    files_for_consensus <- list.files(pattern = "\\_replicated_peaks.bed$")
    bin_size <- rv$bin_size
    consensus_peak_set_function(files_for_consensus,bin_size)
    setwd(initialdir)
  })
  observeEvent(input$user_suppl_peaks_button, {
    output$peaks_input <- renderUI(fileInput('peaks_input', 'Select tab-delimited text file listing peak files with matching IDs and conditions (see instructions):',multiple = FALSE, accept='text/comma-separated-values,text/plain'))
    output$ext_peaks_note <- renderUI("NOTE: A file with replicated peaks for each condition will be created in a subfolder created in the input directory. A consensus peak set (peaks present in all files) will be created in the input directory. Peaks smaller than the input bin size will be discarded.")
  })
  observeEvent(input$peaks_input, {
    updateTabsetPanel(session, "inTabset",selected = "peak_calling")
    if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
    inFile <- input$peaks_input
    ext_peaks_path <- inFile$datapath
    ext_peaks_table <- ext_peaks_import(ext_peaks_path)
    output$ext_peaks_table_title <- renderUI(h3("External peak files:"))
    output$ext_peaks_table <- renderTable(ext_peaks_table)
    replicated_peaks_batch(ext_peaks_table)
    setwd("./replicated_peaks/")
    files_for_consensus <- list.files(pattern = "\\_replicated_peaks.bed$")
    bin_size <- rv$bin_size
    consensus_peak_set_function(files_for_consensus,bin_size)
    setwd("..")
    setwd(initialdir)
  })
  
  # diff enrichment analysis
  observeEvent(input$diff_action, {
    updateTabsetPanel(session, "inTabset",selected = "diff")
    sample_list <- rv$sample_list
    counts <- rv$counts
    edger <- rv$edger
    TMM <- rv$TMM
    coverages <- rv$coverages
    repl_present <- rv$repl_present
    counts_backg_removed <- counts_input_filter(coverages,input$FC_cutoff_all_samples,counts)
    if (is.null(input$bed_input_for_diff)) {
      counts_subset <- counts_backg_removed
      peaks <- rowRanges(counts_subset)
    } else {
      inFile <- input$bed_input_for_diff
      peaks_path <- inFile$datapath
      peaks <- peaks_import(peaks_path)
      counts_subset <- subsetByOverlaps(counts_backg_removed, peaks)
    }
    edger_subset <- edger_subset_keeping_params_all_samples(TMM,counts_subset,sample_list)
    design <- data.frame(rownames(edger_subset$samples),edger_subset$samples$group)
    colnames(design) <- c("Samples","Condition")
    levels(design$Condition)
    design$Condition <- factor(design$Condition, levels = unique(as.character(edger_subset$samples$group)))
    colnames <- levels(design$Condition)
    design <- model.matrix(~0+Condition, data=design)
    colnames(design) <- make.names(colnames)
    contrast_test <- input$contrast_test
    contrast_ref <- input$contrast_ref
    str1 <- make.names(paste(contrast_test,"_treatment",sep=""))
    str2 <- make.names(paste(contrast_ref,"_treatment",sep=""))
    contrast_string <- paste(str1,str2,sep="-")
    assign("contrast_string",contrast_string,envir = .GlobalEnv)
    assign("design",design,envir = .GlobalEnv)
    contrast <- makeContrasts(c(contrast_string), levels=design)
    if (repl_present==TRUE) {
      y <- estimateDisp(edger_subset, design)
      fit <- glmQLFit(y, design, robust=TRUE)
      result <- glmQLFTest(fit, contrast = contrast)
    } else {
      fit <- glmFit(edger_subset, design, dispersion=0.05)
      result <- glmLRT(fit, contrast=contrast)
    }
    olap <- findOverlaps(peaks, rowRanges(counts_subset))
    tabbroad <- combineOverlaps(olap, result$table)
    mcols(peaks) <- tabbroad
    coverages_for_plot <- coverages_function(sample_list,edger_subset)
    binranges <- binranges_fun(coverages_for_plot, sample_list, contrast_test,contrast_ref, counts_subset)
    results <- diff_enrich_tab_fun(peaks,binranges,contrast_test,contrast_ref)
    output$results_plot <- renderPlot(diff_result_plot(results,input$FDR,contrast_test,contrast_ref))  
    output$results_plot_download <- downloadHandler(
      filename <- function() {
        paste ("MA_plot_",contrast_test,"_vs_",contrast_ref,".tiff", sep = "") 
      },
      content <- function(file) {
        tiff(file, width = 20, height = 14, units = "cm", res = 300)
        print(diff_result_plot(results,input$FDR,contrast_test,contrast_ref))
        dev.off()
      }
    )    
    output$boxplot_Q95 <- renderPlot(boxplot_Q95(results,contrast_test,contrast_ref))
    output$boxplot_Q95_download <- downloadHandler(
      filename <- function() {
        paste ("boxplot_",contrast_test,"_vs_",contrast_ref,".tiff", sep = "") 
      },
      content <- function(file) {
        tiff(file, width = 12, height = 10, units = "cm", res = 300)
        print(boxplot_Q95(results,contrast_test,contrast_ref))
        dev.off()
      }
    )
    output$results_table_download <- downloadHandler(
      filename <- function() {
        paste ("diff_enrich_",contrast_test,"_vs_",contrast_ref,".csv", sep = "") 
      },
      content <- function(file) {
        write.csv(results, file, row.names = F)
      }
    ) 
    output$results <- renderTable(results)
  })
  
  # export of coverage files
  observeEvent(input$generate_coverage,{
    if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
    counts <- rv$counts
    edger <- rv$edger
    sample_list <- rv$sample_list
    export_bedgraphs(counts,edger,sample_list)
    setwd(initialdir)
    assign("counts",counts,envir = .GlobalEnv)
    assign("edger",edger,envir = .GlobalEnv)
    assign("sample_list",sample_list,envir = .GlobalEnv)
  })
  
  # heatmaps and metaplots
  observeEvent(input$heatmaps_mode, {
    if (input$heatmaps_mode=="both") {
      output$L <- renderUI(numericInput("L", "Body size (bp):", 500, min = 100))
    } else {
      output$L <- renderUI(NULL)
    }
    output$width_heatmaps <- renderUI(numericInput("width_heatmaps", "Width (cm) of downloaded plot:", 30))
    output$height_heatmaps <- renderUI(numericInput("height_heatmaps", "Height (cm) of downloaded plot:", 15))
    output$width_metaplot <- renderUI(numericInput("width_metaplot", "Width (cm) of downloaded plot:", 20))
    output$height_metaplot <- renderUI(numericInput("height_metaplot", "Height (cm) of downloaded plot:", 15))
  })
  observeEvent(input$heatmaps_button, {
    if (is.null(input$dir)) {} else {setwd(parseDirPath(volumes, directory()))}
	inFile <- input$heatmaps_regions
    heatmaps_regions_path <- inFile$datapath
    sample_list <- rv$sample_list
    regions_string <- heatmaps_regions_path
    binsize <- input$heatmaps_binsize
    length_before <- input$length_before
    length_after <- input$length_after  
    length_all_before <- input$length_before
    length_all_after <- input$length_after  
    if (input$heatmaps_mode=="start") {
      mode <- "TSS"
    } else {
      if (input$heatmaps_mode=="end") {
        mode <- "TTS"
      } else {
        if (input$heatmaps_mode=="both") {
          L <- input$L
          if (length_before==0&length_after==0) {
            mode <- "body"
          } else {
            mode <- "body_plus"
          }
        }
      }
    }
    coverages <- coverage_import_batch(sample_list)
    chrom_sizes <- chrom_sizes_from_coverage_batch(coverages)
    matrices <- matrix_fun(sample_list,regions_string, mode,binsize, length_before,length_after, L, length_all_before,length_all_after,chrom_sizes)
    output$heatmaps_heatmaps <- renderPlot(multiplot_enriched_heatmaps(sample_list, matrices, mode, binsize, length_before,length_after, L,  length_all_before,length_all_after))
    output$heatmaps_metaplot <- renderPlot(metaplot(sample_list, matrices, mode, binsize, length_before,length_after, L,  length_all_before,length_all_after))
    
    output$heatmaps_heatmaps_download <- downloadHandler(
      filename <- function() {
        "heatmaps.tiff" 
      },
      content <- function(file) {
        tiff(file, width = input$width_heatmaps, height = input$height_heatmaps, units = "cm", res = 300)
        print(multiplot_enriched_heatmaps(sample_list, matrices, mode, binsize, length_before,length_after, L,  length_all_before,length_all_after))
        dev.off()
      }
    )
    
    output$heatmaps_metaplot_download <- downloadHandler(
      filename <- function() {
        "metaplot.tiff" 
      },
      content <- function(file) {
        tiff(file, width = input$width_metaplot, height = input$height_metaplot, units = "cm", res = 300)
        print(metaplot(sample_list, matrices, mode, binsize, length_before,length_after, L,  length_all_before,length_all_after))
        dev.off()
      }
    )
    
  })
  
})
