shinyUI(pageWithSidebar(
  
  
  headerPanel("ChIPdig"),
  
  sidebarPanel(
    
    
    
    h2("Select analysis module"),   
    radioButtons("first_choice", label = NULL, choices = list("Align reads to the genome" = 1, "Normalize and compare ChIP-seq data sets" = 2, "Annotate genomic regions" = 3, "Visualize coverage in specific regions using background-normalized coverage files" = 4), selected = FALSE),
    uiOutput("second_choice"),
    
    h3("Select input folder and files"),
    shinyDirButton("dir", "Choose folder with files to be analyzed", "Upload"), br(), br(),
    verbatimTextOutput("dir"),
    
    uiOutput("anno_title"),
    uiOutput("BAM_sample_processing_title"),
    uiOutput("heatmaps_title"),
    uiOutput("align_title"),
    uiOutput("file_input"),
    uiOutput("text7"),
    uiOutput("preprocessing_title"),
    uiOutput("truncateStartBases"),
    uiOutput("truncateEndBases"),
    uiOutput("Lpattern"),
    uiOutput("Rpattern"),
    uiOutput("minLength"),
    uiOutput("nBases"),
    uiOutput("preprocessing_button"), br(),
    uiOutput("preprocessing_no_button"),
    
    uiOutput("heatmaps_regions"),
    uiOutput("heatmaps_mode"),
    uiOutput("heatmaps_binsize"),
    uiOutput("length_before"),
    uiOutput("length_after"),
    uiOutput("L"),
    uiOutput("heatmaps_button"),
    
    uiOutput("mapping_title"),
    uiOutput("genome"),
    uiOutput("qalign_pair"),
    uiOutput("mapping_button"),
    
    uiOutput("TSS_def"),
    uiOutput("anno_button"),
    
    uiOutput("bin_size"),
    uiOutput("strandMode"),
    uiOutput("dupl_reads_remove"),
    uiOutput("read_extension_single"),
    uiOutput("read_extension_numeric"),
    uiOutput("read_extension_decision"),
    uiOutput("max_distance_pairs"),
    uiOutput("process_reads"),
    uiOutput("process_reads_note"),
    
    
    uiOutput("generate_coverage_title"),
    uiOutput("generate_coverage"),
    uiOutput("generate_coverage_note"),
    uiOutput("peak_calling_title"),
    uiOutput("PP"),
    uiOutput("peak_calling"),
    uiOutput("peak_calling_note"), br(),
    uiOutput("user_suppl_peaks_button"), br(),
    
    uiOutput("add_track_line"),
    uiOutput("add_track_line_note"), 
    uiOutput("repl_peaks"),
    uiOutput("repl_peaks_note"),
    uiOutput("consensus_peaks"),
    uiOutput("consensus_peaks_note"),
    

    
    
    uiOutput("peaks_input"),
    uiOutput("diff_title"),
    uiOutput("FDR"),
    uiOutput("FC_cutoff_all_samples"),
    uiOutput("bed_input_for_diff"),
    uiOutput("contrast_test"),
    uiOutput("contrast_ref"),
    uiOutput("diff_action"),
    uiOutput("diff_note")
    
  ),
  
  
  mainPanel(
    tabsetPanel(id = "inTabset",
      tabPanel(title = "Mapping", value = "map", span(h3(textOutput("warning_mapping")) , style="color:red"),uiOutput("text1"),  uiOutput('text2'), tableOutput('singleread_exp_input_example'), uiOutput('text3'),    tableOutput('pairedend_exp_input_example'),uiOutput('text4'),uiOutput("fastq_list_title"), tableOutput('fastq_list'),uiOutput('text5'),uiOutput("mapping_success")),
      tabPanel(title = "Initial processing of mapped reads", value = "initial_proc", uiOutput('text8'), tableOutput('data_frame_processing_example'),uiOutput("sample_list_title"),   span(h3(textOutput("warning1")) , style="color:red"),uiOutput('warning2'), tableOutput('warning3') ,   tableOutput('sample_list'),uiOutput("sample_list_with_parameters_title"),tableOutput('sample_list_with_parameters'),uiOutput("chrom_sizes_title"),tableOutput('chrom_sizes'),uiOutput("MDS_all_title"),plotOutput('MDS_all', width = "65%"),downloadButton('MDS_all_download', 'Download plot')),
      tabPanel(title = "Peak calling", value = "peak_calling", uiOutput("table_with_peaks_title"),tableOutput('table_with_peaks'),uiOutput("ext_peaks_table_title"),tableOutput('ext_peaks_table')),
      tabPanel(title = "Differential enrichment analysis", value = "diff",h3("Mean-difference plot:"),plotOutput('results_plot', width = "65%"),downloadButton('results_plot_download', 'Download plot'),h3("Box-and-whisker plot:"),   plotOutput('boxplot_Q95', width = "50%"),downloadButton('boxplot_Q95_download', 'Download plot')    ,h3("Results table:"),downloadButton('results_table_download', 'Download table'),tableOutput('results')),
      tabPanel(title = "Peak annotation", value = "anno", plotOutput('annotated_pie', width = "50%"),downloadButton('annotated_pie_download', 'Download plot'), tableOutput('annotated_table'),downloadButton('annotated_table_download', 'Download table')),
      tabPanel(title = "Heatmaps and metaplots", value = "heatmaps", span(h3(textOutput("warning1_h")) , style="color:red") ,uiOutput('text6'),tableOutput('bedgraph_tab_example'), tableOutput('sample_list_heatmaps'), h3("Heatmaps:"),plotOutput('heatmaps_heatmaps'),uiOutput("width_heatmaps"),uiOutput("height_heatmaps"),downloadButton('heatmaps_heatmaps_download', 'Download plot'),h3("Metaplot:"),plotOutput('heatmaps_metaplot',width = "65%"),uiOutput("width_metaplot"),uiOutput("height_metaplot"),downloadButton('heatmaps_metaplot_download', 'Download plot'))
   ))    
  
))

