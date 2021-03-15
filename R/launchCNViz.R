#' Launches CNViz, a shiny app to visualize your sample's copy number data.
#'
#' CNViz launches a shiny application to visualize your sample's copy number data.
#' At least one of probe_data, gene_data, or segment_data must be supplied;
#' sample_name, snv_data and meta_data are all optional.
#'
#' @param sample_name A string with the ID/name of your sample.
#' @param probe_data A dataframe containing probe-level data. Column names must include chr, gene, start, end, log2. Optional column: weight. chr column should be formatted as 'chr1' through 'chrX', 'chrY'. start, end and log2 should be numeric.
#' @param gene_data A dataframe containing gene-level data - one row per gene. Column names must include chr, gene, start, end, log2. Optional columns: weight, loh; where loh values are TRUE or FALSE. chr column should be formatted as 'chr1' through 'chrX', 'chrY'. start, end and log2 should be numeric.
#' @param segment_data A dataframe containing segment-level data. Column names must include chr, start, end, log2. Optional column: loh; where loh values are TRUE or FALSE. chr column should be formatted as 'chr1' through 'chrX', 'chrY'. start, end and log2 should be numeric.
#' @param snv_data A dataframe containg SNVs and columns of your choosing. The only required column is gene. Optional columns: start, mutation_id; where start indicates the starting position of the mutation and mutation_id is a string in any format. Additional columns might include depth, allelic_fraction.
#' @param meta_data A dataframe containing your sample's metadata - columns of your choosing. Optional column: ploidy; ploidy will be rounded to the nearest whole number. Additional columns might include purity. This dataframe should only have one row.
#'
#' @return a Shiny application
#'
#' @import shiny
#' @import utils
#' @importFrom grDevices dev.off pdf
#' @importFrom dplyr select filter summarise mutate left_join group_by n between
#' @import stats
#' @importFrom plotly plot_ly add_segments add_trace layout renderPlotly plotlyOutput event_data subplot
#' @importFrom karyoploteR plotKaryotype
#' @importFrom CopyNumberPlots plotCopyNumberCalls
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom magrittr %>%
#' @importFrom cBioPortalData cBioPortal getDataByGenePanel
#' @importFrom DT DTOutput renderDT formatPercentage datatable formatStyle
#' @importFrom scales rescale
#' @importFrom graphics legend
#'
#' @export launchCNViz
#'

launchCNViz <- function(sample_name = "sample", probe_data = data.frame(), gene_data = data.frame(), segment_data = data.frame(), snv_data = data.frame(), meta_data = data.frame()) {

  cbio <- cBioPortal()

  colnames(probe_data) <- tolower(colnames(probe_data))
  colnames(gene_data) <- tolower(colnames(gene_data))
  colnames(segment_data) <- tolower(colnames(segment_data))
  colnames(snv_data) <- tolower(colnames(snv_data))
  colnames(meta_data) <- tolower(colnames(meta_data))

  chromosomes <- c(paste0("chr", "1":"22"), "chrX", "chrY")
  if(nrow(gene_data) > 0){
    genes <- c("", sort(unique(gene_data$gene)))
  } else if(nrow(probe_data) > 0){
    genes <- c("", sort(unique(probe_data$gene)))
  } else{
    genes <- c("No gene information")
  }

  shinyApp(ui <- fluidPage(
    navbarPage("CNViz",
               tabPanel("Patient Data", fluid=TRUE,
                        sidebarLayout(
                          sidebarPanel(
                            width = 3,
                            selectInput(inputId = "chr",
                                        label = "chromosome",
                                        choices = c("all", chromosomes),
                                        selected = "all"),
                            selectizeInput(inputId = "gene",
                                           label = "gene",
                                           choices = genes,
                                           selected = ""),
                            conditionalPanel(condition = 'output.karyotype != null', downloadButton("karyotype", "karyotype")),
                            br(), br(),
                            if(nrow(gene_data) > 0 | nrow(probe_data) > 0){
                              img(src="https://cnviz.s3.amazonaws.com/markers.png", width = "100%")
                            },
                            if(nrow(snv_data) > 0 & (nrow(gene_data) + nrow(probe_data) > 0)){
                              img(src="https://cnviz.s3.amazonaws.com/mutation.png", width = "100%")
                            },
                            if("loh" %in% colnames(gene_data)){
                              img(src="https://cnviz.s3.amazonaws.com/marker_loh.png", width = "100%")
                            },
                            if("loh" %in% colnames(gene_data) & nrow(snv_data) > 0){
                              img(src="https://cnviz.s3.amazonaws.com/mutation_loh.png", width = "100%")
                            },
                            if(nrow(segment_data) > 0){
                              img(src="https://cnviz.s3.amazonaws.com/segment.png", width = "100%")
                            },
                            if("loh" %in% colnames(segment_data)){
                              img(src="https://cnviz.s3.amazonaws.com/segment_loh.png", width = "100%")
                            }
                          ),
                          mainPanel(
                            h3(sample_name),
                            tableOutput("meta"),
                            br(),
                            plotlyOutput("chr_plot", width = "100%"),
                            conditionalPanel("input.chr != 'all'",
                                             column(10, offset = 1,
                                                    h4(textOutput("gene_text")), br(),
                                                    DTOutput("mutations"), br(),
                                                    plotlyOutput("selected_plot"), br(),
                                                    style = 'padding:20px'))
                          )
                        )),
               tabPanel("TCGA Pan-Cancer Atlas 2018 Data", fluid=TRUE,
                        selectizeInput(inputId = "cancer", label = "cancer",
                                       choices = c("", CNViz::cbio_studies$Cancer),
                                       selected = ""),
                        DTOutput("cbioOutput")),
               tabPanel(icon("info-circle", class = NULL, lib = "font-awesome"),
                        fluid=TRUE,
                        h5("Patient Data"),
                        p("The visualization displayed is a representation of your input data. No additional inference has been done. Ploidy, if included, and copy numbers are rounded to the nearest whole number. Only one of probe, gene or segment data is required to launch the application. LOH and SNV data are optional. Please note, probe and gene data may not align with one another if the method used to generate your gene data was adjusted for sample purity and/or tumor ploidy."),
                        h5("TCGA Pan-Cancer Atlas Data"),
                        p("TCGA Pan-Cancer Atlas Data was obtained from cBioPortal's R package cBioPortalData. This same information can be found on cBioPortal's website. Similar information from many additional studies can also be found on their website. The copy number data displayed was generated by the GISTIC algorithm."),
                        p(" - Deep Deletion indicates a deep loss, possibly a homozygous deletion"),
                        p(" - Shallow Deletion indicates a shallow loss, possibley a heterozygous deletion"),
                        p(" - Gain indicates a low-level gain (a few additional copies, often broad)"),
                        p(" - Amplification indicate a high-level amplification (more copies, often focal)"),
                        div(style="display: inline-block;", a("cBioPortal website,", target ="_blank", href = "https://www.cbioportal.org/")),
                        div(style="display: inline-block;", a("cBioPortal FAQ,", target ="_blank", href = "https://docs.cbioportal.org/1.-general/faq")),
                        div(style="display: inline-block;", a("cBioPortalData R Package,", target ="_blank", href = "https://bioconductor.org/packages/release/bioc/html/cBioPortalData.html")),
                        div(style="display: inline-block;", a("GISTIC paper", target ="_blank", href = "https://pubmed.ncbi.nlm.nih.gov/18077431/")))
    )
  ),
  server <- function(input, output, session) {

    output$meta <- renderTable(meta_data, bordered = TRUE)

    green <- "#009E73" # in range gene outline
    white <- "#FFFFFF" # in range gene color
    blue <- "#0072B2" # out of range gene color
    pink <- "#CC79A7" # mutation
    orange <- "#D55E00" # segment
    black <- "#000000" # LOH (gene or segment)

    # add copy number estimate to gene_data
    if(nrow(gene_data) > 0){
      gene_data <- dplyr::filter(gene_data, !(gene %in% c('Antitarget', "", ".")))
      gene_data$m <- (gene_data$start + gene_data$end)/2
      gene_data$cn <- round((2^gene_data$log2)*2)
      gene_data$copies <- ifelse(gene_data$cn == 1, " copy", " copies")
      if("loh" %in% colnames(gene_data)){
        gene_data$loh <- ifelse(is.na(gene_data$loh), FALSE, gene_data$loh)
      }
    }

    if(nrow(probe_data) > 0){
      probe_data <- dplyr::filter(probe_data, !(gene %in% c('Antitarget', "", ".")))
      probe_data$m <- (probe_data$start + probe_data$end)/2
      probe_data$loh <- rep(NULL, nrow(probe_data))
      if("weight" %in% colnames(probe_data)){
        probe_data$total_weight <- scales::rescale(probe_data$weight, c(5,80))
      } else {
        probe_data$total_weight <- rep(10, nrow(probe_data))
      }

    }

    # add copy number estimate to probe_data
    if("weight" %in% colnames(probe_data)){
      probe_by_gene <- probe_data %>% dplyr::group_by(chr, gene) %>%
        dplyr::summarise(s = min(start), e = max(end),
                         mean_log2 = stats::weighted.mean(log2, weight),
                         total_weight = sum(weight)) %>%
        dplyr::mutate(log2 = mean_log2, m = (s+e)/2) %>%
        dplyr::mutate(cn = round((2^log2)*2)) %>%
        dplyr::mutate(blue = as.numeric(log2 < -0.41 | log2 > 0.32),
                      copies = ifelse(cn == 1, " copy", " copies"))
      probe_by_gene$total_weight <- scales::rescale(probe_by_gene$total_weight, c(5,80))
    } else if(nrow(probe_data) > 0) {
      probe_by_gene <- probe_data %>% dplyr::group_by(chr, gene) %>%
        dplyr::summarise(s = min(start), e = max(end),
                         mean_log2 = mean(log2),
                         total_weight = n()) %>%
        dplyr::mutate(log2 = mean_log2, m = (s+e)/2) %>%
        dplyr::mutate(cn = round((2^log2)*2)) %>%
        dplyr::mutate(blue = as.numeric(log2 < -0.41 | log2 > 0.32),
                      copies = ifelse(cn == 1, " copy", " copies"))
      probe_by_gene$loh <- rep(NULL, nrow(probe_by_gene))
      probe_by_gene$total_weight <- scales::rescale(probe_by_gene$total_weight, c(5,80))
    } else probe_by_gene <- data.frame()

    # ensuring log2 values stay in bounds of plot, will not change integer copy number that is displayed
    gene_data$log2 <- ifelse(gene_data$log2 < -2.5, -2.5, gene_data$log2)
    gene_data$log2 <- ifelse(gene_data$log2 > 5, 5, gene_data$log2)
    if(nrow(probe_by_gene) > 0){
      probe_by_gene$log2 <- ifelse(probe_by_gene$log2 < -2.5, -2.5, probe_by_gene$log2)
      probe_by_gene$log2 <- ifelse(probe_by_gene$log2 > 5, 5, probe_by_gene$log2)
    }

    # by gene
    if(nrow(gene_data) > 0){
      by_gene <- gene_data
    } else if (nrow(probe_data) > 0){
      by_gene <- probe_by_gene
    } else {
      by_gene <- NULL
    }

    # if gene data does not have weight, take weight from number of probes
    if(nrow(probe_by_gene) > 0 & nrow(gene_data) > 0 & !("weight" %in% colnames(gene_data))){
      gene_data <- left_join(gene_data, (probe_by_gene %>% dplyr::select(gene, chr, total_weight)), by = c("gene", "chr"))
      gene_data$total_weight <- ifelse(is.na(gene_data$total_weight), 1, gene_data$total_weight)
      gene_data <- gene_data %>% dplyr::mutate(weight = total_weight)
    }

    # highlight segments that are outside of tumor ploidy; if none supplied, use ploidy of 2
    ploidy <- ifelse(is.null(meta_data$ploidy), 2, meta_data$ploidy)
    if(nrow(segment_data) > 0){
      if("loh" %in% colnames(segment_data)){
        segment_data_sig <- filter(segment_data, log2 < log((ploidy-0.5)/2,2) | log2 > log((ploidy+0.5)/2,2) | loh == TRUE)
      } else {
        segment_data_sig <- filter(segment_data, log2 < log((ploidy-0.5)/2,2) | log2 > log((ploidy+0.5)/2,2))
      }
    } else segment_data_sig <- data.frame()

    for(i in c(1:length(chromosomes))){

      if(nrow(gene_data) > 0){
        # use gene data
        if(length(unique(gene_data$gene))<500){sizeref = 0.4} else{sizeref = 0.7}
        if("weight" %in% colnames(gene_data)){
          gene_data$total_weight <- scales::rescale(gene_data$weight, c(5,80))
        } else {
          gene_data$total_weight <- rep(10, nrow(gene_data))
        }
        chr <- filter(gene_data, chr == chromosomes[i])
        if(nrow(snv_data)>0){
          chr <- chr %>% left_join(snv_data %>% dplyr::group_by(gene) %>% dplyr::summarise(mutation_present = TRUE), by = "gene")
          chr$mutation_present <- ifelse(is.na(chr$mutation_present), FALSE, TRUE)
        } else chr$mutation_present <- rep(FALSE, nrow(chr))
        if(!("loh" %in% colnames(gene_data))){
          chr$loh <- rep(FALSE, nrow(chr))
        }
        assign(chromosomes[i], chr)

      } else if(nrow(probe_data) > 0 & nrow(gene_data) == 0){
        # use probe data
        if(length(unique(probe_data$gene))<500){
          probe_sizeref = 0.4
        } else{ probe_sizeref = 0.7 }
        chr <- filter(probe_by_gene, chr == chromosomes[i])
        if(nrow(snv_data)>0){
          chr <- chr %>% left_join(snv_data %>% dplyr::group_by(gene) %>% dplyr::summarise(mutation_present = TRUE), by = "gene")
          chr$mutation_present <- ifelse(is.na(chr$mutation_present), FALSE, TRUE)
        } else chr$mutation_present <- rep(FALSE, nrow(chr))
        chr$loh <- rep(FALSE, nrow(chr))
        assign(chromosomes[i], chr)

      } else{
        #return empty dataframe
        assign(chromosomes[i], data.frame())
      }

      marker_colors <- ifelse(between(get(chromosomes[i])$log2, 0.32, 0.41), white, blue) # make blue if < 1.5 or > 2.5
      marker_colors <- ifelse(get(chromosomes[i])$loh == TRUE, black, marker_colors)
      marker_colors <- ifelse(get(chromosomes[i])$mutation_present == TRUE, pink, marker_colors)
      outline_colors <- ifelse(marker_colors == white, green, marker_colors)
      outline_colors <- ifelse(marker_colors == pink & get(chromosomes[i])$loh == TRUE, black, outline_colors)

      if(nrow(segment_data)>0 & "loh" %in% colnames(segment_data)){
        segment_data_sig$seg_color <- ifelse(segment_data_sig$loh == TRUE, black, orange)
      } else if(nrow(segment_data) > 0) {
        segment_data_sig$seg_color <- orange #error (have to add an else if)
      }

      if(nrow(segment_data_sig) > 0){
        chr_seg <- filter(segment_data_sig, chr == chromosomes[i])
        assign(paste0(chromosomes[i], "_seg"), chr_seg)
      } else assign(paste0(chromosomes[i], "_seg"), data.frame())

      out_of_range <- ifelse(get(chromosomes[i])$cn > 64, " - log-2 outside range of y axis", "") # flagging points that were brough into view

      xmax <- max(filter(CNViz::cytoband_data, chrom == chromosomes[i])$chromEnd)

      plot <- plot_ly(source = "a", type = 'scatter', mode = 'markers') %>%
        add_trace(x = get(chromosomes[i])$m,
                  y = get(chromosomes[i])$log2,
                  text = paste0(get(chromosomes[i])$gene, " (", get(chromosomes[i])$cn, get(chromosomes[i])$copies, ")", out_of_range),
                  hoverinfo = 'text',
                  marker = list(color = marker_colors,
                                line = list(color = outline_colors),
                                size = get(chromosomes[i])$total_weight,
                                sizemode = 'area',
                                sizeref = ifelse(exists("sizeref"), sizeref, 1)),
                  showlegend = FALSE) %>%
        add_segments(x = 0, xend = xmax,
                     y = log((ploidy-0.5)/2, 2), yend = log((ploidy-0.5)/2, 2), line = list(color = "gray", width = 1, dash = "dot"), showlegend = FALSE) %>%
        add_segments(x = 0, xend = xmax,
                     y = log((ploidy+0.5)/2, 2), yend = log((ploidy+0.5)/2, 2), line = list(color = "gray", width = 1, dash = "dot"), showlegend = FALSE) %>%
        layout(yaxis=list(title = "log(2) copy number ratio",
                          titlefont = list(size = 8),
                          range = c(-3.5, 6)),
               xaxis= list(range = c(0,xmax)))

      if(nrow(get(paste0(chromosomes[i], "_seg")))>0){
        for(j in 1:nrow(get(paste0(chromosomes[i], "_seg")))){
          plot <- plot %>% add_segments(x = get(paste0(chromosomes[i], "_seg"))$start[j],
                                        xend = get(paste0(chromosomes[i], "_seg"))$end[j],
                                        y = get(paste0(chromosomes[i], "_seg"))$log2[j],
                                        yend = get(paste0(chromosomes[i], "_seg"))$log2[j],
                                        text = paste0("segment (", get(paste0(chromosomes[i], "_seg"))$start[j], "-", get(paste0(chromosomes[i], "_seg"))$end[j], ")"),
                                        hoverinfo = 'text',
                                        line = list(color = get(paste0(chromosomes[i], "_seg"))$seg_color[j], width = 3), showlegend = FALSE)
        }
      }

      cytoband_chrom <- dplyr::filter(CNViz::cytoband_data, chrom == chromosomes[i])

      subplot <- plot %>% plotly::layout(
        annotations = list(x = 40e6 , y = 6, text = chromosomes[i], showarrow= FALSE),
        xaxis = list(range = c(0, 250e6), dtick = 100e6), yaxis = list(range(-3.5,6)))

      chr_plot <- plot %>%
        add_trace(x = cytoband_chrom$chromStart, y = 6, xaxis = 'x2', showlegend = FALSE, marker = list(size = 0.1), hoverinfo = 'skip') %>%
        plotly::layout(title = gsub("adj_", "", chromosomes[i]),
                       xaxis = list(range = c(0, max(cytoband_chrom$chromEnd)), zeroline = TRUE, showline = TRUE),
                       xaxis2 = list(range = c(0, max(cytoband_chrom$chromEnd)),
                                     ticktext = as.list(cytoband_chrom$name), tickvals = as.list(cytoband_chrom$chromStart),
                                     tickfont = list(size = 8), tickmode = "array", tickangle = 270, side = "top",
                                     overlaying = 'x', zeroline = TRUE, autorange = FALSE, matches = 'x'),
                       margin = list(t = 80))

      assign(paste0(chromosomes[i],"_plot"), chr_plot)
      assign(paste0(chromosomes[i], "_subplot"), subplot)

    }

    all_plot <- plotly::subplot(chr1_subplot, chr2_subplot, chr3_subplot,
                                chr4_subplot, chr5_subplot, chr6_subplot,
                                chr7_subplot, chr8_subplot, chr9_subplot,
                                chr10_subplot, chr11_subplot, chr12_subplot,
                                chr13_subplot, chr14_subplot, chr15_subplot,
                                chr16_subplot, chr17_subplot, chr18_subplot,
                                chr19_subplot, chr20_subplot, chr21_subplot,
                                chr22_subplot, chrX_subplot, chrY_subplot,
                                nrows=9, shareY = TRUE, shareX = TRUE) %>%
      layout(autosize = FALSE, height = 1200)

    plot_todisplay <- reactive({ get(paste0(input$chr, "_plot")) })

    output$chr_plot <- renderPlotly(plot_todisplay())

    observe({
      updateSelectInput(session, "chr", selected = as.character(by_gene[by_gene$gene==input$gene,]$chr[1]))
    })

    d <- reactive ({ event_data(event="plotly_click", source = "a")[[3]] })

    observeEvent(event_data(event = "plotly_click", source = "a"), {
      updateSelectizeInput(session, "gene", selected = by_gene[by_gene$m == d(),]$gene[1])
    })

    observeEvent(input$chr, {
      updateSelectizeInput(session, "gene", selected =
                             ifelse(by_gene[by_gene$gene == input$gene,"chr"][1] == input$chr, input$gene, ""))
    })

    probe_data_select <- eventReactive(input$gene,{
      if(nrow(probe_data)>0){
        filter(probe_data, gene == input$gene)
      } else data.frame()
    })

    snv_data_select <- eventReactive(input$gene,{
      if(nrow(snv_data)>0 & "start" %in% colnames(snv_data)){
        filter(snv_data, gene == input$gene)
      } else data.frame()
    })

    probe_plot_title <- reactive({ paste0("probe data: ", input$gene) })

    probe_data_check <- reactive({ nrow(probe_data_select()) > 0 })

    probe_plot_check <- reactive({ max(!(input$gene %in% c("", "all")), d()) >= 1 })

    probe_sizes <- reactive({
      if(nrow(probe_data_select()) > 1){
        probe_data_select()$total_weight
      } else if(nrow(probe_data_select()) == 1){
        10
      }
    })

    # probe plot
    output$selected_plot <- renderPlotly({
      req(probe_data_check(), probe_plot_check())
      plot_ly(type = 'scatter', mode = 'markers') %>%
        add_trace(x = probe_data_select()$m,
                  y = probe_data_select()$log2,
                  marker = list(
                    color='purple',
                    line = list(color = 'purple'),
                    size = probe_sizes(),
                    sizemode = "area",
                    sizeref = ifelse(exists("probe_sizeref"), probe_sizeref, 1)),
                  showlegend = FALSE) %>%
        add_trace(x = as.numeric(snv_data_select()$start),
                  y = rep(0, nrow(snv_data_select())),
                  marker=list(
                    symbol = 'x',
                    size = 10,
                    color = 'black'),
                  text = snv_data_select()$mutation_id,
                  hoverinfo = 'text',
                  showlegend = FALSE) %>%
        add_segments(x = min(probe_data_select()$start), xend = max(probe_data_select()$end), y = -0.41, yend = -0.41, line = list(color = "gray", width = 1, dash = "dot"), showlegend = FALSE) %>%
        add_segments(x = min(probe_data_select()$start), xend = max(probe_data_select()$end), y = 0.32, yend = 0.32, line = list(color = "gray", width = 1, dash = "dot"), showlegend = FALSE) %>%
        layout(
          title = probe_plot_title(),
          xaxis = list(tickfont = list(size = 6), range = min(probe_data_select()$s), max(probe_data_select()$e)),
          yaxis=list(tickfont = list(size = 6), title = "log(2) copy number ratio",
                     titlefont = list(size = 8),
                     range = c(min(min(probe_data_select()$log2)-1,0), max(max(probe_data_select()$log2)+1),0)),
          paper_bgcolor='#fafafa', plot_bgcolor='#fafafa',margin = list(t = 80))
    })

    gene_snvs <- eventReactive(input$gene,{
      if(nrow(snv_data) > 0){
        return(filter(snv_data, gene == input$gene))
      } else return(data.frame())
    })

    # if no probe data/plot, display text with copy number & loh information
    cn_text <- eventReactive(input$gene,{
      if(nrow(gene_data) > 0 & nrow(probe_data) == 0 & nchar(input$gene) > 1){
        paste0(input$gene, " (", round(filter(gene_data, gene == input$gene)$start/1e6,2), "M", "-", round(filter(gene_data, gene == input$gene)$end/1e6,2), "M", ")", ": ",
               round(filter(gene_data, gene == input$gene)$cn,2),
               ifelse(filter(gene_data, gene == input$gene)$cn == 1, " copy", " copies"))
      } else if(nrow(gene_data) > 0 & nchar(input$gene) > 1){
        paste0(input$gene, ": ",
               round(filter(gene_data, gene == input$gene)$cn,2),
               ifelse(filter(gene_data, gene == input$gene)$cn == 1, " copy", " copies"))
      } else if (nrow(probe_data) > 0 & nchar(input$gene) > 1){
        paste0(input$gene, ": ",
               round(probe_by_gene[probe_by_gene$gene == input$gene,]$cn[1], 2),
               ifelse(round(probe_by_gene[probe_by_gene$gene == input$gene,]$cn[1],2) == 1, " copy", " copies"))
      }
      else{ paste0("") }
    })
    loh_text <- eventReactive(input$gene,{
      if(nrow(gene_data) > 0 & nchar(input$gene) > 1){
        if("loh" %in% colnames(gene_data)){
          ifelse(filter(gene_data, gene == input$gene)$loh == TRUE, ", LOH", "")
        } else paste0("")
      }
    })
    output$gene_text <- reactive({ paste0(cn_text(), loh_text()) })

    # mutation table
    output$mutations <- DT::renderDT({
      if(nrow(gene_snvs()) > 0){
        return(DT::datatable(gene_snvs(),
                             escape = FALSE,
                             rownames = FALSE,
                             options = list(dom = 't', columnDefs = list(list(className = 'dt-center')))) %>%
                 DT::formatStyle(c(1:dim(gene_snvs())[2]), border = '1px solid #ddd'))
      } else return(data.frame())
    })

    output$mutations_title <- renderText({
      if(nrow(gene_snvs()) > 0){
        return("Mutations:")
      } else return("")
    })

    # TCGA data tab
    if(exists("cbio")){
      cbio_studyId <- reactive({CNViz::cbio_studies$studyId[CNViz::cbio_studies$Cancer == input$cancer]})
      cbio_table <- reactive({
        cBioPortalData::getDataByGenePanel(api = cbio,
                                           studyId = cbio_studyId(),
                                           genePanelId = "IMPACT468",
                                           molecularProfileId = paste0(cbio_studyId(), "_gistic"),
                                           sampleListId = paste0(cbio_studyId(), "_cna"))
      })
      cbio_dat <- reactive({ data.frame(cbio_table()[[1]], stringsAsFactors = FALSE) })
      output$cbioOutput <- DT::renderDT({
        req(nchar(input$cancer)>0)
        DT::datatable(cbio_dat() %>% dplyr::group_by(hugoGeneSymbol) %>%
                        dplyr::summarise(Gain = sum(value ==1)/n(),
                                         Amplification = sum(value == 2)/n(),
                                         ShallowDeletion = sum(value == -1)/n(),
                                         DeepDeletion = sum(value == -2)/n()),
                      rownames = FALSE) %>%
          DT::formatPercentage(c("Gain", "Amplification", "ShallowDeletion", "DeepDeletion"), 2)
      })
    }

    # karyotype diagram
    if(nrow(segment_data)>0){
      if("loh" %in% colnames(segment_data)){
        output$karyotype <- downloadHandler(
          filename = paste0(sample_name, "_karyotype.pdf"),
          content = function(file) {
            karyo_data <- segment_data %>% dplyr::select(chr, start, end, log2, loh) %>% dplyr::mutate(cn = round(2^log2*2))
            karyo_data$cn <- ifelse(karyo_data$cn > 6, 6, karyo_data$cn)
            granges <- GenomicRanges::makeGRangesFromDataFrame(karyo_data, keep.extra.columns = TRUE, ignore.strand = TRUE)
            pdf(file)
            kp <- karyoploteR::plotKaryotype("hg38", plot.type = 2)
            CopyNumberPlots::plotCopyNumberCalls(kp, cn.calls = granges, labels = "", label.cex = 0, cn.colors = "red_blue", loh.color = "green")
            graphics::legend(x = "bottomright", fill = c("#EE0000", "#FFC1C1", "#E0E0E0", "#B2DFEE", "#87CEFA", "#1E90FF", "#0000FF", "green"),
                             legend = c("0", "1", "2", "3", "4", "5", "6+", "loh"), title = "copies", bty = "n")
            dev.off()
          })
      } else {
        output$karyotype <- downloadHandler(
          filename = "karyotype.pdf",
          content = function(file) {
            karyo_data <- segment_data %>% dplyr::select(chr, start, end, log2) %>% dplyr::mutate(cn = round(2^log2*2))
            karyo_data$cn <- ifelse(karyo_data$cn > 6, 6, karyo_data$cn)
            granges <- GenomicRanges::makeGRangesFromDataFrame(karyo_data, keep.extra.columns = TRUE, ignore.strand = TRUE)
            pdf(file)
            kp <- karyoploteR::plotKaryotype("hg38", plot.type = 2)
            CopyNumberPlots::plotCopyNumberCalls(kp, cn.calls = granges, labels = "", label.cex = 0, cn.colors = "red_blue")
            graphics::legend(x = "bottomright", fill = c("#EE0000", "#FFC1C1", "#E0E0E0", "#B2DFEE", "#87CEFA", "#1E90FF", "#0000FF"),
                             legend = c("0", "1", "2", "3", "4", "5", "6+"), title = "copies", bty = "n")
            dev.off()
          })
      }
    }

  })

}
