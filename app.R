
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# The following code runs DESeq2 analysis in Rshiny and opens an app which shows DESeq2 analysis for the input data provided by the user.
#
#    http://shiny.rstudio.com/
#
library(DESeq2)
library(shiny)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
#source("helper.R")

# User interface.
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      
      fileInput("file1", "upload counts table",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv", ".tsv", ".rda")
      ),
      fileInput("file2", "upload column data",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      
      selectInput("analysis", "type of analysis:",c("SELECT","single factor" = "SF", "multifactor" = "MF", "PCA" = "PCA")),
      # show this panel if only type of analysis "single factor" is selected.
      downloadButton("downloadData", "Download"),
      conditionalPanel( condition = "input.analysis == 'SF'", selectInput("plot2", "plot type:", c("SELECT", "plotMA" = "Ma","plotcount" = "count", "volcanoplot" = "vplot"))),
      # show this panel if only type of analysis "PCA" is selected.
      conditionalPanel( condition = "input.analysis == 'PCA'", selectInput("PCAplot","plot type:",c("SELECT", "PCA PLOT" = "pcaplot"))),
      # show this panel if only type of analysis "Multi factor" is selected.
      conditionalPanel( condition = "input.analysis == 'MF'", selectInput("Heat","HeatMap type:",c("SELECT", "HeatMap Of count matrix" = "heatmapcount","HeatMap Of sample to sample dists" = "heatmapdist"))),
      # show this panel if only type of analysis "heatmap" is selected.
      conditionalPanel(condition = "input.Heat == 'heatmapcount'", radioButtons("Updown","TYPE:",c("upward"="up","downward"="down"))),
      tags$hr()),
    #selectInput("plot","plot type:", c("MA1"="plotMA","multi" = "multi"))),
   # Main Panel
     mainPanel(
      #tableOutput("contents"),
      #plotOutput("PCAplot", width = "90%", height = "700px" ),
      #plotOutput("plot",width = "90%", height = "700px", click = "plot_click"),
      column(8,conditionalPanel(condition = "input.Heat == 'heatmapcount'", sliderInput("slideheat", "select number of genes:",min = 10, max = 100, value = 1),checkboxInput("Gene","Gene Names",FALSE))),
      #column(8,conditionalPanel(condition = "input.Heat == 'heatmapcount'", checkboxInput("Gene","Gene Names",TRUE))),
      plotOutput("Plotgraph", width = "90%", height = "700px")
      #plotOutput("heat",width = "90%", height = "700px"),
      #dataTableOutput("matrix"),
      #textOutput("test"),
      #verbatimTextOutput("info")
      
    )
  )
  
)

#
server <- function(input, output){
# function will react for each change in input parameters 
# it creates Deseqdataset when the both count data and metadata are available.
# this deseqdataset can be used in the rest of server function to inbuilt functions in DESeq2,
  data <- reactive({
    infile1 <- input$file1
    infile2 <- input$file2
    if (is.null(infile2)){
      return(NULL)}
    #  a <- input$file1
    # data <- dds(infile1, infile2)
    cts <- as.matrix(read.csv(infile1$datapath,sep="\t",row.names="gene_id"))
    coldata <- read.csv(infile2$datapath, row.names=1)
    counts <- cts[, rownames(coldata)]
    x <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,design = ~ condition)
    data <- DESeq(x)
  }) 
  #   output$test <- renderPrint({
  #    data1 <- data()
  #  summary(data1)
  #})
  
  
  
  #####################
# this function reacts if only multifactor analysis is selected.
# create deseqdataset 
  dataMF  <- reactive({
    ddsMF <- data()
    if (is.null(ddsMF)){
      return(NULL)}
    design(ddsMF) <- formula(~ type + condition)
    ddsMF <- DESeq(ddsMF)
  }
  )
# reactive function which reacts for each change in conditional input of Single factor analysis.
# gets the selected input in to "x" variable.
  SF1 <- reactive({
    x <- input$plot2
    # y <- input$pcaplot
    if(is.null(x)){
      return(NULL)
    }
    SF1 <- x
  })
# function for plots in Single factor analysis.
  plot <- reactive({
    data1 <- data()
    if (is.null(data1)){
      return(NULL)}
    a <- SF1()
    if (is.null(a)){
      return(NULL)
    }
    res <- DESeq2::results(data1)
    if(a == "Ma"){
      DESeq2::plotMA(res)
    }
    if(a == "count"){
      plotCounts(data1, gene=which.max(res$padj), intgroup="condition")
    }
    if(a == "vplot"){
      a <- as.data.frame(res)
      a$threshold = as.factor(abs(a$log2FoldChange) >2 & a$padj<0.9)
      x <- ggplot(data = a, aes(res$log2FoldChange, -log10(res$padj), colour=threshold))+ geom_point(alpha = 0.4, size=1.5)+ xlim(c(-4,4)) + ylim(c(0,75))+theme_bw()
      print(x)
      }
  })
  PCAplot <- reactive({
    data1 <- data()
    if(is.null(data1)){
      return(NULL)
    }
  # log transformed values using rlog
    rld <- rlog(data1, blind = FALSE)
    
    if (input$PCAplot == "pcaplot"){
      pcaData <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      x <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
        geom_point(size=3) + theme_bw()+
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()
      print(x)
    }})
  output$info <- renderText({
    xy <- input$plot_click
    if(is.null(xy)){
      return(NULL)
    }
    data1 <- data()
    x <- xy$x
    y <- xy$y
    res <- DESeq2::results(data1)
    #paste(x ,y)
    #a <- grep("xy$y",res$log2FoldChange)
    #rownames(res)[a]
  })
  
  heat <- reactive({
    data1 <- data()
    if(is.null(data1)){
      return(NULL)
    }
    rld <- rlog(data1, blind = FALSE)
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    if(input$Updown == 'up'){
    select <- order(rowMeans(counts(data1,normalized=TRUE)),decreasing=TRUE)[1:input$slideheat]
    }
    else{
    select <- order(rowMeans(counts(data1,normalized=TRUE)),increasing=TRUE)[1:5]
      
    }
    df <- as.data.frame(colData(data1)[,c("condition","type")])
    if (input$Heat == 'heatmapcount'){
      G <- renderText({input$Gene})
      pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames= input$Gene,
               cluster_cols=FALSE, annotation_col=df)}
    if (input$Heat == 'heatmapdist'){
      pheatmap(sampleDistMatrix,
               clustering_distance_rows=sampleDists,
               clustering_distance_cols=sampleDists,
               col=colors)}      
  })
  output$Plotgraph <- renderPlot({
    ptlist <- list(plot(),PCAplot(),heat())
  })
}
shinyApp(ui, server)
