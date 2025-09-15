# app.R
library(shiny)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(circlize)
library(ggrepel)
library(DT)

options(shiny.maxRequestSize = 100*1024^2)  # allow 100 MB uploads

ui <- fluidPage(
  titlePanel("RNA-seq Differential Expression Analysis with DESeq2"),
  sidebarLayout(
    sidebarPanel(
      fileInput("count_file", "Upload Count Matrix", accept = c(".csv",".txt")),
      fileInput("coldata_file", "Upload Sample Info", accept = c(".csv",".txt")),
      actionButton("run_analysis", "Run DESeq2 Analysis"),
      hr(),
      downloadButton("download_all", "Download Annotated Results (All Genes)"),
      downloadButton("download_sig", "Download Significant Genes")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary", tableOutput("summary_table")),
        tabPanel("All Genes", DTOutput("all_table")),
        tabPanel("Significant Genes", DTOutput("sig_table")),
        tabPanel("PCA Plot",
                 plotOutput("pca_plot", height = "600px"),
                 downloadButton("download_pca", "Download PCA Plot")),
        tabPanel("MA Plot",
                 plotOutput("ma_plot", height = "600px"),
                 downloadButton("download_ma", "Download MA Plot")),
        tabPanel("Volcano Plot",
                 plotOutput("volcano_plot", height = "600px"),
                 downloadButton("download_volcano", "Download Volcano Plot")),
        tabPanel("Sample Distance Heatmap",
                 plotOutput("sample_heatmap", height = "800px"),
                 downloadButton("download_sample_heatmap", "Download Sample Distance Heatmap")),
        tabPanel("Top 20 Genes Heatmap",
                 plotOutput("top20_heatmap", height = "900px"),
                 downloadButton("download_top20_heatmap", "Download Top20 Heatmap"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  analysis <- eventReactive(input$run_analysis, {
    req(input$count_file, input$coldata_file)
    
    # ---------- read input ----------
    counts  <- read.csv(input$count_file$datapath, row.names = 1, check.names = FALSE)
    coldata <- read.csv(input$coldata_file$datapath, row.names = 1, check.names = FALSE)
    
    # checks
    if (!all(colnames(counts) == rownames(coldata))) {
      stop("Sample names in counts and coldata do not match. Check column/rownames.")
    }
    
    coldata$disease_state <- factor(coldata$disease_state)
    coldata$cell_type     <- factor(coldata$cell_type)
    
    # ---------- DESeq2 ----------
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ cell_type + disease_state)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds)
    
    # results with explicit contrast
    res_raw <- DESeq2::results(dds, contrast = c("disease_state", "MDS", "healthy donor"), alpha = 0.05)
    
    # LFC shrink
    shrink_type <- if(requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"
    coefs <- resultsNames(dds)
    coef_name <- grep("disease_state.*MDS", coefs, value = TRUE)[1]
    if(is.na(coef_name)) coef_name <- coefs[2]
    res_shr <- tryCatch(lfcShrink(dds, coef = coef_name, type = shrink_type), error=function(e) res_raw)
    
    res_df <- as.data.frame(res_shr)
    res_df$ensembl_id <- sub("\\..*","",rownames(res_df))
    
    # ---------- annotation ----------
    gene_map <- data.frame()
    try({
      mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
      gene_map <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","description"),
                        filters="ensembl_gene_id",
                        values=res_df$ensembl_id,
                        mart=mart)
    }, silent=TRUE)
    
    if(nrow(gene_map) > 0){
      res_annot <- merge(res_df, gene_map, by.x="ensembl_id", by.y="ensembl_gene_id", all.x=TRUE)
    } else {
      res_df$hgnc_symbol <- NA_character_
      res_df$description <- NA_character_
      res_annot <- res_df
    }
    
    res_annot$hgnc_symbol[is.na(res_annot$hgnc_symbol)] <- ""
    rownames(res_annot) <- make.unique(ifelse(res_annot$hgnc_symbol=="", res_annot$ensembl_id, res_annot$hgnc_symbol))
    
    # ---------- Significant genes ----------
    res_annot$regulation <- "Non-significant"
    res_annot$regulation[res_annot$padj < 0.05 & res_annot$log2FoldChange > 1] <- "Upregulated"
    res_annot$regulation[res_annot$padj < 0.05 & res_annot$log2FoldChange < -1] <- "Downregulated"
    res_annot$regulation <- factor(res_annot$regulation, levels=c("Upregulated","Downregulated","Non-significant"))
    
    res_sig <- subset(res_annot, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
    
    summary_stats <- data.frame(
      Total_Genes = nrow(res_annot),
      Significant_Genes = nrow(res_sig),
      Upregulated = sum(res_sig$log2FoldChange >= 1, na.rm=TRUE),
      Downregulated = sum(res_sig$log2FoldChange <= -1, na.rm=TRUE)
    )
    
    # ---------- VST ----------
    vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
    rownames(vsd) <- sub("\\..*","",rownames(vsd))
    
    # ---------- PCA ----------
    pcaData <- plotPCA(vsd, intgroup=c("disease_state","cell_type"), returnData=TRUE)
    percentVar <- round(100*attr(pcaData,"percentVar"))
    pcaData$group <- paste(pcaData$disease_state, pcaData$cell_type, sep="_")
    gg_pca <- ggplot(na.omit(pcaData), aes(PC1, PC2, color=group, shape=group)) +
      geom_point(size=4, alpha=0.9) +
      labs(x=paste0("PC1: ", percentVar[1], "% variance"),
           y=paste0("PC2: ", percentVar[2], "% variance"),
           title="PCA of RNA-seq Samples") +
      theme_bw(base_size=14) + theme(legend.position="bottom")
    
    # ---------- MA Plot ----------
    top_genes <- rbind(head(res_sig[order(-res_sig$log2FoldChange),],10),
                       head(res_sig[order(res_sig$log2FoldChange),],10))
    top_genes$label <- ifelse(top_genes$hgnc_symbol=="", top_genes$ensembl_id, top_genes$hgnc_symbol)
    
    gg_ma <- ggplot(res_annot, aes(x=baseMean, y=log2FoldChange, color=regulation)) +
      geom_point(alpha=0.5, size=1.5) +
      scale_x_log10() +
      scale_color_manual(values=c("Upregulated"="red","Downregulated"="blue","Non-significant"="gray70")) +
      geom_hline(yintercept=0, linetype="dashed") +
      ggrepel::geom_label_repel(data=top_genes, aes(label=label),
                                box.padding=0.4, point.padding=0.3, size=3, max.overlaps=Inf) +
      labs(title="MA Plot (Top up/down genes labeled)", x="Mean expression", y="Log2 Fold Change") +
      theme_bw(base_size=14)
    
    # ---------- Volcano Plot ----------
    vol_df <- res_annot
    vol_df$log10padj <- -log10(vol_df$padj + 1e-300)
    sig_genes_vol <- subset(vol_df, padj < 0.05 & abs(log2FoldChange) >= 1)
    top_up_v   <- head(sig_genes_vol[order(-sig_genes_vol$log2FoldChange), ], 10)
    top_down_v <- head(sig_genes_vol[order(sig_genes_vol$log2FoldChange), ], 10)
    top20_v <- rbind(top_up_v, top_down_v)
    vol_df$label <- ifelse(rownames(vol_df) %in% rownames(top20_v),
                           ifelse(vol_df$hgnc_symbol == "" | is.na(vol_df$hgnc_symbol), vol_df$ensembl_id, vol_df$hgnc_symbol),
                           NA)
    
    gg_volcano <- ggplot(vol_df, aes(x = log2FoldChange, y = log10padj, color = regulation)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "Non-significant" = "gray70")) +
      geom_vline(xintercept = c(-1,1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      ggrepel::geom_text_repel(data = subset(vol_df, !is.na(label)), aes(label = label), size = 3, max.overlaps = Inf) +
      labs(title = "Volcano (Top significant genes labeled)", x = "Log2 Fold Change", y = "-log10(padj)") +
      theme_bw(base_size = 14)
    
    # ---------- Sample distance heatmap ----------
    sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))
    rownames(sampleDistMatrix) <- colnames(vsd)
    colnames(sampleDistMatrix) <- colnames(vsd)
    annotation_df <- coldata[,c("disease_state","cell_type"),drop=FALSE]
    ann_colors <- list(
      disease_state = c("MDS"="blue","healthy donor"="lightgreen"),
      cell_type = c("Mesenchymal stem cells"="orange","Monocytes"="pink")
    )
    
    # ---------- Top20 Heatmap ----------
    top20 <- head(res_sig[order(res_sig$padj),], 20)
    if(nrow(top20) > 0){
      ensembl_rows <- intersect(top20$ensembl_id, rownames(vsd))
      mat <- assay(vsd)[ensembl_rows, , drop=FALSE]
      labels <- ifelse(top20[match(ensembl_rows, top20$ensembl_id), "hgnc_symbol"]=="",
                       ensembl_rows,
                       top20[match(ensembl_rows, top20$ensembl_id),"hgnc_symbol"])
      rownames(mat) <- make.unique(labels)
      mat_to_plot <- t(scale(t(mat)))
      mat_to_plot[is.na(mat_to_plot)] <- 0
    } else mat_to_plot <- matrix(NA, nrow=0, ncol=ncol(vsd))
    
    # ---------- return ----------
    list(
      res_annot=res_annot,
      res_sig=res_sig,
      summary=summary_stats,
      dds=dds,
      vsd=vsd,
      plots=list(
        pca=gg_pca,
        ma=gg_ma,
        volcano=gg_volcano,
        sample_heatmap=list(mat=sampleDistMatrix, ann=annotation_df, ann_colors=ann_colors),
        top20_heatmap=list(mat=mat_to_plot, ann=annotation_df, ann_colors=ann_colors)
      )
    )
  })
  
  # ---------- Tables ----------
  output$summary_table <- renderTable({ analysis()$summary })
  output$all_table <- renderDT({ datatable(analysis()$res_annot, options=list(scrollX=TRUE), rownames=TRUE) })
  output$sig_table <- renderDT({ datatable(analysis()$res_sig, options=list(scrollX=TRUE), rownames=TRUE) })
  
  # ---------- Plots ----------
  output$pca_plot <- renderPlot({ analysis()$plots$pca })
  output$ma_plot <- renderPlot({ analysis()$plots$ma })
  output$volcano_plot <- renderPlot({ analysis()$plots$volcano })
  
  # ---------- Heatmaps ----------
  output$sample_heatmap <- renderPlot({
    ph <- analysis()$plots$sample_heatmap
    pheatmap::pheatmap(ph$mat, annotation_col=ph$ann, annotation_colors=ph$ann_colors,
                       cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE,
                       number_format="%.2f", number_color="black",
                       show_rownames=TRUE, show_colnames=TRUE,
                       main="Sample-to-sample Euclidean distances")
  }, height=700)
  
  output$top20_heatmap <- renderPlot({
    th <- analysis()$plots$top20_heatmap
    if(nrow(th$mat)==0){
      plot.new(); text(0.5,0.5,"No significant top20 genes available", cex=1.2)
    } else {
      pheatmap::pheatmap(th$mat, annotation_col=th$ann, annotation_colors=th$ann_colors,
                         show_rownames=TRUE, fontsize_row=8,
                         cluster_rows=FALSE, cluster_cols=FALSE,
                         display_numbers=TRUE, number_color="black",
                         main="Top 20 DE Genes (row-scaled)")
    }
  }, height=900)
  
  # ---------- Downloads ----------
  output$download_all <- downloadHandler(
    filename="DESeq2_all_annotated.csv",
    content=function(file){ write.csv(analysis()$res_annot, file, row.names=TRUE) }
  )
  output$download_sig <- downloadHandler(
    filename="DESeq2_significant_annotated.csv",
    content=function(file){ write.csv(analysis()$res_sig, file, row.names=TRUE) }
  )
  output$download_pca <- downloadHandler(
    filename="PCA_plot.png",
    content=function(file){ ggsave(file, plot=analysis()$plots$pca, width=8, height=6, dpi=300) }
  )
  output$download_ma <- downloadHandler(
    filename="MA_plot.png",
    content=function(file){ ggsave(file, plot=analysis()$plots$ma, width=8, height=6, dpi=300) }
  )
  output$download_volcano <- downloadHandler(
    filename="Volcano_plot.png",
    content=function(file){ ggsave(file, plot=analysis()$plots$volcano, width=8, height=6, dpi=300) }
  )
}

shinyApp(ui, server)
