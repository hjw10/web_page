#setwd("C:/Users/penguin/Desktop/学习资料/大数据/20240529-Shiny基础使用/myshiny")
library(shiny)
library(shinythemes)

options(shiny.maxRequestSize=60*1024^2)

ui <- fluidPage(

  titlePanel("Bioinformatics Analysis Website"),
  sidebarLayout(
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("datafile", "Upload Data File",
                multiple = FALSE),
      checkboxInput("header", "Header", TRUE),
      
      
      # Input: Select number of rows to display ----
      radioButtons("disp", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head"),
      tags$hr(),
      
      fileInput("meta_data", "Upload Meta Data File",
                multiple = FALSE),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header2", "Header", TRUE),
      
      
      # Input: Select number of rows to display ----
      radioButtons("disp2", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head"),
      tags$hr(),
      
      actionButton("analyze_button", "show diff_gene plot"),
      radioButtons("disp3", "Display",
                   choices = c(Volcano = "volcano",
                               Pheatmap = "pheatmap"),
                   selected = "volcano"),
      
      tags$hr(),
      
      fileInput("clinic_data", "Upload Clinic Data File",
                multiple = FALSE),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header4", "Header", TRUE),
      
      # Input: Select number of rows to display ----
      radioButtons("disp4", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head"),
      textInput("text", "gene to do survival analysis", "BTRC"),
      actionButton("survival_button", "Perform Survival Analysis")
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Input data matrix",tableOutput("gene_exp"),tableOutput("sample_style"),tableOutput("clinical_data")),
        tabPanel("Differential Gene Analysis", plotOutput("diff_gene_plot")),
        tabPanel("Survival Analysis", plotOutput("survival_plot"))
      )
    )
  )
)
library(shiny)
library(DT)
library(limma)
library(survival)

server <- function(input, output) {
  output$gene_exp <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    req(input$datafile)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        df <- read.csv(input$datafile$datapath,
                       header = input$header)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  })
  output$sample_style <- renderTable({

    req(input$meta_data)
    tryCatch(
      {
        df <- read.csv(input$meta_data$datapath,
                       header = input$header2)
      },
      error = function(e) {
        stop(safeError(e))
      }
    )
    
    if(input$disp2 == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  })
  output$clinical_data <- renderTable({
    
    req(input$clinic_data)
    tryCatch(
      {
        df <- read.csv(input$clinic_data$datapath,
                       header = input$header4)
      },
      error = function(e) {
        stop(safeError(e))
      }
    )
    
    if(input$disp4 == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  })

  output$diff_gene_analysis <- ({
    
  }) 
  output$diff_gene_plot <- renderPlot({
    if(input$analyze_button > 0) {
      gene_exp <-read.csv(input$datafile$datapath,
                          header = F)
      sample_type <- read.csv(input$meta_data$datapath,
                              header = T)
      colnames(gene_exp) <- gene_exp[1,]
      gene_exp <- gene_exp[-1,]
      rownames(gene_exp) <- gene_exp[,1]
      gene_exp <- gene_exp[-1]
      
      rownames(sample_type) <- sample_type[,1]
      sample_type <- sample_type[-1]
      sample_type[,1] <- as.factor(sample_type[,1])
      
      forceMatrixToInteger <- function(m){
        apply (m, c (1, 2), function (x) {
          (as.integer(x))
        })
      }
      
      gene_exp <- forceMatrixToInteger(gene_exp)
      
      library(DESeq2)
      dds <- DESeqDataSetFromMatrix(countData=gene_exp, 
                                    colData=sample_type,
                                    design = ~phenotype)
      dds <- DESeq(dds)
      res <- results(dds, contrast = c('phenotype','Case', 'Control'))
      res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
      #      write.table(res1, 'GBM_DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)
      
      #差异分析
      ##筛选差异表达基因
      #首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
      res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
      
      #log2FC≥1 & padj<0.01 标识 up，代表显著上调的基因
      #log2FC≤-1 & padj<0.01 标识 down，代表显著下调的基因
      #其余标识 none，代表非差异的基因
      res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'up'
      res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'down'
      res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'
      
      #输出选择的差异基因总表
      res1_select <- subset(res1, sig %in% c('up', 'down'))
      #      write.table(res1_select, file = 'GBM_DESeq2.select.txt', sep = '\t', col.names = NA, quote = FALSE)
      
      #根据 up 和 down 分开输出
      res1_up <- subset(res1, sig == 'up')
      res1_down <- subset(res1, sig == 'down')
      
      #      write.table(res1_up, file = 'GBM_DESeq2.up.txt', sep = '\t', col.names = NA, quote = FALSE)
      #      write.table(res1_down, file = 'GBM_DESeq2.down.txt', sep = '\t', col.names = NA, quote = FALSE)
      ##ggplot2 差异火山图
      library(ggplot2)
      if(input$disp3 == "volcano"){
      #默认情况下，横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 padj
      p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
        geom_point(size = 1) +  #绘制散点图
        scale_color_manual(values = c('deeppink', 'gray', 'green'), limits = c('up', 'none', 'down')) +  #自定义点的颜色
        labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'control vs case', color = '') +  #坐标轴标题
        theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
              panel.background = element_rect(color = 'black', fill = 'transparent'), 
              legend.key = element_rect(fill = 'transparent')) +
        geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
        geom_hline(yintercept = 2, lty = 3, color = 'black') +
        xlim(-14, 14) + ylim(0, 40)  #定义刻度边界
      
      p
      }
      #热图
      else{
      library(pheatmap)
      
      diff_gene <- intersect(rownames(gene_exp),rownames(head(res1[order(res1[,5]),],100)))
      diff_exp <- gene_exp[rownames(gene_exp) %in% diff_gene,]
      normalization<-function(x){
        return((x-min(x))/(max(x)-min(x)))}#设定normalization函数,将值映射到[0，1]之间
      diff_exp=t(scale(t(diff_exp)))
      
      bk = unique(c(seq(0,1, length=100)))#设定热图参数
      pheatmap(diff_exp,
               color=colorRampPalette(c("royalblue","white","firebrick3"))(100),
               show_colnames = F,
               cluster_rows = T,cluster_cols = T,
               annotation_col = sample_type,
               breaks=bk,cellwidth=4,cellheight=5,
               fontsize=5,
               width = 16,height = 7)#画热图
      }
    }
  },width = 1000,height = 800)

    
  output$survival_plot <- renderPlot({
    if(input$survival_button > 0) {
      gene_exp <-read.csv(input$datafile$datapath,
                          header = F)
      sample_type <- read.csv(input$meta_data$datapath,
                              header = T)
      colnames(gene_exp) <- gene_exp[1,]
      gene_exp <- gene_exp[-1,]
      rownames(gene_exp) <- gene_exp[,1]
      gene_exp <- gene_exp[-1]
      
      rownames(sample_type) <- sample_type[,1]
      sample_type <- sample_type[-1]
      sample_type[,1] <- as.factor(sample_type[,1])
      
      forceMatrixToInteger <- function(m){
        apply (m, c (1, 2), function (x) {
          (as.integer(x))
        })
      }
      
      gene_exp <- forceMatrixToInteger(gene_exp)
      
      library(DESeq2)
      dds <- DESeqDataSetFromMatrix(countData=gene_exp, 
                                    colData=sample_type,
                                    design = ~phenotype)
      dds <- DESeq(dds)
      res <- results(dds, contrast = c('phenotype','Case', 'Control'))
      res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
      #      write.table(res1, 'GBM_DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)
      
      #差异分析
      ##筛选差异表达基因
      #首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
      res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
      
      #log2FC≥1 & padj<0.01 标识 up，代表显著上调的基因
      #log2FC≤-1 & padj<0.01 标识 down，代表显著下调的基因
      #其余标识 none，代表非差异的基因
      res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'up'
      res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'down'
      res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'
      
      #输出选择的差异基因总表
      res1_select <- subset(res1, sig %in% c('up', 'down'))
      #      write.table(res1_select, file = 'GBM_DESeq2.select.txt', sep = '\t', col.names = NA, quote = FALSE)
      
      #根据 up 和 down 分开输出
      res1_up <- subset(res1, sig == 'up')
      res1_down <- subset(res1, sig == 'down')
      
      diff_gene <- intersect(rownames(gene_exp),rownames(head(res1[order(res1[,5]),],100)))
      diff_exp <- gene_exp[rownames(gene_exp) %in% diff_gene,]
      normalization<-function(x){
        return((x-min(x))/(max(x)-min(x)))}#设定normalization函数,将值映射到[0，1]之间
      diff_exp=t(scale(t(diff_exp)))
      
      #生存分析
      GBM_clinic <- read.csv(input$clinic_data$datapath,header = T)
      #dim(GBM_clinic)
      newid<-lapply(strsplit(colnames(gene_exp),'-'),function(i){paste0(i[1:3],collapse = '-')})
      newid<-sapply(1:ncol(gene_exp),function(i){newid[[i]]})
      
      diff_exp<-rbind(diff_exp,newid)
      
      colnames(GBM_clinic)[3]<-'newid'
      df_OS<-merge(GBM_clinic,t(diff_exp),by="newid")
      
      
      alive <- df_OS[which(df_OS$vital_status=="Alive"),]
      dead <- df_OS[which(df_OS$vital_status=="Dead"),]
      alive$time <- alive$days_to_last_follow_up
      dead$time <- dead$days_to_death
      
      all_patients <- rbind(alive,dead)
      
      patient_survive <- cbind(all_patients$newid,all_patients$time,all_patients$vital_status)
      colnames(patient_survive) <- c("newid","time","status")
      
      GBM_survive_matrix <- merge(t(diff_exp),patient_survive,by="newid")
      
      library(limma)
      GBM_survive_matrix= limma::avereps(GBM_survive_matrix,ID=GBM_survive_matrix$newid)
      GBM_survive_matrix[GBM_survive_matrix == "Alive"] = 0
      GBM_survive_matrix[GBM_survive_matrix == "Dead"] = 1
      rownames(GBM_survive_matrix) <- GBM_survive_matrix[,1]
      GBM_survive_matrix<-GBM_survive_matrix[,-1]
      GBM <- as.data.frame(GBM_survive_matrix)
      
      #生存曲线
      library(survival)
      library(survminer)
      
      specific_gene <- c()
      for(i in 1: nrow(GBM)){
        if(GBM[i,which(colnames(GBM) == input$text)] >= 0){
          specific_gene[i] = "high"
        }else specific_gene[i] = "low"
      }
      GBM <- cbind(GBM,specific_gene)
      fit <- survfit(Surv(as.numeric(time),as.numeric(status)) ~ specific_gene , data = GBM)
      ggsurvplot(
        fit,                     
        data = GBM,             
        risk.table = TRUE,       
        pval = TRUE,             
        conf.int = TRUE,         
        xlim = c(0,800),         # 横坐标轴范围，相当于局部放大
        xlab = "Time in days",   # 横坐标标题
        break.time.by = 100,     # 横坐标刻度
        ggtheme = theme_light(), 
        risk.table.y.text.col = T, # risk table文字注释颜色
        risk.table.y.text = FALSE # risk table显示条形而不是文字
      )
    }
  },width = 900,height = 600)
  
}

shinyApp(ui = ui, server = server)








