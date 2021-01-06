packages <- c("shiny","TDAmapper","igraph","Matrix","ggplot2","pcaMethods","plyr","dplyr","missForest","mvdalab","readxl")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    library(x, character.only = TRUE)
  }
)

options(shiny.maxRequestSize=30*1024^2)

# ====================== FUNCTIONS ======================================================
# Mean Value by Columns (excluding NAs)
mean.col.na <- function(df){
  mean.col <- c()
  
  for (i in seq(from=1,to=ncol(df),by=1)){
    if (is.na(mean(df[,i],na.rm = TRUE)) != TRUE){
      mean.col <- c(mean.col, mean(df[,i],na.rm = TRUE))
    }
  }
  
  return(mean.col)
}

# Mean Value by Rows (excluding NAs)
mean.row.na <- function(df){
  mean.df.row <- c()
  
  for (i in seq(from=1,to=nrow(df),by=1)){
    if (is.na(mean(as.numeric(df[i,]),na.rm = TRUE)) != TRUE){
      mean.df.row <- c(mean.df.row, mean(as.numeric(df[i,]),na.rm = TRUE))
    }
  }
  
  return(mean.df.row)
}

# Standard Deviation by Rows (excluding NAs)
sd.row.na <- function(df){
  sd.df.row <- c()
  
  for (i in seq(from=1,to=nrow(df),by=1)){
    if (is.na(sd(as.numeric(df[i,]),na.rm = TRUE)) != TRUE){
      sd.df.row <- c(sd.df.row, sd(as.numeric(df[i,]),na.rm = TRUE))
    }
  }
  
  return(sd.df.row)
}

# Cosine Similarity with NAs
cosine.similarity.na <- function(df){
  cosine.matrix <- matrix(0,nrow(df),nrow(df))
  for (i in 1:nrow(df)){
    for (j in 1:nrow(df)){
      if (i >= j){
        norm.i <- sqrt(sum(df[i,][which(is.na(df[i,]) == FALSE)]^2))
        norm.j <- sqrt(sum(df[j,][which(is.na(df[j,]) == FALSE)]^2))
        cosine.matrix[i,j] <- sum(df[i,]*df[j,], na.rm = TRUE)/(norm.i*norm.j)
        cosine.matrix[j,i] <- cosine.matrix[i,j]
      }
    }
  }
  return(cosine.matrix)
}

# ==================================================================================================



shinyApp(
  ui = tagList(
    navbarPage(
      "OptiMissP", id = "tabs",
      tabPanel("Loading Data",
               sidebarPanel(
                 #useShinyalert(),
                 fileInput("file", "Upload your data (CVS file):"),
                 radioButtons("dataset_structure", "Proteins:",choices = c("Column-wise" = "col", "Row-wise" = "row"), selected = "col"),
                 radioButtons("impdown", "Choose:",choices = c("Upload Imputed Data" = "upimpdata", "Impute Data" = "impdata"), selected = "upimpdata"),
                 conditionalPanel("input.impdown == 'upimpdata'",
                                  fileInput("file2", "Upload Imputed Data:")
                 ),
                 conditionalPanel(condition = "input.impdown == 'impdata'", 
                                  radioButtons("imputationmethod", "Imputation Methods:",choices = c("Lowest Value" = "lvalue", "MissForests" = "missf", "Mice" = "mice", "Probabilistic PCA" = "ppca", "Expectation-Maximization Imputation" = "em")
                                  ),
                                 actionButton("impute", "Impute!"),
                                 downloadButton(outputId = "downloadData.csv", label = "Download Imputed Data")
                 )
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Protein's Missingness",
                            plotOutput(outputId = "hist_prot")
                   ),
                   tabPanel("Patient's Missingness",
                            plotOutput(outputId = "hist_pat")
                   )
                 )
               )
      )
    )
  ),
  server = function(input, output) {
    
    
    # ===== LOADING DATA
    
    # Dataframe 1 = Not Imputed Data
    radiob <- reactiveValues(translation = FALSE)
    observeEvent(input$dataset_structure, {
      # 0 will be coerced to FALSE
      # 1+ will be coerced to TRUE
      if (input$dataset_structure == "row")
        radiob$translation <- TRUE
      else
        radiob$translation <- FALSE
    })
    
    dataframe1 <- reactive({
      if (radiob$translation == FALSE)
      {return(read.csv(input$file$datapath))}
      else
      {df <- read.csv(input$file$datapath)
      row_n <- row.names(df)
      df <- t(df)
      colnames(df) <- row_n
      return(df)}
    })
    
    # Dataframe 2 = Imputed Data
    r <- reactiveValues(doImp = FALSE)
    observeEvent(input$impute, {
      # 0 will be coerced to FALSE
      # 1+ will be coerced to TRUE
      r$doImp <- input$impute
    })
    
    dataframe2 <- reactive({
      
      if (input$impdown == "upimpdata" && is.null(input$file2$datapath) == FALSE){
        if (radiob$translation == FALSE)
        {return(read.csv(input$file2$datapath))}
        else
        {df <- read.csv(input$file2$datapath)
        row_n <- row.names(df)
        df <- t(df)
        colnames(df) <- row_n
        return(df)}
      }
      
      if (input$impdown == "impdata"){
        
        if (r$doImp == FALSE){return(NULL)}
        
        isolate({
          df1 <- read.csv(input$file$datapath)
          
          if(input$imputationmethod == "lvalue"){
            df2 <- df1
            showNotification("Imputing Data with Lowest Value Imputation..", duration = NULL, type = "message", id = "imputwarn")
            # Lowest value imputation
            for (i in 1:nrow(df2)){
              for (j in 1:ncol(df2)){
                if (is.na(df2[i,j])){
                  df2[i,j] = 1
                }
              }
            }
            removeNotification(id = "imputwarn")
          }
          
          if (input$imputationmethod == "missf"){
            df2 <- df1
            showNotification("Imputing Data with MissForest..", duration = NULL, type = "message", id = "imputwarn")
            missForest_imp <- missForest(df2)
            df2 <- missForest_imp$ximp
            removeNotification(id = "imputwarn")
          }
          
          if (input$imputationmethod == "mice"){
            df2 <- df1
            showNotification("Imputing Data with Mice..", duration = NULL, type = "message", id = "imputwarn")
            packages <- c("mice")
            
            package.check <- lapply(
              packages,
              FUN = function(x) {
                if (!require(x, character.only = TRUE)) {
                  install.packages(x, dependencies = TRUE)
                }
                library(x, character.only = TRUE)
              }
            )
            
            mice_imp <- mice(df2,m=10,maxit=100,meth='pmm',seed=500)
            df2 <- complete(mice_imp)
            removeNotification(id = "imputwarn")
          }
          
          if (input$imputationmethod == "ppca"){
            df2 <- df1
            showNotification("Imputing Data with Probabilistic PCA..", duration = NULL, type = "message", id = "imputwarn")
            df2 <- as.data.frame(pca(as.matrix(df1), method="ppca", nPcs=3)@completeObs)
            removeNotification(id = "imputwarn")
          }
          
          if (input$imputationmethod == "em"){
            df2 <- df1
            showNotification("Imputing Data with Expectation Maximization algorithm for imputation..", duration = NULL, type = "message", id = "imputwarn")
            df2 <- imputeEM(df1, scale = FALSE)$Imputed.DataFrames[[1]]
            removeNotification(id = "imputwarn")
          }
          
          
          return(df2)
          
        })
        
      }
      
    })
    
    # Histograms
    output$hist_prot <- renderPlot({
      if (is.null(input$file$datapath) == FALSE){
        df1 <- dataframe1()
        
        # Frequecy on missing values for each protein
        names.col <- colnames(df1)
        number.missing.values.protein <- c()
        for (i in 1:ncol(df1)){
          number.missing.values.protein <- c(number.missing.values.protein,sum(is.na(df1[,i])))
        }
        
        ggplot() + aes(number.missing.values.protein) +
          geom_histogram(breaks=seq(0, nrow(df1), by = 1),
                         col="#00AFBB",
                         fill="#00AFBB",
                         alpha = .2) +
          labs(title="Histogram of Missing Values for each Protein") +
          labs(x="Number of Missing Values", y="Count")
      }
    })
    output$hist_pat <- renderPlot({
      if (is.null(input$file$datapath) == FALSE){
        df1 <- dataframe1()
        
        #showNotification("Loading..", duration = NULL, type = "message", id = "histpat")
        
        # Frequency of missing values for each patient
        names.row <- rownames(df1)
        number.missing.values.patient <- c()
        for (i in 1:nrow(df1)){
          number.missing.values.patient <- c(number.missing.values.patient,sum(is.na(df1[i,])))
        }
        
        ggplot() + aes(number.missing.values.patient) +
          geom_histogram(breaks=seq(0, ncol(df1), by = 1),
                         col="#FC4E07",
                         fill="#FC4E07",
                         alpha = .2) +
          labs(title="Histogram of Missing Values for each Patient") +
          labs(x="Number of Missing Values", y="Count")
        
        #removeNotification(id="h")
      }
    })
    
    # Download Data
    output$downloadData.csv <- downloadHandler(
      filename = function() {
        paste('data-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(dataframe2(), con, row.names = FALSE)
      }
    )

    
    # ===== DATA for DISTRIBUTIONS and TDA
    # Dataframe 2: Uploaded
    
    # CONTROLS
    observeEvent(input$file, {
      observeEvent(input$file2, {
        
        # Dataframes
        df1 <- read.csv(input$file$datapath)
        
        # Imputation/Not Imputation
        df2 <- read.csv(input$file2$datapath)
        
        # Errors
        if(nrow(df1) != nrow(df2)){
          shinyalert("Oops!", "The two datasets have a different number of rows.", type = "error")
        }
        
        if(ncol(df1) != ncol(df2)){
          shinyalert("Oops!", "The two datasets have a different number of columns.", type = "error")
        }
        
        if(sum(colnames(df1) != colnames(df2))==ncol(df1)){
          shinyalert("Oops!", "The two datasets have different column names.", type = "error")
        }
        
      })
    })
    
    
    # ===== DISTRIBUTIONS
    
    # Population Tab
    observeEvent(input$file, {
      observeEvent(dataframe2(), {
        removeTab(inputId = "tabs", target = "Missingness Threshold")
        insertTab(inputId = "tabs",
                  tabPanel("Missingness Threshold",
                           sidebarPanel(
                             sliderInput("slider", "Missingness Threshold:", 1, 100, 60)
                           ),
                           mainPanel(
                             #h4("Density plot"),
                             plotOutput(outputId = "densityPlot"),
                             #h4("Details"),
                             htmlOutput("txtout1",container = tags$li),
                             htmlOutput("txtout4",container = tags$li),
                             htmlOutput("txtout5",container = tags$li),
                             htmlOutput("txtout6",container = tags$li),
                             htmlOutput("txtout7",container = tags$li)
                           )
                  ),
                  target = "Loading Data",
                  position = "after")
      })})
    
    # Protein Tab
    observeEvent(input$file, {
      observeEvent(dataframe2(), {
        removeTab(inputId = "tabs", target = "Protein")
        insertTab(inputId = "tabs",
                  tabPanel("Protein",
                           sidebarPanel(
                             selectInput(inputId = "prot", label = "Select a Protein:", choices = colnames(dataframe1()), multiple = FALSE)
                             ),
                           mainPanel(
                             plotOutput(outputId = "densityPlotProtein"),
                             verbatimTextOutput("textProt"))),
                  target = "Missingness Threshold",
                  position = "after")
      })})
  
    # Text in Tab 1
    # Text in Population Tab
    output$txtout1 <- renderText({
      if (is.null(input$file$datapath) == FALSE && is.null(dataframe2()) == FALSE){
        df1 <- dataframe1()
        df2 <- dataframe2()
        perc <- (100*sum(is.na(df1)))/(nrow(df1)*ncol(df1))
        perc2 <- sum(is.na(df2))/(ncol(df2)*nrow(df2))*100
        number.missing.values <- c()
        for (i in 1:ncol(df1)){
          number.missing.values <- c(number.missing.values, sum(is.na(df1[,i])))
        }
        perc.threshold <- input$slider/100
        indeces <- which(number.missing.values <= (perc.threshold*nrow(df1)))
        if(length(indeces) < 2){
          return()
        }else{
        dataset1 <- df1[,indeces]
        new.percentage.missingness <- (sum(is.na(dataset1))/(nrow(dataset1)*ncol(dataset1)))*100
        
        str1 <- paste("<b>","Percentage of not imputed data in: ","</b>")
        str2 <- paste("- the orginal not imputed dataset: ",as.character(round(perc)),"%")
        str3 <- paste("- the imputed dataset: ",as.character(round(perc2)),"%")
        str4 <- paste("- the subset of the original not imputed dataset at the given missingness threshold: ",as.character(round(new.percentage.missingness)),"%")
        HTML(paste(str1, str2,str3,str4, sep = '<br/>'))
        }
      }
    }) 
    output$txtout4 <- renderText({
      if (is.null(input$file$datapath) == FALSE){
        df1 <- dataframe1()
        number.missing.values <- c()
        for (i in 1:ncol(df1)){
          number.missing.values <- c(number.missing.values, sum(is.na(df1[,i])))
        }
        perc.threshold <- input$slider/100
        indeces <- which(number.missing.values <= (perc.threshold*nrow(df1)))
        if(length(indeces) < 2){
          showNotification("Not enough proteins to compute the operations", duration = NULL, type = "warning", id = "percthreshold")
          return()
        }else{
        dataset1 <- df1[,indeces]
        new.percentage.missingness <- (sum(is.na(dataset1))/(nrow(dataset1)*ncol(dataset1)))*100
        number.of.features.left <- ncol(dataset1)
        HTML(paste("<b>","Number of considered proteins: ","</b>",as.character(number.of.features.left)))
        }
      }
    })
    output$txtout5 <- renderText({
      if (is.null(input$file$datapath) == FALSE && is.null(dataframe2()) == FALSE){
        df1 <- dataframe1()
        df2 <- dataframe2()
        
        perc.threshold <- input$slider/100
        
        perc.missingness <- c()
        
        for (i in 1:ncol(df1)){
          perc.missingness <- c(perc.missingness, sum(is.na(df1[,i])))
        }
        
        threshold <- perc.threshold*nrow(df1)

        quant.indeces <- which(perc.missingness <= threshold)
        
        
        if(length(quant.indeces) < 2){
          return()
        }else{
        df12 <- df1[,quant.indeces]
        
        df22 <- df2[,quant.indeces]
        
        x <- c(as.matrix(df12))
        x <- x[which(is.na(x) == FALSE)]
        y <- c(as.matrix(df22))
        y <- y[which(is.na(y) == FALSE)]
        
        t <- ks.test(x,y)
        
        if(t$p.value < 0.05){
          pval <- "< 0.05"
        }
        else{
          pval <- as.character(signif(t$p.value, digits = 3))
        }
        if(t$p.value < 0.01){
          pval <- "< 0.01"
        }
        
        HTML(paste("<b>",t$method,"'s p-value: ","</b>",pval))
        }
      }
    })
    
    output$txtout6 <- renderText({
      if (is.null(input$file$datapath) == FALSE && is.null(dataframe2()) == FALSE){
        df1 <- dataframe1()
        df2 <- dataframe2()
        
        perc.threshold <- input$slider/100
        
        perc.missingness <- c()
        
        for (i in 1:ncol(df1)){
          perc.missingness <- c(perc.missingness, sum(is.na(df1[,i])))
        }
        
        threshold <- perc.threshold*nrow(df1)
        
        quant.indeces <- which(perc.missingness <= threshold)
        
        if(length(quant.indeces) < 2){
          return()
        }else{
        
        mean.qprot1 <- mean.col.na(df1[,quant.indeces])
        
        mean.qprot2 <- mean.col.na(df2[,quant.indeces])
        
        df <- data.frame("MeanProteinIntensity"=c(mean.qprot2,mean.qprot1), "Data" = c(rep("Imputed", length(mean.qprot1)), rep("Not Imputed", length(mean.qprot2))))
        
        qnt_imp <- quantile(df$MeanProteinIntensity[df$Data == "Imputed"])
        qnt_notimp <- quantile(df$MeanProteinIntensity[df$Data == "Not Imputed"])
        
        str1 <- paste("<b>","Quantiles: ","</b>")
        str2 <- paste("- Imputed data: <b>","0%: ","</b>",signif(qnt_imp[1], digits = 3), "| ", "<b>","25%: ","</b>",signif(qnt_imp[2], digits = 3), "| ", "<b>","50%: ","</b>",signif(qnt_imp[3], digits = 3), "| ", "<b>","75%: ","</b>",signif(qnt_imp[4], digits = 3), "| ", "<b>","100%: ","</b>",signif(qnt_imp[5], digits = 3))
        str3 <- paste("- Not imputed data: <b>","0%: ","</b>",signif(qnt_notimp[1], digits = 3), "| ", "<b>","25%: ","</b>",signif(qnt_notimp[2], digits = 3), "| ", "<b>","50%: ","</b>",signif(qnt_notimp[3], digits = 3), "| ", "<b>","75%: ","</b>",signif(qnt_notimp[4], digits = 3), "| ", "<b>","100%: ","</b>",signif(qnt_notimp[5], digits = 3))
        
        HTML(paste(str1, str3, str2, sep = '<br/>'))
        }
      }
    })
    
    output$txtout7 <- renderText({
      if (is.null(input$file$datapath) == FALSE && is.null(dataframe2()) == FALSE){
        df1 <- dataframe1()
        df2 <- dataframe2()
        
        perc.threshold <- input$slider/100
        
        
        perc.missingness <- c()
        
        for (i in 1:ncol(df1)){
          perc.missingness <- c(perc.missingness, sum(is.na(df1[,i])))
        }
        
        threshold <- perc.threshold*nrow(df1)
        
        quant.indeces <- which(perc.missingness <= threshold)
        
        if(length(quant.indeces) < 2){
          return()
        }else{
        
        mean.qprot1 <- mean.col.na(df1[,quant.indeces])
        
        mean.qprot2 <- mean.col.na(df2[,quant.indeces])
        
        X <- list(NotImputed=mean.qprot1, Imputed=mean.qprot2)
        
        df <- data.frame("MeanProteinIntensity"=c(mean.qprot2,mean.qprot1), "Data" = c(rep("Imputed", length(mean.qprot1)), rep("Not Imputed", length(mean.qprot2))))
        
        
        #mu <- ddply(df, "Data", summarise, grp.mean=median(MeanProteinIntensity))
        peak_imp <- density(df$MeanProteinIntensity[df$Data == "Imputed"])$x[which.max(density(df$MeanProteinIntensity[df$Data == "Imputed"])$y)]
        peak_notimp <- density(df$MeanProteinIntensity[df$Data == "Not Imputed"])$x[which.max(density(df$MeanProteinIntensity[df$Data == "Not Imputed"])$y)]
        diff <- abs(peak_notimp - peak_imp)
        
        str1 <- paste("<b>","Distributions' peaks: ","</b>")
        str2 <- paste("- Not imputed data: ",signif(peak_notimp, digits = 3))
        str3 <- paste("- Imputed data: ",signif(peak_imp, digits = 3))
        str4 <- paste("Distance between peaks: ",signif(diff, digits = 3))
        HTML(paste(str1, str2,str3,str4, sep = '<br/>'))
        }
      }
    })
    
    # Density Plots
    output$densityPlot <- renderPlot({ 
      
      # Dataset 1
      
      if (is.null(input$file$datapath) == FALSE){
        df1 <- dataframe1()
      }
      # Dataset 2
      if (is.null(dataframe2()) == FALSE){
        df2 <- dataframe2()
      }
      
      if (is.null(input$file$datapath) == FALSE && is.null(dataframe2()) == FALSE){
        # Basic function
        perc.threshold <- input$slider/100
        densityplots.funct <- function(df1,df2,perc.threshold){
          
          # ===== Libraries
          
          # Percentage of missing values for every protein
          perc.missingness <- c()
          
          
          for (i in 1:ncol(df1)){
            perc.missingness <- c(perc.missingness, sum(is.na(df1[,i])))
          }
          
          # ggplot() + aes(perc.missingness) +
          #   geom_histogram(breaks=seq(0, nrow(NURTuRE.data1), by = 1),
          #                  col="#00AFBB",
          #                  fill="#00AFBB",
          #                  alpha = .2) +
          #   labs(title="Frequency of Missing Values for each Protein") +
          #   labs(x="Protein", y="Count")
          
          # if(round(perc.threshold*nrow(df1),digits = 0) <= min(perc.missingness)){
          #   threshold <- min(perc.missingness)
          # }
          # else {
          #   threshold <- perc.threshold*nrow(df1)
          # }
          
          threshold <- perc.threshold*nrow(df1)
          quant.indeces <- which(perc.missingness <= threshold)
          
          if(length(quant.indeces) < 2){
            return()
          }else{
            
            quant.indeces <- which(perc.missingness <= threshold)
            # Mean in Nurture.data1
            mean.qprot1 <- mean.col.na(df1[,quant.indeces])
            
            # Mean in Nurture.data2
            mean.qprot2 <- mean.col.na(df2[,quant.indeces])
            
            # Overlap
            X <- list(NotImputed=mean.qprot1, Imputed=mean.qprot2)
            #p_overlap <- overlap(X)$OV
            
            
            df <- data.frame("MeanProteinIntensity"=c(mean.qprot2,mean.qprot1), "Data" = c(rep("Imputed", length(mean.qprot1)), rep("Not Imputed", length(mean.qprot2))))
            
            overlap_plots_sublist <- ggplot(df, aes(x=MeanProteinIntensity, fill=Data)) + geom_density(alpha=0.4)  + xlab("Mean Protein Intensity for each Patient") + ylab("Density") + ggtitle("Distribution of Patients' Mean Protein Intesity")
            return(overlap_plots_sublist)
            
          }
        
        }
        # Output
        densityplots.funct(df1,df2,perc.threshold)
      }
      else {
        return()
      }
      
      
    })
    
    # Density Plots - Single Protein
    output$densityPlotProtein <- renderPlot({ 
      # Dataset 1
      
      if (is.null(input$file$datapath) == FALSE){
        #df1 <- read.csv(input$file$datapath)
        df1 <- dataframe1()
      }
      # Dataset 2
      if (is.null(dataframe2()) == FALSE){
        #     df2 <- read.csv(input$file2$datapath)
        df2 <- dataframe2()
      }


      
      # Output
      if (is.null(input$file$datapath) == FALSE && is.null(dataframe2()) == FALSE){
        # Mean in Nurture.data1
        qprot1 <- df1[,input$prot][which(is.na(df1[,input$prot])==FALSE)]
        
        # Mean in Nurture.data2
        qprot2 <- df2[,input$prot][which(is.na(df2[,input$prot])==FALSE)]
        
        df <- data.frame("ProteinIntensity"=c(qprot2,qprot1), "Data" = c(rep("Imputed", length(qprot1)), rep("Not Imputed", length(qprot2))))
        
        ggplot(df, aes(x=ProteinIntensity, fill=Data)) + geom_density(alpha=0.4)  + xlab("Protein Intensity") + ylab("Density") + ggtitle(paste("Distribution of the Intensities of the", input$prot ,"Protein", sep=" "))
      }
        
    })
    
    # Text in Tab 2
    output$textProt <- renderText({
      if (is.null(input$file$datapath) == FALSE && is.null(dataframe2()) == FALSE){
        df1 <- dataframe1()
        perc <- (sum(is.na(df1[,input$prot]) == TRUE)/nrow(df1))*100
        paste("Percentage of missing values for the selected protein: ",as.character(round(perc)),"%")
      }
    })
    
    # ===== TDA
    
    # 2D TDA Tab
    observeEvent(input$file, {
      observeEvent(dataframe2(), {
        removeTab(inputId = "tabs", target = "2D TDA")
        insertTab(inputId = "tabs",
                  tabPanel("2D TDA",
                           sidebarPanel(
                             h4("Missingness"),
                             sliderInput("slidermiss", "Missingness Threshold:", 1, 100, 100),
                             radioButtons("enrichment", "Enrichment:",
                                          choices = c(None = "none",
                                                      Missingness = "missingness"),selected = "none"),
                             h6("Legend"),
                             h6("Pink: High Missingness"),
                             h6("Green: Low Missingness"),
                             h4("Resolution Parameters"),
                             sliderInput("X", "X Interval:", 3, 20, 5),
                             sliderInput("Y", "Y Interval:", 3, 20, 5),
                             sliderInput("overlap", "Percentage of Overlap:", 30, 60, 50, step = 10),
                             checkboxInput("cluster", "Optimal Clustering", FALSE),
                             h4("Advanced Parameters"),
                             selectInput(inputId = "lens1",
                                         label = "First Lens:",
                                         choices = c("Mean Intensity", "Standard Deviation (Intensity)", "L1 Infinity Centrality", "Mean Distance", "Standard Deviation (Distance)", "PPCA First Component"), selected = "Mean Intensity"),
                             selectInput(inputId = "lens2",
                                         label = "Second Lens:",
                                         choices = c("Mean Intensity", "Standard Deviation (Intensity)", "L1 Infinity Centrality", "Mean Distance", "Standard Deviation (Distance)", "PPCA First Component"), selected = "Standard Deviation (Intensity)"),
                             sliderInput("bins", "Single Linkage Clustering Parameter:", 6, 16, 10),
                             radioButtons("distmatrix", "Distance Matrix",
                                          choices = c(Euclidean = "euclidean",
                                                      Correlation = "correlation"),
                                          selected = "euclidean"),
                             actionButton("run", "Create Topologies!")
                             
                           ),
                           mainPanel(
                             #h4("Topology with Not Imputed Data:"),
                             plotOutput(outputId = "graph1"),
                             #h4("Topology with Imputed Data:"),
                             plotOutput(outputId = "graph2")
                           )
                           
                           
                  ),
                  target = "Protein",
                  position = "after")
      })})
    
    # Button
    v <- reactiveValues(doPlot = FALSE)
    observeEvent(input$run, {
      # 0 will be coerced to FALSE
      # 1+ will be coerced to TRUE
      v$doPlot <- input$run
    })
    
    # Graphs
    output$graph1 <- renderPlot({
      if (v$doPlot == FALSE) return()
      
      isolate({
      if (is.null(input$file$datapath) == FALSE){
        
        # Warning
        showNotification("Working on TDA..", duration = NULL, type = "message", id = "tdawarn")
          
          # Dataframe
          df1 <- read.csv(input$file$datapath)
          # Missingness threshold
          perc.missingness <- c()
          
          for (i in 1:ncol(df1)){
            perc.missingness <- c(perc.missingness, sum(is.na(df1[,i])))
          }
          
          
          perc.threshold <- input$slidermiss/100
          threshold <- perc.threshold*nrow(df1)
          quant.indeces <- which(perc.missingness <= threshold)
          df1 <- df1[,quant.indeces]
          
          # Distance matrix
          if (input$distmatrix == "euclidean"){
            d.matrix <- as.matrix(dist(df1, method = "euclidean"))
          }
          if (input$distmatrix == "correlation") {
            # Warning
            showNotification("Correlation Matrix in progress: this may take a while..", duration = 40, type = "message")
            d.matrix <- cosine.similarity.na(df1)
          }
          # Filter functions
          if (input$lens1 == "Mean Intensity"){
            lens1 <- mean.row.na(df1)
          }
          if (input$lens1 == "Standard Deviation (Intensity)"){
            lens1 <- sd.row.na(df1)
          }
          if (input$lens1 == "L1 Infinity Centrality"){
            if (input$distmatrix == "euclidean"){
              lens1 <- apply(d.matrix, 2, function(x) max(x, na.rm = TRUE))
            }
            if (input$distmatrix == "correlation") {
              d.eucl <- as.matrix(dist(df1, method = "euclidean"))
              lens1 <- apply(d.eucl, 2, function(x) max(x, na.rm = TRUE))
            }
          }
          if (input$lens1 == "Mean Distance"){
            if (input$distmatrix == "euclidean"){
              lens1 <- mean.row.na(d.matrix)
            }
            if (input$distmatrix == "correlation") {
              d.eucl <- as.matrix(dist(df1, method = "euclidean"))
              lens1 <- mean.row.na(d.eucl)
            }
          }
          if (input$lens1 == "Standard Deviation (Distance)"){
            if (input$distmatrix == "euclidean"){
              lens1 <- sd.row.na(d.matrix)
            }
            if (input$distmatrix == "correlation") {
              d.eucl <- as.matrix(dist(df1, method = "euclidean"))
              lens1 <- sd.row.na(d.eucl)
            }
          }
          if (input$lens1 == "PPCA First Component"){
            showNotification("PPCA in progress: this may take a while..", duration = 10, type = "message")
            ppca.components <- ppca(as.matrix(df1), nPCs=2, maxIterations = 100000)
            lens1 <- scores(ppca.components)[,1]
          }
          if (input$lens2 == "Mean Intensity"){
            lens2 <- mean.row.na(df1)
          }
          if (input$lens2 == "Standard Deviation (Intensity)"){
            lens2 <- sd.row.na(df1)
          }
          if (input$lens2 == "L1 Infinity Centrality"){
            if (input$distmatrix == "euclidean"){
              lens2 <- apply(d.matrix, 2, function(x) max(x, na.rm = TRUE))
            }
            if (input$distmatrix == "correlation") {
              d.eucl <- as.matrix(dist(df1, method = "euclidean"))
              lens2 <- apply(d.eucl, 2, function(x) max(x, na.rm = TRUE))
            }
          }
          if (input$lens2 == "Mean Distance"){
            if (input$distmatrix == "euclidean"){
              lens2 <- mean.row.na(d.matrix)
            }
            if (input$distmatrix == "correlation") {
              d.eucl <- as.matrix(dist(df1, method = "euclidean"))
              lens2 <- mean.row.na(d.eucl)
            }
          }
          if (input$lens2 == "Standard Deviation (Distance)"){
            if (input$distmatrix == "euclidean"){
              lens2 <- sd.row.na(d.matrix)
            }
            if (input$distmatrix == "correlation") {
              d.eucl <- as.matrix(dist(df1, method = "euclidean"))
              lens2 <- sd.row.na(d.eucl)
            }
          }
          if (input$lens2 == "PPCA First Component"){
            showNotification("PPCA in progress: this may take a while..", duration = 30, type = "message")
            ppca.components <- ppca(as.matrix(df1), nPCs=2, maxIterations = 100000)
            lens2 <- scores(ppca.components)[,1]
          }
          
          
          map <- mapper2D(distance_matrix = d.matrix, filter_values = list(lens1,lens2), num_intervals = c(input$X, input$Y), percent_overlap = input$overlap, num_bins_when_clustering = input$bins)
          graph <- graph_from_adjacency_matrix(map$adjacency,mode="undirected")
          
          # Vertex dimension
          vertex.attribute <- c()
          
          for (r in 1:length(map$points_in_vertex)){
            vertex.attribute <- c(vertex.attribute,length(map$points_in_vertex[[r]]))
          }
          
          vertex_attr(graph)$label <- vertex.attribute
          edge_attr(graph)$label <- rep(1,gsize(graph))
          
          #Size of each Vertex
          vertex.size <- rep(0,map$num_vertices)
          for (t in 1:map$num_vertices){
            points.in.vertex <-map$points_in_vertex
            vertex.size[t] <- length((map$points_in_vertex[[t]]))
          }
          V(graph)$size <- ( ((vertex.size-min(vertex.size))/(max(vertex.size)-min(vertex.size)) ) *15)+8
          
          # Enrichment
          if (input$enrichment == "missingness"){
            missing.freq <- c()
            for (b in 1:nrow(df1)){
              missing.freq <- c(missing.freq, sum(is.na(df1[b,]))/ncol(df1))
            }
            # 
            colfunc<-colorRampPalette(c("#3bff72","#ff3bc7"))
            # 
            miss.vertex <- c()
            for (v in 1:length(map$points_in_vertex)){
              miss.vertex <- c(miss.vertex, mean(missing.freq[map$points_in_vertex[[v]]]))
            }
            
            enrich.class.df <- data.frame("vertex"=seq(from=1,to=map$num_vertices,by=1), "value"=miss.vertex)
            
            clrlookup<-data.frame(clrs= colfunc(length(unique(enrich.class.df$value))), value=sort(unique(enrich.class.df$value)))
            
            clrMap<-merge(enrich.class.df,clrlookup, by=c("value"), all.x=TRUE)
            clrMap <- clrMap[order(clrMap$vertex),]
            V(graph)$color <- paste(clrMap$clrs)
          }
          
          # Plot
          if(input$cluster == FALSE){
            plot(graph, layout = layout_with_kk(graph)) + title("Topology with Not Imputed Data")
            removeNotification(id = "tdawarn")
          }
          else {
            optClst<-cluster_optimal(graph, weights = NULL)
            NumbOfClustrs<-length(unique(optClst$membership))
            plot(optClst, graph) + title("Topology with Not Imputed Data")
            removeNotification(id = "tdawarn")
          }
          
      }
      })
      
    })
    output$graph2 <- renderPlot({
      if (v$doPlot == FALSE) return()
      
      isolate({
      if (is.null(dataframe2()) == FALSE){
        # Dataframe
        df1 <- dataframe2()
        
        # Dataframe
        df1.notimp <- dataframe1()
        # Missingness threshold
        perc.missingness <- c()
        
        for (i in 1:ncol(df1.notimp)){
          perc.missingness <- c(perc.missingness, sum(is.na(df1.notimp[,i])))
        }
        
        
        perc.threshold <- input$slidermiss/100
        threshold <- perc.threshold*nrow(df1.notimp)
        quant.indeces <- which(perc.missingness <= threshold)
        df1 <- df1[,quant.indeces]
        
        # Distance matrix
        if (input$distmatrix == "euclidean"){
          d.matrix <- as.matrix(dist(df1, method = "euclidean"))
        }
        if (input$distmatrix == "correlation") {
          d.matrix <- cosine.similarity.na(df1)
        }
        # Filter functions
        if (input$lens1 == "Mean Intensity"){
          lens1 <- mean.row.na(df1)
        }
        if (input$lens1 == "Standard Deviation (Intensity)"){
          lens1 <- sd.row.na(df1)
        }
        if (input$lens1 == "L1 Infinity Centrality"){
          if (input$distmatrix == "euclidean"){
            lens1 <- apply(d.matrix, 2, function(x) max(x, na.rm = TRUE))
          }
          if (input$distmatrix == "correlation") {
            d.eucl <- as.matrix(dist(df1, method = "euclidean"))
            lens1 <- apply(d.eucl, 2, function(x) max(x, na.rm = TRUE))
          }
        }
        if (input$lens1 == "Mean Distance"){
          if (input$distmatrix == "euclidean"){
            lens1 <- mean.row.na(d.matrix)
          }
          if (input$distmatrix == "correlation") {
            d.eucl <- as.matrix(dist(df1, method = "euclidean"))
            lens1 <- mean.row.na(d.eucl)
          }
        }
        if (input$lens1 == "Standard Deviation (Distance)"){
          if (input$distmatrix == "euclidean"){
            lens1 <- sd.row.na(d.matrix)
          }
          if (input$distmatrix == "correlation") {
            d.eucl <- as.matrix(dist(df1, method = "euclidean"))
            lens1 <- sd.row.na(d.eucl)
          }
        }
        if (input$lens1 == "PPCA First Component"){
          showNotification("PPCA in progress: this may take a while..", duration = 20, type = "message")
          ppca.components <- ppca(as.matrix(df1), nPCs=2, maxIterations = 100000)
          lens1 <- scores(ppca.components)[,1]
        }
        if (input$lens2 == "Mean Intensity"){
          lens2 <- mean.row.na(df1)
        }
        if (input$lens2 == "Standard Deviation (Intensity)"){
          lens2 <- sd.row.na(df1)
        }
        if (input$lens2 == "L1 Infinity Centrality"){
          if (input$distmatrix == "euclidean"){
            lens2 <- apply(d.matrix, 2, function(x) max(x, na.rm = TRUE))
          }
          if (input$distmatrix == "correlation") {
            d.eucl <- as.matrix(dist(df1, method = "euclidean"))
            lens2 <- apply(d.eucl, 2, function(x) max(x, na.rm = TRUE))
          }
        }
        if (input$lens2 == "Mean Distance"){
          if (input$distmatrix == "euclidean"){
            lens2 <- mean.row.na(d.matrix)
          }
          if (input$distmatrix == "correlation") {
            d.eucl <- as.matrix(dist(df1, method = "euclidean"))
            lens2 <- mean.row.na(d.eucl)
          }
        }
        if (input$lens2 == "Standard Deviation (Distance)"){
          if (input$distmatrix == "euclidean"){
            lens2 <- sd.row.na(d.matrix)
          }
          if (input$distmatrix == "correlation") {
            d.eucl <- as.matrix(dist(df1, method = "euclidean"))
            lens2 <- sd.row.na(d.eucl)
          }
        }
        if (input$lens2 == "PPCA First Component"){
          showNotification("PPCA in progress: this may take a while..", duration = 10, type = "message")
          ppca.components <- ppca(as.matrix(df1), nPCs=2, maxIterations = 100000)
          lens2 <- scores(ppca.components)[,1]
        }
        
        map <- mapper2D(distance_matrix = d.matrix, filter_values = list(lens1,lens2), num_intervals = c(input$X, input$Y), percent_overlap = input$overlap, num_bins_when_clustering = input$bins)
        graph <- graph_from_adjacency_matrix(map$adjacency,mode="undirected")
        
        # Vertex dimension
        vertex.attribute <- c()
        
        for (r in 1:length(map$points_in_vertex)){
          vertex.attribute <- c(vertex.attribute,length(map$points_in_vertex[[r]]))
        }
        
        vertex_attr(graph)$label <- vertex.attribute
        edge_attr(graph)$label <- rep(1,gsize(graph))
        
        #Size of each Vertex
        vertex.size <- rep(0,map$num_vertices)
        for (t in 1:map$num_vertices){
          points.in.vertex <-map$points_in_vertex
          vertex.size[t] <- length((map$points_in_vertex[[t]]))
        }
        V(graph)$size <- ( ((vertex.size-min(vertex.size))/(max(vertex.size)-min(vertex.size)) ) *15)+8
        
        # Enrichment
        if (input$enrichment == "missingness"){
          missing.freq <- c()
          for (b in 1:nrow(df1.notimp)){
            missing.freq <- c(missing.freq, sum(is.na(df1.notimp[b,]))/ncol(df1.notimp))
          }
          # 
          colfunc<-colorRampPalette(c("#3bff72","#ff3bc7"))
          # 
          miss.vertex <- c()
          for (v in 1:length(map$points_in_vertex)){
            miss.vertex <- c(miss.vertex, mean(missing.freq[map$points_in_vertex[[v]]]))
          }
          
          enrich.class.df <- data.frame("vertex"=seq(from=1,to=map$num_vertices,by=1), "value"=miss.vertex)
          
          clrlookup<-data.frame(clrs= colfunc(length(unique(enrich.class.df$value))), value=sort(unique(enrich.class.df$value)))
          
          clrMap<-merge(enrich.class.df,clrlookup, by=c("value"), all.x=TRUE)
          clrMap <- clrMap[order(clrMap$vertex),]
          V(graph)$color <- paste(clrMap$clrs)
        }
        
        # Plot
        if(input$cluster == FALSE){
          plot(graph, layout = layout_with_kk(graph)) + title("Topology with Imputed Data")
        }
        else {
          optClst<-cluster_optimal(graph, weights = NULL)
          NumbOfClustrs<-length(unique(optClst$membership))
          plot(optClst, graph) + title("Topology with Imputed Data")
        }
        
      }
      })
      
    })

  }
)