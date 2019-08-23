# BAT 2.1
# By Mans Thulin
# mans@statistikkonsult.com

####################


library(shiny)

# Run helpers:
source("bat-helpers.R")


# Create UI:
ui <- navbarPage("BAT 2.1",
                 # The different tabs:
                 tabPanel("Data", 
                          sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(HTML("Welcome!<br>Please see <a href=”http://www.mansthulin.se/bat”>this blog post</a> for instructions on how to use BAT.<br>BAT is free software and comes with absolutely no warranty.<br>&nbsp;<br>"),
                                  
                                  # Input: Select a file ----
                                  fileInput("file1", "Choose CSV File",
                                            multiple = TRUE,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv")),
                                  
                                  # Horizontal line ----
                                  tags$hr(),
                                  
                                  # Input: Select separator ----
                                  radioButtons("sep", "Columns in csv file are separated by:",
                                               choices = c(Comma = ",",
                                                           Semicolon = ";",
                                                           Tab = "\t"),
                                               selected = ","),
                                  
                                  # Input: Select encoding ----
                                  radioButtons("encoding", "Encoding for CSV file:",
                                               choices = c(UTF8 = "UTF-8",
                                                           UTF16 = "UTF-16",
                                                           Latin1 = "Latin1",
                                                           Latin2 = "Latin2"),
                                               selected = "UTF-16"),
                                  
                                  # Input: Select decimal point symbol ----
                                  radioButtons("dec", "Decimal separator used in file:",
                                               choices = c(Point = ".",
                                                           Comma = ","),
                                               selected = ".")
                              ),
                              
                              # Main panel for displaying outputs ----
                              mainPanel(tableOutput("contents")) )),
                          
                 tabPanel("Blank wells",  
                          
                          sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                  
                                   uiOutput("blankWells")
                              ),
                              
                              # Main panel for displaying outputs ----
                              #mainPanel(uiOutput(outputId = "blankPlotScaled")
                              mainPanel(uiOutput(outputId = "blankPlot")

                              ))),

                 tabPanel("Reference wells",  
                          sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                  
                                 # sliderInput("mask", "Masking interval:",
                                 #              min = 0, max = 1, value = c(0.02,0.1)),
                                 
                                 radioButtons("blankAdj", "Which method should be used for adjusting using the blanks?",
                                              choices = list("Method 1" = 1, "Method 2" = 2,
                                                             "No adjustment" = 3), selected = 2),
                                 
                                 #actionButton("saveResults", "Save results"),
                                  
                                  # Horizontal line ----
                                  tags$hr(),     
                                  uiOutput("referenceWells")
                                 # uiOutput("emptyWells")
                              ),
                              
                              # Main panel for displaying outputs ----
                              mainPanel(uiOutput(outputId = "referencePlot")
                              ))),
                 
                 tabPanel("Other wells",  
                          sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                  
                                  sliderInput("mask2", "Masking interval:",
                                              min = 0, max = 1, value = c(0.02,0.1)),
                                  
                                  sliderInput("numRep", "Number of replicates:",
                                              min = 1, max = 16, value = 4),
                                  
                                  sliderInput("minR", "Minimum R value:",
                                              min = 0.9, max = 0.9999, value = 0.999, step=0.0001),
                                  
                                  actionButton("saveResults2", "Save results"),
                                  
                                  # Horizontal line ----
                                  tags$hr(),     
                                  uiOutput("remWells")
                              ),
                              
                              # Main panel for displaying outputs ----
                              mainPanel(
                                  # Message while plots loading:
                                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                       tags$div("Plotting growth curves... Please wait.",id="loadmessage")),
                                  
                                  uiOutput(outputId = "remPlot")
                              ))),
                 
                 tabPanel("Manual fitting",
                          sidebarLayout(

                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                  
                                  radioButtons("maskType", "How would you like to chose which points to include for the fit?",
                                               choices = list("Vertical masking interval" = 1, "Horizontal masking interval" = 2),
                                                         selected = 1),

                                  sliderInput("mask3", "Vertical masking interval (log scale):",
                                              min = log(0.0001), max = log(1), value = log(c(0.02,0.1))),

                                  uiOutput("xmask"),
                                  
                                   actionButton("saveResults3", "Save results")

                              ),

                              # Main panel for displaying outputs ----
                              mainPanel(uiOutput(outputId = "manPlotScaled")
                              ))),
                 
                 tabPanel("Results", 
                          sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                  
                                "Scroll right to see the full table of results.",
                                  # Horizontal line ----
                                  tags$hr(),
                                
                                # Download buttons
                                downloadButton("downloadData", "Download results in csv file"),
                                downloadButton("downloadData2", "Download plots in pdf file"),
                                
                                HTML("<br>&nbsp;<br><b>Citations</b><br>Please use the following to cite the use of BAT:<br>&nbsp;<br>Thulin, M. (2018). BAT: an online tool for analysing growth curves. Retrieved from http://www.mansthulin.se/bat/")
                                  
       
                              ),
                              
                              # Main panel for displaying outputs ----
                              mainPanel(tableOutput("calcVal2")) ))
                 
                 )     
                 

# Define server logic:
server <- function(input, output, session) {
    
    # Define variables:
    myVars<- reactiveValues(calculatedValues=data.frame(Well=c(NA),fittedValue=c(NA),R=c(NA),doubTime=c(NA),growthRate=c(NA),pointsUsed=c(NA),groupMeanDoubTime=c(NA),groupMeanGrowthRate=c(NA),groupGrowthRateSD=c(NA),warningMessages=c(NA),mask1=c(NA),mask2=c(NA)),masktype=c(NA))
    
    # Read and clean the uploaded file:
    filedata <- reactive({
        infile <- input$file1
        if (is.null(infile)){
            return(NULL)      
        }
     OD<-read.csv(infile$datapath, header = TRUE, sep = input$sep, fileEncoding=input$encoding, dec=input$dec)
       
        
        
        ########
        
        # Check if first column is empty "Blank" column (new standard)
        if(names(OD)[2]=="Blank") { OD<-OD[,c(1,3:ncol(OD))]; for(i in 2:ncol(OD)) { names(OD)[i]<-gsub("X","Well.",names(OD)[i]) } }
        
        ######
        
        # Remove empty columns:
        OD<-OD[,!is.na(OD[1,])]
        OD<-OD[,!(OD[1,]==0)]
        
        # Fix issues that occur if the Bioscreen runs for more than 24 hours by changing the dates:
        tempTime<-as.vector(OD$Time)
        for(i in 1:length(OD$Time))
        {
            if(substring(OD$Time[i],1,2)<24)
            {
                tempTime[i]<-paste("2000-01-01 ",OD$Time[i],sep="")
            }
            else
            {
                tempTime[i]<-paste("2000-01-02 ",fixTime(OD$Time[i]),substring(OD$Time[i],3,8),sep="")
            }
        }
        
        # Convert strings to time objects
        convertedTimes<-strptime(as.character(tempTime),"%Y-%m-%d %H:%M:%S")
        
        # Calculate time differences in minutes
        diffTimes<-rep(0,length(convertedTimes))
        for(i in 2:length(convertedTimes))
        {
            diffTimes[i]<-diffTimes[i-1]+as.numeric(convertedTimes[i]-convertedTimes[i-1])
        }
        OD$Time<-diffTimes
        
        OD
    })
    
    ############
    
    # After reading the file, create check boxes for blank, reference and empty wells:
    output$blankWells <- renderUI({
        df <- filedata()
        if (is.null(df)) return(NULL)
        items=names(df)[2:ncol(df)]
        names(items)=items
           checkboxGroupInput("blankWells", 
                           h4("Which wells are blank?"), 
                           choices = items,
                           selected = c("Well.101","Well.102"))
    })
    
    output$referenceWells <- renderUI({
        df <- filedata()
        if (is.null(df)) return(NULL)
        items=names(df)[2:ncol(df)]
        names(items)=items
        defaultChecked<-setdiff(items,input$blankWells)[1:4]
        checkboxGroupInput("referenceWells", 
                           h4("Which wells contain the reference strain?"), 
                           choices = items,
                           selected = defaultChecked)
    })
    
    output$remWells <- renderUI({
        df <- filedata()
        if (is.null(df)) return(NULL)
        items=names(df)[2:ncol(df)]
        names(items)=items
        defaultChecked<-setdiff(items,union(input$blankWells,input$referenceWells))
        checkboxGroupInput("remWells", 
                           h4("Which wells should be analysed?"), 
                           choices = defaultChecked,
                           selected = defaultChecked)
    }) 
    
    output$xmask <- renderUI({
        df <- filedata()
        if (is.null(df)) return(NULL)
        time<-df$Time
        sliderInput("xmask", "Horizontal masking interval:",
                    min = min(time), max = max(time), value = quantile(time,c(0.05,0.2)), step=1)
    }) 
    
    ############
    
    # Plot blank wells:
    get_plot_output_list <- function(max_plots, input_n) {
        df <- filedata()
        if (is.null(df)) return(NULL)
        whichWells<-unlist(input$blankWells)
        x<-df$Time
        
        # Insert plot output objects the list
        plot_output_list <- lapply(1:input_n, function(i) {
            plotname <- paste("plot", i, sep="")
        
            if(i %% 2 == 0)
            {
            plot_output_object <- plotOutput(plotname, height = 250, width = 250)
            plot_output_object <- renderPlot({
                whichWell1<-whichWells[i]
                whichWell2<-whichWells[i-1]
                y1<-eval(parse(text=paste("df$",as.character(whichWell1),sep="")))
                y2<-eval(parse(text=paste("df$",as.character(whichWell2),sep="")))
                par(mfrow=c(1,2),cex=1.1)
                plot(x,y2,main=whichWell2,xlab="Time",ylab="OD",type="l",lwd=2)
                plot(x,y1,main=whichWell1,xlab="Time",ylab="OD",type="l",lwd=2)
            })
            }
        })
        
        do.call(tagList, plot_output_list) # needed to display properly.
        
        return(plot_output_list)
    }
    
    observe({
            output$blankPlot <- renderUI({ get_plot_output_list(length(unlist(input$blankWells)), length(unlist(input$blankWells))) })
        })
    
    
    ############
    
    # Plot reference wells:
    
    get_plot_output_list2 <- function(max_plots, input_n) {
        df <- filedata()
        if (is.null(df)) return(NULL)
        whichWells<-unlist(input$referenceWells)
        
        # Adjust using blanks:
        df<-adjustBlanks(df,input$blankAdj,unlist(input$blankWells))
        
        # Insert plot output objects the list
        plot_output_list2 <- lapply(1:input_n, function(i) {
            plotname <- paste("plot", i, sep="")
            
            # Plot growth curves on log scale:
            if(i %% 2 == 0)
            {
            plot_output_object <- plotOutput(plotname, height = 250, width = 250)
            plot_output_object <- renderPlot({
                referencePlotter(whichWells[(i-1):i],df,0.02,0.1,myVars$calculatedValues,2)
            })
            }
        })
        do.call(tagList, plot_output_list2) # needed to display properly.
        
        return(plot_output_list2)
    }
    
    observe({
        output$referencePlot <- renderUI({ get_plot_output_list2(ceiling(length(unlist(input$referenceWells))/2), length(unlist(input$referenceWells))) })
    })
    
    ############
    
    # Plot other wells:
    get_plot_output_list3 <- function(max_plots, input_n) {
        df <- filedata()
        if (is.null(df)) return(NULL)
        whichWells<-unlist(input$remWells)
        
        # Adjust using blanks:
        df<-adjustBlanks(df,input$blankAdj,unlist(input$blankWells))
        
        N<-ceiling(length(whichWells)/2)
        
        # Insert plot output objects the list
         plot_output_list3 <- lapply(1:input_n, function(i) {
            plotname <- paste("plot", i, sep="")
            
            # Plot growth curves on log scale:
            if(i %% 2 == 0)
            {
                plot_output_object <- plotOutput(plotname, height = 250, width = 250)
                plot_output_object <- renderPlot({
                    referencePlotter(whichWells[(i-1):i],df,input$mask2[1],input$mask2[2],myVars$calculatedValues,2)
                })
            }
        })
        do.call(tagList, plot_output_list3) # needed to display properly.
        
        return(plot_output_list3)
    }
    
    observe({
        output$remPlot <- renderUI({ get_plot_output_list3(ceiling(length(unlist(input$remWells))/2), length(unlist(input$remWells))) })
    })

    observeEvent(input$saveResults2, {
        
        df <- filedata()
        df<-adjustBlanks(df,input$blankAdj,unlist(input$blankWells))
        
        myVars$calculatedValues<-referencePlotter2(unlist(input$referenceWells),unlist(input$remWells),df,input$mask2[1],input$mask2[2],input$numRep,myVars$calculatedValues)
        
    })
    
    ############
        
    # Plot wells for manual fitting:
    output$manPlot <- renderPlot({
        df <- filedata()
        if (is.null(df)) return(NULL)
        cvd<-myVars$calculatedValues
        whichWells<-cvd$Well[which(cvd$R<input$minR)]
        
        # Adjust using blanks:
        df<-adjustBlanks(df,input$blankAdj,unlist(input$blankWells))
        
        # Plot growth curves on log scale:
        if(length(whichWells>0)) { referencePlotter(whichWells[1],df,exp(input$mask3[1]),exp(input$mask3[2]),myVars$calculatedValues,1,input$maskType,input$xmask[1],input$xmask[2]) } else { plot(0,0,col="white",main="All wells now have acceptable R-values") }
        
    })
    
    output$manPlotScaled <- renderUI({
        plotOutput("manPlot", height = 500)
    })
    
    observeEvent(input$saveResults3, {
        
        df <- filedata()
        df<-adjustBlanks(df,input$blankAdj,unlist(input$blankWells))
        
        cvd<-myVars$calculatedValues
        whichWells<-cvd$Well[which(cvd$R<input$minR)]
        
        myVars$calculatedValues<-referencePlotter3(whichWells[1],unlist(input$referenceWells),unlist(input$remWells),df,exp(input$mask3[1]),exp(input$mask3[2]),input$numRep,myVars$calculatedValues,input$maskType,input$xmask[1],input$xmask[2])
        
    })

    ############
    
    # Print results:
    output$calcVal2 <- renderTable({
        
        myVars$calculatedValues
    
        
    })
        
    # Print first few rows of the CSV file:
    output$contents <- renderTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file will be shown.
        
        req(input$file1)
        
        df<-filedata()
        if (is.null(df)) return(NULL)
        
        return( head(df) )
        
    })
    
    # Downloadable csv with numerical results:
    output$downloadData <- downloadHandler(
        filename = function() {
            gsub(".csv","-results.csv",input$file1,ignore.case=TRUE)
        },
        content = function(file) {
            write.csv(myVars$calculatedValues, file, row.names = FALSE)
        }
    )
    
    ############
    
    # Download plots:
    output$downloadData2 <- downloadHandler(
        filename = gsub(".csv","-results.pdf", input$file1, ignore.case=TRUE),
        content = function(file) {
            
            # Read data:
            df <- filedata()
            
            # Initiate plot:
            k<-2
            if(input$numRep==3) { k<-3 }
            if(input$numRep==5) { k<-5 }
            pdf(file,width=6*k,height=12)
            par(mfrow=c(2,k),cex=1.1)
            
            # Plot blank wells:
            whichWells<-unlist(input$blankWells)
            x<-df$Time
            
            # Plot growth curves:
            for(i in 1:length(whichWells)){
                whichWell<-whichWells[i]
                y<-eval(parse(text=paste("df$",as.character(whichWell),sep="")))
                
                plot(x,y,main=paste("Blank:",whichWell),xlab="Time",ylab="OD",type="l",lwd=2)
            }
            
            # Which wells should be plotted apart from blanks?
            whichWells<-c(unlist(input$referenceWells),unlist(input$remWells))

            # Adjust using blanks:
            df<-adjustBlanks(df,input$blankAdj,unlist(input$blankWells))

            # Plot non-blank wells:
            finalPlotter(whichWells,df,myVars$calculatedValues)

            dev.off()
        },
        contentType="application/pdf"
    )
    
    session$onSessionEnded(stopApp)
    
}

# Create Shiny app ----
shinyApp(ui, server)
