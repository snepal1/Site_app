########################################################################
library(shiny)
library(shinyWidgets)
library(readr)
library(DT)
library(ggplot2)
library(dplyr)
library(markdown)
# UI code
ui <- fluidPage(
  # Custom CSS styles
  tags$head(
    tags$style(HTML("
      body {
        background-color: #222222; /* Dark background color */
        color: white; /* Text color */
      }

      .btn {
        color: forestgreen; /* Set text color for all buttons */
      }

      .modal-content {
        background-color: #333333; /* Dark background color for modal dialogs */
        color: white; /* Text color for modal dialogs */
      }

      /* Custom styles for fileInput label */
      .custom-file-label::after {
        border-color: forestgreen;
        color: forestgreen;
      }

      /* Custom styles for selectInput dropdown */
      .selectize-input {
        color: forestgreen;
      }

      /* Custom style for the tab labels */
      .nav-tabs .nav-link {
        border-bottom: 4px solid #006400; /* Green line below all tabs */
      }

      /* Remove border from the active tab */
      .nav-tabs .nav-link.active {
        border-bottom: none;
      }

      /* Custom styles for the table */
      .dataTables_wrapper {
        color: white;
      }

      .dataTables_wrapper th {
        color: white;
        font-weight: bold;
      }

      .dataTables_wrapper td {
        color: white;
      }

      /* Custom style for the pagination text */
      .dataTables_wrapper .dataTables_info {
        color: white;
      }
    /* Custom styles for the active tab's bottom border */
    .nav-tabs .active, .nav-tabs .active:hover {
      border-bottom: 2px solid forestgreen;
    }
  " ))),
  
  titlePanel("Site Index Calculation"),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        inputId = "filedata",
        label = tags$span("Upload CSV File", style = "color: forestgreen;"),
        multiple = FALSE,
        accept = c(".csv"),
        buttonLabel = "Choose...",
        placeholder = "No file selected yet"
      ),   
      downloadButton("download_data", "Download Updated File"),
      actionButton("process_data", "Process Data"),
      selectInput(
        inputId = "equation",
        label = tags$span("Select Equation", style = "color: forestgreen;"),
        choices = c("DF_KP1_NC","PP_CR2_MMC","IC_LG1_CA","WF_CR2_CA","RF_KP1_CA",
                    "MC3_CR2_OMC","MC3_CR2_MMC"),
        selected = "DF_KP1_NC"
      ),
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Processing",
                 DTOutput("tb1")
        ),
        tabPanel("Help",
                 includeMarkdown("Help.Rmd")
        ),
        tabPanel("Site Index Plot",
                 plotOutput("site_index_plot")
        )
      ),
    )
  )
)

# Rest of the code remains the same...


# Rest of the code remains the same...


server <- function(input, output, session) {
  
  data <- reactiveVal(NULL)
  
  observeEvent(input$filedata, {
    data(read.csv(input$filedata$datapath, header = TRUE))
  })
  
  processed <- reactiveVal(FALSE)
  
  observeEvent(input$process_data, {
    processed(FALSE)
    showModal(modalDialog(
      title = "Processing",
      "Please wait while the data is being processed...",
      footer = NULL,
      closable = FALSE
    ))
    process_data()
    processed(TRUE)
    removeModal()
  })
  
  process_data <- function() {
    data_df <- data()
    
    equations <- c("DF_KP1_NC","IC_LG1_CA","WF_CR2_CA","RF_KP1_CA",
                   "MC3_CR2_OMC","MC3_CR2_MMC")
    for (equation in equations) {
      data_df[[equation]] <- NA
    }
    
    # Define a data frame to map species to their factors
    species_factors <- data.frame(
      Species = c("DF","DF-C", "RW", "WF", "HL", "SS", "RC", "PC", "IC", "PP", "SP", "LP","WP","RF","YW","MC","DP","KP",
                  "AL","TO","PM","CQ","MA","LA","LO","BO","MH","JU","GS"),
      Factor = c(1.08,1.08,1,1.11,1.2,1.2,1.2,1.2,1.2,1,1.16,1.2,1.2,1.11,1.9,1.8,1.6,1.6,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,
                 1.6,0.9)
    )
    
    for (i in 1:nrow(data_df)) {
      if (data_df$Species[i] == "DF") {
        b1 <- -0.01524
        b2 <- -4.194
        b3 <- 28.35
        L <- log(data_df$HT[i] - 4.5)
        Y <- log(1 - exp(b1 * data_df$Age[i]))
        R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
        data_df$`MC3_CR2_MMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
        
      }   else if (data_df$Species[i] == "DF-C") {
        b1 <- 1.221
        b2 <- -0.8755
        b3 <- 402.8
        R <- ((((data_df$Age[i])^b1)/(data_df$HT[i]-4.5))-b2)/(b3+(data_df$Age[i])^b1)
        #Y <- log(1 - exp(b1 * data_df$Age[i]))
       # R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
        data_df$`DF_KP1_NC`[i] <- 4.5 + ((50^b1) / (b2+(b3*R)+(R*(50)^b1)))
        
      } else if (data_df$Species[i] == "PP") {
        b1 <- -0.01524
        b2 <- -4.194
        b3 <- 28.35
        L <- log(data_df$HT[i] - 4.5)
        Y <- log(1 - exp(b1 * data_df$Age[i]))
        R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
        data_df$`MC3_CR2_MMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
        
      } else if (data_df$Species[i] == "SP") {
        b1 <- -0.01524
        b2 <- -4.194
        b3 <- 28.35
        L <- log(data_df$HT[i] - 4.5)
        Y <- log(1 - exp(b1 * data_df$Age[i]))
        R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
        data_df$`MC3_CR2_MMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
        
     } else if (data_df$Species[i] == "DP") {
        b1 <- -0.01684
        b2 <- -1.255
        b3 <- 12.53
        L <- log(data_df$HT[i] - 4.5)
        Y <- log(1 - exp(b1 * data_df$Age[i]))
        R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
        data_df$`MC3_CR2_OMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
        
         } else if (data_df$Species[i] == "GS") {
         b1 <- -0.01524
        b2 <- -4.194
        b3 <- 28.35
        L <- log(data_df$HT[i] - 4.5)
        Y <- log(1 - exp(b1 * data_df$Age[i]))
        R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
        data_df$`MC3_CR2_MMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
        
        } else if (data_df$Species[i] == "HL") {
        b1 <- -0.01524
        b2 <- -4.194
        b3 <- 28.35
        L <- log(data_df$HT[i] - 4.5)
        Y <- log(1 - exp(b1 * data_df$Age[i]))
        R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
        data_df$`MC3_CR2_MMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
        
         } else if (data_df$Species[i] == "IC") {
        b1 <- 234.1
        b2 <- 3.923
        b3 <- -1.237
        L <- data_df$HT[i] - 4.5
        Y <- exp(b3 * log(data_df$Age[i]))
        R <- ((L - b1) + sqrt((L - b1)^2 + 4 * b2 * Y * L)) / 2
        data_df$`IC_LG1_CA`[i] <- 4.5 + (b1 + R) / (1+((b2/R)*exp(b3*log(50))))
        } 
      
      else if (data_df$Species[i] == "WF") {
        b1 <- -0.02834
        b2 <- -4.336
        b3 <- 31.51
        L <- log(data_df$HT[i] - 4.5)
        Y <- log(1 - exp(b1 * data_df$Age[i]))
        R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
        data_df$`WF_CR2_CA`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
      
     } else if (data_df$Species[i] == "KP") {
           b1 <- -0.01684
           b2 <- -1.255
           b3 <- 12.53
           L <- log(data_df$HT[i] - 4.5)
           Y <- log(1 - exp(b1 * data_df$Age[i]))
           R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
           data_df$`MC3_CR2_OMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
           
         } else if (data_df$Species[i] == "LP") {
           b1 <- -0.01684
           b2 <- -1.255
           b3 <- 12.53
           L <- log(data_df$HT[i] - 4.5)
           Y <- log(1 - exp(b1 * data_df$Age[i]))
           R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
           data_df$`MC3_CR2_OMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
           
         }else if (data_df$Species[i] == "MC") {
           b1 <- -0.01684
           b2 <- -1.255
           b3 <- 12.53
           L <- log(data_df$HT[i] - 4.5)
           Y <- log(1 - exp(b1 * data_df$Age[i]))
           R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
           data_df$`MC3_CR2_OMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
           
         } else if (data_df$Species[i] == "RF") {
        b1 <- 1.741
        b2 <- -110.3
        b3 <- 20100
        R<- (((data_df$Age[i])^b1/(data_df$HT[i]-4.5))-b2)/(b3+(data_df$Age[i])^b1)
        data_df$`RF_KP1_CA`[i] <- 4.5 + ((50^b1) / (b2+(b3*R)+(R*(50)^b1)))
        
        } 
      
      else {
        b1 <- -0.01524
        b2 <- -4.194
        b3 <- 28.35
        L <- log(data_df$HT[i] - 4.5)
        Y <- log(1 - exp(b1 * data_df$Age[i]))
        R <- (((L - b2 * Y)) + sqrt((L - (b2 * Y))^2 - 4 * b3 * Y)) / 2
        data_df$`MC3_CR2_MMC`[i] <- 4.5 + (data_df$HT[i] - 4.5) * ((1 - exp(b1 * 50)) / (1 - exp(b1 * data_df$Age[i])))^(b2 + b3 / R)
        
        }
      
    }
    
    # Combine site index values from all equations into a new column
    data_df$Combined_Site_Index <- apply(data_df[, equations], 1, function(row) {
      non_na_values <- na.omit(row)
      if (length(non_na_values) > 0) {
        paste0(non_na_values, collapse = ", ")
      } else {
        NA
      }
    })
    
    data_df$Combined_Site_Index <- as.numeric(data_df$Combined_Site_Index)
    
    # Calculate adjusted values based on the factor
    data_df$Adjusted <- sapply(1:nrow(data_df), function(i) {
      species <- data_df$Species[i]
      factor_index <- which(species_factors$Species == species)
      if (length(factor_index) > 0) {
        factor_value <- species_factors$Factor[factor_index]
        adjusted_value <- data_df$Combined_Site_Index[i] * factor_value
        return(adjusted_value)
      } else {
        return(NA)
      }
    })
    
    
    
    # Add a new column listing the equation used for each site index
    data_df$Equation_Used <- apply(data_df[, equations], 1, function(row) {
      non_na_values <- na.omit(row)
      if (length(non_na_values) > 0) {
        names(non_na_values)
      } else {
        NA
      }
    })
    
    # Add a new column to count the number of observations within each group based on the first two columns
    data_df$Count <- ave(seq_along(data_df[[1]]), data_df[[1]], FUN = length)
    

    data(data_df)
  }
  
  output$tb1 <- renderDT({
    req(data())
    datatable(data(), options = list(pageLength = 10, dom = 't<"bottom"ip>'))
  })
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("Updated_", input$filedata$name)
    },
    content = function(file) {
      write.csv(data(), file, row.names = FALSE)
    }
  )
  
  
  
  output$site_index_plot <- renderPlot({
    req(data(), input$equation)
    
    avg_data <- data() %>%
      group_by(Age) %>%
      summarise(avg_site_index = mean(get(input$equation), na.rm = TRUE))
    
    ggplot(data = avg_data, aes(x = Age, y = avg_site_index)) +
      geom_point() +
      geom_smooth(method = "auto", se = FALSE) +
      labs(title = "Site Index Plot",
           x = "Age",
           y = "Average Site Index") +
      theme_minimal()
  })
  
  
  
  observeEvent(processed(), {
    if (processed()) {
      shinyjs::disable("process_data")
      shinyjs::enable("download_data")
    } else {
      shinyjs::enable("process_data")
      shinyjs::disable("download_data")
    }
  })
  
  shinyjs::useShinyjs()
  
  shinyjs::runjs(
    "
    $(document).on('shiny:busy', function() {
      $('#process_data').text('Processing...').attr('disabled', true);
    });
  
    $(document).on('shiny:idle', function() {
      $('#process_data').text('Process Data').attr('disabled', false);
    });
    "
  )
}

# Increase maximum upload file size to 50MB
options(shiny.maxRequestSize = 50*1024^2)  # 50MB

shinyApp(ui, server)
