#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(stars)
library(cowplot)
library(cubelyr)

# Read in all of the climate clines
gea_dir <- list.files("outputs/GEA_res/")
all_clines <- lapply(gea_dir,function(dir){
    read.table(paste0("outputs/GEA_res/",dir,"/climate_cline.tsv"),header=T)
})
names(all_clines) <- gea_dir

# Set dir labels
dir_labels <- lapply(gea_dir,function(x) return(x))
names(dir_labels) <- gea_dir

# Get climate names
# climate_labs <- colnames(all_clines[[1]])
# climate_labs <- climate_labs[!(climate_labs %in% c("Lat","Long"))]
climate_vars <- c("Annual Mean Temperature",
                  "Mean Diurnal Range",
                  "Isothermality",
                  "Temperature Seasonality",
                  "Max Temperature of Warmest Month",
                  "Min Temperature of Coldest Month",
                  "Temperature Annual Range",
                  "Mean Temperature of Wettest Quarter",
                  "Mean Temperature of Driest Quarter",
                  "Mean Temperature of Warmest Quarter",
                  "Mean Temperature of Coldest Quarter",
                  "Annual Precipitation",
                  "Precipitation of Wettest Month",
                  "Precipitation of Driest Month",
                  "Precipitation Seasonality",
                  "Precipitation of Wettest Quarter",
                  "Precipitation of Driest Quarter",
                  "Precipitation of Warmest Quarter",
                  "Precipitation of Coldest Quarter")

# Set climate labs
climate_labs <- c("mean_temp",
                  "mean_diurnal",
                  "isothermality",
                  "temp_seasonality",
                  "max_temp_warmest_month",
                  "min_temp_coldest_month",
                  "temp_range",
                  "mean_temp_wet_quarter",
                  "mean_temp_dry_quarter",
                  "mean_temp_warm_quarter",
                  "mean_temp_cold_quarter",
                  "annual_precip",
                  "precip_wet_month",
                  "precip_dry_month",
                  "precip_seasonality",
                  "precip_wet_quarter",
                  "precip_dry_quarter",
                  "precip_warm_quarter",
                  "precip_cold_quarter")
names(climate_vars) <- climate_labs

# Set up temp colours
temp_colors <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))

# Function Library --------------------------------------------------------
plot_dataset_map <- function(cline_input,climate_input){
    
    # Fetch relevant climate info
    tmp_stars <- read_stars(paste0("data/worldclim/wc2-5/bio",climate_input,".bil"))
    if(climate_input<12){
        tmp_stars <- tmp_stars/10
    }
    
    # Set the limits
    long_range <- max(cline_input$Long) - min(cline_input$Long)
    long_min <- min(cline_input$Long) - 0.2*long_range
    long_max <- max(cline_input$Long) + 0.2*long_range

    lat_range <- max(cline_input$Lat) - min(cline_input$Lat)
    lat_min <- min(cline_input$Lat) - 0.2*lat_range
    lat_max <- max(cline_input$Lat) + 0.2*lat_range
    
    # Filter the stars
    tmp_stars_part <- stars:::filter.stars(tmp_stars, x > long_min, x < long_max,y > lat_min,y < lat_max)
    
    # Make the climate plot
    # long_range=abs(max(cline_input$Long)-min(cline_input$Long))
    # lat_range=abs(max(cline_input$Lat)-min(cline_input$Lat))
    
    climate_plot <- ggplot() + 
        geom_stars(data = tmp_stars_part) +
        scale_fill_gradientn(name = climate_vars[climate_input],
                             colors = temp_colors(5),
                             #limits = c(-7, 32),
                             na.value = "white") +
        coord_equal() +
        # scale_x_discrete(expand = c(0, 0)) +
        # scale_y_discrete(expand = c(0, 0)) +
        theme_minimal() +
        theme(legend.position = "top",
              panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank())+
        scale_x_continuous(limits = c(long_min,long_max))+
        scale_y_continuous(limits = c(lat_min,lat_max))
        # scale_x_continuous(limits = c(long_min,long_max))+
        # scale_y_continuous(limits = c(min(cline_input$Lat)-0.2*lat_range,max(cline_input$Lat)+0.2*lat_range))
    
    climate_plot2 <- ggplot() + 
        geom_stars(data = tmp_stars,downsample = 10) +
        scale_fill_gradientn(name = climate_vars[climate_input],
                             colors = temp_colors(5),
                             #limits = c(-7, 32),
                             na.value = "white") +
        coord_equal() +
        # scale_x_discrete(expand = c(0, 0)) +
        # scale_y_discrete(expand = c(0, 0)) +
        theme_minimal() +
        theme(legend.position = "top",
              panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank())+
        geom_point(data=cline_input,aes(x=Long,y=Lat),colour="black",size=2,alpha=1)
    
    # Add all points
    plot_grid(
    climate_plot + geom_point(data=cline_input,aes(x=Long,y=Lat),colour="black",size=2,alpha=1),
    climate_plot2,
    ncol=2)
}

# Plot coefficient of variation
cv_calc <- function(x) {
    sd(x) / mean(x) * 100
}
plot_cv <- function(cline_input){
    
    # Calculate the coefficient of variation for all variables...
    tmp <- cline_input[,names(climate_vars)]
    all_cv <- abs(apply(tmp,2,cv_calc))
    
    # Plot as bars...
    to_plot <- data.frame(climate_var=climate_vars,
                          cv=all_cv)
    to_plot$climate_var_F <- factor(to_plot$climate_var,levels=rev(to_plot$climate_var))
    ggplot(to_plot,aes(y=climate_var_F,x=all_cv))+
        geom_bar(stat="identity")+
        theme_minimal()+
        theme(axis.title.y=element_blank())+
        labs(x="CV")
}


# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("RepAdapt Climate Viewer"),
    
    # Dataset 1 inputs
    fluidRow(
        sidebarLayout(
            sidebarPanel(
                selectInput("data1",
                            h3("Dataset 1:"),
                            choices = dir_labels,
                            selected = gea_dir[1]),
                selectInput("climate1",
                            h3("Climate 1:"),
                            choices = as.list(climate_vars),
                            selected = climate_vars[1]),
            ),
            
            # Show plots in tabs
            mainPanel(
                tabsetPanel(
                    tabPanel("Map", plotOutput("climate_fig1")), 
                    tabPanel("CV", plotOutput("cv_fig1"))
                )
            )
        )
    ),
    # Dataset 2 inputs
    fluidRow(
        sidebarLayout(
            sidebarPanel(
                selectInput("data2",
                            h3("Dataset 2:"),
                            choices = dir_labels,
                            selected = gea_dir[1]),
                selectInput("climate2",
                            h3("Climate 2:"),
                            choices = as.list(climate_vars),
                            selected = climate_vars[1]),
            ),
            
            # Show plots in tabs
            mainPanel(
                tabsetPanel(
                    tabPanel("Map", plotOutput("climate_fig2")), 
                    tabPanel("CV", plotOutput("cv_fig2"))
                )
            )
        )
    )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
    
    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    #     
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white')
    # })
    
    
    # Plot Maps ---------------------------------------------------------------
    output$climate_fig1 <- renderPlot({
        plot_dataset_map(all_clines[[input$data1]],which(climate_vars==input$climate1))
    })
    output$climate_fig2 <- renderPlot({
        plot_dataset_map(all_clines[[input$data2]],which(climate_vars==input$climate2))
    })    
    
    # Plot coefficient of variation for all -----------------------------------
    output$cv_fig1 <- renderPlot({
        plot_cv(all_clines[[input$data1]])
    })
    output$cv_fig2 <- renderPlot({
        plot_cv(all_clines[[input$data2]])
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
