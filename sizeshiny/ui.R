library(shiny)
library(shinythemes)
# Define UI for application that draws a histogram
shinyUI(fluidPage(theme=shinytheme("flatly"),
 #theme = "bootstrap.min.css",title="Depletion",
                
 #Application title-;
  ## headerPanel(
  ##    HTML(
  ##           '
  ##          <img src="http://www.cfrn-rcrp.ca/styles/ccfrn_global/images/cfrn_rcrp_logo.png"  width="399" height="102"  alt="Smiley face" > 
  ##          <img src="http://brand.ubc.ca/files/2012/08/UBCLogo_Black1.png"  width="86" height="102"  alt="Smiley face" align="right" >
  ##           </div>'
  ##       )),

  # Application title
  titlePanel("Size Pyramid"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
      sidebarPanel(
                sliderInput("year",
                            "Year:",
                            sep="",
                  min = 2001,
                  max =2011,
                  value = 2001),

      sliderInput("month",
                  "Month:",
                  min = 1,
                  max =13,
                  value = 13)
    ),
    # Show a plot of the generated distribution
      mainPanel(
          tabsetPanel(type = "tabs",
                      tabPanel("Size Pyramid Count", plotOutput("countPlot")),
                      tabPanel("Size Pyramid Percent", plotOutput("percentPlot"))
    )))))
      
