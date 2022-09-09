library(shiny)
library(shinythemes)
                                        # Define UI for application that draws a histogram
shinyUI(fluidPage(theme=shinytheme("journal"),
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
                  titlePanel("Spawner Recruit"),
                  sidebarLayout(
                      sidebarPanel(
                          radioButtons("SRindex", label = h3("Stock Recruit Relationship"),inline=TRUE,
                                       choices = list("Two SF" = 1, "Single SF" = 2,"Marine Dev" = 3),selected = 3),
                          radioButtons("DevSD", label = h3("Marine Dev SD"),inline=TRUE,
                                       choices = list("2" = 1, "1" = 2,".5" = 3),selected = 2),
                          radioButtons("PriorSD", label = h3("SD on prior"),inline=TRUE,
                                       choices = list("Default" = 1, "Default*.9" = 2,"Default*1.2" = 3),selected = 1),
                          radioButtons("EIS", label = h3("EIS"),inline=TRUE,
                                       choices = list("Alt1" = 1, "Alt2a" = 2, "Alt2b"=3, "Alt3a"=4,"Alt3b"=5,"Alt4"=6,"NAA"=7,"General"=8),selected = 7),
                          radioButtons("Mort", label = h3("Marine Mortality Regime"),inline=TRUE,
                                       choices = list("Sin" = 1, "Rand" = 2,"AR1"=3),selected = 3),
                          ## radioButtons("river", label = h3("River"),inline=TRUE,
                          ##              choices = list("Foster" = 1, "Greenpeter" = 2, "Detroit"=3),selected = 1),
                          radioButtons("ssize", label = h3("Sample Size"),inline=TRUE,
                                       choices = list("1000" = 1, "10000" = 2),selected = 1),
     
     radioButtons("worm", label = h3("Worm Plots or Histograms"),inline=TRUE,
                                       choices = list("Worm" = 1, "Hist" = 2, "Auxiliary"=3),selected = 2),
      
                         sliderInput("dpe",
                                       "DPE:",
                                       min = 0.0,
                                       max =1.0000,
                                      value =1.000),
      
                         sliderInput("dps",
                                       "DPS:",
                                       min = 0.0,
                                       max =1.0000,
                                      value =1.000),
                        sliderInput("marinescale",
                                      "Marine Mortality scaling:",
                                      sep="",
                                      min = 0.0001,
                                      max =2.0,
                                      value =1)

                     ),
                                        # Show a plot of the generated distribution
                      mainPanel(
                          tabsetPanel(type = "tabs",
                                        tabPanel("Fit Plots", plotOutput("fitP")),
                                      tabPanel("Estimates", dataTableOutput("zT")),
                                      tabPanel("Covariance", dataTableOutput("covT")),
                                      tabPanel("Correlation", dataTableOutput("corT")),
                                      tabPanel("Deviates", dataTableOutput("devT")),
                                      tabPanel("Data", dataTableOutput("dataT")),
                                      tabPanel("Ages", dataTableOutput("agesT")),
                                      tabPanel("Forcast Plots Foster",plotOutput("forcastPF")),
                                      tabPanel("Forcast Table Foster",dataTableOutput("forcastTF")),
                                      tabPanel("Forcast Plots Green Peter", plotOutput("forcastPG")),
                                      tabPanel("Forcast Table Green Peter", dataTableOutput("forcastTG")),
                                      tabPanel("Forcast Plots Detroit", plotOutput("forcastPD")),
                                      tabPanel("Forcast Table Detroit", dataTableOutput("forcastTD")),
                                      tabPanel("Forcast Plots Foster and Green Peter", plotOutput("forcastPJoint")),
                                      tabPanel("Forcast Table Foster and Green Peter", dataTableOutput("forcastTJoint"))
                                      ))
       
     )))
