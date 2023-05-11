# Load packages
library(shiny)
library(tidyverse)
library(tibble)
library(EnvStats)         # to load rtri function 
library(truncnorm)        # to load rtruncnorm function
library(jmuOutlier)       # to load rlaplace function
library(formula.tools)    # to load 'rhs'
library(purrr)            # to load 'map'
library(deSolve)          # to load 'ode'
library(rlang)
library(ggplot2)

# Load utility functions
source("UtilityFunctions_dynamic_growth.R")

# Generate database 
BTyper3_input = read.csv("Btyper3_Results.csv")
colnames(BTyper3_input)[1] <- "Isolate.Name"
gp_input = read.csv("simulation_input.csv")
database = cbind(BTyper3_input,gp_input[,3:7])
database <- database %>% 
  separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
  mutate(ANI = gsub("\\)", "", ANI))

# Define ui
ui <- fluidPage(
  titlePanel("BRisk"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n0", "Initial count:", value = 100),  # Numeric input for "initial count"
      numericInput("d", "Storage day:", value = 35),  # Numeric input for "storage day"
      fileInput("file", "Input BTyper3 result for a detected B cereus isolate")  # BTyper3 input for a B cereus isolate
    ),
    
    mainPanel(
      plotOutput("hist")
    )
  )
)

# Define server
server <- function(input, output) {
  
  # Input BTyper3 result for a B cereus isolate 
  data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath)
    
    # Filter the database for rows with the same species as the BTyper3 input
    df$species <- trimws(df$species)
    matching_species_df <- subset(database, species == df$species)
    
    # Simulate 1000 units of HTST milk along the supply chain 
    ## Set up dataframe for modeling 1000 units of HTST milk 
    n_sim = 1000
    ModelData = data.frame(unit_id = rep(seq(1,n_sim)))
    
    ## Set seed
    set.seed(1)
    
    # Stage 1: facility storage 
    ## (a)  Sample the temperature distribution
    ModelData$T_F <- rep(runif(n_sim,min=3.5,max=4.5)) #uniform distribution
    ## (b) Sample the storage time (in days) distribution
    ModelData$t_F <- rep(runif(n_sim,min=1,max=2)) #uniform distribution
    
    # Stage 2: transport from facility to retail store
    ## (a)  Sample the temperature distribution
    ModelData$T_T <- rep(rtri(n_sim,min=1.7,max=10.0,mode=4.4)) #triangular distribution
    ## (b) Sample the transportation time (in days) distribution
    ModelData$t_T <- rep(rtri(n_sim,min=1,max=10,mode=5))
    
    # Stage 3: storage/display at retail store
    ## (a)  Sample the temperature distribution
    ModelData$T_S <- rep(rtruncnorm(n_sim,a=-1.4,b=5.4,mean=2.3,sd=1.8)) #truncated normal distribution
    ## (b) Sample the storage time (in days) distribution
    ModelData$t_S <- rep(rtruncnorm(n_sim,a=0.042,b=10.0, mean=1.821,sd=3.3)) #truncated normal distribution
    
    ## Stage 4: transportation from retail store to home
    ## (a)  Sample the temperature distribution
    ModelData$T_T2 <- rep(rtruncnorm(n_sim,a=0,b=10,mean=8.5,sd=1.0)) #truncated normal distribution
    ## (b) Sample the transportation time (in days) distribution 
    ModelData$t_T2 <- rep(rtruncnorm(n_sim,a=0.01,b=0.24, mean=0.04,sd=0.02)) #truncated normal distribution
    
    ## Stage 5: home storage 
    ## (a)  Sample the temperature distribution
    temps <- rep(NA, n_sim)
    for (i in 1:n_sim){
      number <- rlaplace(1,m=4.06,s=2.31)
      while (number > 15 | number < -1) {
        number <- rlaplace(1,m=4.06,s=2.31) #truncated laplace distribution 
      }
      temps[i] <- number
    }
    ModelData$T_H <- temps
    ## (b) Define t_H as 35 days for all units
    ModelData$t_H <- rep(input$d, each = n_sim)
    
    ## Model temperature profiles of 1000 units HTST milk 
    env_cond_time <- matrix(c(rep(0,1000),
                              ModelData$t_F, 
                              ModelData$t_F+0.001,
                              ModelData$t_F + ModelData$t_T,
                              ModelData$t_F + ModelData$t_T+0.001,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S+0.001,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2+0.001,
                              ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2 + ModelData$t_H), ncol = 10)
    
    env_cond_temp <- matrix(c(ModelData$T_F, 
                              ModelData$T_F,
                              ModelData$T_T,
                              ModelData$T_T,
                              ModelData$T_S,
                              ModelData$T_S,
                              ModelData$T_T2,
                              ModelData$T_T2,
                              ModelData$T_H,
                              ModelData$T_H), ncol = 10)
    
    # Assign serving size to 1000 units of HTST milk 
    serving.size<-sample(x = c(rep(x = 244,50),rep(245,25),rep(488,20),rep(732,5)),size = 1000,replace = TRUE)
    ModelData$serving.size = serving.size
    
    ## Randomly assign isolate names to 1000 units of HTST milk 
    ## All isolates from the same species are equally represented
    isolates <- matching_species_df$Isolate.Name
    samples_per_isolate <- n_sim / length(isolates)
    sampled_isolates <- character()
    for (isolate in isolates) {
      sampled_isolates <- c(sampled_isolates, rep(isolate, samples_per_isolate))
    }
    sampled_isolates <- sample(sampled_isolates)
    ModelData$isolate <- sampled_isolates
    
    ## Generate simulation input 
    ## Assign growth parameters to 1000 units of HTST milk 
    ModelData$index = match(ModelData$isolate, matching_species_df$Isolate.Name)
    ModelData$Q0 = matching_species_df$Q0[ModelData$index]
    ModelData$Nmax = matching_species_df$Nmax[ModelData$index]
    ModelData$b = matching_species_df$b[ModelData$index]
    ModelData$Tmin = matching_species_df$Tmin[ModelData$index]
    ModelData$Clade = matching_species_df$Clade[ModelData$index]
    
    ## Generate N0 from a Poisson distribution 
    set.seed(42)
    N0 = rpois(n = n_sim, lambda = input$n0)
    ModelData$N0 = N0 
    
    ModelData$Topt = sapply(ModelData$Clade, xopt_func)
    ModelData$mu_opt = (ModelData$b*(ModelData$Topt-ModelData$Tmin))^2
    
    # Run simulation
    for (i in 1:nrow(ModelData)){
      my_primary <- list(mu_opt = ModelData$mu_opt[i], Nmax = ModelData$Nmax[i], N0 = ModelData$N0[i], Q0 = ModelData$Q0[i])
      sec_temperature <- list(model = "reducedRatkowsky", xmin = ModelData$Tmin[i], b = ModelData$b[i], clade = ModelData$Clade[i])
      my_secondary <- list(temperature = sec_temperature)
      growth <- predict_dynamic_growth(times = env_cond_time[i,],
                                       env_conditions = tibble(time = env_cond_time[i,],
                                                               temperature = env_cond_temp[i,]),
                                       my_primary,
                                       my_secondary)
      sim <- growth$simulation
      ModelData$conc[i] = tail(sim$logN, 1)
    }
    
    ModelData$realCFU = 10^ModelData$conc
    ModelData$CFU_per_serve = ModelData$realCFU*ModelData$serving.size
    log_CFU_per_serving <- log10(ModelData$CFU_per_serve)
    
    # Return the required data frame
    return(data.frame(ModelData, log_CFU_per_serving))
    })
  
  # Generate a histogram for the distribution of cfu per serving in 1000 units HTST milk 
  output$hist <- renderPlot({
    req(data())
    df <- data()
    ggplot(df, aes(x = log_CFU_per_serving)) +
      geom_histogram(binwidth = 0.1, color = "black", fill = "lightblue") +
      labs(title = "Log CFU per Serving", x = "CFU per Serving (log scale)", y = "Number of Servings")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
