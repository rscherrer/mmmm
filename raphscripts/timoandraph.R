# Initialization

library(truncnorm)



# Run the simulation
out <- mmmm_simul()



source("raphscripts/functions.R")

#### Simulation ####

# Parameters
necol <- 2
nspatial <- 2
nmating <- 2
nindiv <- 10
bounds <- c(0,10)
sd <- 1
dispersal_distance <- 1

# Initialize the population with a function
population <- initialize_population(necol, nspatial, nmating, nindiv, bounds, sd)

eco_dimensions <- 1:necol
geo_dimensions <- (necol + 1):(necol + nspatial)
sex_dimensions <- (necol + nspatial + 1):(necol + nspatial + nmating)

# Initialize the resource mapping function over ecological and geographical space
# Possibly by reading an input file

# Check boundary conditions in reproduction, putain

t <- 0

# Simulation loop
while(t <= tmax) {
  
  # Calculate distance matrices
  eco_distance_matrix <- make_distance_matrix(population, eco_dimensions)
  geo_distance_matrix <- make_distance_matrix(population, geo_dimensions)
  
  # Mating -- create a population of offspring
  # For each individual...
  offspring <- lapply(1:nrow(population), function(individual_id) {
    
    # Find a mate
    mate_id <- find_mate(individual_id)
    
    # Is mating successful? 
    is_mating <- mating_success(individual_id, mate_id)
    
    # If yes, produce offspring
    if(is_mating) offpspring <- produce_offspring(individual_id, mate_id)
    
    return(offspring)
    
  })
  
  # Create a population of survivors
  # For each individual...
  survivors_id <- lapply(1:nrow(population), function(individual_id) {
    
    # Does the focal individual survive?
    does_survive <- survive(individual_id)
    
    return(does_survive)
    
  })
  
  # Extract the survivors from the old population
  survivors <- population[survivors_id,]
  
  # Plug together offspring and survivors
  population <- rbind(survivors, offspring)
  
  # Update the population matrix based on speciation
  population <- update_speciation(population)
  
  # Update the phylogeny
  register <- update_register(register, population)
  
  # Is it time to takea selfie?
  is_selfietime <- t %% selfietime == 0
  
  # Take a selfie
  if(is_selfietime) {
    take_selfie(population, t)
  }
  
  # Advance time
  t <- t + 1
  
}

