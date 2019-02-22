

#set.seed(42)


#### Simulation ####

# Parameters
necol <- 2 # number of ecological dimensions
nspatial <- 2 # number of geographical dimensions
nmating <- 2 # number of sexual dimensions
nindiv <- 10 # initial pop size
bounds <- c(0,10) # bounds of phenotype and spatial space
sd <- 1 # standard deviation of the initial pop in geomorphospace
dispersal_distance <- 1
maxmatedistance <- 1 # max geographical distance of a potential mate
mateslope <- -0.5 # has to be negative, you are more likely to mate with someone at a shorter distance
mateintercept <- 1 # basal mate encounter probability
sexslope <- 0.5 # slope of the relation between the distance in sex space and mating probability
sexintercept <- 0 # basal mating probability
mutation_rate_eco <- 0.01
mutation_rate_sex <- 0.01
mutational_effect_eco <- 0.01 # standard dev of the distribution from where mutations are sampled
mutational_effect_sex <- 0.01
base_survival <- 0.5 # from one gen to the next
eco_cutoff <- 2 # cut-off beyond which there is no competition
geo_cutoff <- 2
niche_width <- 1 # determines the intensity of competition
geo_width <- 1
resource_peaks <- rep(5, necol) # coordinates of the peak, works for multivariate Gaussian for now (single peak)
resource_width <- 3 # width of the resource distribution across ecological space
max_carrying_capacity <- 2.2 # highest point of the resource distribution
speciation_delta <- 0.6 # how much better does a KMeans model with 2 clusters needs to be in order to conclude speciation happened?
selfietime <- 100 # take a screenshot of the population and save it to csv every X generations
minspecsize <- 10 # minimum population size needed for a species to split
tmax <- 200 # time

# Initialize the population with a function
population <- initialize_population(necol, nspatial, nmating, nindiv, bounds, sd)

# Initialize the register with a single species
register <- list()
register[[1]] <- 1

eco_dimensions <- 1:necol
geo_dimensions <- (necol + 1):(necol + nspatial)
sex_dimensions <- (necol + nspatial + 1):(necol + nspatial + nmating)

# Initialize the resource mapping function over ecological and geographical space
# Possibly by reading an input file

# Check boundary conditions in reproduction, putain

t <- 0

popsizes <- 0

# Simulation loop
while(t <= tmax) {
  
  print(t)
  
  # Calculate distance matrices
  eco_distance_matrix <- make_distance_matrix(population, eco_dimensions)
  geo_distance_matrix <- make_distance_matrix(population, geo_dimensions)
  
  # Mating -- create a population of offspring
  # For each individual...
  offspring <- lapply(1:nrow(population), function(individual_id) {
    
    ##print(individual_id)
    
    # Species of the focal individual
    individual_species <- population[individual_id, ncol(population)]
    
    # Find a mate
    mate_id <- find_mate(individual_id, maxmatedistance, geo_distance_matrix, mateslope, mateintercept, population)
    #print(mate_id)
    
    
    # If there is an encounter...
    if(!is.null(mate_id)) {
      
      # Is mating successful? 
      is_mating <- mating_success(individual_id, mate_id, sexslope, sexintercept)
      
    } else {
      
      # Otherwise no mating
      is_mating <- 0
      
    }
    
    # If yes, produce offspring
    if(is_mating) {
      
      offspring <- produce_offspring(individual_id, mate_id, dispersal_distance, eco_dimensions, sex_dimensions, geo_dimensions, mutation_rate_eco, mutation_rate_sex, mutational_effect_eco, mutational_effect_sex, bounds, individual_species)
      
      # Exception tracking
      if(any(is.na(offspring))) stop("NA found!")
    
    } else {
      
      offspring <- NULL
        
    }
      
    return(offspring)
    
  })
  
  offspring <- do.call("rbind", offspring)
  
  # Create a population of survivors
  # For each individual...
  survivors_id <- as.logical(sapply(1:nrow(population), function(individual_id) {
    
    # Does the focal individual survive?
    does_survive <- survive(individual_id, base_survival, eco_dimensions, eco_cutoff, geo_cutoff, eco_distance_matrix, geo_distance_matrix, niche_width, geo_width, resource_peaks, resource_width, max_carrying_capacity)
    
    return(does_survive)
    
  }))
  
  # Extract the survivors from the old population
  survivors <- population[survivors_id,]
  
  # Plug together offspring and survivors
  colnames(offspring) <- colnames(survivors)
  population <- rbind(survivors, offspring)
  
  # Update the population matrix based on speciation
  population <- update_speciation(population, speciation_delta, minspecsize)
  
  # Update the phylogeny
  register <- update_register(register, population, t)
  
  # Is it time to take a selfie?
  is_selfietime <- t %% selfietime == 0
  
  # Take a selfie
  if(is_selfietime) {
    take_selfie(population, t)
  }
  
  popsizes[t + 1] <- nrow(population)
  
  # Advance time
  t <- t + 1
  
}

# Check the results
par(mfrow = c(2,2))
plot(popsizes, pch = 16)
plot(population[,eco_dimensions], col = "darkgreen", pch = 16)
plot(population[,sex_dimensions], xlim = bounds, ylim = bounds, col = "darkred", pch = 16)
plot(population[,geo_dimensions], xlim = bounds, ylim = bounds, col = "darkblue", pch = 16)
par(mfrow = c(1,1))

population
table(population$species)
register
