# Initialization

library(truncnorm)

# Events
# An individual is a list


population <- rbind(c(1,1), c(2,0), c(1,1))

# Distance matrix
make_distance_matrix <- function(population, dimensions) {
  
  population <- population[,dimensions]
  distance_matrix <- dist(population)
  
  return(as.matrix(distance_matrix))
  
}

# Mate finding
find_mate <- function(individual_id, max_distance, distance_matrix, a, b) {
  
  # Extract relevant entries from distance matrix
  distance_vec <- distance_matrix[,individual_id]
  
  # Remove the focal individual
  distance_vec <- distance_vec[-individual_id]
  
  # Remove everything that is too far
  distance_vec <- distance_vec[!distance_vec > max_distance]
  
  # Probability depends on geographical distance (turn this into probability)
  sampling_probs <- a * distance_vec + b
  sampling_probs[sampling_probs < 0] <- 0
  
  # Sample a mate
  mate_id <- sample(names(distance_vec), prob = sampling_probs)
  
  return(mate_id)
}


# Does mating happens?
mating_success <- function(mom_id, dad_id, a, b) {
  
  mom_sex_traits <- population[mom_id, sex_dimensions]
  dad_sex_traits <- population[dad_id, sex_dimensions]
  
  # Distance in mating space
  sex_distance <- abs(mom_sex_traits - dad_sex_traits)
  
  # Mating probability (turn this into a probability)
  mating_prob <- a * sex_distance + b
  if(mating_prob > 1) mating_prob <- 1
  if(mating_prob < 0) mating_prob <- 0
  
  # Sample mating event
  is_mating <- rbinom(1, 1, mating_prob)
  
  return(is_mating)
  
}

# Function to produce offspring
produce_offspring <- function(mom_id, dad_id, dispersal_distance, eco_dimensions, sex_dimensions, geo_dimensions, mutation_rate_eco, mutation_rate_sex, mutational_effect_eco, mutational_effect_sex) {
  
  # Extract parental trait values
  mom_eco_traits <- population[mom_id, eco_dimensions]
  dad_eco_traits <- population[dad_id, eco_dimensions]
  mom_sex_traits <- population[mom_id, sex_dimensions]
  dad_sex_traits <- population[dad_id, sex_dimensions]
  
  # Phenotypic traits of the offspring (regression towards the mean, quantitative genetics assumption)
  off_eco_traits <- (mom_eco_traits + dad_eco_traits) / 2
  off_sex_traits <- (mom_eco_traits + dad_eco_traits) / 2
  
  # Sample mutation events
  ntraits <- length(c(eco_dimensions, sex_dimensions))
  is_mutation_eco <- rbinom(length(eco_dimensions), 1, mutation_rate_eco)
  is_mutation_sex <- rbinom(length(sex_dimensions), 1, mutation_rate_sex)
  
  # Implement mutation
  if(any(is_mutation_eco)) {

    nmutations <- length(which(is_mutation_eco))
    mutations <- sapply(1:nmutations, function(i) rnorm(1, 0, mutational_effect_eco))
    off_eco_traits[is_mutation_eco] <- off_eco_traits[is_mutation_eco] + mutations
    
  }
  
  if(any(is_mutation_sex)) {

    nmutations <- length(which(is_mutation_sex))
    mutations <- sapply(nmutations, function(i) rnorm(1, 0, mutational_effect_sex))
    off_sex_traits[is_mutation_sex] <- off_sex_traits[is_mutation_sex] + mutations
    
  }
  
  # Mom's geographical location
  x_mom <- population[mom_di, geo_dimensions][1]
  y_mom <- population[mom_di, geo_dimensions][2]
  
  # Geographical location of the offspring
  distance_from_mom <- abs(rnorm(1, 0, dispersal_distance))
  angle_from_mom <- runif(1, 0, 2 * pi)
  x_offspring <- x_mom + distance_from_mom * cos(angle_from_mom)
  y_offspring <- y_mom + distance_from_mom * sin(angle_from_mom)
  
  offspring <- c(off_eco_traits, off_sex_traits, x_offspring, y_offspring)
  
  return(offspring)
  
}

# Survival of an individual
survive <- function(individual_id, base_survival, eco_dimensions, eco_cutoff, geo_cutoff, eco_distance_matrix, geo_distance_matrix) {
  
  # Who are the competitors?
  is_competitor <- eco_distance_matrix[focal_id,] <= eco_cutoff & geo_distance_matrix[focal_id,] <= geo_cutoff
  competitors <- population[is_competitor,]
  
  # Competition perceived by the focal individual
  competition <- lapply(competitors, function(competitor_id) {
    
    competition_kernel(individual_id, competitor_id, max_ecol_distance, niche_width, geo_width, eco_distance_matrix, geo_distance_matrix)
    
  })
  
  # Carrying capacity of the focal individual
  carrying_capacity <- get_carrying_capacity(individual_id, population, resource_peaks, resource_width, max_carrying_capacity, eco_dimensions)
  
  # Density dependence factor
  density_dependence <- competition / carrying_capacity
  
  # Survival probability
  survival_prob <- base_survival / density_dependence
  
  # Does it survive?
  is_survivor <- rbinom(1, 1, survival_prob)
  
  return(is_survivor)
  
}

# Competition function
competition_kernel <- function(focal_id, competitor_id, max_ecol_distance, niche_width, geo_width, eco_distance_matrix, geo_distance_matrix) {
  
  ecological_distance <- eco_distance_matrix[focal_id, competitor_id]
  geographical_distance <- geo_distance_matrix[focal_id, competitor_id]
  
  competition_coeff <- exp(- 0.5 * (ecological_distance^2 / niche_width^2 + geographical_distance^2 / geo_width^2))
  
  return(competition_coeff)
  
}

# Carrying capacity function
get_carrying_capacity <- function(individual_id, population, resource_peaks, resource_width, max_carrying_capacity, eco_dimensions) {
  
  # Trait of the individual
  individual_traits <- population[individual_id, eco_dimensions]
  
  # Distances to peak in each dimension
  distances_to_peaks <- individual_traits - resource_peaks
  
  # Take the sum of squares
  ssdistances_to_peaks <- sum((distances_to_peaks)^2)
  
  # Carrying capacity equation
  # Note: multidimensional normal distribution over ecological space, constant across geographical space
  carrying_capacity <- max_carrying_capacity * exp(- 0.5 * ssdistances_to_peaks / resource_width^2)
  
  return(carrying_capacity)
  
}

# Function to initialize the population data frame
initialize_population <- function(necol, nspatial, nmating, nindiv, bounds = c(0,10), sd = 1) {
  
  # Center on the midpoint between the boundaries
  midpoint <- mean(bounds)
  
  # Sample initial values for all individuals and all dimensions
  dimensions <- seq_len(necol + nspatial + nmating)
  population <- lapply(dimensions, function(dim) {
    rtruncnorm(n = nindiv, a = bounds[1], b = bounds[2], mean = midpoint, sd = sd)
  })
  population <- do.call("cbind", population)
  
  return(population)
  
}

#### Simulation ####

# Initialize the population with a function
population <- initialize_population()

# Initialize the resource mapping function over ecological and geographical space
# Possibly by reading an input file

# Check boundary conditions in reproduction, putain

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



# Function to take a snapshot of the population
take_selfie <- function(population, t) {
  
  outfile <- paste0("selfie_", t, ".csv")
  
  # Write the population matrix to CSV
  write.csv(population, outfile)
  
}

#### Post processing ####

#phylogeny <- phylogenize(register) # Newick

#### Speciation related functions ####

# Update the phylogeny
update_register <- function(register, population) {
  
  # Record all species and append them to the register
  allspecies <- levels(population$species)
  register <- append(register, allspecies)
  return(register)
  
}

# Update the population with new species if speciation has happened
update_speciation <- function(population) {
  
  # What are all the species?
  allspecies <- levels(population$species)
  
  # For each species...
  for(i in seq_len(length(allspecies))) {
    
    curr.species <- allspecies[i]
    
    # Check whether the current species has speciated
    species_ids <- check_split(curr.species, population, sex_dimensions, speciation_delta)
    
    # If there is speciation
    if(species_ids != FALSE) {
      
      # Replace the species ids with the new ids
      population$species[population$species == curr.species,] <- species_ids
      
    }
    
  }
  
}

# Function to check if speciation has happened for a given species
check_split <- function(species_id, population, sex_dimensions, speciation_delta) {
  
  # Subset the focal species
  population <- population[population$species == species_id,]
  
  # Extract the mating traits
  all_sex_traits <- population[,sex_dimensions]
  
  # Has the species split? K-means analysis
  mod1 <- kmeans(all_sex_traits, 1)
  mod2 <- kmeans(all_sex_traits, 2)
  
  fit1 <- mod1$betweenss / mod1$totss
  fit2 <- mod2$betweenss / mod2$totss
  
  is_speciation <- fit2 - fit1 > speciation_delta
  
  # Return the species ID of all individuals
  if(is_speciation) {
    return(as.numeric(paste0(species_id, mod2$cluster)))
  } else {
    return(FALSE)
  }
 
}