# 
#' Distance matrix
#'
#' @param population stub
#' @param dimensions stuv
#'
#' @return A distance matrix.
#' @export
make_distance_matrix <- function(population, dimensions) {
  
  population <- population[,dimensions]
  distance_matrix <- dist(population)
  
  return(as.matrix(distance_matrix))
  
}

#' Mate finding
#'
#' @param individual_id stub
#' @param maxmatedistance stub
#' @param geo_distance_matrix stub
#' @param mateslope stub
#' @param mateintercept stub
#' @param population stub
#'
#' @return Numeric indicating mate id.
#' @export
find_mate <- function(individual_id, maxmatedistance, geo_distance_matrix, mateslope, mateintercept, population) {
  
  speciescol <- ncol(population)
  
  # Extract relevant entries from distance matrix
  distance_vec <- geo_distance_matrix[,individual_id]
  
  # Remove the focal individual
  distance_vec <- distance_vec[-individual_id]
  
  # Keep only conspecifics
  conspecifics <- population[-individual_id, speciescol] == population[individual_id, speciescol]
  distance_vec <- distance_vec[conspecifics]
  
  # Remove everything that is too far
  distance_vec <- distance_vec[!distance_vec > maxmatedistance]
  
  # If there are potential mates left
  if (length(distance_vec) > 0) {
    
    # Probability depends on geographical distance (turn this into probability)
    sampling_probs <- mateslope * distance_vec + mateintercept
    sampling_probs[sampling_probs < 0] <- 0
    
    # Sample a mate
    mate_id <- sample(1:length(distance_vec), size = 1, prob = sampling_probs)
    mate_id <- as.numeric(mate_id)
    
    return(mate_id)
    
  } else {
    
    return(NULL)
    
  }
  
}

#' Check mating success
#'
#' @param mom_id 
#' @param dad_id 
#' @param sexslope 
#' @param sexintercept 
#'
#' @return a numeric? int? probability? of mating happening
#' @export
mating_success <- function(mom_id, dad_id, sexslope, sexintercept) {
  
  mom_sex_traits <- t(population[mom_id, sex_dimensions])
  dad_sex_traits <- t(population[dad_id, sex_dimensions])
  
  # Distance in mating space
  sex_distance <- dist(cbind(mom_sex_traits, dad_sex_traits))
  
  # Mating probability (turn this into a probability)
  mating_prob <- sexslope * sex_distance + sexintercept
  if (mating_prob > 1) mating_prob <- 1
  if (mating_prob < 0) mating_prob <- 0
  
  # Sample mating event
  is_mating <- rbinom(1, 1, mating_prob)
  
  return(is_mating)
  
}

#' Produce offspring
#'
#' @param mom_id stub
#' @param dad_id stub
#' @param dispersal_distance stub
#' @param eco_dimensions stub
#' @param sex_dimensions stub
#' @param geo_dimensions stub
#' @param mutation_rate_eco stub
#' @param mutation_rate_sex stub
#' @param mutational_effect_eco stub
#' @param mutational_effect_sex stub
#' @param bounds stuv
#' @param individual_species stub
#'
#' @return Vector of stub.
#' @export
produce_offspring <- function(mom_id, dad_id, dispersal_distance, eco_dimensions, sex_dimensions, geo_dimensions, mutation_rate_eco, mutation_rate_sex, mutational_effect_eco, mutational_effect_sex, bounds, individual_species) {
  
  # Extract parental trait values
  mom_eco_traits <- population[mom_id, eco_dimensions]
  dad_eco_traits <- population[dad_id, eco_dimensions]
  mom_sex_traits <- population[mom_id, sex_dimensions]
  dad_sex_traits <- population[dad_id, sex_dimensions]
  
  # Phenotypic traits of the offspring (regression towards the mean, quantitative genetics assumption)
  off_eco_traits <- (mom_eco_traits + dad_eco_traits) / 2
  off_sex_traits <- (mom_sex_traits + dad_sex_traits) / 2
  
  ##if (any(is.na(c(off_eco_traits, off_eco_traits)))) stop("Here!")
  
  # Sample mutation events
  ntraits <- length(c(eco_dimensions, sex_dimensions))
  is_mutation_eco <- as.logical(rbinom(length(eco_dimensions), 1, mutation_rate_eco))
  is_mutation_sex <- as.logical(rbinom(length(sex_dimensions), 1, mutation_rate_sex))
  
  # Implement mutation
  if (any(is_mutation_eco)) {
    
    # How many dimensions mutate?
    nmutations <- length(which(is_mutation_eco))
    
    # Determine the mutated phenotype
    mutated <- sapply(1:nmutations, function(i) {
      nonmutated <- off_eco_traits[is_mutation_eco]
      truncnorm::rtruncnorm(
        n = 1,
        a = bounds[1],
        b = bounds[2],
        mean = nonmutated,
        sd = mutational_effect_eco
      )
    })
    
    # Update the relevant dimensions for the offspring
    off_eco_traits[is_mutation_eco] <- mutated
    
  }
  
  if (any(is_mutation_sex)) {
    
    # How many dimensions mutate?
    nmutations <- length(which(is_mutation_sex))
    
    # Determine the mutated phenotype
    mutated <- sapply(1:nmutations, function(i) {
      nonmutated <- off_sex_traits[is_mutation_sex]
      truncnorm::rtruncnorm(n = 1, a = bounds[1], b = bounds[2], mean = nonmutated, sd = mutational_effect_sex)
    })
    
    # Update the relevant dimensions for the offspring
    off_sex_traits[is_mutation_sex] <- mutated
    
  }
  
  # Mom's geographical location
  x_mom <- population[mom_id, geo_dimensions][1]
  y_mom <- population[mom_id, geo_dimensions][2]
  
  # Geographical location of the offspring within the bounds
  x_offspring <- truncnorm::rtruncnorm(
    n = 1,
    a = bounds[1],
    b = bounds[2],
    mean = x_mom,
    sd = dispersal_distance
  )
  
  y_offspring <- truncnorm::rtruncnorm(
    n = 1,
    a = bounds[1],
    b = bounds[2],
    mean = y_mom,
    sd = dispersal_distance
  )
  
  ##if (any(is.na(c(x_offspring, y_offspring)))) stop("Here too!")
  
  offspring <- unlist(c(
    off_eco_traits,
    off_sex_traits,
    x_offspring,
    y_offspring,
    individual_species
  ))
  
  return(offspring)
  
}

#' Determines if individual survives
#'
#' @param individual_id stub
#' @param base_survival stub
#' @param eco_dimensions stub
#' @param eco_cutoff stub
#' @param geo_cutoff stub
#' @param eco_distance_matrix stub
#' @param geo_distance_matrix stub
#' @param niche_width stub
#' @param geo_width stub
#' @param resource_peaks stub
#' @param resource_width stub
#' @param max_carrying_capacity stub
#'
#' @return Integer 1 or 0 determining if individual survives (1) or dies (0).
#' @export
survive <- function(individual_id, base_survival, eco_dimensions, eco_cutoff, geo_cutoff, eco_distance_matrix, geo_distance_matrix, niche_width, geo_width, resource_peaks, resource_width, max_carrying_capacity) {
  
  # Who are the competitors?
  is_competitor <- eco_distance_matrix[individual_id,] <= eco_cutoff & geo_distance_matrix[individual_id,] <= geo_cutoff
  competitors <- population[is_competitor,]
  
  # Competition perceived by the focal individual
  competition <- sum(sapply(1:nrow(competitors), function(competitor_id) {
    
    competition_kernel(individual_id, competitor_id, niche_width, geo_width, eco_distance_matrix, geo_distance_matrix)
    
  }))
  
  # Carrying capacity of the focal individual
  carrying_capacity <- get_carrying_capacity(individual_id, population, resource_peaks, resource_width, max_carrying_capacity, eco_dimensions)
  
  # Density dependence factor
  density_dependence <- competition / carrying_capacity
  
  # Survival probability
  survival_prob <- base_survival / density_dependence
  ##print(survival_prob)
  
  if (survival_prob < 0) survival_prob <- 0
  if (survival_prob > 1) survival_prob <- 1
  
  # Does it survive?
  is_survivor <- rbinom(1, 1, survival_prob)
  
  return(is_survivor)
  
}

#' Calculate competition coefficient
#'
#' @param focal_id stub
#' @param competitor_id stub 
#' @param niche_width stub
#' @param geo_width stub
#' @param eco_distance_matrix stub
#' @param geo_distance_matrix stub
#'
#' @return
#' @export numeric with competition coefficient
#'
#' @examples
competition_kernel <- function(focal_id, competitor_id, niche_width, geo_width, eco_distance_matrix, geo_distance_matrix) {
  
  ecological_distance <- eco_distance_matrix[focal_id, competitor_id]
  geographical_distance <- geo_distance_matrix[focal_id, competitor_id]
  
  competition_coeff <- exp(- 0.5 * (ecological_distance^2 / niche_width^2 + geographical_distance^2 / geo_width^2))
  ##print(competition_coeff)
  
  return(competition_coeff)
  
}

# Carrying capacity function
get_carrying_capacity <- function(individual_id, population, resource_peaks, resource_width, max_carrying_capacity, eco_dimensions) {
  
  # Trait of the individual
  individual_traits <- population[individual_id, eco_dimensions]
  
  # Distances to each peak in each dimension
  distances_to_peaks <- individual_traits - resource_peaks
  
  # Take the sum of squares
  ssdistances_to_peaks <- sum(distances_to_peaks^2)
  
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
    truncnorm::rtruncnorm(n = nindiv, a = bounds[1], b = bounds[2], mean = midpoint, sd = sd)
  })
  population <- do.call("cbind", population)
  population <- as.data.frame(population)
  population$species <- 1
  
  return(population)
  
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
update_register <- function(register, population, t) {
  
  # Record all species and append them to the register
  allspecies <- unique(population$species)
  register[[t + 1]] <- allspecies
  return(register)
  
}



# Update the population with new species if speciation has happened
update_speciation <- function(population, speciation_delta, minspecsize) {
  
  # What are all the species?
  allspecies <- unique(population$species)
  
  # For each species...
  for (i in seq_len(length(allspecies))) {
    
    #print(i)
    
    curr.species <- allspecies[i]
    
    # How many individuals?
    n <- length(which(population$species == curr.species))
    
    # If there are enough individuals to look at meaningful clusters...
    if (n >= minspecsize) {
      
      # Check whether the current species has speciated
      species_ids <- check_split(curr.species, population, sex_dimensions, speciation_delta)
      
    } else {
      
      species_ids <- NULL
      
    }
    
    # If there is speciation
    if (length(species_ids) > 0) {
      
      # Replace the species ids with the new ids
      population$species[population$species == curr.species] <- species_ids
      
    }
    
  }
  
  return(population)
  
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
  if (is_speciation) {
    return(as.numeric(paste0(species_id, mod2$cluster)))
  } else {
    return(NULL)
  }
  
}
