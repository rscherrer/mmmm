

#' Stub function to temporarily wrap around loop. needs variables
#'
#' @param necol stub
#' @param nspatial stub
#' @param nmating stub
#' @param nindiv stub
#' @param bounds stub
#' @param sd stub
#' @param tmax stub
#' @param maxmatedistance stub
#' @param mateslope stub
#' @param mateintercept stub
#' @param sexslope stub
#' @param sextintercept stub
#' @param dispersal_distance stub
#' @param mutation_rate_eco stub
#' @param mutation_rate_sex stub
#' @param mutational_effect_eco stub
#' @param mutational_effect_sex stub
#' @param base_survival stub
#' @param eco_cutoff stub
#' @param geo_cutoff stub
#' @param niche_width stub
#' @param geo_width stub
#' @param resource_width stub
#' @param max_carrying_capacity stub
#' @param specialization_delta stub
#' @param minspecsize stub
#' @param selfietime stub
#'
#' @return stub. simulation.
#' @export
#'
#' @author Pedro Neves
do_things <- function(necol,
                      nspatial,
                      nmating,
                      nindiv,
                      bounds,
                      sd,
                      tmax,
                      maxmatedistance,
                      mateslope,
                      mateintercept,
                      sexslope,
                      sextintercept,
                      dispersal_distance,
                      mutation_rate_eco,
                      mutation_rate_sex,
                      mutational_effect_eco,
                      mutational_effect_sex,
                      base_survival,
                      eco_cutoff,
                      geo_cutoff,
                      niche_width,
                      geo_width,
                      resource_width,
                      max_carrying_capacity,
                      specialization_delta,
                      minspecsize,
                      selfietime) {
  # Initialize the population with a function
  population <-
    initialize_population(necol, nspatial, nmating, nindiv, bounds, sd)
  
  # Initialize the register with a single species
  register <- list()
  register[[1]] <- 1
  
  eco_dimensions <- 1:necol
  geo_dimensions <- (necol + 1):(necol + nspatial)
  sex_dimensions <-
    (necol + nspatial + 1):(necol + nspatial + nmating)
  
  # Initialize the resource mapping function
  # over ecological and geographical space
  # Possibly by reading an input file
  
  # Check boundary conditions in reproduction
  
  t <- 0
  
  popsizes <- 0
  
  # Simulation loop
  while (t <= tmax) {
    print(t)
    
    # Calculate distance matrices
    eco_distance_matrix <-
      make_distance_matrix(population, eco_dimensions)
    geo_distance_matrix <-
      make_distance_matrix(population, geo_dimensions)
    
    # Mating -- create a population of offspring
    # For each individual...
    offspring <-
      lapply(1:nrow(population), function(individual_id) {
        # Species of the focal individual
        individual_species <-
          population[individual_id, ncol(population)]
        
        # Find a mate
        mate_id <-
          find_mate(
            individual_id,
            maxmatedistance,
            geo_distance_matrix,
            mateslope,
            mateintercept,
            population
          )
        
        
        # If there is an encounter...
        if (!is.null(mate_id)) {
          # Is mating successful?
          is_mating <-
            mating_success(individual_id, mate_id, sexslope, sexintercept)
          
        } else {
          # Otherwise no mating
          is_mating <- 0
          
        }
        
        # If yes, produce offspring
        if (is_mating) {
          offspring <-
            produce_offspring(
              individual_id,
              mate_id,
              dispersal_distance,
              eco_dimensions,
              sex_dimensions,
              geo_dimensions,
              mutation_rate_eco,
              mutation_rate_sex,
              mutational_effect_eco,
              mutational_effect_sex,
              bounds,
              individual_species
            )
          
          # Exception tracking
          if (any(is.na(offspring)))
            stop("NA found!")
          
        } else {
          offspring <- NULL
          
        }
        
        return(offspring)
        
      })
    
    offspring <- do.call("rbind", offspring)
    
    # Create a population of survivors
    # For each individual...
    survivors_id <- as.logical(
      sapply(1:nrow(population), function(individual_id) {
        
        # Does the focal individual survive?
        does_survive <- survive(
          individual_id,
          base_survival,
          eco_dimensions,
          eco_cutoff,
          geo_cutoff,
          eco_distance_matrix,
          geo_distance_matrix,
          niche_width,
          geo_width,
          resource_peaks,
          resource_width,
          max_carrying_capacity
        )
        return(does_survive)
      }))
    
    # Extract the survivors from the old population
    survivors <- population[survivors_id, ]
    
    # Plug together offspring and survivors
    colnames(offspring) <- colnames(survivors)
    population <- rbind(survivors, offspring)
    
    # Update the population matrix based on speciation
    population <- update_speciation(
      population,
      speciation_delta,
      minspecsize
      )
    
    # Update the phylogeny
    register <- update_register(register, population, t)
    # Is it time to take a selfie?
    is_selfietime <- t %% selfietime == 0
    # Take a selfie
    if (is_selfietime) {
      take_selfie(population, t)
    }
    popsizes[t + 1] <- nrow(population)
    # Advance time
    t <- t + 1
  }
}