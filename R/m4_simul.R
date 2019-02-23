

# Function to simulate the process
#' Simulate the process
#'
#' @param n_geo_dim number of geographical dimensions
#' @param n_eco_dim number of ecological dimensions
#'
#' @return stub
#' @export
#'
m4_simul <- function(n_geo_dim, n_eco_dim) {
  
  # Simulate for a certain time
  
  # Initialization
  # Need to create a multidimensional space
  # What is the unit?
  # A cell in the grid
  # Each cell has a certain density of each species
  # In the end, every species is present in certain cells with certain densities
  # Go asexual first
  
  # Create a world with multiple dimensions
  # Not in matrix format, we want a data frame where each row is a cell
  # So that we can have the densities of multiple species on that row
  geoDimensions <- lapply(seq_len(n_geo_dim), function(i) seq_len(10))
  eco_dimensions <- lapply(seq_len(n_eco_dim), function(i) seq_len(10))
  dimensions <- c(geoDimensions, eco_dimensions)
  #domains <- sapply(dimensions, length)
  #n_cells <- prod(domains)
  n_dim <- length(dimensions) # This variable is not used
  world <- expand.grid(dimensions)
  
  # Initialize with a first species
  # Present in only one cell taken at random
  # Initial density = 10 individuals
  n_cells <- nrow(world)
  starting_cell <- sample(n_cells, 1)
  nSpecies <- 1 # This variable is not used
  species_densities <- as.matrix(rep(0, nrow(world)))
  starting_density <- 10
  species_densities[starting_cell] <- starting_density
  
  # Simulation loop (discrete time)
  # Asexual species
  # One generation life span
  # Dispersal
  # No selection for now
  # Clonal reproduction
  # Local density-dependence in population growth (Ricker model)
  # Death of parents
  # Census
  timestep <- 1 # This variable is not used
  
  # Dispersal
  dispersal_rate <- 0.01
  world_population <- sum(species_densities) # This variable is not used
  
  # Each cell has a density for each species
  # Use these densities as probabilities
  # Record non-empty cell coordinates
  # Create an empty matrix of change in density (deltaDensities)
  # For each non-empty cell
  # Count the number of dispersers from that cell
  # For each disperser
  # Sample a destination
  # The deltaDensity matrix has the same dimensions as the density matrix
  # -1 every time a disperser disperses away
  # +1 every time a destination is chosen
  # After everything has been calculated, add deltaD to D
  # No need to update gradually (goes faster)
  
  delta_densities <- matrix(
    0,
    ncol = ncol(species_densities),
    nrow = nrow(species_densities)
  )
  
  # Count numbers of dispersers across all cells
  n_dispersers <- lapply(seq_len(n_cells), function(i) {
    
    curr.cell <- species_densities[i,]
    
    # Number of dispersers for each species present locally
    n_dispersers <- sapply(
      curr.cell,
      function(n) stats::rbinom(1, n, dispersal_rate)
    )
  })
  
  n_dispersers <- do.call("rbind", n_dispersers)
  if (!all(dim(n_dispersers) == dim(species_densities) &&
           dim(species_densities) == dim(delta_densities))) {
    stop("densities of the matrices do not match")
  }
  
  # Vector of all possible directions
  # Remove the event "no movement"
  all_moves <- expand.grid(lapply(seq_len(n_geo_dim), function(i) moves))
  all_moves <- all_moves[!apply(all_moves, 1, function(x) all(x == 0)),]
  
  # Create a matrix of immigrants
  n_immigrants <- matrix(
    0,
    ncol = ncol(species_densities),
    nrow = nrow(species_densities)
  )
  
  # For each cell, where do the dispersers go?
  # Need a for loop
  for (i in seq_len(nrow(n_dispersers))) {
    for (j in seq_len(ncol(n_dispersers))) {
      if (n_dispersers[i, j] > 0) {
        n_dispersers <- n_dispersers[i, j]
        # For each disperser
        while (n_dispersers > 0) {
          move <- unlist(all_moves[sample(nrow(all_moves), 1), ])
          ii <- i + move[1]
          jj <- j + move[2]
          n_immigrants[ii, jj] <- n_immigrants[ii, jj] + 1
          n_dispersers <- n_dispersers - 1
        }
      }
    }
  }
}
