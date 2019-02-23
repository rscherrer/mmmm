

# Function to simulate the process
#' Simulate the process
#'
#' @param nGeoDim number of geographical dimensions
#' @param nEcoDim number of ecological dimensions
#'
#' @return stub
#' @export
#'
m4_simul <- function(nGeoDim, nEcoDim) {
  
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
  geoDimensions <- lapply(seq_len(nGeoDim), function(i) seq_len(10))
  ecoDimensions <- lapply(seq_len(nEcoDim), function(i) seq_len(10))
  dimensions <- c(geoDimensions, ecoDimensions)
  #domains <- sapply(dimensions, length)
  #nCells <- prod(domains)
  nDim <- length(dimensions) # This variable is not used
  world <- expand.grid(dimensions)
  
  # Initialize with a first species
  # Present in only one cell taken at random
  # Initial density = 10 individuals
  nCells <- nrow(world)
  startingCell <- sample(nCells, 1)
  nSpecies <- 1 # This variable is not used
  speciesDensities <- as.matrix(rep(0, nrow(world)))
  startingDensity <- 10
  speciesDensities[startingCell] <- startingDensity
  
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
  dispersalRate <- 0.01
  worldPopulation <- sum(speciesDensities) # This variable is not used
  
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
  
  deltaDensities <- matrix(0, ncol = ncol(speciesDensities), nrow = nrow(speciesDensities))
  
  # Count numbers of dispersers across all cells
  nDispersers <- lapply(seq_len(nCells), function(i) {
    
    curr.cell <- speciesDensities[i,]
    
    # Number of dispersers for each species present locally
    nDispersers <- sapply(curr.cell, function(n) rbinom(1, n, dispersalRate))
  })
  
  nDispersers <- do.call("rbind", nDispersers)
  if (!all(dim(nDispersers) == dim(speciesDensities) &&
           dim(speciesDensities) == dim(deltaDensities))) {
    stop("densities of the matrices do not match")
  }
  
  # Vector of all possible directions
  # Remove the event "no movement"
  allMoves <- expand.grid(lapply(seq_len(nGeoDim), function(i) moves))
  allMoves <- allMoves[!apply(allMoves, 1, function(x) all(x == 0)),]
  
  # Create a matrix of immigrants
  nImmigrants <- matrix(0, ncol = ncol(speciesDensities), nrow = nrow(speciesDensities))
  
  # For each cell, where do the dispersers go?
  # Need a for loop
  for (i in seq_len(nrow(nDispersers))) {
    for (j in seq_len(ncol(nDispersers))) {
      if (nDispersers[i,j] > 0) {
        ndispersers <- nDispersers[i,j]
        
        # For each disperser
        while (ndispersers > 0) {
          move <- unlist(allMoves[sample(nrow(allMoves), 1),])
          ii <- i + move[1]
          jj <- j + move[2]
          nImmigrants[ii, jj] <- nImmigrants[ii, jj] + 1
          ndispersers <- ndispersers - 1
        }
        
      }
    }
  }
  
  
  
}
