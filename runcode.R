##### File for run code ####
# @Neves-P: Run code must not be stored with functions.
# Added this to Rbuildignore.

#### Parameters for a function? ####
# @Neves-P: From m4_simul. These should be parameters of the function?
nGeoDim <- 2 # number of geographical dimensions
nEcoDim <- 2 # number of ecological dimensions

##### Experiment setup ####
# @Neves-P: from timoandraph. This is one specific case to run?
# These should be parameters to be plugged in function?


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


#### Plot output ####

# @Neves-P: From timoandraph. These should be a seperate function?
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




