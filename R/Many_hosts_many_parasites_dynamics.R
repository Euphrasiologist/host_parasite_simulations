#' Many hosts, many parasites simulation replicates
#' 
#' A simple wrapper function to replicate the simulations of the function manyhost_manyparasite().
#' 
#' @details Model assumes (1) that hosts are fixed in position and perennial, and are unaffected by parasitism.
#'   (2) no effect on parasite of nearby hosts, (3) parasites seed randomly over the field, (4) all parasites survive to end of generation,
#'   (5) parasites do not produce any offspring without a host, and all parasites with a host reproduce and (6) no host independent mortality of parasites
#' @param field.size size of host parasite matrix to play on. Must be square.
#' @param interaction.matrix the selective advantage to a parasite of landing on a host. Supply a matrix of selection coefficients for each parasite on each host. Hosts on rows, parasites on columns! See examples.
#' @param host.number number of hosts to seed the matrix with. Supply a vector of values.
#' @param parasite.number number of parasites to seed the matrix with. Supply a vector of values.
#' @param gens number of generations to run the model for.
#' @param reps how many replicates to run this 
#' @return a list of the parasite dynamics, spatial structure and the host matrix?
#' @import data.table
#' @export
#' @examples
#' 
#' # say I have two host species and two parasite species.
#' int.mat <- matrix(c(1,1,2,2), 2, 2, dimnames = list(c("C.cristatus", "L.corniculatus"), c("E.vigursii", "E.anglica")))
#' 
#' sim <- manyhost_manyparasite_dynamics(field.size = 100^2, 
#'    interaction.matrix = int.mat, 
#'    host.number = c(3000, 3000), 
#'    parasite.number = c(100, 300), 
#'    gens = 100, 
#'    reps = 10)
#' 
#' examples to come...

manyhost_manyparasite_dynamics <- function(field.size, 
                                           interaction.matrix, 
                                           host.number, 
                                           parasite.number, 
                                           gens, 
                                           reps){
  # housekeeping
  if(sqrt(field.size) %% 1 != 0){
    stop("Field size should be square")
  }
  if(!is.matrix(interaction.matrix)){
    stop("Interaction matrix should be a matrix")
  }
  if(dim(interaction.matrix)[1] != length(host.number) | dim(interaction.matrix)[2] != length(parasite.number)){
    stop("Please check that interaction matrix and numbers of each host and parasite are compatible.")
  }
  if(sum(host.number) > field.size | sum(parasite.number) > field.size){
    stop("Too many hosts or parasites for the size of field!")
  }
  if(is.null(dimnames(interaction.matrix))){
    stop("The interaction matrix should have rownames and colnames")
  }
  # add in the replicates to a blank list
  replicates <- list()
  for(i in 1:reps){
    replicates[[i]] <- manyhost_manyparasite(field.size = field.size, 
                                             interaction.matrix = interaction.matrix, 
                                             host.number = host.number, 
                                             parasite.number = parasite.number, 
                                             gens = gens)
    # add replicate number to each iteration
    replicates[[i]]$Reprod.Para.Host$Rep <- as.factor(i)
    replicates[[i]]$Reprod.Para$Rep <- as.factor(i)
    replicates[[i]]$Pop.Size.Para$Rep <- as.factor(i)
  }
  
  Reprod.Para.Host <- lapply(replicates, "[[", 1)
  Reprod.Para <- lapply(replicates, "[[", 2)
  Pop.Size.Para <- lapply(replicates, "[[", 3)
  
  res <- list(Reprod.Para.Host.dyn = rbindlist(Reprod.Para.Host),
              Reprod.Para.dyn = rbindlist(Reprod.Para),
              Pop.Size.Para.dyn = rbindlist(Pop.Size.Para))
  
}