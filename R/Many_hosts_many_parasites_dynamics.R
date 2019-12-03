#' Many hosts, many parasites simulation.
#' 
#' Take a host grid containing a number of hosts with some benefit to the parasite (randomly spaced, but fixed).
#' Then overlay a parasite grid contaning some parasites (randomly spaced).
#' If hosts and parasites co-occur, the parasite will reproduce. Cumulative offspring then randomly overlaid on grid.
#' 
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
#'    parasite.number = c(1, 1), 
#'    gens = 100, 
#'    reps = 10)
#' 
#' examples to come...

manyhost_manyparasite_dynamics <- function(field.size = 100^2, 
                                           interaction.matrix = matrix(c(1,1,2,2), 2, 2, dimnames = list(c("Host 1", "Host 2"), c("Parasite 1", "Parasite 2"))), 
                                           host.number = c(2000, 2000), 
                                           parasite.number = c(1, 1), 
                                           gens = 100, 
                                           reps = 5){
  # add in the replicates to a blank list
  replicates <- list()
  for(i in 1:reps){
    replicates[[i]] <- manyhost_manyparasite(field.size = field.size, 
                                             interaction.matrix = interaction.matrix, 
                                             host.number = host.number, 
                                             parasite.number = parasite.number, 
                                             gens = gens)
  # what was this bit for again?
    if(dim(replicates[[i]])[1] == 0){
      replicates[[i]] <- data.table(Parasite = 0,
                                    Generation = 0,
                                    `Population size` = 0,
                                    rep = i)
    } else 
      replicates[[i]]$rep <- i
  }
  # bind all the lists together!
  res <- rbindlist(l = replicates, fill = TRUE)
  res[Parasite > 0]
}