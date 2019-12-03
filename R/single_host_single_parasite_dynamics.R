#' Single host, single parasite simulation for many iterations.
#' 
#' Take a host grid containing a number of hosts with some benefit to the parasite (randomly spaced, but fixed).
#' Then overlay a parasite grid contaning some parasites (randomly spaced).
#' If hosts and parasites co-occur, the parasite will reproduce. Cumulative offspring then randomly overlaid on grid.
#' 
#' This function is a wrapper for onehost_oneparasite() and simply replicates the number of simulations nicely.
#' 
#' 
#' @details Model assumes (1) that hosts are fixed in position and perennial, and are unaffected by parasitism.
#'   (2) no effect on parasite of nearby hosts, (3) parasites seed randomly over the field, (4) all parasites survive to end of generation,
#'   (5) parasites do not produce any offspring without a host, and all parasites with a host reproduce and (6) no host independent mortality of parasites
#' @seealso onehost_oneparasite_dynamics to replicate many times.
#' @param field.size size of host parasite matrix to play on. Must be square.
#' @param s the selective advantage to a parasite of landing on a host.
#' @param host.number number of hosts to seed the matrix with.
#' @param parasite.number number of parasites to seed the matrix with.
#' @param gens number of generations to run the model for.
#' @param reps number of replicates of each simulation.
#' 
#' @importFrom data.table rbindlist
#' @export
#' @examples
#' 
#' # should print out simulation from default parameters
#' onehost_oneparasite() 


onehost_oneparasite_dynamics <- function(field.size = 10^2, s = 3, host.number = 35, parasite.number = 20, gens = 100, reps = 100){
  # add in the replicates to a blank list
  replicates <- list()
  for(i in 1:reps){
    replicates[i] <- onehost_oneparasite(field.size = field.size, s = s, host.number = host.number, parasite.number = parasite.number, gens = gens)
  }
  # bind all the lists together!
  res <- rbindlist(l = replicates)
  res$rep <- as.factor(rep(1:reps, each = gens))
  
  res
  
}
