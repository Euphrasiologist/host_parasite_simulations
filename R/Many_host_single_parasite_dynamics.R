#' Many hosts, single parasite simulation for many iterations.
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
#' @param s the selective advantage to a parasite of landing on a host. Supply a vector of values.
#' @param host.number number of hosts to seed the matrix with. Supply a vector of values.
#' @param parasite.number number of parasites to seed the matrix with.
#' @param gens number of generations to run the model for.
#' @return a list of the parasite dynamics, spatial structure and the host matrix?
#' @importFrom data.table rbindlist
#' @export
#' @examples
#' # need ggplot
#' library(ggplot2)
#' sim <- manyhost_oneparasite_dynamics()
#' ggplot(sim, aes(x = generation, y = P.size, group = rep)) + geom_line()

manyhost_oneparasite_dynamics <- function(field.size = 100^2, s = c(0,0,0,0,0,5,6,7,8,9), host.number = c(rep(10000/10, 10)), parasite.number = 1, gens = 100, reps = 100){
  # add in the replicates to a blank list
  replicates <- list()
  for(i in 1:reps){
    replicates[i] <- manyhost_oneparasite(field.size = field.size, s = s, host.number = host.number, parasite.number = parasite.number, gens = gens)
  }
  # bind all the lists together
  res <- rbindlist(l = replicates)
  res$rep <- as.factor(rep(1:reps, each = gens))
  res
}
