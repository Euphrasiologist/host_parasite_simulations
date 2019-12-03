#' Single host, single parasite simulation. Take a host grid containing a number of hosts with some benefit to the parasite (randomly spaced, but fixed). Then overlay a parasite grid contaning some parasites (randomly spaced). If hosts and parasites co-occur, the parasite will reproduce. Cumulative offspring then randomly overlaid on grid.
#' 
#' @details Model assumes 1 - that hosts are fixed in position and perennial, and are unaffected by parasitism. 2 - no effect on parasite of nearby hosts, 3 - parasites seed randomly over the field, 4 - all parasites survive to end of generation, 5 - parasites do not produce any offspring without a host, and all parasites with a host reproduce and 6 - no host independent mortality of parasites
#' @param field.size size of host parasite matrix to play on. Must be square.
#' @param s the selective advantage to a parasite of landing on a host.
#' @param host.number number of hosts to seed the matrix with.
#' @param parasite.number number of parasites to seed the matrix with.
#' @param gens number of generations to run the model for.
#' @return a list of the parasite dynamics, spatial structure and the host matrix
#' @export
#' @examples
#' # should print out simulation from default parameters
#' onehost_oneparasite() 

onehost_oneparasite <- function(field.size = 5^2, s = 1, host.number = 5, parasite.number = 5, gens = 100){
  fz <- field.size
  # is whole number function
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  if(!is.wholenumber(sqrt(fz))){
    stop("Field size should be a square!")
  }
  # create a sample of host
  host <- sample(x = c(rep(0,fz-host.number), rep(1,host.number)), size = fz, replace = FALSE)
  # spatial distribution of hosts, multiplied by selective advantage
  H <- matrix(host, sqrt(fz), sqrt(fz))*s
  
  # parasite matrix
  P <- sample(x = c(rep(1, parasite.number), rep(0,fz-parasite.number)), size = fz, replace = FALSE)
  P <- matrix(data = P, sqrt(fz), sqrt(fz), byrow = TRUE)
  
  # i == two empty spaces to fill, one for distribution of parasites
  mats <- list()
  # one for
  sz <- vector(length = gens)
  
  for(i in 2:gens){
    # starting population of Euphrasia
    mats[[1]]<- P
    sz[1] <- sum(P)
    # multiply Euphrasia matrix by host matrix
    #mats[[i]] <- mats[[i-1]]*h1
    # total number of Euphrasia in t+1
    sz[i] <- sum(mats[[i-1]]*H)
    # generate random matrix of Euphrasia with n(t+1) from the sum of elements in mat[[i-1]] (all 1's)
    mats[[i]] <- matrix(data = sample(x = c(rep(1, sz[i]), rep(0, fz- ifelse(test = sz[i] >fz, yes = fz, no = sz[i]))), size = fz, replace = FALSE), 
                        nrow = sqrt(fz), ncol = sqrt(fz), byrow = TRUE)*(H/s)
    sz[i] <- sum(mats[[i]])
  }
  
  output <- list(P.dynamics = data.frame(generation = 1:gens, P.size = sz), 
                 Spatial.structure = mats, 
                 Host.matrix = H)
  output
}
