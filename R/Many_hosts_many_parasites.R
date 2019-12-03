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
#' @return a list of the parasite dynamics, spatial structure and the host matrix?
#' @import data.table
#' @export
#' @examples
#' 
#' # say I have three host species and three parasite species.
#' int.mat <- matrix(c(1,1,1,2,2,2,3,3,3), 3, 3, dimnames = list(c("C.cristatus", "L.corniculatus", "P.lanceolata"), c("E.vigursii", "E.anglica", "E.micrantha")))
#' 
#' manyhost_manyparasite(field.size = 100^2, 
#'    interaction.matrix = int.mat, 
#'    host.number = c(3000, 3000, 3000), 
#'    parasite.number = c(1, 1, 1), 
#'    gens = 100)


manyhost_manyparasite <- function(field.size = 100^2, 
                                  interaction.matrix = matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, dimnames = list(c("Host 1", "Host 2", "Host 3"), c("Parasite 1", "Parasite 2", "Parasite 3"))), 
                                  host.number = c(3000, 3000, 3000), 
                                  parasite.number = c(1, 1, 1), 
                                  gens = 100){
  if(dim(interaction.matrix)[1] != length(host.number | dim(interaction.matrix)[2] != length(parasite.number))){
    stop("Please check that interaction matrix and numbers of each host and parasite are compatible.")
  }
  fz <- field.size
  # create a sample of host, with a number of hosts each with s
  host <- sample(x = c(rep(0,fz-sum(host.number)), rep(1:dim(interaction.matrix)[1], host.number)), size = fz, replace = FALSE)
  # spatial distribution of hosts
  H <- matrix(host, sqrt(fz), sqrt(fz))
  
  # parasite matrix
  P <- sample(x = c(rep(0,fz-sum(parasite.number)), rep(1:dim(interaction.matrix)[2], parasite.number)), size = fz, replace = FALSE)
  P <- matrix(data = P, sqrt(fz), sqrt(fz), byrow = TRUE)
  
  # fill the distribution of parasites at each generation
  mats <- list()
  # a list of the size of the parasite for next generation
  sz <- list()
  # realised population sizes
  pz <- list()
  # table with breakdown of all euphrasia-host interactions
  
  for(j in 2:gens){
    
    # starting population of Euphrasia
    mats[[1]] <- P
    
    # to get the desired interaction values
    select <- matrix(c(mats[[j-1]], H), ncol = 2)
    # we need to be able to identify each host parasite combination...
    Pop.size <- vector(length = dim(select)[1])
    Host.parasite <- vector(length = dim(select)[1])
    
    for(i in 1:dim(select)[1]){
      if(length(interaction.matrix[select[i,1], select[i,2]]) == 0){
        Pop.size[i] <- 0
      } else
        Pop.size[i]<-interaction.matrix[select[i,1], select[i,2]]
      
      Host.parasite[i] <- paste(unlist(dimnames(interaction.matrix[select[i,1], select[i,2], drop = FALSE])), collapse = ", ")
    }
    # bind the results (vec and vec2 from loop just above)
    mat2 <- cbind(select, Pop.size, Host.parasite)
    
    # so here we have a character matrix, which we probably want to convert to a data frame
    mat2 <- as.data.table(mat2)
    mat2 <- mat2[,c("Pop.size", "Host.parasite")]
    mat2 <- mat2[Pop.size > 0]
    mat2$Pop.size <- as.numeric(mat2$Pop.size)
    mat2$gens <- as.factor(j)
    
    mats[[j]] <- matrix(data = sample(x = c(rep(1:dim(interaction.matrix)[2], sum(mat2$Pop.size)), 
                                            rep(0, fz- ifelse(test = sum(mat2$Pop.size) >fz, yes = fz, no = sum(mat2$Pop.size)))), 
                                      size = fz, replace = FALSE), 
                        nrow = sqrt(fz), ncol = sqrt(fz), byrow = TRUE) 
    
    sz[[j]] <- mat2
    pz[[j]] <- as.data.table(table(mats[[j]]))[, rep := j]
  }
  res<-rbindlist(pz)
  colnames(res) <- c("Parasite", "Population size", "Generation")
  list(res[Parasite > 0], sz)
  
  # TODO: would be nice to output parasite pop size on each host... how?
}
