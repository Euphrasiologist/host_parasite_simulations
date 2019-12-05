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


manyhost_manyparasite <- function(field.size = 10^2, 
                                  interaction.matrix = matrix(c(1,1,1,2,2,2,3,3,3), 3, 3, dimnames = list(c("Host 1", "Host 2", "Host 3"), c("Parasite 1", "Parasite 2", "Parasite 3"))), 
                                  host.number = c(10,10,10), 
                                  parasite.number = c(5,5,5), 
                                  gens = 100){
  if(dim(interaction.matrix)[1] != length(host.number | dim(interaction.matrix)[2] != length(parasite.number))){
    stop("Please check that interaction matrix and numbers of each host and parasite are compatible.")
  }
  if(sum(host.number) > field.size | sum(parasite.number) > field.size){
    stop("Too many hosts or parasites for the size of field!")
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
  # a breakdown of the number of parasites on each host per generation
  sz <- list()
  # total population sizes of each parasite per generation
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
    
    # so mat2 contains the reproductive output of parasite from host-parasite interaction
    # get the total population size for different host-parasite interactions.
    # save the host parasite interaction numbers per generation
    mat3 <- mat2[, .(Pop.size.total = sum(Pop.size)), by = .(Host.parasite, gens)]
    # save in sz
    sz[[j]] <- mat3
    # but we also want total parasite number per generation
    mat4 <- mat3[, .(Host.parasite = gsub(".*, ", "", Host.parasite), gens, Pop.size.total)][, .(Total.parasite.number = sum(Pop.size.total)), by = .(Host.parasite, gens)][order(Host.parasite)]
    mat4$Host.parasite <- as.factor(mat4$Host.parasite)
    # save this
    pz[[j]] <- mat4
    
    # early on in iterations a parasite may have died out, and therefore this rep() command below will throw an error
    if(length(colnames(interaction.matrix)) != length(levels(mat4$Host.parasite))){
      # find the difference in levels
      diffs <- setdiff(colnames(interaction.matrix),levels(mat4$Host.parasite))
      add <- data.table(Host.parasite = diffs, gens = j, Total.parasite.number = 0)
      mat4 <- rbind(mat4, add)
    }
    
    # what is this?
    mats[[j]] <- matrix(data = sample(x = c(rep(x = 1:dim(interaction.matrix)[2], times = as.numeric(mat4$Total.parasite.number)), 
                                            rep(0, fz- ifelse(test = sum(mat4$Total.parasite.number) >fz, yes = fz, no = sum(mat4$Total.parasite.number)))), 
                                      size = fz, replace = FALSE), 
                        nrow = sqrt(fz), ncol = sqrt(fz), byrow = TRUE) 
  }
  res<-rbindlist(sz)
  colnames(res) <- c("Host-Parasite", "Generation", "Population Size")
  res2<-rbindlist(pz) 
  colnames(res2) <- c("Parasite", "Generation", "Population Size")
  
  return(list(res, res2))
}
