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

manyhost_manyparasite <- function(field.size, 
                                  interaction.matrix, 
                                  host.number, 
                                  parasite.number, 
                                  gens){
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
  
  fz <- field.size
  # create a sample of host, with a number of hosts each with s
  host <- sample(x = c(rep(0,fz-sum(host.number)), rep(1:dim(interaction.matrix)[1], host.number)), size = fz, replace = FALSE)
  # spatial distribution of hosts
  H <- matrix(host, sqrt(fz), sqrt(fz))
  # and presence/absence of hosts
  Hbin <- as.matrix((H > 0) + 0)
  
  # parasite matrix
  P <- sample(x = c(rep(0,fz-sum(parasite.number)), rep(1:dim(interaction.matrix)[2], parasite.number)), size = fz, replace = FALSE)
  P <- matrix(data = P, sqrt(fz), sqrt(fz), byrow = TRUE)
  
  # these need to be changed.
  
  # distribution of parasites at each generation
  mats <- list()
  
  # count of parasite species (population size of parasites)
  # total population sizes of each parasite per generation
  pop.size <- list()

  # a breakdown of the number of parasites on each host per generation
  # actually the reproductive output of each parasite on each host.
  sz <- list()

  # reproductive output of each parasite regardless of host.
  pz <- list()
  
  for(j in 2:gens){
    
    # starting population of parasite
    mats[[1]] <- P
    
        ###################################
    ####        Population sizes           ####
        ###################################
    
    # get pop sizes from generation 1 onwards
    # multiply by binary host presence matrix??
    gpop <- as.data.table(table(mats[[j-1]]*Hbin))
    # rename the table
    setnames(x = gpop, old = c("V1", "N"), new = c("match", "Population Size"))
    # match on the names from the 
    matches <- data.table(Parasite = colnames(interaction.matrix), Generation = j-1, match = as.character(1:length(colnames(interaction.matrix))))
    fin <- gpop[matches, on = .(match)][,-"match"]
    pop.size[[j-1]] <- fin

        ###################################
    ####    Parasite reproductive output   ####
        ###################################
    
    # to get the desired interaction values for the generation j
    select <- matrix(c(mats[[j-1]], H), ncol = 2)
    # to get the parasites for generation j+1
    Pop.size <- vector(length = dim(select)[1])
    # names of the host-parasite pair for fitness in j+1
    Host.parasite <- vector(length = dim(select)[1])
    
    # need a condition here to deal with when there is only a single parasite
    # otherwise an error that you can't select a 1 column matrix with method
    # below.
    
    # One parasite.
    if(ncol(interaction.matrix) == 1){
      # for switch argument later
      type <- 1 
      
      for(i in 1:dim(select)[1]){
        # # if column 1, row i is zero, popsize is zero. OR if the second column (host presence) is zero.
        if(select[i,1] == 0 | length(interaction.matrix[select[i,2]]) == 0){
          Pop.size[i] <- 0
        } else
          # population size is the value in the interaction matrix. 
          Pop.size[i]<- interaction.matrix[select[i,2]]
        # not totally happy about this. But I think it works.
        Host.parasite[i] <- paste(dimnames(interaction.matrix)[[2]], unlist(dimnames(interaction.matrix))[select[i,2]], sep = ", ")
      }
      
    } else
      
      # One host.
      if(nrow(interaction.matrix) == 1){
        # for switch argument later
        type <- 2
        
        for(i in 1:dim(select)[1]){
          # the next generation parasite contribution is zero if there is no host - parasite overlap.
          if(select[i,1] == 0 | select[i,2] == 0){
            Pop.size[i] <- 0
          } else
            # population size in j+1 is the value in the interaction matrix when host - parasite overlap.
            Pop.size[i]<- interaction.matrix[select[i,1]]
          # not totally happy about this. But I think it works.
          Host.parasite[i] <- paste(dimnames(interaction.matrix)[[2]][select[i,1]], dimnames(interaction.matrix)[[1]], sep = ", ")
        }
    } else
      
      # little worried about whether this bit makes sense?
      # More than one parasite and one host. Can we generalise this to
      if(ncol(interaction.matrix) > 1 & nrow(interaction.matrix) > 1){
        # for switch argument later
        type <- 3
        
      for(i in 1:dim(select)[1]){
        # if there is no host parasite co-incidence
        if(select[i,1] == 0 | select[i,2] == 0){
          Pop.size[i] <- 0
          Host.parasite[i] <- as.character("")
        } else
          # if there are more parasites than hosts or more hosts than parasites
          if(ncol(interaction.matrix) > nrow(interaction.matrix) | ncol(interaction.matrix) < nrow(interaction.matrix)){
            # population size is the value in the interaction matrix
            Pop.size[i]<-interaction.matrix[select[i,2], select[i,1]]
            # name of the host parasite pair
            Host.parasite[i] <- paste(unlist(dimnames(interaction.matrix[select[i,2], select[i,1], drop = FALSE])), collapse = ", ")
          } else
            # if number of hosts is equal to number of parasites
            if(ncol(interaction.matrix) == nrow(interaction.matrix)){
              # population size is the value in the interaction matrix
              Pop.size[i]<-interaction.matrix[select[i,1], select[i,2]]
              # name of the host parasite pair
              Host.parasite[i] <- paste(unlist(dimnames(interaction.matrix[select[i,1], select[i,2], drop = FALSE])), collapse = ", ")
            }
          
      }
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
    # get the total reproductive output for different host-parasite interactions.
    mat3 <- mat2[, .(Pop.size.total = sum(Pop.size)), by = .(Host.parasite, gens)]
    # save in sz
    sz[[j]] <- mat3
    # but we also want total parasite number per generation
    mat4 <- mat3[, .(Host.parasite = switch(type,
                                            gsub(", .*", "", Host.parasite), # for one parasite multiple hosts
                                            gsub(", .*", "", Host.parasite), # for one host multiple parasites, below is multiple of both.
                                            gsub(".*, ", "", Host.parasite)), gens, Pop.size.total)][, .(Total.parasite.number = sum(Pop.size.total)), by = .(Host.parasite, gens)][order(Host.parasite)]
    mat4$Host.parasite <- as.factor(mat4$Host.parasite)
    # save this
    pz[[j]] <- mat4
    
    # early on in iterations a parasite may have died out, and therefore this rep() command below will throw an error
    if(length(colnames(interaction.matrix)) != length(levels(mat4$Host.parasite))){
      diffs <- setdiff(colnames(interaction.matrix),levels(mat4$Host.parasite))
      add <- data.table(Host.parasite = diffs, gens = j, Total.parasite.number = 0)
      mat4 <- rbind(mat4, add)
      }
    
    # seed the next matrix with the numbers of hosts and parasites from the iteration before.
    # having problems with rep()
    # replicate parasites the number of times they appear in the matrix
    # seems to be peculiar behaviour with one-host, many parasites with low fitnesses.
    
    para.int <- c(rep.int(x = 1:dim(interaction.matrix)[2], times = as.integer(mat4$Total.parasite.number)), 
                  rep(0, fz - ifelse(test = sum(mat4$Total.parasite.number) > fz, yes = fz, no = sum(mat4$Total.parasite.number))))
        
    mats[[j]] <- matrix(data = sample(x = para.int, 
                                      size = fz, 
                                      replace = FALSE), 
                        nrow = sqrt(fz), 
                        ncol = sqrt(fz), 
                        byrow = TRUE) 
  }
  res<-rbindlist(sz)
  colnames(res) <- c("Host-Parasite", "Generation", "Population Size")
  res2<-rbindlist(pz) 
  colnames(res2) <- c("Parasite", "Generation", "Population Size")
  res3 <- rbindlist(pop.size)
  
  output <- list(Reprod.Para.Host = res, 
                 Reprod.Para = res2,
                 Pop.Size.Para = res3)
  
  class(output) <- c("HoPaSim")
  output
}
