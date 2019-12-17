#' Many hosts, many parasites simulation.
#' 
#' Take a host grid containing a number of hosts with some benefit to the parasite (randomly spaced, but fixed).
#' Then overlay a parasite grid contaning some parasites (randomly spaced).
#' If hosts and parasites co-occur, the parasite will reproduce. Cumulative offspring then randomly overlaid on grid.
#' 
#' 
#' @details Model assumes:
#'   (1) that hosts are fixed in position and perennial, and are unaffected by parasitism.
#'   (2) no effect on parasite of nearby hosts, 
#'   (3) parasites seed randomly over the field, 
#'   (4) all parasites survive to end of generation,
#'   (5) parasites do not produce any offspring without a host, and all parasites with a host reproduce and (6) no host independent mortality of parasites
#' @param field.size size of host parasite matrix to play on. Must be square.
#' @param interaction.matrix the selective advantage to a parasite of landing on a host. Supply a matrix of selection coefficients for each parasite on each host. Hosts on rows, parasites on columns! See examples.
#' @param host.number number of hosts to seed the matrix with. Supply a vector of values.
#' @param parasite.number number of parasites to seed the matrix with. Supply a vector of values.
#' @param gens number of generations to run the model for.
#' @return a list of (1) reproductive output of parasite-host pair per generation
#'                   (2) reproductive output of parasite per generation and
#'                   (3) population size of each parasite in the current generation
#' @import data.table
#' @export
#' @examples
#' 
#' # say I have three host species and three parasite species.
#' int.mat <- matrix(c(1,1,1,2,2,2,3,3,3), 3, 3, dimnames = list(c("C.cristatus", "L.corniculatus", "P.lanceolata"), c("E.vigursii", "E.anglica", "E.micrantha")))
#' 
#' manyhost_manyparasite(field.size = 25, 
#'    interaction.matrix = int.mat, 
#'    host.number = c(8,8,8), 
#'    parasite.number = c(2, 2, 2), 
#'    gens = 100)

manyhost_manyparasite <- function(field.size,
                                  interaction.matrix,
                                  host.number,
                                  parasite.number,
                                  gens) {
  # housekeeping
  if (sqrt(field.size) %% 1 != 0) {
    stop("Field size should be square")
  }
  if (!is.matrix(interaction.matrix)) {
    stop("Interaction matrix should be a matrix")
  }
  if (dim(interaction.matrix)[1] != length(host.number) |
      dim(interaction.matrix)[2] != length(parasite.number)) {
    stop(
      "Please check that interaction matrix and numbers of each host and parasite are compatible."
    )
  }
  if (sum(host.number) > field.size |
      sum(parasite.number) > field.size) {
    stop("Too many hosts or parasites for the size of field!")
  }
  if (is.null(dimnames(interaction.matrix))) {
    stop("The interaction matrix should have rownames and colnames")
  }
  
  fz <- field.size
  # create a (random) sample of host, with a number of hosts defined by user.
  host <- sample(x = c(rep(0, fz - sum(host.number)), rep(1:dim(interaction.matrix)[1], host.number)),
                 size = fz,
                 replace = FALSE)
  # H is the spatial distribution of hosts
  H <- matrix(host, sqrt(fz), sqrt(fz))
  # Hbin is the presence/absence of hosts
  Hbin <- as.matrix((H > 0) + 0)
  
  # parasite matrix
  P <- sample(x = c(rep(0, fz - sum(parasite.number)), rep(1:dim(interaction.matrix)[2], parasite.number)),
              size = fz,
              replace = FALSE)
  P <- matrix(data = P, sqrt(fz), sqrt(fz), byrow = TRUE)
  
  # distribution of parasites at each generation
  mats <- list()
  
  # count of parasite species (population size of parasites)
  # total population sizes of each parasite per generation
  pop.size <- list()
  
  # a breakdown of the number of parasites on each host per generation
  # actually the reproductive output of each parasite on each host (seeds).
  sz <- list()
  
  # reproductive output of each parasite regardless of host.
  pz <- list()
  
  for (j in 2:gens) {
    
    # starting population of parasite
    mats[[1]] <- P
    
        ###################################
    ####        Population sizes           ####
        ###################################
    
    # get pop sizes from generation 1 onwards
    gpop <- as.data.table(table(mats[[j - 1]] * Hbin))
    # rename the table
    setnames(
      x = gpop,
      old = c("V1", "N"),
      new = c("match", "Population Size")
    )
    # match on the names from the interaction matrix and gpop
    matches <- data.table(Parasite = colnames(interaction.matrix),
                          Generation = j - 1, # is j - 1 right here?
                          match = as.character(1:length(colnames(interaction.matrix))))
    fin <- gpop[matches, on = .(match)][, -"match"][, .(Parasite, Generation, `Population Size`)]
    # save as list in pop.size
    pop.size[[j - 1]] <- fin 
    
        ###################################
    ####    Parasite reproductive output   ####
        ###################################
    
    # to get the desired interaction values for the generation j
    # first column is the parasite species second column is the host species
    select <- matrix(c(mats[[j - 1]], H), ncol = 2)
    
    # convert to data table
    select <- as.data.table(select)
    # change names
    setnames(x = select,
             old = c("V1", "V2"),
             new = c("P", "H"))
    # remove zero elements
    mat2 <- select[P != 0][H != 0]
    # access or subset the interaction matrix for the correct value
    get_int <- function(x, y) {
      interaction.matrix[x, y]
    }
    # get the names of the parasites and hosts
    get_names <- function(x, y) {
      cn <- colnames(interaction.matrix)
      rn <- rownames(interaction.matrix)
      paste(rn[x], cn[y], sep = ", ")
    }
    # append reproductive output to mat2 (2,1 because hosts are on rows of interaction matrix)
    mat2[, `:=`(
      Pop.size = apply(
        X = mat2,
        1,
        FUN = function(x)
          get_int(x[2], x[1]) # apply get_int over all rows of mat2
      ),
      Host.parasite = apply(
        X = mat2,
        1,
        FUN = function(x)
          get_names(x[2], x[1]) # apply get_names over all rows of mat2
      ),
      gens = j - 1
    )]
    
    # so mat2 contains the reproductive output of parasite from host-parasite interaction
    # get the total reproductive output for different host-parasite interactions.
    mat3 <- mat2[, .(Reproductive.output = sum(Pop.size)), by = .(Host.parasite, Generation = gens)]
    # save in sz
    sz[[j-1]] <- mat3
    # but we also want total parasite number per generation
    mat4 <-mat3[, .(Host.parasite = gsub(".*, ", "", Host.parasite),
                    Generation,
                    Reproductive.output)][, .(Reproductive.output = sum(Reproductive.output)), by = .(Host.parasite, Generation)]
    # host parasite to factor
    mat4$Host.parasite <- as.factor(mat4$Host.parasite)
    # save this
    pz[[j]] <- mat4
    
    # v. important that order of mat4 species is the same as that of the interaction matrix
    index <- data.table(Host.parasite = colnames(interaction.matrix),
                        INDEX = 1:length(colnames(interaction.matrix)))
    index <- mat4[index, on = .(Host.parasite)]
    # replace NA values with zero. See https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
    replNA <-  function(DT) {
      for (i in names(DT))
        DT[is.na(get(i)), (i) := 0]
    }
    # excecute
    replNA(index)
    
    # numbers of each parasite for the sample in the next generation
    para.int <-c(rep.int(x = index$INDEX,
                         times = as.integer(index$Reproductive.output)),
                 rep(0, fz - ifelse(test = sum(index$Reproductive.output) > fz,
                                    yes = fz,
                                    no = sum(index$Reproductive.output))))
    # actually sample and save in mats[[]]
    mats[[j]] <- matrix(data = sample(x = para.int,
                                      size = fz,
                                      replace = FALSE),
                        nrow = sqrt(fz),
                        ncol = sqrt(fz),
                        byrow = TRUE)
  } # end loop
  
  # collate results
  res <- rbindlist(sz)
  colnames(res) <- c("Host-Parasite", "Generation", "Seed") # suppose it is technically seed in the current generation
  res2 <- rbindlist(pz)
  colnames(res2) <- c("Parasite", "Generation", "Seed")
  res3 <- rbindlist(pop.size)
  
  output <- list(Reprod.Para.Host = res,
                 Reprod.Para = res2,
                 Pop.Size.Para = res3)
  
  class(output) <- c("HoPaSim")
  output
}
