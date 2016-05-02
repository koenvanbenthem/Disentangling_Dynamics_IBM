
FunctionGvalues <- function(nbLoci=10,nbAlleles=10,dominance=0.5,
                            SDeffects=1,SDalleles=1) {
 
  ## Initialising a matrix that will contain the genotypic effects on trait
  gvalues <- array(data=NA,dim=c(nbAlleles,nbAlleles,nbLoci),
                   dimnames=list(paste("A",1: nbAlleles,sep=""),
                                 paste("A",1: nbAlleles,sep=""),
                                 paste("L",1:nbLoci,sep=""))) 
  
  for (L in 1:nbLoci) {
    ## Setting the effects for the homozygotes [all loci]
    ## alter the locus importance in a realistic way (many small-effect loci, 
    ## few major loci)
    effect <- abs(rnorm(n=1,mean=0,sd=SDeffects))
    diag(gvalues[,,L]) <- 2 * rnorm(n=dim(gvalues)[1],mean=0,
                                 sd=effect * SDalleles)
    ## Setting the effects for the heterozygotes
    ## loop for off-diagonal = heterozygotes (additive and dominance effects)
    for (A in 1:(nbAlleles - 1)) {
      for (D in (A + 1):nbAlleles) {
        d <- dominance * runif(n=1,min=-0.5,max=0.5)
        ## mean of additive effects + dominance, over diagonal
        gvalues[A,D,L] <- (0.5 - d) * gvalues[A,A,L] + (0.5 + d) * 
                          gvalues[D,D,L] 
        ## the same below diagonal    
        gvalues[D,A,L] <- (0.5 - d) * gvalues[A,A,L] + (0.5 + d) * 
                          gvalues[D,D,L] 
      }
    }
  }
  return(gvalues)
}

######################################################################
## Genetic determinism
# dominance <- 1 # for additive effects only, must be 0
# overdominance <- 0 # non-null values generate overdominance
# 
# nbLoci <- 10 # number of loci controling the trait phenotype
# nbAlleles <- 10 # number of existing alleles per loci
# 
# SDZ <- 1 # standard deviation of locus effect on size


## Low heritability
## Genetic determination Z: body size
approxVa <- 1
while (approxVa >= 1) {
  gvaluesZlow <- FunctionGvalues(nbLoci=10, nbAlleles=10, dominance=1,
                                 SDeffects=1, SDalleles=1)
  approxVa <- var(gvaluesZlow)
}

# mean(apply(gvaluesZ,1,FUN = function(x){apply(X = x,MARGIN = 1,FUN = var)}))

## Store the values
gvaluesZ <- gvaluesZlow
save(gvaluesZ,file="GeneticMatrices_lowh2")


## High heritability
## Genetic determination Z: body size
approxVa <- 5
while (approxVa <= 5 | approxVa > 10) {
  gvaluesZhigh <- FunctionGvalues(nbLoci=10, nbAlleles=10, dominance=1,
                                  SDeffects=1,SDalleles=1)
  approxVa <- var(gvaluesZhigh)
}

va <- vector(length=100)
for (i in 1:100) {
  va[i] <- var(FunctionGvalues(nbLoci=10, nbAlleles=10, dominance=1,
               SDeffects=1, SDalleles=1))
}

##hist(va)
##mean(va)
##min(va)

# mean(apply(gvaluesZ,1,FUN = function(x){apply(X = x,MARGIN = 1,FUN = var)}))

## Store the values
gvaluesZ <- gvaluesZhigh
save(gvaluesZ, file="GeneticMatrices_highh2")
