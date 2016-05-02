## Class definition
library("HMMpa")

setClass(Class="Leprechaun",
         representation=representation(ID = "integer",
                                       pID = "integer",
                                       age = "integer",
                                       Birth = "integer",
                                       alive = "logical",
                                       size = "numeric",
                                       energy = "integer",
                                       sex = "character",
                                       DNAZ = "matrix",
                                       DNAR = "matrix",
                                       bvs = "numeric",
                                       ARS = "integer"))

## Definition of the basic methods 
## (for printing to the screen and initialisation)
setMethod("show","Leprechaun",
          function(object){
            cat(object@ID,"\t",object@size,"\t",object@energy,"\t",
                object@age,"\t",object@sex,"\t(",object@pID[1],",",
                object@pID[2],")\t",object@Birth,"\t",object@alive,"\n",
                sep="")
          }
)

## Parent1 is the mother
setMethod("initialize","Leprechaun",function(.Object,parent1,parent2) {
  ## Set empty genome for body size at birth
  .Object@DNAZ <- matrix(NA,nrow=2,ncol=nbLociZ) 
  ## If no parent #1, assign this apart of the genes at random
  if (missing(parent1)) {
    parent1<-NA
    .Object@DNAZ[1,] <- floor(runif(nbLociZ,min=1,max=nbAllelesZ+1))
  ## else, use the genome of the parent to decide what the genome of the
  ## offspring looks like
  } else { 
    # weight1 <- pop[[parent1]]@size
    .Object@DNAZ[1,] <- pop[[parent1]]@DNAZ[cbind(floor(runif(n=nbLociZ, 
                                                              min=1,
                                                              max=3)),
                                                  1:nbAllelesZ)]
  }
  if (missing(parent2)) {
    parent2 <- NA
    .Object@DNAZ[2,] <- floor(runif(n=nbLociZ,min=1,max=nbAllelesZ + 1))
  } else {
    .Object@DNAZ[2,] <- pop[[parent2]]@DNAZ[cbind(floor(runif(n=nbAllelesZ,
                                                              min=1,
                                                              max=3)),
                                                  1:nbAllelesZ)]
  }
  .Object@age <- as.integer(0) # age to zero
  .Object@ID <- CID # assign ID
  .Object@pID <- c(as.integer(parent1),as.integer(parent2)) # parental IDs
  .Object@Birth <- as.integer(YR) # assign birth year
  .Object@alive <- TRUE
  BreedingValueSize <- 0
  ## Take the mean of genetic values
  for (Locus in 1:nbLociZ) {
    BreedingValueSize <- BreedingValueSize + 
                         (gvaluesZ[ .Object@DNAZ[1,Locus],
                                   .Object@DNAZ[2,Locus], Locus] / 
                          nbLociZ)
  }
  .Object@bvs <- BreedingValueSize
  ## Assign birth size based on genotype
  size <- MeanBirthSize + BreedingValueSize 
  
  if (!is.na(parent1)) {
    size <- size + MaternalEffect * pop[[parent1]]@size # add maternal effect
  }
  # sd plasticity birth size
  .Object@size <- abs(rnorm(n=1,mean=size,sd=PlasticityBirthSize)) 
  

  .Object@energy <- as.integer(0)
  .Object@ARS <- as.integer(0) # Annual reproductive success
  
  if (runif(1) > 0.5) {
    .Object@sex <- 'F'
  }  else {
    .Object@sex <- 'M'
  }
  
  CID <<- as.integer(CID + 1)
  
  return(.Object)
})

#######################################################################
## Definition of demographic functions

## Mortality function (bathtub curve)
bathtub <- function(age,shift=4,BaseMortality=0.1) {
  p <- BaseMortality * exp(-(age - shift) / 4) + 
       (-1 + exp((age - shift) * log(2) / (30 - shift)))
  p[p > 1] <- 1
  return(p)
}

## Incorporate the effect of size to the survival function
sizeSurvival <- function(age,size,energy,ARS) {
  page <- 1 - bathtub(age)
  ## because size does not prevent animals of maximal age to die out
  if(page < 1 & age > 0) {
    Philogit <- log(basePhi / (1 - basePhi)) + 
      SurvivalSelection1 * (size - 15) + 
      ARS * SurvivalPenaltyForRepro
    psize <- exp(Philogit) / (1 + exp(Philogit))
    p <- page * psize
    if (energy < 40) {
      p <- p * (energy / 40) ^ (1 / 2)
    }
  }
  if (age == 0) p <- page
  return(p)
}

## Applying the bathtub in a surival function
setGeneric("Surv",function(Object) standardGeneric("Surv"))
setMethod("Surv","Leprechaun",function(Object) {
  p <- sizeSurvival(Object@age,Object@size,Object@energy,
                                         Object@ARS)
  
  if (runif(1) > p) {
    Object@alive <- FALSE
    DEAD <<- c(DEAD,Object@ID)
  }
  return(Object)
  }
)
## Adds 1 to the age
setGeneric("Age",function(Object) standardGeneric("Age"))
setMethod("Age","Leprechaun",function(Object) {
  Object@age <- as.integer(Object@age + 1)
  return(Object)
  }
)

# Growth // Sizes change
setGeneric("Grow",function(Object) standardGeneric("Grow"))
setMethod("Grow","Leprechaun",function(Object) {
  ## Maximal growth converge to 1 with increasing age
  MeanGrowthAge <- (MeanGrowth + Object@age - 1) / Object@age 
  SdGrowthAge <- (MeanGrowthAge - 1) / 2
  ## Gives the energy specific expected growth
  GrowthC <- 1 + (MeanGrowthAge - 1) * 
             (2 / (1 + exp(-GrowthEnergyConversion * (Object@energy))) - 1)
  ## Variation in growth; the animal can always shrink
  Object@size <- Object@size * rnorm(1,mean=GrowthC,sd=SdGrowthAge)
  return(Object)
  }
)

## Retrieving the ID of an individual
setGeneric("IDretrieve",function(Object) standardGeneric("IDretrieve"))
setMethod("IDretrieve","Leprechaun",function(Object) return(Object@ID))

## Retrieving the size of an individual
setGeneric("Size",function(Object) standardGeneric("Size"))
setMethod("Size","Leprechaun",function(Object) return(Object@size))

## Retrieving the sex of an individual
setGeneric("Sex",function(Object) standardGeneric("Sex"))
setMethod("Sex","Leprechaun",function(Object) return(Object@sex))

## Calculating the number of offspring for a females
setGeneric("Num_off",function(Object) standardGeneric("Num_off"))
setMethod("Num_off","Leprechaun",function(Object) { 
  logitV <- Object@age - SexualMaturity
  p <- 1 / (1 + exp(-logitV))
  
  if(rbinom(1,size=1,prob=p) == 1) {
    mu <- exp(log(MeanRepro) + 
                  fertilitySelection1 * (Object@size-15) + 
                  energySelection * ((Object@energy) ^ (1 / 3)) / 15
              )
   
    repro <- rpois(1,mu) + 1 
    Object@ARS <- as.integer(repro)
  } else{
    Object@ARS <- 0L
  } 
  return(Object)
  }
)

## Function for energy attributions 
setGeneric("Food",function(Object,FoodCaptured) standardGeneric("Food"))
setMethod("Food","Leprechaun",function(Object,FoodCaptured) {
  Object@energy <- as.integer(FoodCaptured)
  return(Object)
  }
)
