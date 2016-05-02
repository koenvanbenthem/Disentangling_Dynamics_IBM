########################
####### Type of program
## Written in R
## Object Oriented
#######################

#######################
# Main object "Leprechaun": an individual of our species of interest. 
# The Leprechaun is not very choosy, and mates completely random.
##########################################


## Random seed
#set.seed(12) -- to be used when testing

LeprechaunSimul <- function(DIR= "Data",datafile="dataset",
                          SurvivalSelection1=0.5,
                          fertilitySelection1=0.1,
                          basePhi=0.75,
                          energySelection=0.1,
                          SurvivalPenaltyForRepro=-0.01,
                          MeanResources=15000,
                          SDResources=10,
                          ResourcesTrend=-200,
                          MeanBirthSize=10,
                          MeanGrowth=1.245,
                          GrowthEnergyConversion=0.05,
                          SexualMaturity=1,
                          PlasticityBirthSize=1,
                          PlasticityReproduction=0,
                          MaternalEffect=0.1,
                          MeanRepro=0.5,
                          StudyLength=30,
                          InitialPopSize=200,
                          GeneticMatrices="GeneticMatrices") {
  ## Setting the genotype phenotype map from disk
  ## Loading the genotype-phenotype map in
  load(file=GeneticMatrices)
  ## number of loci targeting Z
  nbLociZ <- dim(gvaluesZ)[3] 
  ## number of existing alleles per locus targeting Z
  nbAllelesZ <- dim(gvaluesZ)[2] 
  
  ## Loads in the class and function definitions 
  ## this is done per run to adapt to changing sizes of genomes
  source("Class.R",local=TRUE) 
  filename <- paste(DIR,"/",datafile,".csv",sep="")
  
  ## Setting counters
  CID <- as.integer(1) # Overall ID
  YR <- 0  # Year
  
  ## Creating an initial population with 100 individuals
  pop <- c(new("Leprechaun"))
  for(i in 2:InitialPopSize){
  	pop <- c(pop,new("Leprechaun"))
  }
  
  ## List of living individuals [their indices], this will save time later,
  ## because dead individuals are not looped over
  ## ALIVE is a bookkeeping variable: it contains the indices of the 
  ## alive individuals. This is used throughout the simulation.
  ALIVE <- 1:length(pop) 
  
  # Write column names to file
  cat("t\tID\tz\tbvs\tC\ts\tARS\tage\tp1\tp2\tphi",
      file=filename,append=FALSE)
  
  #####################################################################
  ## The start of time 
  for (YR in 1:StudyLength) {
     print(YR)

    ## after year 10, food availability increases 
    ExpectedResources <- ifelse(test=  YR > 10, 
                                yes=MeanResources+(YR-10) * ResourcesTrend,
                                no=MeanResources)
    ## resources for year YR
    Resources <- abs(round(rnorm(n=1,
                                 mean=ExpectedResources,
                                 sd=ExpectedResources / SDResources),
                     digits=0))
    
  	## Competition for resources
  	## A randomly assigned hunterquality
    HunterQualities <- runif(n=length(ALIVE),min=0.3,max=1) 
    if(length(HunterQualities) == 1){ HunterQualities <- 1 }
    
    ## the relative hunterquality determines the share of food that an 
    ## individual obtains
    camams <- round(Resources * HunterQualities / sum(HunterQualities),
                    digits=0) 
    
    ## The amount of food is written to the population.
    pop[ALIVE] <- lapply(1:length(ALIVE), function(x) 
                                          Food(pop[[ALIVE[x]]],camams[x])) 
    
    ## Survival
    DEAD<-c()
    ## the surv function is applied on the population, causing the indices
    ## of the 'dead' individuals to be stored in DEAD.
  	pop[ALIVE] <- lapply(pop[ALIVE],Surv) 
  	ALIVE <- ALIVE[!(ALIVE %in% DEAD)]
  	
    if (length(ALIVE) == 0) {
      for(i in DEAD){
        cat("\n",YR,"\t",pop[[i]]@ID,"\t",NA,"\t",
            pop[[i]]@bvs,"\t",
            pop[[i]]@energy,"\t",pop[[i]]@sex,"\t",NA,"\t",
            NA,"\t",pop[[i]]@pID[1],"\t",pop[[i]]@pID[2],"\t",
            0,file=filename,append=TRUE)
      }
      break
    }
    
  	## Age+1 and growth
  	pop[ALIVE] <- lapply(pop[ALIVE],Age)
  	pop[ALIVE] <- lapply(pop[ALIVE],Grow)
  	
  	## Reproduction
  	## Part dedicated to retrieving the indices of all living males and
  	## of all living females
  	## Determine which individuals are males -- logicaly all other
  	## individuals should be females. However, this includes dead ones. 
  	## Simply a list of T,T,F,F,T,F,F,....
  	males <- lapply(pop,Sex) == "M" 
  	females <- which(!males) # Get the indices of the non-males
  	males <- which(males) # Get the indices of the males
  	females <- intersect(females,ALIVE) # Indices of the living(!) females
  	males <- intersect(males,ALIVE) # Indices of the living males
  
  	
  	## Part dedicated to breeding..
  	## We take a female based approach: we determine for each female
  	from <- CID
  	pop[females] <- lapply(pop[females],FUN=Num_off)
  	for (i in females) {
  		Noffs <- pop[[i]]@ARS
  		if (Noffs > 0 & length(males > 0)) {
  			## Determine the father by random sampling
  			fat <- sample(males,1)
  			for (j in 1:Noffs) {
  				# Create the offspring
  				pop <- c(pop,new("Leprechaun",parent1=i,parent2=fat))
  			}
  		}		
  	}
    
  	if (from != CID) {
  	  ALIVE <- c(ALIVE,(from):(CID - 1))
  	}
  	## Write to file
    for (i in ALIVE) {
      cat("\n",YR,"\t",pop[[i]]@ID,"\t",round(pop[[i]]@size,4),"\t",pop[[i]]@bvs,
          "\t",pop[[i]]@energy,
          "\t",pop[[i]]@sex,"\t",pop[[i]]@ARS,"\t",
          pop[[i]]@age,"\t",pop[[i]]@pID[1],"\t",pop[[i]]@pID[2],"\t",
          1,file=filename,append=TRUE)
    }
    for (i in DEAD) {
      cat("\n",YR,"\t",pop[[i]]@ID,"\t",NA,"\t",pop[[i]]@bvs,
          "\t",pop[[i]]@energy,
          "\t",pop[[i]]@sex,"\t",NA,"\t",NA,"\t",
          pop[[i]]@pID[1],"\t",pop[[i]]@pID[2],"\t",0,file=filename,
          append=TRUE)
    }
  }
   return(filename)
}

## End of LeprechaunSimul
