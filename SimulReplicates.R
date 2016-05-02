
## Loading genotype, to be read in later.
## to be run only once!
## source("Genetics.R") # Creates a genotype-phenotype map and saves it to

LeprechaunSimulReplicates <- function(DIR= "Data",n_replicates=100,
                                    GeneticMatrices="GeneticMatrices",
                                    PlasticityBirthSize,
                                    SurvivalSelection1,
                                    fertilitySelection1,
                                    MeanRepro) {
  source("main.R")
  for (i in 1:n_replicates) { # Replicates over one set of parameters
    Datafile_Replicate<-paste("dataset",i,sep="")
    LeprechaunSimul(DIR=DIR,
                    datafile=Datafile_Replicate,
                    StudyLength=50,
                    InitialPopSize=200,
                    GeneticMatrices=GeneticMatrices,
                    MeanRepro=MeanRepro,
                    SurvivalSelection1=SurvivalSelection1,
                    fertilitySelection1=fertilitySelection1,
                    PlasticityBirthSize=PlasticityBirthSize)
  }
}

########################################################################
## Run scenarios
########################################################################

## Low h2, no selection
LeprechaunSimulReplicates(DIR="hs/",n_replicates=100,
                          GeneticMatrices="GeneticMatrices_lowh2", 
                          MeanRepro=0.2,
                          PlasticityBirthSize=1,
                          fertilitySelection1=0,
                          SurvivalSelection1=0)

## Low h2, selection via survival and reproduction
LeprechaunSimulReplicates(DIR="hS/",n_replicates=100,
                          GeneticMatrices="GeneticMatrices_lowh2", 
                          MeanRepro=0.2,
                          PlasticityBirthSize=1,
                          fertilitySelection1=0.04,
                          SurvivalSelection1=0.2)

## High h2, no selection
LeprechaunSimulReplicates(DIR="Hs/",n_replicates=100,
                          GeneticMatrices="GeneticMatrices_highh2", 
                          MeanRepro=0.2,
                          PlasticityBirthSize=0.5,
                          fertilitySelection1=0,
                          SurvivalSelection1=0)

## High h2, selection via survival and reproduction
LeprechaunSimulReplicates(DIR="HS/",n_replicates=100,
                          GeneticMatrices="GeneticMatrices_highh2", 
                          MeanRepro=0.2,
                          PlasticityBirthSize=0.5,
                          fertilitySelection1=0.04,
                          SurvivalSelection1=0.2)
