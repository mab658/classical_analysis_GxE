# Define a custom function ammiModel to fit Additive Main Effect and Multiplicative Interaction Model
# for a trait of interest from the input data

# Define a function to fit AMMI model for any trait of interest
ammiModel <- function(traitName){
  
  geAmmi <- gxeAmmi(TD = setATD, trait = traitName,nPC=2)
  summOutput <- list(Anova=geAmmi$anova,PCA=geAmmi$importance)
  
  #  extract genotypic scores and mean from AMMI2 model
  genScores <- geAmmi$genoScores
  genMean <- geAmmi$genoMean
  
  genMeanScores <- cbind(genMean,genScores) # merge the two files
  
  # Calculate the AMMI stability values
  genStab <- as.data.frame(genMeanScores)
  colnames(genStab)[1] <- "Mean"
  genStab <- tibble::rownames_to_column(genStab, "Genotype")
  
  # extract out PCs SS
  SSIPCA1 <- geAmmi$anova$`Sum Sq`[4] 
  SSIPCA2 <- geAmmi$anova$`Sum Sq`[5]
  
  genStab$ASV <- sqrt(((SSIPCA1/SSIPCA2)*genStab$PC1)^2 + genStab$PC2^2)
  genStab <- arrange(genStab,desc(Mean)) %>%
    mutate(Mean_rank=1:nrow(genStab))
  
  genStab <- arrange(genStab,ASV) %>%
    mutate(ASV_rank=1:nrow(genStab))
  
  
  # calculate genotype selection index
  genStab$GSI <- genStab$Mean_rank + genStab$ASV_rank 
  write.csv(genStab,file="output/AMMI_stability.csv",row.names=FALSE)
  
  #  extract environmental scores and  mean from geAmmi object
  envScores <- geAmmi$envScores
  envMean <- geAmmi$envMean
  
  # merge the two files
  envMeanScores <- cbind(envMean,envScores)
  
  # Calculate the AMMI stability values
  envMeanScores <- as.data.frame(envMeanScores)
  colnames(envMeanScores)[1] <- "Mean"
  envMeanScores <- as.data.frame(envMeanScores)
  envMeanScores <- tibble::rownames_to_column(envMeanScores, "Env")
  
  
  # Create AMMI2 plot with symetric scaling
  #jpeg(file="AMMI2_biplot.jpeg")
  AMMI2_biplot <- plot(geAmmi, plotType = "AMMI2", scale = 0.5, plotConvHull = TRUE,  
                       sizeGeno = 2.8,sizeEnv = 2.5, title = "")
  #dev.off()
  
  ammiOutput <- list(anovaPCA=summOutput, StabilityGenSelIndex=genStab, envMeanScores=envMeanScores)
}