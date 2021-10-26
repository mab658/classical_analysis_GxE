# Define a custom function ggeModel to fit genotype and genotype by environment (GGE) model 
# for a trait of interest from the input data

fitGGE <- gxeGGE(TD = setATD, trait = "fyld")
ggeModel <- function(traitName){
  
  fitGGE <- gxeGGE(TD = setATD, trait = traitName)
  summOutput <- list(Anova=fitGGE$anova,PCA=fitGGE$importance)
  
  # extract genotypic scores from GGE model
  genScores <- fitGGE$genoScores
  
  # Create GGE2 plot
  
  #jpeg(file="GGE2_biplot.jpeg")
  GGE2_biplot <- plot(fitGGE, plotType = "GGE2", scale = 0.5, plotConvHull = TRUE,  
                      sizeGeno = 2.8,sizeEnv = 2.5, title = "")
  #dev.off()
  
  ggeOutput <- list(anovaPCA=summOutput, genScores=genScores)
  
}
