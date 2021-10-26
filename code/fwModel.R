# Define a custom function fwModel to fit Finlay Wilkinson model 
# for a trait of interest from the input data
fwModel <- function(traitName){
  
  # Invoke the function gxeFw from statgenGxE library for model fitting
  geFW <- gxeFw(TD = setATD, trait = traitName)
  
  genSens <-geFW$estimates   # Extract  genotype sensitivities estimate
  
  mostSensGen<- genSens[order(-genSens$sens),] # most sensitive genotypes
  leastSensGen <- genSens[order(genSens$sens),] # least sensitive genotypes
  
  # FW plot
  #plot(geFW, plotType = 'line')
  
  # write the output as a list object
  fwOutput <- list(anovaSumm=geFW, genSens=genSens, mostSensGen=mostSensGen,leastSensGen=leastSensGen)
  
}