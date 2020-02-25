#### All replicates for all possible Dij,
##variance is estimated for Tau^2 using a MAD,
#but sigma^2 uses the variance across replicates
#for each gene and experiment, sigma^2_{ij}.



library(tidyverse)
library(testthat)
library(reshape2)

##a function to calculate likelihood of network state given
#experimental pertubation data
# requires data matrix of expression measurements
#and calculates the -log(likelihood) for each possible
#value of each gene in each experiment using
#the ternary network mixture model
#uniform, normal, uniform for lowered, baseline, or heightened expression


AllLikeAllRepsAllDijThetaOut <- function(z){

  #matrix for variance of all replicates in each gene and experiment
  rep_var_mat <- matrix(nrow = nrow(z),ncol = nrow(z))
  colnames(rep_var_mat) <- unique(colnames(z))
  
    #matrix for mean of all replicates in each gene and experiment
  rep_mean_mat <- matrix(nrow=nrow(z),ncol=nrow(z))
  
  for(i in 1:nrow(z)){ 
    #number of replicates per experiment
    Kj  <-  table(colnames(z))[unique(colnames(z))[i]]

    
    for(j in 1:nrow(z)){
      
      rep_var_mat[i,j] <-  ((Kj-1)/Kj) * var(z[i,which(colnames(z)==names(table(colnames(z))[j]))])
      rep_mean_mat[i,j] <-  mean(z[i,which(colnames(z)==names(table(colnames(z))[j]))])  
    }
  }
  
##Tau squared, the variance of the baseline normal distribution
  #mean will remain fixed, and it is roughly estimated
  #by the MAD for the means across replicates for each
  #gene in each experiment.

  tau2 <- (1.4826*mad(as.vector(rep_mean_mat)))^2
  
  
  
  # will use the genewise max and mins for UMVUE for bounds of the uniform
  deltaNegs <- apply(X = z,MARGIN = 1,FUN = function(x){min(x)+ (min(x)/length(x))})
  deltaPos <- apply(X = z,MARGIN = 1,FUN = function(x){max(x)+ (max(x)/length(x))})
  
  ##list to store data.frame results from three states
  lsrepLikelihoods <- NULL
  
  for(j in 1:3){
    
    #data.frame for experimental likelihoods
    ExperRepLikelihoods <- NULL
    
    for(k in 1:ncol(z)){
      #a vector to hold each experiments likelihoods
      like <- NULL
      
      for(i in 1:nrow(z)){
        
        #sets diagonal to 0, so there is no contribution to the log-likelihood
        if(grepl(rownames(z)[i],colnames(z)[k])){
          like[i] <- 0
          next
        }#end if  
        
        
        #likelihood of gene i measurement if up regulation
        if(j == 1){
          
          state <- 1
          
          fn <- function(x){(deltaPos[i]*sqrt(2*pi*rep_var_mat[i,colnames(z)[k]]))^(-1) * exp(-((z[i,k]-x)^2)/(2*rep_var_mat[i,colnames(z)[k]]))}
          
          like[i] <- -log(integrate(f = fn,lower = 0,upper = deltaPos[i])$value)
          
          
          #likelihood of gene i measurement if baseline regulation
        }else if(j == 2){
          
          state <- 0
          
          fn <- function(x){(2*pi*sqrt(rep_var_mat[i,colnames(z)[k]]*tau2))^(-1) * exp(-((z[i,k]-x)^2)/(2*rep_var_mat[i,colnames(z)[k]])) *
              exp(-(x^2)/(2*tau2))}
          
          like[i] <- -log(integrate(f=fn,lower= -Inf,upper=Inf)$value)
          #likelihood of gene i measurement if down regulation
        }else if(j == 3){
          
          state <- -1
          
          fn <- function(x){(-deltaNegs[i]*sqrt(2*pi*rep_var_mat[i,colnames(z)[k]]))^(-1) *
              exp(-((z[i,k]-x)^2)/(2*rep_var_mat[i,colnames(z)[k]]))}
          
          like[i] <- -log(integrate(f=fn,lower=deltaNegs[i],upper=0)$value)
        }
      }#end replicate loop
      repLikelihoods.df <- data.frame(like,rownames(z),
                                      rep(state,nrow(z)),rep(colnames(z)[k],nrow(z)))
      colnames(repLikelihoods.df) <- c("RepLikelihoods","Gene","State","ExperType")
      
      ExperRepLikelihoods[[k]]  <- repLikelihoods.df
      
    }  #end experiment loop
    
    
    lsrepLikelihoods[[j]] <- ExperRepLikelihoods 
  }#end dij loop
  
  
  df.downreg20 <- bind_rows(lsrepLikelihoods[[3]], .id = "column_label")[,-1]
  
  
  df.basereg20 <- bind_rows(lsrepLikelihoods[[2]], .id = "column_label")[,-1]
  
  
  df.upreg20 <- bind_rows(lsrepLikelihoods[[1]], .id = "column_label")[,-1]
  

  df.ll20 <- rbind(df.downreg20,df.basereg20,df.upreg20)
  #reordering
  df.ll20 <- df.ll20[order(df.ll20$ExperType),]
  df.ll20 <- df.ll20[order(df.ll20$Gene),]
  
  #some null vectors to hold likelihoods
  #and identifying experiments,genes, etc.
  i_exp <- rep(x=seq(1,20)-1,each = 60)
  i_node <- rep(x = seq(1,20)-1,each = 3,times=20)
  outcome <- rep(c(-1,0,1),400)
  value <-  vector(mode="numeric",length=1200)
  
#calculates combined likelihoods
  count <- 1
  for(r in 1:20){
    for(g in 1:20){  
      for(s in 1:3){
        
        score <- df.ll20  %>% filter(Gene == sort(unique(Gene))[g] &  
                                       ExperType == sort(unique(ExperType))[r] & State == unique(State)[s]) %>%
          summarise(sum(RepLikelihoods))
        
        value[count] <- score[1,1]
        
        
        ################1/29/19 perturbs need redone anyway
        #if(startsWith(sort(unique(df.ll20$ExperType))[r],
        #      as.character(sort(unique(df.ll20$Gene))[g]))){ 
        #       perturbation[count] <- 1}else{
        #       perturbation[count] <- 0}
        
        
        count = count + 1
      }  #placement of count update determines if you're doing one long 1200 row set of likelihoods
      #and in what order, of if you're doing 3 columns set up as to compare to original probs graph
      #count = count + 1
    }
  }
  
  is_perturbation <- vector(length = 1200, mode = "integer")
  is_perturbation[which(value == 0)] <- 1
  
  ##This is the scores table for the loglikelihood method####
  scores_ll20 <- data.frame(i_exp,i_node,outcome,value,is_perturbation)
    
  return(scores_ll20)
  
}#end funtion

df.ll20 <- AllLikeAllRepsAllDijThetaOut(z) 
