


# This file contains functions to estimate the largest eigenvalue of 
# a partially known network using hybrid sampling.

# Last modified: 26 September 2023

# Authors: Omar De La Cruz Cabrera and Razan Alsehibani.


# #
library(igraph)
library(Matrix)

nodes = 200
samplings = 100 ## repetitions 

## UK parameters
alpha = 0.0002 # Rate of contact at random
delta = 0.0002 # Rate of contact through network
rao = alpha + delta # Rate of contact at random or through network 
PNK = 1 ## Probability of meeting new unknown individual per unit time (1 means every day we will meet eta individuals)
eta = 1 #Number of individuals being sampled per unit time
AB2<-readMM(file="./matrixData/AdjBAComplex_Power_1_.mtx")

A<-AB2*1 

AAB<- graph.adjacency(A, mode="undirected", weighted=NULL) 
plot(AAB)
truePerronEvalueB = eigen(A)$values[1]
no_of_methods = 5
zeta = (c(0, 0.25, 0.5, 0.75, 1))
Perron_Eigen_Value = array (0, dim = c(nodes, samplings, no_of_methods))

seed = sample( 1:nodes,  size=eta, replace = FALSE,  prob = rep(1,nodes))

for (j in 1:samplings){
  
  k = matrix(0, no_of_methods, nodes)
  Prob_Weights = matrix(1, no_of_methods, nodes)
  P = matrix(1, no_of_methods, nodes)
  
  
  
  Sampled = matrix(0, no_of_methods, nodes)
  
  for (loop in 1: no_of_methods) { 
    Sampled[,1] = seed
    k[,seed] = 1
    P[,seed] = 0
  }
  
  
  for (i in 2:nodes){
    
    
    for (loop in 1:no_of_methods){ 
      Perron_Eigen_Value[i, j, loop]  =  eigen(crossprod(AAB[,Sampled[loop,]]))$val[1]
    }
    
    for (loop in 1:(no_of_methods-1))   {
      
      
      k_temp = (!k[loop,])
      exp_temp = (1-exp(-rao*(A %*% k[loop,])))
      Prob_Temp = k_temp * exp_temp
      Prob_Weights[loop,] = Prob_Temp[,1]
      
      
    }
    
    Prob_Weights[5, ] = !k[5,]
    
    
    for (loop in 1:no_of_methods) {
      Prob_Weights[loop, ] = Prob_Weights[loop, ]/(sum(Prob_Weights[loop, ]))
    }
    
    for (loop in 1: no_of_methods){ 
      
      
      P[loop, ] =   (!k[loop,]/sum(!k[loop,])) * zeta[loop] +  Prob_Weights[loop, ] * (1-zeta[loop])
    }
    
    
    
    for (loop in 1:no_of_methods)   {     
      indivSample = sample( 1:nodes,  size=eta, replace = FALSE,  prob = P[loop,])
      tmp = indivSample
      k[loop, tmp] = 1
      Sampled[loop,i] = tmp
      Prob_Weights[loop, tmp ] = 0
      P[loop,tmp] = 0
    }
    
  } # # i loop
  
}   # # j loop



avge = matrix(0, no_of_methods, nodes)
for (loop in 1:no_of_methods) { 
  avge[loop,] = rowSums(Perron_Eigen_Value[, ,loop])/samplings
}

par( mfrow= c(1,1) )
plot(avge[1,],type = "l", col=1, xlim=c(1,200), ylim = c(1,truePerronEvalueB^2),xlab="Number of sampled nodes",ylab="Largest eigenvalue", main="Complex Erdős–Rényi")
lines(avge[2,],type = "l", col=2)
lines(avge[3,],type = "l", col=3)
lines(avge[4,],type = "l", col=4)
lines(avge[5,],type = "l", col=5)

abline(h=truePerronEvalueB^2,lty=2)
legend("bottomright", legend=c("Z = 0", "Z = 0.25 (HS)", "Z = 0.50 (HS)", "Z = 0.75 (HS)", "Z = 1 (SRS)"),
       col=1:5, lty=1, cex=0.8,box.lty=0, inset=.02)

plot.igraph(AAB, vertex.size=6, vertex.label=NA)

