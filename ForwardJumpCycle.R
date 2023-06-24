#Author: Pasin Marupanthorn
#Co-Author: Rachanai Kaikeaw
#This code contains calculating functions used in the paper 
#Forward Jump Random Walk on a Cycle Graph and Its Hitting Time

##################################################################
# Probability of Arriving State
##################################################################
# - r is the label of arriving state
# - n is the number of walking step
# - d is the number of vertices labelled from 0 to d-1
# - m is the number of the maximum jump from initial node to the next node
#output
# Probability of the arriving state at node r

pmf.arr <- function(r,n,d,m){
  z <- 0+1i
  w <- exp(2*pi*z/d)
  sumt <- 0
  for(k in 1:(d-1)){
    vk <- ((w^k/m)*((w^(m*k)-1)/(w^k - 1)))
    sumt <- sumt + w^(-r*k)*vk^n
  }
  return(Re(sumt*(1/d) + 1/d))
}

##################################################################
# Distribution of Arriving State
##################################################################
#Input
# - n is the number of walking step
# - d is the number of vertices labelled from 0 to d-1
# - m is the number of the maximum jump from initial node to the next node
#output
# Distribution of the arraiving state for nodes in the graph

pmf.arr.dist <- function(n,d,m){
  dist.r <- rep(0,d)
  for(r in 0:(d-1)){
    dist.r[r+1] <- pmf.arr(r,n,d,m)
  }
  return(dist.r)
}



##################################################################
# PMF of Hitting Time 
##################################################################
#Input
# - r is the label of arriving state
# - d is the number of vertices labelled from 0 to d-1
# - m is the number of the maximum jump from initial node to the next node
# - lim.n is maximum  number of walking step shown in density 
#output
# PMF of Hitting Time 




pmf.hit <- function(r,lim.n,d,m){
  if(n == 1){
    hit.re <- pmf.arr(r,n,d,m)
  }else
  {
    hit.re <- rep(0,lim.n-1)
    hit.re[1] <- pmf.arr(r,1,d,m)
    for(k in 2:(lim.n-1)){
      r.re <- rep(0,k-1)
      for(j in 2:k){
        r.re[j-1] <- pmf.arr(0,k-(j-1),d,m) 
      }
      hit.re[k] <- pmf.arr(r,k,d,m) - sum(hit.re[1:(k-1)]*r.re)
    }
  }
    
  return(hit.re)
}




##################################################################
# Expected Value of Hitting Time 
##################################################################
#Input
# - r is the label of arriving state
# - d is the number of vertices labelled from 0 to d-1
# - m is the number of the maximum jump from initial node to the next node
#output
# Expected Value of Hitting Time from node 0 to node r 

ex.hit <- function(r,d,m){
  w <- exp(2*pi*z/d)
  sumt <- 0
  for(k in 1:(d-1)){
    vk <- ((w^k/m)*((w^(m*k)-1)/(w^k - 1)))
    sumt <- sumt + ((1-w^(-r*k))*vk)/(1-vk)
  }
  return(Re(sumt) + d)
}


##################################################################
# Variance of Hitting Time 
##################################################################
#Input
# - r is the label of arriving state
# - d is the number of vertices labelled from 0 to d-1
# - m is the number of the maximum jump from initial node to the next node
#output
# Variance of Hitting Time from node 0 to node r 

var.hit <- function(r,d,m){
  w <- exp(2*pi*z/d)
  sumt.1 <- 0
  sumt.2 <- 0
  sumt.3 <- 0
  for(k in 1:(d-1)){
    vk <- ((w^k/m)*((w^(m*k)-1)/(w^k - 1)))
    sumt.1 <- sumt.1 + ((1-w^(-r*k))*vk)/(1-vk)^2
    sumt.2 <- sumt.2 + ((1-w^(-r*k))*vk)/(1-vk)
    sumt.3 <- sumt.3 + ((1+w^(-r*k))*vk)/(1-vk)
  }
  return(2*Re(sumt.1) + (d + Re(sumt.2))*(d-1+Re(sumt.3)))
}
