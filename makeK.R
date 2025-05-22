# Helper function to generate K matrix for the NIMBLE-based HBKMR model

makeK <- nimbleFunction(
  run = function(Z = double(2), r = double(1), N = integer(0), M = integer(0)){
    
    K <- nimMatrix(value = 0, nrow = N, ncol = N)
    
    for(i in 1:N){
      # Diagonal entries will be exp(0) = 1
      K[i, i] <- 1.0  
      for(j in 1:N){
        # Avoids doing the same computation for the upper triangular 
        if(i < j){
          dist_ij <- 0.0
          for(m in 1:M){
            dist_ij <- dist_ij + r[m] * (Z[i, m] - Z[j, m])^2
          }
          K[i, j] <- exp(-dist_ij)
          K[j, i] <- K[i, j]
        }
      }
    }
    returnType(double(2))
    return(K[1:N, 1:N])
  }
)

compiled_makeK <- compileNimble(makeK)
