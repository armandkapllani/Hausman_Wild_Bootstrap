#==================================================#
# The Wild Bootstrap Hausman Test in Linear Models #
#==================================================#
#============#
# Simulation #
#============#

# Set seed #
#==========#
set.seed(12345)

N = 100  # Set the sample size 
B = 1000 # Set the number of bootstrap draws
M = 2000 # Set the number of repetitions

# Column vector of betas #
#========================#
betas = matrix(0, nrow = 2, ncol = 1)

# Correlation coefficient column vector rho_xz #
#==============================================#
rho_xz = matrix(c(seq(0.1, 0.9, by = 0.2)), ncol = 1)

# Correlation coefficient column vector rho_xu #
#==============================================#
rho_xu = matrix(c(0, rho_xz), ncol = 1)

# All 30 possible pair values (rho_xz, rho_xu) #
#==============================================#
pairs = expand.grid(rho_xz, rho_xu)
pairs = data.matrix(pairs)
pairs = pairs[order(pairs[,1]),]  # sort the pairs in increasing order by rho_xz

# Store p-values #
#================#
p0 = matrix(0,nrow = M, ncol = nrow(pairs))
p = matrix(0,nrow = M, ncol = nrow(pairs))

# Store Simulation results # 
#==========================#
sim = matrix(0,5,6)

# Computing time #
#================#
tic()

# Start the simulation 
#=====================#
for(j in 1:M){
  
  # Column vector of standard normal N x 1 (Instrumental variables) #
  #=================================================================#
  Z = matrix(rnorm(N), ncol = 1)
  
  # Column vector of standard normal N x 1 (xi) #
  #=============================================#
  xi = matrix(rnorm(N), ncol = 1)
  
  # Column vector of standard normal N X 1 #
  #========================================#
  eta = matrix(rnorm(N), ncol = 1)
  
  # Clustering errors generated through mu (N/2 clusters, 2 units for each cluster) #
  #=================================================================================#
  upsilon = matrix(rnorm(50), ncol = 1)
  mu = kronecker(upsilon, c(1,1))
  
  # This nested for loop will create 30 samples (we have 30 pairs of (rho_xz,rho_xu)) for each (Z, xi, eta, upsilon)
  for(k in 1:nrow(pairs)){
    
    # Generate X, V, U, Y [All column vectors N x 1]
    
    X = pairs[k,1]*Z + xi 
    X = cbind(1, X) # Add vector of ones for the intercept 
    V = pairs[k,2]*xi + sqrt(Z^2+1)*eta
    U = mu + V 
    Y = X%*%betas + U
    Z = cbind(1, Z) # Add vector of ones for the intercept 
    
    beta.u = solve(t(X)%*%X)%*%t(X)%*%Y
    beta.r = solve(t(Z)%*%X)%*%t(Z)%*%Y
    
    resid.u = Y - X%*%beta.u
    resid.r = Y - X%*%beta.r
    
    sigma2.beta.u = (sum(resid.u^2)/(N-2))*solve(crossprod(X,X))
    sigma2.beta.r = (sum(resid.r^2)/(N-2))*solve(crossprod(Z,X),crossprod(Z,Z))%*%solve(crossprod(Z,X))
    H = t(beta.u - beta.r)%*%solve(sigma2.beta.r - sigma2.beta.u)%*%(beta.u - beta.r)
    p0[j,k] = pchisq(H,2,lower.tail = FALSE)
    
    # Wild Bootstrap 
    
    dN = matrix(0,B,2)
    beta.u.boot = matrix(0,B,2)
    beta.r.boot = matrix(0,B,2)
    
    for(i in 1:B){
      g = kronecker(rbinom(50,1,0.5)*2 - 1, c(1,1))
      resid.u.boot = resid.u*matrix(g, ncol = 1)
      resid.r.boot = resid.r*matrix(g, ncol = 1)
      beta.u.boot[i,] = t(solve(crossprod(X,X),crossprod(X,resid.u.boot)))
      beta.r.boot[i,] = t(solve(crossprod(Z,X),crossprod(Z,resid.u.boot)))
    }
    
    dN = beta.u.boot - beta.r.boot
    var.dN = var(dN)
    WH = t(beta.u - beta.r)%*%ginv(var.dN)%*%(beta.u - beta.r)
    
    WHb = matrix(0, nrow = B)
    for(i in 1:B){
      WHb[i] = (beta.u.boot[i,] - beta.r.boot[i,])%*%ginv(var.dN)%*%(beta.u.boot[i,] - beta.r.boot[i,])
    }
    
    p[j,k] = sum(ifelse(WHb>WH[1],1,0))/B
    
    Z = Z[,-1] # Drop the column vector of ones
    
  }
}

# Store the simulation results
sim = matrix(0, nrow = 30, ncol = 6)  # Each row represents a pair (rho_xz, rho_xu), each column represents sig.level (0.01, 0.05, 0.10) for Hausman test and Wild Hausman test


for(i in 1:30){
  
  # Size of given level of significance of conventional Hausman's test
  sim[i,1] = sum(ifelse(p0[,i]<0.01,1,0))/M
  sim[i,2] = sum(ifelse(p0[,i]<0.05,1,0))/M
  sim[i,3] = sum(ifelse(p0[,i]<0.10,1,0))/M
    
  # Size of given level of significance of Wild bootstrap Hausman's test
  sim[i,4] = sum(ifelse(p[,i]<0.01,1,0))/M
  sim[i,5] = sum(ifelse(p[,i]<0.05,1,0))/M
  sim[i,6] = sum(ifelse(p[,i]<0.10,1,0))/M

}

write.csv(sim, file = "/Users/armandkapllani/Desktop/Simulation/sim.csv")

toc()
