
library(expm)
library(Rcpp)
library(devtools)
#########--------------------------------------------FUNCTIONS---------------------------------------------############




#the following function is needed in order to find the mu's (i.e. the parameter that defines that policy) for: OMT FDR; marginal locFDR for FDR; OMT mFDR; marginal locFDR for FDR

fdep = function(musmFDR,  musFDR, muspFDR, musmargmFDR,  musmargFDR, musmargpFDR, var1, var2, cv, rhovec, Kvec, prob, block_matrices, block_sum_matrices, block_nones, block_detfacs,block_deltas, block_invdeltas, maxiter=500000){ #SR+SR1

  
  nummus =max(length(musFDR), length(musmFDR),length(muspFDR), length(musmargmFDR),  length(musmargFDR), length(musmargpFDR))
  

  #column order: OMTFDR, OMTpFDR, OMTmFDR, marglocfdrFDR, marglocfdrpFDR, marglocfdrmFDR
  
  
  minprob_mat = matrix(0, nrow =nummus, ncol = 6 )
  lev_mat = matrix(0, nrow =nummus, ncol = 6 ) #the FDR
  pow_mat= matrix(0, nrow =nummus, ncol = 6 )
  ev_mat = matrix(0, nrow =nummus, ncol = 6 )
  er_mat =   matrix(0, nrow =nummus, ncol = 6 )
  #  minprob=lev = pow=ev=maxr=numeric(length(mus))  
  
  
  
  
  for (iter in 1:maxiter){
    
    blockbeta = rbeta(var1, var2, cv, rhovec, Kvec, prob) #SR
    
    #for oracle rule
    out = AandB(blockbeta, Kvec, block_matrices, block_deltas, block_invdeltas, block_sum_matrices,block_nones,block_detfacs) #SR+SR1
    a = out$a
    b=out$b
    
    #for marginal rule
    marglocfdr = NULL
    # Given beta for block of size K compute the locfdr values
    nb = length(blockbeta)
    for (nbi in 1:nb){
      margprob = (1-prob)*dnorm(as.vector(blockbeta[[nbi]]), mean = 0, sd = sqrt(var1))+prob*dnorm(as.vector(blockbeta[[nbi]]), mean = cv, sd = sqrt(var1+var2)) #SR
      marglocfdr =c(marglocfdr, (1-prob)*dnorm(as.vector(blockbeta[[nbi]]), mean = 0, sd = sqrt(var1))/margprob)
    }
    omarglocfdr = sort(marglocfdr)
    amarg = 1-omarglocfdr
    bmarg =  BZCpp(omarglocfdr)
    
    
    
    for (mui in 1:length(musFDR)){

      mu=musFDR[mui]
      
      #OMTFDR
      Rz = a-mu*b
      Dz = DCpp(Rz)
      ind=1
      lev_mat[mui, ind] = lev_mat[mui, ind]+sum(b[Dz==1])
      pow_mat[mui, ind] = pow_mat[mui, ind]+sum(a[Dz==1])
      minprob_mat[mui, ind] = minprob_mat[mui, ind]+Dz[1]
      ev_mat[mui,ind]=ev_mat[mui,ind]+sum(out$olocfdr[Dz==1])
      er_mat[mui,ind]=er_mat[mui,ind]+sum(Dz)
    }
    
    for (mui in 1:length(muspFDR)){
      
      mu=muspFDR[mui]
      
      #OMTpFDR
      Rz = a-mu*b
      Rz[1]=Rz[1]+mu*alpha
      Dz = DCpp(Rz)
      ind=2
      lev_mat[mui, ind] = lev_mat[mui, ind]+sum(b[Dz==1])
      pow_mat[mui, ind] = pow_mat[mui, ind]+sum(a[Dz==1])
      minprob_mat[mui, ind] = minprob_mat[mui, ind]+Dz[1]
      ev_mat[mui,ind]=ev_mat[mui,ind]+sum(out$olocfdr[Dz==1])
      er_mat[mui,ind]=er_mat[mui,ind]+sum(Dz)
    }
    
    for (mui in 1:length(musmargFDR)){
      
      mu=musmargFDR[mui]
      
      
      #marglocfdrFDR
      Rz = amarg-mu*bmarg
      Dz = DCpp(Rz)
      ind = 4
      lev_mat[mui, ind] = lev_mat[mui, ind]+sum(bmarg[Dz==1])
      pow_mat[mui, ind] = pow_mat[mui, ind]+sum(amarg[Dz==1])
      minprob_mat[mui, ind] = minprob_mat[mui, ind]+Dz[1]
      ev_mat[mui,ind]=ev_mat[mui,ind]+sum(omarglocfdr[Dz==1])
      er_mat[mui,ind]=er_mat[mui,ind]+sum(Dz)
    }
    
    for (mui in 1:length(musmargpFDR)){
      
      mu=musmargpFDR[mui]
      
      #marglocfdrpFDR
      Rz =amarg-mu*bmarg
      Rz[1]=Rz[1]+mu*alpha
      Dz = DCpp(Rz)
      ind = 5
      lev_mat[mui, ind] = lev_mat[mui, ind]+sum(bmarg[Dz==1])
      pow_mat[mui, ind] = pow_mat[mui, ind]+sum(amarg[Dz==1])
      minprob_mat[mui, ind] = minprob_mat[mui, ind]+Dz[1]
      ev_mat[mui,ind]=ev_mat[mui,ind]+sum(omarglocfdr[Dz==1])
      er_mat[mui,ind]=er_mat[mui,ind]+sum(Dz)
      
    }
    
    for (mui in 1:length(musmFDR)){
      mu=musmFDR[mui]
      #OMTmFDR
      Rz = a-mu*(out$olocfdr-alpha)
      Dz = DCpp(Rz)
      ind=3
      lev_mat[mui, ind] = lev_mat[mui, ind]+sum(b[Dz==1])
      pow_mat[mui, ind] = pow_mat[mui, ind]+sum(a[Dz==1])
      minprob_mat[mui, ind] = minprob_mat[mui, ind]+Dz[1]
      ev_mat[mui,ind]=ev_mat[mui,ind]+sum(out$olocfdr[Dz==1])
      er_mat[mui,ind]=er_mat[mui,ind]+sum(Dz)
    }  
    
    for (mui in 1:length(musmargmFDR)){
      mu=musmargmFDR[mui]
      
      #marglocfdrmFDR
      Rz = amarg-mu*(omarglocfdr-alpha)
      Dz = DCpp(Rz)
      ind = 6
      lev_mat[mui, ind] = lev_mat[mui, ind]+sum(bmarg[Dz==1])
      pow_mat[mui, ind] = pow_mat[mui, ind]+sum(amarg[Dz==1])
      minprob_mat[mui, ind] = minprob_mat[mui, ind]+Dz[1]
      ev_mat[mui,ind]=ev_mat[mui,ind]+sum(omarglocfdr[Dz==1])
      er_mat[mui,ind]=er_mat[mui,ind]+sum(Dz)
    }      
    
    
    
    
  }#for (iter in 1:maxiter){
  cat(iter,"\n")
  #  cat (mus,"\n")
  cat (lev_mat/iter,"\n")
  cat (pow_mat/iter,"\n")
  cat (minprob_mat/iter,"\n")
  
  return(list(lev_mat = lev_mat/maxiter,pow_mat = pow_mat/maxiter, minprob_mat = minprob_mat/iter, ev_mat = ev_mat/iter, er_mat = er_mat/iter))
}


#the necessary functions for a given policy


rbeta = function(var1, var2, mu, rhovec, Kvec, prob, withh = FALSE){ #sample from a mixture model #SR
  nb = length(Kvec)
  block_beta = list()
  block_h = list()
  vec_h = NULL
  for (b in 1:nb){
    K = Kvec[b]
    rho = rhovec[b]
    
    
    #-----for one block with symmetric correlation
    covmat =  (matrix (rho*var1, nrow=K,ncol=K) + diag(rep((1-rho),K))) *var1
    #    sqrtmat = sqrt(var1)* matrix( (sqrt(1+(K-1)*rho) -sqrt(1-rho))/K,ncol=K,nrow=K) 
    #    diag(sqrtmat) = sqrt(var1)*(sqrt(1+(K-1)*rho) + (K-1)*sqrt(1-rho))/K
    
    
    h = rbinom(K,size=1,p=prob)
    covmat = covmat+diag(h)*var2
    sqrtmat = sqrtm(covmat)
    beta =  sqrtmat%*%rnorm(K)+mu*h #SR
    block_beta[[b]] = beta
    block_h[[b]]= h
  }
  if (!withh){
    return(block_beta)
  }
  if (withh){
    for (b in 1:nb){
      vec_h = c(vec_h, block_h[[b]] ) 
    }
    return(list(block_beta = block_beta, block_h = block_h, vec_h=vec_h))
  }
}  



AandB = function(blockbeta, Kvec, block_matrices, block_deltas, block_invdeltas, block_sum_matrices, block_nones,block_detfacs){#SR+SR1
  locfdr = NULL
  # Given beta for block of size K compute the locfdr values
  nb = length(blockbeta)
  for (nbi in 1:nb){
    Kb = Kvec[nbi]
    matrices = block_matrices[[nbi]]
    sum_matrices = block_sum_matrices[[nbi]]
    nones = block_nones[[nbi]]
    deltas = block_deltas[[nbi]] #SR
    invdeltas = block_invdeltas[[nbi]] #SR
    beta = as.vector(blockbeta[[nbi]])
    detfacs=block_detfacs[[nbi]]
    ##S1 = (matrices[[Kb]]%*%beta)*beta #2^K*K vector #SR-old
    S1 = (matrices[[Kb]]%*%beta)*(beta-2*deltas[[Kb]]) + invdeltas[[Kb]]*deltas[[Kb]] #SR+SR1 - need to make sure works correctly...
    M1 = matrix(S1,nrow=Kb) #turn to a matrix that will sum the columns for each h to get beta'%*%sigmainv_h%*%beta
    pzh = detfacs[[Kb]]*exp(-0.5*apply (M1,2, sum))* prob^nones[[Kb]]*(1-prob)^(Kb-nones[[Kb]]) #joint distribution of Z and h. 
    locfdr = c(locfdr, 1 - (sum_matrices[[Kb]]%*%pzh) /sum(pzh)) # gives vector of length k
  }
  
  
  K = length(locfdr)
  oo = order(locfdr)
  locFDR=locfdr[oo]
  clocFDR = cumsum(locFDR)
  a = 1-locFDR
  b = locFDR/(1:K) - c(0,clocFDR[-K]/(1:(K-1))/(2:K))
  return(list(a=a, b=b, locfdr = locfdr, olocfdr = locFDR))
  
}

cppFunction('IntegerVector DCpp(NumericVector R) {
            int K=R.size();
            IntegerVector D(K,0);
            double cR=0.;
            
            for (int i =0;i<K;i++){
            cR=0.;
            for (int j = i;j<K;j++){
            cR += R[j];
            if (cR>0) {D[i]=1;break;}
            }
            if (D[i]==0) break;
            } 
            return D;
            }')  


cppFunction('NumericVector BZCpp(NumericVector Tz){
            int K=Tz.size();
            NumericMatrix S(K,K);
            std::fill(S.begin(), S.end(), 0.);
            NumericVector b(K,0.);
            S(0,0)=1;
            for (int i=1;i<K;i++){
            S(i,0) = S(i-1,0)*(1-Tz[i-1]);
            for (int v=1;v<=i;v++){
            S(i,v) = S(i-1,v)*(1-Tz[i-1])+S(i-1,v-1)*Tz[i-1];
            }
            }
            
            b[0] = S(0,0)*Tz[0];
            for (int i=1;i<K;i++){
            for (int j = 0; j <= i;j++){
            b[i] += S(i,j)* ((i-j)*Tz[i] - j*(1-Tz[i]))/(i*(i+1));
            }
            }
            
            return b;
            }')


#########--------------------------------------------END FUNCTIONS---------------------------------------------############
