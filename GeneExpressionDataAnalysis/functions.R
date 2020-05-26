

Zrand2sidedalternative = function(K,propvec, cvvec, sigmavec){ #sample from a mixture model
  uvec = runif(K)
  J = length(propvec)
  ret = numeric(K)
  
  Hnull = which(uvec<propvec[1])
  ret[Hnull] = rnorm(length(Hnull), mean = cvvec[1], sd = sigmavec[1]) 
  
  for (j in 2:J){
    if (j<J){
    Hj = which (uvec>= sum(propvec[1:(j-1)]) & uvec<sum(propvec[1:j])) 
    }
    if (j==J){
      Hj = which (uvec>=sum(propvec[1:(j-1)]))
    }
    ret[Hj] = rnorm(length(Hj), mean = cvvec[j], sd = sigmavec[j]) 
  }
  ret
}





marg.dense2sidedalternative = function (z,propvec, cvvec, sigmavec){ #compute the marginal density of a  mixture model 
  J = length(propvec)
  out=rep(0, length(z))
  for (j in 1:J){
    out = out + propvec[j]*dnorm(z, mean = cvvec[j], sd = sigmavec[j])  
  }  
  return(out)
  
}

scaled.null.dens =function(z,propvec, cvvec, sigmavec, alternative = c("two.sided", "less", "greater")){ 
  J = length(propvec)
  out=rep(0, length(z))
  if (alternative =="less"){
  for (j in 1:J){
    if(cvvec[j]>=0){
    out = out + propvec[j]*dnorm(z, mean = cvvec[j], sd = sigmavec[j])}  
  } 
    return(out)
  }  

  if (alternative =="greater"){
    for (j in 1:J){
      if(cvvec[j]<=0){
        out = out + propvec[j]*dnorm(z, mean = cvvec[j], sd = sigmavec[j])}  
    } 
    return(out)
  }  
  if (alternative =="two.sided"){
    return(propvec[1]*dnorm(z, mean = cvvec[1], sd = sigmavec[1]))
  }
}


fZV = function(K,  mus, Rreg, Cpp, propvec, cvvec,sigmavec, alternative = c("two.sided", "less", "greater"), maxiter=500000){#in propvec and cvvec, respectively the mean and  proportion of zero, pos, neg 
  lev = pow=ev=maxr=numeric(length(mus))
  maxK = ifelse(K>200, max(round(K*(1-propvec[1])),200),K)
  for (iter in 1:maxiter){
    z = Zrand2sidedalternative(K, propvec, cvvec, sigmavec) #if sample from the mixture density with parameters propvec and cvvec and sigmavec
    Pz = marg.dense2sidedalternative(z,propvec, cvvec, sigmavec)
    Tz = scaled.null.dens(z,propvec, cvvec, sigmavec, alternative =alternative )/Pz 
    Tz = sort(Tz)
    az = 1-Tz[1:maxK]
    bz = BZCpp(Tz[1:maxK])
    for (mui in 1:length(mus)){
      mu=mus[mui]
      Rz = az-mu*bz
      Dz = DCpp(Rz)
      lev[mui] = lev[mui] + sum(bz[Dz==1])
      pow[mui] = pow[mui] + sum(az[Dz==1])
      #??     ev[mui] = ev[mui]+sum(((1-prop)*(a)/(prop*gU))[D==1])
      maxr[mui] = max(maxr[mui],sum(Dz==1))
    } 
   }
  
  cat(iter,"\n")
  cat (mus,"\n")
  cat (lev/iter,"\n")
  cat (pow/iter,"\n")
  
  
  return(list(lev = lev/maxiter,pow = pow/maxiter))#??, ev = ev/maxiter))
  
}

#install.packages("installr")
#library("installr")
#install.Rtools()
#Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
#Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

library(Rcpp)
library(devtools)
#assignInNamespace("version_info", c(devtools:::version_info, list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), "devtools")

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



findmuOMT = function(K, propvec, cvvec, sigmavec, alternative = c("two.sided", "less", "greater"), alpha = 0.05, mus = seq(4500,15000,5), murange=4000, maxiter= 50000){
  
  out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative, maxiter=maxiter)
  
  mbest = which.min(abs(out$lev-alpha))
  #if mbest equals 1 start again with all mu's below min(mus), if mbest equals length(mus) start again with all mu's above max(mus)
  while(mbest==1 | mbest== length(mus)){
    if (mbest==1){mus = seq(max(0.1, min(mus)-murange), min(mus),10)}
    if (mbest==length(mus)){mus = seq(max(mus), max(mus)+murange,10)}
    out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
    mbest = which.min(abs(out$lev-alpha))
  }
  print(cbind(mbest, mus[mbest]))
  
  mus = seq(max(mus[mbest]-5,0.1),mus[mbest]+5,length.out = 2500)
  out =   fZV(K,  mus, Rreg=FALSE, Cpp=TRUE, propvec= propvec, cvvec=cvvec,sigmavec=sigmavec,alternative =alternative,maxiter=maxiter)
  mbest = which.min(abs(out$lev-alpha))
  print(cbind(mbest, mus[mbest], out$lev[mbest]))
  
  return(mus[mbest])
  
}


#---------------------------------------------------Functions from the library mixFDR by Omkar Muralidharan------------------------------------------------------
`effectSize` <-
  function(z, m, noiseSD = NA){
    e = m
    if(!is.null(m$noiseSD) && is.na(noiseSD)){
      noiseSD = m$noiseSD
    }else{
      if(is.na(noiseSD)) noiseSD = min(m$sig)	
    }
    sig = noiseSD
    if(any(e$sig<sig)) warning(paste("Marginal variance of a group is less than noise variance, this is impossible under the model."))
    gpProb = groupProbs(z,e)
    priorVar = pmax(e$sig^2 - sig^2,0)
    postVars = 1/((1/priorVar) + (1/sig^2)) * (priorVar != 0)
    postMeans =  matrix(z, length(z),length(e$mu),byrow=FALSE)/sig^2
    postMeans = postMeans + matrix((e$mu)/(priorVar),length(z),length(e$mu),byrow=TRUE)
    postMeans = postMeans * matrix(postVars, length(z), length(e$mu), byrow= TRUE)
    postMeans[is.nan(postMeans)] = 0
    postMeans = postMeans + matrix((priorVar ==0)*e$mu, length(z), length(e$mu), byrow = TRUE)
    postMeansSq = postMeans^2 + matrix(postVars, length(z), length(e$mu), byrow = TRUE)
    EdeltZ = rowSums(postMeans * gpProb)
    VdeltZ = rowSums(postMeansSq * gpProb) - EdeltZ^2
    return(cbind(Edz = EdeltZ, Vdz = VdeltZ))
  }

`tailFDRMixModel` <-
  function(z, m, nullGroups = NA){
    e = m
    if(any(is.na(nullGroups))){
      nullGroups = rep(FALSE, length(e$mu))
      nullGroups[1] = TRUE
    }
    # for each group, compute P(|z|>z)
    # P(Z > z) + P(Z < -z)
    # P(Z > (z - mu)/sigma) + P(Z < (-z - mu)/sigma)
    
    A = matrix(z, length(z), length(e$mu))
    ri = A - matrix(e$mu, length(z), length(e$mu), byrow = TRUE)
    ri = ri/matrix(e$sig, length(z), length(e$sig), byrow = TRUE)
    le = -A - matrix(e$mu, length(z), length(e$mu), byrow = TRUE)
    le = le/matrix(e$sig, length(z), length(e$sig), byrow = TRUE)
    pRi = 1 - pnorm(ri)
    pLe = pnorm(le)
    pMix2sd = (pRi + pLe)*matrix(e$pi, length(z), length(e$pi), byrow = TRUE)
    nullProb = rowSums(matrix(pMix2sd[,nullGroups],ncol = sum(nullGroups)))
    overallProb = rowSums(pMix2sd)
    twoSided = ((nullProb/overallProb)*(z>=0) + rep(sum(e$pi[nullGroups]), length(z)) *(z<0))
    ppRi = pRi*matrix(e$pi, length(z), length(e$pi), byrow = TRUE)
    right = rowSums(matrix(ppRi[,nullGroups],ncol = sum(nullGroups)))/rowSums(ppRi)
    ppRi = pnorm(ri)*matrix(e$pi, length(z), length(e$pi), byrow = TRUE)
    left = rowSums(matrix(ppRi[,nullGroups],ncol = sum(nullGroups)))/rowSums(ppRi)
    return(list(FDRTwoSided = twoSided, FDRleft = left, FDRright = right))
  }

`groupProbs` <-
  function(z,m){
    e = m
    mu = e$mu
    pi = c(e$pi)
    sig = e$sig
    phiMat = matDens(z,mu,sig)
    dens = (phiMat %*% pi)[,1]
    gpProbs = phiMat / matrix(dens, nrow(phiMat), ncol(phiMat), byrow = FALSE) * matrix(pi, nrow(phiMat), ncol(phiMat), byrow =TRUE)
  }

`fdrMixModel` <-
  function(z,m, nullGroups = NA){
    e = m
    if(any(is.na(nullGroups))){
      nullGroups = rep(FALSE, length(e$mu))
      nullGroups[1] = TRUE
    }
    G = groupProbs(z,e)
    rowSums(matrix(G[,nullGroups],ncol = sum(nullGroups)))
  }
`dens.mixture` <-
  function(x,m){
    mu = m$mu
    pi = c(m$pi)
    sig = m$sig
    phiMat = matDens(x,mu,sig)
    return((phiMat %*% pi)[,1])
  }

`getStartsNew` <-
  function(x,J, noiseSD = 1){
    # going to have equally
    # one big group in the middle
    pi = rep(1/J,J)
    mid = floor(median(1:J))
    pi[mid] = 0.75
    if(J!=1) pi[-mid] = 0.25/(J-1)
    pi = pi/sum(pi)
    q = quantile(x, c(0,cumsum(pi)))	
    mu = sig = 1:J
    for(i in 1:J){
      dat = x[x<=q[i+1] & x>=q[i]]
      mu[i] = mean(dat)
      sig[i] = noiseSD
    }
    starts = cbind(pi,mu,sig)
    tmp = starts[mid,]
    starts[mid,]= starts[1,]
    starts[1,] = tmp
    return(starts)
  }
`mixModelManyStarts` <-
  function(x, J, N, starts = NA, maxiter = 1000, tol = 0.001, p = NA, muInt = NA, sigInt = 1, theonull = FALSE){
    
    # try different starting points
    # getStartsNew
    # all null
    # random
    # choose highest penalized likelihood
    
    if(any(is.na(starts))){
      starts = array(0, dim=c(5, J, 3))
      sigMAD = median(abs(x))
      sigST = sigInt[1,1]
      if(sigInt[1,1]/sigMAD < 0.25) sigST = sigMAD/2
      starts[1,,] = getStartsNew(x,J, sigST)
      starts[2,,] = cbind(pi = c(0.99, rep(0.01,J-1)/(J-1)), mu = c(0, sample(x,J-1)), sig = rep(sigST,J))
      starts[3,,] = cbind(pi = c(0.99, rep(0.01,J-1)/(J-1)), mu = c(0, sample(x,J-1)), sig = rep(sigST,J))
      for(i in 4:5){
        starts[i,,1] = rep(1,J)/J
        starts[i,,2] = sample(x, J)
        starts[i,,3] = rep(sigST,J)	
      }
    }
    res = array(0, dim = c(dim(starts)[1], J, 3))
    #	print(starts)
    for(i in 1:dim(starts)[1]){
      m = mixModelFitter(x, J, N, starts[i,,], maxiter, tol, p, muInt, sigInt)
      res[i,,1] = m$pi
      res[i,,2] = m$mu
      res[i,,3] = m$sig
    }
    
    penLogLik = function(params, x, N, p){
      pi = params[,1]
      mu = params[,2]
      sig = params[,3]
      ll = sum(log(dens.mixture(x, list(pi = pi, mu = mu, sig = sig))))
      alpha = N*p
      pen = sum((alpha-1)*log(pi))
      return(ll + pen)
    }
    
    plls = apply(res, 1, penLogLik, x, N, p)
    bestInd = which.max(plls)
    return(list(pi = res[bestInd,,1], mu = res[bestInd,,2], sigma = res[bestInd,,3], data = x))
  }


#---------------------mixFDRfunctions-----------------------------------------------#
`mixFdr` <-
  function(x, J = 3, P = NA, noiseSD = NA, theonull = FALSE, calibrate = FALSE, plots= TRUE, nearlyNull = 0.2, starts = NA, p = NA, maxIter = 1000, tol = 0.001, nocheck = FALSE){
    # try lower sigma threshhold	
    # warn if
    # uncalibrated
    # empirical null is substantially different
    # pi0 is small
    
    # return
    # model parameters
    # effect size estimates
    # fdr estimates
    # plot
    
    # calibration code
    
    muInt = sigInt = NA
    sigma = noiseSD
    n = length(x)	
    if(any(is.na(p))) p = c(1, rep(0,J-1))	
    if(any(is.na(muInt))){
      muInt = matrix(0,J,2)
      muInt[,1] = -Inf
      muInt[,2] = Inf		
    }  	
    if(any(is.na(sigInt)) | any(is.null(dim(sigInt)))){
      left = sigma
      if(is.na(sigma)) left = 0.1
      sigInt = matrix(0,J,2)
      sigInt[,1] = left
      sigInt[,2] = Inf
    }
    if(theonull){
      muInt[1,1] = muInt[1,2] = 0
      sigInt[1,1] = sigInt[1,2] = 1	
    }
    if(is.na(P)){
      P = length(x)/5
      if(!calibrate && !theonull){
        warning("Using uncalibrated default penalization P,  which can give misleading results for empirical nulls. Consider rerunning with calibrate = TRUE and using the resulting penalization","\n")
      }
    }
    if(!is.na(sigma)){
      warning(paste("Assuming known noise noiseSD = ",sigma,". If needed rerun with noiseSD = NA to fit noiseSD."))
      
      if(!theonull){
        warning("Note that using known noiseSD constrains the null to have sd at least noiseSD. If underdispersion is suspected, rerun with noiseSD = NA.")
      }
    }
    
    cat("Fitting preliminary model", "\n")
    m =  mixModelManyStarts(x, J, P, starts, maxIter, tol, p, muInt, sigInt, theonull)
    m$noiseSD = min(median(abs(x-m$mu[1]))/0.6745, m$sig)
    if(calibrate == TRUE){
      cat("Calibrating penalization (slow)", "\n")
      cat("To avoid this step next time, note the calibrated penalization parameter and supply it to the function.", "\n")
      # Calibrate
      P = calibrateP(m, B = 25)
      cat("Calibrated penalization is", P, "\n")
    }
    
    if(is.na(noiseSD)){
      cat("Fitting noiseSD\n")
      if(theonull){ 
        cat("Theoretical null, so using noiseSD = 1\n")
        sigma = 1
      }
      sigma = min(median(abs(x-m$mu[1]))/0.6745, m$sig[1])
    }
    
    cat("Fitting final model","\n")
    # using fitted sigma
    sigInt = matrix(0,J,2)
    sigInt[,1] = sigma
    sigInt[,2] = Inf
    if(theonull) sigInt[1,1] = sigInt[1,2] = 1
    m = mixModelManyStarts(x, J, P, starts, maxIter, tol, p, muInt, sigInt, theonull)
    
    # if pi0 is small (<0.85) throw a warning
    if(m$pi[1]<0.85){
      warning("Null proportion pi0 is small. Consider increasing penalization and/or using an empirical null.")
    }
    
    
    if((theonull | !is.na(noiseSD)) && !nocheck){
      # Check if empirical nulls are significantly different
      # if so, throw a warning
      p = c(1, rep(0,J-1))	
      
      muInt = matrix(0,J,2)
      muInt[,1] = -Inf
      muInt[,2] = Inf		
      
      
      left = sigma/2
      sigInt = matrix(0,J,2)
      sigInt[,1] = left
      sigInt[,2] = Inf
      
      
      mEmp = mixModelManyStarts(x, J, P, NA, maxIter, tol, p, muInt, sigInt, theonull = FALSE)
      sigNew = min(median(abs(x-m$mu[1]))/0.6745, mEmp$sig[1])
      
      left = sigNew
      sigInt = matrix(0,J,2)
      sigInt[,1] = left
      sigInt[,2] = Inf
      
      mEmp = mixModelManyStarts(x, J, P, starts, maxIter, tol, p, muInt, sigInt, theonull = FALSE)
      
      diffr =   log( dnorm(2.5, m$mu[1], m$sig[1]) / dnorm(2.5, mEmp$mu[1], mEmp$sig[1])  ) 
      if(abs(diffr) > abs(log(3))) warning("Using an empirical null with a fitted noiseSD gives a substantially different model. Consider rerunning with theonull = FALSE and noiseSD = NA.")  
    }
    
    # calculate junk and return it
    nullGroups = abs(m$mu - m$mu[1])<=nearlyNull
    fdrX = fdrMixModel(x, m, nullGroups)
    FDRX = tailFDRMixModel(x, m, nullGroups)
    EffSize = effectSize(x, m, sigma)
    
    res = list(pi = m$pi, mu = m$mu, sigma = m$sigma, noiseSD = sigma, converged = m$converged, nIter = m$nIter, fdr = fdrX, FDRTwoSided = FDRX$FDRTwoSided, FDRLeft = FDRX$FDRleft, FDRRight = FDRX$FDRright, effectSize = EffSize[,1], effectPostVar = EffSize[,2], call = match.call(), data = x)
    
    # plot
    if(plots){
      par(mfrow=c(1,3))
      xl = paste(J, "Groups;","pi0 = ",round(sum(m$pi[nullGroups]),3),", mu0 = ", round(m$mu[1],3), ", \n sig0 = ", round(m$sig[1],3), ", noiseSD = ", round(sigma,3))
      
      
      z = x
      mi = min(z)
      ma = max(z)
      s = seq(mi - 1, ma + 1, by = 0.001)
      hist(z, pr = TRUE, br = 50, main = "Mixture Model Fit", xlab = xl, ylab = "Density")
      lines(s, dens.mixture(s,res), lwd = 2)
      phiMat = matDens(s, m$mu, m$sig) * matrix(m$pi, length(s), J, byrow = TRUE)
      nullDens = rowSums(matrix(phiMat[,nullGroups], length(s), sum(nullGroups)))
      altDens = rowSums(matrix(phiMat[,!nullGroups], length(s), sum(!nullGroups)))
      lines(s, nullDens, lwd = 2, col = 2)
      lines(s, altDens, lwd = 2, col = 3)
      
      legend("topleft", legend = c("Mixture", "Null", "Alternative"), lwd = 2, col = c(1,2,3))
      o = order(res$data)
      plot(res$data[o], res$effectSize[o], t = 'l', main = "Effect Size Estimates", xlab = "z", ylab = "delta.hat")
      abline(0,1, lty = 3)
      plot(res$data[o], res$fdr[o], t = 'l', main = "fdr/FDR curves", xlab = "z", ylim = c(0,1), ylab = "fdr/FDR", col = 1)
      lines(res$data[o], res$FDRTwoSided[o], col = 2)
      legend("bottomleft", legend = c("fdr", "FDR 2 Sided"), col = c(1,2), lwd = c(1,1))
    }
    cat("\n")
    cat("Fitted Model: J = ",J," groups\n", sep = "")
    cat("----------------------------\n")
    cat("null?\t", format(nullGroups, digits = 4, justify="left"),"\n\n", sep = "\t")
    cat("pi =\t", format(round(res$pi,4), justify="left"),"\n\n", sep = "\t")
    cat("mu = \t", format(round(res$mu,4), digits = 4, justify="left"),"\n\n", sep = "\t")
    cat("sigma = ", format(round(res$sigma,4), digits = 4, justify="left"),"\n\n", sep = "\t")
    cat("noiseSD = ", format(round(res$noiseSD,4),digits=4, justify="left"),"\n\n\n", sep = "\t")
    
    return(res)	
  }

`mixModelFitter` <-
  function(x, J, N, starts, maxiter, tol, p, muInt, sigInt){
    
    
    EMstep = function(x, pi, mu, sig, N, p, muInt, sigInt){
      n = length(x)
      J = length(mu)
      
      # E-step	
      phiMat = matDens(x,mu,sig)
      denominators = (phiMat %*% pi)[,1]
      G = phiMat * ((1/denominators) %o% pi)
      tG = t(G)
      
      # M-step
      numeratorsMu = (tG %*% x)[,1]
      zplusj = rowSums(tG)
      newMu = (numeratorsMu/zplusj)
      newSig = tG %*% (x*x) / zplusj - newMu^2
      
      newMu = pmax(pmin(muInt[,2],newMu),muInt[,1])
      newSig = pmax(pmin(sigInt[,2], sqrt(newSig)), sigInt[,1])
      newPi = (N*p - 1 + zplusj)/(sum(zplusj) + N - J)
      if(N<0) newPi = zplusj/sum(zplusj)
      newPi = (newPi+0.001)
      newPi = pmax(0,pmin(newPi,1))
      newPi = newPi/sum(newPi)
      if(any(is.na(newPi), is.na(newMu), is.na(newSig))){
        print(pi)
        print(mu)
        print(sig)	
      }
      return(list(pi = newPi, mu = newMu, sigma = newSig)) 
    }
    
    # starting points	
    pi = starts[,1]
    mu = starts[,2]
    sigma = starts[,3]
    # that's it, do the loop	
    iter = 0
    converged = FALSE
    while(!converged && iter<maxiter){
      res = EMstep(x, pi, mu, sigma, N, p, muInt, sigInt)
      distance = sum((pi-res$pi)^2)+sum((mu-res$mu)^2) + sum((sigma-res$sigma)^2)
      if(is.na(distance)){
        print(starts)
        print(muInt)
        print(sigInt)	
      }
      if(distance<tol) converged = TRUE
      iter = iter + 1			
      pi = res$pi
      mu = res$mu
      sigma = res$sigma
    }
    #		e = list(pi = pi, mu = mu, sigma = sigma, converged = converged, nIter = iter, data = x)
    #		fdr = fdrMixModel(x, e, abs(e$mu-e$mu[1])<=nearlyNull)
    return(list(pi = pi, mu = mu, sigma = sigma, converged = converged, nIter = iter))
  }


`matDens` <-
  function(x, mu, sig){
    n = length(x)
    J = length(mu)
    res = matrix(x, n, J) - matrix(mu, n, J, byrow = TRUE)
    ss = matrix(1/sig,n,J,byrow=TRUE)
    res = res * ss # divide each column j by sigma[j]
    res = dnorm(0) * exp(-0.5 * res^2) #stupid pi renaming nonsense
    res = res * ss # divide again
    return(res)
    #	return(matDensC(x,mu,sig))
  }



