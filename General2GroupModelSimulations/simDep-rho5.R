

library(expm)
library(Rcpp)
library(devtools)

source("functions.R")

#----------------------------------------for finding the policy-------------------------#
#The parameters of the simulations
var1 = 1 # noise
var2 = 0.01 # extra variance for non-null
rhovec=rep(0.5, 1000) #  #the symmetric correlation in each independent block
Kvec = rep(5, length(rhovec))  #the block sizes
prob = 0.3 #the probability of being nonnull
alpha=0.05
mu=-1.5 #SR
nb = length(Kvec) #the number of blocks
#Saharon's objects
block_matrices = block_sum_matrices = block_nones = block_detfacs = list()
block_deltas=list()#SR
block_invdeltas=list()#SR1
for (b in 1:nb){
  ## RUN THIS ONCE PER BLOCK
  matrices = sum_matrices = nones = detfacs = list()
  deltas = list() #SR
  invdeltas = list() #SR1
  KMAX = Kvec[b]
  rho = rhovec[b]
  for (K in 1:KMAX){
    covmat_base= matrix(rho*var1, nrow=K,ncol=K) + diag(K)*(1-rho)*var1
    matnow= matrix(nrow=0,ncol=K)
    nonenow=NULL
    detfacnow = NULL
    deltnow = NULL #SR
    invdeltnow = NULL #SR1
    bitmask = 2^(1:K)
    for (i in 0:(2^K-1)){
      fn = (i%%bitmask)>= bitmask/2
      covmat = covmat_base + diag(fn*var2,nrow=K,ncol=K)
      invmat = solve(covmat)
      matnow = rbind(matnow, invmat)
      nonenow=c(nonenow, sum(fn))
      detfacnow = c(detfacnow, det(covmat)^(-0.5)) 
      deltnow = c(deltnow, fn*mu) #SR
      invdeltnow = c(invdeltnow, invmat%*%(fn*mu)) #SR1
    }
    matrices[[K]] = matnow
    nones[[K]]=nonenow
    detfacs[[K]]=detfacnow
    deltas[[K]] = deltnow #SR
    invdeltas[[K]] = invdeltnow #SR
    sum_matrix = matrix(nrow=0,ncol=2^K)
    for (i in 1:K){
      cols = ((0:(2^K-1))%%2^i) >= (2^i) / 2
      sum_matrix = rbind(sum_matrix, as.numeric(cols))
    }
    sum_matrices[[K]] = sum_matrix #locations where hj=1
  }
  block_matrices[[b]] = matrices
  block_sum_matrices[[b]] = sum_matrices
  block_nones[[b]] = nones
  block_detfacs[[b]]=detfacs
  block_deltas[[b]] = deltas #SR
  block_invdeltas[[b]] = invdeltas #SR1
}

#----------find muOMTFDR and muOMTpFDR---------

maxiter =5000 #5000

musFDR= seq( 6209.960-25,  6209.960+25,length.out = 1000)
musmargFDR =   seq(3300.766-240, 3300.766+60,length.out = 1000)
muspFDR=    seq(6209.960-25,  6209.960+25,length.out = 1000)
musmargpFDR =   seq( 3397.557-240, 3397.557+60,length.out = 1000)
musmFDR= seq(    15.33934-0.5,   15.33934+0.5,length.out = 1000)
musmargmFDR =   seq( 27.23658-0.06,27.23658+0.06,length.out = 500)

start.time <- Sys.time()
out =fdep(musmFDR,  musFDR, muspFDR, musmargmFDR,  musmargFDR, musmargpFDR, var1, var2, mu, rhovec, Kvec, prob, block_matrices, block_sum_matrices, block_nones, block_detfacs, block_deltas,block_invdeltas, maxiter=maxiter) #SR 
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

mbest = rep(NA,6)
mbest[1]=which.min(abs(out$lev_mat[1:length(musFDR),1]-alpha))
mbest[2]=which.min(abs(out$lev_mat[1:length(muspFDR),2]-alpha*out$minprob_mat[1:length(muspFDR),2]))
mbest[3]=which.min(abs(out$ev_mat[1:length(musmFDR),3]/out$er_mat[1:length(musmFDR),3]-alpha))
mbest[4]=which.min(abs(out$lev_mat[1:length(musmargFDR),4]-alpha))
mbest[5]=which.min(abs(out$lev_mat[1:length(musmargpFDR),5]-alpha*out$minprob_mat[1:length(musmargpFDR),5]))
mbest[6]=which.min(abs(out$ev_mat[1:length(musmargmFDR),6]/out$er_mat[1:length(musmargmFDR),6]-alpha))





colnames = c("OMTFDR","OMTpFDR","OMTmFDR","marglocfdrFDR","marglocfdrpFDR","marglocfdrmFDR")


for (i in 1:6){
  print(colnames[i])
  print(cbind(mbest[i],  out$lev_mat[mbest[i],i], out$ev_mat[mbest[i],i]/out$er_mat[mbest[i],i], out$pow_mat[mbest[i],i], out$ev_mat[mbest[i],i], out$er_mat[mbest[i],i] ))
}

musFDR[mbest[1]]
musmargFDR[mbest[4]]

muspFDR[mbest[2]]
musmargpFDR[mbest[5]]

musmFDR[mbest[3]]
musmargmFDR[mbest[6]]



save("mbest","out", "musmFDR",  "musFDR", "muspFDR", "musmargmFDR",  "musmargFDR", "musmargpFDR","mu", "var1", "var2", "rhovec","Kvec", "prob", "alpha", 
     "block_matrices","block_sum_matrices","block_nones","block_detfacs","block_deltas",file = "D:\\Saharon\\20-04-21\\simul_rho_p5_K_5000_var2_p01_maxiter_5000_cv_1p5.Rdata")


#----------------------------------The simulation including all competitors. ------------------------
load("D:\\Saharon\\20-04-21\\simul_rho_p5_K_5000_var2_p01_maxiter_5000_cv_1p5.Rdata")



#comparison procedures: 
#the procedures for which we know the power: OMTFDR, OMTmFDR, margFDR, margmFDR, 
#procedures that ignore dependence: wrongOMTFDR with muOMT from rho=0,  adaptivemFDR with oracle mixture components, adaptiveBH with oracle mixture components, and BH. 



muOMTFDR=musFDR[mbest[1]]
mumarglocfdrFDR=musmargFDR[mbest[4]]

muOMTpFDR=muspFDR[mbest[2]]
mumarglocfdrpFDR=musmargpFDR[mbest[5]]

muOMTmFDR=musmFDR[mbest[3]]
mumarglocfdrmFDR=musmargmFDR[mbest[6]]

#the optimal mu if ignore dependence
muindOMTFDR= 3277.127
muindOMTpFDR= 3375.401
muindOMTmFDR =27.41241 

print(c(muOMTFDR, mumarglocfdrFDR))
print(c(muOMTpFDR, mumarglocfdrpFDR))
print(c(muOMTmFDR, mumarglocfdrmFDR))
#comparison procedures: 
#the procedures for which we know the power: OMTFDR, OMTmFDR, margFDR, margmFDR, 
#procedures that ignore dependence: wrongOMTFDR with muOMT from rho=0,  adaptivemFDR with oracle mixture components, adaptiveBH with oracle mixture components, and BH. 


maxiter = 5000
#column order: OMTFDR,   marglocfdrFDR,   wrongOMTFDR XXX
V_mat = matrix(0, nrow =maxiter, ncol = 12 )
R_mat =   matrix(0, nrow =maxiter, ncol = 12 )
#  minprob=lev = pow=ev=maxr=numeric(length(mus))
for (iter in 1:maxiter){
  print(iter)
  out = rbeta(var1, var2, mu, rhovec, Kvec, prob, withh=TRUE) #SR
  blockbeta = out$block_beta
  vec_h=out$vec_h
  #for oracle rule
  out = AandB(blockbeta, Kvec, block_matrices, block_deltas, block_invdeltas, block_sum_matrices,block_nones,block_detfacs) #SR+SR1
  a = out$a
  b=out$b
  #for marginal rule
  marglocfdr = NULL
  # Given beta for block of size K compute the locfdr values
  nb = length(blockbeta)
  for (nbi in 1:nb){
    margprob = (1-prob)*dnorm(as.vector(blockbeta[[nbi]]), mean = 0, sd = sqrt(var1))+prob*dnorm(as.vector(blockbeta[[nbi]]), mean = mu, sd = sqrt(var1+var2)) # SR - need to check I used mu correctly
    marglocfdr =c(marglocfdr, (1-prob)*dnorm(as.vector(blockbeta[[nbi]]), mean = 0, sd = sqrt(var1))/margprob)  
  }
  omarglocfdr = sort(marglocfdr)
  amarg = 1-omarglocfdr
  bmarg =  BZCpp(omarglocfdr)
  #need to get the order right
  #OMTFDR
  Rz = a-muOMTFDR*b
  Dz = DCpp(Rz)
  oo = order(out$locfdr)
  origorderDz = Dz[order(oo)]
  ind=1
  V_mat[iter,ind]=sum((1-vec_h)* origorderDz)
  R_mat[iter,ind]=sum(Dz)
  #marglocfdrFDR
  Rz = amarg-mumarglocfdrFDR*bmarg
  Dz = DCpp(Rz)
  oo = order(marglocfdr)
  origorderDz = Dz[order(oo)]
  ind = 2
  V_mat[iter,ind]=sum((1-vec_h)* origorderDz)
  R_mat[iter,ind]=sum(Dz)
  #muindOMTFDR
  Rz = amarg-muindOMTFDR*bmarg
  Dz = DCpp(Rz)
  origorderDz = Dz[order(oo)]
  ind = 3
  V_mat[iter,ind]=sum((1-vec_h)* origorderDz)
  R_mat[iter,ind]=sum(Dz)
  
  
  
  
  #OMTpFDR
  Rz = a-muOMTpFDR*b
  Rz[1]=Rz[1]+muOMTpFDR*alpha
  Dz = DCpp(Rz)
  oo = order(out$locfdr)
  origorderDz = Dz[order(oo)]
  ind=4
  V_mat[iter,ind]=sum((1-vec_h)* origorderDz)
  R_mat[iter,ind]=sum(Dz)
  #marglocfdrpFDR
  Rz = amarg-mumarglocfdrpFDR*bmarg
  Rz[1]=Rz[1]+mumarglocfdrpFDR*alpha
  Dz = DCpp(Rz)
  oo = order(marglocfdr)
  origorderDz = Dz[order(oo)]
  ind = 5
  V_mat[iter,ind]=sum((1-vec_h)* origorderDz)
  R_mat[iter,ind]=sum(Dz)
  #muindOMTpFDR
  Rz = amarg-muindOMTpFDR*bmarg
  Rz[1]=Rz[1]+muindOMTpFDR*alpha
  Dz = DCpp(Rz)
  origorderDz = Dz[order(oo)]
  ind = 6
  V_mat[iter,ind]=sum((1-vec_h)* origorderDz)
  R_mat[iter,ind]=sum(Dz)
  
  
  
  #mFDR comparisons
  
  
  #OMTmFDR
  Rz = a-muOMTmFDR*(out$olocfdr-alpha)
  Dz = DCpp(Rz)
  oo = order(out$locfdr)
  origorderDz = Dz[order(oo)]
  ind=7
  V_mat[iter,ind]=sum((1-vec_h)* origorderDz)
  R_mat[iter,ind]=sum(Dz)
  
  
  
  #marglocfdrmFDR
  Rz = amarg-mumarglocfdrmFDR*(omarglocfdr-alpha)
  Dz = DCpp(Rz)
  oo = order(marglocfdr)
  origorderDz = Dz[order(oo)]
  ind = 8
  V_mat[iter,ind]=sum((1-vec_h)* origorderDz)
  R_mat[iter,ind]=sum(Dz)
  
  #muindOMTmFDR
  Rz = amarg-muindOMTmFDR*(omarglocfdr-alpha)
  Dz = DCpp(Rz)
  origorderDz = Dz[order(oo)]
  ind = 9
  V_mat[iter,ind]=sum((1-vec_h)* origorderDz)
  R_mat[iter,ind]=sum(Dz)
  
  
  
  #---ADD COMPETITORS: 
  #est-OMT-mFDR
  
  
  fortalpha = which(cumsum(omarglocfdr)/(1:length(omarglocfdr))<=alpha)
  if (length(fortalpha)>0){
    talpha =max(omarglocfdr[fortalpha])
    ind=10
    V_mat[iter,ind]=sum(marglocfdr[vec_h==0]<=talpha)
    R_mat[iter,ind]=sum(marglocfdr<=talpha)
  }  
  #adaptive BH
  padjustedBH = p.adjust(pnorm(unlist(blockbeta)), method = "BH") #left sided p-values. 
  ind = 11
  R_mat[iter,ind] = sum(padjustedBH<=(alpha/(1-prob)))
  V_mat[iter,ind] = sum((padjustedBH[vec_h==0])<=(alpha/(1-prob)))
  
  #BH
  ind = 12
  R_mat[iter,ind] = sum(padjustedBH<=(alpha))
  V_mat[iter,ind] = sum((padjustedBH[vec_h==0])<=(alpha))
  
  
  
  
  
  
}



print(round(apply(V_mat/pmax(R_mat,1), 2, mean),3))
print(round(apply(V_mat/pmax(R_mat,1), 2, sd),3)/sqrt(dim(R_mat)[1]))
print(round(apply(R_mat-V_mat, 2, mean)))


boxplot(V_mat/pmax(R_mat,1))

boxplot(R_mat-V_mat)
boxplot(V_mat)


save("muOMTFDR", "mumarglocfdrFDR","muOMTpFDR", "mumarglocfdrpFDR","muOMTmFDR", "mumarglocfdrmFDR","muindOMTFDR","muindOMTpFDR","muindOMTmFDR","mu", "V_mat", "R_mat",  "var1", "var2", "rhovec","Kvec", "prob", "alpha", file = "D:\\Saharon\\20-04-21\\res_rho_p5_K_5000_var2_p01_maxiter_5000_cv_1p5.Rdata")



#> print(round(apply(V_mat/pmax(R_mat,1), 2, mean),3))
#[1] 0.050 0.052 0.060 0.050 0.054 0.061 0.050 0.050 0.050 0.050 0.050 0.035
#> print(round(apply(V_mat/pmax(R_mat,1), 2, sd),3)/sqrt(dim(R_mat)[1]))
#[1] 0.0001838478 0.0011737973 0.0012445079 0.0001838478 0.0013859293 0.0014142136 0.0001555635 0.0002969848 0.0002969848 0.0002969848 0.0003111270 0.0003394113
#> print(round(apply(R_mat-V_mat, 2, mean)))
#[1] 387 170 198 387 164 186 386 121 121 120 123  73
