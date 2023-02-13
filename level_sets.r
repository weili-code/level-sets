# This script is for level curve in regression using slice sampling to locate level curves
# The script used the splines2 package, dbs() to compute derivatives of B-splines

# Author: Wei Li
# 
# ref: Posterior Contraction and Credible Regions for Level Sets
#      Li, W. and Ghosal, S., Electronic Journal of Statistics, 2021.

library("splines")
library("splines2")
#library(rgl)
library(MASS)
library(pracma)
library("gtools")
library(parallel)

source("slice_sampling_2.r")


####################   simulate data ###################

######### -----------dgp4----------------------------------------------------
set.seed(135)
N<-2000
X1<-runif(N,-1,1)
X2<-runif(N,-1,1)
X<-cbind(X1,X2)
#plot(X)

U0 = seq(0, (2)*pi, len = N)   # U0 = seq(0, 2*pi, len = N)  # full circle data  
X_true<-cbind(1/2 * cos(U0), 1/2 * sin(U0))
#plot(X_true)

f_0<-function(x){
	x1<-x[1]
	x2<-x[2]
	u<-sqrt(x1^2+x2^2)
	return( 1+dnorm(u,1/2,.3)^( 1+(cos( atan(x2/x1)))^2    ) )
}

X_list<-lapply(seq_len(nrow(X)), function(i) X[i,])  # create a list consists of row vectors of grid_mat
sigma_0<-.1  # true standard deviation
y<- unname( sapply( X_list, f_0   ))+rnorm(N,0,sd=sigma_0)
#plot3d(X[,1],X[,2],y,col="red", size=3)
nobs<-N
dat<-cbind(X,y)



###### define a function ######

Bdata<-function(nobs,x1bs,x2bs,J1,J2){
  out<-matrix(,0,J1*J2)
  for (i in 1:nobs ){
    out<-rbind(out, kronecker(x1bs[i,], x2bs[i,]))
  } 
  return(out)
}

#########################################
# find the optimal choice of number of J's
#########################################

###############################################################
# parameters (optimal choice)
optJ<-c(9,9)
q1=5; k1=optJ[1]-q1   # k1 is number of interior knots, q1=order
J1=q1+k1;

# parameters (optimal choice)
q2=5; k2=optJ[2]-q2   # k2 is number of interior knots, q2=order
J2=q2+k2;

knots1<-seq(-1.1,1.1,length=k1+2)
bdknots1<-knots1[c(1,length(knots1))]
inknots1<-knots1[2:(length(knots1)-1)]
knots2<-seq(-1.1,1.1,length=k2+2)
bdknots2<-knots2[c(1,length(knots2))]
inknots2<-knots2[2:(length(knots2)-1)]

x1<-dat[,1]
x2<-dat[,2]
x1bs<-bs(x1,knots=inknots1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)
x2bs<-bs(x2,knots=inknots2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)

B<-Bdata(nobs,x1bs,x2bs,J1,J2)  # the capital B matrix which is n by J1*J2


##############################################################################
## constructing function, derivatives, Hessian 
##############################################################################

############# regression function ####################
f_theta<-function(x1,x2,theta) {
# this function evaluates the B-spline approximation given by theta, at point (x1,x2)
# input: x1, x2 should be scalars,theta coefficients ordered in dictionary order
	res1<-bs(x1, knots=inknots1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
	res2<-bs(x2, knots=inknots2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
	#res1<-bbase2(x1, min(x_test),max(x_test),k1+1,q1-1)
	#res2<-bbase2(x2, min(x_test),max(x_test),k2+1,q2-1) 
	out1<- kronecker(res1, res2)%*%as.vector(theta)
	output<-as.numeric(out1)
	return(output)
}


f10<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
	temp1<-dbs(x1, knots=inknots1, derivs=1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
    temp2<-bs(x2, knots=inknots2, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))
 
}	
	
f01<-function(x1,x2,index=NULL,theta) {
#input: x1,x2 are coordinates on plane,
#index is the grid index a 2-vector, say index=c(i,j) means the x1,x2 is in the interval [t_i,t_{i+1}]x[t_j,t_{j+1}]  
    
    temp1<-bs(x1, knots=inknots1, Boundary.knots=bdknots1, degree=q1-1, intercept=TRUE)[1,] 
    temp2<-dbs(x2, knots=inknots2, derivs=1, Boundary.knots=bdknots2, degree=q2-1, intercept=TRUE)[1,] 
    out<- kronecker(temp1, temp2)%*%as.vector(theta)
    return(as.numeric(out))
}




################### priors ###########################

theta_0<-matrix(rep(0,J1*J2))
Lambda_0<-diag(1,J1*J2)
# empirical bayes estiamte for sigma2
sigma2_hat<-nobs^{-1}*t(y-B%*%theta_0)%*%solve(B%*%Lambda_0%*%t(B)+diag(1,nobs))%*%(y-B%*%theta_0)

################## posterior sampling ################

# posterior samples for theta --------------------------------------------
a<-6
b<-2
nsim<-100  # number of posterior samples
sigma2.post<-1/rgamma(nsim, (a+nobs)/2, rate=(b+nobs*sigma2_hat)/2 ) # 1*nsim vector of posterior sample for sigma2
sigma2_hat # as alternative for sigma2.post

BB_lambda<- (solve(Lambda_0)+t(B)%*%B)
BB_lambda_inv<-solve(BB_lambda)
theta_post_EXPmean<- BB_lambda_inv%*%(t(B)%*%y+solve(Lambda_0)%*%theta_0)

E<-matrix(rnorm(nsim*J1*J2,0,sqrt(sigma2.post)), nsim, J1*J2,byrow=FALSE ) 
# E is nsim by J1J2 matrix. each row corresponds to J1*J2 iid simulated normals with sd corresponding to one element of sigma2.post

#using empirical sigma2_hat. E is nsim by J1J2 matrix. each row corresponds to J1*J2 iid simulated normals with sd=sigma2_hat
E<-matrix(rnorm(nsim*J1*J2,0,sqrt(sigma2_hat)), nsim, J1*J2,byrow=FALSE ) 

chol_BB_lambda_inv<-chol(BB_lambda_inv)
theta.post<-t(t(chol_BB_lambda_inv)%*%t(E)+ matrix(theta_post_EXPmean,ncol=nsim, nrow=J1*J2, byrow=FALSE))  # nsim*J1J2 matrix of posterior sample for theta
theta_postmean<-as.matrix(colMeans(theta.post)) #monte carlo mean of the regression function
#thetaMat_postmean<-matrix(theta_postmean,nrow=J1,byrow = TRUE)  # matrix of theta, rows=first index, columns=second index
#as.vector(t(thetaMat_postmean))-theta_postmean=0



#########################################################################
#----------- compute the 1-gamma quantile for GP(0,Sigma) --------------

gamma <- 0.1
stepl<-.02  #.02 gives about 100 points in x1s
x1s<-seq(bdknots1[1],bdknots1[2],by=stepl) 
x2s<-seq(bdknots2[1],bdknots2[2],by=stepl) 

nsim2<-nsim
#nsim2<-100  # number of iteration to compute quantiles.
E_test<-matrix(rnorm(nsim2*J1*J2,0,sqrt(sigma2_hat)), nsim2, J1*J2,byrow=FALSE ) 
MAT<-t(chol_BB_lambda_inv)%*%t(E_test) #a J1J2 by nsim2 matrix. MAT gives nsim2 centered theta from posterior 

compute_quantile<-function(j){
  GP_00_mat<-matrix(, length(x1s),length(x2s)) 
  
  theta_demean<- MAT[,j]
  f_minus_mean<-function(x1,x2){
    out<-f_theta(x1,x2,theta=theta_demean)
    return(abs(out))
  }
  
  GP_00_mat<-outer(x1s,x2s,Vectorize(f_minus_mean))
  
  return(GP_00_mat=GP_00_mat)
}


GP_00_list<-mclapply(1:nsim2,compute_quantile)

#GP_00_list<-lapply(results1, function (x) x[[c("GP_20_mat")]])  # note that results1[[1]][[c("GP_20_mat")]] can subset what we need
#R00<-quantile(sapply(GP_00_list,max) , 1 -gamma ) 

##--------- do a further gradient ascent/descent on each of grid search candidates--------------------

g00<-function(x,index=NULL,theta){
# gradient functions
	x1<-x[1]
	x2<-x[2]
	c(f10(x1,x2,index=NULL,theta),f01(x1,x2,index=NULL,theta))
}	
	
f00_temp<-function(x,index=NULL,theta){
	x1<-x[1]
	x2<-x[2]
	c(f_theta(x1,x2,theta))
}

refine_gd<-function(j,GP_list, fn, gr ){
  # j is index for the j-th simulation 
  # GP_list= the list of demean GP values evaluted on grids
  # fn= the objective function to maxmize/minimize
  # gr= gradient 
	theta_demean<- MAT[,j]
	GP_list<-GP_list
	max_candi<-max(GP_list[[j]])
	indx<-which(GP_list[[j]] == max_candi, arr.ind = TRUE) # the index of the max point from grid search result
	x0<-c(x1s[indx[1]], x1s[indx[2]])
	x0
	box_y<-c( max(x0[2]-stepl,bdknots2[1] ), min(x0[2]+stepl, bdknots2[2]) ) # the lower bound/ upper bound in y coordiante center about x0
	box_x<-c( max(x0[1]-stepl,bdknots1[1] ), min(x0[1]+stepl, bdknots1[2] ) )  # the lower bound/ upper bound in x coordiante center about x0
	local_min<-optim(x0, fn=fn, gr=gr, method = "L-BFGS-B",lower = c(box_x[1],box_y[1] ), upper = c(box_x[2], box_y[2]),index=NULL,theta=theta_demean, control = list(maxit = 2000)) 
	local_max<-optim(x0, fn=fn, gr=gr, method = "L-BFGS-B",lower = c(box_x[1],box_y[1] ), upper = c(box_x[2], box_y[2]),index=NULL,theta=theta_demean,control = list(maxit = 2000,fnscale=-1)) 
	if ( local_min$convergence==0 ){
		max_candi<-max(max_candi,abs(local_min$value) )
		} 
	if ( local_max$convergence==0 ){
		max_candi<-max(max_candi,abs(local_max$value) )
		} 	
	return(max_candi)
}

results<-mclapply(1:nsim2,refine_gd, GP_list=GP_00_list, fn=f00_temp, gr=g00)
results<-do.call(c, results)
R00<-quantile( results , 1 -gamma ) 
R00

rho<-1.2
theta.post_select<-theta.post[results<R00*rho,]

######################################################
#### plot the posterior mean regression function #####
######################################################
x1s<-seq(bdknots1[1],bdknots1[2],length.out=100)
x2s<-seq(bdknots2[1],bdknots2[2],length.out=100)

f_post<-function(x1,x2){
	f_theta(x1,x2,as.vector(theta_postmean))
}
fz_post<-outer(x1s,x2s,Vectorize(f_post))
#write.table(fdenz_post, file="fdenz_post.txt", row.names=FALSE, col.names=FALSE)
#fdenz_post<- matrix(scan("fdenz_post.txt"), nrow=1000, byrow=TRUE)

#contour(x1s,x2s,fz_post,xlab="x1",ylab="x2")
#persp3d(x1s,x2s,fz_post,col="skyblue",aspect = c(1, 1, 1))


######################################################
########  slice sampling #############################
######################################################

#Boltzmann distribution function

Boltzmann_true<-function(k, x1,x2, T, c_level){
	if (k==1){
		ret<- function(x1){
			exp(- (  f_0(c(x1,x2)) - c_level)^2/T )
		}
	}
	if (k==2){
		ret<- function(x2){
			exp(-  (   f_0(c(x1,x2)) - c_level)^2/T )
		}	
	}
	return(ret)
}


Boltzmann<-function(k, x1,x2, theta, T, c_level){
	if (k==1){
		ret<- function(x1){
			exp(- (f_theta(x1,x2,as.vector(theta))- c_level)^2/T )
		}
	}
	if (k==2){
		ret<- function(x2){
			exp(- (f_theta(x1,x2,as.vector(theta))- c_level)^2/T )
		}	
	}
	return(ret)
}

################### Find the true level set of the true function #######################

x_level_true<-matrix(,0,2)

# method 1: find level curve using slice sampling using observations values as starting points

T<-1e-5
c_level<-2.1
niter<-100
tol<-.5

x1s0<-seq(bdknots1[1],bdknots1[2],length.out=1000)
x2s0<-seq(bdknots2[1],bdknots2[2],length.out=1000)
set.seed(135)
#plot(1, type="n", xlab="x1", ylab="x2", xlim=c( min(x1s0), max(x1s0)), ylim=c(min(x2s0), max(x2s0) ))  
for (i in 1:N){
	x1_iter<-X[i,1]
	x2_iter<-X[i,2]
	if(abs(f_0(X[i,]))>tol){
		for (j in 1:niter){
		  #print(paste0( " At iteration ", j ))
		  
		  x1.tmp<-Boltzmann_true(1, x1_iter,x2_iter,  T, c_level)
		  x1_iter<- uni.slice_2(x1_iter, f = x1.tmp, w=1,lower=bdknots1[1], upper=bdknots1[2])

		  x2.tmp<-Boltzmann_true(2, x1_iter,x2_iter, T, c_level)
		  x2_iter<- uni.slice_2(x2_iter, f = x2.tmp, w=1,lower=bdknots2[1], upper=bdknots2[2])
			

		}
		  x_iter<-c(x1_iter, x2_iter)	
		  x_level_true<-rbind(x_level_true, (x_iter))	 
		  
		#points(x_iter[1], x_iter[2], xlab="x1",ylab="x2",col="red",pch=20,cex=.2)  	
	 }
	
}


################### Find and plot the estimated level curve #######################
x_level<-matrix(,0,2)

# method 1: find level curve using slice sampling using observations values as starting points

levels_all_fun<-function(niter,c_level, T, tol, seed){
pdf("J9_uq.pdf",width=6.8,height=6.8)
    levels_all<-list()
	  nsim<-nrow(theta.post_select)
    theta.post_1<-rbind(theta.post_select,t(theta_postmean))
    #contour(x1s,x2s,fz_post, xlab=expression(X[1]), ylab=expression(X[2]))
    plot(1, type="n",  xlab=expression(X[1]), ylab=expression(X[2]), xlim=c( min(x1s), max(x1s)), ylim=c(min(x2s), max(x2s) ))  
    for (k in 1:(nsim+1)){
          theta<-theta.post_1[k,]
          out<-matrix(,0,2)    
          set.seed(seed)
          for (i in 1:N){
            	x1_iter<-X[i,1]
            	x2_iter<-X[i,2]
          	  if(abs(f_theta(x1_iter,x2_iter,theta_postmean ))>tol){
          	        
                		for (j in 1:niter){
                		
                		  # print(paste0( " At iteration ", j ))
                		  
                		  x1.tmp<-Boltzmann(1, x1_iter,x2_iter, theta, T, c_level)
                		  x1_iter<- uni.slice_2(x1_iter, f = x1.tmp, w=1,lower=bdknots1[1], upper=bdknots1[2])
                
                		  x2.tmp<-Boltzmann(2, x1_iter,x2_iter, theta, T, c_level)
                		  x2_iter<- uni.slice_2(x2_iter, f = x2.tmp, w=1,lower=bdknots2[1], upper=bdknots2[2])
                			
                		}
                		  x_iter<-c(x1_iter, x2_iter)	
                		  
                		  if (k==(nsim+1)){ # the filament using posterior mean of parameters
                      points(x_iter[1],x_iter[2] ,xlab="x1",ylab="x2",,col="blue",pch=20,cex=.15)
                		  #points(x_level[i,1],x_level[i,2] ,xlab="x1",ylab="x2",col="blue",pch=4,cex=.3)    
                        } else {
                       points(x_iter[1],x_iter[2] ,xlab="x1",ylab="x2",col=k+2,lwd=.1,lty=2,pch=20,cex=.1)
                      }
                		  
                		  out<-rbind(out, (x_iter))	 
          	   }
          }# end of i  	
          
       levels_all[[k]]<-out
    }#end of k    
    
    for (i in 1:nrow(x_level_true)){
          		  points(x_level_true[i,1],x_level_true[i,2] ,xlab="x1",ylab="x2",col="orange",,lwd=1.3,lty=2,pch=20,cex=.15)
    }
    
dev.off()    
return(levels_all)    
}

# credibility in terms of percentage
nrow(theta.post_select)/nsim2

levels_all_dat<-levels_all_fun(niter=100,c_level=2.1, T=1e-5, tol=.5, seed=135)

