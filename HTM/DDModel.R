# |---------------------------------------------------------------------------|
# | Delay Difference Model
# |---------------------------------------------------------------------------|
setwd("/Users/stevenmartell1/Documents/IPHC/IPHC-project/HTM/")
set.seed(939)



# |---------------------------------------------------------------------------|
# | Model parameters
# |---------------------------------------------------------------------------|
A       <- 8    				#Number of areas
agek    <- 6
bo      <- 100
reck    <- 12
m       <- 0.15
s       <- exp(-m)
alpha   <- 5.7557
rho     <- 0.9048374
wk      <- alpha * (1-rho^agek)/(1-rho)
 

# |----------------------------------------------------------------------------|
# | Movement matrix based on gravity model
# |----------------------------------------------------------------------------|
R       <- 4.5
g       <- rnorm(A,0,1.0)
P       <- matrix(exp(g),nrow=A,ncol=A,byrow=FALSE)
diag(P) <- diag(P) + exp(R)
P       <- t(P)/rowSums(P)
# P       <- diag(1,A)



# |----------------------------------------------------------------------------|
# | Equilibrium unfished conditions.
# |----------------------------------------------------------------------------|
apportion <- runif(A,0.0,1.0)
apportion <- apportion/sum(apportion)

be        <- bo * apportion
ze        <- rep(log(s/(1-s)),A)
se        <- exp(ze %*% P)/(1+exp(ze %*% P))
we        <- -(se*alpha - wk*(se) + wk)/(-1 + (se)*rho)
re        <- -be*(se*alpha+se*rho*we-we)/(we*wk)
ne        <- be / we
ro        <- sum(re)

pe        <- re/sum(re)
so        <- reck * (ro/bo)
beta      <- (reck-1.0)/bo
aj        <- reck * (re/be)
bj        <- (reck-1.0)/be

# |----------------------------------------------------------------------------|
# | Equilibrium fished conditions.
# |----------------------------------------------------------------------------|
# | Given a vector of fishing mortality rates calculate corresponding 
# | equilibrium conditions, including yields (ye)
ue <- rep(0.5,A)
equil <- function(ue)
{	
	ue <- c(1,1,1,1,0.8,0.8,0.8,0.8)*ue
	ze <- rep(log(s/(1-s)),A)
	se <- exp(ze %*% P)/(1+exp(ze %*% P))*(1-ue)
	we <- -(se*alpha - wk*(se) + wk)/(-1 + (se)*rho)
	be <- -(-we+se*alpha+se*rho*we+wk*aj*we)/(bj*(-we+se*alpha+se*rho*we))
	re <- aj*be/(1.0+bj*be)
	for(i in 1:100)
	{
		re <- pe * (so*sum(be)/(1+beta*sum(be)))
		be <- se * be *(alpha/we+rho)+ wk*re 
	}
	ne <- be / we
	ye <- ue * be
	return(ye)
}

YE <- function(log.ue)
{
	ue <- exp(log.ue)
	ye<- equil(ue)
	return(sum(ye))
}


# |----------------------------------------------------------------------------|
# | Dynamic calculations
# |----------------------------------------------------------------------------|
T  <- 100
bt = nt = rt = wt = matrix(NA,nrow=T+1,ncol=A)
bt[1,] = be
nt[1,] = ne
rt[1,] = re #pe*ro 
wt[1,] = be/ne
for(i in 1:T)
{

	if(i <= agek)
	{
		rt[i,] = re#pe*ro 
	}
	else
	{
		st     <- sum(bt[i-agek,])
		rt[i,] = pe*(so*st/(1.+beta*st))
	}

	bt[i+1,] <- (se)*(alpha*nt[i,]+rho*bt[i,]) + wk*rt[i,]
	nt[i+1,] <- (se)*nt[i,] + rt[i,]
	wt[i+1,] <- bt[i+1,]/nt[i+1,]
}

matplot(bt,type="l",las=1)
fe <- seq(0,0.5,by=0.02)
ye <- t(sapply(fe,equil))
matplot(fe,ye,type="o")