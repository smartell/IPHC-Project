# |---------------------------------------------------------------------------|
# | Delay Difference Model
# |---------------------------------------------------------------------------|
setwd("/Users/stevenmartell1/Documents/IPHC/IPHC-project/HTM/")
set.seed(999)



# |---------------------------------------------------------------------------|
# | Model parameters
# |---------------------------------------------------------------------------|
A       <- 4    				#Number of areas
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
R       <- 50.5
g       <- rnorm(A,0,1.0)
P       <- matrix((g),nrow=A,ncol=A,byrow=FALSE)
diag(P) <- diag(P) + (R)
P       <- t(P)/rowSums(P)
# P       <- diag(1,A)



# |----------------------------------------------------------------------------|
# | Equilibrium conditions.
# |----------------------------------------------------------------------------|
apportion <- runif(A,0.0,1.0)
apportion <- apportion/sum(apportion)

be        <- bo * apportion
se        <- exp(rep(s,A) %*% P)/(1+exp(rep(s,A) %*% P))
we        <- -(se*alpha - wk*(se) + wk)/(-1 + (se)*rho)
re        <- -be*(se*alpha+se*rho*we-we)/(we*wk)
ne        <- be / we
ro        <- sum(re)

pe        <- re/sum(re)
so        <- reck * (ro/bo)
beta      <- (reck-1.0)/bo

# wbar      <- (s*alpha+ wk*(1-s)) / (1-rho*s)
# ro        <- bo/wbar*(1-s)
# be        <- bo * apportion
# re        <- ro * apportion 
# se        <- exp(-rep(m,A) %*% P)/(1+exp(-rep(m,A) %*% P))
# we        <- ((se*alpha) + wk*(1-se))/(1-(se*rho))

# be        <- wk*re/(1.-se*(alpha/we+rho))
# be        <- -wk*re*we/(alpha*se+rho*we*se-we)
# bo        <- sum(be)

# print(ne)
# re        <- (be - se*be*(alpha/we+rho))/wk

# ne        <- (be*(-1+se))/(we*(se-1))
# The following is the same as ne=be/we
#ne        <- ((be-wk*re)%*%solve(P)/se - rho*be)/alpha
# print(ne)
# Constraint in that ne >= se*ne
# pe        <- ((ne - se*ne)*wbar)/(bo*(1.-s))


# |----------------------------------------------------------------------------|
# | Dynamic calculations
# |----------------------------------------------------------------------------|
T  <- 100
bt = nt = rt = wt = matrix(NA,nrow=T+1,ncol=A)
bt[1,] = be
nt[1,] = ne
rt[1,] = pe*ro 
wt[1,] = be/ne
for(i in 1:T)
{

	if(i <= agek)
	{
		rt[i,] = pe*ro 
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