# |---------------------------------------------------------------------------|
# | Delay Difference Model
# |---------------------------------------------------------------------------|
set.seed(999)



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
R       <- 2.5
g       <- rnorm(A,0,1.0)
P       <- matrix(exp(g),nrow=A,ncol=A,byrow=FALSE)
diag(P) <- diag(P) + exp(R)
P       <- P/colSums(P)
# P       <- diag(1,A)



# |----------------------------------------------------------------------------|
# | Equilibrium conditions.
# |----------------------------------------------------------------------------|
apportion <- runif(A,1.0,1.0)
apportion <- apportion/sum(apportion)
wbar      <- (s*alpha+ wk*(1-s)) / (1-rho*s)
ro        <- bo/wbar*(1-s)
be        <- bo * apportion 
se        <- rep(s,A)
we        <- ((se*alpha)%*%P + wk*(1-se%*%P))/(1-(se*rho)%*%P)
# we        <- -(se%*%P*alpha - wk*se%*%P + wk)/(-1 + se%*%P*rho)
ne        <- be / we

# ne        <- (be*(-1+se)%*%P)/(we*(se%*%P-1))

pe        <- ((ne - se*ne%*%P)*wbar)/(bo*(1.-s))
so        <- reck * (ro/bo)
beta      <- (reck-1.0)/bo


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

	bt[i+1,] <- s*(alpha*nt[i,]+rho*bt[i,])%*%P + wk*rt[i,]
	nt[i+1,] <- s*nt[i,]%*%P + rt[i,]
	wt[i+1,] <- bt[i+1,]/nt[i+1,]
}
