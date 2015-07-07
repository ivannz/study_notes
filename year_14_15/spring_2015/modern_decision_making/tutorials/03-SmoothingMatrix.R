# Comparison of Linear Regression and Smoothing Splines with incresing Lambda

# Synthetic data
n=100
x=runif(n,-5,10); 
x=sort(x)
y=1+.1*x^2-.1*x^4+.01*x^5+rnorm(n,0,20);
plot(x,y)
h=diff(x)

R=matrix(0,n-2,n-2)
for (j in 2:(n-2)){
  i=j-1
  R[i,i]=(h[j-1]+h[j])/3
  R[i,i+1]=h[j]/6
  R[i+1,i]=h[j]/6
}

Q=matrix(0,n,n-2)
for (j in 2:(n-1)){
  i=j-1
  Q[j-1,i]=1/h[j-1]
  Q[j,i]=-(1/h[j-1]+1/h[j])
  Q[j+1,i]=1/h[j]
}


# SS with Lambda=1
lambda=1;
K=Q%*%solve(R)%*%t(Q);
S1=solve(diag(n)+lambda*K)
eeS1=eigen(S1)
eval1=eeS1$values
evec1=eeS1$vectors

yh1=S1%*%y
lines(x,yh1,col="blue",lwd=2)

# SS with Lambda=1000
lambda=1000;
K=Q%*%solve(R)%*%t(Q);
S1000=solve(diag(n)+lambda*K)
eeS1000=eigen(S1000)
eval1000=eeS1000$values
evec1000=eeS1000$vectors

yh1000=S1000%*%y
lines(x,yh1000,col="purple",lwd=2)

# LR
X=cbind(rep(1,n),x)
H=X%*%solve(t(X)%*%X)%*%t(X)
eeH=eigen(H)
evalH=eeH$values
evecH=eeH$vectors

yh=H%*%y
lines(x,yh,col="magenta",lwd=2)


# Plot eigenvalues
par(mfrow=c(1,3), mar=c(2,2,2,1))
plot(eval1[1:10], pch=19, main="lambda=1")
plot(eval1000[1:10], pch=19, col="magenta", main="lambda=1000")
plot(evalH[1:10], pch=19, col="red", main="linear reg")
dev.off()

# Plot eigenvectors
par(mfrow=c(3,4),mar=c(1,2,1,1))
plot(x,evec1[,2])
plot(x,evec1[,4])
plot(x,evec1[,8])
plot(x,evec1[,10])

plot(x,evec1000[,2])
plot(x,evec1000[,4])
plot(x,evec1000[,8])
plot(x,evec1000[,10])

plot(x,evecH[,2])
plot(x,evecH[,4])
plot(x,evecH[,8])
plot(x,evecH[,10])
dev.off()

