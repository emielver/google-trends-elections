install_github("gabrielrvsc/HDeconometrics")
library(gtrendsR)
library(gcdnet)
library(msaenet)
library(GoFKernel) 
library(Metrics)
library(devtools)
library(HDeconometrics)
library(glmnet)
library(ade4)
library(hdm)
library(covTest)
library(PoSI)

# Dates are 2018-11-06, 2016-11-08, 2012-11-06, 2010-11-02
years = c(2018, 2016, 2012, 2010)
month = "11"
days = c("06", "08", "06", "02")
# ARIZONA, CALIFORNIA, CONNECTICUT, FLORIDA, HAWAII, INDIANA, MARYLAND, MISSOURI, NEVADA, NEW YORK, NORTH DAKOTA, OHIO, PENNSYLVANIA, UTAH, VERMONT, WASHINGTON, WISCONSIN
states = c("US-AZ", "US-CA", "US-CT", "US-FL", "US-HI", "US-IN", "US-MD", "US-MO", "US-NV", "US-NY", "US-ND", "US-OH", "US-PA", "US-UT", "US-VT", "US-WA", "US-WI")
basic.searches = c("unemployment", "nra", "mortgage", "fox news", "msnbc news", "social security", "medicaid", "medicare", "establishment", 
                   "truck", "immigration", "beer", "nfl", "climate change", "bible", "terrorism", "middle east", "israel", "border control", 
                   "mexico", "saudi arabia", "dnc", "gop", "socialism", "supreme court", "constitution", "abortion", "kkk", 
                   "access to education", "college application", "trade agreement", "amazon", "homeland security", "ged", "payday loans", "media")

l <- matrix(data = NA, nrow = length(states) * length(years), ncol = 5 * length(basic.searches))
row = 1
col = 0
# for every state
for (state in states) {
  # for every year
  for (i in 1:length(years)) {
    # for every search term
    for (term in basic.searches) {
      vec = gtrends(keyword = term, geo = state, time = paste(years[i]-1, "-", month, "-", days[i], " ", years[i], "-", month, "-", days[i], sep = ""))$interest_over_time$hits[48:52]
      if (length(vec) == 0) {
        for (j in 1:5) {
          l[row, col + j] <- -1
        }
      }
      else {
        for (j in 1:5) {
          l[row, col + j] <-  vec[j]
        }
      }
      col = col + 5
    }
    row = row + 1
    col = 0
  }
}


# read in and initialise y
#datatable <- read.table("Election results.csv",sep=";", header=TRUE)
#y <- NULL
#for (i in 1:17) {
#  y <- c(y,datatable[i,2], datatable[i,3],datatable[i,4],datatable[i,5])
#}
# find the mean of y for each state
mean <- vector(mode= "numeric", length = 17)
counter <- 1
for (i in seq(1,68,4)) {
  mean[counter] <- (y[i+1] + y[i+2] +y[i+3])/3
  counter = counter + 1
}

#### DEFINE VARIABLES

N <- 17
T <- 4
p <- 180
rows <- N * T

### demeaned y
ydot <- vector(mode = "numeric", length = 68)
counter <- 0
for (i in 1:68) {
  if (i%%4 == 1) {
    counter = counter+1
    ydot[i] <- y[i]
  }
  else {
    ydot[i] <- (y[i] - mean[counter])
  }
}

ydot <- ydot*100
mean <- mean*100
### 
for (i in 1:rows){
  for (j in 1:p){
    if (l[i,j] == "<1") {
      l[i,j] = "0.5"
    }
  }
}

lnum <- matrix(data = 0, nrow = length(states) * length(years), ncol = 5 * length(basic.searches))
for (i in 1:rows){
  for (j in 1:p){
    lnum[i,j] <- as.numeric(l[i,j])
  }
}
####
xdot <- matrix(data = 0, nrow = length(states) * length(years), ncol = 5 * length(basic.searches))
#get xdot so demeaned x
i <- 1
while (i < rows){ #for each state
  for (c in 1:p){ #for each variable
    mean.s <- 0
    for(s in 1:(T-1)){
      mean.s = mean.s + lnum[(i+s),c]
    }
    mean.s = mean.s / 3
    for (s in 0:(T-1)){
      xdot[i+s,c] = lnum[i+s,c] - mean.s
    }
  }
  i=i+4
}

# getting lambda
gamma <- 0.1/log(min(p, rows))
x<- 1-(gamma/(2*p))
c <- 1.1
## inverse of cdf is just quantile function so qnorm
lambda <- 2* c * sqrt(rows) * qnorm(x)
## In and out of sample y
y.out <- rep(0,17)
y.in <- rep(0,51)
outIndex <- 1
inIndex <- 1
for (i in 1:68){
  if (i %% 4 == 1) {
    y.out[outIndex] <- ydot[i]
    outIndex = outIndex + 1
  }
  else {
    y.in[inIndex] <- ydot[i]
    inIndex = inIndex + 1
  }
}
## In and out of sample x
x.out <- matrix(data=0, nrow = 17, ncol=p)
x.in <- matrix(data=0, nrow = 51, ncol=p)
outIndex <- 1
inIndex <- 1
for (i in 1:68){
  if (i %% 4 == 1) {
    x.out[outIndex,] <- xdot[i,]
    outIndex = outIndex + 1
  }
  else {
    x.in[inIndex,] <- xdot[i,]
    inIndex = inIndex + 1
  }
}
## LASSO
lasso <- ic.glmnet(x.in, y.in, crit="bic")
plot(lasso$glmnet, "lambda")
plot(lasso)
## Forecasting
pred.lasso = predict(lasso, newdata=x.out)
pred.lasso <- pred.lasso + mean
plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n")
lines(pred.lasso, col=2, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)

lasso.coef <- lasso$coefficients

for (i in 1:181) {
  if (lasso.coef[i] != 0) {
    index <- ((i-1)/5) + 1
    print(basic.searches[index])
  }
}

weights <- vector(mode = "numeric", length = 180)

for (i in 2:181) {
  if (lasso.coef[i] != 0) {
    weights[i-1] = (1/abs(lasso.coef[i]))^2
  }
  else {
    weights[i-1] = (1/0.000001)^2
  }
}

adalasso.1 = ic.glmnet(x.in, y.in, crit="bic", penalty.factor = weights)
pred.adalasso.1 = predict(adalasso.1, newdata = x.out)

ada.1.coef <- adalasso.1$coefficients

for (i in 1:181) {
  if (ada.1.coef[i] != 0) {
    index <- ((i-1)/5) + 1
    print(basic.searches[index])
  }
}
## Adaptive Lasso

# Initial Phi
phi <- c()
for(j in 1:p){ #for each variable
  phi[j] = 0
  i=1
  while(i < length(x.in[,1])){
    for(t in 0:(T-2)){
      for(tprime in 0:(T-2)){
        
        phi[j] = phi[j] + (x.in[i+t,j] * x.in[i+tprime,j]*y.in[i+t] * y.in[i+tprime])
      }
    }
    i = i + 3
  }
  phi[j] = phi[j] / length(x.in[,1])
  phi[j] = sqrt(phi[j])
}
tau = 0.5
first.step.coef=coef(lasso)[-1]
penalty.factor = abs(first.step.coef)^(-tau)
adalasso = ic.glmnet(x.in, y.in, crit="bic", penalty.factor=phi)
pred.adalasso=predict(adalasso,newdata = x.out)
pred.adalasso <- pred.adalasso + mean
pred.adalasso.1 <- pred.adalasso.1 + mean

plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0, 100), main="Lasso Prediction")
lines(pred.lasso, col=3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(14, 10, legend=c("Actual", "Lasso"),col=c("blue", "green"), lty=1, cex=0.8)
Metrics::mse(y.out, pred.lasso)
Metrics::mse(y.out, pred.adalasso)
Metrics::mse(y.out, pred.adalasso.1)



# Refined Phi
eps.in <- as.vector(y.in - x.in %*% coef(adalasso)[-1])
ref.phi <- c()
for(j in 1:p){ #for each variable
  ref.phi[j] = 0
  i=1
  while(i < length(x.in[,1])){
    for(t in 0:(T-2)){
      for(tprime in 0:(T-2)){
        
        ref.phi[j] = ref.phi[j] + (x.in[i+t,j] * x.in[i+tprime,j]*eps.in[i+t] * eps.in[i+tprime])
      }
    }
    i = i + 3
  }
  ref.phi[j] = ref.phi[j] / length(x.in[,1])
  ref.phi[j] = sqrt(ref.phi[j])
}
adalasso = ic.glmnet(x.in, y.in, crit="bic", penalty.factor=ref.phi)
pred.adalasso=predict(adalasso,newdata = x.out)
pred.adalasso <- pred.adalasso + mean

clust.coef <- adalasso$coefficients
ada.coef <- adalasso.1$coefficients
for (i in 2:181) {
  if (ada.coef[i] != 0) {
    index <- ((i-1)/5) + 1
    print(basic.searches[index])
    print(ada.coef[i])
  }
}

Metrics::mse(y.out, pred.adalasso)
plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Lasso Prediction")
lines(pred.lasso, col=3, type = "b")
lines(mean, col = 2, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(14, 10, legend=c("Actual", "Lasso", "Mean"),col=c("blue", "green", "red"), lty=1, cex=0.8)


plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Adaptive Lasso Prediction")
lines(pred.adalasso.1, col=3, type = "b")
lines(mean, col = 2, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(14, 10, legend=c("Actual", "Adaptive Lasso", "Mean"),col=c("blue", "green", "red"), lty=1, cex=0.8)



plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Cluster Lasso Prediction")
lines(pred.adalasso, col = 3, type = "b")
lines(mean, col = 2, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(14, 10, legend=c("Actual", "Cluster Lasso", "Mean"),col=c("blue", "green", "red"), lty=1, cex=0.8)

#### Deviation from Mean ####

plot(y.out - mean, type="b", col=4, xlab = "", ylab = "Deviation from Mean", xaxt = "n", ylim = c(-20,20), main = "Lasso Prediction")
lines(pred.lasso - mean, col=3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 0, b = 0, col = 1)
legend(14, -15, legend=c("Actual", "Lasso"),col=c("blue", "green"), lty=1, cex=0.8)


plot(y.out - mean, type="b", col=4, xlab = "", ylab = "Deviation from Mean", xaxt = "n", ylim = c(-20,20), main = "Adaptive Lasso Prediction")
lines(pred.adalasso.1 - mean, col=3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 0, b = 0, col = 1)
legend(14, -15, legend=c("Actual", "Adaptive Lasso"),col=c("blue", "green"), lty=1, cex=0.8)



plot(y.out - mean, type="b", col=4, xlab = "", ylab = "Deviation from Mean", xaxt = "n", ylim = c(-20,20), main = "Cluster Lasso Prediction")
lines(pred.adalasso - mean, col = 3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 0, b = 0, col = 1)
legend(14, -15, legend=c("Actual", "Cluster Lasso"),col=c("blue", "green"), lty=1, cex=0.8)



Metrics::mse(y.out, pred.lasso)
Metrics::mse(y.out, pred.adalasso)
Metrics::mse(y.out, pred.adalasso.1)

## Post Adaptive Lasso
# Count number of non-zero coefficients
count <- 0
coefs <- ada.coef
for (i in 2:length(coefs)){
  if (coefs[i] != 0) {
    count = count + 1
  }
}
print(count)
# for every non-zero coefficient, add them to matrix
nonZero.in <- matrix(0, nrow = 51, ncol = count)
nonZero.out <- matrix(0, nrow = 17, ncol = count)
zeroIndex <- 1
for (i in 2:length(coefs)) {
  if (coefs[i] != 0) {
    nonZero.in[,zeroIndex] <- x.in[,i]
    nonZero.out[,zeroIndex] <- x.out[,i]
    zeroIndex = zeroIndex + 1
  }
}
betas <- solve(t(nonZero.in)%*%nonZero.in) %*% t(nonZero.in) %*% y.in

pred.post <- nonZero.out %*% betas
pred.post <- pred.post + mean

## Prediction Plot ##
plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Post Adaptive Lasso Prediction")
lines(pred.post, col=3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(13, 15, legend=c("Actual", "Post Adaptive Lasso"),col=c("blue", "green"), lty=1, cex=0.8)
## Demeaned ##
plot(y.out - mean, type="b", col=4, xlab = "", ylab = "Deviation from Mean", xaxt = "n", ylim = c(-25,25), main = "Post Adaptive Lasso Prediction")
lines(pred.post - mean, col=3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 0, b = 0, col = 1)
legend(13, -15, legend=c("Actual", "Post Adaptive Lasso"),col=c("blue", "green"), lty=1, cex=0.8)


Metrics::mse(y.out, pred.post)

## With Mean ##
plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Post Adaptive Lasso Prediction")
lines(pred.post, col=3, type = "b")
lines(mean, col = 2, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(13, 15, legend=c("Actual", "Post Adaptive Lasso", "Mean"),col=c("blue", "green", "red"), lty=1, cex=0.8)


## Post Cluster Lasso
# Count number of non-zero coefficients
count <- 0
coefs.ada <- coef(adalasso)
for (i in 2:length(coefs)){
  if (coefs.ada[i] != 0) {
    count = count + 1
  }
}
# for every non-zero coefficient, add them to matrix
nonZero.in <- matrix(0, nrow = 51, ncol = count)
nonZero.out <- matrix(0, nrow = 17, ncol = count)
zeroIndex <- 1
for (i in 2:length(coefs.ada)) {
  if (coefs.ada[i] != 0) {
    nonZero.in[,zeroIndex] <- x.in[,i]
    nonZero.out[,zeroIndex] <- x.out[,i]
    zeroIndex = zeroIndex + 1
  }
}

betas <- solve(t(nonZero.in)%*%nonZero.in) %*% t(nonZero.in) %*% y.in

pred.post.clust <- nonZero.out %*% betas
pred.post.clust <- pred.post.clust + mean

plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Post Cluster Lasso Prediction")
lines(pred.post.clust-mean, col=3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(14, 10, legend=c("Actual", "Post Cluster Lasso"),col=c("blue", "green"), lty=1, cex=0.8)

Metrics::mse(y.out, pred.post.clust - mean)


plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Mean Prediction")
lines(mean, col = 3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(13, 20, legend=c("Actual", "Mean"),col=c("blue", "green"), lty=1, cex=0.8)
Metrics::mse(y.out, mean)
Metrics::mse(y.in + mean, mean.in)
Metrics::mse(y.out, pred.adalasso.1)

###### Adaptive Elastic Net
ada.elas.net <- gcdnet(x.in, y.in, nlambda = 10000, method = "ls", lambda2 = 2, pf = weights, pf2 = weights)
pred.ada.elas <- predict(ada.elas.net, x.out)
pred.ada.elas <- pred.ada.elas + mean
plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Adaptive Elastic Net Prediction")
lines(pred.ada.elas[,99], col = 3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(15, 10, legend=c("Actual", "AEN"),col=c("blue", "green"), lty=1, cex=0.8)
Metrics::mse(y.out, pred.ada.elas[,99])


fitted.elas.net = predict(ada.elas.net, x.in)
plot(y.in, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(-20,20), main = "In Sample Lasso Estimation")
lines(fitted.elas.net[,99], col = 3, type= "b")

###### Adaptive Elastic Net 2
ada.net <- aenet(x.in, y.in, family = "gaussian", init = "ridge", tune = "bic", scale = 2)
pred.ada.net <- predict(ada.net, x.out)
pred.ada.net <- pred.ada.net + mean

plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Adaptive Elastic Net Prediction")
lines(pred.ada.net, col = 3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(15, 15, legend=c("Actual", "AEN"),col=c("blue", "green"), lty=1, cex=0.8)

plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Adaptive Elastic Net Prediction")
lines(pred.ada.net, col = 3, type = "b")
lines(mean, col = 2, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(15, 15, legend=c("Actual", "AEN", "Mean"),col=c("blue", "green", "red"), lty=1, cex=0.8)

plot(y.out - mean, type="b", col=4, xlab = "", ylab = "Deviation from Mean", xaxt = "n", ylim = c(-20,20), main = "Adaptive Elastic Net Prediction")
lines(pred.ada.net - mean, col = 3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 0, b = 0, col = 1)
legend(15, -15, legend=c("Actual", "AEN"),col=c("blue", "green"), lty=1, cex=0.8)



Metrics::mse(y.out, pred.ada.net)

fitted.elas.net = predict(ada.net, x.in)
plot(y.in, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(-20,20), main = "In Sample Lasso Estimation")
lines(fitted.elas.net, col = 3, type= "b")
abline(a = 0, b = 0, col = 1)

## In-sample stuff

fitted.cluster=predict(adalasso,newdata = x.in)
fitted.adalasso = predict(adalasso.1, newdata = x.in)
fitted.lasso=predict(lasso,newdata = x.in)
fitted.post = nonZero.in %*% betas
plot(y.in, col=4, type = "b", xlab = "", ylab = "Deviation of Mean", xaxt = "n", ylim = c(-20,20), main = "In sample stuff")
axis(side = 1, at = seq(1,51,3), labels = states, las = 2)
abline(a = 0, b = 0, col = 1)
lines(fitted.post, col= 3, type = "b")
lines(pred.factors, col=3, type = "b")
lines(fitted.adalasso, col = 3, type = "b")
lines(fitted.lasso, col=3, type = "b")
lines(fitted.elas.net, col = 3, type ="b")



Metrics::mse(y.in, fitted.lasso)
Metrics::mse(y.in, fitted.adalasso)
Metrics::mse(y.in, fitted.cluster)
Metrics::mse(y.in, fitted.post)
Metrics::mse(y.in, fitted.elas.net)

mean.in <- rep(0,51)
index <- 1
for (i in seq(1,51,3)){
  mean.in[i] <- mean[index]
  mean.in[i+1] <- mean[index]
  mean.in[i+2] <- mean[index]
  index = index + 1
}


plot(y.in+mean.in, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "In Sample Lasso Estimation")
lines(fitted.lasso + mean.in, col=3, type = "b")
lines(fitted.adalasso + mean.in, col=3, type = "b")
lines(fitted.post+mean.in, col = 3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)

########### Relevant Regressors ###########
# Lasso Relevant
count <- 0
coef.dummy <- lasso$coefficients
for (i in 2:length(coef.dummy)){
  if (coef.dummy[i] != 0) {
    count = count + 1
    #print(basic.searches[ceiling((i-1)/5)])
    print(coef.dummy[[i]])
    #print((i-1) %% 5)
  }
}
print(count)
coef.dummy <- NULL
# Ada Lasso Relevant
count <- 0
coef.dummy <- adalasso.1$coefficients
for (i in 2:length(coef.dummy)){
  if (coef.dummy[i] != 0) {
    count = count + 1
    #print(basic.searches[ceiling((i-1)/5)])
    print(coef.dummy[[i]])
    #print((i-1) %% 5)
  }
}
print(count)
coef.dummy <- NULL
# Cluster Lasso Relevant
count <- 0
coef.dummy <- adalasso$coefficients
for (i in 2:length(coef.dummy)){
  if (coef.dummy[i] != 0) {
    count = count + 1
    #print(basic.searches[ceiling((i-1)/5)])
    print(coef.dummy[[i]])
    #print((i-1) %% 5)
  }
}
print(count)
coef.dummy <- NULL
# Post Adaptive Lasso Relevant
coef.dummy <- betas
for (i in 1:length(coef.dummy)) {
  print(coef.dummy[[i]])
}
coef.dummy <- null
# Adaptive Elastic Net Relevant
count <- 0
coef.dummy <- ada.net$beta
for (i in 2:length(coef.dummy)){
  if (coef.dummy[i] != 0) {
    count = count + 1
    print(basic.searches[ceiling((i-1)/5)])
    #print(coef.dummy[i])
    #print((i) %% 5)
  }
}
print(count)
coef.dummy <- NULL



## Bai and Ng Information Criteria
g1 <- function(p, n) {
  return ((n + p)/(n*p) + log((n*p)/(n+p)))
}

g2 <- function(p,n) {
  return (((n+p)/(n*p)) * log(min(p,n)))
}

g3 <- function(p,n) {
  return(log(min(p,n))/min(p,n))
}

V <- function(r, p, n, x, lambdas, factors) {
  sum <- 0
  for (j in 1:p) {
    for (t in 1:n){
      sum <- sum + (x[t,j] - t(lambdas[,j])%*%(factors[t, 1:r]))^2 
    }
  }
  return (sum/(n*p))
}

V_j <- function(j, r, n, x, lambdas, factors) {
  sum <- 0
  for (t in 1:n){
    sum <- sum + (x[t,j] - t(lambdas[,j])%*%(factors[t, 1:r]))^2 
  }
  return (sum)
}

TSS_j <- function(j, x) {
  sum <- 0
  for (t in 1:n){
    sum <- sum + (x[t,j])^2 
  }
  return (sum)
}
PC <- function(type, k, r, p, n, x, lambdas, factors, sigma) {
  first <- V(r, p, n, x, lambdas, factors)
  g <- 0
  if (type == 1) {
    g <- g1(p, n)
  }
  else if (type == 2) {
    g <- g2(p, n)
  }
  else {
    g <- g3(p, n)
  }
  return (first + (k * sigma * g))
}

IC <- function(type, k, r, p, n, x, lambdas, factors) {
  first <- log(V(r, p, n, x, lambdas, factors))
  g <- 0
  if (type == 1){
    g <- g1(p,n)
  }
  else if (type == 2){
    g <- g2(p,n)
  }
  else {
    g <- g3(p,n)
  }
  return (first + k * g)
}
### FACTOR MODELS ###
pca <- dudi.pca(df = x.in, center = TRUE, scale = TRUE, scannf = FALSE, nf = 25)
f <- pca[["li"]]
f <- as.matrix(f)

x_standardized <- matrix(0, nrow = 51, ncol = p)
for(j in 1:p){
  x_standardized[,j] = scale(x.in[,j])
}

eigen_st <- eigen(x_standardized%*%t(x_standardized))

r_adj <- rep(0,11)
for(k in 1:10){
  r_adj[k+1] = (1/(rows*p))*sum(eigen_st$values[1:k])
}

s <- c(0:9)
plot(s,diff(r_adj), type = "l", main = "Screeplot of factors", xlab = "Number factors", ylab = "Change in Adj R-square")
abline(v=4, col= "red")

getlambdas <- function(r, f) {
  return(solve(t(f[,1:r]) %*% f[,1:r])%*% t(f[,1:r])%*%x.in)
}
sigma <- V(25, p, 51, x.in, getlambdas(25,f), f)

criteria <- matrix(NA, nrow=25, ncol=6)
for (r in 1:25) {
  v <- V(r,p,51,x.in, getlambdas(r,f),f)
  criteria[r,1] <- v + r*sigma*g1(p,51)
  criteria[r,2] <- v + r*sigma*g2(p,51)
  criteria[r,3] <- v + r*sigma*g3(p,51)
  criteria[r,4] <- log(v) + r*g1(p,51)
  criteria[r,5] <- log(v) + r*g2(p,51)
  criteria[r,6] <- log(v) + r*g3(p,51)
}
pc1 <- which.min(criteria[,1])
pc2 <- which.min(criteria[,2])
pc3 <- which.min(criteria[,3])
ic1 <- which.min(criteria[,4])
ic2 <- which.min(criteria[,5])
ic3 <- which.min(criteria[,6])

numberFactors <- 4

factors <- f[,1:numberFactors]
estimation_fact <- solve(t(factors)%*%factors)%*%(t(factors)%*%y.in)
pred.factors <- factors %*% estimation_fact


factors.out <- as.matrix(predict(pca,newdata= x.out, olddata = x.in))
factors.out <- factors.out[,1:numberFactors]
y.predict <- factors.out %*% estimation_fact

y.fitted <- factors %*% estimation_fact

plot(y.out, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Factor Model Prediction (4 Factors)")
lines(y.predict + mean, col=3, type = "b")
lines(mean, col = 2, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
legend(14, 10, legend=c("Actual", "4 Factor Model", "Mean"),col=c("blue", "green", "red"), lty=1, cex=0.8)



plot(y.out - mean, type="b", col=4, xlab = "", ylab = "Deviation from Mean ", xaxt = "n", ylim = c(-20,20), main = "Factor Model Prediction (4 Factors)")
lines(y.predict, col=3, type = "b")
axis(side = 1, at = seq(1,17,1), labels = states, las = 2)
abline(a = 0, b = 0, col = 1)
legend(14, -15, legend=c("Actual", "4 Factor Model"),col=c("blue", "green"), lty=1, cex=0.8)


plot(y.in + mean.in, type="b", col=4, xlab = "", ylab = "Percentage Democrats", xaxt = "n", ylim = c(0,100), main = "Factor Model Prediction (4 Factors)")
lines(y.fitted + mean.in, col=3, type = "b")
axis(side = 1, at = seq(1,51,3), labels = states, las = 2)
abline(a = 50, b = 0, col = 1)
Metrics::mse(y.out, y.predict + mean)
Metrics::mse(y.in , y.fitted)
Metrics::mse(y.out, pred.lasso)
Metrics::mse(y.out, pred.adalasso)

### Post-Double Selection #####

cube <- array(data = NA, dim = c(180,180,3))
for(i in 1:180) {
  new.x <- x.in[, -i]
  inner.est = ic.glmnet(new.x, y.in,crit="bic", penalty.factor = weights[-i])
  tempor <- inner.est$coefficients
  w.est = ic.glmnet(new.x, x.in[,i],crit="bic", penalty.factor = weights[-i])
  temporary <- w.est$coefficients
  drop(inner.est)
  drop(w.est)
  for(j in 1:180) {
    if(i ==j) {
      cube[i,i,1] <- 0
      cube[i,i,2] <- 0
      cube[i,i,3] <- 0
    }
    else if( i < j) {
      cube[i,j,1] <- tempor[(j)]
      cube[i,j,3] <- temporary[(j)]
    }
    else {
      cube[i,j,1] <- tempor[(j+1)]
      cube[i,j,3] <- temporary[(j+1)]
    }
  }
}

for(i in 1:180) {
  for(j in 1:180) {
    if(i != j) {
      if(cube[i,j,1] != 0 || cube[i,j,3] !=0) {
        cube[i,j,2] <- 1
      }
      else {
        cube[i,j,2] <- 0
      }
    }
  }
}

result <-matrix(NA, nrow=180, ncol=4)

for (i in 1:180) {
  xset <- x.in
  for(j in 1:180) {
    if(cube[i,j,2] ==0) {
      xset[,j] <- 0
    }
  }
  OLSreg <- lm(formula = y.in ~ x.in[,i] + xset -1)  ### lm regression without intercept ##
  result[i,] <- summary(OLSreg)$coefficients[1,]  
}
co <- 0
for (i in 1:length(result[,1])) {
  if (result[i,4] < 0.05) {
    co <- co + 1
    print(result[i,4])
    #print(basic.searches[ceiling((i)/5)])
    #print((i) %% 5)
  }
}
print(co)
##### Bonferroni ####
pval <- result[,4]
indx <- seq(1,5)
name <- vector(mode="character",length = 180)
week <- vector(mode="character", length=180)
for(i in 1:36)
{
  for(j in 1:5)
  {
    name[(5*(i-1)+j)] <- basic.searches[i]
    week[(5*(i-1)+j)] <- indx[j]
  }
}
pval<- data.frame(pval,name,indx)
alpha_bonf <- 0.05/180
sign_bonf <- pval[pval[,1] < alpha_bonf,]

##### Holm's #####
alpha_holms <- vector(mode="numeric",length = 180)
alpha_holms[1] <- alpha_bonf
for(i in 2:180)
{
  alpha_holms[i] <- 0.05/(180+1-i)
}
pval_holms <- pval[order(pval$pval),]
pval_holms <- data.frame(pval_holms, alpha_holms)
sign_holms <- pval_holms[pval_holms$pval < pval_holms$alpha_holms,]

################ End of Post double selection ######
