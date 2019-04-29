library(MASS)
data(Boston)
head(Boston)
summary(Boston)

#--------------------------------model selection with BIC, CP, and adjusted R squared ----------------------------


#################################### Best Subset Selection #################################### 
#install.packages("leaps")
library(leaps)
# Display only the best models for each number of subsets
leaps<-regsubsets(medv~.,data=Boston, nbest=1, nvmax = 13)
S<-summary(leaps)
S
#plot results and best choice
par(mfrow =c(2,2))
plot(S$rss ,xlab=" Number of Variables ",ylab=" RSS",type="l",lwd=2)
w<-which.min (S$rss)
points(w,S$rss[w], col ="red",cex =2, pch =20)

plot(S$adjr2 ,xlab =" Number of Variables ",ylab=" Adjusted RSq",type="l",lwd=2)
w<-which.max (S$adjr2)
points(w,S$adjr2[w], col ="red",cex =2, pch =20)

plot(S$cp ,xlab =" Number of Variables ",ylab="Cp",type="l",lwd=2)
w<-which.min (S$cp)
points(w,S$cp[w], col ="red",cex =2, pch =20)

plot(S$bic ,xlab =" Number of Variables ",ylab="BIC",type="l",lwd=2)
w<-which.min (S$bic)
points(w,S$bic[w], col ="red",cex =2, pch =20)

# more graphs
plot(leaps,scale="bic")
plot(leaps,scale="Cp")

# Display coefficients: BIC selection selected 10 coefficient
coef(leaps,10)

#################################### Forward Selection ####################################
leaps.fwd<-regsubsets(medv~.,data=Boston,method ="forward", nvmax = 13)
S<-summary(leaps.fwd)
S
#plotting results and highlight best choice
par(mfrow =c(2,2))
plot(S$rss ,xlab=" Number of Variables ",ylab=" RSS",type="l",lwd=2)
w<-which.min(S$rss)
points(w,S$rss[w], col ="red",cex =2, pch =20)

plot(S$adjr2 ,xlab =" Number of Variables ",ylab=" Adjusted RSq",type="l",lwd=2)
w<-which.max(S$adjr2)
points(w,S$adjr2[w], col ="red",cex =2, pch =20)

plot(S$cp ,xlab =" Number of Variables ",ylab="Cp",type="l",lwd=2)
w<-which.min(S$cp)
points(w,S$cp[w], col ="red",cex =2, pch =20)

plot(S$bic ,xlab =" Number of Variables ",ylab="BIC",type="l",lwd=2)
w<-which.min(S$bic)
points(w,S$bic[w], col ="red",cex =2, pch =20)

# more graphs
par(mfrow =c(1,1))
plot(leaps.fwd,scale="bic")
plot(leaps.fwd,scale="Cp")

# BIC selected 10 coefficients
coef(leaps.fwd, 10)

#################################### Backward Selection ####################################
leaps.bwd<-regsubsets(medv~.,data=Boston,method ="backward", nvmax = 13)
S<-summary(leaps.fwd)
S
#plotting results and highlight best choice
par(mfrow =c(2,2))
plot(S$rss ,xlab=" Number of Variables ",ylab=" RSS",type="l",lwd=2)
w<-which.min(S$rss)
points(w,S$rss[w], col ="red",cex =2, pch =20)

plot(S$adjr2 ,xlab =" Number of Variables ",ylab=" Adjusted RSq",type="l",lwd=2)
w<-which.max(S$adjr2)
points(w,S$adjr2[w], col ="red",cex =2, pch =20)

plot(S$cp ,xlab =" Number of Variables ",ylab="Cp",type="l",lwd=2)
w<-which.min(S$cp)
points(w,S$cp[w], col ="red",cex =2, pch =20)

plot(S$bic ,xlab =" Number of Variables ",ylab="BIC",type="l",lwd=2)
w<-which.min(S$bic)
points(w,S$bic[w], col ="red",cex =2, pch =20)

# more graphs
par(mfrow =c(2,2))
plot(leaps.bwd,scale="bic")
plot(leaps.bwd,scale="Cp")

# BIC selected 10 coefficients
coef(leaps.bwd, 10)

#################################### Stepwise Selection ####################################
leaps.stepwise<-regsubsets(medv~.,data=Boston,method ="seqrep", nvmax = 13)
S<-summary(leaps.stepwise)
plot(leaps.stepwise,scale="bic")
plot(leaps.stepwise,scale="Cp")

#################################### Ridge Regression ####################################
#install.packages("genridge")
library(genridge) 

lambda<-seq(0,100,length.out=500)
L<-ridge(medv~., data=Boston,lambda=lambda)

# Visualize how the coefficients shrink with increasing lambda
traceplot(L, X=c("lambda"))

# Tune the lambda according to GCV and plot
plot(lambda,L$GCV,pch=16,type="b",ylab="GCV",xlab=expression(lambda))
lambda <- L$kGCV
lambda

# Find minimum
abline(v=L$kGCV)

# Display coefficients for optimal lambda
M <- as.matrix(L$GCV)
L$coef[which.min(M),]

# Final result
traceplot(L,X=c("lambda"))
abline(v=L$kGCV,col="black",lwd=3)


#################################### Lasso Regression ####################################
#install.packages("lars")
library(lars)
X <- scale(Boston[,-14])
# computes for entire regularization path
y <- Boston[,14]
m <- lars(X,y,type="lasso")
# plot of ||beta_lasso||_1 / ||\beta_OLS||_1
mlambda<-cv.lars(X,y,K=10,trace=FALSE,type="lasso")
m
plot(m)
# What are the coefficients
coef(m)[which.min(m$Cp),]

#################################### PCR ####################################
#install.packages("pls")
library(pls)
library(ISLR)
set.seed (1)

# Fit PCR model
pcr.fitfull<-pcr(medv~., data=Boston ,scale =TRUE,validation ="LOO")
summary(pcr.fitfull)

# Show Cv plot
plot(RMSEP(pcr.fitfull), legendpos = "topright")

# Fit complete model with 5 components (5 components based on plot above)
pcr.fit<-pcr(medv~., data=Boston, scale =TRUE,ncomp =5)
summary(pcr.fit)

# Plot % explained variances
plot(explvar(pcr.fitfull),type="b",xlab="No. of Components")
explvar(pcr.fitfull)

# Final regression coef.
coef(pcr.fit,intercept=TRUE)

#--------------------------------select model using validation: select model witht the lowest MSE----------------------------
############ best subset, forward, backward selection ############

# split data into training (70%) and testing (30%)
trainIDs = sample(x=1:nrow(Boston), size=nrow(Boston)*0.7)
train = Boston[trainIDs, ]
test = Boston[-trainIDs, ]
test.matrix <- model.matrix(medv~ ., data = test)
train.matrix <- model.matrix(medv~ ., data = train)

# function that finds the model with the lowest MSE (best subset, forward and backward selection)
validation_selection = function(model){
  val.errors <- rep(NA, 13)
  for (i in 1:13) {
    coefi <- coef(model, id = i)
    pred <- test.matrix[, names(coefi)] %*% coefi
    val.errors[i] <- mean((test$medv - pred)^2)
  }
  plot(val.errors, xlab = "Number of predictors", ylab = "Test MSE", pch = 13, type = "b",col="blue")
  number_of_predictors = which.min(val.errors)
  print(paste0("Number of predictors that yield the lowest MSE is ", number_of_predictors))
  coeficients = coef(model, number_of_predictors)
  print("Coeficients: ")
  print(coeficients)
  MSE = min(val.errors)
  print(sprintf("Test MSE is %10.2f", MSE))
}

# best subset selection
best_subset = regsubsets(medv~ ., data = train, nvmax=13)
validation_selection(best_subset)

# forward selection
forward_selection = regsubsets(medv~ ., data = train, nvmax=13, method = "forward")
validation_selection(forward_selection)

# backward selection
backward_selection = regsubsets(medv~ ., data = train, nvmax=13, method = "backward")
validation_selection(backward_selection)

# predict function from chapter 6 from the book
predict.regsubsets <- function(object, newdata, id, ...){
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form, newdata)
  coefi <- coef(object, id=id)
  xvars <- names(coefi)
  mat[,xvars]%*%coefi
}
############ ridge and lasso regression ############
# prepare data
X =model.matrix(medv~., Boston)[, -1] # independent variables
y = Boston$medv # dependent variable
X_train = X[trainIDs, ]
X_test = X[-trainIDs, ]
y_train = y[trainIDs]
y_test = y[-trainIDs]
#install.packages("glmnet")
library(glmnet)

regression = function(alpha){
  regression = glmnet(X_train, y_train, alpha = alpha)
  predict_regression = predict(regression, s = lambda, newx = X_test)
  MSE = mean((predict_regression - y_test)^2)
  print(sprintf("Test MSE is %10.2f", MSE))
  coeficients = predict(regression, s=lambda, type="coefficients")
  print(coeficients)
}
print("Ridge")
ridge = regression(alpha = 0)

print("Lasso")
lasso = regression(alpha =1)

############### OLS and PCR ##################################
# OLS
lm = lm(medv~., data = train)
lm.pred = predict(lm, data.frame(X_test))
OLS.MSE = mean((lm.pred - y_test)^2)
print(sprintf("Test MSE is %10.2f", OLS.MSE))

# PCR
pcr = pcr(medv ~ ., data = train, scale = TRUE, validation = "CV")
validationplot(pcr, val.type = "MSEP")
pcr.predict = predict(pcr, test, ncomp = 3) # ncomp = 3 based on the plot
PCR.MSE = mean((pcr.predict - y_test)^2)
print(sprintf("Test MSE is %10.2f", PCR.MSE))

