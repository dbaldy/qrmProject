## Preamble
library(readr)
library(stats4)
library(VineCopula)
sink(file="/Users/damienbaldy/Documents/Saint-Gall/Quantitative\ Risk\ Management/output.txt")

# Load the data from csv
filepath <- "/Users/damienbaldy/Documents/Saint-Gall/Quantitative\ Risk\ Management/qrm17HSG_assignmentdata.csv"
QRM_DATA <- read_delim(filepath, ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# Clear the empty cells at the end
table_data <- QRM_DATA[1:3197,]

# Parse data and calculate returns, column 1 is X and column 2 is Y
counter <- 1
for (column in table_data[2:ncol(table_data)]) {
  return = c()
  previous <- column[1]
  for (val in column[2:length(column)]) {
    return <- append(return, (val - previous)/previous)
    previous <- val
  }
  if (counter == 1) X=return
  else Y=return
  counter <- counter + 1
}

# Initialize parameters
lambda1 <- 0.3
lambda2 <- 0.7
initial_wealth <- 10000000

# Model 2 MLE function
estimateM2 <- function(Xl, Yl) {
  # We calculate the seeds for the MLE 
  mu1 <- mean(Xl)
  mu2 <- mean(Yl)
  sig1 <- sd(Xl)
  sig2 <- sd(Yl)
  rho <- cor(tail(X, N), tail(Y, N))
  
  # MLE of bivariate Gaussians (Model M2)
  bivariateDensGaussians = function(mu1, mu2, sig1, sig2, rho) -sum(log(1/(2*pi*sig2*sig1*sqrt(1-rho))* 
                                                                          exp(-1/(2*(1-rho^2))*((Xl - mu1)^2/sig1^2+(Yl-mu2)^2/sig2^2-2*rho*(Xl - mu1)*(Yl - mu2)/(sig1*sig2)))))
  estBivariateGaussians <- mle(bivariateDensGaussians, start = list(mu1=mu1, mu2=mu2, sig1=sig1, sig2=sig2, rho=0.5), method = "Nelder-Mead")  
  return(estBivariateGaussians)
}

# Model 4 MLE function
estModel4 <- function(Xl, Yl, mu1, mu2, sig1, sig2) {
  # First we calculate the degrees of freedom of the error terms individually
  Xerror <- (Xl - mu1)/sig1
  Yerror <- (Yl - mu2)/sig2
  XunivariateT = function(df) {
    -sum(log(gamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2))*(1 + Xerror^2/df)^(-(df+1)/2)))
  }
  estDf1 <- mle(XunivariateT, start = list(df=20))
  
  YunivariateT = function(df) {
    -sum(log(gamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2))*(1 + Yerror^2/df)^(-(df+1)/2)))
  }
  estDf2 <- mle(YunivariateT, start = list(df=20))
  
  df1 <- estDf1@coef[1]
  df2 <- estDf2@coef[1]
  
  # We estimate Gaussian Copula parameter
  XerrorT <- pt(Xerror, df1)
  YerrorT <- pt(Yerror, df2)
  estError <- BiCopEst(XerrorT, YerrorT, 1, method="mle")
  # print(paste(c("rho:", estError$par), collapse = " "))
  rhoStudents <- estError$par
  
  # Then we calculate the degrees of freedom of the bivariate t
  bivariateT = function(df) -sum(log(gamma((df+2)/2)*
                                            (1 + (Xerror^2 - 2*rhoStudents*Xerror*Yerror + Yerror^2)/
                                               (df*(1 - rhoStudents^2)))^(-(df+2)/2)
                                          /(gamma(df/2)*pi*df*sqrt(1 - rhoStudents^2))))
  estM4 <- mle(bivariateT, start = list(df=20))
  df <- estM4@coef[1]
  return(c(df, rhoStudents, df1, df2))
}

# General bivariate gaussion simulation 
bivariateGaussianSim <- function(mu1, mu2, sig1, sig2, rhoGaussians) {
  # We define sigma and mu
  sigma <- matrix(c(sig1^2, rhoGaussians * sig1 * sig2, rhoGaussians * sig1 * sig2, sig2^2), nrow=2, ncol=2)
  mu <- matrix(c(mu1, mu2), ncol=1)
  
  # Choleski decomposition
  A <- chol(sigma)
  
  Y <- c()
  for (i in 1:10000) {
    # Generate normally distributed results for Z1 and Z2, normally distributed standard random variables 
    Z <- matrix(rnorm(2), ncol=1)
    
    # Generate Y
    Y <- append(Y, t(mu + t(A) %*% Z))
  }
  # Split Y in a matrix with 2 columns
  Y <- matrix(Y, ncol = 2)
  return(Y)
}

modelM4Sim <- function(rhoStudents, mu, sig1, sig2, df) {
  Z <- bivariateGaussianSim(0, 0, 1, 1, rhoStudents)
  R <- c()
  for (i in 1:10000) {
    # Generate the bivariate t
    R <- append(R, t(mu) + matrix(c(sig1, sig2), ncol=2) * sqrt(df/rchisq(1, df))* Z[i,])
  }
  R <- matrix(R, ncol = 2)
  return (R)
}

# Main function for questions (i) to (iv)
runModels <- function(N)  {
  print(paste(c("We run the model with N=", N, collapse=" ")))
  Xl <- tail(X, N)
  Yl <- tail(Y, N)
  
  estBivariateGaussians <- estimateM2(Xl, Yl)
  mu1 <- estBivariateGaussians@coef[1]
  mu2 <- estBivariateGaussians@coef[2]
  sig1 <- estBivariateGaussians@coef[3]
  sig2 <- estBivariateGaussians@coef[4]
  rhoGaussians <- estBivariateGaussians@coef[5]
  mu <- matrix(c(mu1, mu2), ncol=1)
  print("MLE Model M2")
  print(estBivariateGaussians@coef)
  
  # MLE of Gumbel Copula (VineCopula package, Model M3)
  # We already estimated mu1, mu2, sig1, and sig2, 
  # we need to estimate the parameter of the Gumbel Copula
  if (N == 500) {
    Xlnorm <- pnorm((Xl - mu1)/sig1)
    Ylnorm <- pnorm((Yl - mu2)/sig2)
    estGumb <- BiCopEst(Xlnorm, Ylnorm, 4, method="mle")
    theta <- data.frame(c(estGumb$par))
    colnames(theta) <- c("M3 theta")
    print("Gumbel Copula Model M3")
    print(paste(c("theta:", estGumb$par, collapse=" ")))  
  }
  
  # Model 4 estimation
  model4 <- estModel4(Xl, Yl, mu1, mu2, sig1, sig2)
  df <- unname(model4[1])
  rhoStudents <- unname(model4[2])
  df1 <- unname(model4[3])
  df2 <- unname(model4[4])
  model4Output <- data.frame(df1=df1, df2=df2, df=df, rhoStudents=rhoStudents)
  print("MLE model 4")
  print(model4Output)
  
  # (ii) We define the function to estimate the VaR and the Espected shortfall
  VaREspSh <- function(portfolioReturns, modelNumber) {
    # We calculate the quantiles
    quant90 <- quantile(portfolioReturns, 0.1)
    quant95 <- quantile(portfolioReturns, 0.05)
    quant99 <- quantile(portfolioReturns, 0.01)
    VaRs <- c(-quant90, -quant95, -quant99) * initial_wealth
    
    print(paste(c("Risk measures Model", modelNumber), collapse = " "))
    print("Value at risk")
    print(VaRs)
    espSh <- data.frame(`90%`=-mean(portfolioReturns[portfolioReturns <= quant90]),
                        `95%`=-mean(portfolioReturns[portfolioReturns <= quant95]),
                        `99%`=-mean(portfolioReturns[portfolioReturns <= quant99]),
                        check.names=FALSE) * initial_wealth
    print("Espected shortfalls")
    print(espSh)
  }
  
  # We simulate model M1 portfolio returns
  simulationX <- sample(Xl, 10000, replace = TRUE)
  simulationY <- sample(Yl, 10000, replace = TRUE)
  portfolioReturns <- lambda1*simulationX + lambda2*simulationY
  VaREspSh(portfolioReturns, "M1")
  
  
  # (iii) We simulate model M2
  Y <- bivariateGaussianSim(mu1, mu2, sig1, sig2, rhoGaussians)
  portfolioGaussians <- lambda1 * Y[,1] + lambda2 * Y[,2]
  VaREspSh(portfolioGaussians, "M2")
  
  # We simulate model M4
  R <- modelM4Sim(rhoStudents, mu, sig1, sig2, df)
  portfolioWithTErrors <- lambda1 * R[,1] + lambda2 * R[,2]
  VaREspSh(portfolioWithTErrors, "M4")
}

nbObs <- c(500, 100, 200, 1000)
for(N in nbObs) {
  runModels(N)
}

# (v) Rolling window
getQuantiles <- function(X, Y) {
  # Model M2
  estBivariateGaussians <- estimateM2(X, Y)
  mu1 <- estBivariateGaussians@coef[1]
  mu2 <- estBivariateGaussians@coef[2]
  sig1 <- estBivariateGaussians@coef[3]
  sig2 <- estBivariateGaussians@coef[4]
  rhoGaussians <- estBivariateGaussians@coef[5]
  mu <- matrix(c(mu1, mu2), ncol=1)
  
  # Model M4
  model4 <- estModel4(X, Y, mu1, mu2, sig1, sig2)
  df <- model4[1]
  rhoStudents <- model4[2]
  
  # Simulate the two models
  Y <- bivariateGaussianSim(mu1, mu2, sig1, sig2, rhoGaussians)
  portfolioGaussians <- lambda1 * Y[,1] + lambda2 * Y[,2]
  
  R <- modelM4Sim(rhoStudents, mu, sig1, sig2, df)
  portfolioWithTErrors <- lambda1 * R[,1] + lambda2 * R[,2]
  return(c(quantile(portfolioGaussians, 0.05), quantile(portfolioWithTErrors, 0.05)))
}

nbOfEstimations <- length(X)-201
countViolationsM2 <- 0
countViolationsM4 <- 0
for(i in 1:nbOfEstimations) {
  last <- i + 200
  quant <- getQuantiles(X[i:last], Y[i:last])
  nextValue <- lambda1*X[last+1] + lambda2*Y[last+1]
  if (nextValue < quant[1]) countViolationsM2 <- countViolationsM2 + 1
  if (nextValue < quant[2]) countViolationsM4 <- countViolationsM4 + 1
}

perViolationsM2 <- countViolationsM2 / nbOfEstimations
perViolationsM4 <- countViolationsM4 / nbOfEstimations
violations <- data.frame(`% violations M2`=perViolationsM2, `% violations M4`=perViolationsM4, check.names = FALSE)
print(violations)