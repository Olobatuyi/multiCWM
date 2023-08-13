# # load useful libraries
# library(knitr)
# library(gridExtra);library(grid);require(keras);library(tidyverse);
# library(caret);library(ggplot2);library(SIBER)
# library(factoextra);library(FactoMineR);library(xtable)
# library(MBCbook);library(HDclassif);library(dplyr);library(SIBER)
# require("markdown");require("rattle");require("xtable");require("stringr")
# require("fBasics");require("MASS");require("survival");#require("STAR")
# require("gamlss.dist");require("VGAM");library(rgl);library(ellipse)
# library(threejs);library(plotROC);library(ROCR);library(pROC)
# library(rgl)
#
# library(tidyverse);library(Rlab)
# library(coda);library(VGAM)
# library(mvnfast); library(fastDummies)
# library(Rtsne)
# require("gamlss.dist")
# library(optimization)
# require("ggplot2")


#' Logit Function for depencies

#' The function performs a logit function with centered around the maximum to ignore underflow errors
#' @param x numeric values to be transformed between 0 and 1
#'
#' @export
#' @keywords logit
#'@rdname logw
#'
logw <- function(x) return(exp(x - max(x)) / (1 + sum(exp(x - max(x)))))


#' Logit Function for depencies

#' The function performs a logit function with centered around the maximum to ignore underflow errors
#' @param x numeric values to be transformed between 0 and 1
#'
#' @export
#' @keywords logit
#' @rdname logq
#'
logq <- function(x) return(1 / (1 + sum(exp(x - max(x)))))

#' Confusion Matrix function to plot the Confusion Matrix

#' The function plots the confusion matrix of the classification results
#' @param x The predicted classes returned from the clustering algorithm
#'
#' @export
#' @keywords confusion matrix
#' @rdname conf
#'

conf <- function(x){

  cm_d  <- as.data.frame(x$table)
  cm_st <- data.frame(x$overall)
  cm_st$x.overall <- round(cm_st$x.overall, 2)
  cm_p <- as.data.frame(prop.table(x$table))
  cm_d$perc <- round(cm_p$Freq*100,2)
  cm_d_p <- ggplot(data = cm_d, aes(x = Prediction, y = Reference, fill = Freq)) +
    geom_tile()+
    geom_text(aes(label = paste("", Freq, ",",perc, "%")), col = "red")+
    theme_light()+
    guides(fill = FALSE)
  cm_st_p <- tableGrob(cm_st)

  print(grid.arrange(cm_d_p, cm_st_p, nrow = 1, ncol = 2,
                     top = textGrob("Confusion Matrix and Statistics", gp = gpar(fontsize = 25, font = 1))))
}



################################################################################

#' Simulation Function for fake data

#' Simulated function simulates from the distributions assumed by the model
#' @param n The number of values to simulate. n must be an integer
#' @param G The number of components which must be an integer
#' @param d The dimension number of features x: must be an integer.
#' @param J The dimension number of categorical variable Y: must be an integer
#' @param mean The true parameter mean for each component on each row of the matrix
#' @param sigma The true values of the variance-covariance matrix
#' @param w The true values of coefficients for the multinomial regression Y ~ X
#' @param Pi The mixing probability of the components
#' @param M The repeat trial
#' @export
#' @keywords simulated data
#' @rdname simData
#'


simData <- function(n, G, d = 2, M = 1, J, mean, sigma, w, Pi = NULL){

  if(n %% 1 != 0) stop("STOP!!! n must be an integer")
  if(G %% 1 != 0) stop("STOP!!! G must be an integer")
  if(d %% 1 != 0) stop("STOP!!! d must be an integer")
  if(J %% 1 != 0) stop("STOP!!! J must be an integer")

  c <- J-1;
  a <- pr <- matrix(NA, nrow = n, ncol = c)
  x <- matrix(NA, nrow = n, ncol = d)
  Y <- py <- matrix(NA, nrow = n, ncol = J)

  # P1 parameter
  if(is.null(Pi)) Pi <- rep(1/G,G)

  # Set a seed
  U <- runif(n); z <- NULL

  if(G == 2){

    for (i in 1:n) {

      if(U[i] <= Pi[1]){

        z[i] <- 1
        x[i,] <- mvnfast::rmvn(1, mean[1,], sigma[[1]])
        a[i,]    <- as.matrix(x[i,] %*% t(w[[1]]))
        pr[i,]   <- exp(a[i,]) / (1 + sum(exp(a[i,])))
        pu       <- 1 / (1 + sum(exp(a[i,])))
        py[i,]   <- c(pu, pr[i,])
        Y[i,] <- t(rmultinom(1, M, py[i,]))

      }else{

        z[i]   <- 2
        x[i,]  <- mvnfast::rmvn(1, mean[2,], sigma[[2]])
        a[i,]  <- as.matrix(x[i,] %*% t(w[[2]]))
        pr[i,] <- exp(a[i,]) / (1 + sum(exp(a[i,])))
        pu     <- 1 / (1 + sum(exp(a[i,])))
        py[i,] <- c(pu, pr[i,])
        Y[i,]  <- t(rmultinom(1, M, py[i,]))

      }
    }

  }else{


    for (i in 1:n) {

      if(U[i] <= Pi[1]){

        z[i]   <- 1
        x[i,]  <- mvnfast::rmvn(1, mean[1,], sigma[[1]])
        a[i,]  <- as.matrix(x[i,] %*% t(w[[1]]))
        pr[i,] <- exp(a[i,]) / (1 + sum(exp(a[i,])))
        pu     <- 1 / (1 + sum(exp(a[i,])))
        py[i,] <- c(pu, pr[i,])
        Y[i,]  <- t(rmultinom(1, M, py[i,]))

      }else

        if(U[i] > Pi[1] & U[i] < (Pi[2] + Pi[1])){

          z[i]   <- 2
          x[i,]  <- mvnfast::rmvn(1, mean[2,], sigma[[2]])
          a[i,]  <- as.matrix(x[i,] %*% t(w[[2]]))
          pr[i,] <- exp(a[i,]) / (1 + sum(exp(a[i,])))
          pu     <- 1 / (1 + sum(exp(a[i,])))
          py[i,] <- c(pu, pr[i,])
          Y[i,]  <- t(rmultinom(1, M, py[i,]))

        }else{

          z[i]   <- 3
          x[i,]  <- mvnfast::rmvn(1, mean[3,], sigma[[3]])
          a[i,]  <- as.matrix(x[i,] %*% t(w[[3]]))
          pr[i,] <- exp(a[i,]) / (1 + sum(exp(a[i,])))
          pu     <- 1 / (1 + sum(exp(a[i,])))
          py[i,] <- c(pu, pr[i,])
          Y[i,]  <- t(rmultinom(1, M, py[i,]))

        }
    }

  }

  dat <- list();
  dat$Y <- Y;
  dat$x <- x;
  dat$z <- z
  dat$mean <- mean;
  dat$sigma <- sigma;
  dat$beta <- w

  return(dat)

}




#####################################

#' MultiCwm Function for the clustering of data
#' @param Y   The categorical variables which can be any of the form: vector, one column matrix,
#' or a one-encoding data frame or matrix
#' @param x    The feature or independent variable of the data.
#' @param init This is the initial value to start the clustering algorithm.
#' This can any of the options: random, mclust, or kmeans
#' @param G The number of components which must be an integer
#' @param maxit The maximum number of iterations to terminate the algorithm. This must be an integer.
#' @param tol This is the tolenrance level of convergence dafault to 1e-10.
#' @param show_table This is a logical parameter which is default to FALSE.
#' @export
#' @keywords MultiCwm Clustering Algorithm
#' @rdname MultiCwm
#' @return res which is a list containing the component mean and covariance matrix, the information criterion
#'


MultiCwm <- function(Y, x, G = 1, init = c("random", "mclust", "kmeans"), maxit = 1000, tol = 1e-10, show_table = FALSE){

  if(maxit %% 1 != 0) stop("STOP: maxit must be an integer")
  if(is.data.frame(x)){

    x <- unname(x); x <- as.matrix(x)

  }

  if(is.vector(Y) | ncol(Y) == 1 | is.data.frame(Y)) {

    Y <- unname(Y); Y <- as.matrix(Y)
    Y <- keras::to_categorical(Y)

  }

  eps = sqrt(.Machine$double.eps)
  d  <- ncol(x); J <- ncol(Y); c <- J - 1
  L  <- L1 <- L2 <- ai <- ai3 <- lb <- NULL

  N  <- nrow(x)
  a  <- pr <- matrix(NA, nrow = N, ncol = c)
  py <- matrix(NA, nrow = N, ncol = J)
  dm <- matrix(NA, nrow = N, ncol = G)
  z  <- matrix(0, nrow = N, ncol = G)

  pu  <- NULL
  sig <- array(NA, c(d,d,N))
  #sig <- array(NA, c(N, d))
  et1 <- list()

  w <- w1 <- sigma <- z_val <- p_val <- list()
  m <- (G * c * (1 + d)) + (G * ((d * d - d)/2 + 2 * d + 1))

  if(init == "random"){

    Pi   <- rep(1/G, G)
    mean <- matrix(rnorm(d*G), nrow = G)

    for(l in 1:G){

      sigma[[l]] <- diag(d)
      weights    <- matrix(runif(d*c), nrow = c, ncol = d)
      w[[l]]     <- weights

    }

  }else if(init == "mclust"){

    require(mclust)
    ini_val <- mclust::Mclust(x, G = G, verbose = FALSE)
    Pi      <- ini_val$parameters$pro
    mean    <- t(ini_val$parameters$mean)

    for(l in 1:G){

      sigma[[l]] <- ini_val$parameters$variance$sigma[,,l]
      weights    <- matrix(runif(d*c), nrow = c, ncol = d)
      w[[l]]     <- weights

    }

  }else{

    ini_val <- kmeans(x, centers = G)
    Pi      <- ini_val$size/sum(ini_val$size)
    mean    <- ini_val$centers

    for(l in 1:G){

      sigma[[l]] <- diag(d)
      weights    <- matrix(runif(d*c), nrow = c, ncol = d)
      w[[l]]     <- weights

    }

  }


  count <- 2;

  L <- c(-16000,-15000,-14000);

  ##############################################################################

  # For G = 1

  if(G == 1){

    for(q in 1:2){

      for (i in 1:N) {

        a[i,]   <- x[i,] %*% t(w[[1]])
        pr[i,]  <- logw(a[i,])    #exp(a[i,]) / (1 + sum(exp(a[i,])))
        pu[i]   <- logq(a[i,]); py[i,]  <- c(pu[i], pr[i,])    #1 / (1 + sum(exp(a[i,])));
        dm[i,1] <- dmultinom(Y[i,], sum(Y[i,]), prob = py[i,]) *
          tryCatch(mvnfast::dmvn(x[i,], mean[1,], chol(sigma[[1]]), isChol = TRUE),
                   error = function(e) 1)

      }

      ############ M-Step ##############

      # Posterior distribution
      po <- dm / rowSums(dm)

      # Mean for the group
      mean[1,] <- colSums(x) / N

      L2      <- sum(tryCatch(mvnfast::dmvn(x, mean[1,], chol(sigma[[1]]),
                                            isChol = TRUE, log = TRUE), error = function(e) 1))

      et      <- nnet::multinom(Y ~ x, trace = FALSE)
      w1[[1]] <- summary(et)$coefficients[,-1]

      # Check for NA
      if(!any(is.na(w1[[1]]))) {

        w[[1]]   <- w1[[1]]
        lb    <- as.numeric(logLik(et))
        et1[[1]] <- summary(et)

      }

      lb <- as.numeric(logLik(et))

      for (i in 1:N) sig[,,i] <- outer((x[i,] - mean[1,]),(x[i,] - mean[1,]))
      sigma[[1]]              <- apply(sig, c(1,2), sum)/ N

      # Compute Log_likelihood
      lik <- log1p(1); lik2 <- L2;lik3 <- lb

      # Compute Log_likelihood lki, lk2,
      L[count]     <- lik + lik2 + lik3
      a_k          <- (L[count+1] - L[count]) / (L[count] - L[count-1])
      L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
      dif <- abs(L[count+2] - L[count+1])

      if (show_table) {

        #dif <- abs(L[count+2] - L[count+1])
        out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
        print(knitr::kable(out_table))

      }

    }

  }

  #############################################################################

  # For G > 1

  if(G > 1){

    repeat{

      for(l in 1:G){

        for (i in 1:N) {

          a[i,]   <- x[i,] %*% t(w[[l]])
          pr[i,]  <- logw(a[i,])    #exp(a[i,]) / (1 + sum(exp(a[i,])))
          pu[i]   <- logq(a[i,]); py[i,] <- c(pu[i], pr[i,])    #1 / (1 + sum(exp(a[i,])));
          dm[i,l] <- dmultinom(Y[i,], sum(Y[i,]), prob = py[i,]) *
            tryCatch(mvnfast::dmvn(x[i,], mean[l,], chol(sigma[[l]]), isChol = TRUE),
                     error = function(e) 1) * Pi[l]

        }
      }


      ############ M-Step ##############

      # Posterior distribution

      po <- dm / rowSums(dm)

      # Update Prior
      Pi <- colSums(po) / N

      for(l in 1:G){

        # Mean for the group
        mean[l, ] <- colSums(po[,l] * x) / colSums(po)[l]

        L2[l]  <- sum(po[,l] * tryCatch(mvnfast::dmvn(x, mean[l,], chol(sigma[[l]]),
                                                      isChol = TRUE, log = TRUE), error = function(e) 1))

        et <- nnet::multinom(Y ~ x, weights = po[,l], trace = FALSE)
        w1[[l]] <- summary(et)$coefficients[,-1]

        # Check for NA
        if(!any(is.na(w1[[l]]))) {

          w[[l]]   <- w1[[l]]
          lb[l]    <- as.numeric(logLik(et))
          et1[[l]] <- summary(et)

        }

        for (i in 1:N) sig[,,i] <- po[i,l] * outer((x[i,] - mean[l,]),(x[i,] - mean[l,]))
        sigma[[l]] <- apply(sig, c(1,2), sum)/ sum(po[,l])

      }

      # Compute Log_likelihood
      lik  <- sum(log1p(Pi) * po);lik2 <- sum(L2);lik3 <- sum(lb)

      # Compute Log_likelihood lki, lk2,
      L[count] <- lik + lik2 + lik3
      a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
      L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
      dif <- abs(L[count+2] - L[count+1])

      if (show_table) {

        #dif <- abs(L[count+2] - L[count+1])
        out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
        if(count %% 10 == 0) print(knitr::kable(out_table))

        #if (dif < tol || count == maxit) break;

      }

      if (dif < tol || count == maxit) break;
      count <- count + 1

    }

  }

  Z <- apply(po, 1, which.max)
  Prof <- Poc1 <- po

  for(i in 1:N){
    for(j in 1:G){

      Prof[i,j] <- ifelse(po[i,j] > 0.9, 1, 0)
      if(po[i,j] > 0) Poc1[i,j] <- log1p(po[i,j])

    }

  }
                
  # Calculating z values and p values
  for (l in 1:G) {
    
    # Check the Z-score for the model (wald Z)
    z_val[[l]] <- et1[[l]]$coefficients/et1[[l]]$standard.errors
    p_val[[l]] <- (1 - pnorm(abs(z_val[[l]]), 0, 1)) * 2
    
  }

  ### Information criterion

  ai   <- -2*L[count + 2] + 2*m;
  ai3  <- -2*L[count + 2] + m*log(N)
  AIcc <-  ai - 2*m*(m+1)/(N-m-1)
  AIC3 <- -2*L[count+2] - 3*m
  AICu <-  AIcc - N*log(N/(N-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count+2] - m*(1+log(N))
  AWE  <- -2*L[count+2] - 2*m*(3/2 + log(N))

  res <- list(

    "mean"           = mean,
    "Prob"           = Pi,
    "post"           = po,
    "weights"        = w,
    "classification" = Z,
    "logLik"         = L,
    "AIC"            = ai,
    "BIC"            = ai3,
    "sigma"          = sigma,
    "reg_est"        = et1,
    "z_val"          = z_val,
    "p_val"          = p_val,
    "ICL"            = ICL,
    "AICc"           = AIcc,
    "AIC3"           = AIC3,
    "AICu"           = AICu,
    "Caic"           = Caic,
    "AWE"            = AWE)


  return(res)

}

###############################################

#' zipCwm Function for the clustering of data
#' @param Y The count variables which can be any of the form: vector, one column matrix,
#' @param v The feature or independent variable of the data.
#' @param u The feature or independent variable of the data.
#' @param w The feature or independent variable of the data.
#' @param G The number of components which must be an integer
#' @param maxit The maximum number of iterations to terminate the algorithm. This must be an integer.
#' @param tol This is the tolenrance level of convergence dafault to 1e-10.
#' @param show_table This is a logical parameter which is default to FALSE.
#' @export
#' @keywords zipCwm Clustering Algorithm
#' @rdname zipCwm
#' @return res which is a list containing the component mean and covariance matrix, the information criterion
#'
zipCwm <- function(Y, w, u, v, G, maxit = 1000, tol = 0.1, show_table = FALSE){
  
  x <- cbind(rep(1, nrow(w)), w,u,v);
  if(!is.data.frame(x)) x <- unname(as.matrix(x))

  k <- G
  eps <- sqrt(.Machine$double.eps); 
  n   <- nrow(x)
  W1  <- keras::to_categorical(u); 
  W2  <- keras::to_categorical(v)
  d   <- ncol(x); 
  d1  <- ncol(w); 
  c1  <- ncol(W1); 
  c2  <- ncol(W2)
  L1  <- ai  <- ai3 <- lb <- L <- NULL
  L2  <- lf1 <- lf2 <- lp <- k1 <- matrix(NA, nrow = n, ncol = k)
  dm1 <- lki <- lk2 <- Linf <- al <- NULL
  pr  <- matrix(NA, nrow = n, ncol = (k-1))
  lp1 <- dm <- lam <- matrix(NA, nrow = n, ncol = (k-1))
  b   <- matrix(runif(d*(k-1)), nrow = (k-1))
  
  sigma <- wew <- list()
  mean  <- matrix(runif(d1*(k-1)), nrow = k-1)
  P1i   <- rep(1/k, k)
  P1    <- P1i[1]; Pi <- P1i[-1]
  
  for(l in 1:(k-1))sigma[[l]] <- diag(d1)
  m <- (k-1)*ncol(b) + (2 * (k-1) * d) + (k-1) + (k-1)*(d*d - d)/2
  
  sigma1 <- diag(d1)
  mean1  <-  rnorm(d1)
  
  py1 <- colMeans(W1)
  py2 <- colMeans(W2)
  
  py11 <- matrix(colMeans(W1), nrow = (k-1), ncol = c1, byrow = TRUE)
  py21 <- matrix(colMeans(W2), nrow = (k-1), ncol = c2, byrow = TRUE)
  
  sig1 <- sig <- array(NA, c(d1,d1,n))
  et1  <- list()
  
  count <- 2;
  
  L <- c(-16000,-15000,-14000);
  
  repeat{
    
    a   <- x %*% t(b)
    lam <- exp(a)
    
    for (i in 1:n) {
      
      # Zero count
      dm1[i] <- VGAM::dzipois(Y[i], lambda = 0) * P1 *
        mvnfast::dmvn(w[i,], mean1, sigma1, ncores = 8) *
        dmultinom(W1[i,], 1, prob = py1, log = FALSE) *
        dmultinom(W2[i,], 1, prob = py2, log = FALSE)
      
      for(l in 1:(k-1)){
        
        # Non-zero count
        dm[i,l] <- Pi[l] * mvnfast::dmvn(w[i,], mean[l,], sigma[[l]], ncores = 8) *
          dpois(Y[i], lambda = lam[i,l]) * dmultinom(W1[i,], 1, prob = py11[l,], log = FALSE) *
          dmultinom(W2[i,], 1, prob = py21[l,], log = FALSE) #* y3[i,l+1]
      }
      
    }
    
    po  <- dm / rowSums(cbind(dm1, dm))
    po1 <- 1 - rowSums(po)
    Poc <- cbind(po1, po)
    
    ### For zeroes 
    
    for(i in 1:n){
      
      lf1[i,1]  <- po1[i] * dmultinom(W1[i,], 1, prob = py1, log = FALSE) * dzipois(Y[i], lambda = 0)
      lf2[i,1]  <- po1[i] * dmultinom(W2[i,], 1, prob = py2, log = FALSE) * dzipois(Y[i], lambda = 0)
      
    }
    
    L2[,1]  <- po1 * mvnfast::dmvn(w, mean1, sigma1, ncores = 8,
                     isChol = TRUE, log = FALSE) * dzipois(Y, lambda = 0)
    
    lp[,1] <- po1 * P1 * dzipois(Y, lambda = 0)
    k1[,1] <- lf1[,1] * lf2[,1] * L2[,1] * lp[,1]
    
    ### For non-zeros
    for (l in 1:(k-1)) {
      
      for (i in 1:n) {
        
        lf1[i,l+1] <- po[i,l] * dmultinom(W1[i,], 1, prob = py11[l,], log = FALSE) #* y3[i,l+1]
        lf2[i,l+1] <- po[i,l] * dmultinom(W2[i,], 1, prob = py21[l,], log = FALSE) #* y3[i,l+1]
        
      }
      
      lp[,l+1] <- po[,l ] * dpois(Y, lambda = lam[,l], log = FALSE) #* y3[,l+1]
      lp1[,l] <- po[,l] * Pi[l]
      L2[,l+1]  <- po[,l] * mvnfast::dmvn(w, mean[l,], sigma[[l]], ncores = 8,
                                          isChol = TRUE, log = FALSE) #* y3[,l+1]
      
      k1[,l+1] <- lf1[,l+1] * lf2[,l+1] * lp[,l+1] * L2[,l+1] * lp1[,l]
      
    }
    
    for(l in 1:(k-1)){
      
      LogLike <- function(par) {
        
        lambda <- exp(x %*% par)
        LL <- -sum(po[,l] * dpois(Y, lambda, log = TRUE))
        return(LL)
        
      }
      
      par <- b[l,]
      
      m.like <- optim_sa(fun = LogLike, start = par, 
                         lower = rep(eps,ncol(x)),
                         upper = rep(1, ncol(x)), trace = TRUE,
                         control = list(t0 = 100, nlimit = 550, t_min = 0.1,
                                        dyn_rf = F, rf = 1, r = 0.7))
      b[l,] <- m.like$par
    }
    
    # M-Step Updates
    P1 <- sum(po1) / n
    
    mean1 <- colSums(po1 * w) / sum(po1)
    
    for (i in 1:n) sig1[,,i] <- po1[i] * outer((w[i,] - mean1),(w[i,] - mean1))
    for(l in 1:(k-1)) mean[l,] <- colSums(po[,l] * w) / colSums(po)[l]
    sigma1 <- apply(sig1, c(1, 2), sum)/ sum(po1) + diag(eps, ncol(w))
    
    for (l in 1:(k-1)) {
      
      for (i in 1:nrow(x)) {sig[,,i] <- po[i,l] * outer((w[i,] - mean[l,]),(w[i,] - mean[l,]))}
      sigma[[l]] <- apply(sig, c(1, 2), sum)/ sum(po[,l]) + diag(eps, ncol(w))
      
    }
    
    #if(count %% 12 == 0){
    
    py1 <- (colSums(po1 *  W1) / sum(po1)) #+ 0.2
    py2 <- (colSums(po1 *  W2) / sum(po1)) #+ 0.2
    
    for (l in 1:(k-1)) {
      
      py11[l, ] <- (colSums(po[,l] * W1) / colSums(po)[l])
      py21[l, ] <- (colSums(po[,l] * W2) / colSums(po)[l])
      
    }
    
    #}
    
    Pi <- colSums(po) / n
    
    # Compute Log_likelihood lki, lk2,
    
    L[count] <- sum(log(rowSums(k1)))
    a_k <- (L[count+1] - L[count]) / (L[count] - L[count-1])
    L[count + 2] <- L[count] + ((1-a_k)^-1 * (L[count+1] - L[count]))
    dif <- abs(L[count+2] - L[count+1])
    
    if (show_table) {
      
      #dif <- abs(L[count+2] - L[count+1])
      out_table = data.frame(Iteration = count, Likelihood = L[count+2], difference = dif)
      print(kable(out_table))
      #if (count == maxit || dif < tol) break;
      
    }

    if (count == maxit || dif < tol) break;
    count <- count + 1
    
  }
  
  Z <- apply(Poc, 1, which.max)
  Prof <- Poc1 <- Poc
  
  for(i in 1:n){
    for(j in 1:k){      
      Prof[i,j] <- ifelse(Poc[i,j] > 0.9, 1, 0)
      if(Poc[i,j] > 0){
        Poc1[i,j] <- log(Poc[i,j])
      }
    }
  }
  
  ai   <- -2*L[count + 2] - 2*m;
  ai3  <- -2*L[count + 2] - m*log(n)  
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  AIcc <-  ai - 2*m*(m+1)/(n-m-1)
  AIC3 <- -2*L[count+2] - 3*m
  AICu <-  AIcc - n*log(n/(n-m-1))
  ICL  <-  ai3 + suppressWarnings(sum(rowSums(Prof * Poc1)))
  Caic <- -2*L[count+2] - m*(1+log(n))
  AWE  <- -2*L[count+2] - 2*m*(3/2 + log(n))

  res <- list(
    "mean"           = mean, 
    "mean1"          = mean1, 
    "sigma1"         = sigma1, 
    "Prob1"          = P1,
    "prob2"          = Pi, 
    "post"           = Poc,
    "poiwei"         = b,
    "classification" = Z, 
    "logLik"         = L, 
    "AIC"            = ai, 
    "BIC"            = ai3,
    "sigma"          = sigma, 
    "ICL"            = ICL, 
    "AICc"           = AIcc, 
    "AIC3"           = AIC3, 
    "AICu"           = AICu, 
    "Caic"           = Caic, 
    "AWE"            = AWE)
  
  return(res)
  
}

