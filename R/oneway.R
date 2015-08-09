#' Oneway Anova Table Generating
#' @param z A data frame or a list to be analyzed
#' @return The list of the anova table including the basic information
#' @export
#' @examples
#' oneway(list(c(3, 5, 3, 4, 6, 3), c("A", "A", "B", "A", "B", "B")))
oneway <- function(z, ...) UseMethod("oneway")
oneway.default <- function(z, ...){
  z <- as.list(z)
  x <- z[[1]]
  y <- as.factor(z[[2]])
  #grand mean
  M <- mean(x)
  #group means
  g.m <- by(x, y, mean)
  #group sizes
  g.n <- by(x, y, length)
  #total size
  n=length(x)
  #sum squares of total (here is the x)
  SSY <- sum((x-M)^2)
  #mean squared total
  MSY <- SSY/(n-1)
  #making a list of group means
  g.m <- predict(lm(x~y))
  #sum squres of within groups
  SSE <- sum((x-g.m)^2)
  #df = n-no.groups
  g <- length(levels(y))
  #mean squared within
  MSE <- SSE/(n-g)
  ##Of the total variation SSY, the part not explained by group is
  SSE/SSY
  #sum squares of among groups
  SSF <- SSY-SSE
  SSF/SSY
  #mean squared among groups
  MSF <- SSF/(g-1)
  g.m <- by(x, y, mean)
  #ANOVA via F distribution
  F.F <- MSF/MSE
  boxplot(x~y)
  by(x, y, var)
  #we need to know df for (group) factor and df for error
  dfF <- g-1
  dfE <- n-g
  dfTot <- n-1
  #thus, the P value from F is
  P.value.from.F <- 1-pf(F.F, dfF, dfE)
  Pr <- P.value.from.F
  Pr
  #alpha <- 0.05
  #decision <- function (x) {if (x < alpha) {paste(Pr,"*", collapse = ' ')}}
  #P_r <- decision(Pr)
  anova.tab <- list(group_size=g.n, group_mean=g.m,
                    Df=c(dfF, dfE, dfTot),
                    Sum_Sq=c(SSF, SSE, SSY),
                    Mean_Sq=c(MSF, MSE, MSY),
                    Rsq=c(SSF/SSY, NA, NA),
                    F.val = c(F.F, NA, NA),
                    P.val = c(Pr, NA, NA))
  class(anova.tab) <- c("oneway", class(anova.tab))
  anova.tab
}

#' Oneway Anova Table Generating in a General Way
#' @param x A numeric object
#' @param y A factor object
#' @return calling oneway function and return the same result as oneway
#' @export
oneway.factor <- function(x, y, ...) {
  #x is numerica variable and y is factor variable
  x <- as.numeric(x)
  y <- as.factor(y)
  z <- list(x,y)
  oneway.default(z)
}

#' Oneway Anova Table Generating for R models
#' @param data A list
#' @param formula An oject for the model.frame
#' @return calling oneway.factor and return the same result as oneway
#' @export
oneway.formula <- function(formula, data=list(), ...) {
  mf <- model.frame(formula = formula, data=data)
  x <- model.response(mf)
  y <- as.factor(data[[2]]) ##This because the input data is list
  oneway.factor(x, y) ##Here we call oneway.factor function,
}

#' Print Function to print out the information of oneway anova table
#' @param x A list from oneway function
#' @return the necessory information for a anova table
#' @export
print.oneway <- function(x, ...) {
  cat("Call:\n")
  cat("\nAnalysis of Variance Table\n")
  print(x$Df)
  cat("\nSum Sq :\n")
  print(x$Sum_Sq)
  cat("\nMean Sq :\n")
  print(x$Mean_Sq)
  cat("\nF value :\n")
  print(x$F.val)
  cat("\nPr(>F) :\n")
  print(x$P.val)
}

#' The summary method to create a summary object
#' @param object An object from oneway return
#' @return a table like list includng the anova table information
#' @export
summary.oneway <- function(object, ...) {
  tab <- cbind(df=object$Df,
               SumSq=object$Sum_Sq,
               MeanSq=object$Mean_Sq,
               Fvalue=object$F.val,
               Pvalue=object$P.val)
  res <- list(call=object$call, coefficients=tab)
  class(res) = "summary.oneway"
  res
}

#' The print functional to print anova table
#' @param x An object from oneway return
#' @return an anova table containing the basic information
#' @export
print.summary.oneway <- function(x, ...) {
  cat("Call:\n")
  cat("Analysis of Variance Table: \n")
  print(x$coefficients, row.names=c("treat", "Residuals"))
  cat("\nSignif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n")
}

#' Fisher LSD test
#' @param object An object from oneway function call
#' @param conf.level A number of confidence level
#' @return a list of the mean difference and Fisher LSD results
#' @export
lsmeans.oneway <- function(object, conf.level=0.95) {

  # A function to perform pairwise t-tests using Fisher's LSD procedure
  model <- object

  # Initialize variable
  MSE        <- model$Mean_Sq[[2]]

  # Define matrices
  diff  <- matrix(nrow=length(model$group_size),ncol=length(model$group_size))
  LSD   <- matrix(nrow=length(model$group_size),ncol=length(model$group_size))
  means <- matrix(nrow=length(model$group_size),ncol=length(model$group_size))


  # Loop through all pairs
  for( i in 1:(length(model$group_size)-1) ) {
    for( j in (i+1):length(model$group_size) ) {
      s1 <- names(model$group_size)[i]
      s2 <- names(model$group_size)[j]                     # Defines the second group

      #   s1 <- model$model[[1]][model$model[[2]]==sset1]    # The data in Group 1
      #   s2 <- model$model[[1]][model$model[[2]]==sset2]    # The data in Group 2

      # Calculate the LSD statistic
      t          <- qt(p=1-(1-conf.level)/2, df=model$Df[[2]] )
      LSD[i,j]   <- t*sqrt(MSE * (1/model$group_size[[i]] + 1/model$group_size[[j]]))
      LSD[j,i]   <- t*sqrt(MSE * (1/model$group_size[[j]] + 1/model$group_size[[i]]))

      # Calculate the difference in means between the two groups
      means[i,j] <- model$group_mean[[i]]-model$group_mean[[j]]
      means[j,i] <- model$group_mean[[j]]-model$group_mean[[i]]

      # Determine if the difference in means is statistically significant
      diff[i,j]  <- abs( model$group_mean[[i]]-model$group_mean[[j]]) > LSD[i, j]
      diff[j,i]  <- abs( model$group_mean[[j]]-model$group_mean[[i]]) > LSD[j, i]

    } # End j loop
  }   # End i loop

  # Tidy up the output a bit
  rownames(diff) <- names(model$group_size)
  colnames(diff) <- names(model$group_size )

  rownames(LSD) <- names(model$group_size)
  colnames(LSD) <- names(model$group_size)

  rownames(means) <- names(model$group_size)
  colnames(means) <- names(model$group_size)

  # Define the variable we will return from this function
  item       <- list(diff=diff)
  item$means <- means
  item$LSD   <- LSD

  # Return the variable
  return(item)

}

#' Coagulation Information.
#' A dataset containing the diet groups and the corresponding coag response.
#' @format A data frame with 24 rows and 2 variables:
#' \describe{
#'   \item{coag}{num, measured numbers}
#'   \item{diet}{factor, 4 levels "A", "B", "C", "D"}
#'   ...
#' }
"coagulation"
