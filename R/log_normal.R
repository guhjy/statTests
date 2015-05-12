# Compute the p-value for testing the difference between means of two 
# log-normal random variables, and computes the confidence interval for
# the difference between the means.
#
# x: sample values of variable 1
# y: sample values of variable 2
# conf: desired confidence interval [0,1]
# max.shapiro.p: maximum allowable p-value for normality test of log(x)
#   and log(y), or NA if normality testing is to be skipped
# niter: number of iterations to use when generating the p-value
#
# Reference: Krishnamoorthy, K. and Mathew, T. (2003) "Inferences on the 
# means of lognormal distributions using generalized p-values and 
# generalized confidence intervals." Journal of Statistical Planning and 
# Inference, (2003), 115, 103 – 121
#
# The following code is adapted from the logtwo SAS macro by K. Krishnamoorthy
# (http://www.ucs.louisiana.edu/~kxk4695/statcalc/logtwo.sas). RNG seeding is
# omitted because it is unnecessary in R.
#
# Value: A list with
# * The sizes, means and standard-deviations of log(x) and log(y)
# * The p-values for the Shapiro normality tests of log(x) and log(y)
# * p.value: The p-value of the test for difference between the means of x and y, 
#   to the accuracy allowed by the number of interations specified; if the p-value 
#   was computed to be zero, the returned value is instead 1/niter and the 
#   'limit.reached' element is set to TRUE, indicating that the p-value is between 
#   [0,1/niter]
# * conf.interval: The confidence interval of the difference between the means of
#   x and y
lognorm.test <- function(x, y, conf=0.95, max.shapiro.p=0.1, niter=100000) {
    # log both variables
    logx <- log(x)
    logy <- log(y)
    
    # test for normality
    logx.norm.p <- shapiro.test(logx)$p.value
    logy.norm.p <- shapiro.test(logy)$p.value
    if (!is.na(max.shapiro.p)) {
        if (logx.norm.p < max.shapiro.p) {
            stop("Log(x) does not fit assumption of normality")
        }
        if (logy.norm.p < max.shapiro.p) {
            stop("Log(y) does not fit assumption of normality")
        }
    }
    
    # size
    n.x <- length(x)
    n.y <- length(y)
    # mean
    mean.logx <- mean(logx)
    mean.logy <- mean(logy)
    # variance
    var.logx <- var(logx)
    var.logy <- var(logy)
    # stdev
    sd.logx = sqrt(var.logx)
    sd.logy = sqrt(var.logy)
    
    # degrees of freedom
    df.x <- n.x - 1.0
    df.y <- n.y - 1.0
    
    # compute p-value
    p.value = 0.0
    gv <- rep(0, niter)
    # generate constants and random values ahead of time
    rnorm1 <- rnorm(niter) * (sd.logx / sqrt(n.x)) * sqrt(df.x)
    rnorm2 <- rnorm(niter) * (sd.logy / sqrt(n.y)) * sqrt(df.y)
    rgam1 <- 2.0 * rgamma(niter, df.x * 0.5)
    rgam2 <- 2.0 * rgamma(niter, df.y * 0.5)
    const.x = var.logx * df.x * 0.5
    const.y = var.logy * df.y * 0.5
    
    for (i in 1:niter) {
        tx <- exp(mean.logx + (rnorm1[i] / sqrt(rgam1[i])) + (const.x / rgam1[i]))
        ty <- exp(mean.logy + (rnorm2[i] / sqrt(rgam2[i])) + (const.y / rgam2[i]))
        gv[i] <- tx - ty
        if (gv[i] < 0.0) {
            p.value <- p.value + 1.0
        }
    }
    
    p.value <- p.value / niter
    limit.reached <- F
    if (p.value == 0.0) {
        p.value <- 1/niter
        limit.reached <- T
    }
    
    # confidence interval
    gv <- sort(gv)
    alpha <- 1.0 - conf
    ci.lower <- round(alpha * niter / 2.0)
    ci.upper <- round((1.0 - alpha / 2.0) * niter)
    
    list(n.x=n.x, n.y=n.y, mean.logx=mean.logx, mean.logy=mean.logy, sd.logx=sd.logx, sd.logy=sd.logy,
         shapiro.x.p.value=logx.norm.p, shapiro.y.p.value=logy.norm.p,
         p.value=p.value, niter=niter, limit.reached=limit.reached, conf.level=conf, 
         conf.interval=c(gv[ci.lower], gv[ci.upper]))
}