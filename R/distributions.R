densfuns <- list(cauchy=pcauchy, exponential=pexp, gamma=pgamma, geometric=pgeom, 
    lognormal=function(q, ...) pnorm(q, log.p=T, ...), logistic=plogis, normal=pnorm,
    Poisson=ppois, t=pt, weibull=pweibull, "negative binomial"=pnbinom)

# Tests the probability that `x` was drawn from one or more distributions.
dist.tests <- function(x, dists=c("cauchy", "exponential", "gamma", "geometric", "normal", 
        "Poisson", "t", "weibull", "lognormal", "logistic", "negative binomial"), inf=0.0001) {
    sapply(dists, function(d) {
        tryCatch({
            f <- fitdistr(x, d)
            l <- list(x, densfuns[[d]])
            for (n in names(f$estimate)) {
                l[[n]] <- f$estimate[n]
            }
            k=do.call(ks.test, l)
            list(params=f, D=k$statistic, p=k$p.value)
        }, error=function(e) e$message)
    })
}